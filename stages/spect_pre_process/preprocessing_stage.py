import os
import nibabel as nib
import numpy as np
from scipy.ndimage import zoom
import json
from json_minify import json_minify

class SimindPreprocessStage:
    def __init__(self, config, context):
        self.config = config
        self.context = context

        subdir_name = config["subdir_names"]["spect_preprocessing"]
        output_root = config["output_folder"]["title"]
        self.output_dir = os.path.join(output_root, subdir_name)
        os.makedirs(self.output_dir, exist_ok=True)
 
        self.prefix = config["spect_preprocessing"]["name"]
        self.resize = config["spect_preprocessing"]["xy_dim"]
        self.roi_subset = list(config["spect_preprocessing"]["roi_subset"])
        
        # load total segmentator classes
        self.ts_classes_path = config["spect_preprocessing"]["total_segmentator_classes_path"]
        if not os.path.exists(self.ts_classes_path):
            raise FileNotFoundError(f"Total seg class file not found: {self.ts_classes_path}")
        
        with open(self.ts_classes_path, encoding="utf-8") as f:
            self.ts_classes_json = json.loads(json_minify(f.read()))
            
        self.ts_classes = self.ts_classes_json["total"]
        # Reverse map: class name -> label id (int), more convient method
        self.ts_name_to_id = {name: int(lab) for lab, name in self.ts_classes.items()}


        self.ct_nii_path = context.ct_nii_path
        self.roi_seg_path = context.roi_seg_path
        self.body_seg_dir = getattr(context, "_body_seg_dir", None)

        if self.body_seg_dir is None:
            raise AttributeError("Context is missing _body_seg_dir.")

    # -----------------------------
    # helpers
    # -----------------------------
    @staticmethod
    def _build_class_map(seg_arr, classes):
        class_map = {}
        labels = np.unique(seg_arr.astype(int))
        for lab in labels:
            if lab == 0:
                continue
            name = classes.get(str(lab))
            if name is not None:
                class_map[name] = int(lab)
        return class_map

    @staticmethod
    def _build_label_masks(arr):
        labels = np.unique(arr)
        labels = labels[labels != 0]
        
        if labels.size == 0:
            raise ValueError(
                "Segmentation has no non-zero labels (only background=0). "
                "Segmentation likely failed or ROI subset is empty/mismatched."
                )
        return {int(lab): (arr == lab) for lab in labels}
    

    @staticmethod
    def _hu_to_mu(hu_arr, pixel_size_cm, mu_water=0.1537, mu_bone=0.2234):
        mu_water_pixel = mu_water * pixel_size_cm
        mu_bone_pixel = mu_bone * pixel_size_cm

        mu_map = np.zeros_like(hu_arr, dtype=np.float32)

        soft = hu_arr <= 0
        bone = hu_arr > 0

        mu_map[soft] = mu_water_pixel * (1 + hu_arr[soft] / 1000.0)
        mu_map[bone] = mu_water_pixel + (hu_arr[bone] / 1000.0) * (mu_bone_pixel - mu_water_pixel)

        return mu_map

    def _write_attenuation_bin(self, ct_arr, body_seg_arr, pixel_size_cm,
                               filename="spect_preprocessing_atn_av.bin"):
        mu_map = self._hu_to_mu(np.asarray(ct_arr, dtype=np.float32), pixel_size_cm)
        mu_map *= body_seg_arr

        out_path = os.path.join(self.output_dir, filename)
        mu_map.tofile(out_path)
        return out_path

    def _combine_roi_and_body(self, roi_seg_arr, body_seg_arr):
        body_label = self.ts_name_to_id.get("body")
        if body_label is None:
            raise KeyError("'body' not found in total segmentator classes json")
        roi_mask = roi_seg_arr != 0
        body_mask = body_seg_arr > 0

        out = np.zeros_like(roi_seg_arr, dtype=np.float32)
        out[body_mask] = body_label
        out[roi_mask] = roi_seg_arr[roi_mask]

        masks = self._build_label_masks(out)
        return out, masks

    @staticmethod
    def _to_simind_grid(nii_obj, resize=None, transpose_tuple=(2, 1, 0),zoom_order=0):
        """Convert NIfTI object to SIMIND grid format with optional resizing.
        
        xy will be resized based on resize parameter, z will be scaled accordingly based
        on ratio of xy value to resize parameter.
        assumes xy have shape (almost always true->James)
        """
        arr = np.array(nii_obj.get_fdata(dtype=np.float32))
        arr = np.transpose(arr, transpose_tuple)[:, ::-1, :]

        scale = 1.0
        if resize is not None:
            if arr.shape[1] != arr.shape[2]:
                raise ValueError("Resize parameter requires square in-plane dimensions (x=y).")
            scale = resize / arr.shape[1]
            arr = zoom(arr, (scale, scale, scale), order = zoom_order) 

        return arr, scale

    def _merge_kidneys_if_needed(self,roi_subset, roi_seg_arr):
        has_left = "kidney_left" in roi_subset
        has_right = "kidney_right" in roi_subset

        if has_left and not has_right:
            raise ValueError("PBPK requires both kidneys (missing kidney_right).")
        if has_right and not has_left:
            raise ValueError("PBPK requires both kidneys (missing kidney_left).")

        if has_left and has_right:
            left_label = self.ts_name_to_id.get("kidney_left")
            right_label = self.ts_name_to_id.get("kidney_right")
            merged_label = self.ts_name_to_id.get("kidney")
            new_subset = [r for r in roi_subset if r not in ("kidney_left", "kidney_right")]
            if "kidney" not in new_subset:
                new_subset.append("kidney")

            merged = roi_seg_arr.copy()
            merged[np.isin(roi_seg_arr, (left_label, right_label))] = merged_label
            return merged, new_subset

        return roi_seg_arr, list(roi_subset)

    # -----------------------------
    # main
    # -----------------------------
    def run(self):
        if self.ct_nii_path is None or not os.path.exists(self.ct_nii_path):
            raise FileNotFoundError(f"ct_nii_path not found: {self.ct_nii_path}")

        if self.roi_seg_path is None or not os.path.exists(self.roi_seg_path):
            raise FileNotFoundError(f"roi_seg_path not found: {self.roi_seg_path}")

        if not os.path.isdir(self.body_seg_dir):
            raise NotADirectoryError(f"_body_seg_dir is not a directory: {self.body_seg_dir}")

        body_seg_path = os.path.join(self.body_seg_dir, "body.nii.gz")
        if not os.path.exists(body_seg_path):
            raise FileNotFoundError(f"Body segmentation not found: {body_seg_path}")

        ct_nii = nib.load(self.ct_nii_path)
        roi_nii = nib.load(self.roi_seg_path)
        body_nii = nib.load(body_seg_path)

        ct_arr, scale = self._to_simind_grid(ct_nii, resize=self.resize, zoom_order=1)   # CT: linear
        roi_arr, _    = self._to_simind_grid(roi_nii, resize=self.resize, zoom_order=0)  # mask: nearest
        body_arr, _   = self._to_simind_grid(body_nii, resize=self.resize, zoom_order=0) # mask: nearest

        roi_arr, roi_subset = self._merge_kidneys_if_needed(self.roi_subset, roi_arr)

        roi_body_arr, mask_roi_body = self._combine_roi_and_body(roi_arr, body_arr)
        class_seg = self._build_class_map(roi_body_arr, self.ts_classes)

        # spacing: header zooms are (x,y,z); our array is (z,y,x)
        zooms_mm = np.array(ct_nii.header.get_zooms()[:3], dtype=float) / scale
        zooms_mm = zooms_mm[[2, 1, 0]]
        arr_px_spacing_cm = tuple(float(x) * 0.1 for x in zooms_mm)

        # use in-plane spacing for mu scaling (scalar approximation)
        pixel_size_cm = (arr_px_spacing_cm[1] + arr_px_spacing_cm[2]) / 2.0

        atn_av_path = self._write_attenuation_bin(
            ct_arr,
            body_arr,
            pixel_size_cm=pixel_size_cm,
        )

        roi_bin = os.path.join(self.output_dir, f"{self.prefix}_roi_seg.bin")
        body_bin = os.path.join(self.output_dir, f"{self.prefix}_body_seg.bin")
        roi_body_bin = os.path.join(self.output_dir, f"{self.prefix}_roi_body_seg.bin")

        roi_arr.astype(np.float32).tofile(roi_bin)
        body_arr.astype(np.float32).tofile(body_bin)
        roi_body_arr.astype(np.float32).tofile(roi_body_bin)

        self.context.ct_arr = ct_arr
        self.context.roi_seg_arr = roi_arr
        self.context.body_seg_arr = body_arr
        self.context.roi_body_seg_arr = roi_body_arr

        self.context.mask_roi_body = mask_roi_body
        self.context.class_seg = class_seg
        self.context.roi_subset = roi_subset

        self.context.atn_av_path = atn_av_path
        self.context.roi_seg_bin_path = roi_bin
        self.context.body_seg_bin_path = body_bin
        self.context.roi_body_seg_bin_path = roi_body_bin

        self.context.arr_px_spacing_cm = arr_px_spacing_cm
        
        self.context.extras["preprocess_scale_factor"] = scale
        self.context.extras["mu_pixel_size_cm_used"] = float(pixel_size_cm)


        return self.context
