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

        # ---- ROI subset is now TDT-level names (kidney, liver, ...) ----
        self.roi_subset = config["spect_preprocessing"]["roi_subset"]  
        if isinstance(self.roi_subset, str):  
            self.roi_subset = [self.roi_subset]  
        self.roi_subset = [str(r).strip() for r in self.roi_subset if str(r).strip()]  

        # ---- Load ONLY TDT pipeline label map (name -> id) ----
        self.ts_map_path = "/home/jhubadmin/Theranostic-Virtual-Patient-Pipeline/data/tdt_map.json"
        if not os.path.exists(self.ts_map_path):
            raise FileNotFoundError(f"Class map json not found: {self.ts_map_path}")

        with open(self.ts_map_path, encoding="utf-8") as f:
            ts_map_json = json.loads(json_minify(f.read()))

        # name -> int_label, e.g. {"body":1,"kidney":2,...}
        self.tdt_name2id = {name: int(lab) for lab, name in ts_map_json["TDT_Pipeline"].items()}  

        self.ct_path = context.ct_path
        self.body_ml_path = context.body_ml_path
        self.tdt_roi_seg_path = context.tdt_roi_seg_path  # unified TDT ROI seg

    # -----------------------------
    # helpers
    # -----------------------------
    @staticmethod
    def _build_class_map(seg_arr, id_to_name): 
        """Return {roi_name: label_id} for labels present in seg_arr."""
        class_map = {}
        labels = np.unique(seg_arr.astype(int))
        for lab in labels:
            if lab == 0:
                continue
            name = id_to_name.get(int(lab))  
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
        mu_map *= body_seg_arr  # assumes body_seg_arr is 0/1-ish mask; if labels >0 that's fine too

        out_path = os.path.join(self.output_dir, filename)
        mu_map.tofile(out_path)
        return out_path

    def _filter_to_requested_rois(self, roi_seg_arr):  
        """
        Keep only labels requested in config (plus body if present), zero out the rest.
        This makes downstream masks/class_map reflect *only* user-selected ROIs.
        """
        requested = set(self.roi_subset)  
        keep_ids = set()  

        # Always keep body label if it exists in the unified seg
        if "body" in self.tdt_name2id:  
            keep_ids.add(self.tdt_name2id["body"])  

        for name in requested:  
            lab = self.tdt_name2id.get(name)  
            if lab is None:
                raise ValueError(f"Requested ROI '{name}' not in TDT label map.")  
            keep_ids.add(lab)  

        keep_ids.discard(0)  

        out = roi_seg_arr.copy()  
        out[~np.isin(out, list(keep_ids))] = 0  
        return out  

    @staticmethod
    def _to_simind_grid(nii_obj, resize=None, transpose_tuple=(2, 1, 0), zoom_order=0):
        """
        Convert NIfTI object to SIMIND grid format with optional resizing.

        NOTE: This assumes your SIMIND expects [z,y,x] with y-flip on the middle axis
        as you were doing before.
        """
        arr = np.array(nii_obj.get_fdata(dtype=np.float32))
        arr = np.transpose(arr, transpose_tuple)[:, ::-1, :]

        scale = 1.0
        if resize is not None:
            if arr.shape[1] != arr.shape[2]:
                raise ValueError("Resize parameter requires square in-plane dimensions (x=y).")
            scale = resize / arr.shape[1]
            arr = zoom(arr, (scale, scale, scale), order=zoom_order)

        return arr, scale

    # -----------------------------
    # main
    # -----------------------------
    def run(self):
        # ---- existence checks ----
        if self.ct_path is None or not os.path.exists(self.ct_path):
            raise FileNotFoundError(f"ct_path not found: {self.ct_path}")

        if self.body_ml_path is None or not os.path.exists(self.body_ml_path):
            raise FileNotFoundError(f"Body segmentation not found: {self.body_ml_path}")

        if self.tdt_roi_seg_path is None or not os.path.exists(self.tdt_roi_seg_path):
            raise FileNotFoundError(f"Unified TDT ROI seg not found: {self.tdt_roi_seg_path}")  

        # ---- load nifti ----
        ct_nii = nib.load(self.ct_path)
        body_nii = nib.load(self.body_ml_path)
        roi_nii = nib.load(self.tdt_roi_seg_path)  

        # ---- to simind grid ----
        ct_arr, scale = self._to_simind_grid(ct_nii, resize=self.resize, zoom_order=1)       # CT: linear
        body_arr, _ = self._to_simind_grid(body_nii, resize=self.resize, zoom_order=0)       # mask: body (nearest)
        roi_arr, _ = self._to_simind_grid(roi_nii, resize=self.resize, zoom_order=0)         # mask: unified ROI (nearest)

        # ---- sanity: body should be binary-ish ----
        body_mask = (body_arr > 0).astype(np.float32)  
        roi_arr = roi_arr.astype(np.int16) 

        # ---- keep only requested TDT ROIs ----
        roi_arr = self._filter_to_requested_rois(roi_arr) 

        # ---- build masks + class map from TDT label space ----
        masks = self._build_label_masks(roi_arr)  
        id_to_name = {v: k for k, v in self.tdt_name2id.items()}  
        class_seg = self._build_class_map(roi_arr, id_to_name)  

        # ---- spacing ----
        zooms_mm = np.array(ct_nii.header.get_zooms()[:3], dtype=float) / scale
        zooms_mm = zooms_mm[[2, 1, 0]]  # (z,y,x)
        arr_px_spacing_cm = tuple(float(x) * 0.1 for x in zooms_mm)

        pixel_size_cm = (arr_px_spacing_cm[1] + arr_px_spacing_cm[2]) / 2.0

        atn_av_path = self._write_attenuation_bin(
            ct_arr,
            body_mask,  # use binary mask
            pixel_size_cm=pixel_size_cm,
        )

        # ---- write SIMIND bins ----
        roi_bin = os.path.join(self.output_dir, f"{self.prefix}_roi_seg.bin")
        body_bin = os.path.join(self.output_dir, f"{self.prefix}_body_seg.bin")
        roi_body_bin = os.path.join(self.output_dir, f"{self.prefix}_roi_body_seg.bin")

        roi_arr.astype(np.float32).tofile(roi_bin)
        body_mask.astype(np.float32).tofile(body_bin)  
        roi_arr.astype(np.float32).tofile(roi_body_bin)  
        
        # ---- update context ----
        self.context.ct_arr = ct_arr
        self.context.roi_seg_arr = roi_arr
        self.context.body_seg_arr = body_mask  
        self.context.roi_body_seg_arr = roi_arr 

        self.context.mask_roi_body = masks  
        self.context.class_seg = class_seg
        self.context.roi_subset = list(self.roi_subset)  

        self.context.atn_av_path = atn_av_path
        self.context.roi_seg_bin_path = roi_bin
        self.context.body_seg_bin_path = body_bin
        self.context.roi_body_seg_bin_path = roi_body_bin

        self.context.arr_px_spacing_cm = arr_px_spacing_cm

        self.context.extras["preprocess_scale_factor"] = float(scale)
        self.context.extras["mu_pixel_size_cm_used"] = float(pixel_size_cm)
        self.context.extras["tdt_label_map"] = dict(self.tdt_name2id)  

        return self.context
