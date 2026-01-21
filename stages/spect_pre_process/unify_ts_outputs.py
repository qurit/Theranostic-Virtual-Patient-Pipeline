import os
import json
import numpy as np
import nibabel as nib
from json_minify import json_minify


class TdtRoiUnifyStage:
    """
    Combines TotalSegmentator outputs (body / total / head_glands_cavities)
    into a single multilabel segmentation in TDT label space.

    Output labels: (might have more later)
      0 background
      1 body
      2 kidney
      3 liver
      4 prostate
      5 spleen
      6 heart
      7 salivary_glands
    """

    def __init__(self, config, context):
        self.config = config
        self.context = context

        subdir_name = config["subdir_names"]["spect_preprocessing"]
        output_root = config["output_folder"]["title"]
        self.output_dir = os.path.join(output_root, subdir_name)
        os.makedirs(self.output_dir, exist_ok=True)

        self.prefix = config["spect_preprocessing"]["name"]

        # class map json path
        self.ts_map_path = "/home/jhubadmin/Theranostic-Virtual-Patient-Pipeline/data/tdt_map.json"
        if not os.path.exists(self.ts_map_path):
            raise FileNotFoundError(f"Class map json not found: {self.ts_map_path}")

        with open(self.ts_map_path, encoding="utf-8") as f:
            self.ts_map_json = json.loads(json_minify(f.read()))

        # Convert json maps to name -> int_label for easy lookup
        self.total_name2id = {name: int(lab) for lab, name in self.ts_map_json["total"].items()}
        self.head_name2id  = {name: int(lab) for lab, name in self.ts_map_json["head_glands_cavities"].items()}
        self.tdt_name2id   = {name: int(lab) for lab, name in self.ts_map_json["TDT_Pipeline"].items()}

        self.ct_nii_path = context.ct_nii_path
        self.body_ml_path = context.body_ml_path
        self.total_ml_path = context.total_ml_path
        self.head_ml_path = context.head_glands_cavities_ml_path

        self.plan = context.totseg_plan
        if self.plan is None:
            raise ValueError("Missing context.totseg_plan; run TotalSegmentationStage first.")

    # -----------------------------
    # helpers
    # -----------------------------
    @staticmethod
    def _load_int_seg(path):
        arr = nib.load(path).get_fdata()
        return arr.astype(np.int16) # only need int labels

    def _assert_inputs_exist(self):
        if self.ct_nii_path is None or not os.path.exists(self.ct_nii_path):
            raise FileNotFoundError(f"CT not found: {self.ct_nii_path}")

        if self.body_ml_path is None or not os.path.exists(self.body_ml_path):
            raise FileNotFoundError(f"Body seg not found: {self.body_ml_path}")

        if self.plan.get("run_total", False):
            if self.total_ml_path is None or not os.path.exists(self.total_ml_path):
                raise FileNotFoundError(f"Total seg not found: {self.total_ml_path}")

        if self.plan.get("run_head_glands_cavities", False):
            if self.head_ml_path is None or not os.path.exists(self.head_ml_path):
                raise FileNotFoundError(f"Head seg not found: {self.head_ml_path}")
    
    def _create_roi_unified(self, body_seg, total_seg, head_seg):
        # Start output with background
        roi_unified = np.zeros(body_seg.shape, dtype=np.uint8)

        # Always paint body where body_seg>0
        roi_unified[body_seg > 0] = self.tdt_name2id["body"]

        requested = set(self.plan["tdt_roi_subset"])

        # TOTAL-task mapping
        if total_seg is not None:
            if "kidney" in requested:
                kL = self.total_name2id["kidney_left"]
                kR = self.total_name2id["kidney_right"]
                roi_unified[(total_seg == kL) | (total_seg == kR)] = self.tdt_name2id["kidney"]

            if "liver" in requested:
                roi_unified[total_seg == self.total_name2id["liver"]] = self.tdt_name2id["liver"]

            if "prostate" in requested:
                roi_unified[total_seg == self.total_name2id["prostate"]] = self.tdt_name2id["prostate"]

            if "spleen" in requested:
                roi_unified[total_seg == self.total_name2id["spleen"]] = self.tdt_name2id["spleen"]

            if "heart" in requested:
                roi_unified[total_seg == self.total_name2id["heart"]] = self.tdt_name2id["heart"]

        # HEAD-task mapping (salivary)
        if head_seg is not None and "salivary_glands" in requested:
            pL = self.head_name2id["parotid_gland_left"]
            pR = self.head_name2id["parotid_gland_right"]
            sL = self.head_name2id["submandibular_gland_left"]
            sR = self.head_name2id["submandibular_gland_right"]
            roi_unified[np.isin(head_seg, [pL, pR, sL, sR])] = self.tdt_name2id["salivary_glands"]
        
        return roi_unified
        

    # -----------------------------
    # main
    # -----------------------------
    def run(self):
        self._assert_inputs_exist()

        ct_nii = nib.load(self.ct_nii_path)

        body_seg = self._load_int_seg(self.body_ml_path)
        total_seg = self._load_int_seg(self.total_ml_path) if self.plan["run_total"] else None
        head_seg = self._load_int_seg(self.head_ml_path) if self.plan["run_head_glands_cavities"] else None

        roi_unified = self._create_roi_unified(body_seg, total_seg, head_seg)

        # Save NIfTI aligned to CT
        out_path = os.path.join(self.output_dir, f"{self.prefix}_tdt_roi_seg.nii.gz")
        out_img = nib.Nifti1Image(roi_unified.astype(np.uint8), ct_nii.affine, ct_nii.header)
        out_img.set_data_dtype(np.uint8)
        nib.save(out_img, out_path)

        # Update context
        self.context.tdt_roi_seg_path = out_path

        return self.context
