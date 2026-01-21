import os
import dicom2nifti
import torch
import SimpleITK as sitk
from torch.cuda.amp import GradScaler as _GradScaler
from totalsegmentator.python_api import totalsegmentator

torch.GradScaler = _GradScaler

TDT_ALLOWED_ROIS = {
    "body",
    "kidney",
    "liver",
    "prostate",
    "spleen",
    "heart",
    "salivary_glands",
}

# How each TDT ROI expands into task-specific roi_subset strings
TDT_TO_TOTSEG = {
    "kidney": ("total", ["kidney_left", "kidney_right"]),
    "liver": ("total", ["liver"]),
    "prostate": ("total", ["prostate"]),
    "spleen": ("total", ["spleen"]),
    "heart": ("total", ["heart"]),

    # “salivary_glands” = parotids + submandibulars (your head task)
    "salivary_glands": ("head_glands_cavities", [
        "parotid_gland_left",
        "parotid_gland_right",
        "submandibular_gland_left",
        "submandibular_gland_right",
    ]),

    # “body” lives in task="body" 
    "body": ("body", []),
}

class TotalSegmentationStage:
    def __init__(self, config, context):
        self.config = config
        self.context = context
        
        self.ct_nii_path = os.path.join(self.output_dir, f"{self.prefix}_ct.nii.gz")

        self.ct_input_path = config["ct_input"]["path1"]
        self.roi_subset = config["spect_preprocessing"]["roi_subset"]
        self.ml = True

        subdir_name = config["subdir_names"]["spect_preprocessing"]
        output_root = config["output_folder"]["title"]
        self.output_dir = os.path.join(output_root, subdir_name)
        os.makedirs(self.output_dir, exist_ok=True)

        self.prefix = config["spect_preprocessing"]["name"]

        self.ct_nii_path = None
        self.body_ml_path = None 
        self.head_glands_cavities_ml_path = None
        self.total_ml_path = None

    def _standardize_ct_to_nifti(self):
        if os.path.exists(self.ct_nii_path):
            return

        if os.path.isdir(self.ct_input_path):
            dicom2nifti.dicom_series_to_nifti(
                self.ct_input_path,
                self.ct_nii_path,
                reorient_nifti=True,
            )
        else:
            lower_input = self.ct_input_path.lower()
            if lower_input.endswith((".nii", ".nii.gz")):
                sitk.WriteImage(sitk.ReadImage(self.ct_input_path), self.ct_nii_path, True)
            else:
                raise ValueError(
                    "Unsupported CT input. Provide a DICOM folder or a NIfTI file "
                    f"(.nii/.nii.gz). Got: {self.ct_input_path}"
                )
    def _pre_totalsegmentation_checks(self):
        # -------- roi checks --------
        rois = self.roi_subset
        if isinstance(rois, str):
            rois = [rois]
        rois = [str(r).strip() for r in rois if str(r).strip()]

        if not rois or len(rois) == 0:
            raise ValueError(
                "roi_subset must contain at least one ROI from: "
                f"{sorted(TDT_ALLOWED_ROIS)}"
            )

        # -------- validate --------
        invalid = [r for r in rois if r not in TDT_ALLOWED_ROIS]
        if invalid:
            raise ValueError(
                f"Invalid ROI(s): {invalid}. Allowed: {sorted(TDT_ALLOWED_ROIS)}"
            )

        # -------- build plan --------
        run_body = True  # always run body if any ROI specified within TDT allowed

        total_rois = []
        head_rois = []
        seen_total = set()
        seen_head = set()

        for r in rois:
            task, expanded = TDT_TO_TOTSEG[r]

            if task == "total":
                for x in expanded:
                    if x not in seen_total:
                        total_rois.append(x)
                        seen_total.add(x)

            elif task == "head_glands_cavities":
                for x in expanded:
                    if x not in seen_head:
                        head_rois.append(x)
                        seen_head.add(x)


        plan = {
            "run_body": run_body,
            "run_total": bool(total_rois),
            "run_head_glands_cavities": bool(head_rois),
            "total_roi_subset": total_rois,     # TotSeg names
            "head_roi_subset": head_rois,       # TotSeg names
            "tdt_roi_subset": rois,             # User names
        }
        return plan
    
    def _files_exist(self):
        """all expected outputs"""
        
        # check if expected outputs exist - all files as ML is True
        self.body_ml_path = os.path.join(self.output_dir, f"{self.prefix}_body_ml.nii.gz")
        self.head_glands_cavities_ml_path = os.path.join(self.output_dir, f"{self.prefix}_head_glands_cavities_ml.nii.gz")
        self.total_ml_path  = os.path.join(self.output_dir, f"{self.prefix}_total_ml.nii.gz")

        body_ml_done = os.path.exists(self.body_ml_path)
        head_glands_cavities_ml_done = os.path.exists(self.head_glands_cavities_ml_path)
        total_ml_done = os.path.exists(self.total_ml_path)
        
        return body_ml_done, head_glands_cavities_ml_done, total_ml_done

    def run(self):
        self._standardize_ct_to_nifti()
        
        # --- validate user ROI list + compute which tasks to run ---
        plan = self._pre_totalsegmentation_checks()
        
        body_ml_done, head_glands_cavities_ml_done, total_ml_done = self._files_exist()
        
        if plan["run_body"] and not body_ml_done:
            print("Running TotalSegmentator for task : BODY...")
            totalsegmentator(
                self.ct_nii_path,
                self.body_ml_path,
                ml=True,
                task="body"
                )

        # --- TOTAL (organs etc.) ---
        if plan["run_total"] and not total_ml_done:
            print("Running TotalSegmentator for task : TOTAL...")
            totalsegmentator(
                self.ct_nii_path,
                self.total_ml_path,
                ml=True,
                task="total",
                roi_subset=plan["total_roi_subset"],
            )

        # --- HEAD GLANDS / CAVITIES ---
        if plan["run_head_glands_cavities"] and not head_glands_cavities_ml_done:
            print("Running TotalSegmentator for task : HEAD_GLANDS_CAVITIES...")
            totalsegmentator(
                self.ct_nii_path,
                self.head_glands_cavities_ml_path,
                ml=True,
                task="head_glands_cavities",
            )
            
        # final existence checks (only what was requested)
        body_done, head_done, total_done = self._files_exist()
        if plan["run_body"] and not body_done:
            raise FileNotFoundError(f"Body seg not found: {self.body_ml_path}")
        if plan["run_total"] and not total_done:
            raise FileNotFoundError(f"Total seg not found: {self.total_ml_path}")
        if plan["run_head_glands_cavities"] and not head_done:
            raise FileNotFoundError(f"Head glands seg not found: {self.head_glands_cavities_ml_path}")

        # --- update context ---
        self.context.ct_nii_path = self.ct_nii_path
        self.context.body_ml_path = self.body_ml_path
        self.context.total_ml_path = self.total_ml_path if plan["run_total"] else None
        self.context.head_glands_cavities_ml_path = (
            self.head_glands_cavities_ml_path if plan["run_head_glands_cavities"] else None
        )
        
        return self.context
