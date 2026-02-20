"""
TotalSegmentator-based segmentation stage for the TDT pipeline. 

This stage:
- Standardizes the CT input into a NIfTI file in the stage output directory.
- Runs TotalSegmentator for the required task(s) based on a user-facing ROI list.
- Writes ML output masks (NIfTI) for each task and stores paths + plan in `context`.

Expected Context interface
--------------------------
The incoming `context` object is expected to provide:
- context.ct_input_path : str
- context.config : dict (must include config["spect_preprocessing"]["roi_subset"] and ["name"])
- context.subdir_paths : dict (must include subdir_paths["spect_preprocessing"])

On success, this stage sets:
- context.ct_nii_path : str
- context.body_ml_path : str
- context.total_ml_path : Optional[str]
- context.head_glands_cavities_ml_path : Optional[str]
- context.totseg_plan : dict (plan of tasks/roi subsets executed)

Maintainer / contact: pyazdi@bccrc.ca  
""" 

from __future__ import annotations  

import os
from typing import Any, Dict, List, Literal, Optional, Sequence, Tuple, TypedDict, Union  

import dicom2nifti
import torch
import SimpleITK as sitk
from torch.cuda.amp import GradScaler as _GradScaler
from totalsegmentator.python_api import totalsegmentator

# -------------------------------------------------------------------------
# Compatibility note:
# Some TotalSegmentator/torch combinations expect `torch.GradScaler` to exist.
# This alias helps avoid runtime attribute errors in those environments.
# -------------------------------------------------------------------------
torch.GradScaler = _GradScaler  


# User-facing ROI names supported by this pipeline stage.
TDT_ALLOWED_ROIS = {
    "body",
    "kidney",
    "liver",
    "prostate",
    "spleen",
    "heart",
    "salivary_glands",
}

# How each TDT ROI expands into (TotalSegmentator task, roi_subset entries).
TDT_TO_TOTSEG = {
    "kidney": ("total", ["kidney_left", "kidney_right"]),
    "liver": ("total", ["liver"]),
    "prostate": ("total", ["prostate"]),
    "spleen": ("total", ["spleen"]),
    "heart": ("total", ["heart"]),

    # “salivary_glands” = parotids + submandibulars (head task)
    "salivary_glands": (
        "head_glands_cavities",
        [
            "parotid_gland_left",
            "parotid_gland_right",
            "submandibular_gland_left",
            "submandibular_gland_right",
        ],
    ),

    # “body” lives in task="body"
    "body": ("body", []),
}

CTInputType = Literal["nii", "dicom"] 

class TotSegPlan(TypedDict):  
    """Execution plan describing which TotalSegmentator tasks will run. """  

    run_body: bool
    run_total: bool
    run_head_glands_cavities: bool
    total_roi_subset: List[str]
    head_roi_subset: List[str]
    tdt_roi_subset: List[str]


class TotalSegmentationStage:
    """
    TDT Stage: TotalSegmentator segmentation.

    Parameters
    ----------
    context : Context-like
        Pipeline context object. Must provide `ct_input_path`, `config`, and `subdir_paths`.

    Notes
    -----
    - The stage always runs the `body` task if any ROI is requested (as written in original logic).
    - The `total` and `head_glands_cavities` tasks run only if required by the requested ROIs.
    """  
    
    def __init__(self, context: Any) -> None:  
        self.context = context

        self.ct_input_path: str = context.ct_input_path 

        # `roi_subset` in your config is in TDT ROI-space (e.g., "kidney", "salivary_glands").
        self.roi_subset: Union[str, Sequence[str]] = context.config["spect_preprocessing"]["roi_subset"]  
        self.ml: bool = True 

        self.output_dir: str = context.subdir_paths["spect_preprocessing"]  
        os.makedirs(self.output_dir, exist_ok=True)

        self.prefix: str = context.config["spect_preprocessing"]["name"]  

        self.ct_nii_path: Optional[str] = None  
        self.body_ml_path: Optional[str] = None  
        self.head_glands_cavities_ml_path: Optional[str] = None  
        self.total_ml_path: Optional[str] = None 

    def _standardize_ct_to_nifti(self) -> None:  
        """
        Convert/standardize the CT input into a NIfTI file for downstream tools.

        Output
        ------
        Writes: <output_dir>/<prefix>_ct.nii.gz

        Behavior
        --------
        - If input is a DICOM directory: converts series -> NIfTI (reoriented).
        - If input is an existing NIfTI (.nii/.nii.gz): re-writes into output_dir to standardize.
        - If output already exists, no-op.
        """  
        self.ct_nii_path = os.path.join(self.output_dir, f"{self.prefix}_ct.nii.gz")
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
                # Copy/write into the stage directory for consistent downstream paths.
                sitk.WriteImage(sitk.ReadImage(self.ct_input_path), self.ct_nii_path, True)
            else:
                raise ValueError(
                    "Unsupported CT input. Provide a DICOM folder or a NIfTI file "
                    f"(.nii/.nii.gz). Got: {self.ct_input_path}"
                )

    def _pre_totalsegmentation_checks(self) -> TotSegPlan:  
        """
        Validate the user ROI list and compute which TotalSegmentator tasks to run.

        Returns
        -------
        TotSegPlan
            Plan describing which tasks run and which `roi_subset` items are requested for each.

        Raises
        ------
        ValueError
            If ROI list is empty or contains unsupported ROI names.
        """  
        # -------- normalize rois --------
        rois = self.roi_subset
        if isinstance(rois, str):
            rois = [rois]
        rois = [str(r).strip() for r in rois if str(r).strip()]

        if not rois:
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
        run_body = True # body task runs if any ROI is requested.

        total_rois: List[str] = []
        head_rois: List[str] = []
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

        plan: TotSegPlan = {  
            "run_body": run_body,
            "run_total": bool(total_rois),
            "run_head_glands_cavities": bool(head_rois),
            "total_roi_subset": total_rois,   # TotalSegmentator ROI names
            "head_roi_subset": head_rois,     # TotalSegmentator ROI names
            "tdt_roi_subset": rois,           # User-facing ROI names
        }
        return plan

    def _files_exist(self) -> Tuple[bool, bool, bool]:  
        """
        Check whether expected output mask files already exist.

        Returns
        -------
        tuple[bool, bool, bool]
            (body_ml_done, head_glands_cavities_ml_done, total_ml_done)
        """  
        # All expected outputs exist as ML masks (ml=True).
        self.body_ml_path = os.path.join(self.output_dir, f"{self.prefix}_body_ml.nii.gz")
        self.head_glands_cavities_ml_path = os.path.join(
            self.output_dir, f"{self.prefix}_head_glands_cavities_ml.nii.gz"
        )
        self.total_ml_path = os.path.join(self.output_dir, f"{self.prefix}_total_ml.nii.gz")

        body_ml_done = os.path.exists(self.body_ml_path)
        head_glands_cavities_ml_done = os.path.exists(self.head_glands_cavities_ml_path)
        total_ml_done = os.path.exists(self.total_ml_path)

        return body_ml_done, head_glands_cavities_ml_done, total_ml_done

    def run(self) -> Any:  
        """
        Run the TotalSegmentator stage.

        Returns
        -------
        context : Context-like
            The updated context object (same instance as `self.context`).
        """ 
        self._standardize_ct_to_nifti() # standardize CT input to NIfTI for TotalSegmentator
        assert self.ct_nii_path is not None  # for type checker; will raise if standardization failed.

        # --- validate user ROI list + compute which tasks to run ---
        plan = self._pre_totalsegmentation_checks()

        body_ml_done, head_glands_cavities_ml_done, total_ml_done = self._files_exist()

        # --- BODY ---
        if plan["run_body"] and not body_ml_done:
            print("Running TotalSegmentator for task: BODY...")
            totalsegmentator(
                self.ct_nii_path,
                self.body_ml_path,
                ml=self.ml,  
                task="body",
            )

        # --- TOTAL (organs etc.) ---
        if plan["run_total"] and not total_ml_done:
            print("Running TotalSegmentator for task: TOTAL...")
            totalsegmentator(
                self.ct_nii_path,
                self.total_ml_path,
                ml=self.ml,  
                task="total",
                roi_subset=plan["total_roi_subset"],
            )

        # --- HEAD GLANDS / CAVITIES ---
        if plan["run_head_glands_cavities"] and not head_glands_cavities_ml_done:
            print("Running TotalSegmentator for task: HEAD_GLANDS_CAVITIES...")
            totalsegmentator(
                self.ct_nii_path,
                self.head_glands_cavities_ml_path,
                ml=self.ml, 
                task="head_glands_cavities",
            )

        # Final existence checks (only what was requested)
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
        self.context.totseg_plan = plan

        return self.context