"""
TDT ROI Unification Stage (TotalSegmentator -> TDT multilabel segmentation).  

This stage combines the TotalSegmentator outputs:
- task="body"
- task="total"
- task="head_glands_cavities"

into a single multilabel NIfTI volume in the TDT pipeline label space.

Key behaviors
-------------
- Uses `data/tdt_map.json` to translate TotalSegmentator class IDs -> ROI names -> TDT IDs.
- Paints a single output volume (`*_tdt_roi_seg.nii.gz`) aligned to the CT NIfTI (affine/header).
- Only maps ROIs that were requested in the TotalSegmentator plan (`context.totseg_plan`).

Expected Context interface
--------------------------
The incoming `context` object is expected to provide:
- context.subdir_paths : dict[str, str] with key "spect_preprocessing"
- context.config : dict with:
    - config["spect_preprocessing"]["name"] : str
- context.ct_nii_path : str (path to CT NIfTI created by TotalSegmentationStage)
- context.body_ml_path : str (path to body mask from TotalSegmentator)
- context.total_ml_path : Optional[str] (path to total mask; required if plan.run_total)
- context.head_glands_cavities_ml_path : Optional[str] (path to head mask; required if plan.run_head_glands_cavities)
- context.totseg_plan : dict with at least:
    - "run_total" : bool
    - "run_head_glands_cavities" : bool
    - "tdt_roi_subset" : list[str]

On success, this stage sets:
- context.tdt_roi_seg_path : str (path to unified multilabel NIfTI)

Maintainer / contact: pyazdi@bccrc.ca  
"""  
from __future__ import annotations  

import os
import json
from typing import Any, Dict, List, Optional, Sequence, TypedDict  

import numpy as np
import nibabel as nib
from json_minify import json_minify


class TotSegPlan(TypedDict, total=False):  
    """
    Minimal typing for the TotalSegmentator execution plan.

    Note: `total=False` so we can tolerate additional keys without typing friction.
    """  

    run_total: bool
    run_head_glands_cavities: bool
    tdt_roi_subset: List[str]


class TdtRoiUnifyStage:
    """
    Combine TotalSegmentator outputs into a single multilabel segmentation in TDT label space.

    Notes
    -----
    The label mapping is defined externally in `data/tdt_map.json`:
    - "total": maps TotalSegmentator *label id* -> *ROI name*
    - "head_glands_cavities": maps TotalSegmentator *label id* -> *ROI name*
    - "TDT_Pipeline": maps TDT *label id* -> *TDT ROI name*

    The unified output uses the TDT_Pipeline label IDs (uint8).
    """  

    def __init__(self, context: Any) -> None:  
        self.context = context

        self.output_dir: str = context.subdir_paths["spect_preprocessing"]  
        os.makedirs(self.output_dir, exist_ok=True)

        self.prefix: str = context.config["spect_preprocessing"]["name"]  

        # Locate mapping JSON relative to repository layout.
        # This assumes the stage file lives at:
        #   <repo_root>/stages/spect_pre_process/unify_ts_outputs.py
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")) 
        self.ts_map_path: str = os.path.join(repo_root, "data", "tdt_map.json")  
        if not os.path.exists(self.ts_map_path):
            raise FileNotFoundError(f"Class map json not found: {self.ts_map_path}")

        with open(self.ts_map_path, encoding="utf-8") as f:
            self.ts_map_json: Dict[str, Dict[str, str]] = json.loads(json_minify(f.read()))  

        # Convert json maps to name -> int_label for easy lookup
        # JSON structure: { "<int_id_as_str>": "<roi_name>" }
        self.total_name2id: Dict[str, int] = {  
            name: int(lab) for lab, name in self.ts_map_json["total"].items()
        }
        self.head_name2id: Dict[str, int] = {  
            name: int(lab) for lab, name in self.ts_map_json["head_glands_cavities"].items()
        }
        self.tdt_name2id: Dict[str, int] = {  
            name: int(lab) for lab, name in self.ts_map_json["TDT_Pipeline"].items()
        }

        # Input paths produced by TotalSegmentationStage
        self.ct_nii_path: Optional[str] = context.ct_nii_path  
        self.body_ml_path: Optional[str] = context.body_ml_path 
        self.total_ml_path: Optional[str] = context.total_ml_path  
        self.head_ml_path: Optional[str] = context.head_glands_cavities_ml_path  

        self.plan: TotSegPlan = context.totseg_plan  
        if self.plan is None:
            raise ValueError("Missing context.totseg_plan; run TotalSegmentationStage first.")

    # -----------------------------
    # helpers
    # -----------------------------
    @staticmethod
    def _load_int_seg(path: str) -> np.ndarray:  
        """
        Load a NIfTI segmentation volume and return it as a small integer array.

        Parameters
        ----------
        path : str
            Path to the NIfTI segmentation.

        Returns
        -------
        np.ndarray
            Segmentation array cast to int16 (sufficient for label IDs).
        """  
        arr = nib.load(path).get_fdata()
        return arr.astype(np.int16)  # only need int labels

    def _assert_inputs_exist(self) -> None:  
        """
        Validate that all required inputs exist on disk based on the plan.

        Raises
        ------
        FileNotFoundError
            If the CT or required segmentation masks do not exist.
        """  
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

    def _create_roi_unified(
        self,
        body_seg: np.ndarray,  
        total_seg: Optional[np.ndarray],  
        head_seg: Optional[np.ndarray], 
    ) -> np.ndarray:  
        """
        Create the unified TDT multilabel ROI volume.

        Parameters
        ----------
        body_seg : np.ndarray
            Body segmentation array (TotalSegmentator task="body").
        total_seg : Optional[np.ndarray]
            Total segmentation array (task="total"), or None if not executed.
        head_seg : Optional[np.ndarray]
            Head glands/cavities segmentation array, or None if not executed.

        Returns
        -------
        np.ndarray
            Unified multilabel volume in TDT label space (dtype uint8).
        """  
        # Start output with background
        roi_unified = np.zeros(body_seg.shape, dtype=np.uint8)

        # Always paint body where body_seg>0
        roi_unified[body_seg > 0] = self.tdt_name2id["body"]

        requested = set(self.plan["tdt_roi_subset"])

        # TOTAL-task mapping
        if total_seg is not None:
            if total_seg.shape != body_seg.shape:  
                raise ValueError(  
                    "Shape mismatch between body and total segmentation: "
                    f"body={body_seg.shape}, total={total_seg.shape}"
                )

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

        # HEAD-task mapping (salivary glands)
        if head_seg is not None and "salivary_glands" in requested:
            if head_seg.shape != body_seg.shape:  
                raise ValueError(  
                    "Shape mismatch between body and head segmentation: "
                    f"body={body_seg.shape}, head={head_seg.shape}"
                )

            pL = self.head_name2id["parotid_gland_left"]
            pR = self.head_name2id["parotid_gland_right"]
            sL = self.head_name2id["submandibular_gland_left"]
            sR = self.head_name2id["submandibular_gland_right"]
            roi_unified[np.isin(head_seg, [pL, pR, sL, sR])] = self.tdt_name2id["salivary_glands"]

        return roi_unified

    # -----------------------------
    # main
    # -----------------------------
    def run(self) -> Any:  
        """
        Run ROI unification and write the unified segmentation NIfTI.

        Returns
        -------
        context : Context-like
            Updated context object with `tdt_roi_seg_path` set.
        """  
        self._assert_inputs_exist()

        ct_nii = nib.load(self.ct_nii_path)  

        body_seg = self._load_int_seg(self.body_ml_path)  
        total_seg = self._load_int_seg(self.total_ml_path) if self.plan.get("run_total", False) else None 
        head_seg = (
            self._load_int_seg(self.head_ml_path) if self.plan.get("run_head_glands_cavities", False) else None
        )  

        roi_unified = self._create_roi_unified(body_seg, total_seg, head_seg)

        # Save NIfTI aligned to CT
        out_path = os.path.join(self.output_dir, f"{self.prefix}_tdt_roi_seg.nii.gz")
        out_img = nib.Nifti1Image(roi_unified.astype(np.uint8), ct_nii.affine, ct_nii.header)
        out_img.set_data_dtype(np.uint8)
        nib.save(out_img, out_path)

        # Update context
        self.context.tdt_roi_seg_path = out_path

        return self.context