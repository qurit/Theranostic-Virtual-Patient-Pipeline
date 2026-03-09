"""
Pipeline Context container for the Theranostic Digital Twin (TDT) pipeline.

`Context` is a lightweight, mutable object used to pass configuration, runtime metadata,
and intermediate outputs between pipeline stages.

Maintainer / contact: pyazdi@bccrc.ca
"""

from __future__ import annotations

from typing import Any, Dict, Optional


class Context:
    """
    Shared pipeline state passed between stages.
    """

    def __init__(self, logger: Optional[Any] = None) -> None:
        # Internal / compatibility
        self._logger: Optional[Any] = logger
        self._log_enabled: bool = True

        # Public misc storage for debugging / provenance / stage metadata
        self.extras: Dict[str, Any] = {}

        # ----------------------------- Initial setup -----------------------------
        self.mode: Optional[str] = None
        self.ct_input_path: Optional[str] = None
        self.ct_input_type: Optional[str] = None
        self.ct_indx: Optional[int] = None
        self.output_folder_path: Optional[str] = None
        self.subdir_paths: Optional[Dict[str, str]] = None
        self.subdir_names: Optional[Dict[str, str]] = None  
        self.synthetic_lesions_enabled: Optional[bool] = None  

        # Config fields (snapshots / sections)
        self.config: Dict[str, Any] = {}  # entire config dict for stages to access as needed

        # ----------------------------- Phase 1: Digital Twin -----------------------------  
        # Stage 1.1: TotalSegmentator  
        self.ct_nii_path: Optional[str] = None
        self.body_ml_path: Optional[str] = None
        self.head_glands_cavities_ml_path: Optional[str] = None
        self.total_ml_path: Optional[str] = None
        self.totseg_plan: Optional[Dict[str, Any]] = None

        # Stage 1.2: ROI Unification  
        self.tdt_roi_seg_path: Optional[str] = None

        # Stage 1.3: Synthetic Lesions  
        self.synthetic_lesions_outdir: Optional[str] = None
        self.synthetic_lesions_results: Optional[Dict[str, Any]] = None
        self.synthetic_lesions_backup_seg_path: Optional[str] = None
        self.synthetic_lesions_global_binary_path: Optional[str] = None
        self.synthetic_lesions_global_labels_path: Optional[str] = None

        # ----------------------------- Phase 2: SPECT Simulation ----------------------------- 
        # Stage 2.1: SIMIND Preprocessing  
        self.body_seg_arr: Optional[Any] = None  # np.ndarray (float32 mask)
        self.roi_body_seg_arr: Optional[Any] = None  # np.ndarray (int labels)
        self.mask_roi_body: Optional[Dict[int, Any]] = None  # {label_id: bool mask}
        self.class_seg: Optional[Dict[str, int]] = None  # {roi_name: label_id}
        self.atn_av_path: Optional[str] = None # path to SIMIND-formatted atn map (float32 .bin)
        self.binary_roi_act_map_paths: Optional[Dict[str, str]] = None  # {roi_name: path}
        self.arr_px_spacing_cm: Optional[Any] = None  # tuple[float,float,float] (z,y,x)
        self.arr_shape_new: Optional[Any] = None  # tuple[int,int,int] (z,y,x)

        # Stage 2.2: SIMIND  
        self.spect_sim_output_dir: Optional[str] = None  
        self.simind_projection_paths: Optional[Dict[str, Any]] = None # {roi_name: path to SIMIND projection output} - legacy / optional

        # ----------------------------- Phase 3: SPECT Post-Process -----------------------------  
        # Stage 3.1: PBPK  
        self.pbpk_tacs_by_organ: Optional[Dict[str, Any]] = None  
        self.pbpk_frame_start_times_min: Optional[Any] = None  
        self.pbpk_frame_durations_s: Optional[Any] = None  
        self.pbpk_projection_paths: Optional[Dict[str, Any]] = None  
        self.activity_map_sum: Optional[Any] = None  # np.ndarray (n_frames,) - legacy / optional  
        self.activity_organ_sum: Optional[Dict[str, Any]] = None  # {roi: np.ndarray} - legacy / optional  
        self.activity_map_paths_by_organ: Optional[Any] = None  # list[str] - legacy / optional  
        self.pbpk_height_m: Optional[float] = None
        self.pbpk_weight_kg: Optional[float] = None
        self.pbpk_parameters: Optional[Dict[str, Any]] = None

        # Stage 3.2: Reconstruction  
        self.reconstruction_output_dir: Optional[str] = None  

    def require(self, *names: str) -> None:
        """
        Assert that required Context fields exist and are non-None.

        Parameters
        ----------
        *names : str
            Attribute names that must exist on the Context and be non-None.

        Raises
        ------
        AttributeError
            If one or more required fields are missing or None.
        """
        missing = [n for n in names if not hasattr(self, n) or getattr(self, n) is None]
        if missing:
            raise AttributeError(f"Context missing required fields: {missing}")