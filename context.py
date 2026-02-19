# =========================
# CONTEXT (minimal edits)
# =========================
import logging

class Context:
    def __init__(self, logger=None):
        # internal
        object.__setattr__(self, "_logger", logger or logging.getLogger("TDT_CONFIG_LOGGER"))
        object.__setattr__(self, "_log_enabled", True)

        # public
        self.extras = {}

        # -----------------------------Inital setup -----------------------------
        self.mode = None
        self.ct_input_path = None
        self.ct_input_type = None
        self.ct_indx = None
        self.output_folder_path = None
        self.subdir_paths = None

        # Config fields
        self.config = {}  # entire config dict for stages to access as needed

        self.output_cfg = None  
        self.subdir_names = None  
        self.spect_preprocessing_cfg = None 
        self.pbpk_cfg = None  
        self.spect_simulation_cfg = None  

        # -----------------------------Stages-----------------------------
        # Stage 1: TotalSegmentator
        self.ct_nii_path = None
        self.body_ml_path = None
        self.head_glands_cavities_ml_path = None
        self.total_ml_path = None
        self.totseg_plan = None

        # Stage 2: ROI Unification
        self.tdt_roi_seg_path = None

        # Stage 3: Preprocessing
        self.body_seg_arr = None
        self.roi_body_seg_arr = None
        self.mask_roi_body = None
        self.class_seg = None
        self.atn_av_path = None
        self.arr_px_spacing_cm = None
        self.arr_shape_new = None

        # Stage 4: PBPK
        self.activity_map_sum = None
        self.activity_organ_sum = None
        self.activity_map_paths_by_organ = None
        self.pbpk_height_m = None
        self.pbpk_weight_kg = None
        self.pbpk_parameters = None

        # Stage 5: SIMIND
        self.spect_sim_output_dir = None

        # Stage 6: Reconstruction
        # (add fields here if you have them)

    def require(self, *names):
        missing = [n for n in names if not hasattr(self, n) or getattr(self, n) is None]
        if missing:
            raise AttributeError(f"Context missing required fields: {missing}")

    def attach_logger(self, logger):
        object.__setattr__(self, "_logger", logger)

    def _summarize(self, value):
        # keep logs readable (no huge arrays dumped)
        try:
            import numpy as np
            if isinstance(value, np.ndarray):
                return f"ndarray shape={value.shape} dtype={value.dtype}"
        except Exception:
            pass

        try:
            import torch
            if isinstance(value, torch.Tensor):
                return f"torch.Tensor shape={tuple(value.shape)} dtype={value.dtype}"
        except Exception:
            pass

        if isinstance(value, dict):
            keys = list(value.keys())
            return f"dict n_keys={len(keys)} keys={keys[:8]}"
        if isinstance(value, (list, tuple)):
            return f"{type(value).__name__} len={len(value)}"
        if isinstance(value, str) and len(value) > 220:
            return value[:220] + "..."
        return repr(value)

    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)

        # don't spam logs for internal fields
        if name in {"_logger", "_log_enabled"}:
            return
        if not getattr(self, "_log_enabled", False):
            return

        logger = getattr(self, "_logger", None)
        if logger:
            logger.info("CTX set: %s = %s", name, self._summarize(value))

    def extras_set(self, key, value):
        """Use this instead of context.extras[key] = ... so it gets logged."""
        self.extras[key] = value
        logger = getattr(self, "_logger", None)
        if logger:
            logger.info("CTX extras[%s] = %s", key, self._summarize(value))

    def extras_update(self, key, patch: dict):
        """Convenient for incremental updates without losing the dict."""
        cur = self.extras.get(key, {})
        if not isinstance(cur, dict):
            cur = {}
        cur.update(patch)
        self.extras_set(key, cur)