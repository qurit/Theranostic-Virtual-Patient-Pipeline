class Context:
    def __init__(self):
        # A place for optional debug / convenience info (not required by pipeline)
        self.extras = {}

        # everything below is being actively used by the pipeline stages
        # -----------------------------
        # Stage 1: TotalSegmentator
        # -----------------------------
        self.ct_nii_path = None # path to CT nii file -> from patient 
        self.body_ml_path = None  # path to body seg -> from totalseg
        self.head_glands_cavities_ml_path = None # path to head glands/cavities seg -> from totalseg
        self.total_ml_path = None # path to total seg -> from totalseg
        self.totseg_plan = None  # plan dict from segmentation stage
        
        # -----------------------------
        # Stage 2: ROI Unification
        # -----------------------------
        self.tdt_roi_seg_path = None # path to unified TDT ROI seg -> from unify stage

        # -----------------------------
        # Stage 3: Preprocessing
        # -----------------------------
        self.body_seg_arr = None   # preprocessed body seg array (simind grid)
        self.roi_body_seg_arr = None  # preprocessed unified ROI seg array including body(simind grid)

        self.mask_roi_body = None   # dict of ROI name -> binary mask array (simind grid)
        self.class_seg = None      # class segmentation map array (simind grid)

        self.atn_av_path = None      # path to attenuation map bin file -> from preprocessing stage
        
        self.arr_px_spacing_cm = None # assuming CT and roi seg are same
        self.arr_shape_new = None # assuming CT and roi seg are same

        # -----------------------------
        # Stage 4: PBPK
        # -----------------------------
        self.activity_map_sum = None          # per-frame total activity [MBq], shape (n_frames,)
        self.activity_organ_sum = None        # dict organ -> per-frame activity [MBq]
        self.activity_map_paths_by_organ = None  # list of per-organ FIRST-frame activity maps (for SIMIND)

        # -----------------------------
        # Stage 5: SIMIND
        # -----------------------------
        self.spect_sim_output_dir = None      # where combined frame totals + calib.res live

        # -----------------------------
        # Stage 6: Reconstruction
        # -----------------------------


    def require(self, *names):
        missing = []
        for name in names:
            if not hasattr(self, name) or getattr(self, name) is None:
                missing.append(name)
        if missing:
            raise AttributeError(f"Context missing required fields: {missing}")

    def __repr__(self):
        keys = list(self.__dict__.keys())
        return f"Context(keys={keys})"
