class Context:
    def __init__(self):
        # A place for optional debug / convenience info (not required by pipeline)
        self.extras = {}

        # -----------------------------
        # Stage 1: TotalSegmentator
        # -----------------------------
        self.ct_nii_path = None
        self.roi_seg_path = None
        self.body_seg_dir = None  # directory that contains body.nii.gz

        # -----------------------------
        # Stage 2: Preprocessing
        # -----------------------------
        self.ct_arr = None
        self.roi_seg_arr = None
        self.body_seg_arr = None
        self.roi_body_seg_arr = None

        self.mask_roi_body = None
        self.class_seg = None

        self.atn_av_path = None
        self.roi_seg_bin_path = None
        self.body_seg_bin_path = None
        self.roi_body_seg_bin_path = None

        self.roi_subset = None
        self.arr_px_spacing_cm = None

        # -----------------------------
        # Stage 3: PBPK
        # -----------------------------
        self.activity_map_sum = None          # per-frame total activity [MBq], shape (n_frames,)
        self.activity_organ_sum = None        # dict organ -> per-frame activity [MBq]
        self.activity_map_paths_by_organ = None  # list of per-organ FIRST-frame activity maps (for SIMIND)
        self.activity_map_paths_by_frame = None  # list of per-frame activity maps (optional)

        # -----------------------------
        # Stage 4: SIMIND
        # -----------------------------
        self.spect_sim_output_dir = None      # where combined frame totals + calib.res live

        # -----------------------------
        # Stage 5: Reconstruction
        # -----------------------------
        self.recon_paths = None
        self.recon_atn_img = None
        self.recon_atn_path = None

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
