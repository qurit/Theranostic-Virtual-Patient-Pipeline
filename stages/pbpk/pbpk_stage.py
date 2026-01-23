import os
import numpy as np
import pycno
import nibabel as nib
import pydicom


class PbpkStage:
    def __init__(self, config, context):
        self.config = config
        self.context = context

        subdir_name = config["subdir_names"]["pbpk"]
        output_root = config["output_folder"]["title"]
        self.output_dir = os.path.join(output_root, subdir_name)
        os.makedirs(self.output_dir, exist_ok=True)
        
        self.ct_input_path = config["ct_input"]["path1"]

        self.prefix = config["pbpk"]["name"]

        self.vois_pbpk = config["pbpk"]["VOIs"]
        self.frame_start = config["pbpk"]["FrameStartTimes"]
        self.frame_stop = max(self.frame_start)
        
        self.roi_body_seg_arr = self.context.roi_body_seg_arr
        self.mask_roi_body = self.context.mask_roi_body

    # -----------------------------
    # helpers
    # -----------------------------
    def _remove_background(self, class_seg):
        if "background" in class_seg:  
            class_seg = dict(class_seg)
            del class_seg["background"]  
        return class_seg

    def _voxel_volume_ml(self, arr_px_spacing_cm):
        arr_px_spacing_cm = np.asarray(arr_px_spacing_cm, dtype=float)
        return float(np.prod(arr_px_spacing_cm))  # cm^3 == mL
    
    def _sample_lognormal_from_mean_sd(self, mean, sd):
        """
        Sample LogNormal such that the *resulting* distribution has the requested mean and sd.

        If X ~ LogNormal(mu, sigma^2), then:
          E[X] = exp(mu + sigma^2/2)
          Var[X] = (exp(sigma^2)-1) * exp(2mu + sigma^2)

        Solve:
          sigma^2 = ln(1 + (sd^2 / mean^2))
          mu      = ln(mean) - sigma^2/2
        """
        mean = float(mean)
        sd = float(sd)
        if mean <= 0 or sd <= 0:
            raise ValueError(f"mean and sd must be > 0 for lognormal (got mean={mean}, sd={sd})")

        sigma2 = np.log(1.0 + (sd * sd) / (mean * mean))
        sigma = np.sqrt(sigma2)
        mu = np.log(mean) - 0.5 * sigma2
        return float(np.random.lognormal(mean=mu, sigma=sigma))

    def _extract_height_weight_from_dicom_dir(self, dicom_dir):
        """
        DICOM tags:
          PatientSize   (0010,1020) -> meters
          PatientWeight (0010,1030) -> kg
        Returns (height_m, weight_kg) where either may be None.
        """
        if pydicom is None:
            return None, None
        if not os.path.isdir(dicom_dir):
            return None, None

        # try a handful of files; many DICOMs have no extension
        candidates = []
        for name in sorted(os.listdir(dicom_dir)):
            path = os.path.join(dicom_dir, name)
            if os.path.isfile(path):
                candidates.append(path)
            if len(candidates) >= 50:
                break

        for path in candidates:
            try:
                ds = pydicom.dcmread(path, stop_before_pixels=True, force=True)
                height = getattr(ds, "PatientSize", None)   # meters
                weight = getattr(ds, "PatientWeight", None) # kg

                # sanitize
                height = float(height) if height not in (None, "", " ") else None
                weight = float(weight) if weight not in (None, "", " ") else None

                if height is not None and height <= 0:
                    height = None
                if weight is not None and weight <= 0:
                    weight = None

                # If we found at least one, return (even if the other is None)
                if height is not None or weight is not None:
                    return height, weight
            except Exception:
                continue

        return None, None

    def _parameter_check(self):
        """
        Check inputs and build PBPK parameters.

        Template keys:
          bodyHeight (optional), bodyWeight (optional),
          Rden_Kidney, Rden_SG,
          lambdaRel (kidney), lambdaRel_SG
        """

        # ---- basic validation ----
        if not hasattr(self, "ct_input_path"):
            raise AttributeError("PbpkStage is missing self.ct_input_path")
        if not os.path.exists(self.ct_input_path):
            raise ValueError(f"CT input path does not exist: {self.ct_input_path}")

        if not isinstance(self.vois_pbpk, (list, tuple)) or len(self.vois_pbpk) == 0:
            raise ValueError("PBPK VOIs must be a non-empty list/tuple")

        # allow list/tuple/np array for frame_start
        if not isinstance(self.frame_start, (list, tuple, np.ndarray)) or len(self.frame_start) == 0:
            raise ValueError("Frame start times must be a non-empty list/tuple/array")

        frame_start = np.asarray(self.frame_start, dtype=float)
        if not np.all(np.isfinite(frame_start)):
            raise ValueError("Frame start times contain non-finite values")
        if np.any(frame_start < 0):
            raise ValueError("Frame start times must be >= 0")

        # ---- sample physiological parameters ----
        # (mean, sd) — update if you change literature values
        recep_dens_kidney = self._sample_lognormal_from_mean_sd(30.0, 10.0)     # nmol/L 
        recep_dens_sg     = self._sample_lognormal_from_mean_sd(60.0, 20.0)     # nmol/L
        lambda_rel_kidney = self._sample_lognormal_from_mean_sd(2.88e-4, 0.55e-4)
        lambda_rel_sg     = self._sample_lognormal_from_mean_sd(3.9e-4, 0.63e-4)

        # ---- optional height/weight from DICOM ----
        height = None
        weight = None
        if os.path.isdir(self.ct_input_path):
            height, weight = self._extract_height_weight_from_dicom_dir(self.ct_input_path)

        # ---- build parameters dict (only include height/weight if found) ----
        parameters = {
            "Rden_Kidney": recep_dens_kidney,
            "Rden_SG": recep_dens_sg,
            "lambdaRel_Kidney": lambda_rel_kidney,
            "lambdaRel_SG": lambda_rel_sg,
        }
        if height is not None:
            parameters["bodyHeight"] = height
        if weight is not None:
            parameters["bodyWeight"] = weight

        return parameters

    def _run_psma_model(self):
        """Generate random physiological parameters and run PBPK model."""
        parameters = self._parameter_check()

        model = pycno.Model(
            model_name="PSMA",
            hotamount=10,
            coldamount=100,
            parameters=parameters,
        )

        # PyCNO expects numeric stop/steps; keep your behavior but ensure int steps
        stop = int(self.frame_stop)
        steps = int(np.ceil(stop)) if stop > 0 else 1

        time, tacs = model.simulate(
            stop=stop,
            steps=steps,
            observables=self.vois_pbpk,
        )
        return time, tacs

    def _roi_to_voi(self, roi_name):
        roi_to_voi = {
            "kidney": "Kidney",
            "body": "Rest",
            "liver": "Liver",
            "prostate": "Prostate",
            "heart": "Heart",              
            "spleen": "Spleen",           
            "salivary_glands": "SG",       
        }
        return roi_to_voi.get(roi_name, None)  

    def _save_tac_files(self, roi_name, time, tac_voi, tac_interp):
        roi_tag = roi_name.lower()

        tac_time_f = os.path.join(self.output_dir, f"{self.prefix}_{roi_tag}_TAC_time.bin")
        tac_vals_f = os.path.join(self.output_dir, f"{self.prefix}_{roi_tag}_TAC_values.bin")
        samp_time_f = os.path.join(self.output_dir, f"{self.prefix}_{roi_tag}_sample_times.bin")
        samp_vals_f = os.path.join(self.output_dir, f"{self.prefix}_{roi_tag}_sample_values.bin")

        np.asarray(time, dtype=np.float32).tofile(tac_time_f)
        np.asarray(tac_voi, dtype=np.float32).tofile(tac_vals_f)
        np.asarray(self.frame_start, dtype=np.float32).tofile(samp_time_f)
        np.asarray(tac_interp, dtype=np.float32).tofile(samp_vals_f)

        return {
            "tac_time": tac_time_f,
            "tac_values": tac_vals_f,
            "sample_time": samp_time_f,
            "sample_values": samp_vals_f,
        }

    def _generate_time_activity_arr_roi(
        self,
        roi_name,
        label_value,
        mask_roi_body,
        roi_body_seg_arr,
        voxel_vol_ml,
        time,
        tacs,
        saved_tacs,
    ):
        voi_name = self._roi_to_voi(roi_name)

        if voi_name is None:  
            if "Rest" in self.vois_pbpk:  
                voi_name = "Rest"  
            else:  
                raise ValueError(  
                    f"No VOI mapping for ROI '{roi_name}', and no 'Rest' VOI exists in the PBPK model "
                    f"(supported: {sorted(self.vois_pbpk)})." 
                )

        if voi_name not in self.vois_pbpk:
            if "Rest" in self.vois_pbpk:
                voi_name = "Rest"  # use 'Rest' VOI for unmapped ROIs
            else:
                raise ValueError(
                    f"ROI '{roi_name}' maps to VOI '{voi_name}', but that VOI is not in the PBPK model "
                    f"(supported: {sorted(self.vois_pbpk)}), and no 'Rest' VOI exists."
                )

        voi_index = self.vois_pbpk.index(voi_name)

        mask = mask_roi_body[label_value]
        n_vox = int(np.sum(mask))
        if n_vox == 0:
            raise AssertionError(f"Mask corresponding to {voi_name} is empty")

        tac_voi = tacs[0, :, voi_index]
        tac_interp = np.interp(self.frame_start, time, tac_voi)

        # [frame, z, y, x]
        activity_map_organ = np.zeros((len(self.frame_start), *roi_body_seg_arr.shape), dtype=np.float32)
        activity_map_organ[:, mask] = tac_interp[:, None] / (n_vox * voxel_vol_ml)

        # Save FIRST frame organ map for SIMIND usage
        organ_map_path = os.path.join(self.output_dir, f"{self.prefix}_{roi_name}_act_av.bin")
        activity_map_organ[0].astype(np.float32).tofile(organ_map_path)

        # Organ sum per frame [MBq]
        organ_sum = np.sum(activity_map_organ, axis=(1, 2, 3)) * voxel_vol_ml

        if roi_name not in saved_tacs:  
            saved_tacs[roi_name] = self._save_tac_files(roi_name, time, tac_voi, tac_interp)  

        return activity_map_organ, organ_sum, organ_map_path

    # -----------------------------
    # main
    # -----------------------------
    def run(self):
        for k in ("roi_body_seg_arr", "mask_roi_body", "class_seg", "arr_px_spacing_cm"):  
            if getattr(self.context, k, None) is None:  
                raise AttributeError(f"Context missing required field: {k}")  
            
        class_seg = self._remove_background(self.context.class_seg)
        voxel_vol_ml = self._voxel_volume_ml(self.context.arr_px_spacing_cm)

        time, tacs = self._run_psma_model()

        n_frames = len(self.frame_start)
        activity_map = np.zeros((n_frames, *self.roi_body_seg_arr.shape), dtype=np.float32)

        activity_organ_sum = {}
        organ_paths = []
        saved_tacs = {}

        for roi_name, label_value in class_seg.items():
            out = self._generate_time_activity_arr_roi(
                roi_name=roi_name,
                label_value=label_value,
                mask_roi_body=self.mask_roi_body,
                roi_body_seg_arr=self.roi_body_seg_arr,
                voxel_vol_ml=voxel_vol_ml,
                time=time,
                tacs=tacs,
                saved_tacs=saved_tacs,
            )
            if out is None:
                raise AssertionError(f"Failed to compute activity for ROI '{roi_name}'")

            activity_map_organ, organ_sum, organ_map_path = out
            organ_paths.append(organ_map_path)
            activity_organ_sum[roi_name] = organ_sum

            mask = self.mask_roi_body[label_value]
            activity_map[:, mask] = activity_map_organ[:, mask]

        activity_map_sum = np.sum(activity_map, axis=(1, 2, 3)) * voxel_vol_ml

        # Save full maps per frame (optional, but kept)
        for i, frame in enumerate(activity_map):
            t = self.frame_start[i]
            p = os.path.join(self.output_dir, f"{self.prefix}_{t}_act_av.bin")
            frame.astype(np.float32).tofile(p)

        # Update context
        self.context.activity_map_sum = activity_map_sum
        self.context.activity_organ_sum = activity_organ_sum
        self.context.activity_map_paths_by_organ = organ_paths

        return self.context
