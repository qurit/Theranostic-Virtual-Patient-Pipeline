"""
PBPK Stage (PSMA model via PyCNO) for the TDT pipeline.  

This stage generates time-activity curves (TACs) for user-specified VOIs and converts them
into SIMIND-ready activity maps on the preprocessed SIMIND grid.

Core responsibilities
---------------------
- Validate PBPK configuration inputs (VOIs, frame start times, CT input path).
- Optionally randomize kidney/salivary-gland parameters (lognormal sampling).
- Optionally extract patient height/weight from a DICOM directory (if CT input is DICOM).
- Run the PyCNO PSMA model to obtain TACs (time, tacs).
- For each ROI present in the unified segmentation, map ROI -> VOI and:
    - Interpolate TAC values at configured frame start times.
    - Build a per-frame activity map (MBq/mL) distributed uniformly within the ROI mask.
    - Save first-frame per-organ activity map for SIMIND usage.
    - Save TAC time series + sampled values to .bin files for provenance/analysis.

Expected Context interface
--------------------------
Incoming `context` is expected to provide:
- context.subdir_paths["pbpk"] : str
- context.ct_input_path : str (file path to NIfTI OR directory path to DICOM series)
- context.config["pbpk"] with keys:
    - "name" : str (prefix for outputs)
    - "VOIs" : list[str] (PyCNO observables, e.g., ["Kidney","Liver","Rest",...])
    - "FrameStartTimes" : list[float] (sampling times for frame-wise activity maps)
    - "Randomization_Kidney_SG_Para" : bool
- From preprocessing stage:
    - context.roi_body_seg_arr : np.ndarray (z,y,x int labels)
    - context.mask_roi_body : dict[int, np.ndarray[bool]] (label_id -> mask)
    - context.class_seg : dict[str, int] (roi_name -> label_id)
    - context.arr_px_spacing_cm : tuple[float,float,float] (z,y,x spacing in cm)

On success, this stage sets:
- context.activity_map_sum : np.ndarray (n_frames,) total activity per frame [MBq]
- context.activity_organ_sum : dict[str, np.ndarray] per-ROI activity per frame [MBq]
- context.activity_map_paths_by_organ : list[str] paths to first-frame per-organ maps
- context.pbpk_height_m : Optional[float]
- context.pbpk_weight_kg : Optional[float]
- context.pbpk_parameters : Optional[dict[str, float]]

Maintainer / contact: pyazdi@bccrc.ca  
""" 

from __future__ import annotations  

import os
from typing import Any, Dict, List, Optional, Sequence, Tuple  

import numpy as np
import pycno
import pydicom


class PbpkStage:
    """
    Run PBPK (PyCNO PSMA) and generate Activity Maps for SIMIND.

    Parameters
    ----------
    context : Context-like
        Pipeline context object that carries config, paths, and preprocessed arrays.

    Notes
    -----
    - Activity maps are constructed as uniform concentration within each ROI mask:
        concentration(t) = TAC_interp(t) / (n_voxels * voxel_volume_ml)
      This yields per-voxel activity concentration with units matching your chosen convention
      (commonly MBq/mL if TAC is MBq and voxel_volume_ml is mL).
    """  

    def __init__(self, context: Any) -> None:  
        self.context = context

        self.output_dir: str = context.subdir_paths["pbpk"]  
        os.makedirs(self.output_dir, exist_ok=True)

        self.ct_input_path: str = context.ct_input_path  

        self.prefix: str = context.config["pbpk"]["name"] 
        self.vois_pbpk: List[str] = list(context.config["pbpk"]["VOIs"])  
        self.frame_start: Sequence[float] = context.config["pbpk"]["FrameStartTimes"]  
        self.randomize_kidney_sg_para: bool = context.config["pbpk"]["Randomization_Kidney_SG_Para"]  
        self.frame_stop: float = float(max(self.frame_start)) if self.frame_start else 0.0  

        # Arrays produced by preprocessing stage
        self.roi_body_seg_arr: np.ndarray = self.context.roi_body_seg_arr  
        self.mask_roi_body: Dict[int, np.ndarray] = self.context.mask_roi_body 

        # Optional metadata extracted from DICOM
        self.height: Optional[float] = None  
        self.weight: Optional[float] = None  

        # Parameters passed to PyCNO model
        self.parameters: Dict[str, float] = {}  

    # -----------------------------
    # helpers
    # -----------------------------
    def _remove_background(self, class_seg: Dict[str, int]) -> Dict[str, int]: 
        """
        Remove a "background" entry from a class map if present.

        Parameters
        ----------
        class_seg : dict[str, int]
            ROI name -> label ID mapping.

        Returns
        -------
        dict[str, int]
            The mapping with "background" removed (if it existed).
        """  
        if "background" in class_seg:
            class_seg = dict(class_seg)
            del class_seg["background"]
        return class_seg

    def _voxel_volume_ml(self, arr_px_spacing_cm: Sequence[float]) -> float:  
        """
        Compute voxel volume in mL from pixel spacing in cm.

        Parameters
        ----------
        arr_px_spacing_cm : Sequence[float]
            Spacing in cm in (z, y, x) order.

        Returns
        -------
        float
            Voxel volume in mL (since cm^3 == mL).
        """  
        arr_px_spacing_cm = np.asarray(arr_px_spacing_cm, dtype=float)
        return float(np.prod(arr_px_spacing_cm))  # cm^3 == mL

    def _sample_lognormal_from_mean_sd(self, mean: float, sd: float) -> float:  
        """
        Sample LogNormal so that the resulting distribution has the requested mean and sd.

        If X ~ LogNormal(mu, sigma^2), then:
          E[X] = exp(mu + sigma^2/2)
          Var[X] = (exp(sigma^2)-1) * exp(2mu + sigma^2)

        Solve:
          sigma^2 = ln(1 + (sd^2 / mean^2))
          mu      = ln(mean) - sigma^2/2

        Parameters
        ----------
        mean : float
            Desired mean of the lognormal distribution (must be > 0).
        sd : float
            Desired standard deviation of the lognormal distribution (must be > 0).

        Returns
        -------
        float
            A single sample from the constructed lognormal distribution.

        Raises
        ------
        ValueError
            If mean <= 0 or sd <= 0.
        """  
        mean = float(mean)
        sd = float(sd)
        if mean <= 0 or sd <= 0:
            raise ValueError(f"mean and sd must be > 0 for lognormal (got mean={mean}, sd={sd})")

        sigma2 = np.log(1.0 + (sd * sd) / (mean * mean))
        sigma = np.sqrt(sigma2)
        mu = np.log(mean) - 0.5 * sigma2
        return float(np.random.lognormal(mean=mu, sigma=sigma))

    def _extract_height_weight_from_dicom_dir(self, dicom_dir: str) -> Tuple[Optional[float], Optional[float]]:  
        """
        Try to extract patient height and weight from a DICOM directory.

        DICOM tags:
          - PatientSize   (0010,1020) -> meters
          - PatientWeight (0010,1030) -> kg

        Parameters
        ----------
        dicom_dir : str
            Path to a directory containing DICOM files.

        Returns
        -------
        tuple[Optional[float], Optional[float]]
            (height_m, weight_kg), where either may be None if missing/unreadable.
        """  
        if not os.path.isdir(dicom_dir):
            return None, None

        # Attempts a handful of files; many DICOMs have no extension
        candidates: List[str] = []
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

    def _parameter_check(self) -> Dict[str, float]:  
        """
        Validate PBPK inputs and build the parameters dict passed into PyCNO.

        Returns
        -------
        dict[str, float]
            PBPK parameter overrides (may be empty).

        Raises
        ------
        AttributeError
            If required members are missing.
        ValueError
            If config values are invalid (e.g., empty VOIs or frame times).
        """  
        # ---- basic validation ----
        if not hasattr(self, "ct_input_path"):
            raise AttributeError("PbpkStage is missing self.ct_input_path")
        if not os.path.exists(self.ct_input_path):
            raise ValueError(f"CT input path does not exist: {self.ct_input_path}")

        if not isinstance(self.vois_pbpk, (list, tuple)) or len(self.vois_pbpk) == 0:
            raise ValueError("PBPK VOIs must be a non-empty list/tuple")

        if not isinstance(self.frame_start, (list, tuple, np.ndarray)) or len(self.frame_start) == 0:
            raise ValueError("Frame start times must be a non-empty list/tuple/array")

        frame_start = np.asarray(self.frame_start, dtype=float)
        if not np.all(np.isfinite(frame_start)):
            raise ValueError("Frame start times contain non-finite values")
        if np.any(frame_start < 0):
            raise ValueError("Frame start times must be >= 0")

        parameters: Dict[str, float] = {}

        # ---- Kidney and SG RDen and LamdaRel ----
        if not isinstance(self.randomize_kidney_sg_para, bool):
            raise ValueError("Randomize Parameter must be only True or False in Config")

        vois_set = {str(v).strip().lower() for v in (self.vois_pbpk or [])}
        has_kidney = ("kidney" in vois_set)
        has_sg = ("sg" in vois_set) or any("salivary" in v for v in vois_set)

        if self.randomize_kidney_sg_para:
            if has_sg:
                recep_dens_sg = self._sample_lognormal_from_mean_sd(60.0, 20.0)     # nmol/L
                lambda_rel_sg = self._sample_lognormal_from_mean_sd(3.9e-4, 0.63e-4)
                parameters["Rden_SG"] = recep_dens_sg
                parameters["lambdaRel_SG"] = lambda_rel_sg

            if has_kidney:
                recep_dens_kidney = self._sample_lognormal_from_mean_sd(30.0, 10.0)  # nmol/L
                lambda_rel_kidney = self._sample_lognormal_from_mean_sd(2.88e-4, 0.55e-4)
                parameters["Rden_Kidney"] = recep_dens_kidney
                parameters["lambdaRel_Kidney"] = lambda_rel_kidney

        # ---- Height and Weight (DICOM only) ----
        self.height = None
        self.weight = None
        if os.path.isdir(self.ct_input_path):
            self.height, self.weight = self._extract_height_weight_from_dicom_dir(self.ct_input_path)

        if self.height is not None:
            parameters["bodyHeight"] = float(self.height)
        if self.weight is not None:
            parameters["bodyWeight"] = float(self.weight)

        return parameters

    def _run_psma_model(self) -> Tuple[np.ndarray, np.ndarray]:  
        """
        Generate parameters and run the PyCNO PSMA PBPK model.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            (time, tacs) where:
            - time is shape (T,)
            - tacs is typically shape (1, T, n_vois) in PyCNO
        """  
        self.parameters = self._parameter_check()

        # save provenance early
        self.context.pbpk_height_m = getattr(self, "height", None)
        self.context.pbpk_weight_kg = getattr(self, "weight", None)
        self.context.pbpk_parameters = dict(self.parameters) if self.parameters is not None else None

        model = pycno.Model(
            model_name="PSMA",
            hotamount=10,
            coldamount=100,
            parameters=self.parameters,
        )

        # PyCNO expects numeric stop/steps; keep your behavior but ensure int steps
        stop = int(self.frame_stop)
        steps = int(np.ceil(stop)) if stop > 0 else 1

        time, tacs = model.simulate(
            stop=stop,
            steps=steps,
            observables=self.vois_pbpk,
        )
        return np.asarray(time, dtype=float), np.asarray(tacs, dtype=float)  

    def _roi_to_voi(self, roi_name: str) -> Optional[str]:  
        """
        Map a TDT ROI name (segmentation space) to a PBPK VOI observable name (PyCNO).

        Parameters
        ----------
        roi_name : str
            ROI name from `context.class_seg` (e.g., "kidney", "salivary_glands").

        Returns
        -------
        Optional[str]
            Mapped VOI name, or None if there is no explicit mapping.
        """  
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

    def _save_tac_files(
        self,
        roi_name: str, 
        time: np.ndarray,  
        tac_voi: np.ndarray,  
        tac_interp: np.ndarray,  
    ) -> Dict[str, str]:  
        """
        Save TAC full-resolution and sampled-at-frame-start arrays to binary files.

        Files written (all float32):
        - <prefix>_<roi>_TAC_time.bin
        - <prefix>_<roi>_TAC_values.bin
        - <prefix>_<roi>_sample_times.bin
        - <prefix>_<roi>_sample_values.bin

        Returns
        -------
        dict[str, str]
            Paths to each saved file.
        """  
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
        roi_name: str,  
        label_value: int,  
        mask_roi_body: Dict[int, np.ndarray],  
        roi_body_seg_arr: np.ndarray,  
        voxel_vol_ml: float,  
        time: np.ndarray,  
        tacs: np.ndarray,  
        saved_tacs: Dict[str, Dict[str, str]], 
    ) -> Tuple[np.ndarray, np.ndarray, str]:  
        """
        Build per-frame activity map for a single ROI and save first frame for SIMIND.

        Parameters
        ----------
        roi_name : str
            ROI name (from context.class_seg).
        label_value : int
            Integer label ID in `roi_body_seg_arr`.
        mask_roi_body : dict[int, np.ndarray]
            Mapping label_id -> boolean mask.
        roi_body_seg_arr : np.ndarray
            Multilabel segmentation array (z,y,x).
        voxel_vol_ml : float
            Voxel volume in mL.
        time : np.ndarray
            Model simulation times (T,).
        tacs : np.ndarray
            Model TAC array (typically shape (1, T, n_vois)).
        saved_tacs : dict
            Cache dict used to avoid re-writing TAC files for the same ROI.

        Returns
        -------
        tuple[np.ndarray, np.ndarray, str]
            (activity_map_organ, organ_sum, organ_map_path) where:
            - activity_map_organ: shape (n_frames, z, y, x), float32
            - organ_sum: shape (n_frames,), float64/float32 (MBq)
            - organ_map_path: path to the first-frame organ map binary
        """  
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
        tac_interp = np.interp(np.asarray(self.frame_start, dtype=float), time, tac_voi)

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
    def run(self) -> Any:  
        """
        Execute PBPK simulation and write activity maps.

        Returns
        -------
        context : Context-like
            Updated context with activity summaries and output paths populated.
        """  
        for k in ("roi_body_seg_arr", "mask_roi_body", "class_seg", "arr_px_spacing_cm"):
            if getattr(self.context, k, None) is None:
                raise AttributeError(f"Context missing required field: {k}")

        class_seg = self._remove_background(self.context.class_seg)
        voxel_vol_ml = self._voxel_volume_ml(self.context.arr_px_spacing_cm)

        time, tacs = self._run_psma_model()

        n_frames = len(self.frame_start)
        activity_map = np.zeros((n_frames, *self.roi_body_seg_arr.shape), dtype=np.float32)

        activity_organ_sum: Dict[str, np.ndarray] = {}  
        organ_paths: List[str] = []  
        saved_tacs: Dict[str, Dict[str, str]] = {}  

        for roi_name, label_value in class_seg.items():
            activity_map_organ, organ_sum, organ_map_path = self._generate_time_activity_arr_roi(
                roi_name=roi_name,
                label_value=int(label_value),
                mask_roi_body=self.mask_roi_body,
                roi_body_seg_arr=self.roi_body_seg_arr,
                voxel_vol_ml=voxel_vol_ml,
                time=time,
                tacs=tacs,
                saved_tacs=saved_tacs,
            )

            organ_paths.append(organ_map_path)
            activity_organ_sum[roi_name] = organ_sum

            mask = self.mask_roi_body[int(label_value)]
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