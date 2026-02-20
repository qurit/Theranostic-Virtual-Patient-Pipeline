"""
SIMIND Simulation Stage for the TDT pipeline.  

This stage runs SIMIND to produce SPECT projection totals (energy windows) using:
- An attenuation map binary produced by preprocessing (`context.atn_av_path`)
- Per-organ activity maps produced by PBPK (`context.activity_map_paths_by_organ`)
- SIMIND template files copied into a working directory (`.win`, `.smc`)

High-level workflow
-------------------
1) Validate required context fields (segmentation labels, activity summaries, spacing, files).
2) Configure SIMIND environment variables (SMC_DIR, PATH).
3) Copy SIMIND template files into a per-CT work directory.
4) For each ROI/organ:
   - Run SIMIND in parallel over `num_cores` by using `/rr:<core_id>` random seeds.
   - Aggregate per-core totals into a single per-organ totals file per energy window.
5) Combine per-organ totals into per-frame totals using the PBPK activity values and frame durations.
6) Run a Jaszczak-based calibration (if not already present) to produce `calib.res`.

Expected Context interface
--------------------------
Incoming `context` is expected to provide:
- context.config : dict with sections:
    - config["spect_simulation"] (SIMIND settings)
    - config["pbpk"] (FrameStartTimes, FrameDurations)
- context.subdir_paths["spect_simulation"] : str
- context.mode : str ("DEBUG" or "PRODUCTION")
- context.require(...) method (used for required-field checks)
- context.extras : dict (optional metadata storage)

And the following fields from earlier stages:
- context.class_seg : dict[str, int] (roi_name -> label_id)
- context.activity_organ_sum : dict[str, np.ndarray] (roi_name -> activity per frame [MBq])
- context.activity_map_sum : np.ndarray (total activity per frame [MBq])
- context.arr_shape_new : tuple[int, int, int] (z,y,x)
- context.arr_px_spacing_cm : tuple[float, float, float] (z,y,x spacing in cm)
- context.activity_map_paths_by_organ : list[str] (first-frame per-organ activity maps)
- context.atn_av_path : str (attenuation binary)

On success, this stage sets:
- context.spect_sim_output_dir : str
- context.extras["simind_output_dir"], ["simind_work_dir"], ["simind_num_cores"]

Maintainer / contact: pyazdi@bccrc.ca  
"""  

from __future__ import annotations  

import os
import shutil
import subprocess
from typing import Any, Dict, List, Sequence  

import numpy as np


class SimindSimulationStage:
    """
    Run SIMIND simulations (per-organ, parallel cores) and assemble per-frame totals.

    Parameters
    ----------
    context : Context-like
        Pipeline context containing config, preprocessing outputs, and PBPK activity outputs.
    """  

    def __init__(self, context: Any) -> None:  
        self.context = context
        self.config: Dict[str, Any] = context.config  

        # Repository root is assumed to be two levels above this stage file.
        self.repo_root: str = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

        self.output_dir: str = context.subdir_paths["spect_simulation"]
        os.makedirs(self.output_dir, exist_ok=True)

        # Work directory where SIMIND templates and input binaries are placed.
        self.work_dir: str = os.path.join(self.output_dir, "simind_work")
        os.makedirs(self.work_dir, exist_ok=True)

        self.prefix: str = context.config["spect_simulation"]["name"]

        self.mode: str = self.context.mode

        # Timing
        self.frame_start: Sequence[float] = context.config["pbpk"]["FrameStartTimes"] # min
        self.frame_durations: Sequence[float] = context.config["pbpk"]["FrameDurations"]  # seconds

        # SIMIND configuration
        self.collimator: str = context.config["spect_simulation"]["Collimator"]
        self.isotope: str = context.config["spect_simulation"]["Isotope"]
        self.num_projections: int = context.config["spect_simulation"]["NumProjections"]
        self.detector_distance: float = context.config["spect_simulation"]["DetectorDistance"]
        self.output_img_size: int = context.config["spect_simulation"]["OutputImgSize"]
        self.output_pixel_width: float = context.config["spect_simulation"]["OutputPixelWidth"]
        self.output_slice_width: float = context.config["spect_simulation"]["OutputSliceWidth"]
        self.num_photons: float = context.config["spect_simulation"]["NumPhotons"]
        self.simind_dir: str = context.config["spect_simulation"]["SIMINDDirectory"]
        self.energy_window_width: float = context.config["spect_simulation"]["EnergyWindowWidth"]
        self.detector_width: float = context.config["spect_simulation"]["DetectorWidth"]  # cm
        self.detector_length: float = context.config["spect_simulation"]["DetectorLength"]  # cm (0 -> use CT length)

        # CPU configuration (validated)
        num_cores = context.config["spect_simulation"]["NumCores"]
        max_cores = os.cpu_count() or 1  # guard against None
        if isinstance(num_cores, bool) or not isinstance(num_cores, int) or not (1 <= num_cores <= max_cores):
            self.num_cores = max_cores
        else:
            self.num_cores = num_cores

        # SIMIND executable 
        simind_exe = os.path.join(self.simind_dir, "simind")
        if not os.path.exists(simind_exe) and os.path.exists(simind_exe + ".exe"):  
            simind_exe = simind_exe + ".exe"  
        self.simind_exe: str = simind_exe  

    def _set_simind_environment(self) -> None:  
        """
        Configure environment variables used by SIMIND.

        Sets:
        - SMC_DIR: where SIMIND expects its `smc_dir` resources
        - PATH: ensures the SIMIND executable directory is discoverable
        """  
        os.environ["SMC_DIR"] = os.path.join(self.simind_dir, "smc_dir")  
        os.environ["PATH"] = (self.simind_dir + os.pathsep + os.environ.get("PATH", "")) 

    def _copy_templates(self) -> None:  
        """
        Copy SIMIND template files into the work directory.

        Required templates in <repo_root>/bin:
        - scattwin.win  -> renamed to <prefix>.win
        - smc.smc       -> renamed to <prefix>.smc
        """  
        scatwin_file = os.path.join(self.repo_root, "bin", "scattwin.win")
        shutil.copyfile(scatwin_file, os.path.join(self.work_dir, f"{self.prefix}.win"))

        smc_file = os.path.join(self.repo_root, "bin", "smc.smc")
        shutil.copyfile(smc_file, os.path.join(self.work_dir, f"{self.prefix}.smc"))

    def _organ_totals_exist(self, organ_name: str) -> bool:  
        """
        Check if aggregated per-organ totals exist in the output directory.

        Returns True if all three energy window totals exist for this organ.
        """  
        w1 = os.path.join(self.output_dir, f"{self.prefix}_{organ_name}_tot_w1.a00")
        w2 = os.path.join(self.output_dir, f"{self.prefix}_{organ_name}_tot_w2.a00")
        w3 = os.path.join(self.output_dir, f"{self.prefix}_{organ_name}_tot_w3.a00")
        return os.path.exists(w1) and os.path.exists(w2) and os.path.exists(w3)

    def _frame_totals_exist(self) -> bool:
        """
        Check if per-frame totals exist for all frames in `self.frame_start`.

        These are written as:
          <prefix>_<t>min_tot_w{1,2,3}.a00
        """ 
        for t in self.frame_start:
            w1 = os.path.join(self.output_dir, f"{self.prefix}_{t}min_tot_w1.a00")
            w2 = os.path.join(self.output_dir, f"{self.prefix}_{t}min_tot_w2.a00")
            w3 = os.path.join(self.output_dir, f"{self.prefix}_{t}min_tot_w3.a00")
            if not (os.path.exists(w1) and os.path.exists(w2) and os.path.exists(w3)):
                return False
        return True

    def _calibration_exists(self) -> bool:
        """Return True if `calib.res` exists in the output directory."""
        return os.path.exists(os.path.join(self.output_dir, "calib.res"))

    def _run_jaszczak_calibration(self) -> None:  
        """
        Run the Jaszczak calibration in SIMIND to produce `calib.res`.

        Notes
        -----
        - If `calib.res` already exists, this is a no-op.
        - Requires `jaszak.smc` to exist in <repo_root>/bin.
        """  
        if self._calibration_exists():
            return

        jaszak_file = os.path.join(self.repo_root, "bin", "jaszak.smc")
        shutil.copyfile(jaszak_file, os.path.join(self.output_dir, "jaszak.smc"))

        cmd = (  
            f"{self.simind_exe} jaszak calib"  
            f"/fi:{self.isotope}"
            f"/cc:{self.collimator}"
            f"/29:1"
            f"/15:5"
            f"/fa:11"
            f"/fa:15"
            f"/fa:14"
        )
        subprocess.run(cmd, shell=True, cwd=self.output_dir, stdout=subprocess.DEVNULL)

    def _combine_organs_into_frame_totals(
        self,
        roi_list: Sequence[str],  
        activity_organ_sum: Dict[str, np.ndarray], 
    ) -> None: 
        """
        Combine per-organ totals into per-frame totals using PBPK activity and frame durations.

        Parameters
        ----------
        roi_list : Sequence[str]
            List of organ/ROI names in the same order as activity maps were generated.
        activity_organ_sum : dict[str, np.ndarray]
            Per-organ activity per frame [MBq], shape (n_frames,).
        """  
        for time_index, time in enumerate(self.frame_start):
            xtot_w1 = 0
            xtot_w2 = 0
            xtot_w3 = 0

            for organ in roi_list:
                w1 = np.fromfile(os.path.join(self.output_dir, f"{self.prefix}_{organ}_tot_w1.a00"), dtype=np.float32)
                w2 = np.fromfile(os.path.join(self.output_dir, f"{self.prefix}_{organ}_tot_w2.a00"), dtype=np.float32)
                w3 = np.fromfile(os.path.join(self.output_dir, f"{self.prefix}_{organ}_tot_w3.a00"), dtype=np.float32)

                # Weight totals by activity at this frame and acquisition duration
                xtot_w1 += w1 * activity_organ_sum[organ][time_index] * self.frame_durations[time_index]
                xtot_w2 += w2 * activity_organ_sum[organ][time_index] * self.frame_durations[time_index]
                xtot_w3 += w3 * activity_organ_sum[organ][time_index] * self.frame_durations[time_index]

            np.asarray(xtot_w1, dtype=np.float32).tofile(os.path.join(self.output_dir, f"{self.prefix}_{time}min_tot_w1.a00"))
            np.asarray(xtot_w2, dtype=np.float32).tofile(os.path.join(self.output_dir, f"{self.prefix}_{time}min_tot_w2.a00"))
            np.asarray(xtot_w3, dtype=np.float32).tofile(os.path.join(self.output_dir, f"{self.prefix}_{time}min_tot_w3.a00"))

    def _aggregate_core_totals_for_organ(self, organ_name: str) -> None: 
        """
        Aggregate totals across `num_cores` SIMIND runs for a given organ.

        Reads per-core files written in work_dir:
          <prefix>_<organ>_<j>_tot_w{1,2,3}.a00

        Writes averaged totals to output_dir:
          <prefix>_<organ>_tot_w{1,2,3}.a00

        If mode == "PRODUCTION", per-core files are deleted to save space.
        """  
        xtot_w1 = 0
        xtot_w2 = 0
        xtot_w3 = 0

        for j in range(self.num_cores):
            p1 = os.path.join(self.work_dir, f"{self.prefix}_{organ_name}_{j}_tot_w1.a00")
            p2 = os.path.join(self.work_dir, f"{self.prefix}_{organ_name}_{j}_tot_w2.a00")
            p3 = os.path.join(self.work_dir, f"{self.prefix}_{organ_name}_{j}_tot_w3.a00")

            w1 = np.fromfile(p1, dtype=np.float32)
            w2 = np.fromfile(p2, dtype=np.float32)
            w3 = np.fromfile(p3, dtype=np.float32)

            xtot_w1 += w1
            xtot_w2 += w2
            xtot_w3 += w3

            if self.mode == "PRODUCTION":
                for p in (p1, p2, p3):
                    try:
                        os.remove(p)  # .a00
                    except FileNotFoundError:
                        pass

        # Average across cores
        xtot_w1 /= self.num_cores
        xtot_w2 /= self.num_cores
        xtot_w3 /= self.num_cores

        np.asarray(xtot_w1, dtype=np.float32).tofile(os.path.join(self.output_dir, f"{self.prefix}_{organ_name}_tot_w1.a00"))
        np.asarray(xtot_w2, dtype=np.float32).tofile(os.path.join(self.output_dir, f"{self.prefix}_{organ_name}_tot_w2.a00"))
        np.asarray(xtot_w3, dtype=np.float32).tofile(os.path.join(self.output_dir, f"{self.prefix}_{organ_name}_tot_w3.a00"))

    def _run_simind_for_organ_cores(self, organ_name: str, simind_switches: str) -> None:  
        """
        Run SIMIND for a single organ across multiple cores (parallel processes).

        Parameters
        ----------
        organ_name : str
            ROI name used in output filenames.
        simind_switches : str
            SIMIND switch string (e.g., "/fd:.../fs:.../nn:...").
        """  
        processes: List[subprocess.Popen] = []  
        for j in range(self.num_cores):
            cmd = f"{self.simind_exe} {self.prefix} {self.prefix}_{organ_name}_{j} " + simind_switches + f"/rr:{j}"
            if j == 0:
                p = subprocess.Popen(cmd, shell=True, cwd=self.work_dir)
            else:
                p = subprocess.Popen(cmd, shell=True, cwd=self.work_dir, stdout=subprocess.DEVNULL)
            processes.append(p)

        for p in processes:
            p.wait()

    def _organ_headers_exist(self, organ_name: str) -> bool:  
        """
        check that SIMIND header outputs exist for an organ (core 0).

        This stage checks for:
          <work_dir>/<prefix>_<organ>_0_tot_w2.h00
        """ 
        return os.path.exists(os.path.join(self.work_dir, f"{self.prefix}_{organ_name}_0_tot_w2.h00"))

    def run(self) -> Any: 
        """
        Execute SIMIND simulation for all organs and assemble totals.

        Returns
        -------
        context : Context-like
            Updated context object.
        """  
        # `context.require` is a pipeline utility that should raise if fields are missing.
        self.context.require(
            "class_seg",
            "roi_body_seg_arr",
            "activity_organ_sum",
            "activity_map_sum",
            "arr_px_spacing_cm",
            "activity_map_paths_by_organ",
            "atn_av_path",
        )

        class_seg = self.context.class_seg
        arr_shape = self.context.arr_shape_new
        activity_organ_sum = self.context.activity_organ_sum
        activity_map_sum = self.context.activity_map_sum
        arr_px_spacing_cm = self.context.arr_px_spacing_cm
        organ_act_paths = self.context.activity_map_paths_by_organ
        atn_av_path = self.context.atn_av_path

        if self.num_cores is None or self.num_cores < 1:
            raise RuntimeError("os.cpu_count() returned an invalid value.")

        if not os.path.exists(atn_av_path):
            raise FileNotFoundError(f"Attenuation map not found: {atn_av_path}")

        # IMPORTANT: This stage assumes `organ_act_paths` aligns with `class_seg.keys()` order.
        # The PBPK stage populates both in the same iteration order, so this should remain consistent. 
        roi_list = list(class_seg.keys())

        self._set_simind_environment()
        self._copy_templates()

        # Input geometry derived from preprocessing spacing (cm) and shape (z,y,x)
        input_slice_width = arr_px_spacing_cm[0]
        input_pixel_width = arr_px_spacing_cm[1]
        input_half_length = input_slice_width * arr_shape[0] / 2.0

        output_img_length = input_slice_width * arr_shape[0] / self.output_slice_width

        detector_width_cm = self.detector_width
        if self.detector_length == 0:
            detector_length_cm = arr_shape[0] * input_slice_width
        else:
            detector_length_cm = self.detector_length

        # Ensure attenuation map is present in work_dir (SIMIND reads inputs from cwd)
        atn_name = os.path.basename(atn_av_path)
        atn_work_path = os.path.join(self.work_dir, atn_name)
        if not os.path.exists(atn_work_path):
            shutil.copyfile(atn_av_path, atn_work_path)

        for index, act_path in enumerate(organ_act_paths):
            organ_name = roi_list[index]

            if self._organ_totals_exist(organ_name) and self._organ_headers_exist(organ_name):
                continue
            if not os.path.exists(act_path):
                raise FileNotFoundError(f"Activity map not found: {act_path}")

            ratio_activity_organ = (activity_organ_sum[organ_name][0] / activity_map_sum[0])
            scale_factor = (self.num_photons * ratio_activity_organ / activity_map_sum[0] / self.num_cores)

            act_name = os.path.basename(act_path)
            shutil.copyfile(act_path, os.path.join(self.work_dir, act_name))

            # SIMIND switches
            simind_switches = (
                f"/fd:{atn_name}"
                f"/fs:{act_name}"
                f"/in:x22,3x"
                f"/nn:{scale_factor}"
                f"/cc:{self.collimator}"
                f"/fi:{self.isotope}"
                f"/02:{input_half_length}"
                f"/05:{input_half_length}"
                f"/08:{detector_length_cm:.2f}"
                f"/10:{detector_width_cm:.2f}"
                f"/14:-7"
                f"/15:-7"
                f"/20:{-1*self.energy_window_width}"
                f"/21:{-1*self.energy_window_width}"
                f"/28:{self.output_pixel_width}"
                f"/29:{self.num_projections}"
                f"/31:{input_pixel_width}"
                f"/34:{arr_shape[0]}"
                f"/42:{self.detector_distance}"
                f"/76:{self.output_img_size}"
                f"/77:{output_img_length}"
                f"/78:{arr_shape[1]}"
                f"/79:{arr_shape[2]}"
            )

            self._run_simind_for_organ_cores(organ_name, simind_switches)
            self._aggregate_core_totals_for_organ(organ_name)

        if not self._frame_totals_exist():
            self._combine_organs_into_frame_totals(roi_list, activity_organ_sum)

        self._run_jaszczak_calibration()

        self.context.spect_sim_output_dir = self.output_dir

        # Optional metadata for debugging
        self.context.extras["simind_output_dir"] = self.output_dir
        self.context.extras["simind_work_dir"] = self.work_dir
        self.context.extras["simind_num_cores"] = self.num_cores

        return self.context