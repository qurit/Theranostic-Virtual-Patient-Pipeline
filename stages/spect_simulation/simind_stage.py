import os
import shutil
import subprocess
import numpy as np


class SimindSimulationStage:
    def __init__(self, current_dir_path, config, context):
        self.current_dir_path = current_dir_path
        self.config = config
        self.context = context

        subdir_name = config["subdir_names"]["spect_simulation"]
        output_root = config["output_folder"]["title"]
        self.output_dir = os.path.join(output_root, subdir_name)
        os.makedirs(self.output_dir, exist_ok=True)

        self.work_dir = os.path.join(self.output_dir, "simind_work")
        os.makedirs(self.work_dir, exist_ok=True)

        self.prefix = config["spect_simulation"]["name"]

        self.frame_start = config["pbpk"]["FrameStartTimes"]
        self.frame_durations = config["pbpk"]["FrameDurations"] # seconds, len(frame_start)
        self.collimator = config["spect_simulation"]["Collimator"]
        self.isotope = config["spect_simulation"]["Isotope"]
        self.num_projections = config["spect_simulation"]["NumProjections"]
        self.detector_distance = config["spect_simulation"]["DetectorDistance"]
        self.output_img_size = config["spect_simulation"]["OutputImgSize"]
        self.output_pixel_width = config["spect_simulation"]["OutputPixelWidth"]
        self.output_slice_width = config["spect_simulation"]["OutputSliceWidth"]
        self.num_photons = config["spect_simulation"]["NumPhotons"]
        self.simind_dir = config["spect_simulation"]["SIMINDDirectory"]
        self.energy_window_width = config["spect_simulation"]["EnergyWindowWidth"]
        self.detector_width = config["spect_simulation"]["DetectorWidth"] # cm
        self.detector_length = config["spect_simulation"]["DetectorLength"] # cm (if value is 0 will use length of CT)
        
        num_cores = config["spect_simulation"]["NumCores"]
        max_cores = os.cpu_count() or 1  # guard against None
        if isinstance(num_cores, bool) or not isinstance(num_cores, int) or not (1 <= num_cores <= max_cores):
            self.num_cores = max_cores
        else:
            self.num_cores = num_cores
        
        self.simind_exe = os.path.join(self.simind_dir, "simind")

    def _set_simind_environment(self):
        os.environ["SMC_DIR"] = os.path.join(self.simind_dir, "smc_dir/")
        os.environ["PATH"] = (self.simind_dir + ":" + os.environ.get("PATH", ""))

    def _copy_templates(self):
        scatwin_file = os.path.join(self.current_dir_path, "bin", "scattwin.win")
        shutil.copyfile(scatwin_file, os.path.join(self.work_dir, f"{self.prefix}.win"))

        smc_file = os.path.join(self.current_dir_path, "bin", "smc.smc")
        shutil.copyfile(smc_file, os.path.join(self.work_dir, f"{self.prefix}.smc"))

    def _organ_totals_exist(self, organ_name):
        w1 = os.path.join(self.output_dir, f"{self.prefix}_{organ_name}_tot_w1.a00")
        w2 = os.path.join(self.output_dir, f"{self.prefix}_{organ_name}_tot_w2.a00")
        w3 = os.path.join(self.output_dir, f"{self.prefix}_{organ_name}_tot_w3.a00")
        return os.path.exists(w1) and os.path.exists(w2) and os.path.exists(w3)

    def _frame_totals_exist(self):
        for t in self.frame_start:
            w1 = os.path.join(self.output_dir, f"{self.prefix}_{t}min_tot_w1.a00")
            w2 = os.path.join(self.output_dir, f"{self.prefix}_{t}min_tot_w2.a00")
            w3 = os.path.join(self.output_dir, f"{self.prefix}_{t}min_tot_w3.a00")
            if not (os.path.exists(w1) and os.path.exists(w2) and os.path.exists(w3)):
                return False
        return True

    def _calibration_exists(self):
        return os.path.exists(os.path.join(self.output_dir, "calib.res"))

    def _run_jaszczak_calibration(self):
        if self._calibration_exists():
            return

        jaszak_file = os.path.join(self.current_dir_path, "bin", "jaszak.smc")
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

    def _combine_organs_into_frame_totals(self, roi_list, activity_organ_sum):
        for time_index, time in enumerate(self.frame_start):
            xtot_w1 = 0
            xtot_w2 = 0
            xtot_w3 = 0

            for organ in roi_list:
                w1 = np.fromfile(os.path.join(self.output_dir, f"{self.prefix}_{organ}_tot_w1.a00"), dtype=np.float32)
                w2 = np.fromfile(os.path.join(self.output_dir, f"{self.prefix}_{organ}_tot_w2.a00"), dtype=np.float32)
                w3 = np.fromfile(os.path.join(self.output_dir, f"{self.prefix}_{organ}_tot_w3.a00"), dtype=np.float32)

                xtot_w1 += w1 * activity_organ_sum[organ][time_index] * self.frame_durations[time_index]
                xtot_w2 += w2 * activity_organ_sum[organ][time_index] * self.frame_durations[time_index]
                xtot_w3 += w3 * activity_organ_sum[organ][time_index] * self.frame_durations[time_index]

            np.asarray(xtot_w1, dtype=np.float32).tofile(os.path.join(self.output_dir, f"{self.prefix}_{time}min_tot_w1.a00"))
            np.asarray(xtot_w2, dtype=np.float32).tofile(os.path.join(self.output_dir, f"{self.prefix}_{time}min_tot_w2.a00"))
            np.asarray(xtot_w3, dtype=np.float32).tofile(os.path.join(self.output_dir, f"{self.prefix}_{time}min_tot_w3.a00"))

    def _aggregate_core_totals_for_organ(self, organ_name):
        xtot_w1 = 0
        xtot_w2 = 0
        xtot_w3 = 0

        for j in range(self.num_cores):
            w1 = np.fromfile(os.path.join(self.work_dir, f"{self.prefix}_{organ_name}_{j}_tot_w1.a00"), dtype=np.float32)
            w2 = np.fromfile(os.path.join(self.work_dir, f"{self.prefix}_{organ_name}_{j}_tot_w2.a00"), dtype=np.float32)
            w3 = np.fromfile(os.path.join(self.work_dir, f"{self.prefix}_{organ_name}_{j}_tot_w3.a00"), dtype=np.float32)
            xtot_w1 += w1
            xtot_w2 += w2
            xtot_w3 += w3

        xtot_w1 /= self.num_cores
        xtot_w2 /= self.num_cores
        xtot_w3 /= self.num_cores

        np.asarray(xtot_w1, dtype=np.float32).tofile(os.path.join(self.output_dir, f"{self.prefix}_{organ_name}_tot_w1.a00"))
        np.asarray(xtot_w2, dtype=np.float32).tofile(os.path.join(self.output_dir, f"{self.prefix}_{organ_name}_tot_w2.a00"))
        np.asarray(xtot_w3, dtype=np.float32).tofile(os.path.join(self.output_dir, f"{self.prefix}_{organ_name}_tot_w3.a00"))

    def _run_simind_for_organ_cores(self, organ_name, simind_switches):
        processes = []
        for j in range(self.num_cores):
            cmd = f"{self.simind_exe} {self.prefix} {self.prefix}_{organ_name}_{j} " + simind_switches + f"/rr:{j}"
            if j == 0:
                p = subprocess.Popen(cmd, shell=True, cwd=self.work_dir)
            else:
                p = subprocess.Popen(cmd, shell=True, cwd=self.work_dir, stdout=subprocess.DEVNULL)
            processes.append(p)

        for p in processes:
            p.wait()

    def run(self):
        self.context.require(
            "class_seg",
            "roi_seg_arr",
            "activity_organ_sum",
            "activity_map_sum",
            "arr_px_spacing_cm",
            "activity_map_paths_by_organ",
            "atn_av_path",
        )

        class_seg = self.context.class_seg
        roi_seg_arr = self.context.roi_seg_arr
        activity_organ_sum = self.context.activity_organ_sum
        activity_map_sum = self.context.activity_map_sum
        arr_px_spacing_cm = self.context.arr_px_spacing_cm
        organ_act_paths = self.context.activity_map_paths_by_organ
        atn_av_path = self.context.atn_av_path

        if self.num_cores is None or self.num_cores < 1:
            raise RuntimeError("os.cpu_count() returned an invalid value.")

        if not os.path.exists(atn_av_path):
            raise FileNotFoundError(f"Attenuation map not found: {atn_av_path}")

        roi_list = list(class_seg.keys())

        self._set_simind_environment()
        self._copy_templates()

        input_slice_width = arr_px_spacing_cm[0]
        input_pixel_width = arr_px_spacing_cm[1]
        input_half_length = input_slice_width * roi_seg_arr.shape[0] / 2.0

        output_img_length = input_slice_width * roi_seg_arr.shape[0] / self.output_slice_width

        detector_width_cm = self.detector_width # cm
        if self.detector_length == 0:
            detector_length_cm = roi_seg_arr.shape[0] * input_slice_width # cm
        else:
            detector_length_cm = self.detector_length # cm

        atn_name = os.path.basename(atn_av_path)
        atn_work_path = os.path.join(self.work_dir, atn_name)
        if not os.path.exists(atn_work_path):
            shutil.copyfile(atn_av_path, atn_work_path)

        for index, act_path in enumerate(organ_act_paths):
            organ_name = roi_list[index]

            if self._organ_totals_exist(organ_name):
                continue

            if not os.path.exists(act_path):
                raise FileNotFoundError(f"Activity map not found: {act_path}")

            ratio_activity_organ = (activity_organ_sum[organ_name][0] / activity_map_sum[0])
            scale_factor = (self.num_photons * ratio_activity_organ / activity_map_sum[0] / self.num_cores)

            act_name = os.path.basename(act_path)
            shutil.copyfile(act_path, os.path.join(self.work_dir, act_name))

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
                f"/20:{-1*self.energy_window_width}" # if value of -20 is given, means ±10% around centre
                f"/21:{-1*self.energy_window_width}" # dont know behavior if not given to both index 20 and 21 (gave to both)
                f"/28:{self.output_pixel_width}"
                f"/29:{self.num_projections}"
                f"/31:{input_pixel_width}"
                f"/34:{roi_seg_arr.shape[0]}"
                f"/42:{self.detector_distance}"
                f"/76:{self.output_img_size}"
                f"/77:{output_img_length}"
                f"/78:{roi_seg_arr.shape[1]}"
                f"/79:{roi_seg_arr.shape[2]}"
            )

            self._run_simind_for_organ_cores(organ_name, simind_switches)
            self._aggregate_core_totals_for_organ(organ_name)

        if not self._frame_totals_exist():
            self._combine_organs_into_frame_totals(roi_list, activity_organ_sum)

        self._run_jaszczak_calibration()

        # pipeline-friendly + extras
        self.context.spect_sim_output_dir = self.output_dir
        self.context.extras["simind_output_dir"] = self.output_dir
        self.context.extras["simind_work_dir"] = self.work_dir
        self.context.extras["simind_num_cores"] = self.num_cores

        return self.context
