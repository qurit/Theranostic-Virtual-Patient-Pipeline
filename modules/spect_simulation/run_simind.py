"""
SIMIND SPECT simulation utilities for the TDT pipeline.

This module wraps the SIMIND executable to:

- Set up the required environment variables and template files (.win, .smc).
- Run organ-wise Monte Carlo simulations in parallel across multiple cores.
- Combine organ-specific results into frame-wise projection files.
- Perform a calibration simulation (Jaszczak phantom) to obtain sensitivity.
"""

import logging
import os
import shutil
import subprocess

import numpy as np

logger = logging.getLogger(__name__)


def run_simind(
    base_path,
    class_seg,
    simind_para,
    pbpk_para,
    output_path,
    num_cores,
    segmented_ml_output_arr,
    activity_organ_sum,
    activity_map_sum,
    pixel_spacing_x,
    slice_thickness,
    act_path_all_organ,
    atn_path,
):
    """
    Run SIMIND SPECT simulation using organ-wise activity maps.

    Parameters
    ----------
    base_path : str
        Base path of the project containing the `bin` directory with
        SIMIND template files (e.g., `scattwin.win`, `smc.smc`, `jaszak.smc`).
    class_seg : dict
        Mapping from region names (keys) to integer labels (values). The keys
        define the organ list ordering used in SIMIND runs.
    simind_para : dict
        SIMIND configuration parameters. Expected keys include:
            - "name" : str
            - "SIMINDDirectory" : str
            - "Collimator" : str
            - "Isotope" : str
            - "NumProjections" : int
            - "NumPhotons" : float
            - "DetectorDistance" : float
            - "OutputImgSize" : int
            - "OutputPixelSize" : float
            - "OutputSliceWidth" : float
    pbpk_para : dict
        PBPK configuration parameters. Expected keys include:
            - "FrameStartTimes" : list[float]
            - "FrameDurations" : list[float]
    output_path : str
        Directory to which SIMIND input and output files are written.
    num_cores : int
        Number of parallel SIMIND processes to run for each organ.
    segmented_ml_output_arr : numpy.ndarray
        3D segmentation array (z, y, x), used only to derive volume size.
    activity_organ_sum : dict
        Mapping from organ name (same keys as `class_seg`) to a 1D numpy array
        of activity per frame [MBq].
    activity_map_sum : numpy.ndarray
        1D array of total activity per frame [MBq] across the whole volume.
    pixel_spacing_x : float
        In-plane pixel spacing in cm (used as PixelWidth in SIMIND).
    slice_thickness : float
        Slice thickness in cm (used as SliceWidth in SIMIND).
    act_path_all_organ : list of str
        List of paths to per-organ activity map binaries (frame 0).
    atn_path : str
        Path to the attenuation map binary file.

    Returns
    -------
    int
        Returns 0 on completion.
    """
    output_name = simind_para["name"]

    frames = pbpk_para["FrameStartTimes"]
    frame_durations = pbpk_para["FrameDurations"]

    roi_list = list(class_seg.keys())

    # Set environment for SIMIND
    os.environ["SMC_DIR"] = os.path.join(
        simind_para["SIMINDDirectory"],
        "smc_dir/",
    )
    os.environ["PATH"] = (
        simind_para["SIMINDDirectory"]
        + ":"
        + os.environ.get("PATH", "")
    )

    collimator = simind_para["Collimator"]
    isotope = simind_para["Isotope"]
    num_projections = simind_para["NumProjections"]
    detector_distance = simind_para["DetectorDistance"]
    output_img_size = simind_para["OutputImgSize"]
    output_pixel_size = simind_para["OutputPixelSize"]
    output_slice_width = simind_para["OutputSliceWidth"]
    num_photons = simind_para["NumPhotons"]

    pixel_width = pixel_spacing_x
    slice_width = slice_thickness

    # SIMIND parameters 14 and 15
    index_14 = -7
    index_15 = -7

    simind_exe = os.path.join(simind_para["SIMINDDirectory"], "simind")

    # Copy template files for scatter and system matrix
    scatwin_file = os.path.join(base_path, "bin", "scattwin.win")
    scatwin_file_out = os.path.join(output_path, f"{output_name}.win")
    shutil.copyfile(scatwin_file, scatwin_file_out)

    smc_file = os.path.join(base_path, "bin", "smc.smc")
    smc_file_out = os.path.join(output_path, f"{output_name}.smc")
    shutil.copyfile(smc_file, smc_file_out)

    logger.debug("Scatter and SMC template files copied to: %s", output_path)

    shape = segmented_ml_output_arr.shape

    half_length = slice_width * shape[0] / 2.0
    output_img_length = slice_width * shape[0] / output_slice_width
    
    # --- Detector size based on volume ---
    width_cm = 53.3 # standard detector width cm * Symbia Intevo Scanner
    length_cm = shape[0] * slice_thickness 
    logger.debug(
        "Computed detector size from volume: %d pixels × %f cm = %f cm",
        shape[0],
        slice_thickness,
        length_cm,
    )

    # --- Organ-wise SIMIND simulations (frame 0 activity maps) ---
    for index, act_path in enumerate(act_path_all_organ):
        organ_name = roi_list[index]

        # Ratio of organ activity to total activity in the first frame
        ratio_activity_organ = (
            activity_organ_sum[organ_name][0] / activity_map_sum[0]
        )

        scale_factor = (
            num_photons * ratio_activity_organ / activity_map_sum[0] / num_cores
        )

        logger.debug(
            "Organ %s: scale_factor=%f, activity_ratio=%f",
            organ_name,
            scale_factor,
            ratio_activity_organ,
        )

        atn_name = os.path.basename(atn_path)
        act_name = os.path.basename(act_path)

        # Copy attenuation and activity maps into SIMIND working directory
        shutil.copyfile(atn_path, os.path.join(output_path, atn_name))
        shutil.copyfile(act_path, os.path.join(output_path, act_name))

        simind_switches = ( # ensure has leading zero for single digit numbers
            f"/fd:{atn_name}"             # attenuation file
            f"/fs:{act_name}"             # activity file
            f"/in:x22,3x"
            f"/nn:{scale_factor}"         # number of photons scaled per organ/core
            f"/cc:{collimator}"
            f"/fi:{isotope}"
            f"/02:{half_length}"
            f"/05:{half_length}"
            f"/08:{length_cm:.2f}"       # length of detector
            f"/10:{width_cm:.2f}"      # width of detector 
            f"/14:{index_14}"
            f"/15:{index_15}"
            f"/28:{output_pixel_size}"
            f"/29:{num_projections}"
            f"/31:{pixel_width}"
            f"/34:{shape[0]}"
            f"/42:{-1*detector_distance}" # if 42 is negative, absolute value is patient surface to detector distance
            f"/76:{output_img_size}"
            f"/77:{output_img_length}"
            f"/78:{shape[1]}"
            f"/79:{shape[2]}"
        )

        # Launch SIMIND for this organ in parallel over num_cores
        processes = []
        for j in range(num_cores):
            simind_command = (
                f"{simind_exe} {output_name} {output_name}_{organ_name}_{j} "
            )
            command = simind_command + simind_switches + f"/rr:{j}"  # random seed

            if j == 0:
                p = subprocess.Popen(
                    command,
                    shell=True,
                    cwd=output_path,
                )
            else:
                p = subprocess.Popen(
                    command,
                    shell=True,
                    cwd=output_path,
                    stdout=subprocess.DEVNULL,
                )
            processes.append(p)

        for p in processes:
            p.wait()

        # Combine core-wise results into a single set of TOT files for this organ
        xtot_w1 = 0
        xtot_w2 = 0
        xtot_w3 = 0

        for j in range(num_cores):
            w1 = np.fromfile(
                os.path.join(
                    output_path,
                    f"{output_name}_{organ_name}_{j}_tot_w1.a00",
                ),
                dtype=np.float32,
            )
            w2 = np.fromfile(
                os.path.join(
                    output_path,
                    f"{output_name}_{organ_name}_{j}_tot_w2.a00",
                ),
                dtype=np.float32,
            )
            w3 = np.fromfile(
                os.path.join(
                    output_path,
                    f"{output_name}_{organ_name}_{j}_tot_w3.a00",
                ),
                dtype=np.float32,
            )

            xtot_w1 += w1
            xtot_w2 += w2
            xtot_w3 += w3

        xtot_w1 /= num_cores
        xtot_w2 /= num_cores
        xtot_w3 /= num_cores

        np.asarray(xtot_w1, dtype=np.float32).tofile(
            os.path.join(output_path, f"{output_name}_{organ_name}_tot_w1.a00")
        )
        np.asarray(xtot_w2, dtype=np.float32).tofile(
            os.path.join(output_path, f"{output_name}_{organ_name}_tot_w2.a00")
        )
        np.asarray(xtot_w3, dtype=np.float32).tofile(
            os.path.join(output_path, f"{output_name}_{organ_name}_tot_w3.a00")
        )

        logger.debug(
            "Organ %s: combined TOT files written (w1, w2, w3).",
            organ_name,
        )

    # --- Frame-wise combination across organs ---
    for frame_index, frame_start_time in enumerate(frames):
        xtot_w1 = 0
        xtot_w2 = 0
        xtot_w3 = 0

        for organ in roi_list:
            w1 = np.fromfile(
                os.path.join(
                    output_path,
                    f"{output_name}_{organ}_tot_w1.a00",
                ),
                dtype=np.float32,
            )
            w2 = np.fromfile(
                os.path.join(
                    output_path,
                    f"{output_name}_{organ}_tot_w2.a00",
                ),
                dtype=np.float32,
            )
            w3 = np.fromfile(
                os.path.join(
                    output_path,
                    f"{output_name}_{organ}_tot_w3.a00",
                ),
                dtype=np.float32,
            )

            xtot_w1 += (
                w1
                * activity_organ_sum[organ][frame_index]
                * frame_durations[frame_index]
            )
            xtot_w2 += (
                w2
                * activity_organ_sum[organ][frame_index]
                * frame_durations[frame_index]
            )
            xtot_w3 += (
                w3
                * activity_organ_sum[organ][frame_index]
                * frame_durations[frame_index]
            )

        np.asarray(xtot_w1, dtype=np.float32).tofile(
            os.path.join(
                output_path,
                f"{output_name}_{frame_start_time}min_tot_w1.a00",
            )
        )
        np.asarray(xtot_w2, dtype=np.float32).tofile(
            os.path.join(
                output_path,
                f"{output_name}_{frame_start_time}min_tot_w2.a00",
            )
        )
        np.asarray(xtot_w3, dtype=np.float32).tofile(
            os.path.join(
                output_path,
                f"{output_name}_{frame_start_time}min_tot_w3.a00",
            )
        )

        logger.debug(
            "Frame %s min: combined TOT files written (w1, w2, w3).",
            frame_start_time,
        )

    # --- Calibration (Jaszczak phantom) ---
    jaszak_file = os.path.join(base_path, "bin", "jaszak.smc")
    jaszak_file_out = os.path.join(output_path, "jaszak.smc")
    shutil.copyfile(jaszak_file, jaszak_file_out)

    calibration_command = (
        f"{simind_exe} jaszak calib"
        f"/fi:{isotope}"
        f"/cc:{collimator}"
        f"/29:1"
        f"/15:5"
        f"/fa:11"
        f"/fa:15"
        f"/fa:14"
    )

    subprocess.run(
        calibration_command,
        shell=True,
        cwd=output_path,
        stdout=subprocess.DEVNULL,
    )

    logger.debug("Calibration (Jaszczak phantom) simulation completed.")

    return 0
