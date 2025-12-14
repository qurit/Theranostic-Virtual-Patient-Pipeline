import os
import shutil
import subprocess
import numpy as np

def get_spect_sim_output_path(config):
    subdir_name = config["subdir_names"]["spect_simulation"]
    output_root = config["output_folder"]["title"]
    output_path = os.path.join(output_root, subdir_name)
    os.makedirs(output_path, exist_ok=True)

    spect_sim_file_prefix = config["spect_simulation"]["name"] 
    return output_path, spect_sim_file_prefix


def run_simind(current_dir_path, config, context):
    
    output_path, spect_sim_file_prefix = get_spect_sim_output_path(config)
    
    # find parameters needed from context
    class_seg = context.class_seg
    roi_seg_arr = context.roi_seg_arr
    ActivityOrganSum = context.ActivityOrganSum
    ActivityMapSum = context.ActivityMapSum
    arr_px_spacing_cm = context.arr_px_spacing_cm
    act_path_all_organ = context.act_path_all_organ
    atn_av_path = context.atn_av_path

    # find parameters needed from config
    frame_start = config["pbpk"]["FrameStartTimes"]  # [min]
    frame_durations = config["pbpk"]["FrameDurations"] # same length as frame_start
    collimator = config["spect_simulation"]["Collimator"]
    isotope = config["spect_simulation"]["Isotope"]
    num_projections = config["spect_simulation"]["NumProjections"]
    detector_distance = config["spect_simulation"]["DetectorDistance"]
    output_img_size = config["spect_simulation"]["OutputImgSize"]
    output_pixel_width = config["spect_simulation"]["OutputPixelWidth"]
    output_slice_width = config["spect_simulation"]["OutputSliceWidth"]
    num_photons = config["spect_simulation"]["NumPhotons"]

    # extra SIMIND indices
    num_cores = os.cpu_count()
    roi_list = list(class_seg.keys())

    # _____Set environment for SIMIND_____
    os.environ["SMC_DIR"] = os.path.join(config["spect_simulation"]["SIMINDDirectory"], "smc_dir/")
    os.environ["PATH"] = (config["spect_simulation"]["SIMINDDirectory"]+ ":" + os.environ.get("PATH", ""))
    simind_exe = os.path.join(config["spect_simulation"]["SIMINDDirectory"], "simind")


    # _____Copy template files for scatter and system matrix_____
    scatwin_file = os.path.join(current_dir_path, "bin", "scattwin.win")
    scatwin_file_out = os.path.join(output_path, f"{spect_sim_file_prefix}.win")
    shutil.copyfile(scatwin_file, scatwin_file_out)

    smc_file = os.path.join(current_dir_path, "bin", "smc.smc")
    smc_file_out = os.path.join(output_path, f"{spect_sim_file_prefix}.smc")
    shutil.copyfile(smc_file, smc_file_out)

    # _____Derived parameters_____
    input_slice_width = arr_px_spacing_cm[0] # cm
    input_pixel_width = arr_px_spacing_cm[1] # cm
    input_half_length = input_slice_width * roi_seg_arr.shape[0] / 2.0
    
    output_img_length = input_slice_width * roi_seg_arr.shape[0] / output_slice_width
    
    # --- Detector size based on volume ---
    detector_width_cm = 53.3 # standard detector width cm * Symbia Intevo Scanner [cm]
    detector_length_cm = roi_seg_arr.shape[0] * input_slice_width # [cm]

    # Copy attenuation map into SIMIND working directory
    atn_name = os.path.basename(atn_av_path)
    shutil.copyfile(atn_av_path, os.path.join(output_path, atn_name))
    
    # --- Organ-wise SIMIND simulations (frame 0 activity maps) ---
    for index, act_path in enumerate(act_path_all_organ):
        organ_name = roi_list[index]
        # Ratio of organ activity to total activity in the first frame
        ratio_activity_organ = (ActivityOrganSum[organ_name][0] / ActivityMapSum[0])
        scale_factor = (num_photons * ratio_activity_organ / ActivityMapSum[0] / num_cores)
    
        
        # Copy activity map into SIMIND working directory
        act_name = os.path.basename(act_path)
        shutil.copyfile(act_path, os.path.join(output_path, act_name))

        simind_switches = ( # ensure has leading zero for single digit numbers
            f"/fd:{atn_name}"             # attenuation file
            f"/fs:{act_name}"             # activity file
            f"/in:x22,3x"
            f"/nn:{scale_factor}"         # number of photons scaled per organ/core
            f"/cc:{collimator}"
            f"/fi:{isotope}"
            f"/02:{input_half_length}"
            f"/05:{input_half_length}"
            f"/08:{detector_length_cm:.2f}"       # length of detector
            f"/10:{detector_width_cm:.2f}"      # width of detector 
            f"/14:-7" 
            f"/15:-7"
            f"/28:{output_pixel_width}"
            f"/29:{num_projections}"
            f"/31:{input_pixel_width}"
            f"/34:{roi_seg_arr.shape[0]}"
            f"/42:{-1*detector_distance}" # if 42 is negative, absolute value is patient surface to detector distance
            f"/76:{output_img_size}"
            f"/77:{output_img_length}"
            f"/78:{roi_seg_arr.shape[1]}"
            f"/79:{roi_seg_arr.shape[2]}"
        )

        # Launch SIMIND for this organ in parallel over num_cores
        processes = []
        for j in range(num_cores):
            simind_command = (
                f"{simind_exe} {spect_sim_file_prefix} {spect_sim_file_prefix}_{organ_name}_{j} "
            )
            command = simind_command + simind_switches + f"/rr:{j}"  # /rr:{j} : random seed

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
                    f"{spect_sim_file_prefix}_{organ_name}_{j}_tot_w1.a00",
                ),
                dtype=np.float32,
            )
            w2 = np.fromfile(
                os.path.join(
                    output_path,
                    f"{spect_sim_file_prefix}_{organ_name}_{j}_tot_w2.a00",
                ),
                dtype=np.float32,
            )
            w3 = np.fromfile(
                os.path.join(
                    output_path,
                    f"{spect_sim_file_prefix}_{organ_name}_{j}_tot_w3.a00",
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
            os.path.join(output_path, f"{spect_sim_file_prefix}_{organ_name}_tot_w1.a00")
        )
        np.asarray(xtot_w2, dtype=np.float32).tofile(
            os.path.join(output_path, f"{spect_sim_file_prefix}_{organ_name}_tot_w2.a00")
        )
        np.asarray(xtot_w3, dtype=np.float32).tofile(
            os.path.join(output_path, f"{spect_sim_file_prefix}_{organ_name}_tot_w3.a00")
        )

    # --- Frame-wise combination across organs ---
    for frame_index, frame_start_time in enumerate(frame_start):
        xtot_w1 = 0
        xtot_w2 = 0
        xtot_w3 = 0

        for organ in roi_list:
            w1 = np.fromfile(
                os.path.join(
                    output_path,
                    f"{spect_sim_file_prefix}_{organ}_tot_w1.a00",
                ),
                dtype=np.float32,
            )
            w2 = np.fromfile(
                os.path.join(
                    output_path,
                    f"{spect_sim_file_prefix}_{organ}_tot_w2.a00",
                ),
                dtype=np.float32,
            )
            w3 = np.fromfile(
                os.path.join(
                    output_path,
                    f"{spect_sim_file_prefix}_{organ}_tot_w3.a00",
                ),
                dtype=np.float32,
            )

            xtot_w1 += (
                w1
                * ActivityOrganSum[organ][frame_index]
                * frame_durations[frame_index]
            )
            xtot_w2 += (
                w2
                * ActivityOrganSum[organ][frame_index]
                * frame_durations[frame_index]
            )
            xtot_w3 += (
                w3
                * ActivityOrganSum[organ][frame_index]
                * frame_durations[frame_index]
            )

        np.asarray(xtot_w1, dtype=np.float32).tofile(
            os.path.join(
                output_path,
                f"{spect_sim_file_prefix}_{frame_start_time}min_tot_w1.a00",
            )
        )
        np.asarray(xtot_w2, dtype=np.float32).tofile(
            os.path.join(
                output_path,
                f"{spect_sim_file_prefix}_{frame_start_time}min_tot_w2.a00",
            )
        )
        np.asarray(xtot_w3, dtype=np.float32).tofile(
            os.path.join(
                output_path,
                f"{spect_sim_file_prefix}_{frame_start_time}min_tot_w3.a00",
            )
        )

    # --- Calibration (Jaszczak phantom) ---
    jaszak_file = os.path.join(current_dir_path, "bin", "jaszak.smc")
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


    return context
