import os
import numpy as np
import pycno

def get_pbpk_output_path(config):
    subdir_name = config["subdir_names"]["pbpk"]
    output_root = config["output_folder"]["title"]
    output_path = os.path.join(output_root, subdir_name)
    os.makedirs(output_path, exist_ok=True)

    pbpk_file_prefix = config["pbpk"]["name"] 
    return output_path, pbpk_file_prefix



def class_seg_background_remove(class_seg):
    # Remove background class so only real organs/regions are iterated over. 
    if "Background" in class_seg: 
        del class_seg["Background"] 
    return class_seg

def get_voxel_volume_ml(arr_px_spacing_cm):
    arr_px_spacing_cm = np.asarray(arr_px_spacing_cm, dtype=float)
    voxel_vol_ml = float(np.prod(arr_px_spacing_cm))  # cm^3 == mL
    return voxel_vol_ml


def run_psma_model(config):
    """
    Run the PyCNO PSMA PBPK model and return:
    - list of VOIs
    - frame start times
    - model time grid
    - TAC array
    """
    vois_pbpk = config["pbpk"]["VOIs"]  # e.g. ["Kidney", "Liver", "Rest", ...]
    frame_start = config["pbpk"]["FrameStartTimes"]  # [min] # [min]
    frame_stop = max(frame_start)

    time, TACs = pycno.run_model(
        model_name="PSMA",
        stop=frame_stop,
        observables=vois_pbpk,
    )
    # TACs shape often (1, n_time, n_vois)
    return vois_pbpk, frame_start, time, TACs


def organ_roi_to_voi_name(roi_name):
    roi_to_voi = {
        "kidney": "Kidney",
        "body": "Rest",
        "liver": "Liver",
        "prostate": "Prostate",
    }
    return roi_to_voi.get(roi_name)


def compute_organ_activity_map(roi_name, label_value, mask_roi_body, roi_body_seg_arr, vois_pbpk,
                               frame_start, time, TACs, voxel_vol_ml, output_path, pbpk_file_prefix):
    
    # Map ROI (segmentation name) -> PBPK VOI
    voi_name = organ_roi_to_voi_name(roi_name)

    if voi_name not in vois_pbpk:
        # If VOI isn't in the model observables, treat as 'Rest' (or raise)
        voi_name = "Rest"

    voi_index = vois_pbpk.index(voi_name)

    n_frames = len(frame_start)

    # Number of voxels in this ROI
    mask = mask_roi_body[label_value]
    mask_len_roi = np.sum(mask)

    # TAC for this VOI over the model grid (time steps)
    tac_voi = TACs[0, :, voi_index]

    # Interpolate TAC at the desired frame start times
    tac_voi_interp_time = np.interp(frame_start, time, tac_voi)

    # Activity map for this organ: [frame, z, y, x]
    ActivityMap_Organ = np.zeros(
        (n_frames, *roi_body_seg_arr.shape),
        dtype=np.float32,
    )

    # Voxel-wise activity [MBq/mL]
    # activity_per_voxel = activity_per_voi / (num_voxels * voxel_volume)
    ActivityMap_Organ[:, mask] = (
        tac_voi_interp_time[:, None] / (mask_len_roi * voxel_vol_ml)
    )

    # Save organ-specific activity map (currently only the first frame)
    ActivityMap_Organ_path = os.path.join(
        output_path,
        f"{pbpk_file_prefix}_{roi_name}_act_av.bin",
    )
    ActivityMap_Organ[0].astype(np.float32).tofile(ActivityMap_Organ_path)

    # Summed activity per frame for this organ [MBq]
    ActivityOrganSum = (
        np.sum(ActivityMap_Organ, axis=(1, 2, 3)) * voxel_vol_ml
    )

    # --- Save TAC time-series and sampled TACs ---
    tac_time_f = os.path.join(
        output_path,
        f"{pbpk_file_prefix}_{voi_name}_TAC_time.bin",
    )
    tac_values_f = os.path.join(
        output_path,
        f"{pbpk_file_prefix}_{voi_name}_TAC_values.bin",
    )
    samp_time_f = os.path.join(
        output_path,
        f"{pbpk_file_prefix}_{voi_name}_sample_times.bin",
    )
    samp_values_f = os.path.join(
        output_path,
        f"{pbpk_file_prefix}_{voi_name}_sample_values.bin",
    )

    np.asarray(time, dtype=np.float32).tofile(tac_time_f)
    np.asarray(tac_voi, dtype=np.float32).tofile(tac_values_f)

    np.asarray(frame_start, dtype=np.float32).tofile(samp_time_f)
    np.asarray(tac_voi_interp_time, dtype=np.float32).tofile(samp_values_f)

    tac_paths = {
        "tac_time": tac_time_f,
        "tac_values": tac_values_f,
        "sample_time": samp_time_f,
        "sample_values": samp_values_f,
    }

    return ActivityMap_Organ, ActivityOrganSum, ActivityMap_Organ_path, tac_paths


def run_pbpk(config, context):
    # Output path and naming
    output_path, pbpk_file_prefix = get_pbpk_output_path(config)

    # Geometry / mask info from previous stage
    roi_body_seg_arr = context.roi_body_seg_arr
    mask_roi_body = context.mask_roi_body
    class_seg = class_seg_background_remove(context.class_seg)
    arr_px_spacing_cm = context.arr_px_spacing_cm

    voxel_vol_ml = get_voxel_volume_ml(arr_px_spacing_cm)

    # PBPK model
    vois_pbpk, frame_start, time, TACs = run_psma_model(config)
    n_frames = len(frame_start)

    # Empty sets for outputs
    ActivityMap = np.zeros(
        (n_frames, *roi_body_seg_arr.shape),
        dtype=np.float32,
    )
    ActivityOrganSum = {}
    act_path_all_organ = []
    act_path_all_map = []

    
    for roi_name, label_value in class_seg.items(): # Loop over each ROI class (e.g. 'kidney', 'liver', 'body', ...)

        ActivityMap_Organ, organ_sum, organ_map_path, _tac_paths = compute_organ_activity_map(
            roi_name=roi_name,
            label_value=label_value,
            mask_roi_body=mask_roi_body,
            roi_body_seg_arr=roi_body_seg_arr,
            vois_pbpk=vois_pbpk,
            frame_start=frame_start,
            time=time,
            TACs=TACs,
            voxel_vol_ml=voxel_vol_ml,
            output_path=output_path,
            pbpk_file_prefix=pbpk_file_prefix,
        )

        # Insert this organ's activity into the global ActivityMap
        mask = mask_roi_body[label_value]
        ActivityMap[:, mask] = ActivityMap_Organ[:, mask]

        ActivityOrganSum[roi_name] = organ_sum
        act_path_all_organ.append(organ_map_path)

    # Total activity per frame [MBq] across whole volume
    ActivityMapSum = np.sum(ActivityMap, axis=(1, 2, 3)) * voxel_vol_ml

    # Save full-frame activity maps (one file per frame)
    for i, frame in enumerate(ActivityMap):
        frame_start_time = config["pbpk"]["FrameStartTimes"][i]
        act_path_single = os.path.join(
            output_path,
            f"{pbpk_file_prefix}_{frame_start_time}_act_av.bin",
        )
        frame.astype(np.float32).tofile(act_path_single)
        act_path_all_map.append(act_path_single)

    # --- Update context ---
    context.ActivityMapSum = ActivityMapSum
    context.ActivityOrganSum = ActivityOrganSum
    context.act_path_all_organ = act_path_all_organ
    context.act_path_all_map = act_path_all_map

    return context
