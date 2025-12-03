"""
PBPK modelling utilities for the TDT pipeline.

This module uses the `pycno` PSMA model to generate time–activity curves (TACs)
for a set of physiologically based pharmacokinetic (PBPK) volumes of interest
(VOIs), and then maps those TACs into voxel-wise activity distributions based
on the segmented CT volume.

The main entry point is `run_pbpk`, which:

- Runs the PSMA PBPK model via `pycno.run_model`.
- Interpolates TACs at the requested frame start times.
- Builds a 4D activity map (time, z, y, x) over the entire volume.
- Saves per-organ and per-frame activity maps as binary files.
- Saves TAC time-series and sampled values as binary files for each VOI.
"""

import logging
import os

import numpy as np
import pycno

logger = logging.getLogger(__name__)


def run_pbpk(out_path, pbpk_para, seg_plus_body_arr, masks, class_seg, ct_get_zoom):
    """
    Run the PSMA PBPK model and generate voxel-wise activity maps.

    Parameters
    ----------
    out_path : str
        Directory where PBPK-related output files (activity maps, TAC binaries)
        will be written.
    pbpk_para : dict
        PBPK configuration parameters. Expected keys include:
            - "name" : str
            - "FrameStartTimes" : list[float]
            (Other keys like HotTotalAmount, ColdTotalAmount, LambdaPhys are
             used by upstream configuration/model setup.)
    seg_plus_body_arr : numpy.ndarray
        3D segmentation array (z, y, x) containing combined ROI+body labels.
        Used only for shape; the masks argument encodes region locations.
    masks : dict
        Dictionary mapping integer label values (from `class_seg`) to boolean
        masks (same shape as `seg_plus_body_arr`) indicating voxels belonging
        to each region.
    class_seg : dict
        Mapping from region names (e.g. "liver", "kidney") to integer label
        values used in `seg_plus_body_arr` and `masks`. The key "Background"
        is removed inside this function.
    ct_get_zoom : tuple of float
        Voxel spacing in mm (x, y, z) after any preprocessing scaling. Used to
        convert voxel counts to volume (mL).

    Returns
    -------
    ActivityMapSum : numpy.ndarray
        1D array of length N_frames with total activity (MBq) per frame summed
        over the whole volume.
    ActivityOrganSum : dict
        Mapping from organ name to a 1D array of length N_frames with total
        activity (MBq) per frame for that organ.
    act_path_all_organ : list of str
        List of filesystem paths to per-organ activity map binaries
        (one file per organ, currently storing only the first frame).
    act_path_all_map : list of str
        List of filesystem paths to per-frame whole-volume activity map
        binaries (one file per frame).
    """
    pbpk_name = pbpk_para["name"]

    act_path_all_map = []
    act_path_all_organ = []
    ActivityOrganSum = {}

    # Remove background class so only real organs/regions are iterated over.
    if "Background" in class_seg:
        del class_seg["Background"]

    vois_possible = [
        "Tumor1",
        "Tumor2",
        "Kidney",
        "Heart",
        "SG",
        "Bone",
        "TumorRest",
        "Spleen",
        "Liver",
        "Prostate",
        "GI",
        "Rest",
        "Skin",
        "Muscle",
        "Brain",
        "RedMarrow",
        "Lungs",
        "Adipose",
    ]

    roi_to_voi = {
        # "NA": "Tumor1",
        # "NA": "Tumor2",
        "kidney": "Kidney",
        "heart": "Heart",
        "body": "Rest",
        # "SG": "SG",
        # "Bone": "Bone",
        # "TumorRest": "TumorRest",
        "spleen": "Spleen",
        "liver": "Liver",
        "prostate": "Prostate",
        # "GI": "GI",
        # "Rest": "Rest",
        # "Skin": "Skin",
        # "Muscle": "Muscle",
        "brain": "Brain",
        # "RedMarrow": "RedMarrow",
        # "lungs": "Lungs",
        # "Adipose": "Adipose",
    }

    n_frames = len(pbpk_para["FrameStartTimes"])
    ActivityMap = np.zeros(
        (n_frames, *seg_plus_body_arr.shape),
        dtype=np.float32,
    )

    # Convert voxel sizes (mm) to volume (mL = cm^3)
    pixel_spacing_ml = np.prod(ct_get_zoom) * 0.1**3  # mm^3 * 0.001 -> mL

    # Run PBPK model
    time, TACs = pycno.run_model(
        model_name="PSMA",
        stop=max(pbpk_para["FrameStartTimes"]),
        observables=vois_possible,
    )

    logger.debug("PBPK TACs generated successfully. TAC shape: %s", TACs.shape)

    frame_start = np.asarray(pbpk_para["FrameStartTimes"], float)

    # Populate ActivityMap for each organ/region
    for key, value in class_seg.items():
        if key in roi_to_voi:
            voi = roi_to_voi[key]
        else:
            voi = "Rest"

        voi_index = vois_possible.index(voi)

        # Number of voxels in this ROI
        mask_len_roi = np.sum(masks[value])

        # TAC for this VOI over the model grid (time steps)
        tac_voi = TACs[0, :, voi_index]

        # Interpolate TAC at the desired frame start times
        tac_voi_interp_time = np.interp(frame_start, time, tac_voi)

        # Create an organ-specific activity map: [frame, z, y, x]
        ActivityMap_Organ = np.zeros(
            (n_frames, *seg_plus_body_arr.shape),
            dtype=np.float32,
        )

        # Voxel-wise activity [MBq/mL]:
        # activity_per_voxel = activity_per_voi / (num_voxels * voxel_volume)
        ActivityMap_Organ[:, masks[value]] = (
            tac_voi_interp_time[:, None] / (mask_len_roi * pixel_spacing_ml)
        )

        # Save organ-specific activity map (currently only the first frame)
        ActivityMap_Organ_path = os.path.join(
            out_path,
            f"{pbpk_name}_{key}_act_av.bin",
        )
        ActivityMap_Organ[0].astype(np.float32).tofile(ActivityMap_Organ_path)
        act_path_all_organ.append(ActivityMap_Organ_path)

        # Summed activity per frame for this organ [MBq]
        ActivityOrganSum[key] = (
            np.sum(ActivityMap_Organ, axis=(1, 2, 3)) * pixel_spacing_ml
        )

        # Fill the global ActivityMap with this organ's activity
        ActivityMap[:, masks[value]] = ActivityMap_Organ[:, masks[value]]

        # --- Save TAC time-series to .bin files (float32), no plotting ---
        tac_time_f = os.path.join(out_path, f"{pbpk_name}_{voi}_TAC_time.bin")
        tac_values_f = os.path.join(out_path, f"{pbpk_name}_{voi}_TAC_values.bin")
        samp_time_f = os.path.join(out_path, f"{pbpk_name}_{voi}_sample_times.bin")
        samp_values_f = os.path.join(
            out_path,
            f"{pbpk_name}_{voi}_sample_values.bin",
        )

        # Full model grid (time vs TAC)
        np.asarray(time, dtype=np.float32).tofile(tac_time_f)
        np.asarray(TACs[0, :, voi_index], dtype=np.float32).tofile(tac_values_f)

        # Sampled at frame times
        np.asarray(frame_start, dtype=np.float32).tofile(samp_time_f)
        np.asarray(tac_voi_interp_time, dtype=np.float32).tofile(samp_values_f)

        logger.debug(
            "[PBPK] Saved TAC for VOI %s: %s, %s, %s, %s",
            voi,
            os.path.basename(tac_time_f),
            os.path.basename(tac_values_f),
            os.path.basename(samp_time_f),
            os.path.basename(samp_values_f),
        )

    # Total activity per frame [MBq] across whole volume
    ActivityMapSum = np.sum(ActivityMap, axis=(1, 2, 3)) * pixel_spacing_ml

    logger.debug("ActivityMap shape: %s", ActivityMap.shape)
    logger.debug("ActivityMap total activity per frame (MBq): %s", ActivityMapSum)

    # Save full-frame activity maps (one file per frame)
    for i, frame in enumerate(ActivityMap):
        act_path_single = os.path.join(
            out_path,
            f'{pbpk_name}_{pbpk_para["FrameStartTimes"][i]}_act_av.bin',
        )
        act_path_all_map.append(act_path_single)
        frame.astype(np.float32).tofile(act_path_single)

    return ActivityMapSum, ActivityOrganSum, act_path_all_organ, act_path_all_map

