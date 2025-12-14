import os

import nibabel as nib
import numpy as np
from scipy.ndimage import zoom


def get_preprocessing_output_path(config):
    subdir_name = config["subdir_names"]["spect_preprocessing"]
    output_root = config["output_folder"]["title"]
    output_path = os.path.join(output_root, subdir_name)
    os.makedirs(output_path, exist_ok=True)

    preprocess_file_prefix = config["spect_preprocessing"]["name"] 
    return output_path, preprocess_file_prefix


def segmentation_2_class(nii_arr, classes):
    class_seg = {}

    labels = np.unique(nii_arr.astype(int))
    for n in labels:
        if n == 0:
            class_seg["Background"] = 0
            continue

        name = classes.get(str(n))  # JSON keys are strings
        if name is not None:
            class_seg[name] = n

    return class_seg


def mask_seg(arr):
    labels = np.unique(arr)
    labels = labels[labels != 0]

    return {int(lab): (arr == lab) for lab in labels}


def hu_to_mu(ct_arr, pixel_size_cm, mu_water=0.1537, mu_bone=0.2234):

    mu_water_pixel = mu_water * pixel_size_cm  # pixel/cm
    mu_bone_pixel = mu_bone * pixel_size_cm   # pixel/cm

    mu_map = np.zeros_like(ct_arr, dtype=np.float32)

    soft_tissue_mask = ct_arr <= 0
    bone_mask = ct_arr > 0

    mu_map[soft_tissue_mask] = mu_water_pixel * (
        1 + ct_arr[soft_tissue_mask] / 1000.0
    )
    mu_map[bone_mask] = mu_water_pixel + (ct_arr[bone_mask] / 1000.0) * (
        mu_bone_pixel - mu_water_pixel
    )

    return mu_map


def save_simind_mu_from_hu(
    hu_arr,body_seg_arr,
    out_dir,pixel_size_cm,
    filename="spect_preprocessing_atn_av.bin"):

    hu = np.asarray(hu_arr, dtype=np.float32)

    mu_map = hu_to_mu(hu, pixel_size_cm)
    mu_map = mu_map * body_seg_arr

    bin_path = os.path.join(out_dir, filename)
    mu_map.tofile(bin_path)
    return bin_path


def roi_body_seg(roi_seg_arr, body_seg_arr, body_label=201):
    roi_mask = roi_seg_arr != 0
    body_mask = body_seg_arr > 0

    roi_body_seg = np.zeros_like(roi_seg_arr, dtype=np.float32)

    roi_body_seg[body_mask] = body_label
    roi_body_seg[roi_mask] = roi_seg_arr[roi_mask]

    mask_body_roi = mask_seg(roi_body_seg)
    return roi_body_seg, mask_body_roi


def simind_standard_array(nii_load, transpose_tuple=(2, 1, 0), resize=None):

    load_arr = np.array(nii_load.get_fdata(dtype=np.float32))
    arr_std = np.transpose(load_arr, transpose_tuple)[:, ::-1, :]

    scale_factor = 1.0
    if resize is not None:
        # assuming roughly square in x/y
        scale_factor = resize / arr_std.shape[1]
        arr_std = zoom(
            arr_std,
            (scale_factor, scale_factor, scale_factor),
            order=0,  # safe for label
        )

    return arr_std, scale_factor

def kidney_merge_check(roi_subset, roi_seg_arr):

    left_label = 3
    right_label = 2
    merged_label = 200
    
    has_left = "kidney_left" in roi_subset
    has_right = "kidney_right" in roi_subset
    
    # Error if only one kidney is specified
    if has_left and not has_right:
        raise ValueError(
            "PBPK requires both kidneys. Found 'kidney_left' but missing 'kidney_right' in roi_subset."
        )
    if has_right and not has_left:
        raise ValueError(
            "PBPK requires both kidneys. Found 'kidney_right' but missing 'kidney_left' in roi_subset."
        )

    if has_left and has_right:
        # Build a new roi_subset without left/right, then add unified kidney
        new_roi_subset = []
        for roi in roi_subset:
            if roi not in ("kidney_left", "kidney_right"):
                new_roi_subset.append(roi)
        if "kidney" not in new_roi_subset:
            new_roi_subset.append("kidney")
    
        merged = roi_seg_arr.copy()
        kidney_mask = np.isin(roi_seg_arr, (left_label, right_label))
        merged[kidney_mask] = merged_label
        return merged, new_roi_subset

    # Neither kidney present - return unchanged
    return roi_seg_arr, list(roi_subset)


def preprocess_data_for_simind(config, context):

    output_path, preprocess_file_prefix = get_preprocessing_output_path(config)

    ts_classes = config["total_segmentator_classes"]
    roi_subset = config["spect_preprocessing"]["roi_subset"]
    resize_factor = config["spect_preprocessing"]["resize"]

    ct_nii_path = context.ct_nii_path
    roi_seg_path = context.roi_seg_path
    body_seg_path = os.path.join(context.body_seg_path, "body.nii.gz")

    # --- Load NIfTI images ---
    ct_nii = nib.load(ct_nii_path)
    roi_seg_nii = nib.load(roi_seg_path)
    body_seg_nii = nib.load(body_seg_path)

    # --- Create arrays for SIMIND use ---
    ct_arr, scale_factor = simind_standard_array(ct_nii, resize=resize_factor)
    roi_seg_arr, _ = simind_standard_array(roi_seg_nii, resize=resize_factor)
    roi_seg_arr, roi_subset = kidney_merge_check(roi_subset, roi_seg_arr) # check and merge kidneys if needed
    body_seg_arr, _ = simind_standard_array(body_seg_nii, resize=resize_factor)

    # --- Build combined ROI+body segmentation and masks ---
    roi_body_seg_arr, mask_roi_body = roi_body_seg(roi_seg_arr, body_seg_arr)
    class_seg = segmentation_2_class(roi_body_seg_arr, ts_classes)

    # --- Attenuation map for SIMIND ---
    # Original zooms are in mm; divide by scale_factor (resize), then convert to cm.
    zooms_mm = np.array(ct_nii.header.get_zooms()) / scale_factor
    arr_px_spacing_cm = tuple(float(x) * 0.1 for x in zooms_mm)  # mm -> cm

    atn_av_path = save_simind_mu_from_hu(
        ct_arr,
        body_seg_arr,
        output_path,
        pixel_size_cm=arr_px_spacing_cm[0],
    )

    # --- Save segmentation arrays as binary files ---
    roi_seg_bin_path = os.path.join(
        output_path,
        f"{preprocess_file_prefix}_roi_seg.bin",
    )
    roi_seg_arr.astype(np.float32).tofile(roi_seg_bin_path)

    body_seg_bin_path = os.path.join(
        output_path,
        f"{preprocess_file_prefix}_body_seg.bin",
    )
    body_seg_arr.astype(np.float32).tofile(body_seg_bin_path)

    roi_body_seg_bin_path = os.path.join(
        output_path,
        f"{preprocess_file_prefix}_roi_body_seg.bin",
    )
    roi_body_seg_arr.astype(np.float32).tofile(roi_body_seg_bin_path)

    # --- Update context ---
    context.ct_arr = ct_arr
    context.roi_seg_arr = roi_seg_arr
    context.body_seg_arr = body_seg_arr
    context.roi_body_seg_arr = roi_body_seg_arr

    context.mask_roi_body = mask_roi_body
    context.class_seg = class_seg

    context.atn_av_path = atn_av_path
    context.roi_seg_bin_path = roi_seg_bin_path
    context.body_seg_bin_path = body_seg_bin_path
    context.roi_body_seg_bin_path = roi_body_seg_bin_path
    context.roi_subset = roi_subset

    context.arr_px_spacing_cm = arr_px_spacing_cm

    return context
