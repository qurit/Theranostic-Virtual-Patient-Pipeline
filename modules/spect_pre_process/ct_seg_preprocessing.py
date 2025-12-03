"""
CT and segmentation preprocessing utilities for SIMIND input.

This module provides helper functions to:

- Convert TotalSegmentator segmentation volumes into class maps and
  per-label binary masks.
- Convert CT images from Hounsfield units (HU) to linear attenuation
  coefficients (mu) and save SIMIND-compatible binary attenuation maps.
- Combine ROI and body segmentations into a single label volume with
  an associated set of masks.
- Preprocess CT and segmentation NIfTI files produced by TotalSegmentator
  into arrays with consistent orientation, resized for SIMIND input, and
  saved as binary files required by the SIMIND SPECT simulation.

The main entry point is `preprocess_ct_and_seg_for_simind`, which is called
from the TDT pipeline main script.
"""

import logging
import os

import nibabel as nib
import numpy as np
from scipy.ndimage import zoom

logger = logging.getLogger(__name__)


def segmentation_2_class(nii_arr, classes, roi_subset):
    """
    Build a dictionary mapping class names to integer labels from a segmentation.

    Parameters
    ----------
    nii_arr : numpy.ndarray
        3D segmentation array where each voxel stores an integer label.
    classes : dict
        Mapping used to convert label values to class names. The exact structure
        depends on the TotalSegmentator class map (e.g., label -> name).
    roi_subset : list of str
        List of ROI names requested from TotalSegmentator. Currently used only
        for logging/debugging.

    Returns
    -------
    class_seg : dict
        Dictionary mapping class names (e.g. "liver") to integer label values.
        The background is always stored as:
            {"Background": 0}
    """
    logger.debug("Running segmentation_2_class. ROI subset: %s", roi_subset)
    logger.debug("Unique segmentation values: %s", np.unique(nii_arr))

    class_seg = {}
    for n in np.unique(nii_arr):
        if n == 0:
            class_seg["Background"] = 0
        else:
            name = classes.get(n)
            class_seg[name] = int(n)

    return class_seg


def segmentation_multiple_arrays(seg_arr):
    """
    Construct one binary mask per non-zero label in a segmentation volume.

    Parameters
    ----------
    seg_arr : numpy.ndarray
        3D segmentation array (integer or float). Non-zero values indicate ROIs.

    Returns
    -------
    masks : dict
        Dictionary mapping each integer label (excluding background 0)
        to a boolean mask of the same shape as `seg_arr`:
            {label_value: mask_array}
    """
    arr = np.asarray(seg_arr)
    labels = np.unique(arr)

    # Remove background label
    labels = labels[labels != 0]

    masks = {}
    for lab in labels:
        masks[int(lab)] = arr == lab

    return masks


def hu_to_mu(ct_array, pixel_size, mu_water=0.1537, mu_bone=0.2234):
    """
    Convert CT values from HU to linear attenuation coefficients (mu).

    A bilinear transformation is used to map HU values to mu (1/cm),
    then scaled by the pixel size to obtain per-pixel attenuation.

    Parameters
    ----------
    ct_array : numpy.ndarray
        Input array of HU values.
    pixel_size : float
        Pixel size in cm (cm/pixel).
    mu_water : float, optional
        Linear attenuation coefficient for water at ~140 keV (1/cm).
        Default is 0.1537.
    mu_bone : float, optional
        Linear attenuation coefficient for bone at ~140 keV (1/cm).
        Default is 0.2234.

    Returns
    -------
    mu_map : numpy.ndarray
        Array of the same shape as `ct_array`, containing linear attenuation
        coefficients (per pixel) computed from HU.
    """
    mu_water_pixel = mu_water * pixel_size  # pixel/cm
    mu_bone_pixel = mu_bone * pixel_size   # pixel/cm

    logger.debug("Attenuation of water: %f pixel/cm", mu_water_pixel)
    logger.debug("Attenuation of bone: %f pixel/cm", mu_bone_pixel)

    mu_map = np.zeros_like(ct_array, dtype=np.float32)

    soft_tissue_mask = ct_array <= 0
    bone_mask = ct_array > 0

    mu_map[soft_tissue_mask] = mu_water_pixel * (
        1 + ct_array[soft_tissue_mask] / 1000.0
    )
    mu_map[bone_mask] = mu_water_pixel + (ct_array[bone_mask] / 1000.0) * (
        mu_bone_pixel - mu_water_pixel
    )

    return mu_map


def save_simind_mu_from_hu(
    hu_arr,
    segmentated_body_output_arr,
    out_dir,
    pixel_size,
    filename="TOTSEG_atn_av.bin",
):
    """
    Convert a CT volume from HU to mu and save as a SIMIND voxel phantom.

    The resulting binary file stores linear attenuation coefficients (mu) in a
    format compatible with SIMIND's voxel-based phantom input.

    Parameters
    ----------
    hu_arr : numpy.ndarray
        3D array of CT values in HU.
    segmentated_body_output_arr : numpy.ndarray
        3D body segmentation array used as a mask. Voxels outside the body
        are set to zero in the final mu map.
    out_dir : str
        Directory where the binary attenuation file will be written.
    pixel_size : float
        Pixel size in cm (cm/pixel).
    filename : str, optional
        Name of the output binary file. Default is "TOTSEG_atn_av.bin".

    Returns
    -------
    bin_path : str
        Filesystem path to the saved binary attenuation map.
    """
    logger.debug("Converting HU values to linear attenuation coefficients (mu).")

    hu = np.asarray(hu_arr, dtype=np.float32)

    mu_map = hu_to_mu(hu, pixel_size)
    mu_map = mu_map * segmentated_body_output_arr

    os.makedirs(out_dir, exist_ok=True)
    bin_path = os.path.join(out_dir, filename)
    mu_map.tofile(bin_path)

    logger.debug("Saved SIMIND attenuation map binary to: %s", bin_path)
    return bin_path


def seg_ROI_plus_body(seg_roi_arr, seg_body_arr):
    """
    Combine ROI and body segmentations into a single label volume.

    Parameters
    ----------
    seg_roi_arr : numpy.ndarray
        3D segmentation array with ROI labels (non-zero inside ROIs).
    seg_body_arr : numpy.ndarray
        3D segmentation array representing the body region.

    Returns
    -------
    seg_plus : numpy.ndarray
        3D array where background is 0, body-only voxels are labeled as 201,
        and ROI voxels retain their original labels.
    mask_body_plus_roi : dict
        Dictionary mapping each non-zero label in `seg_plus` to a boolean mask
        of that label:
            {label_value: mask_array}
    """
    roi = np.asarray(seg_roi_arr)
    body = np.asarray(seg_body_arr)

    roi_mask = roi != 0
    body_mask = body > 0

    seg_plus = np.zeros_like(roi, dtype=np.float32)

    # Body-only voxels labeled with a dedicated class (201)
    seg_plus[body_mask] = 201

    # Overwrite body label with ROI labels where applicable
    seg_plus[roi_mask] = roi[roi_mask]

    mask_body_plus_roi = segmentation_multiple_arrays(seg_plus)

    return seg_plus, mask_body_plus_roi


def preprocess_ct_and_seg_for_simind(
    output_path,
    classes,
    simind_para,
    totseg_para,
    ml_file,
    body_file,
):
    """
    Preprocess CT and segmentation NIfTI files for SIMIND SPECT simulation.

    This function performs the following steps:
      1. Loads the CT input and segmentation outputs (multilabel and body).
      2. Reorients and converts them to float32 numpy arrays with
         shape (Z, Y, X).
      3. Optionally merges left/right kidneys into a single label.
      4. Resizes the arrays to the resolution required by SIMIND.
      5. Builds ROI and body+ROI masks and a class map.
      6. Computes an attenuation map (mu) from CT HU values, masked by body.
      7. Saves the attenuation map and segmentation arrays as binary files.

    Parameters
    ----------
    output_path : str
        Directory where CT input and segmentation NIfTI outputs from
        TotalSegmentator are located and where new binaries will be saved.
    classes : dict
        Mapping used by `segmentation_2_class` to convert label values to
        class names.
    simind_para : dict
        SIMIND-related configuration parameters. Must contain:
            - "resize" : int
    totseg_para : dict
        TotalSegmentator-related configuration parameters. Must contain:
            - "name" : str
            - "roi_subset" : list[str]
    ml_file : str
        Filesystem path to the multilabel segmentation NIfTI file.
    body_file : str
        Path used to locate the body segmentation NIfTI file. The body
        segmentation is expected at: f"{body_file}/body.nii.gz".

    Returns
    -------
    ct_input_arr : numpy.ndarray
        Preprocessed CT array (Z, Y, X) in HU.
    segmentated_ml_output_arr : numpy.ndarray
        Preprocessed multilabel segmentation array (Z, Y, X).
    segmentated_body_output_arr : numpy.ndarray
        Preprocessed body segmentation array (Z, Y, X).
    seg_plus_body_arr : numpy.ndarray
        Combined ROI+body segmentation array (Z, Y, X).
    class_seg : dict
        Mapping from class names to integer labels.
    mask_roi : dict
        Per-ROI masks derived from the multilabel segmentation.
    mask_roi_plus_body : dict
        Per-label masks derived from the combined ROI+body segmentation.
    atn_av_path : str
        Filesystem path to the saved attenuation map binary file.
    seg_ml_bin_path : str
        Filesystem path to the saved multilabel segmentation binary file.
    seg_body_bin_path : str
        Filesystem path to the saved body segmentation binary file.
    pixel_spacing_cm : float
        Pixel spacing in cm (assuming isotropic in-plane).
    slice_thickness : float
        Slice thickness in cm (through-plane spacing).
    ct_get_zoom : tuple of float
        Original voxel spacing (mm) divided by the resize scale factor.
    """
    totseg_name = totseg_para["name"]
    roi_subset = totseg_para["roi_subset"]

    seg_ml_file = ml_file
    seg_body_file = f"{body_file}/body.nii.gz"

    # --- Load NIfTI images ---
    ct_input = nib.load(f"{output_path}/{totseg_name}_ct_input.nii.gz")
    segmentated_ml_output = nib.load(seg_ml_file)
    segmentated_body_output = nib.load(seg_body_file)

    logger.debug(
        "CT, multilabel segmentation and body segmentation NIfTI files loaded "
        "in preprocess_ct_and_seg_for_simind."
    )

    # --- Convert to float32 arrays and reorient to (Z, Y, X) ---
    ct_input_arr = np.transpose(
        np.array(ct_input.get_fdata(dtype=np.float32)),
        (2, 1, 0),
    )[:, ::-1, :]
    segmentated_ml_output_arr = np.transpose(
        np.array(segmentated_ml_output.get_fdata(dtype=np.float32)),
        (2, 1, 0),
    )[:, ::-1, :]
    segmentated_body_output_arr = np.transpose(
        np.array(segmentated_body_output.get_fdata(dtype=np.float32)),
        (2, 1, 0),
    )[:, ::-1, :]

    logger.debug("CT array shape after reorientation: %s", ct_input_arr.shape)
    logger.debug(
        "Multilabel segmentation array shape after reorientation: %s",
        segmentated_ml_output_arr.shape,
    )

    # --- Kidney merging logic ---
    if "kidney_left" in roi_subset and "kidney_right" in roi_subset:
        # Merge left/right kidneys into a single label (200).
        logger.debug(
            "Merging left and right kidneys into a single label: 'kidney' (200)."
        )
        kidney_mask = (segmentated_ml_output_arr == 2) | (
            segmentated_ml_output_arr == 3
        )
        segmentated_ml_output_arr[kidney_mask] = 200
    elif (
        "kidney_left" in roi_subset and "kidney_right" not in roi_subset
    ) or (
        "kidney_right" in roi_subset and "kidney_left" not in roi_subset
    ):
        msg = (
            "Sorry, incapable of running program with only one kidney "
            "(left or right) currently. Please segmentate both!"
        )
        logger.critical(msg)
        raise SystemExit(msg)

    # --- Resize arrays for SIMIND input ---
    resize = simind_para["resize"]
    scale_factor = resize / ct_input_arr.shape[1]  # assuming square in x/y

    ct_input_arr = zoom(
        ct_input_arr, (scale_factor, scale_factor, scale_factor), order=0
    )
    segmentated_ml_output_arr = zoom(
        segmentated_ml_output_arr,
        (scale_factor, scale_factor, scale_factor),
        order=0,
    )
    segmentated_body_output_arr = zoom(
        segmentated_body_output_arr,
        (scale_factor, scale_factor, scale_factor),
        order=0,
    )

    logger.debug(
        "Resized CT and segmentation arrays to shape: %s",
        segmentated_ml_output_arr.shape,
    )

    # --- Build combined ROI+body segmentation and masks ---
    seg_plus_body_arr, mask_roi_plus_body = seg_ROI_plus_body(
        segmentated_ml_output_arr,
        segmentated_body_output_arr,
    )
    mask_roi = segmentation_multiple_arrays(segmentated_ml_output_arr)
    class_seg = segmentation_2_class(seg_plus_body_arr, classes, roi_subset)

    # --- Attenuation map for SIMIND ---
    ct_get_zoom = tuple(np.array(ct_input.header.get_zooms()) / scale_factor)
    pixel_spacing_cm = ct_get_zoom[0] * 0.1  # mm -> cm
    slice_thickness = ct_get_zoom[2] * 0.1   # mm -> cm

    logger.debug("Pixel size: %f cm", pixel_spacing_cm)
    logger.debug("Slice thickness: %f cm", slice_thickness)

    atn_av_path = save_simind_mu_from_hu(
        ct_input_arr,
        segmentated_body_output_arr,
        output_path,
        pixel_spacing_cm,
    )

    # --- Save segmentation arrays as binary files ---
    seg_ml_bin_path = os.path.join(output_path, f"{totseg_name}_ml_segmentation.bin")
    segmentated_ml_output_arr.astype(np.float32).tofile(seg_ml_bin_path)
    logger.debug("Saved multilabel segmentation binary to: %s", seg_ml_bin_path)

    seg_body_bin_path = os.path.join(output_path, f"{totseg_name}_body_segmentation.bin")
    segmentated_body_output_arr.astype(np.float32).tofile(seg_body_bin_path)
    logger.debug("Saved body segmentation binary to: %s", seg_body_bin_path)

    seg_body_plus_ml_bin_path = os.path.join(
        output_path,
        f"{totseg_name}_body_plus_ml_segmentation.bin",
    )
    seg_plus_body_arr.astype(np.float32).tofile(seg_body_plus_ml_bin_path)
    logger.debug(
        "Saved body+ROI segmentation binary to: %s", seg_body_plus_ml_bin_path
    )

    return (
        ct_input_arr,
        segmentated_ml_output_arr,
        segmentated_body_output_arr,
        seg_plus_body_arr,
        class_seg,
        mask_roi,
        mask_roi_plus_body,
        atn_av_path,
        seg_ml_bin_path,
        seg_body_bin_path,
        pixel_spacing_cm,
        slice_thickness,
        ct_get_zoom,
    )
