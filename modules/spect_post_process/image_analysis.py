"""
Post-processing functions for SPECT image analysis.

This module provides utilities to:
- Resample reconstructed SPECT images from the SIMIND / PyTomography grid
  back to the grid of the attenuation map / masks used as SIMIND input.
"""

import logging
import os
from typing import Tuple

import numpy as np
import SimpleITK as sitk

logger = logging.getLogger(__name__)


def resample_spect_to_atn_grid(
    spect_image_path,
    atn_shape,
    pixel_spacing_cm,
    slice_thickness_cm,
    output_dir,
    output_basename = "spect_post_process",
):
    print(pixel_spacing_cm, slice_thickness_cm,output_dir)
    # Read the reconstructed SPECT image with its own spacing
    recon_img = sitk.ReadImage(spect_image_path)
    recon_arr = sitk.GetArrayFromImage(recon_img) # shape: (Z, Y, X)
    logger.debug(
        "Loaded SPECT image from %s with shape %s and spacing %s", 
        spect_image_path,
        recon_arr.shape,
        recon_img.GetSpacing(), 
    )

    # Build a reference image representing the attenuation/mask grid
    atn_template_arr = np.zeros(atn_shape, dtype=np.float32)
    atn_ref_img = sitk.GetImageFromArray(atn_template_arr)
    
    # Set spacing to match the atn/mask grid
    pixel_spacing_cm = float(pixel_spacing_cm)
    slice_thickness_cm = float(slice_thickness_cm)

    atn_ref_img.SetSpacing(
        (pixel_spacing_cm, pixel_spacing_cm, slice_thickness_cm)
    )

    # Use the same origin/direction as the recon image (best guess alignment)
    atn_ref_img.SetOrigin(recon_img.GetOrigin())
    atn_ref_img.SetDirection(recon_img.GetDirection())

    # Resample recon into the atn/mask grid
    resampled_img = sitk.Resample(
        recon_img,
        atn_ref_img,
        sitk.Transform(),           # identity transform
        sitk.sitkLinear,            # linear interpolation
        0.0,                        # default value outside FOV
        recon_img.GetPixelID(),     # keep same pixel type
    )

    resampled_arr = sitk.GetArrayFromImage(resampled_img)
    logger.debug(
        "Resampled SPECT to atn grid. New shape: %s (target: %s)",
        resampled_arr.shape,
        atn_shape,
    )

    # Prepare output path
    post_process_dir = os.path.join(output_dir, "post_process")
    os.makedirs(post_process_dir, exist_ok=True)

    resampled_path = os.path.join(post_process_dir, f"{output_basename}.nii")
    sitk.WriteImage(resampled_img, resampled_path, imageIO="NiftiImageIO")

    logger.info(
        "Resampled SPECT image saved to %s with shape %s and spacing %s",
        resampled_path,
        resampled_arr.shape,
        resampled_img.GetSpacing(),
    )

    return resampled_path
