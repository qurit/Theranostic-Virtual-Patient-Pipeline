"""
Organ segmentation utilities using TotalSegmentator for the TDT pipeline.

This module provides a thin wrapper around the TotalSegmentator Python API to:
- Accept either a DICOM folder or NIfTI file as input.
- Convert DICOM to NIfTI (if needed) and standardize the CT input path.
- Run multilabel organ segmentation and body segmentation.
- Return the paths to the resulting NIfTI files.
"""

import logging
import os

import dicom2nifti
import SimpleITK as sitk
import torch
from torch.cuda.amp import GradScaler as _GradScaler
from totalsegmentator.python_api import totalsegmentator

# Ensure GradScaler is available as expected by some environments
torch.GradScaler = _GradScaler

logger = logging.getLogger(__name__)


def run_totseg(input_path, output_path, totseg_para):
    """
    Run TotalSegmentator with TDT-specific convenience options.

    Parameters
    ----------
    input_path : str
        Path to the input data. Can be either:
        - A DICOM series folder, or
        - A NIfTI file (.nii or .nii.gz).
    output_path : str
        Directory where the standardized CT NIfTI and segmentation results
        will be written.
    totseg_para : dict
        Dictionary of TotalSegmentator parameters from the configuration file.
        Expected keys include:
            - "name" : str
            - "device" : str
            - "fast" : bool
            - "roi_subset" : list[str] or None
            - "ml" : bool
            - "statistics" : bool
            - "radiomics" : bool

    Returns
    -------
    ml_file : str or None
        Path to the multilabel segmentation NIfTI file if `ml` mode is enabled,
        otherwise None.
    body_file : str or None
        Path to the body segmentation NIfTI file if `ml` mode is enabled,
        otherwise None.
    """
    out_dir = output_path
    totseg_name = totseg_para["name"]
    device = totseg_para["device"]
    fast = totseg_para["fast"]
    roi_subset = totseg_para["roi_subset"]
    ml = totseg_para["ml"]
    statistics = totseg_para["statistics"]
    radiomics = totseg_para["radiomics"]

    # Radiomics is not supported together with ml mode in TotalSegmentator.
    if radiomics and ml:
        logger.warning(
            "Radiomics is not supported with multilabel (`ml=True`). "
            "Forcing `ml=False` so radiomics can run."
        )
        ml = False

    # ---- Ensure a standardized CT NIfTI in out_dir (e.g. TOTSEG_ct_input.nii.gz) ----
    ct_nifti = os.path.join(out_dir, f"{totseg_name}_ct_input.nii.gz")

    if os.path.isdir(input_path):
        # DICOM folder -> convert to NIfTI
        logger.debug(
            "Detected DICOM folder. Converting DICOM series to NIfTI: %s -> %s",
            input_path,
            ct_nifti,
        )
        dicom2nifti.dicom_series_to_nifti(
            input_path,
            ct_nifti,
            reorient_nifti=True,
        )
    else:
        # NIfTI file (.nii or .nii.gz) -> copy into output directory
        lower_input = input_path.lower()
        if lower_input.endswith((".nii", ".nii.gz")):
            logger.debug(
                "Detected NIfTI file. Copying into output directory: %s -> %s",
                input_path,
                ct_nifti,
            )
            sitk.WriteImage(sitk.ReadImage(input_path), ct_nifti, True)
        else:
            logger.error(
                "Unsupported input type. Provide a DICOM folder or a NIfTI file (.nii/.nii.gz): %s",
                input_path,
            )
            raise ValueError(
                "Unsupported input. Provide a DICOM folder or a NIfTI "
                "file (.nii/.nii.gz)."
            )

    logger.debug("Using CT NIfTI for TotalSegmentator: %s", ct_nifti)

    # Build keyword arguments for TotalSegmentator
    kwargs = {
        "device": device,
    }
    if fast:
        kwargs["fast"] = True
    if roi_subset:
        kwargs["roi_subset"] = roi_subset
    if ml:
        kwargs["ml"] = True
    if statistics:
        kwargs["statistics"] = True
    if radiomics:
        kwargs["radiomics"] = True

    logger.debug("TotalSegmentator keyword arguments: %s", kwargs)

    # If multilabel mode is disabled, we currently skip running segmentation.
    if not ml:
        logger.debug(
            "Multilabel mode (`ml`) is disabled. Skipping organ segmentation. "
            "Returning (None, None)."
        )
        return None, None

    # Run multilabel segmentation
    ml_file = os.path.join(out_dir, f"{totseg_name}_ml_segmentation.nii.gz")
    body_file = os.path.join(out_dir, f"{totseg_name}_body_segmentation.nii.gz")

    logger.debug("Running TotalSegmentator multilabel segmentation: %s", ml_file)
    totalsegmentator(ct_nifti, ml_file, **kwargs)

    # Run body segmentation as a separate task
    logger.debug("Running TotalSegmentator body segmentation: %s", body_file)
    totalsegmentator(
        ct_nifti,
        body_file,
        device=device,
        task="body",
        fast=fast,
    )

    logger.debug("Segmentation complete. Multilabel: %s, Body: %s", ml_file, body_file)

    return ml_file, body_file



