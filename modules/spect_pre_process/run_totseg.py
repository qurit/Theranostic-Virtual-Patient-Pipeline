import os
import dicom2nifti
import torch

import SimpleITK as sitk

from torch.cuda.amp import GradScaler as _GradScaler
from totalsegmentator.python_api import totalsegmentator

# Ensure GradScaler is available as expected by some environments
torch.GradScaler = _GradScaler

def get_preprocessing_output_path(config):

    subdir_name = config["subdir_names"]["spect_preprocessing"]
    output_root = config["output_folder"]["title"]
    output_path = os.path.join(output_root, subdir_name)
    os.makedirs(output_path, exist_ok=True)

    preprocess_file_prefix = config["spect_preprocessing"]["name"]
    return output_path, preprocess_file_prefix

def standardize_ct_to_nifti(ct_input_path, output_path, preprocess_file_prefix):

    ct_nii_path = os.path.join(output_path, f"{preprocess_file_prefix}_ct.nii.gz")

    if os.path.isdir(ct_input_path):
        # DICOM folder -> convert to NIfTI
        dicom2nifti.dicom_series_to_nifti(
            ct_input_path,
            ct_nii_path,
            reorient_nifti=True,
        )
    else:
        # NIfTI file (.nii or .nii.gz) -> copy into output directory
        lower_input = ct_input_path.lower()
        if lower_input.endswith((".nii", ".nii.gz")):
            sitk.WriteImage(sitk.ReadImage(ct_input_path), ct_nii_path, True)
        else:
            raise ValueError(
                "Unsupported CT input. Provide a DICOM folder or a NIfTI "
                f"file (.nii/.nii.gz). Got: {ct_input_path}"
            )

    return ct_nii_path


def run_totseg(config,context):
    
    output_path, preprocess_file_prefix = get_preprocessing_output_path(config)
    
    ct_input_path = config["ct_input"]["path1"]
    device = config["spect_preprocessing"]["device"]
    roi_subset = config["spect_preprocessing"]["roi_subset"]
    ml = config["spect_preprocessing"]["ml"]

    # Standardize CT to a single NIfTI in output_path
    ct_nii_path = standardize_ct_to_nifti(ct_input_path, output_path, preprocess_file_prefix)
    
    # ---- Run TotalSegmentator ----
    # Run multilabel segmentation (roi)
    roi_seg_path = os.path.join(output_path, f"{preprocess_file_prefix}_roi_seg.nii.gz")
    body_seg_path = os.path.join(output_path, f"{preprocess_file_prefix}_body_seg.nii.gz")

    totalsegmentator(
        ct_nii_path, 
        roi_seg_path, 
        device=device,
        ml=ml,
        roi_subset=roi_subset
    )

    # Run body segmentation as a separate task
    totalsegmentator(
        ct_nii_path,
        body_seg_path,
        device=device,
        task="body"
    )

    # Update context
    context.roi_seg_path = roi_seg_path
    context.body_seg_path = body_seg_path
    context.ct_nii_path = ct_nii_path
    
    return context



