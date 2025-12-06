"""
Theranostic Digital Twin (TDT) main entry point.

This script incoperates the various functionalites for the Theranostic Digital Twin
(TDT) pipeline. 

The various capabilties of this pipeline include:
    1. Simulating SPECT images from real CT images via organ segmentation,
         physiologically based pharmacokinetic (PBPK) modelling, Monte Carlo SPECT
         simulation (SIMIND), and image reconstruction (PyTomography).
    2. (Future) Simulating PET images from real CT images.
    3. (Future) Performing lesion insertion and analysis.
    :)

For any questions please contact:
Peter Yazdi  <pyazdi@bccrc.ca>
"""

import logging
import os

import setup.init as setup

from modules.spect_pre_process.ct_seg_preprocessing import preprocess_ct_and_seg_for_simind
from modules.spect_post_process.image_analysis import resample_spect_to_atn_grid
from modules.spect_pre_process.run_totseg import run_totseg
from modules.pbpk.run_pbpk import run_pbpk
from modules.spect_simulation.run_simind import run_simind
from modules.spect_simulation.run_recon import run_recon

logger = logging.getLogger(__name__)

# ----- TDT PIPELINE FUNCTIONS -----
def tdt_ct_to_spect(ct_input_path, patient_num):
    """
    Run the full TDT pipeline: real CT image -> simulated SPECT image.
    
    The pipeline consists of the following stages:
      1. Sets up configuration, logging, and output directories.
      2. Runs organ segmentation using Total Segmentator.
      3. Preprocesses CT and segmentation data for SIMIND.
      4. Generates time–activity curves (TACs) and activity maps via PBPK.
      5. Runs SIMIND Monte Carlo SPECT simulation.
      6. Reconstructs the simulated projection data into a SPECT image.

    Parameters
    ----------
    ct_input_path : str
        Path to the input CT image (single volume) in either NIfTI (.nii.gz)
        or DICOM format.
    patient_num : int
        Integer identifier used for naming the output folder for this patient.

    Returns
    -------
    recon_spect_path : str
        Filesystem path to the reconstructed simulated SPECT image generated
        by the pipeline.
    """
    print("[MAIN] Starting Theranostic Digital Twin (TDT) pipeline...")
    logger.info("Starting TDT pipeline for patient %s", patient_num)
    logger.info("Input CT path: %s", ct_input_path)

    # ----- SETUP -----
    output_folder_title = f"Output_Folder_Testing_{patient_num}"  # User adjustable
    current_dir = os.path.dirname(__file__)

    (
        output_dir_path,
        output_subdir_paths,
        totseg_para,
        pbpk_para,
        simind_para,
        recon_para,
    ) = setup.setup_config(output_folder_title, current_dir)

    # Configure logging to file and console (root logger configuration)
    setup.setup_log(
        output_dir_path,
        output_subdir_paths,
        totseg_para,
        pbpk_para,
        simind_para,
        recon_para,
    )

    logger.info("Configuration and logging initialized.")
    logger.debug("Output directory: %s", output_dir_path)
    logger.debug("Output subdirectories: %s", output_subdir_paths)
    logger.debug("TotalSegmentator parameters: %s", totseg_para)
    logger.debug("PBPK parameters: %s", pbpk_para)
    logger.debug("SIMIND parameters: %s", simind_para)
    logger.debug("Reconstruction parameters: %s", recon_para)

    # ----- STAGE 1: ORGAN SEGMENTATION -----
    stage_name = "Organ segmentation (Total Segmentator)"
    print(f"[MAIN] Stage 1/5: {stage_name}...")
    logger.info("Stage 1/5: %s - started", stage_name)
    
    """
    ml_seg_path, body_seg_path = run_totseg(
        ct_input_path,
        output_subdir_paths["output_ORGAN_SEGMENTATION"],
        totseg_para,
    )
    """
    ml_seg_path = "/home/jhubadmin/Theranostic-Virtual-Patient-Pipeline/Output_Folder_Testing_1/TOTSEG_Outputs/TOTSEG_ml_segmentation.nii.gz"
    body_seg_path = "/home/jhubadmin/Theranostic-Virtual-Patient-Pipeline/Output_Folder_Testing_1/TOTSEG_Outputs/TOTSEG_body_segmentation.nii.gz"

    logger.info("Stage 1/5: %s - completed", stage_name)
    logger.debug("ML segmentation output path: %s", ml_seg_path)
    logger.debug("Body segmentation output path: %s", body_seg_path)

    # ----- STAGE 2: SIMIND PREPROCESSING -----
    stage_name = "SIMIND Preprocessing"
    print(f"[MAIN] Stage 2/5: {stage_name}...")
    logger.info("Stage 2/5: %s - started", stage_name)

    (
        ct_input_array,
        segmented_ml_array,
        segmented_body_array,
        seg_plus_body_array,
        class_seg,
        mask_roi,
        mask_roi_plus_body,
        atn_path,
        seg_ml_bin_path,
        seg_body_bin_path,
        pixel_spacing,
        slice_thickness,
        ct_zoom,
    ) = preprocess_ct_and_seg_for_simind(
        output_subdir_paths["output_ORGAN_SEGMENTATION"],
        simind_para,
        totseg_para,
        ml_seg_path,
        body_seg_path,
    )

    logger.info("Stage 2/5: %s - completed", stage_name)
    logger.debug("Attenuation map path: %s", atn_path)
    logger.debug("ML segmentation binary path: %s", seg_ml_bin_path)
    logger.debug("Body segmentation binary path: %s", seg_body_bin_path)
    logger.debug("Pixel spacing (cm): %s", pixel_spacing)
    logger.debug("Slice thickness (cm): %s", slice_thickness)
    logger.debug("CT zoom factor: %s", ct_zoom)
    logger.debug("CT array shape: %s", getattr(ct_input_array, "shape", None))
    logger.debug(
        "Combined segmentation (seg_plus_body_array) shape: %s",
        getattr(seg_plus_body_array, "shape", None),
    )

    # ----- STAGE 3: PBPK MODELLING -----
    stage_name = "PBPK modelling"
    print(f"[MAIN] Stage 3/5: {stage_name}...")
    logger.info("Stage 3/5: %s - started", stage_name)

    (
        activity_map_sum,
        activity_organ_sum,
        act_path_all_organ,
        act_path_all_map,
    ) = run_pbpk(
        output_subdir_paths["output_PBPK"],
        pbpk_para,
        seg_plus_body_array,
        mask_roi_plus_body,
        class_seg,
        ct_zoom,
    )

    logger.info("Stage 3/5: %s - completed", stage_name)
    logger.debug("Activity map (summed) path: %s", act_path_all_map)
    logger.debug("Per-organ activity paths: %s", act_path_all_organ)

    # ----- STAGE 4: SIMIND SIMULATION -----
    stage_name = "SPECT simulation"
    print(f"[MAIN] Stage 4/5: {stage_name}...")
    logger.info("Stage 4/5: %s - started", stage_name)
    """
    run_simind(
        current_dir,
        class_seg,
        simind_para,
        pbpk_para,
        output_subdir_paths["output_SIMIND"],
        os.cpu_count(),
        seg_plus_body_array,
        activity_organ_sum,
        activity_map_sum,
        pixel_spacing,
        slice_thickness,
        act_path_all_organ,
        atn_path,
    )

    logger.info("Stage 4/5: %s - completed", stage_name)
    logger.debug(
        "SIMIND output directory: %s",
        output_subdir_paths["output_SIMIND"],
    )

    # ----- STAGE 5: RECONSTRUCTION -----
    stage_name = "SPECT reconstruction (PyTomography)"
    print(f"[MAIN] Stage 5/5: {stage_name}...")
    logger.info("Stage 5/5: %s - started", stage_name)

    recon_spect_path = run_recon(
        recon_para,
        pbpk_para,
        simind_para,
        output_subdir_paths["output_SIMIND"],
        output_subdir_paths["output_RECON"],
        class_seg,
        activity_map_sum,
    )
    """
    recon_spect_path = "/home/jhubadmin/Theranostic-Virtual-Patient-Pipeline/Output_Folder_Testing_1/RECON_Outputs/RECON_frame0.nii"
    logger.info("Stage 5/5: %s - completed", stage_name)
    logger.info("Reconstructed SPECT image path: %s", recon_spect_path)
    logger.debug(
        "Reconstruction output directory: %s",
        output_subdir_paths["output_RECON"],
    )
    
    logger.info("Post-processing started: Resampling reconstructed SPECT to attenuation map grid.")
    resampled_recon_path = resample_spect_to_atn_grid(
        recon_spect_path,
        seg_plus_body_array.shape,
        pixel_spacing,
        slice_thickness,
        output_subdir_paths["output_RECON"],
    )
    

    # ----- PIPELINE FINISH -----
    print("[MAIN] TDT pipeline finished.")
    logger.info("TDT pipeline finished successfully for patient %s", patient_num)

    return recon_spect_path

def tdt_dosimetry_simulation():
    """
    Placeholder for future dosimetry simulation pipeline function.
    """
    pass

def tdt_ct_to_pet(ct_input_path, patient_num):
    """
    Placeholder for future CT to PET simulation pipeline function.
    """
    pass

def tdt_lesion_stuff():
    """
    Placeholder for future lesion insertion and analysis pipeline function.
    """
    pass


if __name__ == "__main__":
    ct_input_dir = "/home/jhubadmin/Theranostic-Virtual-Patient-Pipeline/input_dir/CT"
    tdt_ct_to_spect(ct_input_dir, patient_num=1)

