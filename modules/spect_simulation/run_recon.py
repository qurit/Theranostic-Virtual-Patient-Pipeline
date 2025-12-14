import os

import numpy as np
import SimpleITK as sitk
import torch
import pytomography
from pytomography.algorithms import OSEM
from pytomography.io.SPECT import simind
from pytomography.io.shared import get_header_value
from pytomography.likelihoods import PoissonLogLikelihood
from pytomography.projectors.SPECT import SPECTSystemMatrix
from pytomography.transforms.SPECT import (
    SPECTAttenuationTransform,
    SPECTPSFTransform,
)

def get_spect_sim_output_path(config):
    subdir_name = config["subdir_names"]["spect_simulation"]
    output_root = config["output_folder"]["title"]
    output_path = os.path.join(output_root, subdir_name)
    os.makedirs(output_path, exist_ok=True)

    spect_sim_file_prefix = config["spect_simulation"]["name"] 
    return output_path, spect_sim_file_prefix

def get_sensitvity_from_calibration_file(calibration_file):
    sensitvity_line = 70
    with open(calibration_file, "r") as file:
        lines = file.readlines()
        sensitivity_line = lines[sensitvity_line].strip()
        # Expect something like: "Sensitivity Cps/MBq: 10.6609"
        sensitivity = float(sensitivity_line.split(":")[-1].strip().split()[0])
    return sensitivity

def get_object_and_proj_metadata(photopeak_path, cor_path):
    object_meta, proj_meta = simind.get_metadata(photopeak_path,cor_path)
    return object_meta, proj_meta

def get_metadata_from_header(photopeak_path, lower_path, upper_path):
    with open(photopeak_path, "r") as f:
        headerdata = np.array(f.readlines())
    proj_dim1 = get_header_value(headerdata, "matrix size [1]", int)
    proj_dim2 = get_header_value(headerdata, "matrix size [2]", int)
    num_proj = get_header_value(headerdata, "total number of images", int)

    ww_peak, ww_lower, ww_upper = [simind.get_energy_window_width(path)
        for path in (photopeak_path, lower_path, upper_path)]
    
    return proj_dim1, proj_dim2, num_proj, ww_peak, ww_lower, ww_upper

def get_cor_data(cor_path):
    cor_path_new = np.loadtxt(cor_path).astype(float)
    if cor_path_new.ndim == 2:
        cor_path_new = cor_path_new[:, 0]
        np.savetxt(f"{cor_path}", cor_path_new)
    return cor_path_new
    
def reshape_projection_data(projection, proj_dim1, proj_dim2, num_proj):
    # Reshape to (num_proj, y, x), flip y, then permute to (num_proj, x, y)
    projection_reshaped = np.transpose(
        projection.reshape((num_proj, proj_dim2, proj_dim1))[:, ::-1], (0, 2, 1))
    return projection_reshaped

def convert_counts_to_mbq_per_ml(reconstructed_image, sensitivity, frame_duration,
                                 output_pixel_width, output_slice_width):
    # Convert to MBq/mL: 
    # counts / (counts/s/MBq) / s / cm^2 / cm -> MBq/mL
    recon_img_arr = (
        reconstructed_image.cpu().T
        / sensitivity
        / frame_duration
        / (output_pixel_width**2)
        / output_slice_width
        )
    return recon_img_arr
    
def get_recon_img(recon_img_arr, output_tuple, likelihood, iterations, subsets,
                    sensitivity, frame_duration, output_pixel_width, output_slice_width):
    
    recon_algorithm = OSEM(likelihood)
    reconstructed_image = recon_algorithm(
            n_iters=iterations,
            n_subsets=subsets,
        )
    
    # Convert to MBq/mL:
    # counts / (counts/s/MBq) / s / cm^2 / cm -> MBq/mL
    recon_img_arr = convert_counts_to_mbq_per_ml(
        reconstructed_image, sensitivity,
        frame_duration * 60,  # min to sec
        output_pixel_width, output_slice_width,
    )
    
    recon_img = sitk.GetImageFromArray(recon_img_arr)
    recon_img.SetSpacing(output_tuple)
    return recon_img


def run_recon(config, context):

    class_seg = context.class_seg
    ActivityMapSum = context.ActivityMapSum

    output_path, spect_sim_file_prefix = get_spect_sim_output_path(config)
    
    frame_start = config["pbpk"]["FrameStartTimes"]  # [min]  
    frame_duration = config["pbpk"]["FrameDurations"][0]  # using first value
    iterations = config["spect_simulation"]["Iterations"]
    subsets = config["spect_simulation"]["Subsets"]
    output_pixel_width = config["spect_simulation"]["OutputPixelWidth"]
    output_slice_width = config["spect_simulation"]["OutputSliceWidth"]

    roi_list = list(class_seg.keys())
    output_tuple = (output_pixel_width, output_pixel_width, output_slice_width)

    # ----- Calibration sensitivity -----
    calibration_file = os.path.join(output_path, "calib.res")
    sensitivity = get_sensitvity_from_calibration_file(calibration_file) # in cps/MBq


    # ----- Projection metadata -----
    
    # Load one of the projection files to get metadata
    photopeak_path = os.path.join(output_path,
                                  f"{spect_sim_file_prefix}_{roi_list[0]}_0_tot_w2.h00")
    lower_path = os.path.join(output_path,
                              f"{spect_sim_file_prefix}_{roi_list[0]}_0_tot_w1.h00")
    upper_path = os.path.join(output_path,
                              f"{spect_sim_file_prefix}_{roi_list[0]}_0_tot_w3.h00")
    
    
    # load cor file 
    cor_path = os.path.join(output_path, 
                            f"{spect_sim_file_prefix}_{roi_list[0]}_0.cor")
    cor_path_new = get_cor_data(cor_path)
    
    object_meta, proj_meta = get_object_and_proj_metadata(photopeak_path, cor_path)
    proj_dim1, proj_dim2, num_proj, ww_peak, ww_lower, ww_upper = get_metadata_from_header(photopeak_path, lower_path, upper_path)

    # ----- Frame loop -----
    for idx, time in enumerate(frame_start):
        # Load SIMIND projection data (counts)
        lower = np.fromfile(os.path.join(output_path,
                                        f"{spect_sim_file_prefix}_{time}min_tot_w1.a00"),
                                        dtype=np.float32)
        photopeak = np.fromfile(os.path.join(output_path,
                                            f"{spect_sim_file_prefix}_{time}min_tot_w2.a00"),
                                            dtype=np.float32)
        upper = np.fromfile(os.path.join(output_path,
                                        f"{spect_sim_file_prefix}_{time}min_tot_w3.a00"),
                                        dtype=np.float32)

        # Reshape to (num_proj, y, x), flip y, then permute to (num_proj, x, y)
        lower = reshape_projection_data(lower, proj_dim1, proj_dim2, num_proj)
        photopeak = reshape_projection_data(photopeak, proj_dim1, proj_dim2, num_proj)
        upper = reshape_projection_data(upper, proj_dim1, proj_dim2, num_proj)

        # Move data to GPU
        lower = torch.tensor(lower.copy()).to(pytomography.device)
        photopeak = torch.tensor(photopeak.copy()).to(pytomography.device)
        upper = torch.tensor(upper.copy()).to(pytomography.device)

        # Poisson realizations of projections (counts)
        photopeak_realization = torch.poisson(photopeak)
        lower_realization = torch.poisson(lower)
        upper_realization = torch.poisson(upper)

        # TEW scatter estimate
        scatter_estimate_tew = simind.compute_EW_scatter(
            lower_realization,
            upper_realization,
            ww_lower,
            ww_upper,
            ww_peak,
        )

        # Attenuation map for reconstruction
        path_amap = os.path.join( output_path, 
                                f"{spect_sim_file_prefix}_{roi_list[0]}_0.hct")
        
        amap = simind.get_attenuation_map(path_amap)
        att_transform = SPECTAttenuationTransform(amap)

        # PSF (resolution) transform
        psf_meta = simind.get_psfmeta_from_header(photopeak_path)
        psf_transform = SPECTPSFTransform(psf_meta)

        # System matrix
        system_matrix = SPECTSystemMatrix(
            obj2obj_transforms=[att_transform, psf_transform],
            proj2proj_transforms=[],
            object_meta=object_meta,
            proj_meta=proj_meta)

        # Likelihood and reconstruction algorithm
        likelihood = PoissonLogLikelihood(
            system_matrix=system_matrix,
            projections=photopeak_realization,
            additive_term=scatter_estimate_tew)

        recon_output_path = os.path.join( output_path,
                                         f"{spect_sim_file_prefix}_{time}min.nii")
        recon_img = get_recon_img(
            recon_img_arr, output_tuple, likelihood, iterations, subsets,
            sensitivity, frame_duration, output_pixel_width, output_slice_width)
        
        sitk.WriteImage(recon_img, recon_output_path, imageIO="NiftiImageIO")



    # Save attenuation image as NIfTI (same spacing as reconstruction)
    atn_img = sitk.GetImageFromArray(amap.cpu().T)
    atn_img.SetSpacing(output_tuple)
    atn_path = os.path.join(output_path, f"{spect_sim_file_prefix}_atn_img.nii")
    sitk.WriteImage(atn_img, atn_path, imageIO="NiftiImageIO")

    return context
