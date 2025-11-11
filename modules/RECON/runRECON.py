import os
import torch
import numpy as np
from pytomography.io.SPECT import simind
from pytomography.projectors.SPECT import SPECTSystemMatrix
from pytomography.transforms.SPECT import SPECTAttenuationTransform, SPECTPSFTransform
from pytomography.algorithms import OSEM
from pytomography.likelihoods import PoissonLogLikelihood
import SimpleITK as sitk
import logging as log
from pytomography.io.shared import get_header_value
import pytomography
import torch


def runRECON(recon_para,pbpk_para,simind_para,out_paths, class_seg, ActivityMapSum):
    output_name_simind = simind_para['name']
    output_name_recon = recon_para['name']
    output_path_simind = out_paths['output_SIMIND']
    output_path_recon = out_paths['output_PyTomography']
    
    Iterations = recon_para["Iterations"]
    Subsets = recon_para["Subsets"]
    OutputPixelSize = simind_para["OutputPixelSize"]
    OutputSliceWidth = simind_para["OutputSliceWidth"]
    FrameDurations = pbpk_para["FrameDurations"][0]
    
    roi_list = list(class_seg.keys())

    CalibrationFile = os.path.join(output_path_simind, 'calib.res')
    with open(CalibrationFile, 'r') as file:
        Sensitivity = float(file.readlines()[70].split(' ')[3]) #sensitivity in counts/s/MBq

    photopeak_path = os.path.join(output_path_simind, f'{output_name_simind}_{roi_list[0]}_0_tot_w2.h00')
    lower_path = os.path.join(output_path_simind, f'{output_name_simind}_{roi_list[0]}_0_tot_w1.h00')
    upper_path = os.path.join(output_path_simind, f'{output_name_simind}_{roi_list[0]}_0_tot_w3.h00')

    object_meta, proj_meta = simind.get_metadata(photopeak_path)

    with open(photopeak_path) as f:
        headerdata = f.readlines()
        
    headerdata = np.array(headerdata)
    proj_dim1 = get_header_value(headerdata, 'matrix size [1]', int)
    proj_dim2 = get_header_value(headerdata, 'matrix size [2]', int)

    num_proj = get_header_value(headerdata, 'total number of images', int)
    
    ww_peak, ww_lower, ww_upper = [simind.get_energy_window_width(path) for path in [photopeak_path, lower_path, upper_path]]
    
    for i, fr in enumerate(pbpk_para["FrameStartTimes"]):
        log.info(f"Reconstructing frame {i}")
        print(f'[RECON] Reconstructing frame {i}...')
        
        lower = np.fromfile(os.path.join(output_path_simind, f'{output_name_simind}_{fr}min_tot_w1.a00'),dtype=np.float32)
        photopeak = np.fromfile(os.path.join(output_path_simind, f'{output_name_simind}_{fr}min_tot_w2.a00'),dtype=np.float32)
        upper = np.fromfile(os.path.join(output_path_simind, f'{output_name_simind}_{fr}min_tot_w3.a00'),dtype=np.float32)

        lower = np.transpose(lower.reshape((num_proj,proj_dim2,proj_dim1))[:,::-1], (0,2,1))
        photopeak = np.transpose(photopeak.reshape((num_proj,proj_dim2,proj_dim1))[:,::-1], (0,2,1))
        upper = np.transpose(upper.reshape((num_proj,proj_dim2,proj_dim1))[:,::-1], (0,2,1))

        lower = torch.tensor(lower.copy()).to(pytomography.device)
        photopeak = torch.tensor(photopeak.copy()).to(pytomography.device)
        upper = torch.tensor(upper.copy()).to(pytomography.device)

        photopeak_realization = torch.poisson(photopeak) #counts
        lower_realization = torch.poisson(lower)
        upper_realization = torch.poisson(upper)

        scatter_estimate_TEW = simind.compute_EW_scatter(lower_realization, upper_realization , ww_lower, ww_upper, ww_peak)

        path_amap = os.path.join(output_path_simind, f'{output_name_simind}_{roi_list[0]}_0.hct') #change this line later
        amap = simind.get_attenuation_map(path_amap)
        att_transform = SPECTAttenuationTransform(amap)

        psf_meta = simind.get_psfmeta_from_header(photopeak_path)
        psf_transform = SPECTPSFTransform(psf_meta)

        system_matrix = SPECTSystemMatrix(
            obj2obj_transforms = [att_transform, psf_transform],
            proj2proj_transforms = [],
            object_meta = object_meta,
            proj_meta = proj_meta
        )

        likelihood = PoissonLogLikelihood(
        system_matrix = system_matrix,
        projections = photopeak_realization,
        additive_term = scatter_estimate_TEW
        )

        recon_algorithm = OSEM(likelihood)

        reconstructed_image = recon_algorithm(
        n_iters=Iterations,
        n_subsets=Subsets,
        )

        recon_img = sitk.GetImageFromArray(reconstructed_image.cpu().T/Sensitivity/FrameDurations/OutputPixelSize**2/OutputSliceWidth) #counts/(counts/s/MBq)/s/ -> MBq/ml (dont want depend on image size)
        recon_img.SetSpacing((OutputPixelSize,OutputPixelSize,OutputSliceWidth))
        sitk.WriteImage(recon_img, os.path.join(output_path_recon, f'{output_name_recon}_frame{i}.nii'), imageIO='NiftiImageIO')

    atn_img = sitk.GetImageFromArray(amap.cpu().T)
    atn_img.SetSpacing((OutputPixelSize,OutputPixelSize,OutputSliceWidth))
    sitk.WriteImage(atn_img, os.path.join(output_path_recon, f'{output_name_recon}_atn_img.nii'), imageIO='NiftiImageIO')

    return None