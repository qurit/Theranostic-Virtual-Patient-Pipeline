import os
import torch
from pytomography.io.SPECT import simind
from pytomography.projectors.SPECT import SPECTSystemMatrix
from pytomography.transforms.SPECT import SPECTAttenuationTransform, SPECTPSFTransform
from pytomography.algorithms import OSEM
from pytomography.likelihoods import PoissonLogLikelihood
import SimpleITK as sitk


def runRECON(config, output_name, output_path, FrameAct):

    Iterations = config["Recon"]["Iterations"]
    Subsets = config["Recon"]["Subsets"]
    OutputPixelSize = config["SIMIND"]["OutputPixelSize"]
    OutputSliceWidth = config["SIMIND"]["OutputSliceWidth"]
    FrameDurations = config["PBPK"]["FrameDurations"]

    CalibrationFile = os.path.join(output_path, 'calib.res')
    with open(CalibrationFile, 'r') as file:
        Sensitivity = float(file.readlines()[70].split(' ')[3])

    for i in range(len(config["PBPK"]["FrameStartTimes"])):
        print(f'Reconstructing frame {i}...')

        photopeak_path = os.path.join(output_path, f'{output_name}_frame{i}_tot_w2.h00')
        lower_path = os.path.join(output_path, f'{output_name}_frame{i}_tot_w1.h00')
        upper_path = os.path.join(output_path, f'{output_name}_frame{i}_tot_w3.h00')

        object_meta, proj_meta = simind.get_metadata(photopeak_path)

        photopeak = simind.get_projections(photopeak_path)
        lower = simind.get_projections(lower_path)
        upper = simind.get_projections(upper_path)

        photopeak_realization = torch.poisson(photopeak * FrameAct[i] * FrameDurations[i])
        lower_realization = torch.poisson(lower * FrameAct[i] * FrameDurations[i])
        upper_realization = torch.poisson(upper * FrameAct[i] * FrameDurations[i])

        ww_peak, ww_lower, ww_upper = [simind.get_energy_window_width(path) for path in [photopeak_path, lower_path, upper_path]]
        scatter_estimate_TEW = simind.compute_EW_scatter(lower_realization, upper_realization , ww_lower, ww_upper, ww_peak)

        path_amap = os.path.join(output_path, f'{output_name}_frame0.hct')
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

        recon_img = sitk.GetImageFromArray(reconstructed_image.cpu().T/Sensitivity/FrameDurations[i]/OutputPixelSize**2/OutputSliceWidth)
        recon_img.SetSpacing((OutputPixelSize,OutputPixelSize,OutputSliceWidth))
        sitk.WriteImage(recon_img, os.path.join(output_path, f'recon_frame{i}.nii'), imageIO='NiftiImageIO')

    atn_img = sitk.GetImageFromArray(amap.cpu().T)
    atn_img.SetSpacing((OutputPixelSize,OutputPixelSize,OutputSliceWidth))
    sitk.WriteImage(atn_img, os.path.join(output_path, f'atn_img.nii'), imageIO='NiftiImageIO')

    return None