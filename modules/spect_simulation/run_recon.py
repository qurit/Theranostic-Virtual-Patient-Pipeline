"""
SPECT image reconstruction utilities for the TDT pipeline using PyTomography.

This module provides a wrapper around PyTomography to:

- Read SIMIND-generated projection data (photopeak and scatter windows).
- Generate noisy Poisson realizations of the projection data.
- Estimate scatter using the triple-energy window (TEW) method.
- Build a SPECT system matrix with attenuation and PSF transforms.
- Run OSEM reconstruction for each time frame.
- Save reconstructed SPECT images and attenuation maps as NIfTI files.
"""

import logging
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

logger = logging.getLogger(__name__)


def run_recon(
    recon_para,
    pbpk_para,
    simind_para,
    out_paths_simind,
    out_paths_recon,
    class_seg,
    activity_map_sum,
):
    """
    Run OSEM reconstruction on SIMIND projections for all time frames.

    Parameters
    ----------
    recon_para : dict
        Reconstruction-related parameters. Expected keys include:
            - "name" : str
            - "Iterations" : int
            - "Subsets" : int
    pbpk_para : dict
        PBPK-related parameters. Expected keys include:
            - "FrameStartTimes" : list[float]
            - "FrameDurations" : list[float]
              (only the first element is currently used to normalize
               reconstructed images by acquisition duration)
    simind_para : dict
        SIMIND-related parameters. Expected keys include:
            - "name" : str
            - "OutputPixelSize" : float
            - "OutputSliceWidth" : float
    out_paths_simind : str
        Directory containing SIMIND outputs (projections, attenuation, calib).
    out_paths_recon : str
        Directory where reconstructed images and attenuation maps are saved.
    class_seg : dict
        Mapping from region names to labels. Only the keys are used here to
        select example organ files for metadata and attenuation.
    activity_map_sum : numpy.ndarray
        Total activity per frame [MBq] across the whole volume. Currently not
        used inside this function but kept for interface compatibility.

    Returns
    -------
    recon_last_frame_path : str
        Filesystem path to the reconstructed image of the last frame.
    """
    output_name_simind = simind_para["name"]
    output_name_recon = recon_para["name"]
    output_path_simind = out_paths_simind
    output_path_recon = out_paths_recon

    iterations = recon_para["Iterations"]
    subsets = recon_para["Subsets"]
    output_pixel_size = simind_para["OutputPixelSize"]
    output_slice_width = simind_para["OutputSliceWidth"]
    frame_duration = pbpk_para["FrameDurations"][0]  # using first value

    roi_list = list(class_seg.keys())

    # ----- Calibration sensitivity -----
    calibration_file = os.path.join(output_path_simind, "calib.res")
    with open(calibration_file, "r") as file:
        lines = file.readlines()
        sensitivity_line = lines[70].strip()
        # Expect something like: "Sensitivity Cps/MBq: 10.6609"
        sensitivity = float(sensitivity_line.split(":")[-1].strip().split()[0])

    logger.debug("Parsed sensitivity from calibration file: %f Cps/MBq", sensitivity)

    # ----- Projection metadata -----
    photopeak_path = os.path.join(
        output_path_simind,
        f"{output_name_simind}_{roi_list[0]}_0_tot_w2.h00",
    )
    lower_path = os.path.join(
        output_path_simind,
        f"{output_name_simind}_{roi_list[0]}_0_tot_w1.h00",
    )
    upper_path = os.path.join(
        output_path_simind,
        f"{output_name_simind}_{roi_list[0]}_0_tot_w3.h00",
    )
    
    # ----- Cor metadata -----
    cor_path = os.path.join(
        output_path_simind,
        f"{output_name_simind}_{roi_list[0]}_0.cor",
    )
    cor_path_new = np.loadtxt(cor_path)[:, 0].astype(float)
    np.savetxt(f"{cor_path}", cor_path_new)

    object_meta, proj_meta = simind.get_metadata(photopeak_path,cor_path)

    with open(photopeak_path, "r") as f:
        headerdata = np.array(f.readlines())

    proj_dim1 = get_header_value(headerdata, "matrix size [1]", int)
    proj_dim2 = get_header_value(headerdata, "matrix size [2]", int)
    num_proj = get_header_value(headerdata, "total number of images", int)

    ww_peak, ww_lower, ww_upper = [
        simind.get_energy_window_width(path)
        for path in (photopeak_path, lower_path, upper_path)
    ]

    logger.debug(
        "Projection dimensions: (%d, %d), num_proj=%d",
        proj_dim1,
        proj_dim2,
        num_proj,
    )
    logger.debug(
        "Energy window widths (peak, lower, upper): (%s, %s, %s)",
        ww_peak,
        ww_lower,
        ww_upper,
    )

    recon_last_frame_path = None

    # ----- Frame loop -----
    for i, frame_start_time in enumerate(pbpk_para["FrameStartTimes"]):
        logger.info("Reconstructing frame %d (start time = %s min)", i, frame_start_time)

        # Load SIMIND projection data (counts)
        lower = np.fromfile(
            os.path.join(
                output_path_simind,
                f"{output_name_simind}_{frame_start_time}min_tot_w1.a00",
            ),
            dtype=np.float32,
        )
        photopeak = np.fromfile(
            os.path.join(
                output_path_simind,
                f"{output_name_simind}_{frame_start_time}min_tot_w2.a00",
            ),
            dtype=np.float32,
        )
        upper = np.fromfile(
            os.path.join(
                output_path_simind,
                f"{output_name_simind}_{frame_start_time}min_tot_w3.a00",
            ),
            dtype=np.float32,
        )

        # Reshape to (num_proj, y, x), flip y, then permute to (num_proj, x, y)
        lower = np.transpose(
            lower.reshape((num_proj, proj_dim2, proj_dim1))[:, ::-1],
            (0, 2, 1),
        )
        photopeak = np.transpose(
            photopeak.reshape((num_proj, proj_dim2, proj_dim1))[:, ::-1],
            (0, 2, 1),
        )
        upper = np.transpose(
            upper.reshape((num_proj, proj_dim2, proj_dim1))[:, ::-1],
            (0, 2, 1),
        )

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

        # Attenuation map (from a representative organ file)
        path_amap = os.path.join(
            output_path_simind,
            f"{output_name_simind}_{roi_list[0]}_0.hct",  # TODO: generalize if needed
        )
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
            proj_meta=proj_meta,
        )

        # Likelihood and reconstruction algorithm
        likelihood = PoissonLogLikelihood(
            system_matrix=system_matrix,
            projections=photopeak_realization,
            additive_term=scatter_estimate_tew,
        )

        recon_algorithm = OSEM(likelihood)

        reconstructed_image = recon_algorithm(
            n_iters=iterations,
            n_subsets=subsets,
        )

        # Convert to MBq/mL:
        # counts / (counts/s/MBq) / s / cm^2 / cm -> MBq/mL
        recon_img_arr = (
            reconstructed_image.cpu().T
            / sensitivity
            / frame_duration
            / (output_pixel_size**2)
            / output_slice_width
        )

        recon_img = sitk.GetImageFromArray(recon_img_arr)
        recon_img.SetSpacing(
            (output_pixel_size, output_pixel_size, output_slice_width)
        )

        recon_frame_path = os.path.join(
            output_path_recon,
            f"{output_name_recon}_frame{i}.nii",
        )
        sitk.WriteImage(
            recon_img,
            recon_frame_path,
            imageIO="NiftiImageIO",
        )

        logger.debug("Saved reconstructed frame %d to: %s", i, recon_frame_path)
        recon_last_frame_path = recon_frame_path

    # Save attenuation image as NIfTI (same spacing as reconstruction)
    atn_img = sitk.GetImageFromArray(amap.cpu().T)
    atn_img.SetSpacing(
        (output_pixel_size, output_pixel_size, output_slice_width)
    )
    atn_path = os.path.join(
        output_path_recon,
        f"{output_name_recon}_atn_img.nii",
    )
    sitk.WriteImage(
        atn_img,
        atn_path,
        imageIO="NiftiImageIO",
    )
    logger.debug("Saved attenuation image to: %s", atn_path)

    return recon_last_frame_path
