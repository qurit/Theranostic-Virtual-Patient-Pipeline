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
from pytomography.transforms.SPECT import SPECTAttenuationTransform, SPECTPSFTransform


class SpectReconstructionStage:
    def __init__(self, config, context):
        self.config = config
        self.context = context

        # output_dir: where frame totals + calib.res are 
        self.output_dir = getattr(context, "spect_sim_output_dir", None)
        if self.output_dir is None:
            self.output_dir = context.extras.get("simind_output_dir", None)

        # fallback (if user runs recon standalone)
        if self.output_dir is None:
            subdir_name = config["subdir_names"]["spect_simulation"]
            output_root = config["output_folder"]["title"]
            self.output_dir = os.path.join(output_root, subdir_name)

        os.makedirs(self.output_dir, exist_ok=True)

        # header_dir: where SIMIND wrote .h00/.cor/.hct (work_dir)
        self.header_dir = context.extras.get("simind_work_dir", self.output_dir)

        self.prefix = config["spect_simulation"]["name"]

        self.frame_start = config["pbpk"]["FrameStartTimes"]
        self.frame_duration = config["pbpk"]["FrameDurations"][0]

        self.iterations = config["spect_simulation"]["Iterations"]
        self.subsets = config["spect_simulation"]["Subsets"]
        self.output_pixel_width = config["spect_simulation"]["OutputPixelWidth"]
        self.output_slice_width = config["spect_simulation"]["OutputSliceWidth"]

        self.output_tuple = (self.output_pixel_width, self.output_pixel_width, self.output_slice_width)

    # -----------------------------
    # helpers
    # -----------------------------
    def _get_sensitivity_from_calibration_file(self, calibration_file):
        sensitvity_line = 70
        with open(calibration_file, "r") as file:
            lines = file.readlines()
            sensitivity_line = lines[sensitvity_line].strip()
            sensitivity = float(sensitivity_line.split(":")[-1].strip().split()[0])
        return sensitivity

    def _get_object_and_proj_metadata(self, photopeak_path, cor_path):
        object_meta, proj_meta = simind.get_metadata(photopeak_path, cor_path)
        return object_meta, proj_meta

    def _get_metadata_from_header(self, photopeak_path, lower_path, upper_path):
        with open(photopeak_path, "r") as f:
            headerdata = np.array(f.readlines())

        proj_dim1 = get_header_value(headerdata, "matrix size [1]", int)
        proj_dim2 = get_header_value(headerdata, "matrix size [2]", int)
        num_proj = get_header_value(headerdata, "total number of images", int)

        ww_peak, ww_lower, ww_upper = [
            simind.get_energy_window_width(path) for path in (photopeak_path, lower_path, upper_path)
        ]
        return proj_dim1, proj_dim2, num_proj, ww_peak, ww_lower, ww_upper

    def _get_cor_data(self, cor_path):
        cor_path_new = np.loadtxt(cor_path).astype(float)
        if cor_path_new.ndim == 2:
            cor_path_new = cor_path_new[:, 0]
            np.savetxt(cor_path, cor_path_new)
        return cor_path_new

    def _expected_proj_len(self, proj_dim1, proj_dim2, num_proj):
        return int(num_proj) * int(proj_dim1) * int(proj_dim2)

    def _reshape_projection_data(self, projection, proj_dim1, proj_dim2, num_proj):
        return np.transpose(
            projection.reshape((num_proj, proj_dim2, proj_dim1))[:, ::-1],
            (0, 2, 1),
        )

    def _convert_counts_to_mbq_per_ml(self, reconstructed_image, sensitivity, frame_duration,
                                     output_pixel_width, output_slice_width):
        return (
            reconstructed_image.cpu().T
            / sensitivity
            / frame_duration
            / (output_pixel_width ** 2)
            / output_slice_width
        )

    def _get_recon_img(self, likelihood, sensitivity):
        recon_algorithm = OSEM(likelihood)
        reconstructed_image = recon_algorithm(
            n_iters=self.iterations,
            n_subsets=self.subsets,
        )

        recon_img_arr = self._convert_counts_to_mbq_per_ml(
            reconstructed_image,
            sensitivity,
            self.frame_duration * 60,  # min -> sec
            self.output_pixel_width,
            self.output_slice_width,
        )

        recon_img = sitk.GetImageFromArray(recon_img_arr)
        recon_img.SetSpacing(self.output_tuple)
        return recon_img

    def _write_recon_atn_img(self, amap):
        recon_atn_path = os.path.join(self.output_dir, f"{self.prefix}_atn_img.nii")

        if os.path.exists(recon_atn_path):
            recon_atn_img = sitk.ReadImage(recon_atn_path)
            return recon_atn_img, recon_atn_path

        recon_atn_img = sitk.GetImageFromArray(amap.cpu().T)
        recon_atn_img.SetSpacing(self.output_tuple)
        sitk.WriteImage(recon_atn_img, recon_atn_path, imageIO="NiftiImageIO")
        return recon_atn_img, recon_atn_path

    # -----------------------------
    # main
    # -----------------------------
    def run(self):
        self.context.require("class_seg")

        roi_list = list(self.context.class_seg.keys())
        
        # sensitivity comes from output_dir (calibration runs in output_dir)
        calibration_file = os.path.join(self.output_dir, "calib.res")
        if not os.path.exists(calibration_file):
            raise FileNotFoundError(f"Calibration file not found: {calibration_file}")
        sensitivity = self._get_sensitivity_from_calibration_file(calibration_file)

        # headers come from header_dir (work_dir)
        photopeak_h = os.path.join(self.header_dir, f"{self.prefix}_{roi_list[0]}_0_tot_w2.h00")
        lower_h = os.path.join(self.header_dir, f"{self.prefix}_{roi_list[0]}_0_tot_w1.h00")
        upper_h = os.path.join(self.header_dir, f"{self.prefix}_{roi_list[0]}_0_tot_w3.h00")

        for p in (photopeak_h, lower_h, upper_h):
            if not os.path.exists(p):
                raise FileNotFoundError(f"Missing SIMIND header file: {p}")

        cor_path = os.path.join(self.header_dir, f"{self.prefix}_{roi_list[0]}_0.cor")
        if not os.path.exists(cor_path):
            raise FileNotFoundError(f"Missing COR file: {cor_path}")
        _ = self._get_cor_data(cor_path)

        object_meta, proj_meta = self._get_object_and_proj_metadata(photopeak_h, cor_path)

        proj_dim1, proj_dim2, num_proj, ww_peak, ww_lower, ww_upper = self._get_metadata_from_header(
            photopeak_h, lower_h, upper_h
        )
        expected_len = self._expected_proj_len(proj_dim1, proj_dim2, num_proj)

        path_amap = os.path.join(self.header_dir, f"{self.prefix}_{roi_list[0]}_0.hct")
        if not os.path.exists(path_amap):
            raise FileNotFoundError(f"Missing attenuation map header (.hct): {path_amap}")

        amap = simind.get_attenuation_map(path_amap)

        att_transform = SPECTAttenuationTransform(amap)
        psf_meta = simind.get_psfmeta_from_header(photopeak_h)
        psf_transform = SPECTPSFTransform(psf_meta)

        system_matrix = SPECTSystemMatrix(
            obj2obj_transforms=[att_transform, psf_transform],
            proj2proj_transforms=[],
            object_meta=object_meta,
            proj_meta=proj_meta,
        )

        recon_paths = []
        for time in self.frame_start:
            recon_output_path = os.path.join(self.output_dir, f"{self.prefix}_{time}min.nii")

            if os.path.exists(recon_output_path):
                recon_paths.append(recon_output_path)
                continue

            # frame totals are in output_dir
            lower_a00 = os.path.join(self.output_dir, f"{self.prefix}_{time}min_tot_w1.a00")
            photopeak_a00 = os.path.join(self.output_dir, f"{self.prefix}_{time}min_tot_w2.a00")
            upper_a00 = os.path.join(self.output_dir, f"{self.prefix}_{time}min_tot_w3.a00")

            for p in (lower_a00, photopeak_a00, upper_a00):
                if not os.path.exists(p):
                    raise FileNotFoundError(f"Missing projection file: {p}")

            lower = np.fromfile(lower_a00, dtype=np.float32)
            photopeak = np.fromfile(photopeak_a00, dtype=np.float32)
            upper = np.fromfile(upper_a00, dtype=np.float32)

            if lower.size != expected_len:
                raise ValueError(f"Projection size mismatch for lower: got {lower.size}, expected {expected_len}")
            if photopeak.size != expected_len:
                raise ValueError(f"Projection size mismatch for photopeak: got {photopeak.size}, expected {expected_len}")
            if upper.size != expected_len:
                raise ValueError(f"Projection size mismatch for upper: got {upper.size}, expected {expected_len}")

            lower = self._reshape_projection_data(lower, proj_dim1, proj_dim2, num_proj)
            photopeak = self._reshape_projection_data(photopeak, proj_dim1, proj_dim2, num_proj)
            upper = self._reshape_projection_data(upper, proj_dim1, proj_dim2, num_proj)

            lower = torch.tensor(lower.copy()).to(pytomography.device)
            photopeak = torch.tensor(photopeak.copy()).to(pytomography.device)
            upper = torch.tensor(upper.copy()).to(pytomography.device)

            photopeak_realization = torch.poisson(photopeak)
            lower_realization = torch.poisson(lower)
            upper_realization = torch.poisson(upper)

            scatter_estimate_tew = simind.compute_EW_scatter(
                lower_realization,
                upper_realization,
                ww_lower,
                ww_upper,
                ww_peak,
            )

            likelihood = PoissonLogLikelihood(
                system_matrix=system_matrix,
                projections=photopeak_realization,
                additive_term=scatter_estimate_tew,
            )

            recon_img = self._get_recon_img(likelihood, sensitivity)
            sitk.WriteImage(recon_img, recon_output_path, imageIO="NiftiImageIO")
            recon_paths.append(recon_output_path)

        recon_atn_img, recon_atn_path = self._write_recon_atn_img(amap)

        # Update context
        self.context.spect_sim_output_dir = self.output_dir
        self.context.recon_paths = recon_paths
        self.context.recon_atn_img = recon_atn_img
        self.context.recon_atn_path = recon_atn_path

        self.context.extras["recon_output_dir"] = self.output_dir
        self.context.extras["recon_header_dir_used"] = self.header_dir

        return self.context
