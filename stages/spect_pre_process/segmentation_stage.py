import os
import dicom2nifti
import torch
import SimpleITK as sitk
from torch.cuda.amp import GradScaler as _GradScaler
from totalsegmentator.python_api import totalsegmentator

torch.GradScaler = _GradScaler


class TotalSegmentationStage:
    def __init__(self, config, context):
        self.config = config
        self.context = context

        self.ct_input_path = config["ct_input"]["path1"]
        self.device = config["spect_preprocessing"]["device"]
        self.roi_subset = config["spect_preprocessing"]["roi_subset"]
        self.ml = config["spect_preprocessing"]["ml"]

        subdir_name = config["subdir_names"]["spect_preprocessing"]
        output_root = config["output_folder"]["title"]
        self.output_dir = os.path.join(output_root, subdir_name)
        os.makedirs(self.output_dir, exist_ok=True)

        self.prefix = config["spect_preprocessing"]["name"]

        self.ct_nii_path = None
        self.roi_seg_path = None
        self.body_seg_dir = None

    def _standardize_ct_to_nifti(self):
        self.ct_nii_path = os.path.join(self.output_dir, f"{self.prefix}_ct.nii.gz")

        if os.path.exists(self.ct_nii_path):
            return

        if os.path.isdir(self.ct_input_path):
            dicom2nifti.dicom_series_to_nifti(
                self.ct_input_path,
                self.ct_nii_path,
                reorient_nifti=True,
            )
        else:
            lower_input = self.ct_input_path.lower()
            if lower_input.endswith((".nii", ".nii.gz")):
                sitk.WriteImage(sitk.ReadImage(self.ct_input_path), self.ct_nii_path, True)
            else:
                raise ValueError(
                    "Unsupported CT input. Provide a DICOM folder or a NIfTI file "
                    f"(.nii/.nii.gz). Got: {self.ct_input_path}"
                )

    def run(self):
        if not self.ml:
            raise ValueError(
                "This pipeline expects a single multilabel ROI seg NIfTI. "
                "Set config['spect_preprocessing']['ml'] = true."
            )

        self._standardize_ct_to_nifti()

        self.roi_seg_path = os.path.join(self.output_dir, f"{self.prefix}_roi_seg.nii.gz")
        self.body_seg_dir = os.path.join(self.output_dir, f"{self.prefix}_body_seg_dir")
        os.makedirs(self.body_seg_dir, exist_ok=True)

        roi_done = os.path.exists(self.roi_seg_path)
        body_done = os.path.exists(os.path.join(self.body_seg_dir, "body.nii.gz"))

        if not roi_done:
            totalsegmentator(
                self.ct_nii_path,
                self.roi_seg_path,
                device=self.device,
                ml=self.ml,
                roi_subset=self.roi_subset,
            )

        if not body_done:
            totalsegmentator(
                self.ct_nii_path,
                self.body_seg_dir,
                device=self.device,
                task="body",
            )

        if not os.path.exists(self.roi_seg_path):
            raise FileNotFoundError(f"ROI seg not found: {self.roi_seg_path}")

        if not os.path.exists(os.path.join(self.body_seg_dir, "body.nii.gz")):
            raise FileNotFoundError(f"Body seg not found in: {self.body_seg_dir}")

        self.context.ct_nii_path = self.ct_nii_path
        self.context.roi_seg_path = self.roi_seg_path
        self.context._body_seg_dir = self.body_seg_dir
        
        self.context.extras["totseg_output_dir"] = self.output_dir
        self.context.extras["totseg_cached"] = {"ct": os.path.exists(self.ct_nii_path),
                                                "roi": os.path.exists(self.roi_seg_path),
                                                "body": os.path.exists(os.path.join(self.body_seg_dir, "body.nii.gz"))}


        return self.context
