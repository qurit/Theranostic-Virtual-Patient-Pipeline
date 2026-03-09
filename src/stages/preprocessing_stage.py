"""
SIMIND Preprocessing Stage for the TDT pipeline.  

This stage prepares inputs needed by SIMIND by:
- Converting CT + segmentation NIfTIs into the SIMIND grid convention (z, y, x with y-flip).
- Optionally resizing to a cubic grid using a single isotropic scaling factor.
- Building ROI masks and a label->name class map from a unified TDT multilabel segmentation.
- Writing binary files used by SIMIND (attenuation map, ROI segmentation, body mask).

Expected Context interface
--------------------------
Incoming `context` is expected to provide:
- context.subdir_paths["spect_preprocessing"] : str
- context.config["spect_preprocessing"] with keys:
    - "name" : str (prefix for outputs)
    - "xy_dim" : int | None (target dimension for x/y; scaling is isotropic)
    - "roi_subset" : str | list[str] (TDT ROI names, e.g., "kidney", "liver", ...)
- context.ct_nii_path : str
- context.body_ml_path : str
- context.tdt_roi_seg_path : str  (unified TDT ROI segmentation NIfTI)

On success, this stage sets:
- context.body_seg_arr : np.ndarray (float32 mask in SIMIND grid)
- context.roi_body_seg_arr : np.ndarray (int16 labels in SIMIND grid; filtered to requested ROIs)
- context.mask_roi_body : dict[int, np.ndarray] (label_id -> boolean mask)
- context.class_seg : dict[str, int] (roi_name -> label_id present in the filtered seg)
- context.atn_av_path : str (path to attenuation binary)
- context.arr_px_spacing_cm : tuple[float, float, float] (z,y,x spacing in cm)
- context.arr_shape_new : tuple[int, int, int] (z,y,x array shape)

Maintainer / contact: pyazdi@bccrc.ca  
""" 

from __future__ import annotations  

import os
import json
from typing import Any, Dict, Optional, Sequence, Tuple, Union  

import nibabel as nib
import numpy as np
from scipy.ndimage import zoom
from json_minify import json_minify


class SimindPreprocessStage:
    """
    Prepare CT + ROI masks for SIMIND simulation.

    Notes
    -----
    - Grid convention: arrays are converted to (z, y, x) and flipped in `y` to match SIMIND's expected orientation.
    - If `xy_dim` is provided, an isotropic scaling is applied to (z, y, x). This keeps
      voxels cubic *in index space* after scaling (not necessarily physical isotropy).
    """  

    def __init__(self, context: Any) -> None:  
        self.context = context

        self.output_dir: str = context.subdir_paths["spect_preprocessing"]  
        os.makedirs(self.output_dir, exist_ok=True)

        self.prefix: str = context.config["spect_preprocessing"]["name"]  
        self.resize: Optional[int] = context.config["spect_preprocessing"]["xy_dim"]  

        # ---- ROI subset is now TDT-level names ----
        roi_subset = context.config["spect_preprocessing"]["roi_subset"]
        if isinstance(roi_subset, str):
            roi_subset = [roi_subset]
        self.roi_subset: Sequence[str] = [str(r).strip() for r in roi_subset if str(r).strip()]  

        # ---- Load ONLY TDT pipeline label map (name -> id) ----
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
        self.ts_map_path: str = os.path.join(repo_root, "data", "tdt_map.json")  
        if not os.path.exists(self.ts_map_path):
            raise FileNotFoundError(f"Class map json not found: {self.ts_map_path}")

        with open(self.ts_map_path, encoding="utf-8") as f:
            ts_map_json = json.loads(json_minify(f.read()))

        # name -> int_label, e.g. {"body":1,"kidney":2,...}
        self.tdt_name2id: Dict[str, int] = {name: int(lab) for lab, name in ts_map_json["TDT_Pipeline"].items()}  

        self.ct_nii_path: Optional[str] = context.ct_nii_path 
        self.body_ml_path: Optional[str] = context.body_ml_path  
        self.tdt_roi_seg_path: Optional[str] = context.tdt_roi_seg_path  

    # -----------------------------
    # helpers
    # -----------------------------
    @staticmethod
    def _build_class_map(seg_arr: np.ndarray, id_to_name: Dict[int, str]) -> Dict[str, int]:  
        """
        Return {roi_name: label_id} for labels present in `seg_arr`.

        Parameters
        ----------
        seg_arr : np.ndarray
            Multilabel segmentation array (labels as integers).
        id_to_name : dict[int, str]
            Mapping from label_id -> roi_name.

        Returns
        -------
        dict[str, int]
            ROI names present in `seg_arr` mapped to their label IDs.
        """  
        class_map: Dict[str, int] = {}
        labels = np.unique(seg_arr.astype(int))
        for lab in labels:
            if lab == 0:
                continue
            name = id_to_name.get(int(lab))
            if name is not None:
                class_map[name] = int(lab)
        return class_map

    @staticmethod
    def _build_label_masks(arr: np.ndarray) -> Dict[int, np.ndarray]:  
        """
        Build boolean masks for each non-zero label in a multilabel segmentation.

        Parameters
        ----------
        arr : np.ndarray
            Segmentation array with integer labels (0 is background).

        Returns
        -------
        dict[int, np.ndarray]
            Mapping label_id -> boolean mask.
        """  
        labels = np.unique(arr)
        labels = labels[labels != 0]

        if labels.size == 0:
            raise ValueError(
                "Segmentation has no non-zero labels (only background=0). "
                "Segmentation likely failed or ROI subset is empty/mismatched."
            )
        return {int(lab): (arr == lab) for lab in labels}

    @staticmethod
    def _hu_to_mu(
        hu_arr: np.ndarray,
        pixel_size_cm: float,
        mu_water: float = 0.1537,
        mu_bone: float = 0.2234,
    ) -> np.ndarray: 
        """
        Convert HU (Hounsfield Units) CT values to a linear attenuation map (mu) per pixel.

        Parameters
        ----------
        hu_arr : np.ndarray
            CT array in HU (float32 recommended).
        pixel_size_cm : float
            Effective pixel size (cm) used to scale mu to per-pixel units.
        mu_water : float, default=0.1537
            Linear attenuation coefficient of water (1/cm).
        mu_bone : float, default=0.2234
            Linear attenuation coefficient of bone (1/cm).

        Returns
        -------
        np.ndarray
            Mu map in "per pixel" units (scaled by `pixel_size_cm`), float32.
        """  
        mu_water_pixel = mu_water * pixel_size_cm
        mu_bone_pixel = mu_bone * pixel_size_cm

        mu_map = np.zeros_like(hu_arr, dtype=np.float32)

        soft = hu_arr <= 0
        bone = hu_arr > 0

        mu_map[soft] = mu_water_pixel * (1 + hu_arr[soft] / 1000.0)
        mu_map[bone] = mu_water_pixel + (hu_arr[bone] / 1000.0) * (mu_bone_pixel - mu_water_pixel)

        return mu_map

    def _write_attenuation_bin(
        self,
        ct_arr: np.ndarray,
        body_seg_arr: np.ndarray,
        pixel_size_cm: float,
        filename: str = "spect_preprocessing_atn_av.bin",
    ) -> str:  
        """
        Write the attenuation map binary used by SIMIND.

        Parameters
        ----------
        ct_arr : np.ndarray
            CT array on SIMIND grid (float32).
        body_seg_arr : np.ndarray
            Body mask on SIMIND grid (float32 or bool-like).
        pixel_size_cm : float
            Effective pixel size (cm) used for HU->mu conversion.
        filename : str
            Output filename inside `self.output_dir`.

        Returns
        -------
        str
            Absolute output path to the written `.bin` file.
        """  
        mu_map = self._hu_to_mu(np.asarray(ct_arr, dtype=np.float32), pixel_size_cm)
        mu_map *= body_seg_arr  # assumes mask is 0/1

        out_path = os.path.join(self.output_dir, filename)
        mu_map.tofile(out_path)
        return out_path

    def _filter_to_requested_rois(self, roi_seg_arr: np.ndarray) -> np.ndarray: 
        """
        Keep only labels requested in config (plus body if present); zero out the rest.

        This ensures downstream masks/class_map reflect *only* user-selected ROIs.

        Parameters
        ----------
        roi_seg_arr : np.ndarray
            Unified TDT ROI segmentation (integer labels).

        Returns
        -------
        np.ndarray
            Filtered segmentation array (same shape), with non-requested labels set to 0.
        """  
        requested = set(self.roi_subset)
        keep_ids = set()

        # Always keep body label if it exists in the unified seg
        if "body" in self.tdt_name2id:
            keep_ids.add(self.tdt_name2id["body"])

        for name in requested:
            lab = self.tdt_name2id.get(name)
            if lab is None:
                raise ValueError(f"Requested ROI '{name}' not in TDT label map.")
            keep_ids.add(lab)

        keep_ids.discard(0)

        out = roi_seg_arr.copy()
        out[~np.isin(out, list(keep_ids))] = 0
        return out

    @staticmethod
    def _to_simind_grid(
        nii_obj: nib.Nifti1Image,
        resize: Optional[int] = None,
        transpose_tuple: Tuple[int, int, int] = (2, 1, 0),
        zoom_order: int = 0,
    ) -> Tuple[np.ndarray, float]: 
        """
        Convert NIfTI object to SIMIND grid format with optional resizing.

        Parameters
        ----------
        nii_obj : nib.Nifti1Image
            Loaded NIfTI image object.
        resize : Optional[int]
            Target in-plane dimension (expects square y=x after transpose). If provided,
            applies isotropic zoom to (z,y,x) using a single scale factor.
        transpose_tuple : tuple[int, int, int]
            Permutation applied to convert input array to (z,y,x).
        zoom_order : int
            Interpolation order for `scipy.ndimage.zoom`:
            - 0: nearest neighbor (seg masks)
            - 1: linear (CT intensities)

        Returns
        -------
        tuple[np.ndarray, float]
            - Array in SIMIND grid convention (z,y,x) with y-flip applied.
            - Scale factor applied (1.0 if no resize).
        """  
        arr = np.array(nii_obj.get_fdata(dtype=np.float32))
        arr = np.transpose(arr, transpose_tuple)[:, ::-1, :]

        scale = 1.0
        if resize is not None:
            if arr.shape[1] != arr.shape[2]:
                raise ValueError("Resize parameter requires square in-plane dimensions (x=y).")
            scale = resize / arr.shape[1]
            arr = zoom(arr, (scale, scale, scale), order=zoom_order)

        return arr, scale

    # -----------------------------
    # main
    # -----------------------------
    def run(self) -> Any:  
        """
        Execute preprocessing and write SIMIND-ready binaries.

        Returns
        -------
        context : Context-like
            Updated context object with SIMIND arrays, masks, and file paths populated.
        """  
        # ---- existence checks ----
        if self.ct_nii_path is None or not os.path.exists(self.ct_nii_path):
            raise FileNotFoundError(f"ct_nii_path not found: {self.ct_nii_path}")

        if self.body_ml_path is None or not os.path.exists(self.body_ml_path):
            raise FileNotFoundError(f"Body segmentation not found: {self.body_ml_path}")

        if self.tdt_roi_seg_path is None or not os.path.exists(self.tdt_roi_seg_path):
            raise FileNotFoundError(f"Unified TDT ROI seg not found: {self.tdt_roi_seg_path}")

        # ---- load nifti ----
        ct_nii = nib.load(self.ct_nii_path)
        body_nii = nib.load(self.body_ml_path)
        roi_nii = nib.load(self.tdt_roi_seg_path)

        # ---- to simind grid ----
        ct_arr, scale = self._to_simind_grid(ct_nii, resize=self.resize, zoom_order=1)  # CT: linear
        body_arr, _ = self._to_simind_grid(body_nii, resize=self.resize, zoom_order=0)  # mask: body (nearest)
        roi_arr, _ = self._to_simind_grid(roi_nii, resize=self.resize, zoom_order=0)  # mask: unified ROI (nearest)

        # ---- sanity: body should be binary-ish ----
        body_mask = (body_arr > 0).astype(np.float32)
        roi_arr = roi_arr.astype(np.int16)

        # ---- keep only requested TDT ROIs ----
        roi_arr = self._filter_to_requested_rois(roi_arr)

        # ---- build masks + class map from TDT label space ----
        masks = self._build_label_masks(roi_arr)
        id_to_name = {v: k for k, v in self.tdt_name2id.items()}
        class_seg = self._build_class_map(roi_arr, id_to_name)

        # ---- spacing ----
        # `ct_nii.header.get_zooms()` is in mm. After zoom, spacing shrinks by `scale`.
        zooms_mm = np.array(ct_nii.header.get_zooms()[:3], dtype=float) / scale
        zooms_mm = zooms_mm[[2, 1, 0]]  # (z,y,x)
        arr_px_spacing_cm = tuple(float(x) * 0.1 for x in zooms_mm)

        # Pixel size for HU->mu conversion: mean of in-plane (y,x) spacing in cm
        pixel_size_cm = (arr_px_spacing_cm[1] + arr_px_spacing_cm[2]) / 2.0

        atn_av_path = self._write_attenuation_bin(
            ct_arr,
            body_mask,  # use binary mask
            pixel_size_cm=pixel_size_cm,
        )

        # ---- write SIMIND bins ----
        roi_bin = os.path.join(self.output_dir, f"{self.prefix}_roi_seg.bin")
        body_bin = os.path.join(self.output_dir, f"{self.prefix}_body_seg.bin")
        roi_body_bin = os.path.join(self.output_dir, f"{self.prefix}_roi_body_seg.bin")

        # ROI labels masked to body to prevent non-body label artifacts in SIMIND inputs
        roi_arr.astype(np.float32).tofile(roi_bin)
        body_mask.astype(np.float32).tofile(body_bin)
        (roi_arr * body_mask).astype(np.float32).tofile(roi_body_bin)

        # ---- update context ----
        self.context.body_seg_arr = body_mask
        self.context.roi_body_seg_arr = roi_arr * body_mask

        self.context.mask_roi_body = masks
        self.context.class_seg = class_seg

        self.context.atn_av_path = atn_av_path

        self.context.arr_px_spacing_cm = arr_px_spacing_cm
        self.context.arr_shape_new = ct_arr.shape

        return self.context