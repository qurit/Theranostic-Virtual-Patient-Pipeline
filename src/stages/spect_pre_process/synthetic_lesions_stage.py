"""
Synthetic Lesions Stage for the TDT pipeline.

Goal
----
Generate synthetic spherical lesions inside user-specified organ ROIs (from the unified TDT ROI seg),
then insert them into the unified segmentation as the "synthetic_lesion" label (from tdt_map.json).

Key behavior (matches your notebook logic)
------------------------------------------
Constraints:
- Center must be inside ROI.
- Sphere must remain inside ROI (enforced using distance transform boundary constraint).
- Lesions must not overlap (physical distance in mm).
Sampling:
- prob="uniform": uniform candidate sampling
- prob="gaussian": Gaussian weights centered at ROI centroid
- prob="user_defined": user provides centers_zyx explicitly (validated)
- prob="tom": TODO (not implemented; raises)

Outputs
-------
Writes into:
  <spect_preprocessing_outputs>/<synthetic_lesions_name>_outputs/

- Backup of unified seg BEFORE lesions:
    tdt_roi_seg_pre_lesions.nii.gz
- Global lesion masks:
    all_lesions_binary.nii.gz   (uint8 0/1)
    all_lesions_labels.nii.gz   (int16 0=bg, 1..K=lesion id across ALL ROIs)
- Per-ROI outputs (for QC):
    <roi>/<roi>_lesions_labels.nii.gz
    <roi>/<roi>_lesions_binary.nii.gz
    <roi>/<roi>_organ_minus_lesions.nii.gz
    <roi>/<roi>_lesion_metadata.json

Most important side-effect
--------------------------
OVERWRITES `context.tdt_roi_seg_path` on disk so that:
- organ voxels remain their organ label
- lesion voxels become label = TDT_Pipeline["synthetic_lesion"] (e.g. 8)

Expected Context interface
--------------------------
Incoming `context` must provide:
- context.subdir_paths["spect_preprocessing"]
- context.config["synthetic_lesions"] with:
    - "name": str
    - "specs": dict | None
- context.tdt_roi_seg_path: str (unified multilabel seg produced by TdtRoiUnifyStage)

Maintainer / contact: pyazdi@bccrc.ca
"""

from __future__ import annotations

import os
import json
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import nibabel as nib
from scipy.ndimage import distance_transform_edt
from json_minify import json_minify


class SyntheticLesionsStage:
    def __init__(self, context: Any) -> None:
        self.context = context

        self.output_dir: str = context.subdir_paths["spect_preprocessing"]
        os.makedirs(self.output_dir, exist_ok=True)

        self.cfg: Dict[str, Any] = context.config.get("synthetic_lesions", {})
        self.prefix: str = str(self.cfg.get("name", "synthetic_lesions"))
        self.specs: Optional[Dict[str, Dict[str, Any]]] = self.cfg.get("specs", None)

        self.tdt_roi_seg_path: Optional[str] = getattr(context, "tdt_roi_seg_path", None)
        
        self.roi_subset: List[str] = self.context.config["spect_preprocessing"]["roi_subset"]

        # Load label map from src/data/tdt_map.json (same approach as your other stages)
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
        self.ts_map_path: str = os.path.join(repo_root, "data", "tdt_map.json")
        if not os.path.exists(self.ts_map_path):
            raise FileNotFoundError(f"Class map json not found: {self.ts_map_path}")

        with open(self.ts_map_path, encoding="utf-8") as f:
            ts_map_json: Dict[str, Dict[str, str]] = json.loads(json_minify(f.read()))

        # Build name -> id map for TDT pipeline labels
        self.tdt_name2id: Dict[str, int] = {
            name: int(lab) for lab, name in ts_map_json["TDT_Pipeline"].items()
        }

        if "synthetic_lesion" not in self.tdt_name2id:
            raise ValueError(
                "tdt_map.json missing 'synthetic_lesion' in TDT_Pipeline. "
                "Add it (e.g. \"8\": \"synthetic_lesion\")."
            )
        self.synthetic_lesion_id: int = int(self.tdt_name2id["synthetic_lesion"])

        self.lesions_outdir = os.path.join(self.output_dir, f"{self.prefix}_outputs")

    @staticmethod
    def _xyz_to_zyx(arr_xyz: np.ndarray) -> np.ndarray:
        # NIfTI arrays are typically (X,Y,Z). We operate in (Z,Y,X).
        return np.transpose(arr_xyz, (2, 1, 0))

    @staticmethod
    def _zyx_to_xyz(arr_zyx: np.ndarray) -> np.ndarray:
        """Transpose back from (Z,Y,X) to (X,Y,Z) for NIfTI saving."""
        return np.transpose(arr_zyx, (2, 1, 0))

    @staticmethod
    def _spacing_zyx_from_nifti(nii: nib.Nifti1Image) -> np.ndarray:
        """Extract voxel spacing in (Z,Y,X) order from NIfTI header."""
        sx, sy, sz = nii.header.get_zooms()[:3]
        return np.array([sz, sy, sx], dtype=np.float64)

    @staticmethod
    def _phys_dist_mm(
        idx1_zyx: Sequence[int],
        idx2_zyx: Sequence[int],
        spacing_zyx: np.ndarray,
    ) -> float:
        """Compute physical distance in mm between two voxel indices."""
        d = (np.array(idx1_zyx, dtype=np.float64) - np.array(idx2_zyx, dtype=np.float64)) * spacing_zyx
        return float(np.sqrt(np.sum(d * d)))

    @staticmethod
    def _candidate_weights(
        mask_zyx: np.ndarray,
        cand_zyx: np.ndarray,
        spacing_zyx: np.ndarray,
        prob: str,
        sigma_mm: float | None = None,
        tom_map_zyx: np.ndarray | None = None,
    ) -> np.ndarray:
        """
Returns weight for each candidate center voxel based on the specified probability scheme.
        """
        prob = prob.lower()
        n = cand_zyx.shape[0]

        if prob == "uniform":
            return np.ones(n, dtype=np.float64)

        if prob == "gaussian":
            if sigma_mm is None:
                raise ValueError("sigma_mm required for prob='gaussian'")
            roi_pts = np.argwhere(mask_zyx)
            if roi_pts.size == 0:
                raise ValueError("ROI mask is empty; cannot compute gaussian centroid.")
            mu = roi_pts.mean(axis=0)  # (z,y,x) centroid in voxel coords

            dz = (cand_zyx[:, 0] - mu[0]) * spacing_zyx[0]
            dy = (cand_zyx[:, 1] - mu[1]) * spacing_zyx[1]
            dx = (cand_zyx[:, 2] - mu[2]) * spacing_zyx[2]
            r2 = dx * dx + dy * dy + dz * dz
            w = np.exp(-0.5 * r2 / (sigma_mm**2))
            return w.astype(np.float64)

        if prob == "tom":
            raise NotImplementedError(
                "prob='tom' selected, but TOM integration is TODO. "
                "Use prob='uniform' or prob='gaussian' for now."
            )

        if prob == "user_defined":
            return np.ones(n, dtype=np.float64)

        raise ValueError(f"Unknown prob choice: {prob}")

    @staticmethod
    def _place_lesion_centers(
        mask_zyx: np.ndarray,
        radii_mm: List[float],
        spacing_zyx: np.ndarray,
        prob: str = "uniform",
        sigma_mm: float | None = None,
        margin_mm: float = 1.0,
        seed: int = 0,
        max_attempts_per_lesion: int = 4000,
        tom_map_zyx: np.ndarray | None = None,
        user_centers_zyx: List[Tuple[int, int, int]] | None = None,
    ) -> Tuple[List[Tuple[int, int, int]], List[float], np.ndarray]:
        """
        Returns:
          centers_zyx: list of (z,y,x) ints
          radii_mm:    list of radii (same order)
          dist_mm:     distance-to-boundary map (mm) inside ROI
        """
        rng = np.random.default_rng(seed)
        dist_mm = distance_transform_edt(mask_zyx.astype(np.uint8), sampling=spacing_zyx)

        centers_zyx: List[Tuple[int, int, int]] = []
        placed_r: List[float] = []
        prob = prob.lower()

        # user defined: validate and skip sampling
        if prob == "user_defined":
            if user_centers_zyx is None:
                raise ValueError("prob='user_defined' but user_centers_zyx=None")
            if len(user_centers_zyx) != len(radii_mm):
                raise ValueError("user_centers_zyx length must match radii_mm length")

            for c, r in zip(user_centers_zyx, radii_mm):
                c = tuple(map(int, c))
                if not mask_zyx[c]:
                    raise ValueError(f"User center {c} not inside ROI")
                if dist_mm[c] < (r + margin_mm):
                    raise ValueError(f"User center {c} too close to ROI boundary for radius {r} mm")
                for cj, rj in zip(centers_zyx, placed_r):
                    if SyntheticLesionsStage._phys_dist_mm(c, cj, spacing_zyx) < (r + rj + margin_mm):
                        raise ValueError(f"User center {c} overlaps existing lesion at {cj}")
                centers_zyx.append(c)
                placed_r.append(r)

            return centers_zyx, placed_r, dist_mm

        # sample lesions sequentially (same as notebook)
        for i, r in enumerate(radii_mm, start=1):
            cand_mask = dist_mm >= (r + margin_mm)
            cand = np.argwhere(cand_mask)
            if cand.shape[0] == 0:
                raise RuntimeError(
                    f"No valid candidate centers for radius={r} mm. "
                    "Try smaller radius or smaller margin."
                )

            w = SyntheticLesionsStage._candidate_weights(
                mask_zyx, cand, spacing_zyx, prob, sigma_mm=sigma_mm, tom_map_zyx=tom_map_zyx
            )
            w = np.maximum(w, 0)
            w_sum = w.sum()
            if w_sum <= 0:
                raise RuntimeError("All candidate weights are zero. Check gaussian sigma / probability map.")
            p = w / w_sum

            placed = False
            for _ in range(max_attempts_per_lesion):
                k = rng.choice(cand.shape[0], p=p)
                c = tuple(map(int, cand[k]))

                ok = True
                for cj, rj in zip(centers_zyx, placed_r):
                    if SyntheticLesionsStage._phys_dist_mm(c, cj, spacing_zyx) < (r + rj + margin_mm):
                        ok = False
                        break

                if ok:
                    centers_zyx.append(c)
                    placed_r.append(r)
                    placed = True
                    break

            if not placed:
                raise RuntimeError(
                    f"Failed to place lesion {i}/{len(radii_mm)} (r={r}mm) after {max_attempts_per_lesion} attempts. "
                    "Try reducing radii, margin_mm, or switching prob='uniform'."
                )

        return centers_zyx, placed_r, dist_mm

    @staticmethod
    def _build_lesion_labelmap(
        mask_zyx: np.ndarray,
        centers_zyx: List[Tuple[int, int, int]],
        radii_mm: List[float],
        spacing_zyx: np.ndarray,
    ) -> np.ndarray:
        """
        Returns int16 labelmap (Z,Y,X):
          0 = background
          1..N = lesion id
        """
        Z, Y, X = mask_zyx.shape
        labels = np.zeros((Z, Y, X), dtype=np.int16)

        for lbl, (c, r) in enumerate(zip(centers_zyx, radii_mm), start=1):
            z0, y0, x0 = c

            rz = int(np.ceil(r / spacing_zyx[0]))
            ry = int(np.ceil(r / spacing_zyx[1]))
            rx = int(np.ceil(r / spacing_zyx[2]))

            zmin, zmax = max(0, z0 - rz), min(Z, z0 + rz + 1)
            ymin, ymax = max(0, y0 - ry), min(Y, y0 + ry + 1)
            xmin, xmax = max(0, x0 - rx), min(X, x0 + rx + 1)

            zz, yy, xx = np.ogrid[zmin:zmax, ymin:ymax, xmin:xmax]
            dz = (zz - z0) * spacing_zyx[0]
            dy = (yy - y0) * spacing_zyx[1]
            dx = (xx - x0) * spacing_zyx[2]

            sphere = (dx * dx + dy * dy + dz * dz) <= (r * r)
            sphere &= mask_zyx[zmin:zmax, ymin:ymax, xmin:xmax]

            labels[zmin:zmax, ymin:ymax, xmin:xmax][sphere] = lbl

        return labels

    @staticmethod
    def _write_json(path: str, payload: Dict[str, Any]) -> None:
        """Write JSON with indentation for readability (for QC metadata files)."""
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)

    def _assert_inputs_exist(self) -> None:
        """Check that required inputs exist before running. Raises informative errors if not."""
        if self.specs is None or len(self.specs) == 0:
            return
        if self.tdt_roi_seg_path is None or not os.path.exists(self.tdt_roi_seg_path):
            raise FileNotFoundError(f"Unified TDT ROI seg not found: {self.tdt_roi_seg_path}")

    def run(self) -> Any:
        """
        Run Synthetic Lesions Stage.

        Returns
        -------
        context : Context-like
            Updated context object with `tdt_roi_seg_path` set.
        """
        self._assert_inputs_exist() # early check to fail fast if inputs are missing

        # no-op if disabled
        if self.specs is None or len(self.specs) == 0:
            self.context.synthetic_lesions_outdir = None
            self.context.synthetic_lesions_results = {}
            return self.context

        os.makedirs(self.lesions_outdir, exist_ok=True)

        # IMPORTANT:
        # Preprocessing stage filters labels based on config["spect_preprocessing"]["roi_subset"].
        # We do NOT want lesion voxels to be dropped to 0, so ensure "synthetic_lesion" is included
        # AFTER segmentation/unify has already run (so TotalSegmentationStage doesn't see it).
        self.roi_subset = [str(r).strip() for r in self.roi_subset if str(r).strip()] # remove empty/whitespace entries
        if "synthetic_lesion" not in self.roi_subset:
            self.roi_subset.append("synthetic_lesion")
            self.context.config["spect_preprocessing"]["roi_subset"] = self.roi_subset

        # Load unified seg
        seg_nii = nib.load(self.tdt_roi_seg_path)
        seg_xyz = seg_nii.get_fdata(dtype=np.float32)  # load as float32 to avoid issues with large int16 lesion labels
        seg_zyx = self._xyz_to_zyx(seg_xyz)
        spacing_zyx = self._spacing_zyx_from_nifti(seg_nii)

        # Backup pre-lesion seg (for QC)
        backup_path = os.path.join(self.lesions_outdir, "tdt_roi_seg_pre_lesions.nii.gz")
        backup_img = nib.Nifti1Image(seg_xyz.astype(np.uint8), seg_nii.affine, seg_nii.header.copy())
        backup_img.set_data_dtype(np.uint8)
        nib.save(backup_img, backup_path)

        results: Dict[str, Any] = {}

        global_lesion_binary_zyx = np.zeros(seg_zyx.shape, dtype=np.uint8) # binary mask of all lesions across all ROIs (for QC)
        global_lesion_labels_zyx = np.zeros(seg_zyx.shape, dtype=np.int16) # labelmap of all lesions across all ROIs, with unique lesion ids (for QC)
        global_next_id = 1 # start lesion IDs at 1 in global labelmap (0 is background). This is separate from the "synthetic_lesion" label in the unified seg, which remains constant (e.g. 8) for all lesion voxels.

        for roi_name, spec in self.specs.items():
            if roi_name not in self.tdt_name2id:
                raise ValueError(f"ROI '{roi_name}' not found in TDT_Pipeline label map.")

            if roi_name == "synthetic_lesion":
                raise ValueError("Do not specify lesions inside ROI='synthetic_lesion'.")

            roi_id = int(self.tdt_name2id[roi_name]) # get integer label for this ROI from tdt_map.json
            organ_mask_zyx = (seg_zyx == roi_id) # binary mask of the current organ/ROI in (Z,Y,X) order
            if organ_mask_zyx.sum() == 0:
                raise ValueError(f"[{roi_name}] mask is empty in unified segmentation.")

            n_lesions = int(spec.get("n_lesions", 0))
            radii_mm = [float(r) for r in spec.get("radii_mm", [])]
            prob = str(spec.get("prob", "uniform"))
            sigma_mm = spec.get("sigma_mm", None)
            margin_mm = float(spec.get("margin_mm", 1.0))
            seed = int(spec.get("seed", 0))
            user_centers_raw = spec.get("user_centers_zyx", None)

            if n_lesions <= 0:
                raise ValueError(f"[{roi_name}] n_lesions must be > 0")
            if len(radii_mm) != n_lesions:
                raise ValueError(f"[{roi_name}] radii_mm length must match n_lesions")

            user_centers_zyx: Optional[List[Tuple[int, int, int]]] = None
            if prob.lower() == "user_defined":
                if user_centers_raw is None:
                    raise ValueError(f"[{roi_name}] prob='user_defined' requires user_centers_zyx in config.")
                user_centers_zyx = [tuple(map(int, c)) for c in user_centers_raw]

            centers_zyx, placed_radii, dist_mm = self._place_lesion_centers(
                mask_zyx=organ_mask_zyx,
                radii_mm=radii_mm,
                spacing_zyx=spacing_zyx,
                prob=prob,
                sigma_mm=float(sigma_mm) if sigma_mm is not None else None,
                margin_mm=margin_mm,
                seed=seed,
                max_attempts_per_lesion=4000, # hard coded
                tom_map_zyx=None,
                user_centers_zyx=user_centers_zyx,
            )

            roi_lesion_labels_zyx = self._build_lesion_labelmap(
                mask_zyx=organ_mask_zyx,
                centers_zyx=centers_zyx,
                radii_mm=placed_radii,
                spacing_zyx=spacing_zyx,
            )
            roi_lesion_binary_zyx = (roi_lesion_labels_zyx > 0).astype(np.uint8)
            roi_organ_minus_lesions_zyx = (organ_mask_zyx & (roi_lesion_binary_zyx == 0)).astype(np.uint8)

            # Update global lesion masks
            global_lesion_binary_zyx |= roi_lesion_binary_zyx
            for local_id in range(1, int(roi_lesion_labels_zyx.max()) + 1):
                global_lesion_labels_zyx[roi_lesion_labels_zyx == local_id] = global_next_id
                global_next_id += 1

            # Save per-ROI QC outputs
            roi_dir = os.path.join(self.lesions_outdir, roi_name)
            os.makedirs(roi_dir, exist_ok=True)

            roi_labels_xyz = self._zyx_to_xyz(roi_lesion_labels_zyx).astype(np.int16)
            roi_bin_xyz = self._zyx_to_xyz(roi_lesion_binary_zyx).astype(np.uint8)
            roi_minus_xyz = self._zyx_to_xyz(roi_organ_minus_lesions_zyx).astype(np.uint8)

            labels_path = os.path.join(roi_dir, f"{roi_name}_lesions_labels.nii.gz")
            bin_path = os.path.join(roi_dir, f"{roi_name}_lesions_binary.nii.gz")
            minus_path = os.path.join(roi_dir, f"{roi_name}_organ_minus_lesions.nii.gz")
            meta_path = os.path.join(roi_dir, f"{roi_name}_lesion_metadata.json")

            out_labels = nib.Nifti1Image(roi_labels_xyz, seg_nii.affine, seg_nii.header.copy())
            out_labels.set_data_dtype(np.int16)
            nib.save(out_labels, labels_path)

            out_bin = nib.Nifti1Image(roi_bin_xyz, seg_nii.affine, seg_nii.header.copy())
            out_bin.set_data_dtype(np.uint8)
            nib.save(out_bin, bin_path)

            out_minus = nib.Nifti1Image(roi_minus_xyz, seg_nii.affine, seg_nii.header.copy())
            out_minus.set_data_dtype(np.uint8)
            nib.save(out_minus, minus_path)

            meta = {
                "roi": roi_name,
                "roi_id": roi_id,
                "synthetic_lesion_id": self.synthetic_lesion_id,
                "prob": prob,
                "sigma_mm": float(sigma_mm) if sigma_mm is not None else None,
                "margin_mm": margin_mm,
                "seed": seed,
                "spacing_zyx_mm": spacing_zyx.tolist(),
                "centers_zyx": [list(c) for c in centers_zyx],
                "radii_mm": [float(r) for r in placed_radii],
                "dist_to_boundary_mm": [float(dist_mm[c]) for c in centers_zyx],
                "paths": {
                    "lesions_labels": labels_path,
                    "lesions_binary": bin_path,
                    "organ_minus_lesions": minus_path,
                },
            }
            self._write_json(meta_path, meta)
            results[roi_name] = meta

        # Save global lesion masks (QC)
        global_bin_xyz = self._zyx_to_xyz(global_lesion_binary_zyx).astype(np.uint8)
        global_lbl_xyz = self._zyx_to_xyz(global_lesion_labels_zyx).astype(np.int16)

        global_bin_path = os.path.join(self.lesions_outdir, "all_lesions_binary.nii.gz")
        global_lbl_path = os.path.join(self.lesions_outdir, "all_lesions_labels.nii.gz")

        out_gb = nib.Nifti1Image(global_bin_xyz, seg_nii.affine, seg_nii.header.copy())
        out_gb.set_data_dtype(np.uint8)
        nib.save(out_gb, global_bin_path)

        out_gl = nib.Nifti1Image(global_lbl_xyz, seg_nii.affine, seg_nii.header.copy())
        out_gl.set_data_dtype(np.int16)
        nib.save(out_gl, global_lbl_path)

        # ---------------------------------------------------------
        # Insert lesions into unified seg as synthetic_lesion label
        # ---------------------------------------------------------
        seg_zyx_mod = seg_zyx.copy()
        seg_zyx_mod[global_lesion_binary_zyx > 0] = np.uint8(self.synthetic_lesion_id)
        seg_xyz_mod = self._zyx_to_xyz(seg_zyx_mod).astype(np.uint8)

        # OVERWRITE original unified seg on disk
        out_seg = nib.Nifti1Image(seg_xyz_mod, seg_nii.affine, seg_nii.header.copy())
        out_seg.set_data_dtype(np.uint8)
        nib.save(out_seg, self.tdt_roi_seg_path)

        # Update context
        self.context.synthetic_lesions_outdir = self.lesions_outdir
        self.context.synthetic_lesions_results = results
        self.context.synthetic_lesions_backup_seg_path = backup_path
        self.context.synthetic_lesions_global_binary_path = global_bin_path
        self.context.synthetic_lesions_global_labels_path = global_lbl_path
        self.context.tdt_roi_seg_path = self.tdt_roi_seg_path  # same path, overwritten contents

        return self.context