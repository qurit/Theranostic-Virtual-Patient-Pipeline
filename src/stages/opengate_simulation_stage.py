"""
OpenGATE dosimetry stage for the Theranostic Digital Twin (TDT) pipeline.

This stage performs voxel-source Monte Carlo dose calculations directly on the
native phase-1 CT grid using OpenGATE.

High-level workflow
-------------------
1) Read the clean phase-1 CT handoff (`context.ct_nii_path`).
2) Read the clean phase-1 digital twin label map (`context.tdt_roi_seg_path`).
3) Build binary source masks for the ROIs currently tracked in
   `context.downstream_roi_subset`.
4) Run one OpenGATE voxel-source simulation per ROI.
5) Save raw dose NIfTI outputs in the phase-4 root directory.
6) Save work-dir metadata and optional OpenGATE MHD outputs.

Important note on thread handling
---------------------------------
In the user workflow this stage is based on, OpenGATE source history counts were
observed to scale with the thread count when `source.n` was set directly.
To keep the requested total histories under control, this stage converts the
requested total histories into a per-thread history count before assigning it to
OpenGATE sources.

Expected Context interface
--------------------------
Incoming `context` is expected to provide:
- context.subdir_paths["phase_4"] : str
- context.config["phase_4"]["opengate_simulation_stage"] : dict  # 
- context.ct_nii_path : str
- context.tdt_roi_seg_path : str
- context.downstream_roi_subset : list[str]
- context.config["phase_1"]["unification_stage"]["label_map_path"] : str

On success, this stage sets:
- context.dosimetry_output_dir : str
- context.dosimetry_stage_output_dir : str
- context.dosimetry_work_dir : str
- context.dosimetry_metadata_path : str
- context.dosimetry_mask_paths : dict[str, str]
- context.dosimetry_raw_dose_paths : dict[str, str]
- context.dosimetry_raw_uncertainty_paths : dict[str, str]
- context.dosimetry_sum_dose_path : str | None
- context.dosimetry_material_label_path : str | None

Maintainer / contact: pyazdi@bccrc.ca
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import json
import os

from json_minify import json_minify
import numpy as np
import SimpleITK as sitk
import opengate as gate


class DosimetryStage:
    """Run OpenGATE voxel-source dosimetry on the native CT grid."""

    def __init__(self, context: Any) -> None:
        context.require(
            "subdir_paths",
            "config",
            "ct_nii_path",
            "tdt_roi_seg_path",
            "downstream_roi_subset",
        )
        self.context = context

        self.phase_output_dir: str = context.subdir_paths["phase_4"]
        self.stage_cfg: Dict[str, Any] = context.config["phase_4"]["opengate_simulation_stage"]  
        self.stage_output_dir: str = os.path.join(
            self.phase_output_dir,
            self.stage_cfg.get("sub_dir_name", "opengate_simulation"), 
        )
        self.output_dir: str = self.stage_output_dir  
        self.work_dir: str = os.path.join(self.stage_output_dir, "work_dir")
        self.source_mask_dir: str = os.path.join(self.stage_output_dir, "source_masks")
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.stage_output_dir, exist_ok=True)
        os.makedirs(self.work_dir, exist_ok=True)
        os.makedirs(self.source_mask_dir, exist_ok=True)

        self.prefix: str = str(self.stage_cfg.get("file_prefix", "dosimetry"))
        self.metadata_path: str = os.path.join(self.work_dir, f"{self.prefix}_metadata.json")

        self.save_per_roi_dose_maps: bool = bool(self.stage_cfg.get("save_per_roi_dose_maps", True))
        self.save_summed_dose_map: bool = bool(self.stage_cfg.get("save_summed_dose_map", True))
        self.save_uncertainty_map: bool = bool(self.stage_cfg.get("save_uncertainty_map", True))
        self.save_material_label_image: bool = bool(self.stage_cfg.get("save_material_label_image", True))
        self.write_mhd_outputs: bool = bool(self.stage_cfg.get("write_mhd_outputs", False))

        gate_cfg = self.stage_cfg.get("gate", {})
        self.requested_total_histories: int = int(gate_cfg.get("total_histories", 100000))
        self.requested_num_threads: int = int(gate_cfg.get("num_threads", 1))
        self.start_new_process: bool = bool(gate_cfg.get("start_new_process", True))
        self.random_seed: Any = gate_cfg.get("random_seed", "auto")

        if self.requested_total_histories <= 0:
            raise ValueError("phase_4.opengate_simulation_stage.gate.total_histories must be > 0")
        if self.requested_num_threads <= 0:
            raise ValueError("phase_4.opengate_simulation_stage.gate.num_threads must be > 0")

        self.num_threads: int = min(self.requested_num_threads, self.requested_total_histories)
        if self.num_threads <= 0:
            raise ValueError("Effective OpenGATE thread count must be > 0")

        self.histories_per_thread_total: int = self.requested_total_histories // self.num_threads
        if self.histories_per_thread_total <= 0:
            raise ValueError(
                "Requested total histories are too small for the requested number of threads. "
                "Increase total_histories or reduce num_threads."
            )

        self.actual_total_histories: int = self.histories_per_thread_total * self.num_threads
        self.history_rounding_loss: int = self.requested_total_histories - self.actual_total_histories

        source_cfg = self.stage_cfg.get("source", {})
        self.source_components: List[Dict[str, Any]] = self._parse_source_components(source_cfg)
        self.component_histories_per_thread: List[int] = self._allocate_component_histories_per_thread(
            self.histories_per_thread_total,
            self.source_components,
        )

        physics_cfg = self.stage_cfg.get("physics", {})
        self.density_tolerance_gcm3: float = float(physics_cfg.get("density_tolerance_gcm3", 0.05))
        self.world_margin_scale: float = float(physics_cfg.get("world_margin_scale", 1.4))
        self.world_min_size_mm: float = float(physics_cfg.get("world_min_size_mm", 400.0))

        self.ct_nii_path: Path = Path(context.ct_nii_path)
        self.tdt_roi_seg_path: Path = Path(context.tdt_roi_seg_path)

        label_map_path = context.config["phase_1"]["unification_stage"]["label_map_path"]
        self.label_map_path: Path = Path(label_map_path)
        self.tdt_name2id: Dict[str, int] = self._load_tdt_label_map(self.label_map_path)

        roi_subset = getattr(context, "downstream_roi_subset", None)
        if isinstance(roi_subset, str):
            roi_subset = [roi_subset]
        if roi_subset is None:
            raise ValueError("context.downstream_roi_subset must be provided for opengate_simulation_stage")
        self.requested_roi_subset: List[str] = self._normalize_roi_subset(roi_subset)
        if not self.requested_roi_subset:
            raise ValueError("No ROI names were found in context.downstream_roi_subset")

    # -----------------------------
    # helpers
    # -----------------------------
    @staticmethod
    def _normalize_roi_subset(roi_subset: Sequence[str]) -> List[str]:
        """Return a de-duplicated ROI list while preserving order."""
        out: List[str] = []
        seen = set()
        for roi_name in roi_subset:
            name = str(roi_name).strip()
            if not name or name in seen:
                continue
            seen.add(name)
            out.append(name)
        return out

    @staticmethod
    def _load_tdt_label_map(label_map_path: Path) -> Dict[str, int]:
        """Load the TDT name->label_id mapping from the phase-1 label map JSON."""
        if not label_map_path.exists():
            raise FileNotFoundError(f"TDT label map not found: {label_map_path}")

        with open(label_map_path, encoding="utf-8") as f:
            map_json = json.loads(json_minify(f.read()))

        if "TDT_Pipeline" not in map_json:
            raise KeyError("Label map JSON is missing top-level key 'TDT_Pipeline'")

        return {name: int(label_id) for label_id, name in map_json["TDT_Pipeline"].items()}

    @staticmethod
    def _parse_source_components(source_cfg: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Validate and normalize the configured monoenergetic source components."""
        components_raw = source_cfg.get("components")
        if components_raw is None:
            particle = source_cfg.get("particle")
            energy_kev = source_cfg.get("energy_kev")
            if particle is None or energy_kev is None:
                raise ValueError(
                    "phase_4.opengate_simulation_stage.source must define 'components' or legacy "
                    "'particle' + 'energy_kev' fields"
                )
            components_raw = [{
                "particle": particle,
                "energy_kev": energy_kev,
                "relative_weight": source_cfg.get("relative_weight", 1.0),
            }]

        if not isinstance(components_raw, list) or len(components_raw) == 0:
            raise ValueError("phase_4.opengate_simulation_stage.source.components must be a non-empty list")

        components: List[Dict[str, Any]] = []
        for idx, item in enumerate(components_raw):
            if not isinstance(item, dict):
                raise TypeError("Each source component must be a dictionary")

            particle = str(item.get("particle", "")).strip()
            if not particle:
                raise ValueError(f"Source component {idx} is missing 'particle'")

            if "energy_kev" not in item:
                raise ValueError(f"Source component {idx} is missing 'energy_kev'")
            energy_kev = float(item["energy_kev"])
            if energy_kev <= 0:
                raise ValueError(f"Source component {idx} has non-positive energy_kev={energy_kev}")

            relative_weight = float(item.get("relative_weight", 1.0))
            if relative_weight <= 0:
                raise ValueError(
                    f"Source component {idx} has non-positive relative_weight={relative_weight}"
                )

            components.append(
                {
                    "particle": particle,
                    "energy_kev": energy_kev,
                    "relative_weight": relative_weight,
                }
            )
        return components

    @staticmethod
    def _allocate_component_histories_per_thread(
        histories_per_thread_total: int,
        components: Sequence[Dict[str, Any]],
    ) -> List[int]:
        """Split the per-thread history budget across source components by weight."""
        weights = np.asarray([float(c["relative_weight"]) for c in components], dtype=np.float64)
        weights = weights / np.sum(weights)

        raw_counts = weights * int(histories_per_thread_total)
        counts = np.floor(raw_counts).astype(int)
        remainder = int(histories_per_thread_total) - int(np.sum(counts))

        if remainder > 0:
            fractional = raw_counts - counts
            order = np.argsort(-fractional)
            for idx in order[:remainder]:
                counts[idx] += 1

        return counts.astype(int).tolist()

    @staticmethod
    def _images_match_geometry(reference_img: sitk.Image, other_img: sitk.Image, atol: float = 1e-6) -> bool:
        """Return True if two SimpleITK images share size, spacing, origin, and direction."""
        if reference_img.GetSize() != other_img.GetSize():
            return False
        if not np.allclose(reference_img.GetSpacing(), other_img.GetSpacing(), atol=atol):
            return False
        if not np.allclose(reference_img.GetOrigin(), other_img.GetOrigin(), atol=atol):
            return False
        if not np.allclose(reference_img.GetDirection(), other_img.GetDirection(), atol=atol):
            return False
        return True

    @staticmethod
    def _save_array_like_reference(reference_img: sitk.Image, arr_zyx: np.ndarray, out_path: Path) -> str:
        """Write a numpy array to NIfTI while copying geometry from a reference image."""
        out_path = Path(out_path)
        out_img = sitk.GetImageFromArray(np.asarray(arr_zyx))
        out_img.CopyInformation(reference_img)
        sitk.WriteImage(out_img, str(out_path), imageIO="NiftiImageIO")
        return str(out_path)

    @staticmethod
    def _cleanup_mhd_sidecars(mhd_path: Path) -> None:
        """Delete an .mhd file and common paired raw files if they exist."""
        mhd_path = Path(mhd_path)
        for path in [mhd_path, mhd_path.with_suffix(".raw"), mhd_path.with_suffix(".zraw")]:
            try:
                if path.exists():
                    path.unlink()
            except Exception:
                pass

    @staticmethod
    def _find_gate_hu_material_tables() -> Tuple[Path, Path]:
        """Locate OpenGATE's Schneider HU conversion tables from the installed package."""
        gate_pkg_dir = Path(gate.__file__).resolve().parent
        candidate_roots = [
            gate_pkg_dir,
            gate_pkg_dir.parent,
            *list(gate_pkg_dir.parents[:4]),
        ]
        relative_pairs = [
            (
                Path("tests/data/Schneider2000MaterialsTable.txt"),
                Path("tests/data/Schneider2000DensitiesTable.txt"),
            ),
            (
                Path("opengate/tests/data/Schneider2000MaterialsTable.txt"),
                Path("opengate/tests/data/Schneider2000DensitiesTable.txt"),
            ),
        ]

        for root in candidate_roots:
            for mat_rel, den_rel in relative_pairs:
                mat_path = root / mat_rel
                den_path = root / den_rel
                if mat_path.exists() and den_path.exists():
                    return mat_path, den_path

        material_hits: List[Path] = []
        density_hits: List[Path] = []
        for root in candidate_roots[:2]:
            try:
                material_hits.extend(root.rglob("Schneider2000MaterialsTable.txt"))
                density_hits.extend(root.rglob("Schneider2000DensitiesTable.txt"))
            except Exception:
                pass

        if material_hits and density_hits:
            return material_hits[0], density_hits[0]

        raise FileNotFoundError(
            "Could not automatically locate OpenGATE's HU conversion tables. "
            "Expected Schneider2000MaterialsTable.txt and Schneider2000DensitiesTable.txt "
            "to be present inside the installed OpenGATE package."
        )

    def _get_gate_hu_voxel_materials(self, sim: Any) -> Tuple[Any, Any, Path, Path]:
        """Build the OpenGATE voxel-material mapping from installed HU tables."""
        gcm3 = gate.g4_units.g_cm3
        mat_table_path, den_table_path = self._find_gate_hu_material_tables()
        voxel_materials, generated_materials = gate.geometry.materials.HounsfieldUnit_to_material(
            sim,
            float(self.density_tolerance_gcm3) * gcm3,
            str(mat_table_path),
            str(den_table_path),
        )
        return voxel_materials, generated_materials, mat_table_path, den_table_path

    def _build_source_masks(
        self,
        ct_img: sitk.Image,
        seg_img: sitk.Image,
    ) -> Tuple[Dict[str, str], Dict[str, int], List[str]]:
        """Create binary NIfTI source masks for all requested ROIs present in the digital twin."""
        if not self._images_match_geometry(ct_img, seg_img):
            raise ValueError(
                "CT image and digital twin segmentation do not share the same geometry. "
                "Dosimetry stage requires native-space alignment."
            )

        seg_arr = sitk.GetArrayFromImage(seg_img).astype(np.int32)
        mask_paths: Dict[str, str] = {}
        roi_voxel_counts: Dict[str, int] = {}
        simulated_roi_names: List[str] = []

        for roi_name in self.requested_roi_subset:
            label_id = self.tdt_name2id.get(roi_name)
            if label_id is None:
                raise ValueError(f"Requested ROI '{roi_name}' not present in TDT label map")

            mask_arr = (seg_arr == int(label_id)).astype(np.uint8)
            voxel_count = int(np.count_nonzero(mask_arr))
            if voxel_count == 0:
                continue

            out_path = Path(self.source_mask_dir) / f"{self.prefix}_{roi_name}_source_mask.nii.gz"
            self._save_array_like_reference(ct_img, mask_arr, out_path)

            mask_paths[roi_name] = str(out_path)
            roi_voxel_counts[roi_name] = voxel_count
            simulated_roi_names.append(roi_name)

        if not simulated_roi_names:
            raise ValueError(
                "None of the requested ROIs had voxels in the digital twin segmentation. "
                "Check context.downstream_roi_subset and the unified TDT label map."
            )

        return mask_paths, roi_voxel_counts, simulated_roi_names

    def _get_roi_output_paths(self, roi_name: str) -> Dict[str, str]:
        """Return the final phase-4 NIfTI outputs for a single ROI."""
        return {
            "dose": os.path.join(self.output_dir, f"{self.prefix}_dose_raw_{roi_name}.nii.gz"),
            "uncertainty": os.path.join(
                self.output_dir,
                f"{self.prefix}_dose_uncertainty_raw_{roi_name}.nii.gz",
            ),
        }

    def _run_single_roi(
        self,
        roi_name: str,
        mask_path: str,
        save_material_labels_this_run: bool,
    ) -> Dict[str, Any]:
        """Run one OpenGATE voxel-source simulation for a single ROI mask."""
        roi_work_dir = Path(self.work_dir) / roi_name
        roi_work_dir.mkdir(exist_ok=True, parents=True)

        ct_img = sitk.ReadImage(str(self.ct_nii_path))
        size_xyz = ct_img.GetSize()
        spacing_xyz = ct_img.GetSpacing()
        physical_size_mm = np.asarray(size_xyz, dtype=np.float64) * np.asarray(spacing_xyz, dtype=np.float64)

        mm = gate.g4_units.mm
        keV = gate.g4_units.keV

        sim = gate.Simulation()
        sim.g4_verbose = False
        sim.visu = False
        sim.number_of_threads = int(self.num_threads)
        sim.random_seed = self.random_seed
        sim.output_dir = roi_work_dir

        world_size_mm = np.maximum(physical_size_mm * float(self.world_margin_scale), float(self.world_min_size_mm))
        sim.world.size = [world_size_mm[0] * mm, world_size_mm[1] * mm, world_size_mm[2] * mm]
        sim.world.material = "G4_AIR"

        patient = sim.add_volume("Image", "patient")
        patient.image = str(self.ct_nii_path)
        patient.material = "G4_AIR"

        voxel_materials_auto, generated_materials, used_mat_table, used_den_table = self._get_gate_hu_voxel_materials(sim)
        patient.voxel_materials = voxel_materials_auto

        label_image_mhd_path: Optional[Path] = None
        if save_material_labels_this_run and self.save_material_label_image:
            patient.dump_label_image = "patient_material_labels.mhd"
            label_image_mhd_path = roi_work_dir / "patient_material_labels.mhd"

        translation = gate.image.get_translation_between_images_center(
            str(self.ct_nii_path),
            str(mask_path),
        )

        active_components: List[Dict[str, Any]] = []
        for idx, (component, histories_per_thread) in enumerate(
            zip(self.source_components, self.component_histories_per_thread)
        ):
            if int(histories_per_thread) <= 0:
                continue

            src = sim.add_source("VoxelSource", f"src_{idx}")
            src.particle = str(component["particle"])
            src.n = int(histories_per_thread)
            src.image = str(mask_path)
            src.attached_to = patient
            src.direction.type = "iso"
            src.energy.type = "mono"
            src.energy.mono = float(component["energy_kev"]) * keV
            src.position.translation = translation

            active_components.append(
                {
                    "particle": str(component["particle"]),
                    "energy_kev": float(component["energy_kev"]),
                    "relative_weight": float(component["relative_weight"]),
                    "histories_per_thread": int(histories_per_thread),
                    "actual_total_histories": int(histories_per_thread) * int(self.num_threads),
                }
            )

        if not active_components:
            raise ValueError(
                f"ROI '{roi_name}' has no active source components after history allocation. "
                "Increase total_histories or adjust component weights."
            )

        stats = sim.add_actor("SimulationStatisticsActor", "Stats")

        dose = sim.add_actor("DoseActor", "dose")
        dose.attached_to = patient
        dose.size = list(size_xyz)
        dose.spacing = [float(s) * mm for s in spacing_xyz]
        dose.dose.active = True
        dose.dose.output_filename = "dose_map.mhd"
        dose.dose_uncertainty.active = bool(self.save_uncertainty_map)
        if self.save_uncertainty_map:
            dose.dose_uncertainty.output_filename = "dose_uncertainty.mhd"

        sim.run(start_new_process=self.start_new_process)
        
        if label_image_mhd_path is not None and not label_image_mhd_path.exists():
            label_image_mhd_path = None

        dose_mhd_path = roi_work_dir / "dose_map.mhd"
        if not dose_mhd_path.exists():
            raise FileNotFoundError(f"OpenGATE dose output not found: {dose_mhd_path}")

        dose_img = sitk.ReadImage(str(dose_mhd_path))
        dose_arr = sitk.GetArrayFromImage(dose_img).astype(np.float32)

        unc_arr: Optional[np.ndarray] = None
        dose_uncertainty_mhd_path: Optional[Path] = None
        if self.save_uncertainty_map:
            dose_uncertainty_mhd_path = roi_work_dir / "dose_uncertainty.mhd"
            if not dose_uncertainty_mhd_path.exists():
                raise FileNotFoundError(
                    f"OpenGATE uncertainty output not found: {dose_uncertainty_mhd_path}"
                )
            unc_img = sitk.ReadImage(str(dose_uncertainty_mhd_path))
            unc_arr = sitk.GetArrayFromImage(unc_img).astype(np.float32)

        stats_events = None
        try:
            stats_events = int(stats.counts.events)
        except Exception:
            stats_events = None

        return {
            "roi_name": roi_name,
            "roi_work_dir": str(roi_work_dir),
            "dose_mhd_path": str(dose_mhd_path),
            "dose_uncertainty_mhd_path": None if dose_uncertainty_mhd_path is None else str(dose_uncertainty_mhd_path),
            "label_image_mhd_path": None if label_image_mhd_path is None else str(label_image_mhd_path),
            "dose_arr": dose_arr,
            "unc_arr": unc_arr,
            "size_xyz": tuple(int(v) for v in size_xyz),
            "spacing_xyz": tuple(float(v) for v in spacing_xyz),
            "stats_str": str(stats),
            "stats_events": stats_events,
            "hu_materials_table_path_used": str(used_mat_table),
            "hu_densities_table_path_used": str(used_den_table),
            "generated_material_count": len(generated_materials),
            "source_components": active_components,
        }

    def _save_stage_metadata(
        self,
        simulated_roi_names: Sequence[str],
        roi_voxel_counts: Dict[str, int],
        roi_run_metadata: Dict[str, Any],
        mask_paths: Dict[str, str],
        dose_paths: Dict[str, str],
        uncertainty_paths: Dict[str, str],
        sum_dose_path: Optional[str],
        material_label_path: Optional[str],
    ) -> None:
        """Save dosimetry-stage metadata for debugging and provenance."""
        metadata: Dict[str, Any] = {
            "stage": "dosimetry_stage",
            "ct_nii_path": str(self.ct_nii_path),
            "tdt_roi_seg_path": str(self.tdt_roi_seg_path),
            "label_map_path": str(self.label_map_path),
            "requested_roi_subset": list(self.requested_roi_subset),
            "simulated_roi_names": list(simulated_roi_names),
            "roi_voxel_counts": {k: int(v) for k, v in roi_voxel_counts.items()},
            "gate": {
                "requested_total_histories": int(self.requested_total_histories),
                "requested_num_threads": int(self.requested_num_threads),
                "effective_num_threads": int(self.num_threads),
                "histories_per_thread_total": int(self.histories_per_thread_total),
                "actual_total_histories": int(self.actual_total_histories),
                "history_rounding_loss": int(self.history_rounding_loss),
                "random_seed": self.random_seed,
                "start_new_process": bool(self.start_new_process),
            },
            "source_components_config": self.source_components,
            "source_component_histories_per_thread": self.component_histories_per_thread,
            "physics": {
                "density_tolerance_gcm3": float(self.density_tolerance_gcm3),
                "world_margin_scale": float(self.world_margin_scale),
                "world_min_size_mm": float(self.world_min_size_mm),
            },
            "save_flags": {
                "save_per_roi_dose_maps": bool(self.save_per_roi_dose_maps),
                "save_summed_dose_map": bool(self.save_summed_dose_map),
                "save_uncertainty_map": bool(self.save_uncertainty_map),
                "save_material_label_image": bool(self.save_material_label_image),
                "write_mhd_outputs": bool(self.write_mhd_outputs),
            },
            "mask_paths": mask_paths,
            "dose_paths": dose_paths,
            "uncertainty_paths": uncertainty_paths,
            "sum_dose_path": sum_dose_path,
            "material_label_path": material_label_path,
            "roi_run_metadata": roi_run_metadata,
        }

        with open(self.metadata_path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=4)

    # -----------------------------
    # public API
    # -----------------------------
    def run(self) -> Any:
        """Execute the OpenGATE dosimetry stage and update the pipeline context."""
        ct_img = sitk.ReadImage(str(self.ct_nii_path))
        seg_img = sitk.ReadImage(str(self.tdt_roi_seg_path))

        mask_paths, roi_voxel_counts, simulated_roi_names = self._build_source_masks(ct_img, seg_img)

        dose_paths: Dict[str, str] = {}
        uncertainty_paths: Dict[str, str] = {}
        roi_run_metadata: Dict[str, Any] = {}
        sum_dose_arr: Optional[np.ndarray] = None
        material_label_path: Optional[str] = None

        for idx, roi_name in enumerate(simulated_roi_names):
            run_res = self._run_single_roi(
                roi_name=roi_name,
                mask_path=mask_paths[roi_name],
                save_material_labels_this_run=(idx == 0),
            )

            dose_arr = np.asarray(run_res["dose_arr"], dtype=np.float64)
            if sum_dose_arr is None:
                sum_dose_arr = dose_arr.copy()
            else:
                sum_dose_arr += dose_arr

            if self.save_per_roi_dose_maps:
                roi_paths = self._get_roi_output_paths(roi_name)
                dose_paths[roi_name] = self._save_array_like_reference(
                    ct_img,
                    np.asarray(run_res["dose_arr"], dtype=np.float32),
                    Path(roi_paths["dose"]),
                )

            if self.save_uncertainty_map and run_res["unc_arr"] is not None:
                roi_paths = self._get_roi_output_paths(roi_name)
                uncertainty_paths[roi_name] = self._save_array_like_reference(
                    ct_img,
                    np.asarray(run_res["unc_arr"], dtype=np.float32),
                    Path(roi_paths["uncertainty"]),
                )   
                
            if (
                material_label_path is None
                and self.save_material_label_image
                and run_res["label_image_mhd_path"] is not None
            ):
                label_path = run_res["label_image_mhd_path"]
                if label_path is not None and os.path.exists(label_path):
                    label_img = sitk.ReadImage(str(label_path))
                    label_arr = sitk.GetArrayFromImage(label_img).astype(np.int16)
                    out_path = Path(self.output_dir) / f"{self.prefix}_material_labels.nii.gz"
                    material_label_path = self._save_array_like_reference(ct_img, label_arr, out_path)
                else:
                    material_label_path = None

            if not self.write_mhd_outputs:
                self._cleanup_mhd_sidecars(Path(run_res["dose_mhd_path"]))
                if run_res["dose_uncertainty_mhd_path"] is not None:
                    self._cleanup_mhd_sidecars(Path(run_res["dose_uncertainty_mhd_path"]))
                if run_res["label_image_mhd_path"] is not None:
                    self._cleanup_mhd_sidecars(Path(run_res["label_image_mhd_path"]))

            roi_run_metadata[roi_name] = {
                "roi_work_dir": run_res["roi_work_dir"],
                "stats_str": run_res["stats_str"],
                "stats_events": run_res["stats_events"],
                "size_xyz": list(run_res["size_xyz"]),
                "spacing_xyz": list(run_res["spacing_xyz"]),
                "hu_materials_table_path_used": run_res["hu_materials_table_path_used"],
                "hu_densities_table_path_used": run_res["hu_densities_table_path_used"],
                "generated_material_count": int(run_res["generated_material_count"]),
                "source_components": run_res["source_components"],
            }

        sum_dose_path: Optional[str] = None
        if self.save_summed_dose_map and sum_dose_arr is not None:
            out_path = Path(self.output_dir) / f"{self.prefix}_dose_raw_sum.nii.gz"
            sum_dose_path = self._save_array_like_reference(
                ct_img,
                np.asarray(sum_dose_arr, dtype=np.float32),
                out_path,
            )

        self._save_stage_metadata(
            simulated_roi_names=simulated_roi_names,
            roi_voxel_counts=roi_voxel_counts,
            roi_run_metadata=roi_run_metadata,
            mask_paths=mask_paths,
            dose_paths=dose_paths,
            uncertainty_paths=uncertainty_paths,
            sum_dose_path=sum_dose_path,
            material_label_path=material_label_path,
        )

        self.context.dosimetry_output_dir = self.output_dir
        self.context.dosimetry_stage_output_dir = self.stage_output_dir
        self.context.dosimetry_work_dir = self.work_dir
        self.context.dosimetry_metadata_path = self.metadata_path
        self.context.dosimetry_mask_paths = mask_paths
        self.context.dosimetry_raw_dose_paths = dose_paths
        self.context.dosimetry_raw_uncertainty_paths = uncertainty_paths
        self.context.dosimetry_sum_dose_path = sum_dose_path
        self.context.dosimetry_material_label_path = material_label_path
        self.context.extras["opengate_simulation_stage"] = {
            "simulated_roi_names": list(simulated_roi_names),
            "requested_total_histories": int(self.requested_total_histories),
            "effective_num_threads": int(self.num_threads),
            "actual_total_histories": int(self.actual_total_histories),
            "history_rounding_loss": int(self.history_rounding_loss),
            "metadata_path": self.metadata_path,
        }
        return self.context
