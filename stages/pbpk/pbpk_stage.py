import os
import numpy as np
import pycno


class PbpkStage:
    def __init__(self, config, context):
        self.config = config
        self.context = context

        subdir_name = config["subdir_names"]["pbpk"]
        output_root = config["output_folder"]["title"]
        self.output_dir = os.path.join(output_root, subdir_name)
        os.makedirs(self.output_dir, exist_ok=True)

        self.prefix = config["pbpk"]["name"]

        self.vois_pbpk = config["pbpk"]["VOIs"]
        self.frame_start = config["pbpk"]["FrameStartTimes"]
        self.frame_stop = max(self.frame_start)

    # -----------------------------
    # helpers
    # -----------------------------
    def _remove_background(self, class_seg):
        if "Background" in class_seg:
            class_seg = dict(class_seg)
            del class_seg["Background"]
        return class_seg

    def _voxel_volume_ml(self, arr_px_spacing_cm):
        arr_px_spacing_cm = np.asarray(arr_px_spacing_cm, dtype=float)
        return float(np.prod(arr_px_spacing_cm))  # cm^3 == mL

    def _run_psma_model(self):
        time, tacs = pycno.run_model(
            model_name="PSMA",
            stop=self.frame_stop,
            observables=self.vois_pbpk,
        )
        return time, tacs

    def _roi_to_voi(self, roi_name):
        roi_to_voi = {
            "kidney": "Kidney",
            "body": "Rest",
            "liver": "Liver",
            "prostate": "Prostate",
        }
        return roi_to_voi.get(roi_name)

    def _save_tac_files(self, roi_name, time, tac_voi, tac_interp):
        roi_tag = roi_name.lower()

        tac_time_f = os.path.join(self.output_dir, f"{self.prefix}_{roi_tag}_TAC_time.bin")
        tac_vals_f = os.path.join(self.output_dir, f"{self.prefix}_{roi_tag}_TAC_values.bin")
        samp_time_f = os.path.join(self.output_dir, f"{self.prefix}_{roi_tag}_sample_times.bin")
        samp_vals_f = os.path.join(self.output_dir, f"{self.prefix}_{roi_tag}_sample_values.bin")

        np.asarray(time, dtype=np.float32).tofile(tac_time_f)
        np.asarray(tac_voi, dtype=np.float32).tofile(tac_vals_f)
        np.asarray(self.frame_start, dtype=np.float32).tofile(samp_time_f)
        np.asarray(tac_interp, dtype=np.float32).tofile(samp_vals_f)

        return {
            "tac_time": tac_time_f,
            "tac_values": tac_vals_f,
            "sample_time": samp_time_f,
            "sample_values": samp_vals_f,
        }
    def _compute_activity_for_roi(
        self,
        roi_name,
        label_value,
        mask_roi_body,
        roi_body_seg_arr,
        voxel_vol_ml,
        time,
        tacs,
        saved_tacs,
    ):
        voi_name = self._roi_to_voi(roi_name)

        if voi_name not in self.vois_pbpk:
            if "Rest" in self.vois_pbpk:
                voi_name = "Rest" # use 'Rest' VOI for unmapped ROIs
            else:
                raise ValueError(
                    f"ROI '{roi_name}' maps to VOI '{voi_name}', but that VOI is not in the PBPK model "
                    f"(supported: {sorted(self.vois_pbpk)}), and no 'Rest' VOI exists."
                )

        voi_index = self.vois_pbpk.index(voi_name)

        mask = mask_roi_body[label_value]
        n_vox = int(np.sum(mask))
        if n_vox == 0:
            raise AssertionError(f"Mask corresponding to {voi_name} is empty")

        tac_voi = tacs[0, :, voi_index]
        tac_interp = np.interp(self.frame_start, time, tac_voi)

        # [frame, z, y, x]
        activity_map_organ = np.zeros((len(self.frame_start), *roi_body_seg_arr.shape), dtype=np.float32)
        activity_map_organ[:, mask] = tac_interp[:, None] / (n_vox * voxel_vol_ml)

        # Save FIRST frame organ map for SIMIND usage
        organ_map_path = os.path.join(self.output_dir, f"{self.prefix}_{roi_name}_act_av.bin")
        activity_map_organ[0].astype(np.float32).tofile(organ_map_path)

        # Organ sum per frame [MBq]
        organ_sum = np.sum(activity_map_organ, axis=(1, 2, 3)) * voxel_vol_ml

        if voi_name not in saved_tacs:
            saved_tacs[roi_name] = self._save_tac_files(roi_name, time, tac_voi, tac_interp)

        return activity_map_organ, organ_sum, organ_map_path

    # -----------------------------
    # main
    # -----------------------------
    def run(self):
        self.context.require("roi_body_seg_arr", "mask_roi_body", "class_seg", "arr_px_spacing_cm")

        roi_body_seg_arr = self.context.roi_body_seg_arr
        mask_roi_body = self.context.mask_roi_body
        class_seg = self._remove_background(self.context.class_seg)
        voxel_vol_ml = self._voxel_volume_ml(self.context.arr_px_spacing_cm)

        print("arr_px_spacing_cm:", self.context.arr_px_spacing_cm)
        print("voxel_vol_ml:", voxel_vol_ml)
        
        time, tacs = self._run_psma_model()

        n_frames = len(self.frame_start)
        activity_map = np.zeros((n_frames, *roi_body_seg_arr.shape), dtype=np.float32)

        activity_organ_sum = {}
        organ_paths = []
        frame_paths = []
        saved_tacs = {}

        for roi_name, label_value in class_seg.items():
            out = self._compute_activity_for_roi(
                roi_name=roi_name,
                label_value=label_value,
                mask_roi_body=mask_roi_body,
                roi_body_seg_arr=roi_body_seg_arr,
                voxel_vol_ml=voxel_vol_ml,
                time=time,
                tacs=tacs,
                saved_tacs=saved_tacs,
            )
            if out is None:
                continue

            activity_map_organ, organ_sum, organ_map_path = out
            organ_paths.append(organ_map_path)
            activity_organ_sum[roi_name] = organ_sum

            mask = mask_roi_body[label_value]
            activity_map[:, mask] = activity_map_organ[:, mask]

        activity_map_sum = np.sum(activity_map, axis=(1, 2, 3)) * voxel_vol_ml

        # Save full maps per frame (optional, but kept)
        for i, frame in enumerate(activity_map):
            t = self.frame_start[i]
            p = os.path.join(self.output_dir, f"{self.prefix}_{t}_act_av.bin")
            frame.astype(np.float32).tofile(p)
            frame_paths.append(p)

        # Update context 
        self.context.activity_map_sum = activity_map_sum
        self.context.activity_organ_sum = activity_organ_sum
        self.context.activity_map_paths_by_organ = organ_paths
        self.context.activity_map_paths_by_frame = frame_paths

        # extras
        self.context.extras["pbpk_output_dir"] = self.output_dir
        self.context.extras["pbpk_saved_tacs"] = saved_tacs

        return self.context