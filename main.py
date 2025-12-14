import os
import json
from json_minify import json_minify

from context import Context

from modules.spect_pre_process.segmentation_stage import TotalSegmentationStage
from modules.spect_pre_process.preprocessing_stage import SimindPreprocessStage

from modules.pbpk.pbpk_stage import PbpkStage
from modules.spect_simulation.simind_stage import SimindSimulationStage
from modules.spect_simulation.reconstruction_stage import SpectReconstructionStage


class TdtPipeline:
    def __init__(self, config_path):
        self.config_path = config_path
        self.current_dir_path = os.path.dirname(__file__)
        self.setup_config(config_path)

    def setup_config(self, config_path):
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        with open(config_path, encoding="utf-8") as f:
            self.config = json.loads(json_minify(f.read()))

        output_folder_title = self.config["output_folder"]["title"]
        output_folder_path = os.path.join(self.current_dir_path, output_folder_title)
        os.makedirs(output_folder_path, exist_ok=True)

        subdir_dict = self.config["subdir_names"]
        for _, name in subdir_dict.items():
            subdir_path = os.path.join(output_folder_path, name)
            os.makedirs(subdir_path, exist_ok=True)

    def run(self):
        context = Context()

        # -----------------------------
        # Stage 1: TotalSegmentator
        # -----------------------------
        context = TotalSegmentationStage(self.config, context).run()

        # -----------------------------
        # Stage 2: Preprocess for SIMIND
        # -----------------------------
        context = SimindPreprocessStage(self.config, context).run()

        # -----------------------------
        # Stage 3: PBPK (function for now)
        # -----------------------------
        context = PbpkStage(self.config, context).run()

        # -----------------------------
        # Stage 4: SIMIND (function for now)
        # -----------------------------
        context = SimindSimulationStage(self.current_dir_path, self.config, context).run()

        # -----------------------------
        # Stage 5: Recon (function for now)
        # -----------------------------
        context = SpectReconstructionStage(self.config, context).run()

        # -----------------------------
        # Stage 6: Post-process (optional)
        # -----------------------------
        # context = resample_spect_to_atn_grid(self.config, context)

        print(context)
        return context


if __name__ == "__main__":
    config_path = "setup/config.json"
    pipeline1 = TdtPipeline(config_path)
    pipeline1.run()
