import os
import json
from json_minify import json_minify

from context import Context

from stages.spect_pre_process.segmentation_stage import TotalSegmentationStage
from stages.spect_pre_process.preprocessing_stage import SimindPreprocessStage

from stages.pbpk.pbpk_stage import PbpkStage
from stages.spect_simulation.simind_stage import SimindSimulationStage
from stages.spect_simulation.reconstruction_stage import SpectReconstructionStage


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
        print("Loading TDT Pipeline configuration and setting up context manager...")
        context = Context()
        print("Starting TDT Pipeline...")
        # -----------------------------
        # Stage 1: TotalSegmentator
        # -----------------------------
        print("Running TotalSegmentator Stage...")
        context = TotalSegmentationStage(self.config, context).run()
        print("TotalSegmentator Stage completed.")
        # -----------------------------
        # Stage 2: Preprocess for SIMIND
        # -----------------------------
        print("Running SIMIND Preprocessing Stage...")
        context = SimindPreprocessStage(self.config, context).run()
        print("SIMIND Preprocessing Stage completed.")

        # -----------------------------
        # Stage 3: PBPK 
        # -----------------------------
        print("Running PBPK Stage...")
        context = PbpkStage(self.config, context).run()
        print("PBPK Stage completed.")

        # -----------------------------
        # Stage 4: SIMIND 
        # -----------------------------
        print("Running SIMIND Simulation Stage...")
        context = SimindSimulationStage(self.current_dir_path, self.config, context).run()
        print("SIMIND Simulation Stage completed.")
        
        # -----------------------------
        # Stage 5: Recon 
        # -----------------------------
        print("Running SPECT Reconstruction Stage...")
        context = SpectReconstructionStage(self.config, context).run()
        print("SPECT Reconstruction Stage completed.")

        print("TDT Pipeline completed successfully.")
        
        print("printing context for debugging...")
        print(context)
        return context


if __name__ == "__main__":
    config_path = "inputs/config.json"
    pipeline1 = TdtPipeline(config_path)
    pipeline1.run()
