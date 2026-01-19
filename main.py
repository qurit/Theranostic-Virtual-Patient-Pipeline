import os
import json
from json_minify import json_minify
import logging
import time

from context import Context

from stages.spect_pre_process.segmentation_stage import TotalSegmentationStage
from stages.spect_pre_process.unify_ts_outputs import TdtRoiUnifyStage
from stages.spect_pre_process.preprocessing_stage import SimindPreprocessStage

from stages.pbpk.pbpk_stage import PbpkStage
from stages.spect_simulation.simind_stage import SimindSimulationStage
from stages.spect_simulation.reconstruction_stage import SpectReconstructionStage


class TdtPipeline:
    def __init__(self, config_path):
        self.config_path = config_path
        self.current_dir_path = os.path.dirname(__file__)
        self.setup_config(config_path)
        self.mode = self.config["mode"]["mode"]
    
    def _log_config_to_file(self):
        '''
        For now only saved Config, will add on more logging later on
        '''
        log_path = os.path.join(self.output_folder_path, "logging_file.log")

        logger = logging.getLogger("TDT_CONFIG_LOGGER")
        logger.setLevel(logging.INFO)
        logger.propagate = False  # prevent duplicate output if root logger is set elsewhere
        
        
        # Avoid adding multiple handlers if pipeline is constructed multiple times
        if not any(isinstance(h, logging.FileHandler) and getattr(h, "baseFilename", "") == log_path for h in logger.handlers):
            fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
            fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
            logger.addHandler(fh)
            
        logger.info("Log started") 
            
        config = json.dumps(self.config, indent=2, sort_keys=True)
        logger.info("CONFIG SOURCE: %s", os.path.abspath(self.config_path))
        logger.info("CONFIG CONTENTS:\n%s", config)

    def setup_config(self, config_path):
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        with open(config_path, encoding="utf-8") as f:
            self.config = json.loads(json_minify(f.read()))

        output_folder_title = self.config["output_folder"]["title"]
        output_folder_path = os.path.join(self.current_dir_path, output_folder_title)
        os.makedirs(output_folder_path, exist_ok=True)
        
        self.output_folder_path = output_folder_path

        subdir_dict = self.config["subdir_names"]
        for _, name in subdir_dict.items():
            subdir_path = os.path.join(output_folder_path, name)
            os.makedirs(subdir_path, exist_ok=True)
        
        self._log_config_to_file()

    def run(self):
        logger = logging.getLogger("TDT_CONFIG_LOGGER")
        t_pipeline = time.perf_counter()                 

        logger.info("Pipeline start | mode=%s", self.mode)  
        
        print("Loading TDT Pipeline configuration and setting up context manager...")
        context = Context()
        print("Starting TDT Pipeline in mode:", self.mode) # TODO: mode will be used later, determine what gets saved
        
        # -----------------------------
        # Stage 1: TotalSegmentator
        # -----------------------------
        
        logger.info("Stage start: TotalSegmentator")  
        t_stage = time.perf_counter()                 
        
        print("Running TotalSegmentator Stage...")
        context = TotalSegmentationStage(self.config, context).run()
        print("TotalSegmentator Stage completed.")
        
        logger.info("Stage end: TotalSegmentator | elapsed=%.2fs", time.perf_counter() - t_stage)  
        
        # -----------------------------
        # Stage 2: Unification of TS outputs to TDT ROIs
        # -----------------------------
        logger.info("Stage start: TDT ROI Unification")
        t_stage = time.perf_counter()
        
        print("Running TDT ROI Unification Stage...")
        context = TdtRoiUnifyStage(self.config, context).run()
        print("TDT ROI Unification Stage completed.")
        
        logger.info("Stage end: TDT ROI Unification | elapsed=%.2fs", time.perf_counter() - t_stage)

        # -----------------------------
        # Stage 3: Preprocess for SIMIND
        # -----------------------------
        logger.info("Stage start: SIMIND Preprocessing") 
        t_stage = time.perf_counter()               
               
        print("Running SIMIND Preprocessing Stage...")
        context = SimindPreprocessStage(self.config, context).run()
        print("SIMIND Preprocessing Stage completed.")
        
        logger.info("Stage end: SIMIND Preprocessing | elapsed=%.2fs", time.perf_counter() - t_stage)  
        
        # -----------------------------
        # Stage 4: PBPK 
        # -----------------------------
        
        logger.info("Stage start: PBPK")  
        t_stage = time.perf_counter()     
        
        print("Running PBPK Stage...")
        context = PbpkStage(self.config, context).run()
        print("PBPK Stage completed.")
        
        logger.info("Stage end: PBPK | elapsed=%.2fs", time.perf_counter() - t_stage)  

        # -----------------------------
        # Stage 5: SIMIND 
        # -----------------------------
        logger.info("Stage start: SIMIND Simulation")  
        t_stage = time.perf_counter()                 
        
        print("Running SIMIND Simulation Stage...")
        context = SimindSimulationStage(self.current_dir_path, self.config, context).run()
        print("SIMIND Simulation Stage completed.")
        
        logger.info("Stage end: SIMIND Simulation | elapsed=%.2fs", time.perf_counter() - t_stage)  
        
        # -----------------------------
        # Stage 6: Recon 
        # -----------------------------
        
        logger.info("Stage start: SPECT Reconstruction")  
        t_stage = time.perf_counter()                      
        
        print("Running SPECT Reconstruction Stage...")
        context = SpectReconstructionStage(self.config, context).run()
        print("SPECT Reconstruction Stage completed.")

        logger.info("Stage end: SPECT Reconstruction | elapsed=%.2fs", time.perf_counter() - t_stage)  

        logger.info("Pipeline end | total_elapsed=%.2fs", time.perf_counter() - t_pipeline)  

        print("TDT Pipeline completed successfully.")
        return context


if __name__ == "__main__":
    config_path = "inputs/config.json"
    pipeline1 = TdtPipeline(config_path)
    pipeline1.run()
