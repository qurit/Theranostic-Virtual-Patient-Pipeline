import os
import json
from json_minify import json_minify
import logging
import time
import shutil

from context import Context

from stages.spect_pre_process.segmentation_stage import TotalSegmentationStage
from stages.spect_pre_process.unify_ts_outputs import TdtRoiUnifyStage
from stages.spect_pre_process.preprocessing_stage import SimindPreprocessStage

from stages.pbpk.pbpk_stage import PbpkStage
from stages.spect_simulation.simind_stage import SimindSimulationStage
from stages.spect_simulation.reconstruction_stage import SpectReconstructionStage


class TdtPipeline:
    def __init__(self, config_path,dicom_path,ct_indx):
        self.config_path = config_path
        self.dicom_path = dicom_path 
        self.ct_indx = ct_indx 
        self.current_dir_path = os.path.dirname(__file__)
        self.setup_config(config_path)
        self.mode = self.config["mode"]["mode"]
    
    def _copy_dicom_into_output_root(self):
        """
        Copy the *entire* DICOM folder (self.dicom_path) into the CT output root
        (same level as the log file), preserving the folder name.
        """
        logger = logging.getLogger(f"TDT_CONFIG_LOGGER_CT_{self.ct_indx}")

        src = self.dicom_path
        if not os.path.isdir(src):
            raise FileNotFoundError(f"DICOM path is not a directory: {src}")

        dst = os.path.join(self.output_folder_path, os.path.basename(os.path.normpath(src)))

        # If already copied, skip
        if os.path.exists(dst):
            logger.info("DICOM copy skipped (already exists): %s", dst)
            return dst

        logger.info("Copying DICOM folder to output root: %s -> %s", src, dst)

        # Full copy (can be large!)
        shutil.copytree(src, dst, copy_function=shutil.copy2)

        logger.info("DICOM copy done: %s", dst)
        return dst
    def _copy_config_into_output_root(self):
        """
        Copy the config JSON used for this CT into the CT output root (beside the log),
        preserving the original filename.
        """
        logger = logging.getLogger(f"TDT_CONFIG_LOGGER_CT_{self.ct_indx}")

        src = os.path.abspath(self.config_path)
        if not os.path.isfile(src):
            raise FileNotFoundError(f"Config path is not a file: {src}")

        dst = os.path.join(self.output_folder_path, os.path.basename(src))

        # overwrite to guarantee it matches this run
        shutil.copy2(src, dst)
        logger.info("Copied config into output root: %s -> %s", src, dst)
        return dst
        
    def _log_config_to_file(self):
        '''
        For now only saved Config, will add on more logging later on
        '''
        log_path = os.path.join(self.output_folder_path, f"logging_file_CT_{self.ct_indx}.log")

        logger = logging.getLogger(f"TDT_CONFIG_LOGGER_CT_{self.ct_indx}")
        logger.setLevel(logging.INFO)
        logger.propagate = False  # prevent duplicate output if root logger is set elsewhere
        
        
        # Avoid adding multiple handlers if pipeline is constructed multiple times
        if not any(isinstance(h, logging.FileHandler) and getattr(h, "baseFilename", "") == log_path for h in logger.handlers):
            fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
            fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
            logger.addHandler(fh)
            
        logger.info("----Log started----") 
        
        logger.info("CT INDEX: %d", self.ct_indx)  # for validation study
        logger.info("DICOM PATH: %s", self.dicom_path)  # for validation study
            
        config = json.dumps(self.config, indent=2, sort_keys=True)
        logger.info("CONFIG SOURCE: %s", os.path.abspath(self.config_path))
        logger.info("CONFIG CONTENTS:\n%s", config)
        

    def setup_config(self, config_path):
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        with open(config_path, encoding="utf-8") as f:
            self.config = json.loads(json_minify(f.read()))

        output_folder_title = f"{self.config['output_folder']['title']}_CT_{self.ct_indx}"
        output_folder_path = os.path.join(self.current_dir_path, output_folder_title)
        os.makedirs(output_folder_path, exist_ok=True)

        # ✅ set attribute FIRST
        self.output_folder_path = output_folder_path

        # ✅ IMPORTANT: point stages to the CT-specific output root
        self.config["output_folder"]["title"] = self.output_folder_path

        # ✅ IMPORTANT: point stages to THIS CT’s DICOM folder
        self.config["ct_input"]["path1"] = self.dicom_path

        # create subdirs
        for _, name in self.config["subdir_names"].items():
            os.makedirs(os.path.join(self.output_folder_path, name), exist_ok=True)

        self._log_config_to_file()
        self._copy_config_into_output_root()
        self._copy_dicom_into_output_root()

    def run(self):
        logger = logging.getLogger(f"TDT_CONFIG_LOGGER_CT_{self.ct_indx}")
        t_pipeline = time.perf_counter()                 

        logger.info("Pipeline start | mode=%s", self.mode)  
        
        print("Loading TDT Pipeline configuration and setting up context manager...")
        logger = logging.getLogger(f"TDT_CONFIG_LOGGER_CT_{self.ct_indx}")
        context = Context(logger=logger)
        self.context = context
        print("TDT Pipeline configuration loaded and context manager set up.")
        print("Starting TDT Pipeline in mode:", self.mode) 
        
        # if final recon outputs exist, skip whole CT in production
        if str(self.mode).strip().lower() == "production":
            recon_done = all(
                os.path.exists(os.path.join(
                    self.output_folder_path,
                    self.config["subdir_names"]["spect_simulation"],
                    f'{self.config["spect_simulation"]["name"]}_{t}min.nii'
                ))
                for t in self.config["pbpk"]["FrameStartTimes"]
            )
            if recon_done:
                logger.info("Skipping CT_%d: recon outputs already exist.", self.ct_indx)
                return self.context
        
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
    config_path = "inputs/config_validation_study.json"
    parent_dir = '/home/jhubadmin/Theranostic-Virtual-Patient-Pipeline/inputs/validation_study_dataset/CT_Peter'
    
    
    ct_indx = 0
    for patient_folder in sorted(os.listdir(parent_dir)):
        pf = os.path.join(parent_dir, patient_folder)
        if not os.path.isdir(pf):
            continue
        for dicom_folder in sorted(os.listdir(pf)):
            dicom_path = os.path.join(pf, dicom_folder)
            if os.path.isdir(dicom_path):
                ct_indx += 1
                print(f"Processing CT {ct_indx}: {dicom_path}")
                pipeline = TdtPipeline(config_path, dicom_path, ct_indx)
                pipeline.run()
