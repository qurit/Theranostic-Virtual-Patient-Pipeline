"""
Theranostic Digital Twin (TDT) Pipeline Runner.

This module provides:
- `TdtPipeline`: an orchestrator that runs all pipeline stages for a single CT input.
- A CLI entrypoint that iterates through a directory of CT inputs and runs the pipeline.

Notes:
- A CT input may be either a NIfTI file (.nii / .nii.gz) or a DICOM directory.
- The pipeline writes outputs into an output folder derived from the config and CT index.

For any questions or issues, please contact: pyazdi@bccrc.ca
"""

from __future__ import annotations

# -----------------------------
# Standard library imports
# -----------------------------
import os
import json
from json_minify import json_minify
import logging
import time
import shutil
import argparse
import copy
from typing import Any, Dict, Literal, Optional

# -----------------------------
# Local imports
# -----------------------------
from src.io.context import Context

from src.stages.spect_pre_process.segmentation_stage import TotalSegmentationStage
from src.stages.spect_pre_process.unify_ts_outputs import TdtRoiUnifyStage
from src.stages.spect_pre_process.preprocessing_stage import SimindPreprocessStage
from src.stages.synthetic_lesions.synthetic_lesions_stage import SyntheticLesionsStage

from src.stages.pbpk.pbpk_stage import PbpkStage
from src.stages.spect_simulation.simind_stage import SimindSimulationStage
from src.stages.spect_simulation.reconstruction_stage import SpectReconstructionStage


CTInputType = Literal["nii", "dicom"]


class TdtPipeline:
    """
    Orchestrates the full TDT pipeline for a single CT input.

    Parameters
    ----------
    config_path : str
        Path to the JSON config file (comments allowed via `json_minify`).
    ct_input : str
        Path to a CT input (either a .nii/.nii.gz file OR a DICOM directory).
    ct_indx : int
        Index used for naming (e.g., output folder suffix "_CT_{ct_indx}").
    logging_on : bool, default=True
        If True, writes a per-CT log file into the CT output folder.
    save_ct_scan : bool, default=False
        If True, copies the CT input into the output folder for provenance.
    save_config : bool, default=False
        If True, saves a copy of the config JSON into the output folder. 
    mode : {"DEBUG", "PRODUCTION"}, default="PRODUCTION"
        Affects logging verbosity.
    synthetic_lesions : bool, default=False
        If True, runs the synthetic lesions stage to generate lesions in the CT scan.

    Attributes
    ----------
    config : dict[str, Any]
        Parsed configuration loaded from `config_path`.
    output_folder_path : str
        Absolute path to the CT-specific output folder.
    ct_input_type : {"nii", "dicom"}
        Determined input type for `ct_input`.
    context : Context
        Runtime context object passed through pipeline stages.
    logger : logging.Logger
        Per-CT logger, optionally writing to a file handler.
    """

    def __init__(
        self,
        config_path: str,
        ct_input: str,
        ct_indx: int,
        logging_on: bool = True,
        save_ct_scan: bool = False,
        save_config: bool = False,  
        mode: Literal["DEBUG", "PRODUCTION"] = "PRODUCTION",
        synthetic_lesions: bool = False,
    ) -> None:
        self.config_path: str = config_path
        self.ct_input: str = ct_input
        self.ct_indx: int = ct_indx
        self.current_dir_path: str = os.path.abspath(os.path.dirname(__file__))

        self.logging_on: bool = logging_on  # if True, enables file logging in the output folder
        self.save_ct_scan: bool = save_ct_scan  # if True, saves a copy of the CT input in the output folder for provenance
        self.save_config: bool = save_config  # if True, saves config json into the CT output folder  
        self.mode: Literal["DEBUG", "PRODUCTION"] = mode
        self.synthetic_lesions: bool = synthetic_lesions  # if True, runs the synthetic lesions stage to generate lesions in the CT scan

        self.config: Dict[str, Any] = {}  # will be populated in _config_setup()
        self.output_folder_path: str = ""  # will be set in _config_setup()
        self.ct_input_type: CTInputType = "dicom"  # default, will be set properly in _config_setup() after validation

        self._config_setup(config_path)

        self.logger: logging.Logger = logging.getLogger(f"TDT_CONFIG_LOGGER_CT_{self.ct_indx}")
        self.logger.setLevel(logging.DEBUG if self.mode == "DEBUG" else logging.INFO)
        self.logger.propagate = False

        if self.logging_on:
            self._log_setup()
        else:
            self.logger.disabled = True

        self.context: Context
        self._context_setup()
        if not self.logging_on:
            self.context._log_enabled = False
            
        self.run_synthetic_lesions = self.synthetic_lesions 
        if self.run_synthetic_lesions and self.config["synthetic_lesions"]['specs'] is not None:
            self.run_synthetic_lesions = True
        elif self.run_synthetic_lesions and self.config["synthetic_lesions"]['specs'] is None:
            self.logger.info("Synthetic lesions stage is disabled based on missing specs in config.")
            self.run_synthetic_lesions = False

    def _save_ct_input_copy(self) -> None:
        """
        Copy the CT input into the output folder under `ct_input_copy/`.

        Behavior
        --------
        - If `self.ct_input` is a file: copies the file into `ct_input_copy/`
        - If `self.ct_input` is a directory: copies the entire directory into
          `ct_input_copy/<dirname>/`
        """
        dst_root = os.path.join(self.output_folder_path, "ct_input_copy")
        os.makedirs(dst_root, exist_ok=True)

        if os.path.isfile(self.ct_input):  # nii file
            shutil.copy2(self.ct_input, os.path.join(dst_root, os.path.basename(self.ct_input)))
        else:  # dicom dir
            dst = os.path.join(dst_root, os.path.basename(os.path.normpath(self.ct_input)))
            if not os.path.exists(dst):
                shutil.copytree(self.ct_input, dst)

    def _save_config_copy(self, config_path: str) -> None:  
        """
        Save a copy of the config JSON into the output folder.

        Writes:
            <output_folder_path>/config.json
        """
        dst = os.path.join(self.output_folder_path, "config.json")
        shutil.copy2(config_path, dst)

    def _log_setup(self) -> logging.Logger:
        """
        Configure a per-CT log file handler writing to:
            <output_folder_path>/logging_file_CT_<ct_indx>.log

        Returns
        -------
        logging.Logger
            The configured per-CT logger.
        """
        log_path = os.path.join(self.output_folder_path, f"logging_file_CT_{self.ct_indx}.log")
        logger = self.logger

        # Avoid adding multiple handlers if pipeline is constructed multiple times
        if not any(
            isinstance(h, logging.FileHandler) and getattr(h, "baseFilename", "") == log_path
            for h in logger.handlers
        ):
            fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
            fh.setFormatter(
                logging.Formatter(
                    "%(asctime)s | %(levelname)s | %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S",
                )
            )
            logger.addHandler(fh)

        logger.info("----Log started----")
        logger.info("Output folder: %s", self.output_folder_path)
        return logger

    def _config_setup(self, config_path: str) -> None:
        """
        Load config from disk and prepare output folder + subdirectories.

        Parameters
        ----------
        config_path : str
            Path to JSON config. (May include // comments, stripped via json_minify.)

        Raises
        ------
        FileNotFoundError
            If the config file does not exist or the CT input path is invalid.
        ValueError
            If a file input is provided but is not a supported NIfTI extension.
        """
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        with open(config_path, encoding="utf-8") as f:
            self.config = json.loads(json_minify(f.read()))

        output_folder_title = f"{self.config['output_folder']['title']}_CT_{self.ct_indx}"
        self.output_folder_path = os.path.join(self.current_dir_path, output_folder_title)
        os.makedirs(self.output_folder_path, exist_ok=True)

        # Validate + tag ct input (file must be nifti; dir assumed dicom)
        if os.path.isfile(self.ct_input):
            ct_lower = self.ct_input.lower()
            if not (ct_lower.endswith(".nii") or ct_lower.endswith(".nii.gz")):
                raise ValueError(f"CT input file must be .nii or .nii.gz, got: {self.ct_input}")
            self.ct_input_type = "nii"
        elif os.path.isdir(self.ct_input):
            self.ct_input_type = "dicom"
        else:
            raise FileNotFoundError(f"CT input not found: {self.ct_input}")

        # Optionally save a copy of the CT input
        if self.save_ct_scan:
            self._save_ct_input_copy()

        # Optionally save a copy of the config
        if self.save_config:  
            self._save_config_copy(config_path)  

        # Create subdirs
        for _, name in self.config["subdir_names"].items():
            os.makedirs(os.path.join(self.output_folder_path, name), exist_ok=True)

    def _context_setup(self) -> None:
        """
        Create and populate the Context object.

        The Context stores:
        - A deep-copied snapshot of relevant config sections
        - Runtime metadata (mode, CT path/type/index, output folder)
        - Computed subdir paths used across stages
        """
        context = Context(logger=self.logger)
        self.context = context

        # Dump config snapshot into context FIRST (so subdir_names exists)
        context.config = copy.deepcopy(self.config)
        context.output_cfg = copy.deepcopy(self.config.get("output_folder", {}))
        context.subdir_names = copy.deepcopy(self.config.get("subdir_names", {}))
        context.spect_preprocessing_cfg = copy.deepcopy(self.config.get("spect_preprocessing", {}))
        context.pbpk_cfg = copy.deepcopy(self.config.get("pbpk", {}))
        context.spect_simulation_cfg = copy.deepcopy(self.config.get("spect_simulation", {}))

        # Initial setup for Context (runtime)
        context.mode = self.mode
        context.ct_input_path = self.ct_input
        context.ct_input_type = self.ct_input_type
        context.ct_indx = self.ct_indx
        context.output_folder_path = self.output_folder_path

        # Computed paths (now subdir_names exists)
        context.subdir_paths = {
            k: os.path.join(self.output_folder_path, name)
            for k, name in context.subdir_names.items()
        }

        self.logger.debug("Context initialized for CT_%s", self.ct_indx)

    def run(self) -> Context:
        """
        Execute the pipeline stages sequentially for this CT.

        Stages:
        1. TotalSegmentator
        2. Unification of TS outputs to TDT ROIs
        2.5. Generate synthetic lesions (if enabled)
        3. Preprocess for SIMIND
        4. PBPK
        5. SIMIND Simulation
        6. SPECT Reconstruction

        Returns
        -------
        Context
            The updated context after all stages complete.

        Notes
        -----
        Each stage is expected to implement `.run()` and return an updated Context.
        """
        logger = self.logger
        t_pipeline = time.perf_counter()

        print(f"Starting TDT Pipeline (CT_{self.ct_indx}) | mode={self.mode} | input={self.ct_input}")

        context = self.context

        logger.info("Pipeline start | mode=%s", self.mode)
        logger.info("CT input | path=%s | type=%s", self.ct_input, self.ct_input_type)

        # -----------------------------
        # Stage 1: TotalSegmentator
        # -----------------------------
        logger.info("Stage start: TotalSegmentator")
        t_stage = time.perf_counter()

        print("Running TotalSegmentator Stage...")
        context = TotalSegmentationStage(context).run()
        print("TotalSegmentator Stage completed.")

        logger.info("Stage end: TotalSegmentator | elapsed=%.2fs", time.perf_counter() - t_stage)

        # -----------------------------
        # Stage 2: Unification of TS outputs to TDT ROIs
        # -----------------------------
        logger.info("Stage start: TDT ROI Unification")
        t_stage = time.perf_counter()

        print("Running TDT ROI Unification Stage...")
        context = TdtRoiUnifyStage(context).run()
        print("TDT ROI Unification Stage completed.")

        logger.info("Stage end: TDT ROI Unification | elapsed=%.2fs", time.perf_counter() - t_stage)
        
        
        # -----------------------------
        # Stage 2.5: Generate synthetic lesions (if enabled)
        # -----------------------------
        if self.run_synthetic_lesions:
            logger.info("Stage start: Synthetic Lesions Generation")
            t_stage = time.perf_counter()

            print("Running Synthetic Lesions Generation Stage...")
            context = SyntheticLesionsStage(context).run()
            print("Synthetic Lesions Generation Stage completed.")

            logger.info("Stage end: Synthetic Lesions Generation | elapsed=%.2fs", time.perf_counter() - t_stage)

        # -----------------------------
        # Stage 3: Preprocess for SIMIND
        # -----------------------------
        logger.info("Stage start: SIMIND Preprocessing")
        t_stage = time.perf_counter()

        print("Running SIMIND Preprocessing Stage...")
        context = SimindPreprocessStage(context).run()
        print("SIMIND Preprocessing Stage completed.")

        logger.info("Stage end: SIMIND Preprocessing | elapsed=%.2fs", time.perf_counter() - t_stage)
        
        # -----------------------------
        # Stage 4: PBPK
        # -----------------------------
        logger.info("Stage start: PBPK")
        t_stage = time.perf_counter()

        print("Running PBPK Stage...")
        context = PbpkStage(context).run()
        print("PBPK Stage completed.")

        logger.info("Stage end: PBPK | elapsed=%.2fs", time.perf_counter() - t_stage)
        break # 
        # -----------------------------
        # Stage 5: SIMIND
        # -----------------------------
        logger.info("Stage start: SIMIND Simulation")
        t_stage = time.perf_counter()

        print("Running SIMIND Simulation Stage...")
        context = SimindSimulationStage(context).run()
        print("SIMIND Simulation Stage completed.")

        logger.info("Stage end: SIMIND Simulation | elapsed=%.2fs", time.perf_counter() - t_stage)

        # -----------------------------
        # Stage 6: Recon
        # -----------------------------
        logger.info("Stage start: SPECT Reconstruction")
        t_stage = time.perf_counter()

        print("Running SPECT Reconstruction Stage...")
        context = SpectReconstructionStage(context).run()
        print("SPECT Reconstruction Stage completed.")

        logger.info("Stage end: SPECT Reconstruction | elapsed=%.2fs", time.perf_counter() - t_stage)

        logger.info("Pipeline end | total_elapsed=%.2fs", time.perf_counter() - t_pipeline)

        print("TDT Pipeline completed successfully.")
        return context


def build_arg_parser() -> argparse.ArgumentParser:
    """
    Build the CLI argument parser.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser for the pipeline CLI.
    """
    parser = argparse.ArgumentParser(
        description="Theranostic Digital Twin (TDT) Pipeline Runner"
    )

    # Required arguments
    parser.add_argument("--config_file", required=True, type=str)
    parser.add_argument("--input_ct_dir", required=True, type=str)

    # Optional flags
    parser.add_argument(
        "--logging_on",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable file logging. Use --logging_on / --no-logging_on. Default: enabled",
    )
    parser.add_argument(
        "--save_ct_scan",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Copy CT input into output folder. Use --save_ct_scan / --no-save_ct_scan. Default: disabled",
    )
    parser.add_argument(  
        "--save_config",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Copy the config JSON into each CT output folder. Use --save_config / --no-save_config. Default: disabled",
    )  

    parser.add_argument(
        "--mode",
        default="PRODUCTION",
        choices=["DEBUG", "PRODUCTION"],
    )
    
    parser.add_argument(
        "--synthetic_lesions",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Generate synthetic lesions. Use --synthetic_lesions / --no-synthetic_lesions. Default: disabled",
    )

    return parser


def main() -> int:
    """
    CLI entrypoint.

    Returns
    -------
    int
        Process exit code (0 = success).
    """
    parser = build_arg_parser()
    args = parser.parse_args()

    # Require a directory and iterate directly
    ct_inputs_dir = os.path.abspath(args.input_ct_dir)
    if not os.path.isdir(ct_inputs_dir):
        raise NotADirectoryError(f"input_ct_dir must be a directory: {ct_inputs_dir}")

    # Filter hidden files/directories and keep deterministic ordering for repeatability
    items = [n for n in sorted(os.listdir(ct_inputs_dir)) if not n.startswith(".")]

    print(f"Discovered {len(items)} CT item(s) in: {ct_inputs_dir}")

    for idx, name in enumerate(items):
        ct_path = os.path.join(ct_inputs_dir, name)

        try:
            pipeline = TdtPipeline(
                config_path=args.config_file,  # required by user
                ct_input=ct_path,
                ct_indx=idx,
                logging_on=args.logging_on,
                save_ct_scan=args.save_ct_scan,
                save_config=args.save_config, 
                mode=args.mode,
                synthetic_lesions=args.synthetic_lesions,
            )
            pipeline.run()
        except Exception as e:
            # Keep console-friendly failure reporting per CT input
            print(f"[ERROR] CT index {idx} failed for input: {ct_path}\n{e}")

    return 0  # 0 indicates successful completion of the CLI process


if __name__ == "__main__":
    raise SystemExit(main())