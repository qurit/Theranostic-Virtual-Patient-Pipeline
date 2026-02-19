#imports
import os
import json
from json_minify import json_minify
import logging
import time
import shutil
import argparse
from distutils.util import strtobool  
import copy

from context import Context

from stages.spect_pre_process.segmentation_stage import TotalSegmentationStage
from stages.spect_pre_process.unify_ts_outputs import TdtRoiUnifyStage
from stages.spect_pre_process.preprocessing_stage import SimindPreprocessStage

from stages.pbpk.pbpk_stage import PbpkStage
from stages.spect_simulation.simind_stage import SimindSimulationStage
from stages.spect_simulation.reconstruction_stage import SpectReconstructionStage


class TdtPipeline:
    """
    Orchestrates the full TDT pipeline for a single CT input.

    Args:
        config_path (str): Path to the JSON config file (comments allowed via json_minify).
        ct_input (str): Path to a CT input (either a .nii/.nii.gz file OR a DICOM directory).
        ct_indx (int): Index used for naming (e.g., output folder suffix "_CT_{ct_indx}").
        logging_on (bool): If True, writes a per-CT log file into the CT output folder.
        save_ct_scan (bool): If True, copies the CT input into the output folder for provenance.
        mode (str): "DEBUG" or "PRODUCTION" (affects logging verbosity).
    """

    def __init__(self, config_path, ct_input, ct_indx, logging_on=True, save_ct_scan=False, mode='PRODUCTION'):
        self.config_path = config_path
        self.ct_input = ct_input
        self.ct_indx = ct_indx
        self.current_dir_path = os.path.abspath(os.path.dirname(__file__))

        self.logging_on = logging_on
        self.save_ct_scan = save_ct_scan
        self.mode = mode

        self._config_setup(config_path)

        self.logger = logging.getLogger(f"TDT_CONFIG_LOGGER_CT_{self.ct_indx}")
        self.logger.setLevel(logging.DEBUG if self.mode == "DEBUG" else logging.INFO)
        self.logger.propagate = False

        if self.logging_on:
            self._log_setup()
        else:
            self.logger.disabled = True

        self._context_setup()
        if not self.logging_on:
            self.context._log_enabled = False

    def _save_ct_input_copy(self):
        """
        Copy the CT input into the output folder under `ct_input_copy/`.

        Accepts:
            - If self.ct_input is a file: copies file into ct_input_copy/
            - If self.ct_input is a directory: copies entire directory into ct_input_copy/<dirname>/
        """
        dst_root = os.path.join(self.output_folder_path, "ct_input_copy")
        os.makedirs(dst_root, exist_ok=True)

        if os.path.isfile(self.ct_input): # nii file
            shutil.copy2(self.ct_input, os.path.join(dst_root, os.path.basename(self.ct_input)))
        else: # dicom dir
            dst = os.path.join(dst_root, os.path.basename(os.path.normpath(self.ct_input)))
            if not os.path.exists(dst):
                shutil.copytree(self.ct_input, dst)

    def _log_setup(self):
        """
        Configure a per-CT log file handler writing to:
            <output_folder_path>/logging_file_CT_<ct_indx>.log
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

    def _config_setup(self, config_path):
        """
        Load config from disk and prepare output folder + subdirectories.

        Args:
            config_path (str): Path to JSON config. (May include // comments, stripped via json_minify.)
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

        # Create subdirs
        for _, name in self.config["subdir_names"].items():
            os.makedirs(os.path.join(self.output_folder_path, name), exist_ok=True)

    def _context_setup(self):
        """
        Create and populate the Context object.

        This stores:
          - A snapshot of the user config (deep-copied)
          - Runtime metadata (mode, ct_input_path, ct_input_type, ct_indx, output_folder_path)
          - Computed subdir paths
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

    def run(self):
        """
        Execute the pipeline stages sequentially for this CT.

        Returns:
            Context: The updated context after all stages complete.
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


def build_arg_parser():
    """
    Build the CLI argument parser.

    Returns:
        argparse.ArgumentParser: Configured parser for main.py
    """
    parser = argparse.ArgumentParser(
        description="Theranostic Digital Twin (TDT) Pipeline Runner"
    )

    # REQUIRED
    parser.add_argument(
        "--config_file",
        required=True,
        type=str,
        help="Path to the pipeline config JSON, typically /inputs/config.json",
    )
    parser.add_argument(
        "--input_ct_dir",
        required=True,
        type=str,
        help="Path to CT inputs directory, typically /inputs/ct_inputs/",
    )

    # Extra
    parser.add_argument(
        "--logging_on",
        default=True,
        type=lambda x: bool(strtobool(x)),
        help="Enable file logging (true/false). Default: true",
    )
    parser.add_argument(
        "--save_ct_scan",
        default=False,
        type=lambda x: bool(strtobool(x)),
        help="Copy CT input into output folder for provenance (true/false). Default: false",
    )
    parser.add_argument(
        "--mode",
        default='PRODUCTION',
        choices=["DEBUG", "PRODUCTION"],
        help=(
            "DEBUG is better for single runs while PRODUCTION saves space for "
            "multi-patient runs. Default PRODUCTION"
        ),
    )

    return parser


if __name__ == "__main__":
    parser = build_arg_parser()
    args = parser.parse_args()

    # Require a directory and iterate directly
    ct_inputs_dir = os.path.abspath(args.input_ct_dir)
    if not os.path.isdir(ct_inputs_dir):
        raise NotADirectoryError(f"input_ct_dir must be a directory: {ct_inputs_dir}")

    # so indx isnt thrown off by hidden files; also sort for consistency
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
                mode=args.mode,
            )
            pipeline.run()
        except Exception as e:
            print(f"[ERROR] CT index {idx} failed for input: {ct_path}\n{e}")
