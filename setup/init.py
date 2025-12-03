"""
TDT initialization utilities.

This module provides helper functions to:

- Load and parse the pipeline configuration JSON file.
- Create parent and per-module output directories.
- Configure logging to a per-run log file and record a summary
  of the current configuration.
"""

import json
import logging
import os

from json_minify import json_minify

logger = logging.getLogger(__name__)


def setup_config(output_folder_title, path):
    """
    Set up output folders and load configuration parameters.

    Parameters
    ----------
    output_folder_title : str
        Name of the parent output folder for this run.
    path : str
        Base path of the project (directory containing the setup folder).

    Returns
    -------
    output_dir_path : str
        Path to the parent output directory.
    output_subdir_paths : dict
        Mapping from module key -> subdirectory path.
    totseg_para : dict
        TotalSegmentator-related configuration parameters.
    pbpk_para : dict
        PBPK-related configuration parameters.
    simind_para : dict
        SIMIND-related configuration parameters.
    recon_para : dict
        Reconstruction-related configuration parameters.
    """
    config_file = os.path.join(path, "setup", "config.json")

    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Configuration file not found: {config_file}")

    # Load configuration JSON (json_minify removes comments and extra whitespace).
    with open(config_file, encoding="utf-8") as f:
        minified_json = json_minify(f.read())
        config = json.loads(minified_json)

    logger.info("Loading configuration from %s", config_file)

    # Create parent output folder
    output_dir_path = os.path.join(path, output_folder_title)
    os.makedirs(output_dir_path, exist_ok=True)
    logger.info("Created/verified parent output directory: %s", output_dir_path)

    # Create subfolder for each module
    output_names = config.get("OutputNames", {})
    output_subdir_paths = {}

    for key, name in output_names.items():
        subdir_path = os.path.join(output_dir_path, name)
        output_subdir_paths[key] = subdir_path
        os.makedirs(subdir_path, exist_ok=True)

    logger.debug("Per-module output directories: %s", output_subdir_paths)

    # TOTAL SEGMENTATOR PARAMETERS
    totseg_para = dict(config.get("TotalSegmentator", {}))

    # PBPK PARAMETERS
    pbpk_para = dict(config.get("PBPK", {}))

    # SIMIND PARAMETERS
    simind_para = dict(config.get("SIMIND", {}))

    # RECON PARAMETERS
    recon_para = dict(config.get("RECON", {}))

    logger.debug("TotalSegmentator config: %s", totseg_para)
    logger.debug("PBPK config: %s", pbpk_para)
    logger.debug("SIMIND config: %s", simind_para)
    logger.debug("Reconstruction config: %s", recon_para)

    return (
        output_dir_path,
        output_subdir_paths,
        totseg_para,
        pbpk_para,
        simind_para,
        recon_para,
    )


def setup_log(
    output_dir_path,
    output_subdir_paths,
    totseg_para,
    pbpk_para,
    simind_para,
    recon_para,
):
    """
    Set up logging configuration and write initial configuration summary.

    This function configures the root logger to write to a per-run log file and
    records a summary of the current TDT configuration. Other modules should
    obtain loggers via ``logging.getLogger(__name__)`` and will inherit this
    configuration.

    Parameters
    ----------
    output_dir_path : str
        Path to the parent output directory.
    output_subdir_paths : dict
        Mapping from module key -> subdirectory path.
    totseg_para : dict
        TotalSegmentator-related configuration parameters.
    pbpk_para : dict
        PBPK-related configuration parameters.
    simind_para : dict
        SIMIND-related configuration parameters.
    recon_para : dict
        Reconstruction-related configuration parameters.

    Returns
    -------
    None
    """
    # Create logging directory and file
    logs_dir = os.path.join(output_dir_path, "LOG_folder")
    os.makedirs(logs_dir, exist_ok=True)
    log_file = os.path.join(logs_dir, "TDT_Framework_Log.log")

    # Root logging configuration
    logging.basicConfig(
        level=logging.DEBUG,
        filename=log_file,
        encoding="utf-8",
        filemode="a",
        format="{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logger.info("=== TDT CT-to-SPECT configuration summary ===")

    logger.info("Output module paths (see DEBUG for details).")
    for key, subdir_path in output_subdir_paths.items():
        logger.debug("Output path for %s: %s", key, subdir_path)

    logger.info("TotalSegmentator parameters (see DEBUG for details).")
    for key, value in totseg_para.items():
        logger.debug("TotalSegmentator - %s: %s", key, value)

    logger.info("PBPK parameters (see DEBUG for details).")
    for key, value in pbpk_para.items():
        logger.debug("PBPK - %s: %s", key, value)

    logger.info("SIMIND parameters (see DEBUG for details).")
    for key, value in simind_para.items():
        logger.debug("SIMIND - %s: %s", key, value)

    logger.info("Reconstruction parameters (see DEBUG for details).")
    for key, value in recon_para.items():
        logger.debug("Reconstruction - %s: %s", key, value)
