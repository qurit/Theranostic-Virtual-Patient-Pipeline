'''
Main TDT Framework

Need to add description
'''
import json, os
from json_minify import json_minify
import logging as log
from modules.TOTSEG.runTOTSEG import runTOTSEG
import modules.TOTSEG.totalseg_classmap as ts_cl
from modules.TOTSEG.nii_processing import NII_PROCCESSING
from modules.SIMIND_SPECT.runSIMIND import runSIMIND
from modules.PBPK.runPBPK import runPBPK
from modules.RECON.runRECON import runRECON


output_folder_title = 'Output_Folder_CTvsSPECT_Comparison'

def setup():
    path = os.path.dirname(__file__) #finds where this folder is
    config_file = os.path.join(path, 'Configure_Files/config.json')
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Configuration file not found: {config_file}")
    with open(config_file) as f:
        minified_json = json_minify(f.read())
        config = json.loads(minified_json)
        
    #Output Folder 
    output_path = f"{path}/{output_folder_title}"
    os.makedirs(output_path, exist_ok =True)
    
    # --- Create subfolders from OutputNames and return paths ---
    output_names = config.get("OutputNames", {})
    out_paths = {}
    for key, name in output_names.items():
        out_paths[key] = f"{output_path}/{name}"
        os.makedirs(out_paths[key], exist_ok =True)
        
    # --- Create subfolders from Input and return paths ---
    input_names = config.get("InputPaths", {})
    input_paths = {}
    for key, name in input_names.items():
        input_paths[key] = name
        
    # --- Create subfolders from Total Seg parameters ---
    totseg_names = config.get("TotalSegmentator", {})
    totseg_para = {}
    for key, name in totseg_names.items():
        totseg_para[key] = name
    
    # --- Create subfolders from PBPK parameters ---
    pbpk_names = config.get("PBPK", {})
    pbpk_para = {}
    for key, name in pbpk_names.items():
        pbpk_para[key] = name
        
    # --- Create subfolders from SIMIND parameters ---
    simind_names = config.get("SIMIND", {})
    simind_para = {}
    for key, name in simind_names.items():
        simind_para[key] = name
        
    # --- Create subfolders from RECON parameters ---
    recon_names = config.get("RECON", {})
    recon_para = {}
    for key, name in recon_names.items():
        recon_para[key] = name
    

    return path,config,out_paths,input_paths,totseg_para,pbpk_para,simind_para,recon_para, output_path

path, config, out_paths, input_paths, totseg_para , pbpk_para, simind_para, recon_para, output_path = setup()


logs_dir = os.path.join(output_path, "logs")
os.makedirs(logs_dir, exist_ok=True)
log_file = os.path.join(logs_dir, "TDT_Framework_Log.log")

#logging stuff
log.basicConfig(
    level=log.DEBUG,
    filename=f"{log_file}",
    encoding="utf-8",
    filemode="a",
    format="{asctime} - {levelname} - {message}",
    style="{", 
    datefmt="%Y-%m-%d %H:%M:%S")

def log_config(out_paths, input_paths, totseg_para, pbpk_para, simind_para, recon_para):
    """Log all configuration parameters at startup"""
    log.info("=== TDT Framework Configuration ===")
    
    log.info("Output Paths:")
    for key, path in out_paths.items():
        log.info(f"  {key}: {path}")
    
    log.info("Input Paths:")
    for key, path in input_paths.items():
        log.info(f"  {key}: {path}")
        
    log.info("Total Segmentator Parameters:")
    for key, value in totseg_para.items():
        log.info(f"  {key}: {value}")
        
    log.info("PBPK Parameters:")
    for key, value in pbpk_para.items():
        log.info(f"  {key}: {value}")
        
    log.info("SIMIND Parameters:")
    for key, value in simind_para.items():
        log.info(f"  {key}: {value}")
        
    log.info("Reconstruction Parameters:")
    for key, value in recon_para.items():
        log.info(f"  {key}: {value}")
    
    log.info("=== End Configuration ===\n")

log_config(out_paths, input_paths, totseg_para, pbpk_para, simind_para, recon_para)


def simulate():
    """
    need to fill
    """
    #####TOTAL SEGMENENTATION#####
    log.info("Beginning TDT Program")
    print("[MAIN] Beginning TDT Program...")
    
    log.debug(f"Input file:{input_paths['ct_spect_input_dicom']}")
    
    log.info("Beginning Segmentation using Total Segmentator")
    print("[MAIN] Beginning Segmentation using Total Segmentator...")

    #ml_file,body_file = runTOTSEG(input_paths['ct_spect_input_dicom'],out_paths['output_total_seg'], totseg_para)

    ml_file= '/home/jhubadmin/Theranostic-Virtual-Patient-Pipeline/Output_Folder_CTvsSPECT_Comparison/TOTSEG_Outputs/TOTSEG_ml_segmentation.nii.gz'
    body_file = '/home/jhubadmin/Theranostic-Virtual-Patient-Pipeline/Output_Folder_CTvsSPECT_Comparison/TOTSEG_Outputs/TOTSEG_body_segmentation.nii.gz'

    print("[MAIN] Segmentation Complete")

    log.debug(f"CT input and Segmentated Output can be found:{out_paths['output_total_seg']}")

    
    #####NII PROCCESSING#####
    log.info("Beginning Nifti file processing")
    print("[MAIN] Beginning Nifti file processing...")
    ts_classes = ts_cl.class_map["total"]
    ct_input_arr,segmentated_ml_output_arr, segmentated_body_output_arr, class_seg, masks, atn_path, seg_ml_bin_path, seg_body_bin_path, pixel_spacing,slice_thickness,ct_get_zoom = NII_PROCCESSING(out_paths['output_total_seg'],ts_classes,simind_para,totseg_para,ml_file,body_file)
    log.debug(f"CT seg. processed, atten bin can be found:{atn_path}")
    
    
    ##### PBPK #####
    log.info("Beginning PBPK modelling")
    print("[MAIN] Beginning PBPK modelling...")
    ActivityMapSum,ActivityOrganSum, act_path_all_organ, act_path_all_map = runPBPK(out_paths,pbpk_para,segmentated_ml_output_arr,masks,class_seg,ct_get_zoom)
    log.debug(f"PBPK TAC created, activity bin can be found: {act_path_all_map}")


    ######SIMIND########
    log.info("Beginning SIMIND simulation")
    print("[MAIN] Beginning SIMIND simulation...")
    runSIMIND(path, class_seg, simind_para, pbpk_para, out_paths['output_SIMIND'],os.cpu_count(), 
              segmentated_ml_output_arr, ActivityOrganSum, ActivityMapSum, pixel_spacing, slice_thickness, 
              act_path_all_organ, atn_path)
    
    
    ######RECON########
    log.info("Beginning RECON using PyTomography")
    print("[MAIN] Beginning RECON using PyTomography...")
    runRECON(recon_para,pbpk_para,simind_para,out_paths, class_seg, ActivityMapSum)
    
    log.info("TDT Framework finished")
    print("[MAIN] TDT Framework finished")
    return 1



if __name__ == "__main__":
    simulate()