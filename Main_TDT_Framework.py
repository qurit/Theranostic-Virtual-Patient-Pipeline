'''
Main TDT Framework

Need to add description
'''
import json
import os
import numpy as np
import matplotlib.pyplot as plt

#modules
'''
from modules.GATE_PET.runGATE import runGATE
from modules.RECON.runRECON import runRECON
'''
import json, datetime
from json_minify import json_minify
from modules.TOTSEG.runTOTSEG import runTOTSEG

import modules.TOTSEG.totalseg_classmap as ts_cl
from modules.TOTSEG.nii_processing import NII_PROCCESSING

from modules.SIMIND_SPECT.runSIMIND import runSIMIND

from modules.PBPK.runPBPK import runPBPK


#Comment out if not using Mac :
# -------------------- Environment & worker limits (set BEFORE heavy imports) --------------------
import multiprocessing as mp

# Keep parallelism tame on laptops (adjust if you have more cores/RAM)
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "4")
os.environ.setdefault("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS", "2")

# Reduce internal workers used by TotalSegmentator/MONAI to avoid CPU thrash on macOS
os.environ.setdefault("TOTALSEGMENTATOR_NUM_WORKERS", "1")
os.environ.setdefault("MONAI_NUM_WORKERS", "0")

# Allow Torch to fall back gracefully when MPS ops are missing
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

# -------------------- End of comment--------------------

def setup():
    path = os.path.dirname(__file__) #finds where this folder is
    config_file = os.path.join(path, 'Configure_Files/config.json')
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Configuration file not found: {config_file}")
    with open(config_file) as f:
        minified_json = json_minify(f.read())
        config = json.loads(minified_json)
        
    #Output Folder 
    output_path = f"{path}/Output_{datetime.date.today()}"
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
        
    # --- Create subfolders from SIMIND parameters ---
    simind_names = config.get("SIMIND", {})
    simind_para = {}
    for key, name in simind_names.items():
        simind_para[key] = name
    
    # --- Create subfolders from PBPK parameters ---
    pbpk_names = config.get("PBPK", {})
    pbpk_para = {}
    for key, name in pbpk_names.items():
        pbpk_para[key] = name


        
    return path,config,out_paths,input_paths,totseg_para,simind_para,pbpk_para


path, config, out_paths, input_paths, totseg_para , simind_para, pbpk_para = setup()

def simulate():

    #####TOTAL SEGMENENTATION#####
    # Comment out if not using Mac :
    # IMPORTANT on macOS: spawn method prevents re-exec storms
    try:
        mp.set_start_method("spawn", force=True)
    except RuntimeError:
        pass
    # end of comment out
    output_file_name_totseg = f"{datetime.date.today()}_Segmenation"
    
    """
    runTOTSEG(
    input_path=input_paths['ct_input_dicom'],
    output_path=out_paths['output_total_seg'],
    output_name = output_file_name_totseg,
    device=totseg_para['device'], 
    fast=totseg_para['fast'],
    roi_subset=totseg_para['roi_subset'],          
    ml=totseg_para['ml'],                 
    statistics=totseg_para['statistics'],         
    radiomics=totseg_para['radiomics'],          
    )
    """
    
    print(f"[MAIN] Segemenation Complete, ct input and segmentated output can be found:\n[MAIN]{out_paths['output_total_seg']}/{output_file_name_totseg}")
    
    #####NII PROCCESSING#####
    ts_classes = ts_cl.class_map["total"]
    ct_input_arr,segmentated_ml_output_arr, class_seg, masks, atn_path, pixel_spacing_x,slice_thickness = NII_PROCCESSING(out_paths['output_total_seg'],output_file_name_totseg,ts_classes)
    print(f"[MAIN] CT seg. processed, atten bin can be found:\n[MAIN] {atn_path}")

    ##### PBPK #####
    FrameAct,act_path = runPBPK(out_paths,pbpk_para,totseg_para,ct_input_arr,pixel_spacing_x,slice_thickness,segmentated_ml_output_arr,masks,class_seg)
    print(f"[MAIN] PBPK TAC created, activity bin can be found:\n[MAIN] {act_path}")
    quit()

    ######SIMIND########

    runSIMIND(path=path, simind_para = simind_para, pbpk_para=pbpk_para,
              output_path=out_paths['output_SIMIND'],num_cores=1, 
              InputImageSize=np.shape(ct_input_arr)[0], 
              ImgLength=np.shape(segmentated_ml_output_arr)[2],
              ImgSum=FrameAct, 
              pixel_spacing_x=pixel_spacing_x, 
              slice_thickness=slice_thickness, act_path=act_path, 
              atn_path=atn_path)
    

    
    
    return

if __name__ == "__main__":
    simulate()