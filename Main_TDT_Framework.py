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

from modules.RECON.runRECON import runRECON


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
    

    return path,config,out_paths,input_paths,totseg_para,pbpk_para,simind_para,recon_para

path, config, out_paths, input_paths, totseg_para , pbpk_para, simind_para, recon_para = setup()

def simulate():

    #####TOTAL SEGMENENTATION#####
    
    
    # Comment out if not using Mac :
    # IMPORTANT on macOS: spawn method prevents re-exec storms
    try:
        mp.set_start_method("spawn", force=True)
    except RuntimeError:
        pass
    # end of comment out
    
    print("[MAIN] Beginning TDT Program")
    print(f"[MAIN] Input file:{input_paths['ct_input_dicom']}")
    
    print(f"[MAIN] Beginning Segmentation using Total Segmentator")
    """
    runTOTSEG(input_paths['ct_input_dicom'],out_paths['output_total_seg'], totseg_para)
    """
    print(f"[MAIN] Segemenation Complete, ct input and segmentated output can be found:\n[MAIN]{out_paths['output_total_seg']}")
    
    
    #####NII PROCCESSING#####
    #had to resize to reduce time for the runPBPK to work 
    #(need to end up using full resolution at the end, for now this will work)
    
    resize_tuple=(128,128,128)
    print(f"[MAIN] Beginning Niifty file processing")
    print(f"[MAIN] Will need to resize ct arrays to {resize_tuple} for efficiancy")
    
    ts_classes = ts_cl.class_map["total"]
    ct_input_arr,segmentated_ml_output_arr, ct_input_arr_resize, segmentated_ml_output_arr_resize, class_seg, masks, atn_path, pixel_spacing,slice_thickness = NII_PROCCESSING(out_paths['output_total_seg'],ts_classes,resize_tuple,simind_para,totseg_para)
    print(f"[MAIN] CT seg. processed, atten bin can be found:\n[MAIN] {atn_path}")

    ##### PBPK #####
    FrameAct,act_path_all = runPBPK(out_paths,pbpk_para,totseg_para,ct_input_arr_resize,
                                    pixel_spacing,slice_thickness,resize_tuple,
                                    segmentated_ml_output_arr_resize,masks,class_seg)
    print(f"[MAIN] PBPK TAC created, activity bin can be found:\n[MAIN] {act_path_all}")
    
    ######SIMIND########
    runSIMIND(path, simind_para, pbpk_para, out_paths['output_SIMIND'],os.cpu_count(), 
              np.shape(ct_input_arr_resize)[0], 
              np.shape(segmentated_ml_output_arr_resize)[2],
              FrameAct, pixel_spacing, slice_thickness, 
              act_path_all, atn_path)
    
    print(f"[MAIN] SIMIND process finished, beginning RECON using PyTomography")
    
    ######RECON########
    runRECON(recon_para,pbpk_para,simind_para,out_paths,FrameAct)
    
    print(f"[MAIN] RECON process finished")
    
    return 0



if __name__ == "__main__":
    simulate()