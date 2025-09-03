'''
Main TDT Framework

Need to add description
'''
import json
import os
import numpy as np

#modules
'''
from modules.XCAT.runXCAT import runXCAT # shouldnt be running

from modules.PBPK.runPBPK import runPBPK
from modules.SIMIND_SPECT.runSIMIND import runSIMIND
from modules.GATE_PET.runGATE import runGATE
from modules.RECON.runRECON import runRECON
'''
import json, datetime
from json_minify import json_minify
from modules.TOTSEG.runTOTSEG import runTOTSEG

import modules.TOTSEG.totalseg_classmap as ts_cl
from modules.TOTSEG.nii_processing import NII_PROCCESSING

from modules.SIMIND_SPECT.runSIMIND import runSIMIND


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

# -------------------- End of comment (make sure to import os) --------------------



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

        
    return config,out_paths,input_paths,totseg_para


config, out_paths, input_paths, totseg_para = setup()

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
    
    
    '''
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
    '''
    print(f"[MAIN] Segemenation Complete, ct input and segmentated output can be found:\n[MAIN]{out_paths['output_total_seg']}/{output_file_name_totseg}")
    
    
    #####NII PROCCESSING#####
    ts_classes = ts_cl.class_map["total"]
    ct_input_arr,segmentated_ml_output_arr, class_seg, masks, atn_av_path = NII_PROCCESSING(out_paths['output_total_seg'],output_file_name_totseg,ts_classes)
    #print(f"[MAIN] CT seg. processed, atten bin can be found:\n[MAIN]{atn_av_path}")
    
    act_path = '/Users/peteryazdi/Desktop/BC_Cancer/TDT/Output_2025-09-03/TotalSegmentator_Outputs/2025-09-03_Segmenation/_act_av.bin'
    atn_path = '/Users/peteryazdi/Desktop/BC_Cancer/TDT/Output_2025-09-03/TotalSegmentator_Outputs/2025-09-03_Segmenation/_atn_av.bin'
    
    ImgSum=sum(segmentated_ml_output_arr)
    ImgLength= np.shape(segmentated_ml_output_arr)[2]
    num_cores = os.cpu_count()
    output_name = "fake_pbpk_run1"
    path = os.path.dirname(__file__)
    output_path = "Output_2025-09-03/TotalSegmentator_Outputs"
    runSIMIND(path, config, output_name, output_path, num_cores, ImgLength, ImgSum,act_path,atn_path)
    

    
    
    return

if __name__ == "__main__":
    simulate()