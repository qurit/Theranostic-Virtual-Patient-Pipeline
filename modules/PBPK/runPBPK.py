import numpy as np
import pycno, os
import logging as log
# ensure a non-GUI backend so savefig works on headless servers
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage
import math


def runPBPK(out_paths, pbpk_para,seg_plus_body_arr,masks,class_seg,ct_get_zoom):
    pbpk_name = pbpk_para['name'] 
    act_path_all_map = []
    act_path_all_organ = []
    ActivityOrganSum = {}
    del class_seg['Background']
    VOIs_possible = ['Tumor1', 'Tumor2', 'Kidney', 'Heart', 'SG', 'Bone', 'TumorRest', 'Spleen', 'Liver', 'Prostate', 'GI', 'Rest', 'Skin', 'Muscle', 'Brain', 'RedMarrow', 'Lungs', 'Adipose']
    roi_2_VOI = {
#        "NA":'Tumor1', 
#        "NA":'Tumor2', 
        "kidney":'Kidney', 
        "heart":'Heart', 
        "body" : 'Rest',
        #'SG',
        #'Bone',
        #'TumorRest', 
        'spleen':'Spleen', 
        'liver':'Liver', 
        "prostate":'Prostate',
        #'GI', 
        #'Rest',
        #'Skin',
        #'Muscle',
        'brain':'Brain'
        #'RedMarrow',
        #'lungs':'Lungs',
        #'Adipose'
        }
    ActivityMap = np.zeros((len(pbpk_para["FrameStartTimes"]),*seg_plus_body_arr.shape),dtype=np.float32) # 4D array (time, 3D)
    pixel_spacing_ml = np.prod(ct_get_zoom)*0.1**3 #multiplies 3 different spacing -> ml
    
    

    time, TACs = pycno.run_model(
        model_name="PSMA",
        stop=max(pbpk_para['FrameStartTimes']),
        observables=VOIs_possible
    )
    
    #TAC shape (1, 100, 18) (patient, time step, VOIs)
    log.debug("TAC created successfully")

    for key,value in class_seg.items(): # populate ActivityMap
        if key in roi_2_VOI:
            VOI = roi_2_VOI[key]
        else :
            VOI = 'Rest'
        VOI_index = VOIs_possible.index(VOI)
        
        mask_len_ROI = np.sum(masks[value])  # number of voxels in that ROI mask
        
        TAC_VOI = TACs[0, :, VOI_index]           # 1-D vector of the TAC for that VOI across the model grid, length 100
        frame_start = np.asarray(pbpk_para["FrameStartTimes"], float)  # 1-D vector of your frame start times length 3 [240, 1440]
        TAC_VOI_interp_time  = np.interp(frame_start, time, TAC_VOI)  #  1-D vector of interpolated activities at each frame time

        ActivityMap_Organ = np.zeros((len(pbpk_para["FrameStartTimes"]),*seg_plus_body_arr.shape),dtype=np.float32) # reset to zero for next organ
        ActivityMap_Organ[:, masks[value]] = TAC_VOI_interp_time[:, None]/(mask_len_ROI*pixel_spacing_ml)
        ActivityMap_Organ_path = os.path.join(out_paths['output_PBPK'], f'{pbpk_name}_{key}_act_av.bin')
        ActivityMap_Organ = ActivityMap_Organ.astype(np.float32)
        
        # #to be taken out 
        # i = 0.93
        # print(ActivityMap_Organ.shape)
        # x= math.floor((ActivityMap_Organ.shape[1]-ActivityMap_Organ.shape[1]*i)/2)
        # y=math.ceil((ActivityMap_Organ.shape[2]-ActivityMap_Organ.shape[2]*i)/2)
        # z=math.ceil((ActivityMap_Organ.shape[3]-ActivityMap_Organ.shape[3]*i)/2)
        # ActivityMap_Organ = ndimage.zoom(ActivityMap_Organ,(1,i,i,i),order=0)

        # ActivityMap_Organ = np.pad(ActivityMap_Organ,pad_width=((0,0),(x,x),(y,y),(z,z)))

        
        ActivityMap_Organ[0].tofile(ActivityMap_Organ_path) # saves only first frame for checking
        act_path_all_organ.append(ActivityMap_Organ_path)

        ActivityOrganSum[key] = np.sum(ActivityMap_Organ, axis=(1,2,3))*pixel_spacing_ml  # MBq
        ActivityMap[:, masks[value]] = ActivityMap_Organ[:, masks[value]]  # MBq/(vox*ml/vox) -> MBq/ml ### extra step to create full ActivityMap


        ### SAVING TAC FILES ###
        
        # --- save TAC time-series to .bin files (float32), no plotting ---
        tac_time_f    = os.path.join(out_paths['output_PBPK'], f"{pbpk_name}_{VOI}_TAC_time.bin")
        tac_values_f  = os.path.join(out_paths['output_PBPK'], f"{pbpk_name}_{VOI}_TAC_values.bin")
        samp_time_f   = os.path.join(out_paths['output_PBPK'], f"{pbpk_name}_{VOI}_sample_times.bin")
        samp_values_f = os.path.join(out_paths['output_PBPK'], f"{pbpk_name}_{VOI}_sample_values.bin")

        # Save: full model grid (time vs TAC)
        np.asarray(time, dtype=np.float32).tofile(tac_time_f)
        np.asarray(TACs[0, :, VOI_index], dtype=np.float32).tofile(tac_values_f)

        # Save: your frame sampling (frame_start vs interpolated TAC)
        np.asarray(frame_start, dtype=np.float32).tofile(samp_time_f)
        np.asarray(TAC_VOI_interp_time, dtype=np.float32).tofile(samp_values_f)

        log.debug(f"[PBPK] Saved TAC for {VOI}: "
                  f"{os.path.basename(tac_time_f)}, {os.path.basename(tac_values_f)}, "
                  f"{os.path.basename(samp_time_f)}, {os.path.basename(samp_values_f)}")
    
    
    ActivityMapSum = np.sum(ActivityMap, axis=(1,2,3))*pixel_spacing_ml #MBq
    
    log.debug(f"Activity Map Shape: {ActivityMap.shape}")
    log.debug(f"Activity Map Sum: {ActivityMapSum}")
    

    for i, frame in enumerate(ActivityMap):
        act_path_single = os.path.join(out_paths['output_PBPK'], f'{pbpk_name}_{pbpk_para["FrameStartTimes"][i]}_act_av.bin')
        act_path_all_map.append(act_path_single)
        frame = frame.astype(np.float32)
        frame.tofile(act_path_single)
    
    return ActivityMapSum, ActivityOrganSum, act_path_all_organ, act_path_all_map
