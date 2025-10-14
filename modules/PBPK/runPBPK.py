import numpy as np
import matplotlib.pyplot as plt
import pycno, os
import logging as log

def runPBPK(out_paths, pbpk_para,segmentated_ml_output_arr,masks,class_seg):
    pbpk_name = pbpk_para['name']
    
    VOIs_possible = ['Tumor1', 'Tumor2', 'Kidney', 'Heart', 'SG', 'Bone', 'TumorRest', 'Spleen', 'Liver', 'Prostate', 'GI', 'Rest', 'Skin', 'Muscle', 'Brain', 'RedMarrow', 'Lungs', 'Adipose']
    roi_2_VOI = {
#        "NA":'Tumor1', 
#        "NA":'Tumor2', 
        "kidney_right":'Kidney', 
        "kidney_left":'Kidney',
        "heart":'Heart', 
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
        #'Lungs',
        #'Adipose'
        }

    time, TACs = pycno.run_model(model_name="PSMA", stop=pbpk_para['StopTime']) #TAC shape (1, 100, 18) (patient, time step, VOIs)
    log.debug("TAC created successfully")

    del class_seg['Background']

    ActivityMap = np.zeros((len(pbpk_para["FrameStartTimes"]),*segmentated_ml_output_arr.shape),dtype=np.float32) # 4D array (time, 3D)

    for key,value in class_seg.items(): # populate ActivityMap 
        if key in roi_2_VOI:
            VOI = roi_2_VOI[key]
        else :
            VOI = 'Rest'
        VOI_index = VOIs_possible.index(VOI)
        
        plt.plot(time,TACs[0,:,VOI_index],label = VOI)

        
        TAC_VOI = TACs[0, :, VOI_index]                           # 1-D vector of the TAC for that VOI across the model grid, length 100
        frame_start = np.asarray(pbpk_para["FrameStartTimes"], float)  # 1-D vector of your frame start times length 3
        TAC_VOI_interp_time  = np.interp(frame_start, time, TAC_VOI)                         #  1-D vector of interpolated activities at each frame time
        ActivityMap[:, masks[value]] = TAC_VOI_interp_time[:, None]
    
    log.info("Created Activity map using TACs and Segmentated Mask")
    print("Created Activity map using TACs and Segmentated Mask")
    
    
    ActivityMapSum = np.sum(ActivityMap, axis=(1,2,3))
    
    log.debug(f"Activity Map Shape: {ActivityMap.shape}")
    log.debug(f"Activity Map Sum: {ActivityMapSum}")
    

    act_path_all = []
    for i, frame in enumerate(ActivityMap):
        act_path_single = os.path.join(out_paths['output_PBPK'], f'{pbpk_name}_{pbpk_para["FrameStartTimes"][i]}_act_av.bin')
        act_path_all.append(act_path_single)
        frame = frame.astype(np.float32)
        frame.tofile(act_path_single)


    return ActivityMapSum, act_path_all
