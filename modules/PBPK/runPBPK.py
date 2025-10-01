import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pycno
import os
from skimage.transform import resize

def runPBPK(out_paths, pbpk_para, totseg_para, ct_input_arr, pixel_spacing_x,slice_thickness,resize_tuple,segmentated_ml_output_arr_resize,masks,class_seg):
    pbpk_name = pbpk_para['name']
    
    print("[PBPK] Beginning creating Time Activity Curves")
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

    time, TACs = pycno.run_model(model_name="PSMA", stop=pbpk_para['StopTime'])
    
    print("[PBPK] TAC created!")

    del class_seg['Background']
 

    for key,value in masks.items():
        value_resize = resize(value, resize_tuple)
        masks[key] = value_resize.astype(bool)
        
    #for now have all ones, need to turn back to zeros once know inside and outside body
    FrameArrs = np.ones((len(pbpk_para["FrameStartTimes"]),*segmentated_ml_output_arr_resize.shape),dtype=np.float32) # 4D array (time, 3D)

    for key,value in class_seg.items(): #{'Background': 0, 'liver': 5, 'pancreas': 7, 'esophagus': 15}
        if key in roi_2_VOI:
            VOI = roi_2_VOI[key]
        else :
            VOI = 'Rest'
        VOI_index = VOIs_possible.index(VOI)
        print(key,VOI,VOI_index, key, value)
        
        ActInterpolate = interp1d(time, TACs[0,:,VOI_index], kind='linear')
        for i, activity in enumerate(ActInterpolate(pbpk_para["FrameStartTimes"])):
            FrameArrs[i][masks[value]] = activity

    print("got the TAC correct!!")
    FrameAct = np.sum(FrameArrs, axis=(1,2,3))

    #turn to act_path
    print("[PBPK] Creating Activity Bin")
    act_path_all = []
    for i, frame in enumerate(FrameArrs):
        act_path_single = os.path.join(out_paths['output_PBPK'], f'{pbpk_name}_{pbpk_para["FrameStartTimes"][i]}_act_av.bin')
        act_path_all.append(act_path_single)
        frame = frame.astype(np.float32)
        frame.tofile(act_path_single)
    
    print(FrameAct,act_path_all)

    return FrameAct, act_path_all
