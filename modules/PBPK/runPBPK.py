import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pycno
import os

def runPBPK(out_paths, pbpk_para, totseg_para, ct_input_arr, pixel_spacing_x,slice_thickness,segmentated_ml_output_arr,masks,class_seg):
    
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
        "prostate":'Prostate'
        #'GI', 
        #'Rest',
        #'Skin',
        #'Muscle',
        #'Brain',
        #'RedMarrow',
        #'Lungs',
        #'Adipose'
        }

    roi_list = totseg_para['roi_subset']
    

    time, TACs = pycno.run_model(model_name="PSMA", stop=pbpk_para['StopTime'])
    
    print("[PBPK] TAC created!")

    del class_seg['Background']
    #print(class_seg)
    
    FrameArrs = np.zeros((len(pbpk_para["FrameStartTimes"]),*segmentated_ml_output_arr.shape),dtype=np.float32) # 4D array (time, 3D)
    #print(FrameArrs.shape)

    for key,value in class_seg.items(): #{'Background': 0, 'liver': 5, 'pancreas': 7, 'esophagus': 15}
        if key in roi_2_VOI:
            VOI = roi_2_VOI[key]
        else :
            VOI = 'Rest'
        VOI_index = VOIs_possible.index(VOI)
        print(key,VOI,VOI_index)
        
        ActInterpolate = interp1d(time, TACs[0,:,VOI_index], kind='linear')
        
        for i, activity in enumerate(ActInterpolate(pbpk_para["FrameStartTimes"])):
            print(1)
            FrameArrs[i][masks[value]] = activity

        #FrameArrs[:, masks[value]] = ActInterpolate(pbpk_para["FrameStartTimes"])

    print("got the TAC correct!!")

    FrameSum = np.sum(FrameArrs, axis=0)
    
    '''
    for i in range(ImgAct.shape[0]):
        mask = ImgAct[i] == 0
        ImgAct[i, mask] = ImgAct[i, 2]
        ImgAct[i, 0] = 0
    '''
    #turn to act_path
    
    print("[PBPK] Creating Activity Bin")
    
    """
    ROICounts = np.zeros(np.size(ImgAct, 1))
    for i in range(np.size(ImgAct, 1)):
        ROICounts[i] = np.sum(ROIMap == i)
    ROICounts[ImgAct[0] == ImgAct[0][2]] = np.sum(ROICounts[ImgAct[0] == ImgAct[0][2]])
    ROIVolumes = ROICounts * pixel_spacing_x**2 * slice_thickness

    with np.errstate(divide='ignore', invalid='ignore'):
        ImgAct = np.divide(ImgAct,ROIVolumes[np.newaxis, :])
    ImgAct[np.isnan(ImgAct)] = 0
    ImgAct[np.isinf(ImgAct)] = 0
    """

    
    act_path_all = []
    for i, frame in enumerate(FrameArrs):
        act_path_single = os.path.join(out_paths['output_PBPK'], f'PBPK_frame_{pbpk_para["FrameStartTimes"][i]}_act_av.bin')
        act_path_all.append(act_path_single)
        
        frame = frame.astype(np.float32)
        frame.tofile(act_path_single)
    
    print(FrameSum,act_path_all)

    return FrameSum, act_path_all
