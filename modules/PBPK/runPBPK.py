import numpy as np
import pycno, os
import logging as log
# ensure a non-GUI backend so savefig works on headless servers
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def runPBPK(out_paths, pbpk_para,segmentated_ml_output_arr,masks,class_seg,ct_get_zoom):
    pbpk_name = pbpk_para['name']
    
    VOIs_possible = ['Tumor1', 'Tumor2', 'Kidney', 'Heart', 'SG', 'Bone', 'TumorRest', 'Spleen', 'Liver', 'Prostate', 'GI', 'Rest', 'Skin', 'Muscle', 'Brain', 'RedMarrow', 'Lungs', 'Adipose']
    roi_2_VOI = {
#        "NA":'Tumor1', 
#        "NA":'Tumor2', 
        "kidney":'Kidney', 
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
        "stomach"
        }

 
    time, TACs = pycno.run_model(
        model_name="PSMA",
        stop=max(pbpk_para['FrameStartTimes']),
        observables=VOIs_possible
    )
    #TAC shape (1, 100, 18) (patient, time step, VOIs)
    log.debug("TAC created successfully")

    del class_seg['Background']

    ActivityMap = np.zeros((len(pbpk_para["FrameStartTimes"]),*segmentated_ml_output_arr.shape),dtype=np.float32) # 4D array (time, 3D)
    print(class_seg)

    pixel_spacing_ml = np.prod(ct_get_zoom)*0.1**3 #multiplies 3 different spacing -> ml

    for key,value in class_seg.items(): # populate ActivityMap
        if key in roi_2_VOI:
            VOI = roi_2_VOI[key]
        else :
            VOI = 'Rest'
        VOI_index = VOIs_possible.index(VOI)
        
        mask_len_ROI = np.sum(masks[value])  # number of voxels in that ROI mask
        print(f"mask_len_ROI {key} : {mask_len_ROI}")
        print(pixel_spacing_ml)
        
        TAC_VOI = TACs[0, :, VOI_index]           # 1-D vector of the TAC for that VOI across the model grid, length 100
        frame_start = np.asarray(pbpk_para["FrameStartTimes"], float)  # 1-D vector of your frame start times length 3 [240, 1440]
        TAC_VOI_interp_time  = np.interp(frame_start, time, TAC_VOI)            #  1-D vector of interpolated activities at each frame time
        ActivityMap[:, masks[value]] = TAC_VOI_interp_time[:, None]/(mask_len_ROI*pixel_spacing_ml) # MBq/(vox*ml/vox) -> MBq/ml
        
        
        # create a figure, plot, save and close it
        fig, ax = plt.subplots()
        ax.plot(time, TACs[0, :, VOI_index], 'b-', label='TAC')
        
        # Plot interpolated points and vertical lines
        ax.plot(frame_start, TAC_VOI_interp_time, 'ro', label='Sample points')  # red dots at interpolated points
        for ft, interp_val in zip(frame_start, TAC_VOI_interp_time):
            ax.axvline(ft, color='red', linestyle='--', alpha=0.3, linewidth=0.9)
            ax.text(ft, interp_val, f't={int(ft)}s\n{interp_val:.2f}MBq', 
                   rotation=0, va='bottom', ha='right', color='red', fontsize=8)
        
        ax.set_title(f"TAC for {VOI}")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Activity (MBq)")
        ax.legend()
        
        fig.tight_layout()
        out_file = os.path.join(out_paths['output_PBPK'], f"{pbpk_name}_{VOI}_TAC.png")
        fig.savefig(out_file, bbox_inches='tight', dpi=150)
        plt.close(fig)
        
    
    log.info("Created Activity map using TACs and Segmentated Mask")
    print("Created Activity map using TACs and Segmentated Mask")
    
    
    ActivityMapSum = np.sum(ActivityMap, axis=(1,2,3))*pixel_spacing_ml #MBq
    
    log.debug(f"Activity Map Shape: {ActivityMap.shape}")
    log.debug(f"Activity Map Sum: {ActivityMapSum}")
    

    act_path_all = []
    for i, frame in enumerate(ActivityMap):
        act_path_single = os.path.join(out_paths['output_PBPK'], f'{pbpk_name}_{pbpk_para["FrameStartTimes"][i]}_act_av.bin')
        act_path_all.append(act_path_single)
        frame = frame.astype(np.float32)
        frame.tofile(act_path_single)


    return ActivityMapSum, act_path_all 
