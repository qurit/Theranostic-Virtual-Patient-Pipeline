import nibabel as nib
import numpy as np
import os
import logging as log
from scipy.ndimage import zoom


def segmentation_2_class(nii_arr,classes,roi_subset):
    """
    The segmentation stores each ROI as a certain value
    
    The Total Seg Classes are saved as a Python file in dict form
    
    ie {'Background': 0, 'liver': 5, 'pancreas': 7, 'esophagus': 15}
    
    """
    log.debug(f"Running segmentation_2_class, ROI subset: {roi_subset}")
    log.debug(f"Unique Segmentation Values: {np.unique(nii_arr)}")

    class_seg = {}
    for n in np.unique(nii_arr):
        if n == 0:
            class_seg["Background"] = 0
        else:
            name = classes.get(n)
            class_seg[name] = int(n)
            
    return class_seg

def segmentation_multiple_arrays(seg_arr):
    """
    Build one binary mask per unique label in the segmentation volume.

    seg_arr            : 3D segmentation array (can be float or int).

    Returns: dict {bool : mask_array}, where mask_array has same shape as seg_arr.
    """
    arr = seg_arr
    labels = np.unique(arr)
    
    labels = labels[labels != 0] # gets rid of background 

    masks = {}
    for lab in labels:
        m = (arr == lab)
        masks[int(lab)] = m

    return masks

def hu_to_mu(CT,pixel_size, mu_water=0.1537,mu_bone=0.2234):
    """
    Convert a CT image from HU to μ (linear attenuation coefficient in cm^-1) using a bilinear transformation.
    Args: Values taken using XCAT phantom at 140keV
        CT (numpy array): Input tensor of HU values 
        mu_water (float): μ for water at 140.5 keV (default 0.1537 1/cm)
        mu_bone (float): μ for bone at 140.5 keV (default 0.2234 1/cm)
        pixel_size (float) : pixel size set in cm (cm/pixel)
        
    Return:
    mu_map: tensor of the CT in linear attentuion coefficant 

    """
    mu_water_pixel = mu_water*pixel_size #pixel/cm
    mu_bone_pixel = mu_bone*pixel_size #pixel cm
    
    log.debug(f"Attenuation of Water : {mu_water_pixel} pixel/cm")
    log.debug(f"Attenuation of Bone : {mu_bone_pixel} pixel/cm")
    
    mu_map = np.zeros_like(CT, dtype=np.float32)

    soft_tissue_mask = CT <= 0
    bone_mask = CT > 0
 
    mu_map[soft_tissue_mask] = mu_water_pixel * (1 + CT[soft_tissue_mask] / 1000)
    mu_map[bone_mask] = mu_water_pixel + (CT[bone_mask] / 1000) * (mu_bone_pixel - mu_water_pixel)
 
    return mu_map

def save_simind_mu_from_hu(hu_arr,segmentated_body_output_arr,out_dir,pixel_size,filename="TOTSEG_atn_av.bin"):
    """
    Convert HU -> mu (g/cc)
    then write a SIMIND voxel-based phantom file with 
    values = linear attn 

    hu_arr shape is assumed (X, Y, Z) when axis_order='XYZ'.
    """
    
    log.info("Converting HU values to Linear Attenutaion")
    print("[NII_PROCCESSING] Converting CT HU values to Linear Attenutaion...")
    hu = np.asarray(hu_arr, dtype=np.float32)
    
    mu_map = hu_to_mu(hu,pixel_size)
    mu_map = mu_map*segmentated_body_output_arr

    
    os.makedirs(out_dir, exist_ok=True)
    bin_path = os.path.join(out_dir, filename)
    mu_map.tofile(bin_path)

    return bin_path


def seg_ROI_plus_body(seg_roi_arr, seg_body_arr):
    # Ensure arrays are numpy arrays
    roi = np.asarray(seg_roi_arr)
    body = np.asarray(seg_body_arr)

    # Masks
    roi_mask  = roi != 0
    body_mask = body > 0

    # Start with zeros (background)
    seg_plus = np.zeros_like(roi, dtype=np.float32)

    # Body-only voxels labeled in totalseg map
    seg_plus[body_mask] = 201

    seg_plus[roi_mask] = roi[roi_mask]

    # Union mask 
    mask_body_plus_roi = segmentation_multiple_arrays(seg_plus)
    
    return seg_plus, mask_body_plus_roi
    



def NII_PROCCESSING(output_path,classes,simind_para,totseg_para,ml_file,body_file):
    """
    Processing Nifti File of segmentated CT Scan. 

    asigns 0 to anything outside of body
    """
    totseg_name = totseg_para['name']
    roi_subset = totseg_para['roi_subset']
    seg_ml_file = ml_file
    seg_body_file = f"{body_file}/body.nii.gz"
    
    
    ############### LOAD IMAGE ###############
    
    ct_input = nib.load(f"{output_path}/{totseg_name}_ct_input.nii.gz")
    segmentated_ml_output  = nib.load(f"{seg_ml_file}")
    segmentated_body_output = nib.load(f"{seg_body_file}")

    log.debug("CT, Segmentated ROI's and Segmentated Body successfully loaded into nii_processing.py")
    
    ############### CONVERT TO float32 ARRAYS + REORIENT  ###############
    
    print("[NII_PROCCESSING] Obtaining data from Nifti files (storing as numpy array) ")

    

    ct_input_arr = np.transpose(np.array(ct_input.get_fdata(dtype=np.float32)),(2,1,0))[:,::-1,:]
    segmentated_ml_output_arr  = np.transpose(np.array(segmentated_ml_output.get_fdata(dtype=np.float32)),(2,1,0))[:,::-1,:]
    segmentated_body_output_arr  = np.transpose(np.array(segmentated_body_output.get_fdata(dtype=np.float32)),(2,1,0))[:,::-1,:]
    
    log.debug("CT, Segmentated ROI's and Segmentated Body successfully reorientated and converted to float32 arrays")

    
    log.debug(f"CT input matrix has shape: {ct_input_arr.shape}")
    log.debug(f"Segmentated ROIs matrix has shape: {segmentated_ml_output_arr.shape}")

    if "kidney_left" in roi_subset and "kidney_right" in roi_subset: # converts left and right kidney into one
        log.debug("Turning kidney left and kidney right into one set : kidney")
        kidney_mask = (segmentated_ml_output_arr == 2) | (segmentated_ml_output_arr == 3)
        segmentated_ml_output_arr[kidney_mask] = 200
    elif ("kidney_left" in roi_subset and 'kidney_right' not in roi_subset) or ("kidney_right" in roi_subset and 'kidney_right' not in roi_subset):
        print("Sorry, incapable of running program with only one of the kidney (left or right) currently. Please segmentate both!")
        log.critical("Sorry, incapable of running program with only one of the kidney (left or right) currently. Please segmentate both!")
        quit()

    ############### RESIZE ###############
    resize = simind_para['resize']
    scale_factor = resize / ct_input_arr.shape[1]  # assuming square in x and y dimensions
    ct_input_arr = zoom(ct_input_arr, (scale_factor, scale_factor, scale_factor), order=0) # didnt resize z, should make functional so in config
    segmentated_ml_output_arr = zoom(segmentated_ml_output_arr, (scale_factor, scale_factor, scale_factor), order=0)
    segmentated_body_output_arr = zoom(segmentated_body_output_arr, (scale_factor, scale_factor, scale_factor), order=0)

    log.info(f"Resized CT and Segmentated arrays to {segmentated_ml_output_arr.shape} for SIMIND input")
    print(f"[NII_PROCCESSING] Resized CT and Segmentated arrays to {segmentated_ml_output_arr.shape} for SIMIND input")
    
    ############### CLASS SEG. & MASK SEG. ###############
    seg_plus_body_arr,mask_roi_plus_body = seg_ROI_plus_body(segmentated_ml_output_arr,segmentated_body_output_arr)
    
    mask_roi = segmentation_multiple_arrays(segmentated_ml_output_arr)
    
    class_seg = segmentation_2_class(seg_plus_body_arr,classes,roi_subset)

    ############### FOR SIMIND: DENSITY / ATTENUTATION MAP ###############
    
    #Linear Atten Path(will save atn_av as a bin float32 file )
    ct_get_zoom = tuple(np.array(ct_input.header.get_zooms()) / scale_factor) #[x,y,z] mm
    pixel_spacing_cm = ct_get_zoom[0] * 0.1 #cm
    slice_thickness = ct_get_zoom[2] * 0.1   # Slice thickness (indxK-direction spacing) (cm)

    log.debug(f"Pixel Size : {pixel_spacing_cm} cm")
    log.debug(f"Slice Thickness  : {slice_thickness} cm")
    
    log.info("Saving CT scan as a Linear Attenutation matrix from a HU matrix.")
    print("[NII_PROCCESSING] Saving CT scan as a Linear Attenutation matrix from a HU matrix...")
    atn_av_path = save_simind_mu_from_hu(ct_input_arr, segmentated_body_output_arr, output_path,pixel_spacing_cm)

    ############### SAVE SEG ROI & BODY AS BINARY ###############
    seg_ml_bin_path = os.path.join(output_path, f"{totseg_name}_ml_segmentation.bin")
    segmentated_ml_output_arr.astype(np.float32).tofile(seg_ml_bin_path)
    log.info(f"Saved Segmentated ROI's as binary file at: {seg_ml_bin_path}")
    print(f"[NII_PROCCESSING] Saved Segmentated ROI's as binary file at: {seg_ml_bin_path}")    
    
    seg_body_bin_path = os.path.join(output_path, f"{totseg_name}_body_segmentation.bin")
    segmentated_body_output_arr.astype(np.float32).tofile(seg_body_bin_path)
    log.info(f"Saved Segmentated Body as binary file at: {seg_body_bin_path}")
    print(f"[NII_PROCCESSING] Saved Segmentated Body as binary file at: {seg_body_bin_path}")
    
    seg_body_plus_ml_bin_path = os.path.join(output_path, f"{totseg_name}_body_plus_ml_segmentation.bin")
    seg_plus_body_arr.astype(np.float32).tofile(seg_body_plus_ml_bin_path)
    log.info(f"Saved Segmentated Body as binary file at: {seg_body_plus_ml_bin_path}")
    print(f"[NII_PROCCESSING] Saved Segmentated Body as binary file at: {seg_body_plus_ml_bin_path}")    
    
    
    return ct_input_arr,segmentated_ml_output_arr,segmentated_body_output_arr,seg_plus_body_arr,class_seg, mask_roi, mask_roi_plus_body,atn_av_path , seg_ml_bin_path, seg_body_bin_path, pixel_spacing_cm, slice_thickness,ct_get_zoom 
