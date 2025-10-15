import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter
import os
import logging as log


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

def hu_to_rho(nii_arr): # not used
    """
    note that the nii.arrays have indicies referring to coord. 
    CT Hounsfield unit (HU) values quantify X-ray attenuation (density) 
    with water set at 0 HU, air at -1000 HU, 
    and other tissues having positive or negative values accordingly
    
    HU isn’t strictly −1000…+1000. 
    Air is ~−1000 HU and water is 0, 
    but dense bone is often >+1000, and metal can be 2000–3000+ HU. 
    Many CTs therefore span roughly −1024 to ~3071 
    (that’s the 12-bit stored range 0…4095 with an intercept of −1024 applied).
    The value −2048 usually means “padding/outside FOV/invalid” 
    introduced by the scanner or by the DICOM→NIfTI converter 
    (often used as a sentinel outside the body). It’s not a real tissue value.
    
    HU = 1000*(u-u_water)/u_water 
    u = linear attenuation coefficient of material
    u_water = linear attenuation coefficient of  water
    
    Attenuation depends on both density and composition 
    (via the mass attenuation coefficient μ/ρ)
    and on the effective X-ray energy of your scan
    
    Conversion based on Schneider et al. 2000 (using GATE's material db example)
    Convert a CT array, in HU into a density map in g/cc
    """
    # Define the bin edges for HU values
    bins = np.array(
        [
            -1050,
            -950,
            -852.884,
            -755.769,
            -658.653,
            -561.538,
            -464.422,
            -367.306,
            -270.191,
            -173.075,
            -120,
            -82,
            -52,
            -22,
            8,
            19,
            80,
            120,
            200,
            300,
            400,
            500,
            600,
            700,
            800,
            900,
            1000,
            1100,
            1200,
            1300,
            1400,
            1500,
            1640,
            1807.5,
            1975.01,
            2142.51,
            2300,
            2467.5,
            2635.01,
            2802.51,
            2970.02,
            3000,
        ]
    )

    # Define the corresponding density values for each bin
    values = np.array(
        [
            0.00121,
            0.102695,
            0.202695,
            0.302695,
            0.402695,
            0.502695,
            0.602695,
            0.702695,
            0.802695,
            0.880021,
            0.926911,
            0.957382,
            0.984277,
            1.01117,
            1.02955,
            1.0616,
            1.1199,
            1.11115,
            1.16447,
            1.22371,
            1.28295,
            1.34219,
            1.40142,
            1.46066,
            1.5199,
            1.57914,
            1.63838,
            1.69762,
            1.75686,
            1.8161,
            1.87534,
            1.94643,
            2.03808,
            2.13808,
            2.23808,
            2.33509,
            2.4321,
            2.5321,
            2.6321,
            2.7321,
            2.79105,
            2.9,
        ]
    )

    # Clip the HU array values to be within the range of defined bins
    hu_clipped = np.clip(nii_arr, bins[0], bins[-1])

    # Apply Median filter to remove a bit of remaining noise
    hu_clipped = median_filter(hu_clipped, size=2)

    # Find the corresponding bin for each HU value
    bin_indices = np.digitize(hu_clipped, bins, right=True)

    # Map each bin index to the corresponding density value
    rho = values[bin_indices - 1]

    return rho

def save_simind_density_bin_from_hu(hu_arr,out_dir,filename="d1000.dmi"): #not used
    """
    Convert HU -> density (g/cc) ( assume in XYZ order)
    then write a SIMIND voxel-based phantom file with 
    values = density * 1000 as 16-bit integers.

    """
    hu = np.asarray(hu_arr, dtype=np.float32)
    hu = np.where(hu <= -2000.0, -1000.0, hu)  # padding -> air


    print("[NII_PROCCESSING] Converting HU values to mass density")
    rho = hu_to_rho(hu).astype(np.float32)     
    scaled = rho * 1000.0                 # g/cc -> “density×1000”
    scaled = np.clip(scaled, 0, 30000)    # keep in [0, 30000] to avoid nonsense/overflow
    scaled = np.rint(scaled)              # unbiased rounding to nearest int
    d1000  = scaled.astype(np.int16)      # compact, SIMIND-friendly type   
    
    os.makedirs(out_dir, exist_ok=True)
    bin_path = os.path.join(out_dir, filename)
    d1000.tofile(bin_path)
    print("[NII_PROCCESSING] Mass density scaled and stored as 16 bit int(dmi file)")

    return str(bin_path)

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
=
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
    
     ############### CLASS SEG. & MASK SEG. ###############
    class_seg = segmentation_2_class(segmentated_ml_output_arr,classes,roi_subset)
    
    masks = segmentation_multiple_arrays(segmentated_ml_output_arr)


    ############### FOR SIMIND: DENSITY / ATTENUTATION MAP ###############
    
    '''
    #Density Path
    #will save d1000 ( density*1000 in 16 bit) bin file
    #hu_arr shape is assumed (X, Y, Z) when axis_order='XYZ'.
    #Set Index-32 in your .smc to match `index32_orientation`
    
    print("[NII_PROCCESSING] Saving CT scan as a density matrix from a HU matrix.")
    dmi_file_path = save_simind_density_bin_from_hu(ct_input_arr, f"{output_path}/{output_file_name}")
    print(f"[NII_PROCCESSING] Saved: {dmi_file_path}")
    '''
    
    #Linear Atten Path(will save atn_av as a bin float32 file )
    ct_get_zoom = ct_input.header.get_zooms() #[x,y,z] mm 
    pixel_spacing_cm = ct_get_zoom[0] * 0.1 #cm
    slice_thickness = ct_get_zoom[2] * 0.1  # Slice thickness (indxK-direction spacing) (cm)

    log.debug(f"Pixel Size :{pixel_spacing_cm} cm")
    log.debug(f"Slice Thickness  :{slice_thickness} cm")
    
    log.info("Saving CT scan as a Linear Attenutation matrix from a HU matrix.")
    print("[NII_PROCCESSING] Saving CT scan as a Linear Attenutation matrix from a HU matrix...")
    atn_av_path = save_simind_mu_from_hu(ct_input_arr, segmentated_body_output_arr, output_path,pixel_spacing_cm)

    return ct_input_arr,segmentated_ml_output_arr,segmentated_body_output_arr, class_seg, masks, atn_av_path , pixel_spacing_cm, slice_thickness,ct_get_zoom
