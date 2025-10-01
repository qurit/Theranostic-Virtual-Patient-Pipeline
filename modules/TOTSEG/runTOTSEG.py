#Comment out if not using Mac :
# -------------------- Environment & worker limits (set BEFORE heavy imports) --------------------
import os, multiprocessing as mp

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

import torch
from torch.cuda.amp import GradScaler as _GradScaler
torch.GradScaler = _GradScaler
from totalsegmentator.python_api import totalsegmentator
import SimpleITK as sitk
import dicom2nifti
import datetime




def runTOTSEG(input_path,output_path, totseg_para):
    """
    Run TotalSegmentator with common convenience options.

    Notes:
      - If you run on CPU, use --fast and/or --roi_subset for big speedups.
    """
    out_dir = output_path
    os.makedirs(out_dir, exist_ok=True)
    
    totseg_name = totseg_para['name']
    device = totseg_para['device']
    fast = totseg_para['fast']
    roi_subset = totseg_para['roi_subset']          
    ml = totseg_para['ml']
    statistics = totseg_para['statistics']
    radiomics = totseg_para['radiomics']

    print("[TOTSEG]=== TotalSegmentator Run ===")
    print(f"[TOTSEG] Input:   {input_path}")
    print(f"[TOTSEG] Output:  {out_dir}")
    print(f"[TOTSEG] Device:  {device}")
    print(f"[TOTSEG] fast:  {fast}")
    print(f"[TOTSEG] roi_subset:  {roi_subset}")
    print(f"[TOTSEG] ml:  {ml}")
    print(f"[TOTSEG] statistics:  {statistics}")
    print(f"[TOTSEG] radiomics:  {radiomics}")
    
    

    # If radiomics is on with ml, force ml=False (TS limitation).
    if radiomics==True and ml==True:
        print("[TOTSEG] Radiomics not supported with --ml; forcing ml=False so radiomics can run.")
        ml = False
    
    # ---- Guarantee a NIfTI in out_dir (ct_input.nii.gz) ----
    ct_nifti = os.path.join(out_dir, f"{totseg_name}_ct_input.nii.gz")

    if os.path.isdir(input_path):  # DICOM folder
        print("[TOTSEG] DICOM folder -> converting to NIfTI (ct_input.nii.gz)")

        dicom2nifti.dicom_series_to_nifti(input_path, ct_nifti, reorient_nifti=True)
    else:  # NIfTI file (.nii or .nii.gz) -> write a copy to output
        low = input_path.lower()
        if low.endswith((".nii", ".nii.gz")):
            print("[TOTSEG] NIfTI -> copying to ct_input.nii.gz")
            sitk.WriteImage(sitk.ReadImage(input_path), ct_nifti, True)
        else:
            raise ValueError("[TOTSEG] Unsupported input. Provide a DICOM folder or a NIfTI (.nii/.nii.gz).")

    print(f"[TOTSEG] Using NIfTI: {ct_nifti}")

    kwargs = {
        "device": device,
        # "task": task,  # keep commented unless you explicitly need non-default task (e.g., "total_mr")
    }
    if fast: 
        kwargs["fast"] = True
        print("[TOTSEG] Mode:    FAST (3mm model)")
    if roi_subset: 
        kwargs["roi_subset"] = roi_subset 
        print(f"[TOTSEG] ROI subset: {roi_subset}")
    if ml: 
        kwargs["ml"] = True
        print("[TOTSEG] ml mode ON")
    if statistics: 
        kwargs["statistics"] = True
        print("[TOTSEG] stats mode ON")
    if radiomics: 
        kwargs["radiomics"] = True
        print("[TOTSEG] radiomics mode ON")

    # Run TS
    if ml:
        # write a single multilabel NIfTI inside your output folder
        ml_file = os.path.join(out_dir, f"{totseg_name}_ml_segmentation.nii.gz") 
        totalsegmentator(ct_nifti, ml_file, **kwargs)
        print("[TOTSEG] Segmentation complete (multilabel).")
        print(f"[TOTSEG] Multilabel file: {ml_file}")
    else:
        # write the usual per-class directory structure
        totalsegmentator(ct_nifti, out_dir, **kwargs)
        print("[TOTSEG] Segmentation complete.")

 