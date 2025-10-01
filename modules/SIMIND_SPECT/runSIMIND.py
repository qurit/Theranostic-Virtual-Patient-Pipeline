import os
import subprocess
import shutil
import numpy as np

def runSIMIND(path,simind_para, pbpk_para, output_path, num_cores, InputImageSize, ImgLength, ImgSum,pixel_spacing_x, slice_thickness,act_path_all,atn_path,index_14=-7,index_15=-7):
    print("[SIMIND] Beginning SIMIND Simulation")
    output_name = simind_para['name']
    os.environ['SMC_DIR'] = os.path.join(simind_para["SIMINDDirectory"], 'smc_dir/')
    os.environ['PATH'] = simind_para["SIMINDDirectory"] + ':' + os.environ.get('PATH', '') 
    InputImgSize = InputImageSize
    Collimator = simind_para["Collimator"]
    Isotope = simind_para["Isotope"]
    NumProjections = simind_para["NumProjections"]
    DetectorDistance = simind_para["DetectorDistance"]
    OutputImgSize = simind_para["OutputImgSize"]
    PixelWidth = pixel_spacing_x
    SliceWidth = slice_thickness
    OutputPixelSize = simind_para["OutputPixelSize"]
    OutputSliceWidth = simind_para["OutputSliceWidth"]
    OutputImgLength = SliceWidth*ImgLength/OutputSliceWidth 
    NumPhotons = simind_para["NumPhotons"]
    index_32 = simind_para["Index_32"]
    
    print("[SIMIND] Inital Parameters Stored")
    
    ScatwinFile = os.path.join(path, 'bin/scattwin.win')
    ScatwinFileOut = os.path.join(output_path, f'{output_name}.win')
    shutil.copyfile(ScatwinFile, ScatwinFileOut)

    SmcFile = os.path.join(path, 'bin/smc.smc')
    SmcFileOut = os.path.join(output_path, f'{output_name}.smc')
    shutil.copyfile(SmcFile, SmcFileOut)
    
    print("[SIMIND] Scat and Smc files stored")

    HalfLength = SliceWidth*ImgLength/2

    for index,value in enumerate((pbpk_para["FrameStartTimes"])):

        ScaleFactor = NumPhotons/ImgSum[index]/num_cores

        simind_exe = os.path.join(simind_para["SIMINDDirectory"], 'simind')
        
        atn_name = os.path.basename(atn_path)                    # e.g., "test_atn_av.bin"
        act_name = os.path.basename(act_path_all[index])         # e.g., "PBPK_frame_100_act_av.bin"

        shutil.copyfile(atn_path, os.path.join(output_path, atn_name))
        shutil.copyfile(act_path_all[index], os.path.join(output_path, act_name))
        
        print(f"[SIMIND] saved in SIMIND dir, inputed atn file :{atn_name}")
        print(f"[SIMIND] saved in SIMIND dir, inputed act file :{act_name}")
        simind_switches = (
            f'/fd:{atn_name}'              
            f'/fs:{act_name}'              
            f'/14:{index_14}'             # -7
            f'/15:{index_15}'             # -7
            f'/cc:{Collimator}'
            f'/fi:{Isotope}'
            f'/nn:{ScaleFactor}'
            f'/in:x22,6x'               
            f'/02:{HalfLength}'
            f'/05:{HalfLength}'
            f'/12:{DetectorDistance}'
            f'/28:{OutputPixelSize}'
            f'/29:{NumProjections}'
            f'/31:{PixelWidth}'
            f'/32:{index_32}'
            f'/34:{ImgLength}'
            f'/76:{OutputImgSize}'
            f'/77:{OutputImgLength}'
            f'/78:{InputImgSize}'
            f'/79:{InputImgSize}'
        )

        processes = []
        for j in range(num_cores):
            simind_command = f'{simind_exe} {output_name} {output_name}_frame{index}_{j} '
            command = simind_command + simind_switches
            if j == 0:
                p = subprocess.Popen(command, shell=True, cwd=output_path)
            else:
                p = subprocess.Popen(command, shell=True, cwd=output_path, stdout=subprocess.DEVNULL)
            processes.append(p)

        for p in processes:
            p.wait()
        

        xtot_w1 = 0
        xtot_w2 = 0
        xtot_w3 = 0
        for j in range(0, num_cores):
            w1 = np.fromfile(os.path.join(output_path, f'{output_name}_frame{index}_{j}_tot_w1.a00'), dtype=np.float32)
            w2 = np.fromfile(os.path.join(output_path, f'{output_name}_frame{index}_{j}_tot_w2.a00'), dtype=np.float32)
            w3 = np.fromfile(os.path.join(output_path, f'{output_name}_frame{index}_{j}_tot_w3.a00'), dtype=np.float32)
            xtot_w1+=w1
            xtot_w2+=w2
            xtot_w3+=w3
        xtot_w1 = xtot_w1 / num_cores
        xtot_w2 = xtot_w2 / num_cores
        xtot_w3 = xtot_w3 / num_cores

        xtot_w1.tofile(os.path.join(output_path, f'{output_name}_frame{index}_tot_w1.a00'))
        xtot_w2.tofile(os.path.join(output_path, f'{output_name}_frame{index}_tot_w2.a00'))
        xtot_w3.tofile(os.path.join(output_path, f'{output_name}_frame{index}_tot_w3.a00'))

        for h in range(1, 4):
            HeaderFile = os.path.join(output_path, f'{output_name}_frame{index}_0_tot_w{h}.h00')
            HeaderFileOut = os.path.join(output_path, f'{output_name}_frame{index}_tot_w{h}.h00')
            shutil.copyfile(HeaderFile, HeaderFileOut)

        hctFile = os.path.join(output_path, f'{output_name}_frame{index}_0.hct')
        hctFileOut = os.path.join(output_path, f'{output_name}_frame{index}.hct')
        shutil.copyfile(hctFile, hctFileOut)

        ictFile = os.path.join(output_path, f'{output_name}_frame{index}_0.ict')
        ictFileOut = os.path.join(output_path, f'{output_name}_frame{index}.ict')
        shutil.copyfile(ictFile, ictFileOut)

    JaszakFile = os.path.join(path, 'bin/jaszak.smc')
    JaszakFileOut = os.path.join(output_path, f'jaszak.smc')
    shutil.copyfile(JaszakFile, JaszakFileOut)

    calibration_command = f'{simind_exe} jaszak calib/fi:{Isotope}/cc:{Collimator}/29:1/15:5/fa:11/fa:15/fa:14'
    subprocess.run(calibration_command, shell=True, cwd=output_path, stdout=subprocess.DEVNULL)
    print("[SIMIND] Files Saved")
    return
