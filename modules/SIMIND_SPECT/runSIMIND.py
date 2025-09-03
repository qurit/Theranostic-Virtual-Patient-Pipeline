import os
import subprocess
import shutil
import numpy as np
       
def runSIMIND(path, config, output_name, output_path, num_cores, ImgLength, ImgSum,act_path,atn_path):
    os.environ['SMC_DIR'] = os.path.join(config["SIMIND"]["SIMINDDirectory"], 'smc_dir/')
    os.environ['PATH'] = os.path.join(config["SIMIND"]["SIMINDDirectory"], ':$PATH') + ':' + os.environ['PATH']
    InputImgSize = config["XCAT"]["InputImgSize"]
    Collimator = config["SIMIND"]["Collimator"]
    Isotope = config["SIMIND"]["Isotope"]
    NumProjections = config["SIMIND"]["NumProjections"]
    DetectorDistance = config["SIMIND"]["DetectorDistance"]
    OutputImgSize = config["SIMIND"]["OutputImgSize"]
    PixelWidth = config["XCAT"]["PixelWidth"]
    SliceWidth = config["XCAT"]["SliceWidth"]
    OutputPixelSize = config["SIMIND"]["OutputPixelSize"]
    OutputSliceWidth = config["SIMIND"]["OutputSliceWidth"]
    OutputImgLength = SliceWidth*ImgLength/OutputSliceWidth # what is OutputSliceWidth purpose??
    NumPhotons = config["SIMIND"]["NumPhotons"]

    ScatwinFile = os.path.join(path, 'bin/scattwin.win')
    ScatwinFileOut = os.path.join(output_path, f'{output_name}.win')
    shutil.copyfile(ScatwinFile, ScatwinFileOut)

    SmcFile = os.path.join(path, 'bin/smc.smc')
    SmcFileOut = os.path.join(output_path, f'{output_name}.smc')
    shutil.copyfile(SmcFile, SmcFileOut)

    HalfLength = SliceWidth*ImgLength/2

    for frame in range(len(config["PBPK"]["FrameStartTimes"])):

        ScaleFactor = NumPhotons/ImgSum[frame]/num_cores

        simind_exe = os.path.join(config["SIMIND"]["SIMINDDirectory"], 'simind')

        simind_switches = (
            f'/fd:{atn_path}' 
            f'/fs:{act_path}'
            f'/cc:{Collimator}'
            f'/fi:{Isotope}'
            f'/nn:{ScaleFactor}'
            f'/in:x22,3x' #set output mu map
            f'/02:{HalfLength}'
            f'/05:{HalfLength}'
            f'/12:{DetectorDistance}'
            f'/28:{OutputPixelSize}'
            f'/29:{NumProjections}'
            f'/31:{PixelWidth}' #density pixel size
            f'/34:{ImgLength}' #num density images
            f'/76:{OutputImgSize}'
            f'/77:{OutputImgLength}'
            f'/78:{InputImgSize}' #density image size
            f'/79:{InputImgSize}' #source image size
        )

        processes = []
        for j in range(num_cores):
            simind_command = f'{simind_exe} {output_name} {output_name}_frame{frame}_{j}'
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
            w1 = np.fromfile(os.path.join(output_path, f'{output_name}_frame{frame}_{j}_tot_w1.a00'), dtype=np.float32)
            w2 = np.fromfile(os.path.join(output_path, f'{output_name}_frame{frame}_{j}_tot_w2.a00'), dtype=np.float32)
            w3 = np.fromfile(os.path.join(output_path, f'{output_name}_frame{frame}_{j}_tot_w3.a00'), dtype=np.float32)
            xtot_w1+=w1
            xtot_w2+=w2
            xtot_w3+=w3
        xtot_w1 = xtot_w1 / num_cores
        xtot_w2 = xtot_w2 / num_cores
        xtot_w3 = xtot_w3 / num_cores

        xtot_w1.tofile(os.path.join(output_path, f'{output_name}_frame{frame}_tot_w1.a00'))
        xtot_w2.tofile(os.path.join(output_path, f'{output_name}_frame{frame}_tot_w2.a00'))
        xtot_w3.tofile(os.path.join(output_path, f'{output_name}_frame{frame}_tot_w3.a00'))

        for h in range(1, 4):
            HeaderFile = os.path.join(output_path, f'{output_name}_frame{frame}_0_tot_w{h}.h00')
            HeaderFileOut = os.path.join(output_path, f'{output_name}_frame{frame}_tot_w{h}.h00')
            shutil.copyfile(HeaderFile, HeaderFileOut)

        hctFile = os.path.join(output_path, f'{output_name}_frame{frame}_0.hct')
        hctFileOut = os.path.join(output_path, f'{output_name}_frame{frame}.hct')
        shutil.copyfile(hctFile, hctFileOut)

        ictFile = os.path.join(output_path, f'{output_name}_frame{frame}_0.ict')
        ictFileOut = os.path.join(output_path, f'{output_name}_frame{frame}.ict')
        shutil.copyfile(ictFile, ictFileOut)

    JaszakFile = os.path.join(path, 'bin/jaszak.smc')
    JaszakFileOut = os.path.join(output_path, f'jaszak.smc')
    shutil.copyfile(JaszakFile, JaszakFileOut)

    calibration_command = f'{simind_exe} jaszak calib/fi:{Isotope}/cc:{Collimator}/29:1/15:5/fa:11/fa:15/fa:14'
    subprocess.run(calibration_command, shell=True, cwd=output_path, stdout=subprocess.DEVNULL)

    return