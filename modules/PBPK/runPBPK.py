import matlab.engine
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def runPBPK(config):
    print("Modeling PBPK...")
    eng = matlab.engine.start_matlab()
    eng.addpath('modules/matlab', nargout=0)

    eng.workspace['HotTotalAmount'] = config["PBPK"]["HotTotalAmount"]
    eng.workspace['ColdTotalAmount'] = config["PBPK"]["ColdTotalAmount"]
    eng.workspace['LambdaPhys'] = config["PBPK"]["LambdaPhys"]
    eng.workspace['StopTime'] = config["PBPK"]["StopTime"]

    eng.runPBPK(nargout=0)
    
    time = np.array(eng.workspace['time']).flatten()
    TACs = np.array(eng.workspace['TACs'])
    VOI_List = ["Tumor1", "Tumor2", "Background", "SG", "Liver", "Prostate", "Spleen", "RedMarrow", "Bone", "Muscle", "Heart", "Lungs", "Brain", "Skin", "GI", "Kidney", "Arteries", "Veins"]

    eng.quit()

    ImgAct = np.zeros((len(config["PBPK"]["FrameStartTimes"]), len(VOI_List)))

    for VOI in config["PBPK"]["VOIs"]:
        VOI_index = VOI_List.index(VOI)
        ActInterpolate = interp1d(time, TACs[:, VOI_index], kind='linear')

        ImgAct[:, VOI_index] = ActInterpolate(config["PBPK"]["FrameStartTimes"])

    FrameAct = np.sum(ImgAct, axis=1)
    
    for i in range(ImgAct.shape[0]):
        mask = ImgAct[i] == 0
        ImgAct[i, mask] = ImgAct[i, 2]
        ImgAct[i, 0] = 0

    return ImgAct, FrameAct

if __name__ == "__main__":
    runPBPK()