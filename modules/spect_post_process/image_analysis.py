

import numpy as np


def tac_calculation_spect(spect_image_path,img_shape,mask,reshape_ratio):
    mask_resize = mask * reshape_ratio
    arr = np.fromfile(spect_image_path, dtype=np.float32).reshape(img_shape)
    
    