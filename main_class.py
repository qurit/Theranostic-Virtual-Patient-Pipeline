"""
Theranostic Digital Twin (TDT) main entry point.

This script incoperates the various functionalites for the Theranostic Digital Twin
(TDT) pipeline. 

The various capabilties of this pipeline include:
    1. Simulating SPECT images from real CT images via organ segmentation,
         physiologically based pharmacokinetic (PBPK) modelling, Monte Carlo SPECT
         simulation (SIMIND), and image reconstruction (PyTomography).
    2. (Future) Simulating PET images from real CT images.
    3. (Future) Performing lesion insertion and analysis.
    :)

For any questions please contact:
Peter Yazdi  <pyazdi@bccrc.ca>
"""

import os
import json

from json_minify import json_minify

from modules.spect_pre_process.data_preprocessing import preprocess_data_for_simind
from modules.spect_post_process.image_analysis import resample_spect_to_atn_grid
from modules.spect_pre_process.run_totseg import run_totseg
from modules.pbpk.run_pbpk import run_pbpk
from modules.spect_simulation.run_simind import run_simind
from modules.spect_simulation.run_recon import run_recon

class TdtPipeline:
    def __init__(self, config_path):
        self.config_path = config_path
        self.current_dir_path = os.path.dirname(__file__)
        self.setup_config(config_path)
    
    def setup_config(self, config_path):
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        with open(config_path, encoding="utf-8") as f:
            self.config = json.loads(json_minify(f.read()))

        output_folder_title = self.config["output_folder"]["title"]
        output_folder_path = os.path.join(self.current_dir_path, output_folder_title)
        os.makedirs(output_folder_path, exist_ok=True)

        subdir_dict = self.config["subdir_names"]
        for _, name in subdir_dict.items():
            subdir_path = os.path.join(output_folder_path, name)
            os.makedirs(subdir_path, exist_ok=True)
    
    def run(self):
        context = Context()
        # stage 1 : tdt segementation 
        #context = run_totseg(self.config,context)
        context.ct_nii_path = '/home/jhubadmin/Theranostic-Virtual-Patient-Pipeline/TDT_Output_classes_test3/spect_preprocessing_outputs/spect_preprocessing_ct.nii.gz'
        context.roi_seg_path = '/home/jhubadmin/Theranostic-Virtual-Patient-Pipeline/TDT_Output_classes_test3/spect_preprocessing_outputs/spect_preprocessing_roi_seg.nii.gz'
        context.body_seg_path = '/home/jhubadmin/Theranostic-Virtual-Patient-Pipeline/TDT_Output_classes_test3/spect_preprocessing_outputs/spect_preprocessing_body_seg.nii.gz'
        
        # stage 2 : preprocess ct and seg for simind
        context = preprocess_data_for_simind(self.config, context)
        
        # stage 3 : run pbpk
        context = run_pbpk(self.config,context)

        # stage 4 : run simind
        #context = run_simind(self.current_dir_path, self.config, context)

        # stage 5 : run recon
        context = run_recon(self.config, context)
        
        for i in context.__dict__:
            print(f"\n{i}:\n{getattr(context, i)}")
        
        
        return context
    
class Context:
    def __init__(self):
        pass
        

if __name__ == "__main__":
    config_path = "setup/config.json"
    pipeline1 = TdtPipeline(config_path)
    pipeline1.run()
    