import os
OZYPATH = os.path.realpath(__file__).split('/__init__.py')[0]
AMRPATH = os.path.join(OZYPATH,'amr')
PARTPATH = os.path.join(OZYPATH,'part')
VISPATH = os.path.join(OZYPATH,'visualisation')
import sys
sys.path.append(OZYPATH)
sys.path.append(AMRPATH)
sys.path.append(PARTPATH)
sys.path.append(VISPATH)

from ozy.loader import load
from ozy.main import Snapshot, CosmoSnapshot
from ozy.driver import drive

import importlib.util
from ozy.variables_settings import plotting_dictionary,circle_dictionary,basic_conv
def load_custom_settings():
    """Load custom settings if 'ozy_settings.py' exists in the current directory."""
    global plotting_dictionary,circle_dictionary,basic_conv  # Declare the dictionaries as global to override them
    
    current_dir = os.getcwd()
    ozy_settings_path = os.path.join(current_dir, 'ozy_settings.py')
    
    if os.path.isfile(ozy_settings_path):
        # Load the ozy_settings module dynamically
        spec = importlib.util.spec_from_file_location("ozy_settings", ozy_settings_path)
        ozy_settings = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(ozy_settings)
        
        # Override the dictionaries
        plotting_dictionary = getattr(ozy_settings, 'plotting_dictionary', plotting_dictionary)
        circle_dictionary = getattr(ozy_settings, 'circle_dictionary', circle_dictionary)
        basic_conv = getattr(ozy_settings, 'basic_conv', basic_conv)

# Call the function to load custom settings at package initialization
load_custom_settings()