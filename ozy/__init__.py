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
import variables_settings