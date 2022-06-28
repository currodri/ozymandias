import sys
import os
from ozy.utils import get_SNevents_log,sn_data_hdf5


# 1. Check main reading function
# print('Reading logfile '+str(sys.argv[1]))
# sn_data = get_SNevents_log(sys.argv[1])
# print(sn_data)

# 2. Check full catalogue function
os.chdir('/mnt/extraspace/currodri/NUT/cosmoNUTcrmhd')
# os.chdir('/mnt/extraspace/currodri/NUT/cosmoNUThd')
sn_data = sn_data_hdf5(['cosmoNUTcrmhd.o797249','cosmoNUTcrmhd.o798029'],have_crs=True,filename='sn_test.hdf5')
# sn_data = sn_data_hdf5(['cosmoNUThd.o872195'],have_crs=True)
# sn_data = sn_data_hdf5(['cosmoNUTcrmhd.o797249','cosmoNUTcrmhd.o798029','cosmoNUTcrmhd.o803999'],have_crs=True)
# sn_data = sn_data_hdf5(['cosmoNUTcrmhd.o798029'],have_crs=True)