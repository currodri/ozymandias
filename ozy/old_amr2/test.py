import numpy as np
import _amr2prof
profiles = np.zeros((1,1,3,100),order='F')
xprofiles = np.zeros((1,100),order='F')
args = '/mnt/extraspace/currodri/NUT/cosmoNUThd/output_00044',1,1,3,['r_cyl'],['v_rot'],['mass','volume','metallicity'],100,'cylinder',np.array([0.70336995, 0.33568012, 0.341516,0.0,0.0028299868106842043,-0.0028299868106842043,0.0028299868106842043,-0.35676084, 0.40252779, 0.84302614,28.26907349, -20.34032822, -25.46469116],order='F'),10,False,profiles,xprofiles

_amr2prof.f90wrap_onedprofile(*args)
print(profiles,xprofiles)