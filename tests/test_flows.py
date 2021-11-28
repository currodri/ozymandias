import numpy as np
import ozy
from ozy.outflow_inflow import compute_flows

obj = ozy.load('test_00010.hdf5')
masses = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(masses)
gal = obj.galaxies[progind]

compute_flows(gal,'test_00010.hdf5','outflow',rmin=(0.19,'rvir'),rmax=(0.21,'rvir'),save=True,recompute=True)
compute_flows(gal,'test_00010.hdf5','inflow',rmin=(0.49,'rvir'),rmax=(0.51,'rvir'),save=True,recompute=True)