import numpy as np
import ozy
from ozy.profiles import compute_profile

obj = ozy.load('test_00010.hdf5')
masses = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(masses)
gal = obj.galaxies[progind]

prof = compute_profile(gal,'test_00010.hdf5', 'r_sphere', ['gas/density','star/v_sphere_r','dm/v_sphere_r'],
                        ['gas/volume','star/mass','dm/mass'],save=True,recompute=True,nbins=10,rmax=(15,'kpc'))
print(prof.ydata)