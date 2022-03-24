import numpy as np
import ozy
from ozy.profiles import compute_profile

obj = ozy.load('test_00010.hdf5')
masses = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(masses)
gal = obj.galaxies[progind]

prof = compute_profile(gal,'test_00010.hdf5', 'r_cyl', ['gas/mass','star/mass','dm/mass'],
                        ['gas/cumulative','star/cumulative','dm/cumulative'],region_type='sphere',save=False,recompute=True,nbins=50,rmax=(15,'kpc'))
print(prof.ydata)