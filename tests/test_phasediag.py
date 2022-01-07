import os
import numpy as np
import ozy
from ozy.phase_diagrams import compute_phase_diagram,plot_single_phase_diagram

obj = ozy.load('test_00010.hdf5')
virial_mass = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(virial_mass)

gal = obj.galaxies[progind]

pd = compute_phase_diagram(gal,'test_00010.hdf5', 'density','temperature', ['gas/mass'],
                            ['gas/cumulative','gas/mass'],save=True,recompute=False,rmax=(20,'kpc'))

plot_single_phase_diagram(pd,'mass','NUT_00010_pd_mass',stats='mean',gent=True)