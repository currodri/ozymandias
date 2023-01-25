import os
import numpy as np
import ozy
from ozy.phase_diagrams import compute_phase_diagram,plot_single_phase_diagram

obj = ozy.load('test_00035.hdf5')
virial_mass = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(virial_mass)

gal = obj.galaxies[progind]

pd = compute_phase_diagram(gal,'test_00035.hdf5', 'density','temperature', ['gas/mass'],
                            ['gas/cumulative','gas/mass'],save=True,recompute=True,rmax=(0.2,'rvir'),
                            nbins=[500,500])#,filter_conds='eff_FKmag/>=/0.01/dimensionless',
                            #filter_name='starforming')

print('Tot_mass: ',np.sum(pd.zdata['hydro'][0][:,:,0][:,:,0]))
plot_single_phase_diagram(pd,'mass','NUT_00035_pd_mass',stats='mean',gent=True)