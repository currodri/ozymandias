import os
import numpy as np
import ozy
from ozy.phase_diagrams import compute_phase_diagram,plot_single_phase_diagram

obj = ozy.load('test_00035.hdf5')
virial_mass = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(virial_mass)

gal = obj.galaxies[progind]

pd = compute_phase_diagram(gal,'test_00035.hdf5', 'density','temperature', 
                            ['gas/mass','gas/alfvendiff_ratio','gas/stheatcooling_ratio','gas/streaming_heating','gas/net_cooling','gas/total_coolingtime'],
                            ['gas/cumulative','gas/mass','gas/density','gas/volume'],save=True,recompute=False,rmax=(0.2,'rvir'),
                            nbins=[100,100],filter_name='galaxy',cr_st=True,cr_heat=True)#,filter_conds='eff_FKmag/>=/0.01/dimensionless',
                            #filter_name='starforming')

print('Tot_mass: ',np.sum(pd.zdata['hydro'][0][:,:,0][:,:,0]))
plot_single_phase_diagram(pd,'mass','NUT_00035_pd_mass',stats='mean',weightvar='cumulative',gent=True)
plot_single_phase_diagram(pd,'alfvendiff_ratio','NUT_00035_pd_alfvendiff_ratio',stats='mean',weightvar='volume',gent=True)
plot_single_phase_diagram(pd,'stheatcooling_ratio','NUT_00035_pd_stheatcooling_ratio',stats='mean',weightvar='volume',gent=True)
plot_single_phase_diagram(pd,'streaming_heating','NUT_00035_pd_streaming_heating',stats='mean',weightvar='volume',gent=True)
plot_single_phase_diagram(pd,'net_cooling','NUT_00035_pd_net_cooling',stats='mean',weightvar='volume',gent=True)
plot_single_phase_diagram(pd,'total_coolingtime','NUT_00035_pd_total_coolingtime',stats='mean',weightvar='volume',gent=True)

pd = compute_phase_diagram(gal,'test_00035.hdf5', 'density','temperature', 
                            ['gas/mass','gas/alfvendiff_ratio','gas/stheatcooling_ratio','gas/streaming_heating','gas/net_cooling','gas/total_coolingtime'],
                            ['gas/cumulative','gas/mass','gas/density','gas/volume'],save=True,recompute=False,rmax=(1.0,'rvir'),
                            rmin=(0.2,'rvir'),nbins=[100,100],filter_name='halo',cr_st=True,cr_heat=True)

print('Tot_mass: ',np.sum(pd.zdata['hydro'][0][:,:,0][:,:,0]))
plot_single_phase_diagram(pd,'mass','halo_00035_pd_mass',stats='mean',weightvar='cumulative',gent=True)
plot_single_phase_diagram(pd,'alfvendiff_ratio','halo_00035_pd_alfvendiff_ratio',stats='mean',weightvar='volume',gent=True)
plot_single_phase_diagram(pd,'stheatcooling_ratio','halo_00035_pd_stheatcooling_ratio',stats='mean',weightvar='volume',gent=True)
plot_single_phase_diagram(pd,'streaming_heating','halo_00035_pd_streaming_heating',stats='mean',weightvar='volume',gent=True)
plot_single_phase_diagram(pd,'net_cooling','halo_00035_pd_net_cooling',stats='mean',weightvar='volume',gent=True)
plot_single_phase_diagram(pd,'total_coolingtime','halo_00035_pd_total_coolingtime',stats='mean',weightvar='volume',gent=True)