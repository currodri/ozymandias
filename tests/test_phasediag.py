import os
import numpy as np
import ozy
from ozy.phase_diagrams import compute_phase_diagram,plot_single_phase_diagram

obj = ozy.load('test_00010.hdf5')
virial_mass = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(virial_mass)

gal = obj.galaxies[progind]

minrho = obj.quantity(1e-30,'g/cm**3')
maxrho = obj.quantity(1e-20,'g/cm**3')
minT = obj.quantity(1e0,'K')
maxT = obj.quantity(1e8,'K')
pds = compute_phase_diagram(gal,'test_00010.hdf5', 'density','temperature', 
                            ['gas/mass','gas/alfvendiff_ratio','gas/stheatcooling_ratio','gas/streaming_heating','gas/net_cooling','gas/entropy_specific','gas/absv_sphere_r'],
                            ['gas/cumulative','gas/mass','gas/density','gas/volume'],
                            [minrho,minT],[maxrho,maxT],
                            save=True,recompute=True,rmax=(0.2,'rvir'),
                            nbins=[100,100],filter_name=['all','cold','warm','hot'],
                            filter_conds=['none','entropy_specific/</4.4e+8/erg*K**-1*g**-1',
                                      ['entropy_specific/>=/4.4e+8/erg*K**-1*g**-1','entropy_specific/<=/23.2e+8/erg*K**-1*g**-1'],
                                      'entropy_specific/>/23.2e+8/erg*K**-1*g**-1'],
                            cr_st=True,cr_heat=True,remove_subs=False)

# print('Tot_mass: ',np.sum(pd.zdata['hydro'][0][:,:,0][:,:,0]))
plot_single_phase_diagram(pds[0],'mass','NUT_00035_pd_mass',stats='mean',weightvar='cumulative',gent=True)
plot_single_phase_diagram(pds[1],'mass','NUT_00035_pd_mass_cold',stats='mean',weightvar='cumulative',gent=True)
plot_single_phase_diagram(pds[2],'mass','NUT_00035_pd_mass_warm',stats='mean',weightvar='cumulative',gent=True)
plot_single_phase_diagram(pds[3],'mass','NUT_00035_pd_mass_hot',stats='mean',weightvar='cumulative',gent=True)
plot_single_phase_diagram(pds[0],'absv_sphere_r','NUT_00035_pd_absv_sphere_r',stats='mean',weightvar='mass',gent=True)
# plot_single_phase_diagram(pd,'alfvendiff_ratio','NUT_00035_pd_alfvendiff_ratio',stats='mean',weightvar='volume',gent=True)
# plot_single_phase_diagram(pd,'stheatcooling_ratio','NUT_00035_pd_stheatcooling_ratio',stats='mean',weightvar='volume',gent=True)
# plot_single_phase_diagram(pd,'streaming_heating','NUT_00035_pd_streaming_heating',stats='mean',weightvar='volume',gent=True)
# plot_single_phase_diagram(pd,'net_cooling','NUT_00035_pd_net_cooling',stats='mean',weightvar='volume',gent=True)
# plot_single_phase_diagram(pd,'entropy_specific','NUT_00035_pd_entropy_specific',stats='mean',weightvar='density',gent=True)

# pd = compute_phase_diagram(gal,'test_00035.hdf5', 'density','temperature', 
#                             ['gas/mass','gas/alfvendiff_ratio','gas/stheatcooling_ratio','gas/streaming_heating','gas/net_cooling','gas/total_coolingtime'],
#                             ['gas/cumulative','gas/mass','gas/density','gas/volume'],save=True,recompute=True,rmax=(1.0,'rvir'),
#                             rmin=(0.2,'rvir'),nbins=[100,100],filter_name='halo',cr_st=True,cr_heat=True)

# print('Tot_mass: ',np.sum(pd.zdata['hydro'][0][:,:,0][:,:,0]))
# plot_single_phase_diagram(pd,'mass','halo_00035_pd_mass',stats='mean',weightvar='cumulative',gent=True)
# plot_single_phase_diagram(pd,'alfvendiff_ratio','halo_00035_pd_alfvendiff_ratio',stats='mean',weightvar='volume',gent=True)
# plot_single_phase_diagram(pd,'stheatcooling_ratio','halo_00035_pd_stheatcooling_ratio',stats='mean',weightvar='volume',gent=True)
# plot_single_phase_diagram(pd,'streaming_heating','halo_00035_pd_streaming_heating',stats='mean',weightvar='volume',gent=True)
# plot_single_phase_diagram(pd,'net_cooling','halo_00035_pd_net_cooling',stats='mean',weightvar='volume',gent=True)
# plot_single_phase_diagram(pd,'total_coolingtime','halo_00035_pd_total_coolingtime',stats='mean',weightvar='volume',gent=True)