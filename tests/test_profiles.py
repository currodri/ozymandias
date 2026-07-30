import numpy as np
import ozy
from ozy.profiles import compute_profile
from ozy.utils import init_region
import matplotlib.pyplot as plt
import sys
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
from amr2 import io_ramses
from amr2 import filtering,amr_integrator
io_ramses.verbose = True
obj = ozy.load('fake_00035.hdf5')
masses = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(masses)
gal = obj.galaxies[progind]

print('Center of the galaxy: ', gal.position[0].d, gal.position[1].d, gal.position[2].d)
radius = obj.quantity(1.0,'kpc')
print('Radius of the sphere: ', radius.to('code_length').d)
norm_L = gal.angular_mom['total']/np.linalg.norm(gal.angular_mom['total'])
print('Angular momentum direction: ',norm_L[0].d, norm_L[1].d, norm_L[2].d, np.linalg.norm(norm_L))
velocity = gal.velocity.in_units('code_velocity').d
print('NUT velocity: ', velocity)
zmin = obj.quantity(-4.0,'kpc')
zmax = obj.quantity(4.0, 'kpc')
linthresh = obj.quantity(0.01,'kpc')
profs = compute_profile(gal,'fake_00035.hdf5', 'z', ['gas/density','gas/pseudo_entropy','gas/grad_entropy_rsphere',
                                                     'gas/mass','star/mass','dm/mass'],
                        ['gas/cumulative','gas/volume','gas/density','star/cumulative','dm/cumulative'],
                        zmin, zmax,linthresh=linthresh,
                        region_type='custom_cylinder',save=True,recompute=True,nbins=80,
                        rmin=(0.0,'kpc'),rmax=(4.0,'kpc'),zmin=(-4.0,'kpc'),zmax=(4.0,'kpc'),
                        logscale=True,myaxis=norm_L,mycentre=(gal.position.d,gal.position.units),
                        filter_conds=['none','entropy_specific/</4.4e+8/erg*K**-1*g**-1',
                                      ['entropy_specific/>=/4.4e+8/erg*K**-1*g**-1','entropy_specific/<=/23.2e+8/erg*K**-1*g**-1'],
                                      'entropy_specific/>/23.2e+8/erg*K**-1*g**-1',
                                      'v_sphere_r/<=/-100/km*s**-1'],
                        filter_name=['all','cold','warm','hot','inflow'],
                        regime_type='galaxy',force_neigh=True)

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(6,4), dpi=100, facecolor='w', edgecolor='k')
        
ax.set_ylabel(r'$\rho$ [g cm$^{-3}$]', fontsize=16)
ax.set_xlabel(r'$z$ [kpc]', fontsize=16)
ax.set_yscale('log')
ax.set_xscale('symlog')
ax.tick_params(labelsize=12,direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.minorticks_on()
ax.tick_params(which='both',axis="both",direction="in")

new_x = 0.5*(profs[0].xdata[0][:-1]+profs[0].xdata[0][1:])
ind = profs[0].yvars['hydro'].index('density')
ax.plot(new_x.in_units('kpc'),profs[0].ydata['hydro'][ind][:,1,1].in_units('g*cm**-3'), label='all')
ax.plot(new_x.in_units('kpc'),profs[1].ydata['hydro'][ind][:,1,1].in_units('g*cm**-3'), label='cold')
ax.plot(new_x.in_units('kpc'),profs[2].ydata['hydro'][ind][:,1,1].in_units('g*cm**-3'), label='warm')
ax.plot(new_x.in_units('kpc'),profs[3].ydata['hydro'][ind][:,1,1].in_units('g*cm**-3'), label='hot')
ax.plot(new_x.in_units('kpc'),profs[4].ydata['hydro'][ind][:,1,1].in_units('g*cm**-3'), label='warm inflow')
ax.legend(loc='best',frameon=False)
fig.savefig('NUT_00035_z_density.png',dpi=200,format='png')

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(6,4), dpi=100, facecolor='w', edgecolor='k')
        
ax.set_ylabel(r'$K$ [KeV cm$^{2}$]', fontsize=16)
ax.set_xlabel(r'$z$ [kpc]', fontsize=16)
ax.set_yscale('log')
ax.set_xscale('symlog')
ax.tick_params(labelsize=12,direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.minorticks_on()
ax.tick_params(which='both',axis="both",direction="in")

new_x = 0.5*(profs[0].xdata[0][:-1]+profs[0].xdata[0][1:])
ind = profs[0].yvars['hydro'].index('pseudo_entropy')
ax.plot(new_x.in_units('kpc'),profs[0].ydata['hydro'][ind][:,1,1].in_units('keV*cm**2'), label='all')
ax.plot(new_x.in_units('kpc'),profs[1].ydata['hydro'][ind][:,1,1].in_units('keV*cm**2'), label='cold')
ax.plot(new_x.in_units('kpc'),profs[2].ydata['hydro'][ind][:,1,1].in_units('keV*cm**2'), label='warm')
ax.plot(new_x.in_units('kpc'),profs[3].ydata['hydro'][ind][:,1,1].in_units('keV*cm**2'), label='hot')
ax.plot(new_x.in_units('kpc'),profs[4].ydata['hydro'][ind][:,1,1].in_units('keV*cm**2'), label='warm inflow')
ax.legend(loc='best',frameon=False)
fig.savefig('NUT_00035_z_pseudo_entropy.png',dpi=200,format='png')

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(6,4), dpi=100, facecolor='w', edgecolor='k')
        
ax.set_ylabel(r'$\nabla K$ [KeV cm]', fontsize=16)
ax.set_xlabel(r'$z$ [kpc]', fontsize=16)
ax.set_yscale('symlog')
ax.set_xscale('symlog')
ax.tick_params(labelsize=12,direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.minorticks_on()
ax.tick_params(which='both',axis="both",direction="in")

new_x = 0.5*(profs[0].xdata[0][:-1]+profs[0].xdata[0][1:])
ind = profs[0].yvars['hydro'].index('grad_entropy_rsphere')
ax.plot(new_x.in_units('kpc'),profs[0].ydata['hydro'][ind][:,1,1].in_units('keV*cm'), label='all')
ax.plot(new_x.in_units('kpc'),profs[1].ydata['hydro'][ind][:,1,1].in_units('keV*cm'), label='cold')
ax.plot(new_x.in_units('kpc'),profs[2].ydata['hydro'][ind][:,1,1].in_units('keV*cm'), label='warm')
ax.plot(new_x.in_units('kpc'),profs[3].ydata['hydro'][ind][:,1,1].in_units('keV*cm'), label='hot')
ax.plot(new_x.in_units('kpc'),profs[4].ydata['hydro'][ind][:,1,1].in_units('keV*cm'), label='warm inflow')
ax.legend(loc='best',frameon=False)
fig.savefig('NUT_00035_z_grad_entropy.png',dpi=200,format='png')

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(6,4), dpi=100, facecolor='w', edgecolor='k')
        
ax.set_ylabel(r'Mass [$M_{\odot}$]', fontsize=16)
ax.set_xlabel(r'$z$ [kpc]', fontsize=16)
ax.set_yscale('symlog')
ax.set_xscale('symlog')
ax.tick_params(labelsize=12,direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.minorticks_on()
ax.tick_params(which='both',axis="both",direction="in")

new_x = 0.5*(profs[0].xdata[0][:-1]+profs[0].xdata[0][1:])
ind = profs[0].yvars['hydro'].index('mass')
ax.plot(new_x.in_units('kpc'),profs[0].ydata['hydro'][ind][:,0,0].in_units('Msun'), label='gas')
print(profs[0].ydata['star'][0][:,0,0].in_units('Msun'))
ax.plot(new_x.in_units('kpc'),profs[0].ydata['star'][0][:,0,0].in_units('Msun'), label='stars')
print(profs[0].ydata['dm'][0][:,0,0].in_units('Msun'))
ax.plot(new_x.in_units('kpc'),profs[0].ydata['dm'][0][:,0,0].in_units('Msun'), label='DM')
ax.legend(loc='best',frameon=False)
fig.savefig('NUT_00035_z_mass.png',dpi=200,format='png')
