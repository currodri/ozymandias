import numpy as np
import ozy
from ozy.profiles import compute_profile
from ozy.utils import init_region
import matplotlib.pyplot as plt
import sys
import yt
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
from amr2 import io_ramses
from amr2 import filtering,amr_integrator

obj = ozy.load('test_00035.hdf5')
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

prof = compute_profile(gal,'test_00035.hdf5', 'r_sphere', ['gas/density'],
                        ['gas/cumulative','gas/volume','gas/density'],
                        region_type='sphere',save=True,recompute=True,nbins=80,rmax=(4.0,'kpc'),
                        logscale=True)

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(6,4), dpi=100, facecolor='w', edgecolor='k')
        
ax.set_ylabel(r'$\rho$ [g cm$^{-3}$]', fontsize=16)
ax.set_xlabel(r'$r$ [kpc]', fontsize=16)
ax.set_yscale('log')
ax.set_xscale('log')
ax.tick_params(labelsize=12,direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.minorticks_on()
ax.tick_params(which='both',axis="both",direction="in")

print(prof.ydata['hydro'][0].shape)

new_x = 0.5*(prof.xdata[0][:-1]+prof.xdata[0][1:])
new_x = np.concatenate(([0.5*prof.xdata[0][0]], new_x),axis=0)
new_x = obj.array(new_x,'code_length')
print(new_x.in_units('kpc'),prof.ydata['hydro'][0][:,1,0].in_units('g*cm**-3'))
ax.plot(new_x.in_units('kpc'),prof.ydata['hydro'][0][:,1,0].in_units('g*cm**-3'), label='ozymandias')

# prof = compute_profile(gal,'test_00035.hdf5', 'r_sphere', ['gas/density'],
#                         ['gas/cumulative','gas/volume','gas/density'],
#                         region_type='sphere',save=True,recompute=True,nbins=20,rmax=(1.0,'kpc'))

# new_x = 0.5*(prof.xdata[0][:-1]+prof.xdata[0][1:])
# new_x = np.concatenate(([0.5*prof.xdata[0][0]], new_x),axis=0)
# new_x = obj.array(new_x,'code_length')
# print(new_x.in_units('kpc'),prof.ydata['hydro'][0][:,1,0].in_units('g*cm**-3'))
# ax.plot(new_x.in_units('kpc'),prof.ydata['hydro'][0][:,1,0].in_units('g*cm**-3'), label='ozymandias 1.0',alpha=0.4)

# prof = compute_profile(gal,'test_00035.hdf5', 'r_sphere', ['gas/density'],
#                         ['gas/cumulative','gas/volume','gas/density'],
#                         region_type='sphere',save=True,recompute=True,nbins=40,rmax=(2.0,'kpc'))

# new_x = 0.5*(prof.xdata[0][:-1]+prof.xdata[0][1:])
# new_x = np.concatenate(([0.5*prof.xdata[0][0]], new_x),axis=0)
# new_x = obj.array(new_x,'code_length')
# print(new_x.in_units('kpc'),prof.ydata['hydro'][0][:,1,0].in_units('g*cm**-3'))
# ax.plot(new_x.in_units('kpc'),prof.ydata['hydro'][0][:,1,0].in_units('g*cm**-3'), label='ozymandias 2.0',alpha=0.4)

# prof = compute_profile(gal,'test_00035.hdf5', 'r_sphere', ['gas/density'],
#                         ['gas/cumulative','gas/volume','gas/density'],
#                         region_type='sphere',save=True,recompute=True,nbins=80,rmax=(4.0,'kpc'))

# new_x = 0.5*(prof.xdata[0][:-1]+prof.xdata[0][1:])
# new_x = np.concatenate(([0.5*prof.xdata[0][0]], new_x),axis=0)
# new_x = obj.array(new_x,'code_length')
# print(new_x.in_units('kpc'),prof.ydata['hydro'][0][:,1,0].in_units('g*cm**-3'))
# ax.plot(new_x.in_units('kpc'),prof.ydata['hydro'][0][:,1,0].in_units('g*cm**-3'), label='ozymandias 4.0',alpha=0.4)

# # Doing it with yt
ds = yt.load('/mnt/extraspace/currodri/NUT/cosmoNUTcrmhd/output_00035/info_00035.txt')
sphere = ds.sphere(center=(gal.position[0].d, gal.position[1].d, gal.position[2].d), radius=(4.0, "kpc"))
rp0 = yt.create_profile(
    sphere,
    ("index", "radius"),
    ("gas", "density"),
    weight_field=("index", "cell_volume"),
    n_bins=80
)
print(rp0.x.in_units('kpc'))
ax.plot(rp0.x.in_units('kpc'),rp0['gas','density'].in_units('g*cm**-3'), label='yT')

ax.legend(loc='best')
fig.savefig('NUT_00035_radial_density.png',dpi=200,format='png')

# Initialise region
# selected_reg = init_region(gal,'sphere',rmin=(0.0,'kpc'),rmax=(10.0,'kpc'))
# all_filt = filtering.filter()

# gas_attrs = amr_integrator.amr_region_attrs()
# gas_attrs.nvars = 1
# gas_attrs.nwvars = 1
# gas_attrs.nfilter = 1
# amr_integrator.allocate_amr_regions_attrs(gas_attrs)
# gas_attrs.varnames.T.view('S128')[0] = 'mass'.ljust(128)
# gas_attrs.wvarnames.T.view('S128')[0] = 'cumulative'.ljust(128)

# gas_attrs.filters[0] = all_filt

# # Begin integration
# amr_integrator.integrate_region(obj.simulation.fullpath,selected_reg,gas_attrs)

# tot_gas = obj.quantity(gas_attrs.data[0,0,0,0],'code_mass')

# print('Comparing masses: ',tot_gas,np.sum(prof.ydata['hydro'][0][:,0,0]))