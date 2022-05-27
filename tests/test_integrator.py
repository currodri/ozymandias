import numpy as np
import ozy
import sys
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
from amr2 import amr_integrator
from part2 import part_integrator
from ozy.utils import init_region,init_filter

obj = ozy.load('test_00035.hdf5')
virial_mass = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(virial_mass)

gal = obj.galaxies[progind]
output_path = obj.simulation.fullpath

# Initialise region
selected_reg = init_region(gal,'sphere',rmax=(0.12,'rvir'))

nvar = 0
quantity_names = ['mass','density','temperature']
weight_names = ['cumulative','mass','volume']
        
# Initialise Fortran derived type with attributes
# This object hold the following attributes:
# - nvars: number of variables
# - nwvars:  number of variables for weighting
# - varnames: names of variables
# - wvarnames:  names of variables for weighting
# - data: organised in numpy array of shape (nvars,nwvars,4)
#           each of those 4 values are (final, min, max, sum of weights)
glob_attrs = amr_integrator.amr_region_attrs()
glob_attrs.nvars = len(quantity_names)
glob_attrs.nwvars = len(weight_names)
glob_attrs.nfilter = 3
amr_integrator.allocate_amr_regions_attrs(glob_attrs)
for i in range(0, len(quantity_names)):
    glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
for i in range(0, len(weight_names)):
    glob_attrs.wvarnames.T.view('S128')[i] = weight_names[i].ljust(128)

glob_attrs.filters[0] = init_filter(cond_strs=['entropy_specific/</4.4e+8/erg*K**-1*g**-1'],name='cold',group=gal)
glob_attrs.filters[1] = init_filter(cond_strs=['entropy_specific/>=/4.4e+8/erg*K**-1*g**-1','entropy_specific/<=/23.2e+8/erg*K**-1*g**-1'],name='warm',group=gal)
glob_attrs.filters[2] = init_filter(cond_strs=['entropy_specific/>/23.2e+8/erg*K**-1*g**-1'],name='hot',group=gal)

# Begin integration
amr_integrator.integrate_region(output_path,selected_reg,glob_attrs)
mass_names = ['COLD','WARM','HOT']
tot_mass = 0
for i in range(0,3):
    mass = obj.quantity(glob_attrs.data[i,0,0,0],'code_mass')
    tot_mass += mass
    print(mass_names[i]+': ',mass.to('Msun'))
print('Total mass: ', tot_mass.to('Msun'))

glob_attrs = part_integrator.part_region_attrs()
glob_attrs.nvars = 11
glob_attrs.nwvars = 2
part_integrator.allocate_part_regions_attrs(glob_attrs)
glob_attrs.varnames.T.view('S128')[0] = b'star/mass'.ljust(128)
# TODO: SFR indicators should be taken as arguments
glob_attrs.varnames.T.view('S128')[1] = b'star/sfr_10'.ljust(128)
glob_attrs.varnames.T.view('S128')[2] = b'star/sfr_100'.ljust(128)
glob_attrs.varnames.T.view('S128')[3] = b'star/metallicity'.ljust(128)

glob_attrs.varnames.T.view('S128')[4] = b'star/ang_momentum_x'.ljust(128)
glob_attrs.varnames.T.view('S128')[5] = b'star/ang_momentum_y'.ljust(128)
glob_attrs.varnames.T.view('S128')[6] = b'star/ang_momentum_z'.ljust(128)

glob_attrs.varnames.T.view('S128')[7] = b'dm/mass'.ljust(128)
glob_attrs.varnames.T.view('S128')[8] = b'dm/ang_momentum_x'.ljust(128)
glob_attrs.varnames.T.view('S128')[9] = b'dm/ang_momentum_y'.ljust(128)
glob_attrs.varnames.T.view('S128')[10] = b'dm/ang_momentum_z'.ljust(128)


glob_attrs.wvarnames.T.view('S128')[0] = b'cumulative'.ljust(128)
glob_attrs.wvarnames.T.view('S128')[1] = b'mass'.ljust(128)

# Begin integration
filt = init_filter(cond_strs=['none'],name='none',group=gal)
part_integrator.integrate_region(output_path,selected_reg,filt,glob_attrs,False,True)

