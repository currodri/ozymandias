import numpy as np
import ozy
import sys
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
from amr2 import amr_integrator
from ozy.utils import init_region,init_filter

obj = ozy.load('test_00010.hdf5')
virial_mass = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(virial_mass)

gal = obj.galaxies[progind]
output_path = obj.simulation.fullpath

# Initialise region
selected_reg = init_region(gal,'sphere')

filt = init_filter(cond_strs=['entropy_specific/</3.7e+8/erg*K**-1*g**-1'],name='cold',group=gal)

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
amr_integrator.allocate_amr_regions_attrs(glob_attrs)
for i in range(0, len(quantity_names)):
    glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
for i in range(0, len(weight_names)):
    glob_attrs.wvarnames.T.view('S128')[i] = weight_names[i].ljust(128)

# Begin integration
amr_integrator.integrate_region(output_path,selected_reg,filt,glob_attrs)
print('COLD')
print(glob_attrs.data)

filt = init_filter(cond_strs=['entropy_specific/>=/3.7e+8/erg*K**-1*g**-1','entropy_specific/<=/23.2e+8/erg*K**-1*g**-1'],name='warm',group=gal)

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
amr_integrator.allocate_amr_regions_attrs(glob_attrs)
for i in range(0, len(quantity_names)):
    glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
for i in range(0, len(weight_names)):
    glob_attrs.wvarnames.T.view('S128')[i] = weight_names[i].ljust(128)

# Begin integration
amr_integrator.integrate_region(output_path,selected_reg,filt,glob_attrs)
print('WARM')
print(glob_attrs.data)

filt = init_filter(cond_strs=['entropy_specific/>/23.2e+8/erg*K**-1*g**-1'],name='hot',group=gal)

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
amr_integrator.allocate_amr_regions_attrs(glob_attrs)
for i in range(0, len(quantity_names)):
    glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
for i in range(0, len(weight_names)):
    glob_attrs.wvarnames.T.view('S128')[i] = weight_names[i].ljust(128)

# Begin integration
amr_integrator.integrate_region(output_path,selected_reg,filt,glob_attrs)
print('HOT')
print(glob_attrs.data)