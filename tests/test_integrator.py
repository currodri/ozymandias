import numpy as np
import ozy
import sys
from amr2 import amr_integrator,stats_utils
from part2 import part_integrator
from ozy.utils import init_region,init_filter,\
                        structure_regions,get_code_bins,\
                        pdf_handler_to_stats
from ozy.dict_variables import get_code_units

obj = ozy.load('test_00035.hdf5')
virial_mass = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(virial_mass)

gal = obj.galaxies[progind]
output_path = obj.simulation.fullpath

# Initialise region
selected_reg = init_region(gal,'sphere',rmax=(0.2,'rvir'))

# subs = structure_regions(gal, add_substructure=True, add_neighbours=False,
#                                     tidal_method='BT87_simple')
# nsubs = len(subs)
subs=[];nsubs=0

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
glob_attrs.nfilter = 4
glob_attrs.nsubs = nsubs
amr_integrator.allocate_amr_regions_attrs(glob_attrs)
for i in range(0, len(quantity_names)):
    glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
    glob_attrs.result[i].nbins = 200
    glob_attrs.result[i].nfilter = 4
    glob_attrs.result[i].nwvars = len(weight_names)
    glob_attrs.result[i].varname = quantity_names[i]
    mybins = get_code_bins(obj,'gas/'+quantity_names[i],200)
    glob_attrs.result[i].scaletype = mybins[1]
    stats_utils.allocate_pdf(glob_attrs.result[i])
    glob_attrs.result[i].bins = mybins[0]
    for j in range(0, len(weight_names)):
       glob_attrs.result[i].wvarnames.T.view('S128')[j] = weight_names[j].ljust(128)
       
for i in range(0,nsubs):
    glob_attrs.subs[i] = subs[i]

glob_attrs.filters[0] = init_filter(cond_strs=['entropy_specific/</4.4e+8/erg*K**-1*g**-1'],name='cold',group=gal)
glob_attrs.filters[1] = init_filter(cond_strs=['entropy_specific/>=/4.4e+8/erg*K**-1*g**-1','entropy_specific/<=/23.2e+8/erg*K**-1*g**-1'],name='warm',group=gal)
glob_attrs.filters[2] = init_filter(cond_strs=['entropy_specific/>/23.2e+8/erg*K**-1*g**-1'],name='hot',group=gal)
glob_attrs.filters[3] = init_filter(cond_strs=['none'],name='all',group=gal)
# Begin integration

amr_integrator.integrate_region(output_path,selected_reg,False,glob_attrs)
mass_names = ['COLD','WARM','HOT']
tot_mass = 0
for i in range(0,3):
    varname = str(glob_attrs.result[0].wvarnames.T.view('S128')[1][0].decode("utf-8")).rstrip()
    unit = get_code_units(varname)
    mass = glob_attrs.result[0].totweights[i,0]
    mass = obj.quantity(mass,unit)
    tot_mass += mass
    print(mass_names[i]+': ',mass.to('Msun'))
print('Total mass: ', tot_mass.to('Msun'))
mass = glob_attrs.result[0].totweights[3,1]
mass = obj.quantity(mass,'code_mass')
print('Total mass no filter: ',mass.to('Msun'))
print('From Ozymandias: ',gal.mass['gas'].to('Msun'))
print('Ratio filter: ',tot_mass.to('Msun')/gal.mass['gas'].to('Msun'))
print('Ratio no filter: ',mass.to('Msun')/gal.mass['gas'].to('Msun'))
print(pdf_handler_to_stats(obj,glob_attrs.result[0],3).to('Msun'))

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

