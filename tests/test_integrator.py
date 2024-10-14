import numpy as np
import ozy
import sys
from amr2 import amr_integrator,stats_utils
from amr2 import io_ramses
from part2 import part_integrator
from ozy.utils import init_region,init_filter,\
                        structure_regions,get_code_bins,\
                        pdf_handler_to_stats
from ozy.dict_variables import get_code_units

io_ramses.verbose = True
obj = ozy.load('/mnt/extraspace/currodri/NUT/cosmoNUTmhd/Groups/ozy_00016.hdf5')#('test_00035.hdf5')
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

# quantity_names = ['mass','density','temperature',
#                             'momentum_sphere_r','v_sphere_r','v_tangential',
#                             'thermal_energy','thermal_energy_specific',
#                             'grav_therpfrsphere','grav_therpfrspherepos',
#                             'grav_therpfrsphereneg']
# do_binning = [True,True,True,
#                         True,True,True,
#                         True,True,
#                         True,True,
#                         True]
# weight_names = ['cumulative','mass','volume']
        
# # Initialise Fortran derived type with attributes
# # This object hold the following attributes:
# # - nvars: number of variables
# # - nwvars:  number of variables for weighting
# # - varnames: names of variables
# # - wvarnames:  names of variables for weighting
# # - data: organised in numpy array of shape (nvars,nwvars,4)
# #           each of those 4 values are (final, min, max, sum of weights)
# glob_attrs = amr_integrator.amr_region_attrs()
# glob_attrs.nvars = len(quantity_names)
# glob_attrs.nwvars = len(weight_names)
# glob_attrs.nfilter = 4
# glob_attrs.nsubs = nsubs
# amr_integrator.allocate_amr_regions_attrs(glob_attrs)
# for i in range(0,len(quantity_names)):
#     glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
# glob_attrs.result.nbins = 100
# glob_attrs.result.nvars = len(quantity_names)
# glob_attrs.result.nwvars = len(weight_names)
# glob_attrs.result.nfilter = 4
# stats_utils.allocate_pdf(glob_attrs.result)
# for i in range(0, len(quantity_names)):
#     mybins = get_code_bins(obj,'gas/'+quantity_names[i],100)
#     print(quantity_names[i],mybins)
#     glob_attrs.result.varname.T.view('S128')[i] = quantity_names[i].ljust(128)
#     glob_attrs.result.scaletype.T.view('S128')[i] = mybins[1].ljust(128)
#     glob_attrs.result.bins[:,i] = mybins[0]
#     glob_attrs.result.do_binning[i] = do_binning[i]
#     glob_attrs.result.zero_index[i] = mybins[2]
#     glob_attrs.result.linthresh[i] = mybins[3]
# for i in range(0, len(weight_names)):
#     glob_attrs.result.wvarnames.T.view('S128')[i] = weight_names[i].ljust(128)
       
# for i in range(0,nsubs):
#     glob_attrs.subs[i] = subs[i]

# glob_attrs.filters[0] = init_filter(cond_strs=['entropy_specific/</4.4e+8/erg*K**-1*g**-1'],name='cold',group=gal)
# glob_attrs.filters[1] = init_filter(cond_strs=['entropy_specific/>=/4.4e+8/erg*K**-1*g**-1','entropy_specific/<=/23.2e+8/erg*K**-1*g**-1'],name='warm',group=gal)
# glob_attrs.filters[2] = init_filter(cond_strs=['entropy_specific/>/23.2e+8/erg*K**-1*g**-1'],name='hot',group=gal)
# glob_attrs.filters[3] = init_filter(cond_strs=['none'],name='all',group=gal)
# # Begin integration

# amr_integrator.integrate_region(output_path,selected_reg,True,True,glob_attrs)
# mass_names = ['COLD','WARM','HOT']
# tot_mass = 0
# for i in range(0,3):
#     varname = str(glob_attrs.result.varname.T.view('S128')[0][0].decode()).rstrip()
#     unit = get_code_units(varname)
#     mass = glob_attrs.result.totweights[0,i,0]
#     mass = obj.quantity(mass,unit)
#     tot_mass += mass
#     print(mass_names[i]+': ',mass.to('Msun'))
# print('Total mass: ', tot_mass.to('Msun'))
# mass = glob_attrs.result.totweights[0,3,1]
# mass = obj.quantity(mass,'code_mass')
# print('Total mass no filter: ',mass.to('Msun'))
# print('From Ozymandias: ',gal.mass['gas'].to('Msun'))
# print('Ratio filter: ',tot_mass.to('Msun')/gal.mass['gas'].to('Msun'))
# print('Ratio no filter: ',mass.to('Msun')/gal.mass['gas'].to('Msun'))
# #print(pdf_handler_to_stats(obj,glob_attrs.result[0],3).to('Msun'))
# #print(pdf_handler_to_stats(obj,glob_attrs.result[4],3))
# #print(glob_attrs.result[4])
# print(pdf_handler_to_stats(obj,glob_attrs.result,3,3))

quantity_names = ['star/mass','star/sfr_10','star/sfr_100','star/metallicity',
                  'star/ang_momentum_x','star/ang_momentum_y','star/ang_momentum_z',
                  'dm/ang_momentum_x','dm/ang_momentum_y','dm/ang_momentum_z']
do_binning = [False,False,False,
              False,False,False,
              False,False,False]
weight_names = ['cumulative']

glob_attrs = part_integrator.part_region_attrs()
glob_attrs.nvars = len(quantity_names)
glob_attrs.nwvars = len(weight_names)
glob_attrs.nfilter = 1
glob_attrs.nsubs = nsubs
part_integrator.allocate_part_regions_attrs(glob_attrs)
for i in range(0,len(quantity_names)):
    glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
glob_attrs.result.nbins = 100
glob_attrs.result.nvars = len(quantity_names)
glob_attrs.result.nwvars = len(weight_names)
glob_attrs.result.nfilter = 4
stats_utils.allocate_pdf(glob_attrs.result)

for i in range(0, len(quantity_names)):
    mybins = get_code_bins(obj,quantity_names[i],100)
    print(quantity_names[i],mybins)
    glob_attrs.result.varname.T.view('S128')[i] = quantity_names[i].ljust(128)
    glob_attrs.result.scaletype.T.view('S128')[i] = mybins[1].ljust(128)
    glob_attrs.result.bins[:,i] = mybins[0]
    glob_attrs.result.do_binning[i] = do_binning[i]
    glob_attrs.result.zero_index[i] = mybins[2]
    glob_attrs.result.linthresh[i] = mybins[3]
for i in range(0, len(weight_names)):
    glob_attrs.result.wvarnames.T.view('S128')[i] = weight_names[i].ljust(128)

for i in range(0,nsubs):
    glob_attrs.subs[i] = subs[i]

glob_attrs.filters[0] = init_filter(cond_strs=['none'],name='all',group=gal)

glob_attrs.varnames.T.view('S128')[0] = b'star/mass'.ljust(128)

# Begin integration
part_integrator.integrate_region(output_path,selected_reg,glob_attrs,False,True)

