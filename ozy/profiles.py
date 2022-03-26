import os
import sys

import h5py
import numpy as np
from unyt import unyt_array,unyt_quantity
# TODO: Allow for parallel computation of profiles.
import ozy
from ozy.utils import init_region,init_filter
from ozy.dict_variables import (common_variables, get_code_units,
                                grid_variables, particle_variables)

sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/part')
from amr2 import amr_profiles as amrprofmod
from amr2 import geometrical_regions as geo
from part2 import part_profiles as partprofmod

blacklist = [
    'yvars','weightvars','data','xdata','ydata'
]
class Profile(object):

    def __init__(self,group):
        self.group = group
        self.nbins = 0
        self.xvar = None
        self.region = None
        self.filter = None
        self.lmax = 0
        self.yvars = {}
        self.weightvars = {}
        self.xdata = None
        self.ydata = {}
    
    def _serialise(self, hdd):
        """This makes possible to save the group profile attrs as dataset attributes of an HDF5 group."""
        for k,v in self.__dict__.items():
            if k in blacklist:
                continue
            if isinstance(v, unyt_array):
                hdd.attrs.create(k, v.d)
            elif isinstance(v, (int, float, bool, np.number)):
                hdd.attrs.create(k, v)
            elif isinstance(v, str):
                hdd.attrs.create(k, v.encode('utf8'))
            elif isinstance(v,dict):
                for kd,vd in v.items():
                    if isinstance(vd, unyt_array):
                        hdd.attrs.create(kd, vd.d)
                    elif isinstance(vd, (int, float, bool, np.number)):
                        hdd.attrs.create(kd, vd)
                    elif isinstance(vd, str):
                        hdd.attrs.create(kd, vd.encode('utf8'))
                    elif isinstance(vd, list):
                        hdd.create_dataset('conditions', data=vd, compression=1)
    def _get_python_region(self,reg):
        """Save the Fortran derived type as a dictionary inside the Profile class (only the necessary info)."""
        self.region = {}
        self.region['type'] = reg.name.decode().split(' ')[0]
        self.region['centre'] = self.group.obj.array([reg.centre.x, reg.centre.y, reg.centre.z], 'code_length')
        self.region['axis'] = self.group.obj.array([reg.axis.x, reg.axis.y, reg.axis.z], 'dimensionless')
    
    def _get_python_filter(self,filt):
        """Save the Fortran derived type as a dictionary inside the Profile class (only the necessary info)."""
        self.filter = {}
        self.filter['name'] = filt.name.decode().split(' ')[0]
        self.filter['conditions'] = []
        if self.filter['name'] != 'none':
            for i in range(0, filt.ncond):
                cond_var = filt.cond_vars.T.view('S128')[i][0].decode().split(' ')[0]
                cond_op = filt.cond_ops.T.view('S2')[i][0].decode().split(' ')[0]
                cond_units = get_code_units(cond_var)
                cond_value = self.group.obj.quantity(filt.cond_vals[i], str(cond_units))
                cond_str = cond_var+'/'+cond_op+'/'+str(cond_value.d)+'/'+cond_units
                self.filter['conditions'].append(cond_str)



def compute_profile(group,ozy_file,xvar,yvars,weightvars,lmax=0,nbins=100,region_type='sphere',filter_conds='none',
                    filter_name='none',recompute=False,save=False,logscale=False,
                    rmin=(0.0,'rvir'), rmax=(0.2,'rvir'), zmin=(0.0,'rvir'), zmax=(0.2,'rvir')):
    """Function which computes a 1D profile for a given group object."""

    if not isinstance(xvar, str):
        print('Single x variable 1D profile supported!')
        exit
    prof = Profile(group)
    prof.nbins = nbins
    prof.xvar = xvar
    prof.yvars = dict(hydro = [],star = [], dm = [])
    prof.weightvars = dict(hydro = [],star = [], dm = [])
    # Begin by checking the combination of variables is correct
    # Variables should be given in the form "type"/"variables":
    # e.g. gas/density, or star/ang_momentum_x
    for var in yvars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            if var_name in common_variables or var_name in grid_variables:
                prof.yvars['hydro'].append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!: %s', var)
        elif var_type == 'star':
            if var_name in common_variables or var_name in particle_variables:
                prof.yvars['star'].append(var_name)
            else:
                raise KeyError('This star variable is not supported. Please check!')
        elif var_type == 'dm':
            if var_name in common_variables or var_name in particle_variables:
                prof.yvars['dm'].append(var_name)
            else:
                raise KeyError('This DM variable is not supported. Please check!')
    for var in weightvars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            if var_name in common_variables or var_name in grid_variables:
                prof.weightvars['hydro'].append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!')
        elif var_type == 'star':
            if var_name in common_variables or var_name in particle_variables:
                prof.weightvars['star'].append(var_name)
            else:
                raise KeyError('This star variable is not supported. Please check!')
        elif var_type == 'dm':
            if var_name in common_variables or var_name in particle_variables:
                prof.weightvars['dm'].append(var_name)
            else:
                raise KeyError('This DM variable is not supported. Please check!')
    
    # Check that we do not have any inconsistency...
    if xvar in grid_variables and len(prof.yvars['star'])>0 or xvar in grid_variables and len(prof.yvars['dm'])>0:
        raise KeyError("Having grid vs particle 1D profiles is not well-defined.")
    elif xvar in particle_variables and len(prof.yvars['hydro'])>0:
        raise KeyError("Having particle vs grid 1D profiles is not well-defined.")
    
    # Now create region
    if isinstance(region_type, geo.region):
        selected_reg = region_type
    else:
        selected_reg = init_region(group,region_type,rmin=rmin,rmax=rmax,zmin=zmin,zmax=zmax)
    
    # Save region details to profile object
    prof._get_python_region(selected_reg)

    # Now create filter, if any conditions have been given…
    filt = init_filter(filter_conds, filter_name, group)

    # And save to profile object
    prof._get_python_filter(filt)

    # Check if profile data is already present and if it coincides with the new one
    f = h5py.File(ozy_file, 'r+')
    prof_present,prof_key = check_if_same_profile(f, prof)
    if prof_present and recompute:
        del f[str(prof.group.type)+'_data/profiles/'+str(group._index)+'/'+str(prof_key)]
        print('Overwriting profile data in %s_data'%group.type)
    elif prof_present and not recompute:
        print('Profile data with same details already present for galaxy %s. No overwritting!'%group._index)
        group._init_profiles()
        for i,p in enumerate(group.profiles):
            if p.key == prof_key:
                selected_prof = i
                break
        return group.profiles[selected_prof]
    elif save and recompute:
        print('Writing profile data in %s_data'%group.type)
    f.close()
    
    # Initialise hydro profile data object
    if len(prof.yvars['hydro'])>0 and len(prof.weightvars['hydro'])>0:
        hydro_data = amrprofmod.profile_handler()
        hydro_data.profdim = 1
        hydro_data.xvarname = xvar
        hydro_data.nyvar = len(prof.yvars['hydro'])
        hydro_data.nwvar = len(prof.weightvars['hydro'])
        hydro_data.nbins = nbins

        amrprofmod.allocate_profile_handler(hydro_data)
        for i in range(0, len(prof.yvars['hydro'])):
            hydro_data.yvarnames.T.view('S128')[i] = prof.yvars['hydro'][i].ljust(128)
        for i in range(0, len(prof.weightvars['hydro'])):
            hydro_data.wvarnames.T.view('S128')[i] = prof.weightvars['hydro'][i].ljust(128)
        
        # And now, compute hydro data profiles!
        if hydro_data.nyvar > 0 and hydro_data.nwvar > 0:
            amrprofmod.onedprofile(group.obj.simulation.fullpath,selected_reg,filt,hydro_data,lmax,logscale)
    else:
        hydro_data = None
    
    # Initialise particles profile data object
    if len(prof.yvars['star'])>0 and len(prof.weightvars['star']):
        star_data = partprofmod.profile_handler()
        star_data.profdim = 1
        star_data.xvarname = 'star/'+xvar
        star_data.nyvar = len(prof.yvars['star'])
        star_data.nwvar = len(prof.weightvars['star'])
        star_data.nbins = nbins

        partprofmod.allocate_profile_handler(star_data)
        for i in range(0, len(prof.yvars['star'])):
            tempstr = 'star/'+prof.yvars['star'][i]
            star_data.yvarnames.T.view('S128')[i] = tempstr.ljust(128)
        for i in range(0, len(prof.weightvars['star'])):
            tempstr = 'star/'+prof.weightvars['star'][i]
            star_data.wvarnames.T.view('S128')[i] = tempstr.ljust(128)

        # And now, compute star data profiles!
        if star_data.nyvar > 0 and star_data.nwvar > 0:
            partprofmod.onedprofile(group.obj.simulation.fullpath,selected_reg,filt,star_data,lmax)
    else:
        star_data = None

    
    # And the same for dm particles
    if len(prof.yvars['dm'])>0 and len(prof.yvars['dm'])>0:
        dm_data = partprofmod.profile_handler()
        dm_data.profdim = 1
        dm_data.xvarname = 'dm/'+xvar
        dm_data.nyvar = len(prof.yvars['dm'])
        dm_data.nwvar = len(prof.yvars['dm'])
        dm_data.nbins = nbins

        partprofmod.allocate_profile_handler(dm_data)
        for i in range(0, len(prof.yvars['dm'])):
            tempstr = 'dm/'+prof.yvars['dm'][i]
            dm_data.yvarnames.T.view('S128')[i] = tempstr.ljust(128)
        for i in range(0, len(prof.weightvars['dm'])):
            tempstr = 'dm/'+prof.weightvars['dm'][i]
            dm_data.wvarnames.T.view('S128')[i] = tempstr.ljust(128)

        # And now, compute dm data profiles!
        if dm_data.nyvar > 0 and dm_data.nwvar > 0:
            partprofmod.onedprofile(group.obj.simulation.fullpath,selected_reg,filt,dm_data,lmax)
    else:
        dm_data = None

    

    # Organise everything in the Profile object
    xdata = np.zeros((3,prof.nbins))
    if hydro_data != None:
        xdata[0,:] = hydro_data.xdata[1:]
    if star_data != None:
        xdata[1,:] = star_data.xdata[1:]
    if dm_data != None:
        xdata[2,:] = dm_data.xdata[1:]
    prof.xdata = group.obj.array(xdata, get_code_units(prof.xvar))
    # Save hydro y data
    if hydro_data != None:
        prof.ydata['hydro'] = []
        for v,var in enumerate(prof.yvars['hydro']):
            arr = np.copy(hydro_data.ydata[:,v,:,0:2])
            prof.ydata['hydro'].append(group.obj.array(arr, get_code_units(prof.yvars['hydro'][v])))
    # Save star y data
    if star_data != None:
        prof.ydata['star'] = []
        for v,var in enumerate(prof.yvars['star']):
            arr = np.copy(star_data.ydata[:,v,:,0:2])
            prof.ydata['star'].append(group.obj.array(arr, get_code_units(prof.yvars['star'][v])))
    # Save dm y data
    if dm_data != None:
        prof.ydata['dm'] = []
        for v,var in enumerate(prof.yvars['dm']):
            arr = np.copy(dm_data.ydata[:,v,:,0:2])
            prof.ydata['dm'].append(group.obj.array(arr, get_code_units(prof.yvars['dm'][v])))
    if save:
        write_profiles(group.obj, ozy_file, hydro_data, star_data, dm_data, prof)
    return prof

def check_if_same_profile(hd, profile):
    """This function checks if a profile for an object already exists with the same attributes."""
    if not str(profile.group.type)+'_data/profiles/'+str(profile.group._index) in hd:
        return False, 'none'
    for p in hd[str(profile.group.type)+'_data/profiles/'+str(profile.group._index)].keys():
        check_xvar = (p.split('|')[0] == profile.xvar)
        check_filtername = (p.split('|')[1] == profile.filter['name'])
        check_regiontype = (p.split('|')[2] == profile.region['type'])
        if check_xvar and check_filtername and check_regiontype:
            return True, p
    return False, 'none'

def get_profile_name(profiles_group,prof):
    """Create an individual profile identifier name."""
    name = str(prof.xvar)
    name += '|'+str(prof.filter['name'])
    name += '|'+str(prof.region['type'])
    if name in profiles_group:
        name += '|new'
    print(name)
    return name
def write_profiles(obj, ozy_file, hydro, star, dm, prof):
    """This function writes the resulting profile data for this group to the original OZY HDF5 file."""

    f = h5py.File(ozy_file, 'r+')

    # Create group in HDF5 file
    try:
        profiles = f.create_group(str(prof.group.type)+'_data/profiles/'+str(prof.group._index))
    except:
        profiles = f[str(prof.group.type)+'_data/profiles/'+str(prof.group._index)]
    
    # Clean data and save to dataset
    prof_name = get_profile_name(profiles, prof)
    hdprof = profiles.create_group(prof_name)
    prof._serialise(hdprof)
    # Save x data
    xdata = np.zeros((3,prof.nbins))
    if hydro != None:
        xdata[0,:] = hydro.xdata[1:]
    if star != None:
        xdata[1,:] = star.xdata[1:]
    if dm != None:
        xdata[2,:] = dm.xdata[1:]
    hdprof.create_dataset('xdata', data=xdata)
    hdprof['xdata'].attrs.create('units', get_code_units(prof.xvar))
    # Save hydro y data
    if hydro != None:
        clean_hydro = hdprof.create_group('hydro')
        for v,var in enumerate(prof.yvars['hydro']):
            clean_hydro.create_dataset(var, data=hydro.ydata[:,v,:,0:2])
            clean_hydro[var].attrs.create('units', get_code_units(prof.yvars['hydro'][v]))
            clean_hydro[var].attrs.create('weightvars', prof.weightvars['hydro'][:])
    # Save star y data
    if star != None:
        clean_star = hdprof.create_group('star')
        for v,var in enumerate(prof.yvars['star']):
            clean_star.create_dataset(var, data=star.ydata[:,v,:,0:2])
            clean_star[var].attrs.create('units', get_code_units(prof.yvars['star'][v]))
            clean_star[var].attrs.create('weightvars', prof.weightvars['star'][:])
    # Save dm y data
    if dm != None:
        clean_dm = hdprof.create_group('dm')
        for v,var in enumerate(prof.yvars['dm']):
            clean_dm.create_dataset(var, data=dm.ydata[:,v,:,0:2])
            clean_dm[var].attrs.create('units', get_code_units(prof.yvars['dm'][v]))
            clean_dm[var].attrs.create('weightvars', prof.weightvars['dm'][:])
    f.close()
    return


    
