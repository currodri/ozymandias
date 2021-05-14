import numpy as np
import h5py
import os
import ozy
from ozy.dict_variables import common_variables,grid_variables,particle_variables,get_code_units
# TODO: Allow for parallel computation of profiles.
from joblib import Parallel, delayed
import sys
sys.path.append('/mnt/extraspace/currodri/Codes/ozymandias/ozy/amr')
sys.path.append('/mnt/extraspace/currodri/Codes/ozymandias/ozy/part')
from amr2 import vectors
from amr2 import geometrical_regions as geo
from amr2 import filtering
from amr2 import amr_profiles as amrprofmod
from part2 import part_profiles as partprofmod
from ozy.saver import _write_attrib

blacklist = [
    'yvars','weightvars','data'
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
    
    def _serialise(self, hdd):
        """This makes possible to save the group profile attrs as dataset attributes of an HDF5 group."""
        from yt import YTArray
        for k,v in self.__dict__.items():
            if k in blacklist:
                continue
            if isinstance(v, YTArray):
                hdd.attrs.create(k, v.d)
            elif isinstance(v, (int, float, bool, np.number)):
                hdd.attrs.create(k, v)
            elif isinstance(v, str):
                hdd.attrs.create(k, v.encode('utf8'))
            elif isinstance(v,dict):
                for kd,vd in v.items():
                    if isinstance(vd, YTArray):
                        hdd.attrs.create(kd, vd.d)
                    elif isinstance(vd, (int, float, bool, np.number)):
                        hdd.attrs.create(kd, vd)
                    elif isinstance(vd, str):
                        hdd.attrs.create(kd, vd.encode('utf8'))
                    elif isinstance(vd, list):
                        hdd.create_dataset('conditions', data=vd, compression=1)
    def _get_python_region(self,reg):
        """Save the Fortran derived type as a dictionary inside the Profile class (only the necessary info)."""
        from yt import YTArray
        self.region = {}
        self.region['type'] = reg.name.decode().split(' ')[0]
        self.region['centre'] = YTArray([reg.centre.x, reg.centre.y, reg.centre.z], 'code_length', registry=self.group.obj.unit_registry)
        self.region['axis'] = YTArray([reg.axis.x, reg.axis.y, reg.axis.z], 'dimensionless', registry=self.group.obj.unit_registry)
    
    def _get_python_filter(self,filt):
        """Save the Fortran derived type as a dictionary inside the Profile class (only the necessary info)."""
        from yt import YTQuantity
        self.filter = {}
        self.filter['name'] = filt.name.decode().split(' ')[0]
        self.filter['conditions'] = []
        if self.filter['name'] != 'none':
            for i in range(0, filt.ncond):
                cond_var = filt.cond_vars.T.view('S128')[i][0].decode().split(' ')[0]
                cond_op = filt.cond_ops.T.view('S2')[i][0].decode().split(' ')[0]
                cond_units = get_code_units(cond_var)
                cond_value = YTQuantity(filt.cond_vals[i], str(cond_units), registry=self.group.obj.unit_registry)
                cond_str = cond_var+'/'+cond_op+'/'+str(cond_value.d)+'/'+cond_units
                self.filter['conditions'].append(cond_str)

def init_region(group, region_type):
    """Initialise region Fortran derived type with details of group."""
    reg = geo.region()

    if region_type == 'sphere':
        reg.name = 'sphere'
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0], group.position[1], group.position[2]
        reg.centre = centre
        axis = vectors.vector()
        norm_L = group.angular_mom['total']/np.linalg.norm(group.angular_mom['total'])
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        reg.axis = axis
        bulk = vectors.vector()
        bulk.x, bulk.y, bulk.z = group.velocity[0], group.velocity[1], group.velocity[2]
        reg.bulk_velocity = bulk
        reg.rmin = 0.0
        # Basic configuration: 0.2 of the virial radius of the host halo
        reg.rmax = 0.2*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
    else:
        raise KeyError('Region type not supported. Please check!')
    return reg

def init_filter(cond_strs, name, group):
    """Initialise filter Fortran derived type with the condition strings provided."""
    from yt import YTQuantity
    if isinstance(cond_strs, str):
        cond_strs = [cond_strs]
    filt = filtering.filter()
    if cond_strs[0] == None:
        filt.ncond = 0
        filt.name = 'none'
        return filt
    elif name != 'none':
        filt.ncond = len(cond_strs)
        filt.name = name
        filtering.allocate_filter(filt)
        for i in range(0, filt.ncond):
            # Variable name
            filt.cond_vars.T.view('S128')[i] = cond_strs[i].split('/')[0].ljust(128)
            # Expresion operator
            filt.cond_ops.T.view('S2')[i] = cond_strs[i].split('/')[1].ljust(2)
            # Value transformed to code units
            value = YTQuantity(float(cond_strs[i].split('/')[2]), cond_strs[i].split('/')[3], registry=group.obj.unit_registry)
            filt.cond_vals[i] = value.in_units(get_code_units(cond_strs[i].split('/')[0])).d
        return filt
    else:
        raise ValueError("Condition strings are given, but a name for the filter. Please set!")

def compute_profile(group,ozy_file,xvar,yvars,weightvars,lmax=0,nbins=100,region_type='sphere',filter_conds='none',
                    filter_name='none',recompute=False,save=False):
    """Function which computes a 1D profile for a given group object."""

    if not isinstance(xvar, str):
        print('Single x variable 1D profile supported!')
        exit
    prof = Profile(group)
    prof.nbins = nbins
    prof.xvar = xvar
    prof.yvars = dict(for_hydro = [],for_star = [], for_dm = [])
    prof.weightvars = dict(for_hydro = [],for_star = [], for_dm = [])
    # Begin by checking the combination of variables is correct
    # Variables should be given in the form "type"/"variables":
    # e.g. gas/density, or star/ang_momentum_x
    for var in yvars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            if var_name in common_variables or var_name in grid_variables:
                prof.yvars['for_hydro'].append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!: %s', var)
        elif var_type == 'star':
            if var_name in common_variables or var_name in particle_variables:
                prof.yvars['for_star'].append(var_name)
            else:
                raise KeyError('This star variable is not supported. Please check!')
        elif var_type == 'dm':
            if var_name in common_variables or var_name in particle_variables:
                prof.yvars['for_dm'].append(var_name)
            else:
                raise KeyError('This DM variable is not supported. Please check!')
    for var in weightvars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            if var_name in common_variables or var_name in grid_variables:
                prof.weightvars['for_hydro'].append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!')
        elif var_type == 'star':
            if var_name in common_variables or var_name in particle_variables:
                prof.weightvars['for_star'].append(var_name)
            else:
                raise KeyError('This star variable is not supported. Please check!')
        elif var_type == 'dm':
            if var_name in common_variables or var_name in particle_variables:
                prof.weightvars['for_dm'].append(var_name)
            else:
                raise KeyError('This DM variable is not supported. Please check!')
    
    # Check that we do not have any inconsistency...
    if xvar in grid_variables and len(prof.yvars['for_star'])>0 or xvar in grid_variables and len(prof.yvars['for_dm'])>0:
        raise KeyError("Having grid vs particle 1D profiles is not well-defined.")
    elif xvar in particle_variables and len(prof.yvars['for_hydro'])>0:
        raise KeyError("Having particle vs grid 1D profiles is not well-defined.")
    
    # Now create region
    if isinstance(region_type, geo.region):
        selected_reg = region_type
    else:
        selected_reg = init_region(group,region_type)
    
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
        del f[str(prof.group.obj_type)+'_data/profiles/'+str(group._index)+'/'+str(prof_key)]
        print('Overwriting profile data in %s_data'%group.obj_type)
    elif prof_present and not recompute:
        print('Profile data with same details already present for galaxy %s. No overwritting!'%group._index)
        group._init_profiles()
        for i,p in enumerate(group.profiles):
            if p.key == prof_key:
                selected_prof = i
                break
        return group.profiles[selected_prof]
    else:
        print('Writing profile data in %s_data'%group.obj_type)
    f.close()
    
    # Initialise hydro profile data object
    hydro_data = amrprofmod.profile_handler()
    hydro_data.profdim = 1
    hydro_data.xvarname = xvar
    hydro_data.nyvar = len(prof.yvars['for_hydro'])
    hydro_data.nwvar = len(prof.weightvars['for_hydro'])
    hydro_data.nbins = 100

    amrprofmod.allocate_profile_handler(hydro_data)
    for i in range(0, len(prof.yvars['for_hydro'])):
        hydro_data.yvarnames.T.view('S128')[i] = prof.yvars['for_hydro'][i].ljust(128)
    for i in range(0, len(prof.weightvars['for_hydro'])):
        hydro_data.wvarnames.T.view('S128')[i] = prof.weightvars['for_hydro'][i].ljust(128)
    
    # And now, compute hydro data profiles!
    if hydro_data.nyvar > 0 and hydro_data.nwvar > 0:
        amrprofmod.onedprofile(group.obj.simulation.fullpath,selected_reg,filt,hydro_data,lmax)

    # Initialise particles profile data object
    star_data = partprofmod.profile_handler()
    star_data.profdim = 1
    star_data.xvarname = xvar
    star_data.nyvar = len(prof.yvars['for_star'])
    star_data.nwvar = len(prof.weightvars['for_star'])
    star_data.nbins = nbins

    partprofmod.allocate_profile_handler(star_data)
    for i in range(0, len(prof.yvars['for_star'])):
        tempstr = 'star/'+prof.yvars['for_star'][i]
        star_data.yvarnames.T.view('S128')[i] = tempstr.ljust(128)
    for i in range(0, len(prof.weightvars['for_star'])):
        tempstr = 'star/'+prof.weightvars['for_star'][i]
        star_data.wvarnames.T.view('S128')[i] = tempstr.ljust(128)

    # And now, compute star data profiles!
    if star_data.nyvar > 0 and star_data.nwvar > 0:
        partprofmod.onedprofile(group.obj.simulation.fullpath,selected_reg,filt,star_data,lmax)

    # And the same for dm particles
    dm_data = partprofmod.profile_handler()
    dm_data.profdim = 1
    dm_data.xvarname = xvar
    dm_data.nyvar = len(prof.yvars['for_dm'])
    dm_data.nwvar = len(prof.weightvars['for_dm'])
    dm_data.nbins = nbins

    partprofmod.allocate_profile_handler(dm_data)
    for i in range(0, len(prof.yvars['for_dm'])):
        tempstr = 'dn/'+prof.yvars['for_dm'][i]
        dm_data.yvarnames.T.view('S128')[i] = tempstr.ljust(128)
    for i in range(0, len(prof.weightvars['for_dm'])):
        tempstr = 'star/'+prof.weightvars['for_dm'][i]
        dm_data.wvarnames.T.view('S128')[i] = tempstr.ljust(128)

    # And now, compute star data profiles!
    if dm_data.nyvar > 0 and dm_data.nwvar > 0:
        partprofmod.onedprofile(group.obj.simulation.fullpath,selected_reg,filt,dm_data,lmax)

    if save:
        write_profiles(group.obj, ozy_file, hydro_data, star_data, dm_data, prof)
    
    return hydro_data, star_data, dm_data

def check_if_same_profile(hd, profile):
    """This function checks if a profile for an object already exists with the same attributes."""
    if not str(profile.group.obj_type)+'_data/profiles/'+str(profile.group._index) in hd:
        return False, 'none'
    for p in hd[str(profile.group.obj_type)+'_data/profiles/'+str(profile.group._index)].keys():
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
        profiles = f.create_group(str(prof.group.obj_type)+'_data/profiles/'+str(prof.group._index))
    except:
        profiles = f[str(prof.group.obj_type)+'_data/profiles/'+str(prof.group._index)]
    
    # Clean data and save to dataset
    prof_name = get_profile_name(profiles, prof)
    hdprof = profiles.create_group(prof_name)
    prof._serialise(hdprof)
    # Save x data
    xdata = np.zeros((3,prof.nbins))
    if hydro != None:
        xdata[0,:] = hydro.xdata
    if star != None:
        xdata[1,:] = star.xdata
    if dm != None:
        xdata[2,:] = dm.xdata
    hdprof.create_dataset('xdata', data=xdata)
    hdprof['xdata'].attrs.create('units', get_code_units(prof.xvar))
    # Save hydro y data
    if hydro != None:
        clean_hydro = hdprof.create_group('hydro')
        for v,var in enumerate(prof.yvars['for_hydro']):
            clean_hydro.create_dataset(var, data=hydro.ydata[:,v,:,0:2])
            clean_hydro[var].attrs.create('units', get_code_units(prof.yvars['for_hydro'][v]))
            clean_hydro[var].attrs.create('weightvars', prof.weightvars['for_hydro'][:])
    # Save star y data
    if star != None:
        clean_star = hdprof.create_group('star')
        for v,var in enumerate(prof.yvars['for_star']):
            clean_star.create_dataset(var, data=star.ydata[:,v,:,0:2])
            clean_star[var].attrs.create('units', get_code_units(prof.yvars['for_star'][v]))
            clean_star[var].attrs.create('weightvars', prof.weightvars['for_star'][:])
    # Save dm y data
    if dm != None:
        clean_dm = hdprof.create_group('dm')
        for v,var in enumerate(prof.yvars['for_dm']):
            clean_dm.create_dataset(var, data=dm.ydata[:,v,:,0:2])
            clean_dm[var].attrs.create('units', get_code_units(prof.yvars['for_dm'][v]))
            clean_dm[var].attrs.create('weightvars', prof.weightvars['for_dm'][:])
    f.close()
    return


    