import os
import sys

import h5py
import numpy as np
from unyt import unyt_array,unyt_quantity
# TODO: Allow for parallel computation of profiles.
import ozy
from ozy.utils import init_region,init_filter,pdf_handler_to_stats
from ozy.dict_variables import (common_variables, get_code_units,
                                grid_variables, particle_variables)
from amr2 import amr_profiles as amrprofmod
from amr2 import stats_utils
from amr2 import geometrical_regions as geo
from part2 import part_profiles as partprofmod

blacklist = [
    'yvars','weightvars','data','xdata','ydata'
]
class Profile(object):

    def __init__(self,group):
        self.obj = group.obj
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
        self.rm_subs = False
    
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
        """Save the Fortran derived type as a dictionary inside the GalacticFlow class (only the necessary info)."""
        self.region = {}
        self.region['type'] = reg.name.decode().split(' ')[0]
        self.region['centre'] = self.obj.array([reg.centre.x, reg.centre.y, reg.centre.z], 'code_length')
        self.region['axis'] = self.obj.array([reg.axis.x, reg.axis.y, reg.axis.z], 'dimensionless')
        self.region['rmin'] = self.obj.quantity(reg.rmin, 'code_length')
        self.region['rmax'] = self.obj.quantity(reg.rmax, 'code_length')
    
    def _get_python_filter(self,filt):
        """Save the Fortran derived type as a dictionary inside the Profile class (only the necessary info)."""
        self.filter = {}
        self.filter['name'] = filt.name.decode().split(' ')[0]
        self.filter['conditions'] = []
        if self.filter['name'] != 'none':
            for i in range(0, filt.ncond):
                if filt.use_var[i] == 1:
                    cond_var = filt.cond_vars.T.view('S128')[i][0].decode().split(' ')[0]
                    cond_op = filt.cond_ops.T.view('S2')[i][0].decode().split(' ')[0]
                    cond_vars_comp = filt.cond_vars_comp.T.view('S128')[i][0].decode().split(' ')[0]
                    cond_vars_factor = filt.cond_vals[i]
                    cond_str = cond_var+'/'+cond_op+'/'+cond_vars_comp+'/'+str(cond_vars_factor)
                    self.filter['conditions'].append(cond_str)
                else:
                    cond_var = filt.cond_vars.T.view('S128')[i][0].decode().split(' ')[0]
                    cond_op = filt.cond_ops.T.view('S2')[i][0].decode().split(' ')[0]
                    cond_units = get_code_units(cond_var)
                    cond_value = self.group.obj.quantity(filt.cond_vals[i], str(cond_units))
                    cond_str = cond_var+'/'+cond_op+'/'+str(cond_value.d)+'/'+cond_units
                    self.filter['conditions'].append(cond_str)



def compute_profile(group,ozy_file,xvar,yvars,weightvars,minval,maxval,linthresh=None,lmax=0,nbins=100,
                    region_type='sphere',filter_conds=['none'],
                    filter_name=['none'],recompute=False,save=False,logscale=False,
                    pdf_bins=100,regime_type='',
                    rmin=(0.0,'rvir'), rmax=(0.2,'rvir'), zmin=(0.0,'rvir'), zmax=(0.2,'rvir'),
                    mycentre=([0.5,0.5,0.5],'rvir'), myaxis=np.array([1.,0.,0.]),
                    remove_subs=False,cr_st=False,cr_heat=False,Dcr=0.0,
                    verbose=False,force_neigh=False,force_read_gravity=False,
                    do_binning=True):
    """Function which computes a 1D profile for a given group object."""
    from ozy.plot_settings import plotting_dictionary
    from ozy.utils import structure_regions,get_code_bins
    from ozy.dict_variables import check_need_neighbours
    if not isinstance(xvar, str):
        if verbose: print('Single x variable 1D profile supported!')
        exit
    nfilter = len(filter_name)
    profs = []
    for i in range(0, nfilter):
        prof = Profile(group)
        prof.nbins = nbins
        remove_all = False
        if remove_subs =='all':
            remove_all = True
            remove_subs = True
        use_neigh=False
        prof.rm_subs = remove_subs
        prof.xvar = xvar
        prof.yvars = dict(hydro = [],star = [], dm = [])
        prof.weightvars = dict(hydro = [],star = [], dm = [])
        profs.append(prof)
    
    # Begin by checking the combination of variables is correct
    # Variables should be given in the form "type"/"variables":
    # e.g. gas/density, or star/ang_momentum_x
    for var in yvars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            use_neigh = check_need_neighbours(var_name) or use_neigh
            if var_name in common_variables or var_name in grid_variables:
                for i in range(0, nfilter):
                    profs[i].yvars['hydro'].append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!',var_name)
        elif var_type == 'star':
            if var_name in common_variables or var_name in particle_variables:
                for i in range(0, nfilter):
                    profs[i].yvars['star'].append(var_name)
            else:
                raise KeyError('This star variable is not supported. Please check!',var_name)
        elif var_type == 'dm':
            if var_name in common_variables or var_name in particle_variables:
                for i in range(0, nfilter):
                    profs[i].yvars['dm'].append(var_name)
            else:
                raise KeyError('This DM variable is not supported. Please check!',var_name)
    for var in weightvars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            use_neigh = check_need_neighbours(var_name) or use_neigh
            if var_name in common_variables or var_name in grid_variables:
                for i in range(0, nfilter):
                    profs[i].weightvars['hydro'].append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!',var_name)
        elif var_type == 'star':
            if var_name in common_variables or var_name in particle_variables:
                for i in range(0, nfilter):
                    profs[i].weightvars['star'].append(var_name)
            else:
                raise KeyError('This star variable is not supported. Please check!',var_name)
        elif var_type == 'dm':
            if var_name in common_variables or var_name in particle_variables:
                for i in range(0, nfilter):
                    profs[i].weightvars['dm'].append(var_name)
            else:
                raise KeyError('This DM variable is not supported. Please check!',var_name)
    
    # Check if the xvar needs neighbours
    use_neigh = check_need_neighbours(xvar) or use_neigh or force_neigh
    
    if use_neigh:
        if verbose: print('At least one variable needs neighbours!')
    
    # Check that we do not have any inconsistency...
    if xvar in grid_variables and len(profs[0].yvars['star'])>0 or xvar in grid_variables and len(profs[0].yvars['dm'])>0:
        raise KeyError("Having grid vs particle 1D profiles is not well-defined.")
    elif xvar in particle_variables and len(profs[0].yvars['hydro'])>0:
        raise KeyError("Having particle vs grid 1D profiles is not well-defined.")
    
    # Check that the xaxis min and max quantities have the units expected for that variable
    try:
        minval = minval.to(get_code_units(xvar))
        maxval = maxval.to(get_code_units(xvar))
    except:
        raise ValueError(f"It seems the dimensions of your bins min ({minval.units}) and max ({minval.units}) values do not agree with the dimensions of the chosen xvar ({xvar},{get_code_units(xvar)})")
    
    # Now create region
    
    if isinstance(region_type, geo.region):
        selected_reg = region_type
        enclosing_sphere_r = rmax
        enclosing_sphere_p = group.position
        if not np.array_equal(mycentre,group.position):
            enclosing_sphere_p = mycentre
    else:
        if not np.array_equal(mycentre, group.position) and not np.array_equal(myaxis,group.angular_mom['total']):
            selected_reg,enclosing_sphere_p,enclosing_sphere_r = init_region(group,region_type,rmin=rmin,
                                                                         rmax=rmax,zmin=zmin,zmax=zmax,
                                                                         mycentre=mycentre,myaxis=myaxis,
                                                                         return_enclosing_sphere=True)
        elif not np.array_equal(mycentre,group.position):
            selected_reg,enclosing_sphere_p,enclosing_sphere_r = init_region(group,region_type,rmin=rmin,
                                                                         rmax=rmax,zmin=zmin,zmax=zmax,
                                                                         mycentre=mycentre,
                                                                         return_enclosing_sphere=True)
        elif not np.array_equal(myaxis,group.angular_mom['total']):
            selected_reg,enclosing_sphere_p,enclosing_sphere_r = init_region(group,region_type,rmin=rmin,
                                                                         rmax=rmax,zmin=zmin,zmax=zmax,
                                                                         myaxis=myaxis,
                                                                         return_enclosing_sphere=True)
        else:
            selected_reg,enclosing_sphere_p,enclosing_sphere_r = init_region(group,region_type,rmin=rmin,
                                                                         rmax=rmax,zmin=zmin,zmax=zmax,
                                                                         return_enclosing_sphere=True)
        
    
    
    # Now create filters and regions for each profile
    filts = []
    for i in range(0,nfilter):
        filt = init_filter(filter_conds[i], filter_name[i], group)
        filts.append(filt)
        # Save region details to profile object
        prof = profs[i]
        prof._get_python_region(selected_reg)
        # And save to profile object
        prof._get_python_filter(filt)

    # Check if profile data is already present and if it coincides with the new one
    f = h5py.File(ozy_file, 'r+')
    profs_fr = []
    for i in range(0,nfilter):
        prof  = profs[i]
        prof_present,prof_key = check_if_same_profile(f, prof)
        if prof_present and recompute:
            if remove_subs:
                del f[str(prof.group.type)+'_data/profiles_nosubs/'+str(group._index)+'/'+str(prof_key)]
            else:
                del f[str(prof.group.type)+'_data/profiles/'+str(group._index)+'/'+str(prof_key)]
            profs_fr.append(True)
            if verbose: print('Overwriting profile data in %s_data'%group.type)
        elif prof_present and not recompute:
            if verbose: print('Profile data with same details already present for galaxy %s. No overwritting!'%group._index)
            group._init_profiles()
            if remove_subs:
                if verbose: print('Removing substructure!')
                for j,p in enumerate(group.profiles_nosubs):
                    if p.key == prof_key:
                        selected_prof = j
                        break
                profs_fr.append(group.profiles_nosubs[selected_prof])
            else:
                for j,p in enumerate(group.profiles):
                    if p.key == prof_key:
                        selected_prof = j
                        break
                profs_fr.append(group.profiles[selected_prof])
        elif save:
            profs_fr.append(True)
            if verbose: print('Writing profile data in %s_data'%group.type)
    f.close()
    
    nfilter_real = profs_fr.count(True)
    if nfilter_real == 0:
        if nfilter > 1:
            return profs_fr
        else:
            return profs_fr[0]
    
    # If substructre is removed, obtain regions
    remove_all = False
    if remove_subs == 'all':
        remove_all = True
        remove_subs = True
    if remove_all:
        subs = structure_regions(group, add_substructure=True, add_neighbours=False,
                                    add_intersections=True,position=enclosing_sphere_p,
                                    radius=enclosing_sphere_r,tidal_method='BT87_simple')
        nsubs = len(subs)
    elif remove_subs:
        if verbose: print('Removing substructure!')
        subs = structure_regions(group, add_substructure=True, add_neighbours=False,
                                    tidal_method='BT87_simple')
        nsubs = len(subs)
    else:
        nsubs = 0
        
    if regime_type != '':
        regime_type = '_'+regime_type
    
    # Initialise hydro profile data object
    if len(profs[0].yvars['hydro'])>0 and len(profs[0].weightvars['hydro'])>0:
        hydro_data = amrprofmod.profile_handler()
        hydro_data.profdim = 1
        hydro_data.xvarname = xvar
        hydro_data.nfilter = nfilter_real
        hydro_data.nyvar = len(profs[0].yvars['hydro'])
        hydro_data.nwvar = len(profs[0].weightvars['hydro'])
        hydro_data.nbins = nbins
        hydro_data.nsubs = nsubs
        hydro_data.cr_st = cr_st
        hydro_data.cr_heat = cr_heat
        hydro_data.Dcr = Dcr

        amrprofmod.allocate_profile_handler(hydro_data)
        for i in range(0, len(profs[0].yvars['hydro'])):
            hydro_data.yvarnames.T.view('S128')[i] = profs[0].yvars['hydro'][i].ljust(128)
        for i in range(0, len(profs[0].weightvars['hydro'])):
            hydro_data.wvarnames.T.view('S128')[i] = profs[0].weightvars['hydro'][i].ljust(128)
        
        # Add the scaletype for the xaxis and the pre-computed bin edges
        bin_edges, stype, zero_index, lint = get_code_bins(group.obj,'gas/'+xvar,nbins=nbins,logscale=logscale,
                                                    minval=minval,maxval=maxval,linthresh=linthresh)
        hydro_data.xdata = bin_edges
        hydro_data.scaletype = stype.ljust(128)
        hydro_data.linthresh = lint
        hydro_data.zero_index = zero_index
        
        if remove_subs and nsubs>0:
            for i in range(0,nsubs):
                hydro_data.subs[i] = subs[i]
        
        counter = 0
        for i in range(0,nfilter):
            if profs_fr[i] == True:
                hydro_data.filters[counter] = filts[i]
                counter += 1

        # Now add the PDF handler for each quantity and each xbin
        mybins = []
        for i in range(0, len(profs[0].yvars['hydro'])):
            plot_def = plotting_dictionary[profs[0].yvars['hydro'][i]]
            minv = group.obj.quantity(plot_def['vmin'+regime_type],plot_def['units'])
            maxv = group.obj.quantity(plot_def['vmax'+regime_type],plot_def['units'])
            mybins.append(get_code_bins(group.obj,'gas/'+profs[0].yvars['hydro'][i],pdf_bins,
                                        minval=minv,maxval=maxv))
            
        for j in range(0,nbins):
            hydro_data.ydata[j].nbins = pdf_bins
            hydro_data.ydata[j].nfilter = nfilter_real
            hydro_data.ydata[j].nvars = len(profs[0].yvars['hydro'])
            hydro_data.ydata[j].nwvars = len(profs[0].weightvars['hydro'])
            stats_utils.allocate_pdf(hydro_data.ydata[j])
            for k in range(0, len(profs[0].yvars['hydro'])):
                hydro_data.ydata[j].varname.T.view('S128')[k] = profs[0].yvars['hydro'][k].ljust(128)
                hydro_data.ydata[j].scaletype.T.view('S128')[k] = mybins[k][1].ljust(128)
                hydro_data.ydata[j].bins[:,k] = mybins[k][0]
                hydro_data.ydata[j].do_binning[k] = do_binning
                hydro_data.ydata[j].zero_index[k] = mybins[k][2]
                hydro_data.ydata[j].linthresh[k] = mybins[k][3]
            for k in range(0, len(profs[0].weightvars['hydro'])):
                hydro_data.ydata[j].wvarnames.T.view('S128')[k] = profs[0].weightvars['hydro'][k].ljust(128)
            
        # And now, compute hydro data profiles!
        read_gravity = force_read_gravity
        if hydro_data.nyvar > 0 and hydro_data.nwvar > 0:
            amrprofmod.onedprofile(group.obj.simulation.fullpath,selected_reg,hydro_data,lmax,use_neigh,read_gravity)
    else:
        hydro_data = None
    
    # Initialise particles profile data object
    if len(profs[0].yvars['star'])>0 and len(profs[0].weightvars['star']):
        star_data = partprofmod.profile_handler()
        star_data.profdim = 1
        star_data.xvarname = 'star/'+xvar
        star_data.nfilter = nfilter_real
        star_data.nyvar = len(profs[0].yvars['star'])
        star_data.nwvar = len(profs[0].weightvars['star'])
        star_data.nbins = nbins
        star_data.nsubs = nsubs

        partprofmod.allocate_profile_handler(star_data)
        for i in range(0, len(profs[0].yvars['star'])):
            tempstr = 'star/'+profs[0].yvars['star'][i]
            star_data.yvarnames.T.view('S128')[i] = tempstr.ljust(128)
        for i in range(0, len(profs[0].weightvars['star'])):
            tempstr = 'star/'+profs[0].weightvars['star'][i]
            star_data.wvarnames.T.view('S128')[i] = tempstr.ljust(128)
            
        # Add the scaletype for the xaxis and the pre-computed bin edges
        bin_edges, stype, zero_index, lint = get_code_bins(group.obj,'star/'+xvar,nbins=nbins,logscale=logscale,
                                         minval=minval,maxval=maxval,linthresh=linthresh)
        star_data.xdata = bin_edges
        star_data.scaletype = stype.ljust(128)
        star_data.linthresh = lint
        star_data.zero_index = zero_index
        
        if remove_subs and nsubs>0:
            for i in range(0,nsubs):
                star_data.subs[i] = subs[i]
        
        counter = 0
        for i in range(0,nfilter):
            if profs_fr[i] == True:
                star_data.filters[counter] = filts[i]
                counter += 1
        
        # Now add the PDF handler for each quantity and each xbin
        mybins = []
        for i in range(0, len(profs[0].yvars['star'])):
            tempstr = 'star_'+profs[0].yvars['star'][i]
            plot_def = plotting_dictionary[tempstr]
            minv = group.obj.quantity(plot_def['vmin'+regime_type],plot_def['units'])
            maxv = group.obj.quantity(plot_def['vmax'+regime_type],plot_def['units'])
            mybins.append(get_code_bins(group.obj,'star/'+profs[0].yvars['star'][i],pdf_bins,
                                        minval=minv,maxval=maxv))
            
        for j in range(0,nbins):
            star_data.ydata[j].nbins = pdf_bins
            star_data.ydata[j].nfilter = nfilter_real
            star_data.ydata[j].nvars = len(profs[0].yvars['star'])
            star_data.ydata[j].nwvars = len(profs[0].weightvars['star'])
            stats_utils.allocate_pdf(star_data.ydata[j])
            for k in range(0, len(profs[0].yvars['star'])):
                tempstr = 'star/'+profs[0].yvars['star'][k]
                star_data.ydata[j].varname.T.view('S128')[k] = tempstr.ljust(128)
                star_data.ydata[j].scaletype.T.view('S128')[k] = mybins[k][1].ljust(128)
                star_data.ydata[j].bins[:,k] = mybins[k][0]
                star_data.ydata[j].do_binning[k] = do_binning
                star_data.ydata[j].zero_index[k] = mybins[k][2]
                star_data.ydata[j].linthresh[k] = mybins[k][3]
            for k in range(0, len(profs[0].weightvars['star'])):
                tempstr = 'star/'+profs[0].weightvars['star'][k]
                star_data.ydata[j].wvarnames.T.view('S128')[k] = tempstr.ljust(128)
        
        # And now, compute star data profiles!
        if star_data.nyvar > 0 and star_data.nwvar > 0:
            partprofmod.onedprofile(group.obj.simulation.fullpath,selected_reg,star_data,lmax)
    else:
        star_data = None

    
    # And the same for dm particles
    if len(profs[0].yvars['dm'])>0 and len(profs[0].yvars['dm'])>0:
        dm_data = partprofmod.profile_handler()
        dm_data.profdim = 1
        dm_data.xvarname = 'dm/'+xvar
        dm_data.nfilter = nfilter_real
        dm_data.nyvar = len(profs[0].yvars['dm'])
        dm_data.nwvar = len(profs[0].yvars['dm'])
        dm_data.nbins = nbins
        dm_data.nsubs = nsubs

        partprofmod.allocate_profile_handler(dm_data)
        for i in range(0, len(profs[0].yvars['dm'])):
            tempstr = 'dm/'+profs[0].yvars['dm'][i]
            dm_data.yvarnames.T.view('S128')[i] = tempstr.ljust(128)
        for i in range(0, len(profs[0].weightvars['dm'])):
            tempstr = 'dm/'+profs[0].weightvars['dm'][i]
            dm_data.wvarnames.T.view('S128')[i] = tempstr.ljust(128)
            
        # Add the scaletype for the xaxis and the preo-computed bin edges
        bin_edges, stype, zero_index, lint = get_code_bins(group.obj,'dm/'+xvar,nbins=nbins,logscale=logscale,
                                         minval=minval,maxval=maxval,linthresh=linthresh)
        dm_data.xdata = bin_edges
        dm_data.scaletype = stype.ljust(128)
        dm_data.linthresh = lint
        dm_data.zero_index = zero_index
        
        if remove_subs and nsubs>0:
            for i in range(0,nsubs):
                dm_data.subs[i] = subs[i]
        
        counter = 0
        for i in range(0,nfilter):
            if profs_fr[i] == True:
                dm_data.filters[counter] = filts[i]
                counter += 1
                
        # Now add the PDF handler for each quantity and each xbin
        mybins = []
        for i in range(0, len(profs[0].yvars['dm'])):
            tempstr = 'dm_'+profs[0].yvars['dm'][i]
            plot_def = plotting_dictionary[tempstr]
            minv = group.obj.quantity(plot_def['vmin'+regime_type],plot_def['units'])
            maxv = group.obj.quantity(plot_def['vmax'+regime_type],plot_def['units'])
            mybins.append(get_code_bins(group.obj,'dm/'+profs[0].yvars['dm'][i],pdf_bins,
                                        minval=minv,maxval=maxv))
            
        for j in range(0,nbins):
            dm_data.ydata[j].nbins = pdf_bins
            dm_data.ydata[j].nfilter = nfilter_real
            dm_data.ydata[j].nvars = len(profs[0].yvars['dm'])
            dm_data.ydata[j].nwvars = len(profs[0].weightvars['dm'])
            stats_utils.allocate_pdf(dm_data.ydata[j])
            for k in range(0, len(profs[0].yvars['dm'])):
                tempstr = 'dm/'+profs[0].yvars['dm'][k]
                dm_data.ydata[j].varname.T.view('S128')[k] = tempstr.ljust(128)
                dm_data.ydata[j].scaletype.T.view('S128')[k] = mybins[k][1].ljust(128)
                dm_data.ydata[j].bins[:,k] = mybins[k][0]
                dm_data.ydata[j].do_binning[k] = do_binning
                dm_data.ydata[j].zero_index[k] = mybins[k][2]
                dm_data.ydata[j].linthresh[k] = mybins[k][3]
            for k in range(0, len(profs[0].weightvars['dm'])):
                tempstr = 'dm/'+profs[0].weightvars['dm'][k]
                dm_data.ydata[j].wvarnames.T.view('S128')[k] = tempstr.ljust(128)

        # And now, compute dm data profiles!
        if dm_data.nyvar > 0 and dm_data.nwvar > 0:
            partprofmod.onedprofile(group.obj.simulation.fullpath,selected_reg,dm_data,lmax)
    else:
        dm_data = None
    
    # Organise everything in the Profile object
    xdata = np.zeros((3,profs[0].nbins+1))
    if hydro_data != None:
        xdata[0,:] = hydro_data.xdata
    if star_data != None:
        xdata[1,:] = star_data.xdata
    if dm_data != None:
        xdata[2,:] = dm_data.xdata
    for i in range(0,nfilter):
        profs[i].xdata = group.obj.array(xdata, get_code_units(profs[i].xvar))
    # Save hydro y data
    if hydro_data != None:
        counter = 0
        for i in range(0,nfilter):
            if profs_fr[i] == True:
                prof = profs[i]
                prof.ydata['hydro'] = []
                for v,var in enumerate(prof.yvars['hydro']):
                    mydata = np.zeros((nbins,len(prof.weightvars['hydro']),7))
                    for j in range(0,nbins):
                        if hydro_data.ydata[j].do_binning[v]:
                            mydata[j,:,:] = pdf_handler_to_stats(group.obj,hydro_data.ydata[j],v,counter,verbose=True)
                        else:
                            mydata[j,:,:] = np.full((len(prof.weightvars['hydro']),7),hydro_data.ydata[j].total[v,counter,:,0])
                    prof.ydata['hydro'].append(group.obj.array(mydata, get_code_units(prof.yvars['hydro'][v])))
                counter += 1
    # Save star y data
    if star_data != None:
        counter = 0
        for i in range(0,nfilter):
            if profs_fr[i] == True:
                prof = profs[i]
                prof.ydata['star'] = []
                for v,var in enumerate(prof.yvars['star']):
                    mydata = np.zeros((nbins,len(prof.weightvars['star']),7))
                    for j in range(0,nbins):
                        if star_data.ydata[j].do_binning[v]:
                            mydata[j,:,:] = pdf_handler_to_stats(group.obj,star_data.ydata[j],v,counter,verbose=True)
                        else:
                            mydata[j,:,:] = np.full((len(prof.weightvars['star']),7),star_data.ydata[j].total[v,counter,:,0])
                    prof.ydata['star'].append(group.obj.array(mydata, get_code_units(prof.yvars['star'][v])))
                counter += 1
    # Save dm y data
    if dm_data != None:
        counter = 0
        for i in range(0,nfilter):
            if profs_fr[i] == True:
                prof = profs[i]
                prof.ydata['dm'] = []
                for v,var in enumerate(prof.yvars['dm']):
                    mydata = np.zeros((nbins,len(prof.weightvars['dm']),7))
                    for j in range(0,nbins):
                        if dm_data.ydata[j].do_binning[v]:
                            mydata[j,:,:] = pdf_handler_to_stats(group.obj,dm_data.ydata[j],v,counter,verbose=True)
                        else:
                            mydata[j,:,:] = np.full((len(prof.weightvars['dm']),7),dm_data.ydata[j].total[v,counter,:,0])
                    prof.ydata['dm'].append(group.obj.array(mydata, get_code_units(prof.yvars['dm'][v])))
                counter += 1
    if save:
        profs_to_save = [profs[i] for i in range(0,nfilter) if profs_fr[i] == True]
        write_profiles(group.obj, nfilter_real, ozy_file, hydro_data, star_data, dm_data, profs_to_save)
    
    for index, porig in enumerate(profs_fr):
        if porig != True:
            profs[index] = profs_fr[index]
    if nfilter > 1:
        return profs
    else:
        return profs[0]

def check_if_same_profile(hd, profile):
    """This function checks if a profile for an object already exists with the same attributes."""
    if profile.rm_subs:
        prof_key = '_data/profiles_nosubs/'
    else:
        prof_key = '_data/profiles/'
    if not str(profile.group.type)+prof_key+str(profile.group._index) in hd:
        return False, 'none'
    for p in hd[str(profile.group.type)+prof_key+str(profile.group._index)].keys():
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
    return name
def write_profiles(obj, nfilter, ozy_file, hydro, star, dm, profs):
    """This function writes the resulting profile data for this group to the original OZY HDF5 file."""

    f = h5py.File(ozy_file, 'r+')
    if profs[0].rm_subs:
        prof_key = '_data/profiles_nosubs/'
    else:
        prof_key = '_data/profiles/'
        
    # Create group in HDF5 file
    try:
        profiles = f.create_group(str(profs[0].group.type)+prof_key+str(profs[0].group._index))
    except:
        profiles = f[str(profs[0].group.type)+prof_key+str(profs[0].group._index)]
    for i in range(0,nfilter):
        prof = profs[i]
        # Clean data and save to dataset
        prof_name = get_profile_name(profiles, prof)
        hdprof = profiles.create_group(prof_name)
        prof._serialise(hdprof)
        # Save x data
        xdata = np.zeros((3,prof.nbins+1))
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
            for v,var in enumerate(prof.yvars['hydro']):
                mydata = np.zeros((prof.nbins,len(prof.weightvars['hydro']),7))
                for j in range(0,prof.nbins):
                    mydata[j,:,:] = pdf_handler_to_stats(obj,hydro.ydata[j],v,i)
                    if hydro.ydata[j].do_binning[v]:
                        mydata[j,:,:] = pdf_handler_to_stats(obj,hydro.ydata[j],v,i)
                    else:
                        mydata[j,:,:] = np.full((len(prof.weightvars['hydro']),7),hydro.ydata[j].total[v,i,:,0])
                clean_hydro.create_dataset(var, data=mydata)
                clean_hydro[var].attrs.create('units', get_code_units(prof.yvars['hydro'][v]))
                clean_hydro[var].attrs.create('weightvars', prof.weightvars['hydro'][:])
        # Save star y data
        if star != None:
            clean_star = hdprof.create_group('star')
            for v,var in enumerate(prof.yvars['star']):
                mydata = np.zeros((prof.nbins,len(prof.weightvars['star']),7))
                for j in range(0,prof.nbins):
                    mydata[j,:,:] = pdf_handler_to_stats(obj,star.ydata[j],v,i)
                    if star.ydata[j].do_binning[v]:
                        mydata[j,:,:] = pdf_handler_to_stats(obj,star.ydata[j],v,i)
                    else:
                        mydata[j,:,:] = np.full((len(prof.weightvars['star']),7),star.ydata[j].total[v,i,:,0])
                clean_star.create_dataset(var, data=mydata)
                clean_star[var].attrs.create('units', get_code_units(prof.yvars['star'][v]))
                clean_star[var].attrs.create('weightvars', prof.weightvars['star'][:])
        # Save dm y data
        if dm != None:
            clean_dm = hdprof.create_group('dm')
            for v,var in enumerate(prof.yvars['dm']):
                mydata = np.zeros((prof.nbins,len(prof.weightvars['dm']),7))
                for j in range(0,prof.nbins):
                    mydata[j,:,:] = pdf_handler_to_stats(obj,dm.ydata[j],v,i)
                    if dm.ydata[j].do_binning[v]:
                        mydata[j,:,:] = pdf_handler_to_stats(obj,dm.ydata[j],v,i)
                    else:
                        mydata[j,:,:] = np.full((len(prof.weightvars['dm']),7),dm.ydata[j].total[v,i,:,0])
                clean_dm.create_dataset(var, data=mydata)
                clean_dm[var].attrs.create('units', get_code_units(prof.yvars['dm'][v]))
                clean_dm[var].attrs.create('weightvars', prof.weightvars['dm'][:])
    f.close()
    return

def find_bin_pos(value, bins, zero_index, zero_eps=0.1, scaletype='symlog'):
    nbins = len(bins) - 1

    # Handle different scale types
    if scaletype == 'log_even':
        if value <= 0:
            raise ValueError("log_even scale type expects positive values.")
        ibin = int(nbins * (np.log10(value) - np.log10(bins[0])) / (np.log10(bins[-1]) - np.log10(bins[0])))
    elif scaletype == 'linear_even':
        ibin = int(nbins * (value - bins[0]) / (bins[-1] - bins[0]))
    elif scaletype == 'symlog':
        if value < -zero_eps:
            # Negative logarithmic region
            value = np.log10(-value)
            ibin = -int((zero_index-1)*(value - np.log10(-bins[0])) / (np.log10(-bins[0]) - np.log10(-bins[zero_index - 1]))) + 1
        elif value > zero_eps:
            # Positive logarithmic region
            value = np.log10(value)
            ibin = int((nbins - zero_index + 1) * (value - np.log10(bins[zero_index])) / (np.log10(bins[-1]) - np.log10(bins[zero_index]))) + zero_index
        else:
            # Linear region around zero
            ibin = zero_index  # Assuming zero_epsilon is located at zero_index

    return ibin