import numpy as np
import h5py
import os
import ozy
from unyt import unyt_array,unyt_quantity
from .utils import init_region,init_filter_hydro,\
                    get_code_units, get_plotting_def,\
                    check_need_gravity,check_need_rt,\
                    check_need_neighbours,pdf_handler_to_stats
from .variables_settings import geometrical_variables,raw_gas_variables,\
                                derived_gas_variables,gravity_variables,\
                                raw_part_variables,derived_part_variables,\
                                star_variables

from .amr.amr2_pkg import amr_profiles as amrprofmod
from .amr.amr2_pkg import stats_utils
from .amr.amr2_pkg import geometrical_regions as geo
from .amr.amr2_pkg import io_ramses
from .part.part2_pkg import part_profiles as partprofmod
from .part.part2_pkg import filtering_part
from .amr.amr2_pkg import filtering_hydro

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
        self.yvars = []
        self.weightvars = []
        self.xdata = None
        self.ydata = []
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
        self.region['zmin'] = self.obj.quantity(reg.zmin, 'code_length')
        self.region['zmax'] = self.obj.quantity(reg.zmax, 'code_length')
    
    def _get_python_filter(self,filt):
        """Save the Fortran derived type as a dictionary inside the PhaseDiagram class (only the necessary info)."""
        if isinstance(filt,filtering_hydro.filter_hydro):
            self.filter = dict()
            self.filter['name'] = filt.name.decode().split(' ')[0]
            self.filter['conditions'] = []
            if self.filter['name'] != 'none':
                for i in range(0, filt.ncond):
                    cond_var = filt.cond_vars_name.T.view('S128')[i][0].decode().split(' ')[0]
                    cond_op = filt.cond_ops.T.view('S2')[i][0].decode().split(' ')[0]
                    cond_units = get_code_units(cond_var,'gas')
                    cond_value = self.obj.quantity(filt.cond_vals[i], str(cond_units))
                    cond_str = cond_var+'/'+cond_op+'/'+str(cond_value.d)+'/'+cond_units
                    self.filter['conditions'].append(cond_str)
        elif isinstance(filt,filtering_part.filter_part):
            self.filter = dict()
            self.filter['name'] = filt.name.decode().split(' ')[0]
            self.filter['conditions'] = []
            if self.filter['name'] != 'none':
                for i in range(0, filt.ncond):
                    cond_var = filt.cond_vars_name.T.view('S128')[i][0].decode().split(' ')[0]
                    cond_op = filt.cond_ops.T.view('S2')[i][0].decode().split(' ')[0]
                    cond_units = get_code_units(cond_var,'part')
                    if get_part_vartype(cond_var) == 1:
                        value = filt.cond_vals_d[i]
                    elif get_part_vartype(cond_var) == 2:
                        value = filt.cond_vals_i[i]
                    elif get_part_vartype(cond_var) == 3:
                        value = filt.cond_vals_b[i]
                    cond_value = self.obj.quantity(value, str(cond_units))
                    cond_str = cond_var+'/'+cond_op+'/'+str(cond_value.d)+'/'+cond_units
                    self.filter['conditions'].append(cond_str)



def compute_profile_hydro(group,ozy_file,xvar,yvars,weightvars,minval,maxval,linthresh=None,lmax=0,nbins=100,
                    region_type='sphere',filter_conds=['none'],
                    filter_name=['none'],recompute=False,save=False,scaletype='log_even',
                    pdf_bins=100,regime_type='',
                    rmin=(0.0,'rvir'), rmax=(0.2,'rvir'), zmin=(0.0,'rvir'), zmax=(0.2,'rvir'),
                    mycentre=([0.5,0.5,0.5],'rvir'), myaxis=np.array([1.,0.,0.]),
                    remove_subs=False,cr_st=False,cr_heat=False,Dcr=0.0,
                    verbose=False,do_binning=True):
    """Function which computes a 1D profile for a given group object."""
    from .utils import structure_regions,get_code_bins
    # 1. Determine if we are handling a snapshot or a catalogue OZY object
    if isinstance(group,ozy.Snapshot):
        from ozy.group import Group
        obj = group
        group = Group(obj)
        use_snapshot = True
    else:
        obj = group.obj
        use_snapshot = False

    # 2. Activate the Fortran90 tools verbose options
    if verbose:
        io_ramses.activate_verbose()

    # 3. Check for initial compliance
    if not isinstance(xvar, str):
        if verbose: print('Single x variable 1D profile supported!')
        exit

    # 4. Setup the Profile object for each requested filter
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
        prof.yvars = []
        prof.weightvars = []
        profs.append(prof)
    
    # 5. Loop over the yvars and weightvars to fill the Profile objects
    use_gravity = False
    use_rt = False
    use_neigh = False
    for var in yvars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        
        if var_type == 'gas':
            use_gravity = check_need_gravity(var_name,var_type) or use_gravity
            use_rt = check_need_rt(var_name,var_type) or use_rt
            use_neigh = check_need_neighbours(var_name,var_type) or use_neigh
            if var_name in geometrical_variables or var_name in raw_gas_variables \
                or var_name in derived_gas_variables or var_name in gravity_variables:
                for i in range(0, nfilter):
                    profs[i].yvars.append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!',var_name)
        elif var_type == 'part':
            raise KeyError('Particle variables are not supported in hydro profiles. Please check!',var_name)
    for var in weightvars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            if var_name != 'cumulative' and var_name != 'count':
                use_gravity = check_need_gravity(var_name,var_type) or use_gravity
                use_rt = check_need_rt(var_name,var_type) or use_rt
                use_neigh = check_need_neighbours(var_name,var_type) or use_neigh
            if var_name in geometrical_variables or var_name in raw_gas_variables \
                or var_name in derived_gas_variables or var_name in gravity_variables:
                for i in range(0, nfilter):
                    profs[i].weightvars.append(var_name)
            elif var_name == 'cumulative' or var_name == 'count':
                for i in range(0, nfilter):
                    profs[i].weightvars.append(var_name)
            else:
                raise KeyError('This gas weight variable is not supported. Please check!',var_name)
        elif var_type == 'part':
            raise KeyError('Particle weight variables are not supported in hydro profiles. Please check!',var_name)
    
    # 6. Check if the xvar needs gravity,rt,neighbours
    use_gravity = check_need_gravity(xvar,'gas') or use_gravity
    use_neigh = check_need_neighbours(xvar,'gas') or use_neigh
    use_rt = check_need_rt(xvar,'gas') or use_rt
    
    # 7. Check that the xaxis min and max quantities have the units expected for that variable
    try:
        minval = minval.to(get_code_units(xvar,'gas'))
        maxval = maxval.to(get_code_units(xvar,'gas'))
    except:
        raise ValueError(f"It seems the dimensions of your bins min \
                            ({minval.units}) and max ({minval.units}) \
                            values do not agree with the dimensions of \
                            the chosen xvar ({xvar},{get_code_units(xvar,'gas')})")
    
    # 8. Now create region
    if use_snapshot:
        # This is for the case of not including an OZY catalogue
        group.position = obj.array(mycentre[0],mycentre[1])
        group.angular_mom['total'] = np.array([0.,0.,1.])
        group.velocity = obj.array([0.,0.,0.],'code_velocity')
    
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
        elif not np.array_equal(mycentre,group.position) and np.array_equal(myaxis,group.angular_mom['total']):
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

    # 9. Now create filters and regions for each profile
    filts = []
    for i in range(0,nfilter):
        if isinstance(filter_conds[i],list):
            cond_var = filter_conds[i][0].split('/')[0]
        else:
            cond_var = filter_conds[i].split('/')[0]
        if cond_var in geometrical_variables or cond_var in raw_gas_variables \
            or cond_var in derived_gas_variables or cond_var in gravity_variables:
            f = init_filter_hydro(filter_conds[i],filter_name[i],obj)
            print('Initialized hydro filter %s with conditions %s'%(filter_name[i],filter_conds[i]))
        else:
            # When a filter asks for a variable not existent in the common_variables
            # or the grid_variables dictionaries just ignore it and set it to blank
            f = init_filter_hydro('none','none',obj)
        filts.append(f)
        # Save region details to profile object
        prof = profs[i]
        prof._get_python_region(selected_reg)
        # And save to profile object
        prof._get_python_filter(f)

    # 10. Check if profile data is already present and if it coincides with the new one
    if not ozy_file is None:
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
    else:
        nfilter_real = nfilter
        profs_fr = [True]*nfilter
    
    # 11. If substructre is removed, obtain regions
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
    
    # 12. Initialise hydro profile data object
    hydro_data = amrprofmod.profile_handler()
    hydro_data.profdim = 1
    hydro_data.xvarname = xvar
    hydro_data.nfilter = nfilter_real
    hydro_data.nyvar = len(profs[0].yvars)
    hydro_data.nwvar = len(profs[0].weightvars)
    hydro_data.nbins = nbins
    hydro_data.nsubs = nsubs
    hydro_data.cr_st = cr_st
    hydro_data.cr_heat = cr_heat
    hydro_data.Dcr = Dcr

    amrprofmod.allocate_profile_handler(hydro_data)
    for i in range(0, len(profs[0].yvars)):
        hydro_data.yvarnames.T.view('S128')[i] = profs[0].yvars[i].ljust(128)
    for i in range(0, len(profs[0].weightvars)):
        hydro_data.wvarnames.T.view('S128')[i] = profs[0].weightvars[i].ljust(128)
    
    # 13. Add the scaletype for the xaxis and the pre-computed bin edges
    bin_edges, stype, zero_index, lint = get_code_bins(group.obj,'gas',xvar,nbins=nbins,logscale=scaletype,
                                                minval=minval,maxval=maxval,linthresh=linthresh)
    hydro_data.xdata = bin_edges
    hydro_data.scaletype = stype.ljust(128)
    hydro_data.linthresh = lint
    hydro_data.zero_index = zero_index
    
    if remove_subs and nsubs>0:
        for i in range(0,nsubs):
            hydro_data.subs[i] = subs[i]
    hydro_data.use_gravity = use_gravity
    hydro_data.use_rt = use_rt
    hydro_data.use_neigh = use_neigh
    
    counter = 0
    for i in range(0,nfilter):
        if profs_fr[i] == True:
            hydro_data.filters[counter] = filts[i]
            counter += 1

    # 14. Now add the PDF handler for each quantity and each xbin
    mybins = []
    for i in range(0, len(profs[0].yvars)):
        plot_def = get_plotting_def(profs[0].yvars[i],'gas')
        minv = group.obj.quantity(plot_def['bin_min'+regime_type],plot_def['units'])
        maxv = group.obj.quantity(plot_def['bin_max'+regime_type],plot_def['units'])
        mybins.append(get_code_bins(group.obj,'gas',profs[0].yvars[i],pdf_bins,
                                    minval=minv,maxval=maxv))
        
    for j in range(0,nbins):
        hydro_data.ydata[j].nbins = pdf_bins
        hydro_data.ydata[j].nfilter = nfilter_real
        hydro_data.ydata[j].nvars = len(profs[0].yvars)
        hydro_data.ydata[j].nwvars = len(profs[0].weightvars)
        stats_utils.allocate_pdf(hydro_data.ydata[j])
        for k in range(0, len(profs[0].yvars)):
            hydro_data.ydata[j].varname.T.view('S128')[k] = profs[0].yvars[k].ljust(128)
            hydro_data.ydata[j].scaletype.T.view('S128')[k] = mybins[k][1].ljust(128)
            hydro_data.ydata[j].bins[:,k] = mybins[k][0]
            hydro_data.ydata[j].do_binning[k] = do_binning
            hydro_data.ydata[j].zero_index[k] = mybins[k][2]
            hydro_data.ydata[j].linthresh[k] = mybins[k][3]
        for k in range(0, len(profs[0].weightvars)):
            hydro_data.ydata[j].wvarnames.T.view('S128')[k] = profs[0].weightvars[k].ljust(128)
        
    # 15. And now, compute hydro data profiles!
    if hydro_data.nyvar > 0 and hydro_data.nwvar > 0:
        if obj.use_vardict:
            amrprofmod.onedprofile(group.obj.simulation.fullpath,selected_reg,hydro_data,lmax,obj.vardict)
        else:
            amrprofmod.onedprofile(group.obj.simulation.fullpath,selected_reg,hydro_data,lmax)
    
    # 16. Organise everything in the Profile object
    counter = 0
    for i in range(0,nfilter):
        profs[i].xdata = group.obj.array(hydro_data.xdata, get_code_units(profs[i].xvar,'gas'))
        if profs_fr[i] == True:
            prof = profs[i]
            prof.ydata = []
            for v,var in enumerate(prof.yvars):
                mydata = np.zeros((nbins,len(prof.weightvars),7))
                for j in range(0,nbins):
                    if hydro_data.ydata[j].do_binning[v]:
                        mydata[j,:,:] = pdf_handler_to_stats(group.obj,'gas',hydro_data.ydata[j],v,counter,verbose=True)
                    else:
                        mydata[j,:,:] = np.full((len(prof.weightvars),7),hydro_data.ydata[j].total[v,counter,:,0])
                prof.ydata.append(group.obj.array(mydata, get_code_units(prof.yvars[v],'gas')))
            counter += 1
    if save:
        profs_to_save = [profs[i] for i in range(0,nfilter) if profs_fr[i] == True]
        write_profiles(group.obj, nfilter_real, ozy_file, hydro_data, profs_to_save)
    
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