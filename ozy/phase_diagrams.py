import numpy as np
import h5py
import os
import ozy
from unyt import unyt_array,unyt_quantity
from ozy.utils import init_region,init_filter,gent_curve_T,\
                    check_need_neighbours, get_code_units, \
                    get_plotting_def
from variables_settings import geometrical_variables,raw_gas_variables,\
    raw_star_variables,raw_dm_variables,derived_gas_variables,\
    derived_star_variables,derived_dm_variables,gravity_variables,\
    basic_conv,circle_dictionary
# TODO: Allow for parallel computation of phase diagrams.
from joblib import Parallel, delayed
from amr2 import vectors
from amr2 import geometrical_regions as geo
from amr2 import amr_profiles as amrprofmod

blacklist = [
    'zvars','weightvars','xdata','ydata','zdata'
]

class PhaseDiagram(object):

    def __init__(self,group):
        self.group = group
        self.obj = group.obj
        self.nbins = [0]
        self.xvar = None
        self.yvar = None
        self.region = None
        self.filter = None
        self.rm_subs = False
        self.lmax = 0
        self.zvars = {}
        self.weightvars = {}
        self.xdata = []
        self.ydata = []
        self.zdata = {}

    def _serialise(self, hdd):
        """This makes possible to save the group phase diagram attrs as dataset attributes of an HDF5 group."""
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
        """Save the Fortran derived type as a dictionary inside the PhaseDiagram class (only the necessary info)."""
        self.region = {}
        self.region['type'] = reg.name.decode().split(' ')[0]
        self.region['centre'] = self.obj.array([reg.centre.x, reg.centre.y, reg.centre.z], 'code_length')
        self.region['axis'] = self.obj.array([reg.axis.x, reg.axis.y, reg.axis.z], 'dimensionless')
    
    def _get_python_filter(self,filt):
        """Save the Fortran derived type as a dictionary inside the PhaseDiagram class (only the necessary info)."""
        self.filter = {}
        self.filter['name'] = filt.name.decode().split(' ')[0]
        self.filter['conditions'] = []
        if self.filter['name'] != 'none':
            for i in range(0, filt.ncond):
                cond_var = filt.cond_vars.T.view('S128')[i][0].decode().split(' ')[0]
                cond_op = filt.cond_ops.T.view('S2')[i][0].decode().split(' ')[0]
                cond_units = get_code_units(cond_var)
                cond_value = self.obj.quantity(filt.cond_vals[i], str(cond_units))
                cond_str = cond_var+'/'+cond_op+'/'+str(cond_value.d)+'/'+cond_units
                self.filter['conditions'].append(cond_str)

def compute_phase_diagram(group,ozy_file,xvar,yvar,zvars,weightvars,minval,maxval,linthresh=None,
                            lmax=0,nbins=[100,100],region_type='sphere',
                            filter_conds=['none'],filter_name=['none'],scaletype=['log_even','log_even'],
                            recompute=False,save=False,
                            rmin=(0.0,'rvir'), rmax=(0.2,'rvir'), zmin=(0.0,'rvir'), zmax=(0.2,'rvir'),
                            mycentre=([0.5,0.5,0.5],'code_length'), myaxis=np.array([1.,0.,0.]),
                            remove_subs=False,cr_st=False,cr_heat=False,Dcr=0.0,verbose=False):
    """Function which computes a phase diagram (2D profile) for a given group object."""
    from ozy.utils import structure_regions,get_code_bins

    if isinstance(group,ozy.Snapshot):
        obj = group
        group = ozy.group.Group(obj)
        use_snapshot = True
    else:
        obj = group.obj
        use_snapshot = False
    
    if not isinstance(xvar,str) or not isinstance(yvar,str):
        print('Only single x and y variable supported!')
        exit
    if len(nbins) != 2 or not isinstance(nbins,list):
        print('Number of bins should be given for both x and y axis. No more, no less!')
        exit
    remove_all = False
    if remove_subs =='all':
        remove_all = True
        remove_subs = True
    nfilter = len(filter_name)
    pds = []
    for i in range(0,nfilter):
        pd = PhaseDiagram(group)
        pd.rm_subs = remove_subs
        pd.nbins = nbins
        pd.xvar = xvar
        pd.yvar = yvar
        pd.zvars = dict(hydro = [],for_star = [], for_dm = [])
        pd.weightvars = dict(hydro = [],for_star = [], for_dm = [])
        pds.append(pd)
    # Begin by checking the combination of variables is correct
    # Variables should be given in the form "type"/"variable":
    # e.g. gas/density or star/ang_momentum_x
    for var in zvars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            if var_name in common_variables or var_name in grid_variables:
                for i in range(0, nfilter):
                    pds[i].zvars['hydro'].append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!: %s', var)
        elif var_type == 'star':
            if var_name in common_variables or var_name in particle_variables:
                for i in range(0, nfilter):
                    pds[i].zvars['for_star'].append(var_name)
            else:
                raise KeyError('This star variable is not supported. Please check!')
        elif var_type == 'dm':
            if var_name in common_variables or var_name in particle_variables:
                for i in range(0, nfilter):
                    pds[i].zvars['for_dm'].append(var_name)
            else:
                raise KeyError('This DM variable is not supported. Please check!')
    for var in weightvars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            if var_name in common_variables or var_name in grid_variables:
                for i in range(0, nfilter):
                    pds[i].weightvars['hydro'].append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!')
        elif var_type == 'star':
            if var_name in common_variables or var_name in particle_variables:
                for i in range(0, nfilter):
                    pds[i].weightvars['for_star'].append(var_name)
            else:
                raise KeyError('This star variable is not supported. Please check!')
        elif var_type == 'dm':
            if var_name in common_variables or var_name in particle_variables:
                for i in range(0, nfilter):
                    pds[i].weightvars['for_dm'].append(var_name)
            else:
                raise KeyError('This DM variable is not supported. Please check!')

    # Check that we do not have any inconsistency...
    if xvar in grid_variables and len(pds[0].zvars['for_star'])>0 or xvar in grid_variables and len(pds[0].zvars['for_dm'])>0:
        raise KeyError("Having grid vs particle phase diagrams is not well-defined.")
    elif xvar in particle_variables and len(pds[0].zvars['hydro'])>0:
        raise KeyError("Having particle vs grid phase diagrams is not well-defined.")
    
    # Check that the xaxis min and max quantities have the units expected for that variable
    try:
        minval[0] = minval[0].to(get_code_units(xvar))
        maxval[0] = maxval[0].to(get_code_units(xvar))
    except:
        raise ValueError(f"It seems the dimensions of your bins min ({minval[0].units}) and max ({maxval[0].units}) values do not agree with the dimensions of the chosen xvar ({xvar},{get_code_units(xvar)})")
    # Check that the yaxis min and max quantities have the units expected for that variable
    try:
        minval[1] = minval[1].to(get_code_units(yvar))
        maxval[1] = maxval[1].to(get_code_units(yvar))
    except:
        raise ValueError(f"It seems the dimensions of your bins min ({minval[1].units}) and max ({maxval[1].units}) values do not agree with the dimensions of the chosen yvar ({yvar},{get_code_units(yvar)})")

    # Now create region
    if use_snapshot:
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
    # Now create filters and regions for each phase diagram
    filts = []
    for i in range(0,nfilter):
        filt = init_filter(filter_conds[i], filter_name[i], group)
        filts.append(filt)
        # Save region details to phase diagram object
        pd = pds[i]
        pd._get_python_region(selected_reg)
        # And save to phase diagram object
        pd._get_python_filter(filt)

    # Check if phase data is already present and if it coincides with the new one
    if not ozy_file is None:
        f = h5py.File(ozy_file, 'r+')
        pds_fr = []
        for i in range(0,nfilter):
            pd = pds[i]
            pd_present, pd_key = check_if_same_phasediag(f,pd)
            if pd_present:
                group._init_phase_diagrams()
                if remove_subs:
                    for i,p in enumerate(group.phase_diagrams_nosubs):
                        if p.key == pd_key:
                            selected_pd = i
                            break
                    pd_temp = group.phase_diagrams_nosubs[selected_pd]
                    # TODO: This should be done for all type of variables, not just hydro
                    for pd_field in pd.zvars['hydro']:
                        if not pd_field in pd_temp.zvars['hydro']:
                            if verbose: print('Recomputing phase-diagram because of missing requested variables!')
                            recompute = True
                            pds_fr.append(True)
                            break
                else:
                    for i,p in enumerate(group.phase_diagrams):
                        if p.key == pd_key:
                            selected_pd = i
                            break
                    pd_temp = group.phase_diagrams[selected_pd]
                    # TODO: This should be done for all type of variables, not just hydro
                    for pd_field in pd.zvars['hydro']:
                        if not pd_field in pd_temp.zvars['hydro']:
                            if verbose: print('Recomputing phase-diagram because of missing requested variables!')
                            recompute = True
                            pds_fr.append(True)
                            break
            if pd_present and recompute:
                if remove_subs:
                    del f[str(pd.group.type)+'_data/phase_diagrams_nosubs/'+str(group._index)+'/'+str(pd_key)]
                else:
                    del f[str(pd.group.type)+'_data/phase_diagrams/'+str(group._index)+'/'+str(pd_key)]
                pds_fr.append(True)
                if verbose: print('Overwriting phase diagram data in %s_data'%group.type)
            elif pd_present and not recompute:
                if verbose: print('Phase diagram data with same details already present for galaxy %s. No overwritting!'%group._index)
                group._init_phase_diagrams()
                if remove_subs:
                    for i,p in enumerate(group.phase_diagrams_nosubs):
                        if p.key == pd_key:
                            selected_pd = i
                            break
                    pds_fr.append(group.phase_diagrams_nosubs[selected_pd])
                else:
                    for i,p in enumerate(group.phase_diagrams):
                        if p.key == pd_key:
                            selected_pd = i
                            break
                    pds_fr.append(group.phase_diagrams[selected_pd])
            else:
                pds_fr.append(True)
                if verbose: print('Writing phase diagram data in %s_data'%group.type)
        f.close()
    
    nfilter_real = pds_fr.count(True)
    if nfilter_real == 0:
        if nfilter > 1:
            return pds_fr
        else:
            return pds_fr[0]
        
    # If substructre is removed, obtain regions
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

    # Initialise hydro phase diagram data object
    hydro_data = amrprofmod.profile_handler_twod()
    hydro_data.profdim = 2
    hydro_data.xvarname = xvar
    hydro_data.yvarname = yvar
    hydro_data.nfilter = nfilter_real
    hydro_data.nzvar = len(pd.zvars['hydro'])
    hydro_data.nwvar = len(pd.weightvars['hydro'])
    hydro_data.nbins = np.asarray(nbins,order='F')
    hydro_data.nsubs = nsubs
    hydro_data.cr_st = cr_st
    hydro_data.cr_heat = cr_heat
    hydro_data.Dcr = Dcr
    amrprofmod.allocate_profile_handler_twod(hydro_data)
    for i in range(0, len(pd.zvars['hydro'])):
        hydro_data.zvarnames.T.view('S128')[i] = pd.zvars['hydro'][i].ljust(128)
    for i in range(0, len(pd.weightvars['hydro'])):
        hydro_data.wvarnames.T.view('S128')[i] = pd.weightvars['hydro'][i].ljust(128)
    
    if remove_subs and nsubs>0:
        for i in range(0,nsubs):
            hydro_data.subs[i] = subs[i]
    
    # Add the scaletype for the xaxis and the preo-computed bin edges
    bin_edges, stype, zero_index, linthresh = get_code_bins(group.obj,'gas/'+xvar,nbins=nbins[0],logscale=scaletype[0],
                                        minval=minval[0],maxval=maxval[0],linthresh=linthresh)
    hydro_data.xdata = bin_edges
    hydro_data.scaletype.T.view('S128')[0] = stype.ljust(128)
    hydro_data.linthresh[0] = linthresh
    hydro_data.zero_index[0] = zero_index
    # Add the scaletype for the yaxis and the preo-computed bin edges
    bin_edges, stype, zero_index, linthresh = get_code_bins(group.obj,'gas/'+yvar,nbins=nbins[1],logscale=scaletype[1],
                                        minval=minval[1],maxval=maxval[1],linthresh=linthresh)
    hydro_data.ydata = bin_edges
    hydro_data.scaletype.T.view('S128')[1] = stype.ljust(128)
    hydro_data.linthresh[1] = linthresh
    hydro_data.zero_index[1] = zero_index
    
    counter = 0
    for i in range(0,nfilter):
        if pds_fr[i] == True:
            hydro_data.filters[counter] = filts[i]
            counter += 1
            
    # And now, compute hydro data phase diagrams!
    if hydro_data.nzvar > 0 and hydro_data.nwvar > 0:
        if obj.use_vardict:
            amrprofmod.twodprofile(obj.simulation.fullpath,selected_reg,filt,hydro_data,lmax,obj.vardict)
        else:
            amrprofmod.twodprofile(obj.simulation.fullpath,selected_reg,filt,hydro_data,lmax)
    
    xdata = 0.5*(hydro_data.xdata[:-1]+hydro_data.xdata[1:])
    ydata = 0.5*(hydro_data.ydata[:-1]+hydro_data.ydata[1:])
    
    counter = 0
    for i in range(0,nfilter):
        if pds_fr[i] == True:
            pd  = pds[i]
            pd.xdata.append(group.obj.array(xdata, get_code_units(pd.xvar)))
            pd.ydata.append(group.obj.array(ydata, get_code_units(pd.yvar)))
            pd.zdata['hydro'] = []
            
            for j in range(0, len(pd.zvars['hydro'])):
                code_units = get_code_units(pd.zvars['hydro'][j])
                copy_data = np.copy(hydro_data.zdata[counter,:,:,j,:,::2])
                pd.zdata['hydro'].append(group.obj.array(copy_data, code_units))
            counter += 1

    # TODO: Add phase diagram for particles
    star_data = None
    dm_data = None
    if save:
        pds_to_save = [pds[i] for i in range(0,nfilter) if pds_fr[i] == True]
        write_phasediag(group.obj, nfilter_real, ozy_file, hydro_data, star_data, dm_data, pds_to_save)
    
    for index, porig in enumerate(pds_fr):
        if porig != True:
            pds[index] = pds_fr[index]
    if nfilter > 1:
        return pds
    else:
        return pds[0]

def get_phasediag_name(phasediag_group,pd):
    """Create an individual phase diagram identifier name."""
    name = str(pd.xvar)+'_'+str(pd.yvar)
    name += '|'+str(pd.filter['name'])
    name += '|'+str(pd.region['type'])
    if name in phasediag_group:
        name += '|new'
    return name

def check_if_same_phasediag(hd,pd):
    """This function checks if a phase diagram for a group already exists with the same attributes."""
    if pd.rm_subs:
        prof_key = '_data/phase_diagrams_nosubs/'
    else:
        prof_key = '_data/phase_diagrams/'
    if not str(pd.group.type)+prof_key+str(pd.group._index) in hd:
        return False, 'none'
    for p in hd[str(pd.group.type)+prof_key+str(pd.group._index)].keys():
        check_xvar = (p.split('|')[0] == pd.xvar+'_'+pd.yvar)
        check_filtername = (p.split('|')[1] == pd.filter['name'])
        check_regiontype = (p.split('|')[2] == pd.region['type'])
        if check_xvar and check_filtername and check_regiontype:
            return True, p
    return False, 'none'

def write_phasediag(obj,nfilter,ozy_file,hydro,star,dm,pds):
    """This function writes the resulting phase diagram for this group to the original OZY HDF5 file."""

    f = h5py.File(ozy_file, 'r+')
    
    if pds[0].rm_subs:
        prof_key = '_data/phase_diagrams_nosubs/'
    else:
        prof_key = '_data/phase_diagrams/'

    # Create group in HDF5 file
    try:
        phase_diagrams = f.create_group(str(pds[0].group.type)+prof_key+str(pds[0].group._index))
    except:
        phase_diagrams = f[str(pds[0].group.type)+prof_key+str(pds[0].group._index)]
    
    for i in range(0,nfilter):
        pd = pds[i]
        # Clean data and save to dataset
        pd_name = get_phasediag_name(phase_diagrams,pd)
        hdpd = phase_diagrams.create_group(pd_name)
        pd._serialise(hdpd)

        # Save x and y data
        xdata = np.zeros((3,pd.nbins[0]))
        if hydro != None:
            xdata[0,:] = hydro.xdata[1:]
        if star != None:
            xdata[1,:] = star.xdata[1:]
        if dm != None:
            xdata[2,:] = dm.xdata[1:]
        hdpd.create_dataset('xdata', data=xdata)
        hdpd['xdata'].attrs.create('units', get_code_units(pd.xvar))

        ydata = np.zeros((3,pd.nbins[1]))
        if hydro != None:
            ydata[0,:] = hydro.ydata[1:]
        if star != None:
            ydata[1,:] = star.ydata[1:]
        if dm != None:
            ydata[2,:] = dm.ydata[1:]
        hdpd.create_dataset('ydata', data=ydata)
        hdpd['ydata'].attrs.create('units', get_code_units(pd.yvar))

        # Save hydro z data
        if hydro != None:
            clean_hydro = hdpd.create_group('hydro')
            for v,var in enumerate(pd.zvars['hydro']):
                clean_hydro.create_dataset(var, data=hydro.zdata[i,:,:,v,:,::2])
                clean_hydro[var].attrs.create('units', get_code_units(pd.zvars['hydro'][v]))
            clean_hydro.create_dataset('weightvars', data=pd.weightvars['hydro'][:])
        # TODO: Save particle data
    f.close()
    return

def plot_single_phase_diagram(pd,field,name,weightvar='cumulative',logscale=True,redshift=True,stats='none',extra_labels='none',gent=False,powell=False):
    """This function uses the information from PhaseDiagram following the OZY format and
        plots following the OZY standards."""
    
    # Make required imports
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from matplotlib.colors import LogNorm,SymLogNorm
    import seaborn as sns
    sns.set(style="white")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    hfont = {'fontname':'Helvetica'}
    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', **{'family' : "serif"})

    # Check that the required field is actually in the PhaseDiagram
    if field not in pd.zvars['hydro']:
        print('The field %s is not included in this PhaseDiagram object. Aborting!'%field)
        exit
    else:
        for f in range(0,len(pd.zvars['hydro'])):
            if pd.zvars['hydro'][f] == field:
                field_index = f
                break
        if weightvar not in pd.weightvars['hydro']:
            print('The weight field %s is not included in this PhaseDiagram object. Aborting!'%weightvar)
            exit
        else:
            for w in range(0, len(pd.weightvars['hydro'])):
                if pd.weightvars['hydro'][w] == weightvar:
                    weight_index = w
                    break

    # Since everything appears fine, we begin plotting...
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    plotting_x = get_plotting_def(pd.xvar)
    plotting_y = get_plotting_def(pd.yvar)
    plotting_z = get_plotting_def(field)
    ax.set_xlabel(plotting_x['label'],fontsize=18)
    ax.set_ylabel(plotting_y['label'],fontsize=18)

    ax.tick_params(labelsize=14,direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.minorticks_on()
    ax.tick_params(which='both',axis="both",direction="in")
    ax.set_xlim([1e-30,1e-20])
    ax.set_ylim([15,3e+8])
    ax.set_xscale('log')
    ax.set_yscale('log')
    code_units_x = get_code_units(pd.xvar)
    x = pd.obj.array(pd.xdata[0].d,code_units_x)
    x = x.in_units(plotting_x['units'])
    code_units_y = get_code_units(pd.yvar)
    y = pd.obj.array(pd.ydata[0].d,code_units_y)
    y = y.in_units(plotting_y['units'])
    code_units_z = get_code_units(field)
    z = np.array(pd.zdata['hydro'][field_index][:,:,weight_index].d,order='F')
    z = pd.obj.array(z,code_units_z)
    print(field,np.nanmin(z).in_units(plotting_z['units']),np.nanmax(z).in_units(plotting_z['units']))
    sim_z = pd.obj.simulation.redshift
    if logscale:
        if plotting_z['symlog']:
            plot = ax.pcolormesh(x,y,
                                z[:,:,0].in_units(plotting_z['units']).T.d,
                                shading='auto',
                                cmap=plotting_z['cmap'],
                                norm=SymLogNorm(linthresh=plotting_z['linthresh'],
                                linscale=plotting_z['linscale'],
                                vmin=plotting_z['vmin_galaxy'],
                                vmax=plotting_z['vmax_galaxy'],
                                base=10))
        else:
            plot = ax.pcolormesh(x,y,
                                z[:,:,0].in_units(plotting_z['units']).T.d,
                                shading='auto',
                                cmap=plotting_z['cmap'],
                                norm=LogNorm(vmin=plotting_z['vmin_galaxy'],
                                vmax=plotting_z['vmax_galaxy']))
    else:
        plot = ax.pcolormesh(x,y,
                            z[:,:,0].in_units(plotting_z['units']).T.d,
                            shading='auto',
                            cmap=plotting_z['cmap'],
                            vmin=plotting_z['vmin'],
                            vmax=plotting_z['vmax'])
                            
    axcb = fig.colorbar(plot, orientation="vertical", pad=0.0)
    axcb.set_label(plotting_z['label'],fontsize=18)
    if redshift:
        ax.text(0.80, 0.9, 'z = '+str(round(sim_z, 2)),
                    transform=ax.transAxes, fontsize=20,verticalalignment='top',
                    color='black')
    if stats == 'mean':
        y_mean = np.zeros(len(x))
        for i in range(0, len(x)):
            a = np.nansum(y * z[:,i].T)
            b = np.nansum(z[:,i].T)
            y_mean[i] = a/b
        ax.plot(x,y_mean,color='r',linewidth=2)

    if gent:
        XX,YY = np.meshgrid(x,y)
        ax.fill_between([1e-30,8e-20], [gent_curve_T('cold',1e-30),gent_curve_T('cold',8e-20)],
                            [1,1],color='b', zorder=-10,alpha=0.2)
        ax.fill_between([1e-30,8e-20], [gent_curve_T('cold',1e-30),gent_curve_T('cold',8e-20)],
                            [gent_curve_T('hot',1e-30),gent_curve_T('hot',8e-20)],color='orange', 
                            zorder=-10,alpha=0.2)
        ax.fill_between([1e-30,8e-20], [gent_curve_T('hot',1e-30),gent_curve_T('hot',8e-20)],
                        [1e8,1e8],color='r', zorder=-10,alpha=0.2)
        if field == 'mass':
            z_cold = np.sum(z[:,:,0].T[YY<gent_curve_T('cold',XX)])
            z_warm = np.sum(z[:,:,0].T[(YY>gent_curve_T('cold',XX)) & (YY<gent_curve_T('hot',XX))])
            z_hot = np.sum(z[:,:,0].T[YY>gent_curve_T('hot',XX)])
            z_tot = np.sum(z[:,:,0])
            print(z_cold,z_warm,z_hot)
            print('Distribution of masses in the Gent phases:')
            print('Cold: %.3f, %.3e'%(100*z_cold/z_tot,z_cold))
            print('Warm: %.3f, %.3e'%(100*z_warm/z_tot,z_warm))
            print('Hot: %.3f, %.3e'%(100*z_hot/z_tot, z_hot))
            
    fig.subplots_adjust(top=0.97,bottom=0.1,left=0.1,right=0.95)
    fig.savefig(name+'.png',format='png',dpi=300)

def plot_compare_phase_diagram(pds,field,name,weightvar='cumulative',
                                scaletype='log_even',redshift=True,
                                stats='none',extra_labels='none',
                                gent=False,powell=False,doflows=False,
                                do_sf=False,layout='compact',
                                returnfig=False):

    # Make required imports
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from matplotlib.colors import LogNorm,SymLogNorm
    import seaborn as sns
    sns.set(style="white")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    hfont = {'fontname':'Helvetica'}
    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', **{'family' : "serif"})

    # How many phase diagram datasets have been provided
    if isinstance(pds,list):
        npds = len(pds)
    else:
        raise TypeError('You need to provide all phase diagram datasets in a list!')
    
    # Check that the required field is actually in the PDs provided
    field_indexes = []
    weight_indexes = []
    for i,pd in enumerate(pds):
        if doflows or do_sf:
            if field[i] not in pd[0].zvars['hydro']:
                raise KeyError('The field %s is not included in this PhaseDiagram object. Aborting!'%field[i])
            else:
                for f in range(0,len(pd[0].zvars['hydro'])):
                    if pd[0].zvars['hydro'][f] == field[i]:
                        field_indexes.append(f)
                        break
                if weightvar not in pd[0].weightvars['hydro']:
                    raise KeyError('The weight field %s is not included in this PhaseDiagram object. Aborting!'%weightvar)
                    exit
                else:
                    for w in range(0, len(pd[0].weightvars['hydro'])):
                        if pd[0].weightvars['hydro'][w] == weightvar:
                            weight_indexes.append(w)
                            break
        else:
            if field[i] not in pd.zvars['hydro']:
                raise KeyError('The field %s is not included in this PhaseDiagram object. Aborting!'%field[i])
            else:
                for f in range(0,len(pd.zvars['hydro'])):
                    if pd.zvars['hydro'][f] == field[i]:
                        field_indexes.append(f)
                        break
                if weightvar not in pd.weightvars['hydro']:
                    raise KeyError('The weight field %s is not included in this PhaseDiagram object. Aborting!'%weightvar)
                    exit
                else:
                    for w in range(0, len(pd.weightvars['hydro'])):
                        if pd.weightvars['hydro'][w] == weightvar:
                            weight_indexes.append(w)
                            break
    if len(field_indexes) != npds or len(weight_indexes) != npds:
        print('You should check the fields and weights available, not all your phase diagrams have them!')
        exit
    
    # With everything fine, we begin plotting…
    if layout == 'compact':
        ncol = 2
        nrow = int(len(pds)/ncol)
    elif layout == 'extended':
        nrow = 1
        ncol = len(pds)
    else:
        print('The layout asked is not allowed. Please check!')
        exit
    print(layout)
    if layout == 'compact':
        figsize = plt.figaspect(float((5.0 * nrow) / (5.0 * ncol)))
        fig = plt.figure(figsize=2*figsize, facecolor='w',edgecolor='k')
    elif layout == 'extended':
        figsize = plt.figaspect(float((6.0 * nrow) / (5.0 * ncol)))
        fig = plt.figure(figsize=figsize, facecolor='w',edgecolor='k')
    plot_grid = fig.add_gridspec(nrow, ncol, wspace=0, hspace=0)#,  hspace=0,left=0,right=1, bottom=0, top=1)
    ax = []
    for i in range(0,nrow):
        ax.append([])
        for j in range(0,ncol):
            ax[i].append(fig.add_subplot(plot_grid[i,j]))
    ax = np.asarray(ax)
    for i in range(0, ax.shape[0]):
        for j in range(0, ax.shape[1]):
            ipd = i*ax.shape[1] + j
            if ipd >= npds:
                # Clear that extra panel
                ax[i,j].get_xaxis().set_visible(False)
                ax[i,j].get_yaxis().set_visible(False)
                break
            if doflows or do_sf:
                pd = pds[ipd][0]
            else:
                pd = pds[ipd]
            plotting_x = get_plotting_def(pd.xvar)
            plotting_y = get_plotting_def(pd.yvar)
            plotting_z = get_plotting_def(field)
            print(field,plotting_z)
            ax[i,j].set_xlabel(plotting_x['label'],fontsize=20)
            # Get rid of the y-axis labels for the panels in the middle
            # TODO: The ticks should not be only for density and temperature!
            if ipd%2 == 0 and layout == 'compact':
                ax[i,j].set_ylabel(plotting_y['label'],fontsize=18)
            elif layout == 'compact':
                ax[i,j].axes.yaxis.set_visible(False)
                ax[i,j].set_xticks([1e-28,1e-26,1e-24,1e-22,1e-20])
            elif ipd != 0 and layout == 'extended':
                ax[i,j].axes.yaxis.set_visible(False)
                ax[i,j].set_xticks([1e-28,1e-26,1e-24,1e-22,1e-20])
            else:
                ax[i,j].set_ylabel(plotting_y['label'],fontsize=18)

            ax[i,j].tick_params(labelsize=14,direction='in')
            ax[i,j].xaxis.set_ticks_position('both')
            ax[i,j].yaxis.set_ticks_position('both')
            ax[i,j].minorticks_on()
            ax[i,j].tick_params(which='major',axis="both",direction="in")
            ax[i,j].set_xlim([1e-30,8e-20])
            ax[i,j].set_ylim([15,1e+8])
            
            code_units_x = get_code_units(pd.xvar)
            if scaletype=='log_even':
                ax[i,j].set_xscale('log')
                ax[i,j].set_yscale('log')
                x = pd.obj.array(10**(pd.xdata[0].d),code_units_x)
            else:
                x = pd.obj.array(pd.xdata[0].d,code_units_x)
            x = x.in_units(plotting_x['units'])
            code_units_y = get_code_units(pd.yvar)
            
            if scaletype=='log_even':
                y = pd.obj.array(10**(pd.ydata[0].d),code_units_y)
            else:
                y = pd.obj.array(pd.ydata[0].d,code_units_y)
            y = y.in_units(plotting_y['units'])
            code_units_z = get_code_units(field[i])
            z = np.array(pd.zdata['hydro'][field_indexes[ipd]][:,:,weight_indexes[ipd],0].d,order='F')
            z = pd.obj.array(z,code_units_z)
            print(field[i],np.nanmin(z).to(plotting_z['units']),np.nanmax(z).to(plotting_z['units']),pd.obj.simulation.redshift)
            sim_z = pd.obj.simulation.redshift

            if gent:
                XX,YY = np.meshgrid(x,y)
                ax[i,j].fill_between([1e-30,8e-20], [gent_curve_T('cold',1e-30),gent_curve_T('cold',8e-20)],
                                    [1,1],color='b', zorder=1,alpha=0.2)
                ax[i,j].fill_between([1e-30,8e-20], [gent_curve_T('cold',1e-30),gent_curve_T('cold',8e-20)],
                                    [gent_curve_T('hot',1e-30),gent_curve_T('hot',8e-20)],color='orange', 
                                    zorder=1,alpha=0.2)
                ax[i,j].fill_between([1e-30,8e-20], [gent_curve_T('hot',1e-30),gent_curve_T('hot',8e-20)],
                                [1e8,1e8],color='r', zorder=1,alpha=0.2)
                if field == 'mass':
                    z_cold = np.sum(z.T[YY<gent_curve_T('cold',XX)])
                    z_warm = np.sum(z.T[(YY>gent_curve_T('cold',XX)) & (YY<gent_curve_T('hot',XX))])
                    z_hot = np.sum(z.T[YY>gent_curve_T('hot',XX)])
                    z_tot = np.sum(z)
                    print(z_cold,z_warm,z_hot)
                    print('Distribution of masses in the Gent phases in %s:'%extra_labels[ipd])
                    print('Cold: %.3f, %.3e'%(100*z_cold/z_tot,z_cold))
                    print('Warm: %.3f, %.3e'%(100*z_warm/z_tot,z_warm))
                    print('Hot: %.3f, %.3e'%(100*z_hot/z_tot, z_hot))

            if plotting_z['symlog']:
                plot = ax[i,j].pcolormesh(x,y,
                                    z.in_units(plotting_z['units']).T,
                                    shading='auto',
                                    cmap=plotting_z['cmap'],
                                    norm=SymLogNorm(linthresh=plotting_z['linthresh'],
                                    linscale=plotting_z['linscale'],
                                    vmin=plotting_z['vmin_galaxy'],
                                    vmax=plotting_z['vmax_galaxy'],
                                    base=10))
            else:
                z[z<0.0] = np.min(z[z>0.0])
                z = np.nan_to_num(z, np.min(z))
                plot = ax[i,j].pcolormesh(x,y,
                                    z.in_units(plotting_z['units']).T,
                                    shading='auto',
                                    cmap=plotting_z['cmap'],
                                    norm=LogNorm(vmin=plotting_z['vmin_galaxy'],
                                    vmax=plotting_z['vmax_galaxy']))
            if redshift:
                ax[i,j].text(0.05, 0.1, r'$z = %s$'%str(round(sim_z, 2)),
                            transform=ax[i,j].transAxes, fontsize=20,verticalalignment='top',
                            color='black')
            if isinstance(extra_labels,list):
                ax[i,j].text(0.40, 0.95, extra_labels[ipd],
                            transform=ax[i,j].transAxes, fontsize=14,verticalalignment='top',
                            color='black')

            if powell:
                ax[i,j].text(1e-29, 1e6, 'HD', fontsize=16,verticalalignment='top',
                            color='black')
                ax[i,j].text(1e-29, 5e4, 'WD', fontsize=16,verticalalignment='top',
                            color='black')
                ax[i,j].text(1e-28, 1e3, 'CD', fontsize=16,verticalalignment='top',
                            color='black')
                ax[i,j].text(1e-24, 3e2, 'F', fontsize=16,verticalalignment='top',
                            color='white')
                ax[i,j].text(4e-22, 1e6, 'CL', fontsize=16,verticalalignment='top',
                            color='black')
                ax[i,j].plot([1e-30,1e-23],[2e5,2e5],color='k',linewidth=1)
                ax[i,j].plot([1e-30,1e-23],[2e4,2e4],color='k',linewidth=1)
                ax[i,j].plot([1e-25,1e-25],[0,2e4],color='k',linewidth=1)
                ax[i,j].plot([1e-23,1e-23],[0,1e8],color='k',linewidth=1)
            
            

            if stats == 'mean':
                y_mean = np.zeros(len(x))
                z = z.T
                for k in range(0, len(x)):
                    a = np.nansum(y * z[:,k])
                    b = np.nansum(z[:,k])
                    y_mean[k] = a/b
                ax[i,j].plot(x,y_mean,color='w',linewidth=2, linestyle='--')
            elif stats == 'median':
                y_median = np.zeros(len(x))
                z = z.T
                for k in range(0, len(x)):
                    y_median[k] = np.nanmedian(y * z[:,k])
                ax[i,j].plot(x,y_median,color='w',linewidth=2, linestyle='--')
            
            if ipd==0:
                if layout == 'compact':
                    cbaxes = inset_axes(ax[i,j], width="200%", height="5%", loc='upper left',
                                        bbox_to_anchor=(0.0, 0., 1.0, 1.05),
                                        bbox_transform=ax[i,j].transAxes,borderpad=0)
                elif layout == 'extended':
                    width_ex = str(int(100*ncol))
                    cbaxes = inset_axes(ax[i,j], width=width_ex+"%", height="5%", loc='upper left',
                                        bbox_to_anchor=(0.0, 0., 1.0, 1.05),
                                        bbox_transform=ax[i,j].transAxes,borderpad=0)
                cbar = fig.colorbar(plot, cax=cbaxes, orientation='horizontal')
                cbar.set_label(plotting_z['label'],fontsize=20)
                cbar.ax.tick_params(labelsize=10)
                cbaxes.xaxis.set_label_position('top')
                cbaxes.xaxis.set_ticks_position('top')

            # If we want to include countor for outflows and inflows
            if doflows:
                # Outflow
                outpd = pds[ipd][1]
                zoutflow = np.array(outpd.zdata['hydro'][field_indexes[i]][:,:,weight_indexes[i]].d,order='F')
                zoutflow = outpd.obj.array(zoutflow,code_units_z)
                ztot = np.sum(zoutflow)
                n = 1000
                t = np.linspace(0, zoutflow.max(), n)
                integral = ((zoutflow >= t[:, None, None]) * zoutflow).sum(axis=(1,2))
                from scipy import interpolate
                f = interpolate.interp1d(integral, t)
                t_contours = f(np.array([0.8*ztot]))
                ax[i,j].contour(x, y, zoutflow.T, t_contours, colors='darkorange', linewidths=2)
                # Inflow
                inpd = pds[ipd][2]
                zinflow = np.array(inpd.zdata['hydro'][field_indexes[i]][:,:,weight_indexes[i]].d,order='F')
                zinflow = inpd.obj.array(zinflow,code_units_z)
                ztot = np.sum(zinflow)
                t = np.linspace(0, ztot, n)
                integral = ((zinflow >= t[:, None, None]) * zinflow).sum(axis=(1,2))
                f = interpolate.interp1d(integral, t)
                t_contours = f(np.array([0.8*ztot]))
                ax[i,j].contour(x, y, zinflow.T, t_contours, colors='darkblue', linewidths=2)
                # Escaping
                espd = pds[ipd][3]
                zescape = np.array(espd.zdata['hydro'][field_indexes[i]][:,:,weight_indexes[i]].d,order='F')
                zescape = espd.obj.array(zescape,code_units_z)
                ztot = np.sum(zescape)
                t = np.linspace(0, ztot, n)
                integral = ((zescape >= t[:, None, None]) * zescape).sum(axis=(1,2))
                f = interpolate.interp1d(integral, t)
                t_contours = f(np.array([0.8*ztot]))
                ax[i,j].contour(x, y, zescape.T, t_contours, colors='darkred', linewidths=2)
            if do_sf:
                from scipy import interpolate
                n = 1000
                expd_index = 1
                if doflows:
                    expd_index += 3
                # Star forming
                stf = pds[ipd][expd_index].zdata['hydro'][field_indexes[i]][:,:,weight_indexes[i],0]
                ztot_sf = np.sum(stf)
                t = np.linspace(0, ztot_sf, n)
                integral = ((stf >= t[:, None, None]) * stf).sum(axis=(1,2))
                f = interpolate.interp1d(integral, t)
                t_contours = f(np.array([0.8*ztot_sf]))
                ax[i,j].contour(x, y, stf.T, t_contours, colors='black', linewidths=2)
                print('Star forming: %.3f, %.3f'%(100*ztot_sf/z_tot, ztot_sf))
            
        
        if layout == 'compact':
            fig.subplots_adjust(top=0.92,bottom=0.05,left=0.1,right=0.95)
        elif layout == 'extended':
            fig.subplots_adjust(top=0.85,bottom=0.12,left=0.07,right=0.98)
        if returnfig:
            return fig,ax
        else:
            fig.savefig(name+'.png',format='png',dpi=300)

def plot_compare_stacked_pd(pds,weights,field,name,weightvar='cumulative',
                            scaletype='log_even',redshift=True,stats='none',
                            extra_labels='none',gent=False,powell=False,
                            doflows=False,do_sf=False,layout='compact',
                            returnfig=False):

    # Make required imports
    from astropy.cosmology import FlatLambdaCDM
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import matplotlib.patheffects as pe
    from matplotlib.colors import LogNorm
    import seaborn as sns
    from ozy.variables_settings import plotting_dictionary
    sns.set(style="white")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    hfont = {'fontname':'Helvetica'}
    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', **{'family' : "serif"})
    
    npds = len(pds)
    # Check that the required field is actually in the PDs provided
    field_indexes = []
    weight_indexes = []
    for i,pd in enumerate(pds):
        field_indexes.append([])
        weight_indexes.append([])
        if doflows or do_sf:
            for pd2 in pd:
                if field[i] not in pd2[0].zvars['hydro']:
                    raise KeyError('The field %s is not included in this PhaseDiagram object. Aborting!'%field[i])
                else:
                    for f in range(0,len(pd2[0].zvars['hydro'])):
                        if pd2[0].zvars['hydro'][f] == field[i]:
                            field_indexes[-1].append(f)
                            break
                    print(pd2[0].weightvars['hydro'])
                    if weightvar not in pd2[0].weightvars['hydro']:
                        raise KeyError('The weight field %s is not included in this PhaseDiagram object. Aborting!'%weightvar)
                    else:
                        for w in range(0, len(pd2[0].weightvars['hydro'])):
                            if pd2[0].weightvars['hydro'][w] == weightvar:
                                weight_indexes[-1].append(w)
                                break
        else:
            for pd2 in pd:
                if field[i] not in pd2.zvars['hydro']:
                    raise KeyError('The field %s is not included in this PhaseDiagram object. Aborting!'%field[i])
                else:
                    for f in range(0,len(pd2.zvars['hydro'])):
                        if pd2.zvars['hydro'][f] == field[i]:
                            field_indexes[-1].append(f)
                            break
                    if weightvar not in pd2.weightvars['hydro']:
                        raise KeyError('The weight field %s is not included in this PhaseDiagram object. Aborting!'%weightvar)
                    else:
                        for w in range(0, len(pd2.weightvars['hydro'])):
                            if pd2.weightvars['hydro'][w] == weightvar:
                                weight_indexes[-1].append(w)
                                break
    if len(field_indexes) != npds or len(weight_indexes) != npds:
        print('You should check the fields and weights available, not all your phase diagrams have them!')
        exit
    # With everything fine, we begin plotting…
    if layout == 'compact':
        ncol = 2
        nrow = int(len(pds)/ncol)
    elif layout == 'extended':
        nrow = 1
        ncol = len(pds)
    else:
        print('The layout asked is not allowed. Please check!')
        exit
    
    if layout == 'compact':
        figsize = plt.figaspect(float((5.0 * nrow) / (5.0 * ncol)))
        fig = plt.figure(figsize=2*figsize, facecolor='w',edgecolor='k')
    elif layout == 'extended':
        figsize = plt.figaspect(float((6.0 * nrow) / (5.0 * ncol)))
        fig = plt.figure(figsize=figsize, facecolor='w',edgecolor='k')
    plot_grid = fig.add_gridspec(nrow, ncol, wspace=0, hspace=0)#,  hspace=0,left=0,right=1, bottom=0, top=1)
    ax = []
    for i in range(0,nrow):
        ax.append([])
        for j in range(0,ncol):
            ax[i].append(fig.add_subplot(plot_grid[i,j]))
    ax = np.asarray(ax)
    for i in range(0, ax.shape[0]):
        for j in range(0, ax.shape[1]):
            ipd = i*ax.shape[1] + j
            if ipd >= npds:
                # Clear that extra panel
                ax[i,j].get_xaxis().set_visible(False)
                ax[i,j].get_yaxis().set_visible(False)
                break
            pd = pds[ipd]
            pd_weight = weights[ipd]
            if doflows or do_sf:
                plotting_x = plotting_dictionary[pd[0][0].xvar]
                plotting_y = plotting_dictionary[pd[0][0].yvar]
            else:
                plotting_x = plotting_dictionary[pd[0].xvar]
                plotting_y = plotting_dictionary[pd[0].yvar]
            plotting_z = plotting_dictionary[field[ipd]]
            ax[i,j].set_xlabel(plotting_x['label'],fontsize=18)
            ax[i,j].xaxis.set_ticks_position('both')
            ax[i,j].yaxis.set_ticks_position('left')
            ax[i,j].minorticks_on()
            ax[i,j].set_xlim([1e-30,8e-20])
            ax[i,j].set_ylim([15,1e+8])
            ax[i,j].set_xscale('log')
            ax[i,j].set_yscale('log')
            # Get rid of the y-axis labels for the panels in the middle
            # TODO: The ticks should not be only for density and temperature!
            if ipd%2 == 0 and layout == 'compact':
                ax[i,j].tick_params(labelsize=14)
                ax[i,j].set_ylabel(plotting_y['label'],fontsize=18)
            elif layout == 'compact':
                ax[i,j].tick_params(labelsize=14)
                ax[i,j].axes.get_yaxis().set_ticks([])
                ax[i,j].set_xticks([1e-28,1e-26,1e-24,1e-22,1e-20])
            elif ipd != 0 and layout == 'extended':
                ax[i,j].set_yticks([1e2,1e3,1e4,1e5,1e6,1e7,1e8])
                ax[i,j].set_yticklabels([r"",r"$10^3$",r"$10^4$",r"$10^5$",r"$10^6$",r"$10^7$",""])
                ax[i,j].tick_params(labelsize=14,which='major',axis="y",direction="out")
                ax[i,j].tick_params(which='minor',axis="y",direction="out")
                ax[i,j].set_xticks([1e-28,1e-26,1e-24,1e-22,1e-20])
            else:
                ax[i,j].set_ylabel(plotting_y['label'],fontsize=18)
                ax[i,j].tick_params(labelsize=14,which='major',axis="y",direction="out")
                ax[i,j].tick_params(which='minor',axis="y",direction="out")

            
            if doflows or do_sf:
                code_units_x = get_code_units(pd[0][0].xvar)
            else:
                code_units_x = get_code_units(pd[0].xvar)

            z = np.zeros((100,100))
            sim_z = np.zeros(len(pd_weight))
            sim_t = np.zeros(len(pd_weight))
            tot_weight = 0
            for w in range(0, len(pd_weight)):
                if doflows or do_sf:
                    temp_pd = pd[w][0]
                    temp_weight = pd_weight[w][0]
                else:
                    temp_pd = pd[w]
                    temp_weight = pd_weight[w]
                x = temp_pd.xdata[0]
                print(x)
                x = x.in_units(plotting_x['units'])
                code_units_y = get_code_units(temp_pd.yvar)
                y = temp_pd.ydata[0]
                y = y.in_units(plotting_y['units'])
                code_units_z = get_code_units(field[ipd])
                XX,YY = np.meshgrid(x,y)
                ztemp = np.array(temp_pd.zdata['hydro'][field_indexes[ipd][w]][:,:,weight_indexes[ipd][w],0].d,order='F')
                ztemp = np.nan_to_num(ztemp, nan=0)
                ztemp = temp_pd.obj.array(ztemp,code_units_z).in_units(plotting_z['units'])
                z = z + ztemp.d*temp_weight
                tot_weight = tot_weight + temp_weight
                sim_z[w] = temp_pd.obj.simulation.redshift
                h = temp_pd.obj.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=temp_pd.obj.simulation.hubble_constant, Om0=temp_pd.obj.simulation.omega_matter, 
                                        Ob0=temp_pd.obj.simulation.omega_baryon,Tcmb0=2.73)
                sim_t[w] = cosmo.age(sim_z[w]).value
            z = z / tot_weight
            XX,YY = np.meshgrid(x,y)
            delta_t = sim_t.max() - sim_t.min()
            sim_z = np.mean(sim_z)
            if gent:
                XX,YY = np.meshgrid(x,y)
                ax[i,j].fill_between([1e-30,8e-20], [gent_curve_T('cold',1e-30),gent_curve_T('cold',8e-20)],
                                    [1,1],color='b', zorder=1,alpha=0.2)
                ax[i,j].fill_between([1e-30,8e-20], [gent_curve_T('cold',1e-30),gent_curve_T('cold',8e-20)],
                                    [gent_curve_T('hot',1e-30),gent_curve_T('hot',8e-20)],color='orange', 
                                    zorder=1,alpha=0.2)
                ax[i,j].fill_between([1e-30,8e-20], [gent_curve_T('hot',1e-30),gent_curve_T('hot',8e-20)],
                                [1e8,1e8],color='r', zorder=1,alpha=0.2)
                if field == 'mass':
                    z_cold = np.nansum(z.T[YY<gent_curve_T('cold',XX)])
                    z_hot = np.nansum(z.T[YY>gent_curve_T('hot',XX)])
                    z_warm = np.nansum(z.T[(YY>gent_curve_T('cold',XX)) & (YY<gent_curve_T('hot',XX))])
                    z_tot = np.nansum(z)
                    print(z_cold,z_warm,z_hot)
                    print('Distribution of masses in the Gent phases in %s:'%extra_labels[ipd])
                    print('Cold: %.3f, %.3e'%(100*z_cold/z_tot,z_cold))
                    print('Warm: %.3f, %.3e'%(100*z_warm/z_tot,z_warm))
                    print('Hot: %.3f, %.3e'%(100*z_hot/z_tot, z_hot))
            if scaletype=='log_even':
                plot = ax[i,j].pcolormesh(x,y,
                                    z.T,
                                    shading='auto',
                                    cmap=sns.color_palette("Spectral", as_cmap=True),
                                    norm=LogNorm(vmin=plotting_z['vmin_galaxy'],
                                    vmax=plotting_z['vmax_galaxy']),
                                    rasterized=True,
                                    antialiased=True)
            else:
                plot = ax[i,j].pcolormesh(x,y,
                                    z.T,
                                    shading='auto',
                                    cmap=plotting_z['cmap'],
                                    vmin=plotting_z['vmin'],
                                    vmax=plotting_z['vmax'],
                                    rasterized=True,
                                    antialiased=True)
            if redshift:
                ax[i,j].text(0.05, 0.2, r'$\langle z \rangle = %s$'%str(round(sim_z, 2))+'\n'+
                                        r'$\delta t = %s$ Gyr'%str(round(delta_t, 3)),
                            transform=ax[i,j].transAxes, fontsize=16,verticalalignment='top',
                            color='black',path_effects=[pe.withStroke(linewidth=2, foreground="white")])
            if isinstance(extra_labels,list):
                ax[i,j].text(0.7, 0.9, extra_labels[ipd],
                            transform=ax[i,j].transAxes, fontsize=14,verticalalignment='top',
                            color='black')
            if powell:
                ax[i,j].text(1e-29, 1e6, 'HD', fontsize=16,verticalalignment='top',
                            color='black')
                ax[i,j].text(1e-29, 5e4, 'WD', fontsize=16,verticalalignment='top',
                            color='black')
                ax[i,j].text(1e-28, 1e3, 'CD', fontsize=16,verticalalignment='top',
                            color='black')
                ax[i,j].text(1e-24, 3e2, 'F', fontsize=16,verticalalignment='top',
                            color='white')
                ax[i,j].text(4e-22, 1e6, 'CL', fontsize=16,verticalalignment='top',
                            color='black')
                ax[i,j].plot([1e-30,1e-23],[2e5,2e5],color='k',linewidth=1)
                ax[i,j].plot([1e-30,1e-23],[2e4,2e4],color='k',linewidth=1)
                ax[i,j].plot([1e-25,1e-25],[0,2e4],color='k',linewidth=1)
                ax[i,j].plot([1e-23,1e-23],[0,1e8],color='k',linewidth=1)

            if stats == 'mean':
                z_tot = np.sum(z)
                y_mean = np.zeros(len(x))
                n_points = np.zeros(len(x))
                z = z.T
                for k in range(0, len(x)):
                    a = np.nansum(y * z[:,k])
                    b = np.nansum(z[:,k])
                    y_mean[k] = a/b
                    n_points[k] = np.nansum(z[:,k])/z_tot
                ax[i,j].plot(x[n_points>1e-5],y_mean[n_points>1e-5],color='w',linewidth=2, linestyle='--')
            
            if ipd==0:
                if layout == 'compact':
                    cbaxes = inset_axes(ax[i,j], width="200%", height="5%", loc='upper left',
                                        bbox_to_anchor=(0.0, 0., 1.0, 1.05),
                                        bbox_transform=ax[i,j].transAxes,borderpad=0)
                elif layout == 'extended':
                    width_ex = str(int(100*ncol))
                    cbaxes = inset_axes(ax[i,j], width=width_ex+"%", height="5%", loc='upper left',
                                        bbox_to_anchor=(0.0, 0., 1.0, 1.05),
                                        bbox_transform=ax[i,j].transAxes,borderpad=0)
                cbar = fig.colorbar(plot, cax=cbaxes, orientation='horizontal')
                cbar.set_label(plotting_z['label'],fontsize=20)
                cbar.ax.tick_params(labelsize=14)
                cbaxes.xaxis.set_label_position('top')
                cbaxes.xaxis.set_ticks_position('top')
            # If we want to include countor for outflows and inflows
            if doflows:
                from scipy import interpolate
                # Outflow
                zoutflow = np.zeros((100,100))
                tot_weight = 0
                for w in range(0, len(weights[ipd])):
                    outpd = pds[ipd][w][1]
                    zoutflow_temp = np.array(outpd.zdata['hydro'][field_indexes[ipd][w]][:,:,weight_indexes[ipd][w],0].d,order='F')
                    zoutflow = zoutflow + outpd.obj.array(zoutflow_temp,code_units_z).in_units(plotting_z['units']).d*weights[ipd][w][1]
                    tot_weight = tot_weight + weights[ipd][w][1]
                zoutflow = zoutflow / tot_weight
                ztot = np.sum(zoutflow)
                n = 1000
                t = np.linspace(0, zoutflow.max(), n)
                integral = ((zoutflow >= t[:, None, None]) * zoutflow).sum(axis=(1,2))
                f = interpolate.interp1d(integral, t)
                t_contours = f(np.array([0.8*ztot]))
                ax[i,j].contour(x, y, zoutflow.T, t_contours, colors='darkorange', linewidths=2)
                # Inflow
                zinflow = np.zeros((100,100))
                tot_weight = 0
                for w in range(0, len(weights[ipd])):
                    inpd = pds[ipd][w][2]
                    zinflow_temp = np.array(inpd.zdata['hydro'][field_indexes[ipd][w]][:,:,weight_indexes[ipd][w],0].d,order='F')
                    zinflow = zinflow + inpd.obj.array(zinflow_temp,code_units_z).in_units(plotting_z['units']).d*weights[ipd][w][2]
                    tot_weight = tot_weight + weights[ipd][w][2]
                zinflow = zinflow / tot_weight
                ztot = np.sum(zinflow)
                t = np.linspace(0, ztot, n)
                integral = ((zinflow >= t[:, None, None]) * zinflow).sum(axis=(1,2))
                f = interpolate.interp1d(integral, t)
                t_contours = f(np.array([0.8*ztot]))
                ax[i,j].contour(x, y, zinflow.T, t_contours, colors='darkblue', linewidths=2)
                # Escaping
                zescape = np.zeros((100,100))
                tot_weight = 0
                for w in range(0, len(weights[ipd])):
                    espd = pds[ipd][w][3]
                    zescape_temp = np.array(espd.zdata['hydro'][field_indexes[ipd][w]][:,:,weight_indexes[ipd][w],0].d,order='F')
                    zescape = zescape + espd.obj.array(zescape_temp,code_units_z).in_units(plotting_z['units']).d*weights[ipd][w][3]
                    tot_weight = tot_weight + weights[ipd][w][3]
                zescape = zescape / tot_weight
                ztot = np.sum(zescape)
                t = np.linspace(0, ztot, n)
                integral = ((zescape >= t[:, None, None]) * zescape).sum(axis=(1,2))
                f = interpolate.interp1d(integral, t)
                t_contours = f(np.array([0.8*ztot]))
                ax[i,j].contour(x, y, zescape.T, t_contours, colors='darkred', linewidths=2)
            if do_sf:
                from scipy import interpolate
                n = 1000
                expd_index = 1
                if doflows:
                    expd_index += 3
                # Star forming
                zsf = np.zeros((100,100))
                tot_weight = 0
                for w in range(0, len(weights[ipd])):
                    espd = pds[ipd][w][expd_index]
                    zsf_temp = np.array(espd.zdata['hydro'][field_indexes[ipd][w]][:,:,weight_indexes[ipd][w],0].d,order='F')
                    zsf = zsf + espd.obj.array(zsf_temp,code_units_z).in_units(plotting_z['units']).d*weights[ipd][w][expd_index]
                    tot_weight = tot_weight + weights[ipd][w][expd_index]
                zsf = zsf / tot_weight
                ztot_sf = np.sum(zsf)
                t = np.linspace(0, ztot_sf, n)
                integral = ((zsf >= t[:, None, None]) * zsf).sum(axis=(1,2))
                f = interpolate.interp1d(integral, t)
                t_contours = f(np.array([0.8*ztot_sf]))
                ax[i,j].contour(x, y, zsf.T, t_contours, colors='black', linewidths=2)
                print('Star forming: %.3f, %.3f'%(100*ztot_sf/z_tot, ztot_sf))
                
        if layout == 'compact':
            fig.subplots_adjust(top=0.92,bottom=0.05,left=0.1,right=0.95)
        elif layout == 'extended':
            fig.subplots_adjust(top=0.83,bottom=0.13,left=0.07,right=0.98)
        if returnfig:
            return fig, ax
        else:
            fig.savefig(name+'.pdf',format='pdf',dpi=300)
