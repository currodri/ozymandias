import numpy as np
import h5py
import os
from yt import YTArray
import ozy
from ozy.dict_variables import common_variables,grid_variables,particle_variables,get_code_units
# TODO: Allow for parallel computation of phase diagrams.
from joblib import Parallel, delayed
import sys
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/part')
from amr2 import vectors
from amr2 import geometrical_regions as geo
from amr2 import filtering
from amr2 import amr_profiles as amrprofmod
from ozy.saver import _write_attrib
blacklist = [
    'zvars','weightvars','xdata','ydata','zdata'
]

class PhaseDiagram(object):

    def __init__(self,group):
        self.group = group
        self.obj = group.obj
        self.nbins = [0,0]
        self.xvar = None
        self.yvar = None
        self.region = None
        self.filter = None
        self.lmax = 0
        self.zvars = {}
        self.weightvars = {}
        self.xdata = []
        self.ydata = []
        self.zdata = {}

    def _serialise(self, hdd):
        """This makes possible to save the group phase diagram attrs as dataset attributes of an HDF5 group."""
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
        """Save the Fortran derived type as a dictionary inside the PhaseDiagram class (only the necessary info)."""
        from yt import YTArray
        self.region = {}
        self.region['type'] = reg.name.decode().split(' ')[0]
        self.region['centre'] = YTArray([reg.centre.x, reg.centre.y, reg.centre.z], 'code_length', registry=self.obj.unit_registry)
        self.region['axis'] = YTArray([reg.axis.x, reg.axis.y, reg.axis.z], 'dimensionless', registry=self.obj.unit_registry)
    
    def _get_python_filter(self,filt):
        """Save the Fortran derived type as a dictionary inside the PhaseDiagram class (only the necessary info)."""
        from yt import YTQuantity
        self.filter = {}
        self.filter['name'] = filt.name.decode().split(' ')[0]
        self.filter['conditions'] = []
        if self.filter['name'] != 'none':
            for i in range(0, filt.ncond):
                cond_var = filt.cond_vars.T.view('S128')[i][0].decode().split(' ')[0]
                cond_op = filt.cond_ops.T.view('S2')[i][0].decode().split(' ')[0]
                cond_units = get_code_units(cond_var)
                cond_value = YTQuantity(filt.cond_vals[i], str(cond_units), registry=self.obj.unit_registry)
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
    if cond_strs[0] == 'none':
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
        raise ValueError("Condition strings are given, but not a name for the filter. Please set!")

def compute_phase_diagram(group,ozy_file,xvar,yvar,zvars,weightvars,lmax=0,nbins=[100,100],region_type='sphere',
                            filter_conds='none',filter_name='none',logscale=True,recompute=False,save=False):
    """Function which computes a phase diagram (2D profile) for a given group object."""

    if not isinstance(xvar,str) or not isinstance(yvar,str):
        print('Only single x and y variable supported!')
        exit
    if len(nbins) != 2 or not isinstance(nbins,list):
        print('Number of bins should be given for both x and y axis. No more, no less!')
        exit
    pd = PhaseDiagram(group)
    pd.nbins = nbins
    pd.xvar = xvar
    pd.yvar = yvar
    pd.zvars = dict(hydro = [],for_star = [], for_dm = [])
    pd.weightvars = dict(hydro = [],for_star = [], for_dm = [])
    # Begin by checking the combination of variables is correct
    # Variables should be given in the form "type"/"variable":
    # e.g. gas/density or star/ang_momentum_x
    for var in zvars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            if var_name in common_variables or var_name in grid_variables:
                pd.zvars['hydro'].append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!: %s', var)
        elif var_type == 'star':
            if var_name in common_variables or var_name in particle_variables:
                pd.zvars['for_star'].append(var_name)
            else:
                raise KeyError('This star variable is not supported. Please check!')
        elif var_type == 'dm':
            if var_name in common_variables or var_name in particle_variables:
                pd.zvars['for_dm'].append(var_name)
            else:
                raise KeyError('This DM variable is not supported. Please check!')
    for var in weightvars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            if var_name in common_variables or var_name in grid_variables:
                pd.weightvars['hydro'].append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!')
        elif var_type == 'star':
            if var_name in common_variables or var_name in particle_variables:
                pd.weightvars['for_star'].append(var_name)
            else:
                raise KeyError('This star variable is not supported. Please check!')
        elif var_type == 'dm':
            if var_name in common_variables or var_name in particle_variables:
                pd.weightvars['for_dm'].append(var_name)
            else:
                raise KeyError('This DM variable is not supported. Please check!')

    # Check that we do not have any inconsistency...
    if xvar in grid_variables and len(pd.yvars['for_star'])>0 or xvar in grid_variables and len(pd.yvars['for_dm'])>0:
        raise KeyError("Having grid vs particle phase diagrams is not well-defined.")
    elif xvar in particle_variables and len(pd.yvars['hydro'])>0:
        raise KeyError("Having particle vs grid phase diagrams is not well-defined.")

    # Now create region
    if isinstance(region_type, geo.region):
        selected_reg = region_type
    else:
        selected_reg = init_region(group,region_type)

    # Save region details to PhaseDiagram object
    pd._get_python_region(selected_reg)

    # Now create filter, if any conditions have been given…
    filt = init_filter(filter_conds, filter_name, group)
    # ...and save to PhaseDiagram object
    pd._get_python_filter(filt)

    # Check if phase data is already present and if it coincides with the new one
    f = h5py.File(ozy_file, 'r+')
    pd_present, pd_key = check_if_same_phasediag(f,pd)
    if pd_present and recompute:
        del f[str(pd.group.obj_type)+'_data/phase_diagrams/'+str(group._index)+'/'+str(pd_key)]
        print('Overwriting phase diagram data in %s_data'%group.obj_type)
    elif pd_present and not recompute:
        print('Phase diagram data with same details already present for galaxy %s. No overwritting!'%group._index)
        group._init_phase_diagrams()
        for i,p in enumerate(group.phase_diagrams):
            if p.key == pd_key:
                selected_pd = i
                break
        return group.phase_diagrams[selected_pd]
    else:
        print('Writing phase diagram data in %s_data'%group.obj_type)
    f.close()

    # Initialise hydro phase diagram data object
    hydro_data = amrprofmod.profile_handler_twod()
    hydro_data.profdim = 2
    hydro_data.xvarname = xvar
    hydro_data.yvarname = yvar
    hydro_data.nzvar = len(pd.zvars['hydro'])
    hydro_data.nwvar = len(pd.weightvars['hydro'])
    hydro_data.nbins = np.asarray(nbins,order='F')

    amrprofmod.allocate_profile_handler_twod(hydro_data)
    for i in range(0, len(pd.zvars['hydro'])):
        hydro_data.zvarnames.T.view('S128')[i] = pd.zvars['hydro'][i].ljust(128)
    for i in range(0, len(pd.weightvars['hydro'])):
        hydro_data.wvarnames.T.view('S128')[i] = pd.weightvars['hydro'][i].ljust(128)
    
    # And now, compute hydro data phase diagrams!
    if hydro_data.nzvar > 0 and hydro_data.nwvar > 0:
        amrprofmod.twodprofile(group.obj.simulation.fullpath,selected_reg,filt,hydro_data,lmax,logscale)
    
    copy_data = np.copy(hydro_data.xdata)
    code_units = get_code_units(pd.xvar)
    pd.xdata.append(YTArray(copy_data, code_units, registry=group.obj.unit_registry))

    code_units = get_code_units(pd.yvar)
    copy_data = np.copy(hydro_data.ydata)
    pd.ydata.append(YTArray(copy_data, code_units, registry=group.obj.unit_registry))
    pd.zdata['hydro'] = []
    for i in range(0, len(pd.zvars['hydro'])):
        code_units = get_code_units(pd.zvars['hydro'][i])
        copy_data = np.copy(hydro_data.zdata[:,:,i,:,0:2])
        pd.zdata['hydro'].append(YTArray(copy_data, code_units, registry=group.obj.unit_registry))

    # TODO: Add phase diagram for particles
    star_data = None
    dm_data = None
    if save:
        write_phasediag(group.obj, ozy_file, hydro_data, star_data, dm_data, pd)
    
    return pd

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
    if not str(pd.group.obj_type)+'_data/phase_diagrams/'+str(pd.group._index) in hd:
        return False, 'none'
    for p in hd[str(pd.group.obj_type)+'_data/phase_diagrams/'+str(pd.group._index)].keys():
        check_xvar = (p.split('|')[0] == pd.xvar+'_'+pd.yvar)
        check_filtername = (p.split('|')[1] == pd.filter['name'])
        check_regiontype = (p.split('|')[2] == pd.region['type'])
        if check_xvar and check_filtername and check_regiontype:
            return True, p
    return False, 'none'

def write_phasediag(obj,ozy_file,hydro,star,dm,pd):
    """This function writes the resulting phase diagram for this group to the original OZY HDF5 file."""

    f = h5py.File(ozy_file, 'r+')

    # Create group in HDF5 file
    try:
        phase_diagrams = f.create_group(str(pd.group.obj_type)+'_data/phase_diagrams/'+str(pd.group._index))
    except:
        phase_diagrams = f[str(pd.group.obj_type)+'_data/phase_diagrams/'+str(pd.group._index)]
    
    # Clean data and save to dataset
    pd_name = get_phasediag_name(phase_diagrams,pd)
    hdpd = phase_diagrams.create_group(pd_name)
    pd._serialise(hdpd)

    # Save x and y data
    xdata = np.zeros((3,pd.nbins[0]))
    if hydro != None:
        xdata[0,:] = hydro.xdata
    if star != None:
        xdata[1,:] = star.xdata
    if dm != None:
        xdata[2,:] = dm.xdata
    hdpd.create_dataset('xdata', data=xdata)
    hdpd['xdata'].attrs.create('units', get_code_units(pd.xvar))

    ydata = np.zeros((3,pd.nbins[1]))
    if hydro != None:
        ydata[0,:] = hydro.ydata
    if star != None:
        ydata[1,:] = star.ydata
    if dm != None:
        ydata[2,:] = dm.ydata
    hdpd.create_dataset('ydata', data=ydata)
    hdpd['ydata'].attrs.create('units', get_code_units(pd.yvar))

    # Save hydro z data
    if hydro != None:
        clean_hydro = hdpd.create_group('hydro')
        for v,var in enumerate(pd.zvars['hydro']):
            clean_hydro.create_dataset(var, data=hydro.zdata[:,:,v,:,0:2])
            clean_hydro[var].attrs.create('units', get_code_units(pd.zvars['hydro'][v]))
        clean_hydro.create_dataset('weightvars', data=pd.weightvars['hydro'][:])
    # TODO: Save particle data
    f.close()
    return

def plot_single_phase_diagram(pd,field,name,weightvar='cumulative',logscale=True,redshift=True,stats='none'):
    """This function uses the information from PhaseDiagram following the OZY format and
        plots following the OZY standards."""
    
    # Make required imports
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from matplotlib.colors import LogNorm
    import seaborn as sns
    from ozy.plot_settings import plotting_dictionary
    sns.set(style="white")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    hfont = {'fontname':'Helvetica'}
    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', **{'family' : "serif"})
    params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    matplotlib.rcParams.update(params)

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
    plotting_x = plotting_dictionary[pd.xvar]
    plotting_y = plotting_dictionary[pd.yvar]
    plotting_z = plotting_dictionary[field]
    ax.set_xlabel(plotting_x['label'],fontsize=18)
    ax.set_ylabel(plotting_y['label'],fontsize=18)

    ax.tick_params(labelsize=14,direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.minorticks_on()
    ax.tick_params(which='both',axis="both",direction="in")

    ax.set_xscale('log')
    ax.set_yscale('log')
    code_units_x = get_code_units(pd.xvar)
    x = YTArray(10**(pd.xdata[0].d),code_units_x,
                registry=pd.obj.unit_registry)
    x = x.in_units(plotting_x['units'])
    code_units_y = get_code_units(pd.yvar)
    y = YTArray(10**(pd.ydata[0].d),code_units_y,
                registry=pd.obj.unit_registry)
    y = y.in_units(plotting_y['units'])
    code_units_z = get_code_units(field)
    z = np.array(pd.zdata['hydro'][field_index][:,:,weight_index,0].d,order='F')
    z = YTArray(z,code_units_z,
            registry=pd.obj.unit_registry)
    sim_z = pd.obj.simulation.redshift
    if logscale:
        plot = ax.pcolormesh(x,y,
                            z.in_units(plotting_z['units']).T,
                            shading='auto',
                            cmap=plotting_z['cmap'],
                            norm=LogNorm(vmin=plotting_z['vmin'],
                            vmax=plotting_z['vmax']))
    else:
        plot = ax.pcolormesh(x,y,
                            z.in_units(plotting_z['units']).T,
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
        z = z.T
        for i in range(0, len(x)):
            a = np.nansum(y * z[:,i])
            b = np.nansum(z[:,i])
            y_mean[i] = a/b
        ax.plot(x,y_mean,color='r',linewidth=2)
    fig.subplots_adjust(top=0.97,bottom=0.1,left=0.1,right=0.99)
    fig.savefig(name+'.png',format='png',dpi=300)

<<<<<<< HEAD
def plot_compare_phase_diagram(pds,field,name,weightvar='cumulative',logscale=True,redshift=True,powell=False,gent=False,stats='none',extra_labels='none'):
=======
def plot_compare_phase_diagram(pds,field,name,weightvar='cumulative',logscale=True,redshift=True,stats='none',extra_labels='none'):
>>>>>>> 2fe29f6cd2a2b2f6a393973835311a85476523be

    # Make required imports
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from matplotlib.colors import LogNorm
    import seaborn as sns
    from ozy.plot_settings import plotting_dictionary
    sns.set(style="white")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    hfont = {'fontname':'Helvetica'}
    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', **{'family' : "serif"})
    params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    matplotlib.rcParams.update(params)

    # How many phase diagram datasets have been provided
    if isinstance(pds,list):
        npds = len(pds)
    else:
        raise TypeError('You need to provide all phase diagram datasets in a list!')
    
    # Check that the required field is actually in the PDs provided
    field_indexes = []
    weight_indexes = []
    for pd in pds:
        if field not in pd.zvars['hydro']:
            raise KeyError('The field %s is not included in this PhaseDiagram object. Aborting!'%field)
        else:
            for f in range(0,len(pd.zvars['hydro'])):
                if pd.zvars['hydro'][f] == field:
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
    figsize = plt.figaspect(5.0 / float(5 * npds))
    fig = plt.figure(figsize=figsize, facecolor='w',edgecolor='k')
    plot_grid = fig.add_gridspec(1, npds, wspace=0)#,  hspace=0,left=0,right=1, bottom=0, top=1)
    ax = []
    for j in range(0,npds):
        ax.append(fig.add_subplot(plot_grid[j]))
    ax = np.asarray(ax)
    for i in range(0, npds):
        pd = pds[i]
        plotting_x = plotting_dictionary[pd.xvar]
        plotting_y = plotting_dictionary[pd.yvar]
        plotting_z = plotting_dictionary[field]
        ax[i].set_xlabel(plotting_x['label'],fontsize=18)
        if i == 0:
            ax[i].set_ylabel(plotting_y['label'],fontsize=18)
        else:
            ax[i].axes.yaxis.set_visible(False)

        ax[i].tick_params(labelsize=14,direction='in')
        ax[i].xaxis.set_ticks_position('both')
        ax[i].yaxis.set_ticks_position('both')
        ax[i].minorticks_on()
        ax[i].tick_params(which='major',axis="both",direction="in")

        ax[i].set_xscale('log')
        ax[i].set_yscale('log')
        code_units_x = get_code_units(pd.xvar)

        x = YTArray(10**(pd.xdata[0].d),code_units_x,
                registry=pd.obj.unit_registry)
        x = x.in_units(plotting_x['units'])
        code_units_y = get_code_units(pd.yvar)
        y = YTArray(10**(pd.ydata[0].d),code_units_y,
                    registry=pd.obj.unit_registry)
        y = y.in_units(plotting_y['units'])
        code_units_z = get_code_units(field)
        z = np.array(pd.zdata['hydro'][field_indexes[i]][:,:,weight_indexes[i],0].d,order='F')
        z = YTArray(z,code_units_z,
                    registry=pd.obj.unit_registry)
        sim_z = pd.obj.simulation.redshift

        if logscale:
            plot = ax[i].pcolormesh(x,y,
                                z.in_units(plotting_z['units']).T,
                                shading='auto',
                                cmap=plotting_z['cmap'],
                                norm=LogNorm(vmin=plotting_z['vmin'],
                                vmax=plotting_z['vmax']))
        else:
            plot = ax[i].pcolormesh(x,y,
                                z.in_units(plotting_z['units']).T,
                                shading='auto',
                                cmap=plotting_z['cmap'],
                                vmin=plotting_z['vmin'],
                                vmax=plotting_z['vmax'])
        if redshift:
            ax[i].text(0.05, 0.1, 'z = '+str(round(sim_z, 2)),
                        transform=ax[i].transAxes, fontsize=20,verticalalignment='top',
                        color='black')
        if isinstance(extra_labels,list):
            ax[i].text(0.5, 0.9, extra_labels[i],
                        transform=ax[i].transAxes, fontsize=20,verticalalignment='top',
                        color='black')

        if powell:
            ax[i].text(1e-29, 1e6, 'HD', fontsize=16,verticalalignment='top',
                        color='black')
            ax[i].text(1e-29, 5e4, 'WD', fontsize=16,verticalalignment='top',
                        color='black')
            ax[i].text(1e-28, 1e3, 'CD', fontsize=16,verticalalignment='top',
                        color='black')
            ax[i].text(1e-24, 3e2, 'F', fontsize=16,verticalalignment='top',
                        color='white')
            ax[i].text(4e-22, 1e6, 'CL', fontsize=16,verticalalignment='top',
                        color='black')
            ax[i].plot([1e-30,1e-23],[2e5,2e5],color='k',linewidth=1)
            ax[i].plot([1e-30,1e-23],[2e4,2e4],color='k',linewidth=1)
            ax[i].plot([1e-25,1e-25],[0,2e4],color='k',linewidth=1)
            ax[i].plot([1e-23,1e-23],[0,1e8],color='k',linewidth=1)
        
        if gent:
            ax[i].plot([1e-30,1e-23],[2e5,2e5],color='k',linewidth=1)
            ax[i].plot([1e-30,1e-23],[2e4,2e4],color='k',linewidth=1)

        if stats == 'mean':
            y_mean = np.zeros(len(x))
            z = z.T
            for k in range(0, len(x)):
                a = np.nansum(y * z[:,k])
                b = np.nansum(z[:,k])
                y_mean[k] = a/b
            ax[i].plot(x,y_mean,color='r',linewidth=2)
        
        if i==npds-1:
            cbaxes = inset_axes(ax[i], width="5%", height="100%", loc='lower left',
                                bbox_to_anchor=(1.05, 0., 1, 1),
                                bbox_transform=ax[i].transAxes,borderpad=0)
            cbar = fig.colorbar(plot, cax=cbaxes, orientation='vertical')
            cbar.set_label(plotting_z['label'],fontsize=20)
            cbar.ax.tick_params(labelsize=14)
            
        
        fig.subplots_adjust(top=0.97,bottom=0.12,left=0.07,right=0.88)
<<<<<<< HEAD
        fig.savefig(name+'.png',format='png',dpi=300)
=======
        fig.savefig(name+'.png',format='png',dpi=300)
>>>>>>> 2fe29f6cd2a2b2f6a393973835311a85476523be