import numpy as np
import h5py
import os
import ozy
from unyt import unyt_array,unyt_quantity
from ozy.utils import init_region,init_filter
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

def compute_phase_diagram(group,ozy_file,xvar,yvar,zvars,weightvars,lmax=0,nbins=[100,100],region_type='sphere',
                            filter_conds='none',filter_name='none',logscale=True,recompute=False,save=False,
                            rmin=(0.0,'rvir'), rmax=(0.2,'rvir'), zmin=(0.0,'rvir'), zmax=(0.2,'rvir')):
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
        selected_reg = init_region(group,region_type,rmin=rmin,rmax=rmax,zmin=zmin,zmax=zmax)

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
        del f[str(pd.group.type)+'_data/phase_diagrams/'+str(group._index)+'/'+str(pd_key)]
        print('Overwriting phase diagram data in %s_data'%group.type)
    elif pd_present and not recompute:
        print('Phase diagram data with same details already present for galaxy %s. No overwritting!'%group._index)
        group._init_phase_diagrams()
        for i,p in enumerate(group.phase_diagrams):
            if p.key == pd_key:
                selected_pd = i
                break
        return group.phase_diagrams[selected_pd]
    else:
        print('Writing phase diagram data in %s_data'%group.type)
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
    
    copy_data = np.copy(hydro_data.xdata[1:])
    code_units = get_code_units(pd.xvar)
    pd.xdata.append(group.obj.array(copy_data, code_units))

    code_units = get_code_units(pd.yvar)
    copy_data = np.copy(hydro_data.ydata[1:])
    pd.ydata.append(group.obj.array(copy_data, code_units))
    pd.zdata['hydro'] = []
    for i in range(0, len(pd.zvars['hydro'])):
        code_units = get_code_units(pd.zvars['hydro'][i])
        copy_data = np.copy(hydro_data.zdata[:,:,i,::2])
        pd.zdata['hydro'].append(group.obj.array(copy_data, code_units))

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
    if not str(pd.group.type)+'_data/phase_diagrams/'+str(pd.group._index) in hd:
        return False, 'none'
    for p in hd[str(pd.group.type)+'_data/phase_diagrams/'+str(pd.group._index)].keys():
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
        phase_diagrams = f.create_group(str(pd.group.type)+'_data/phase_diagrams/'+str(pd.group._index))
    except:
        phase_diagrams = f[str(pd.group.type)+'_data/phase_diagrams/'+str(pd.group._index)]
    
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
            clean_hydro.create_dataset(var, data=hydro.zdata[:,:,v,::2])
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
    ax.set_xlim([1e-30,1e-20])
    ax.set_ylim([2,1e+8])
    ax.set_xscale('log')
    ax.set_yscale('log')
    code_units_x = get_code_units(pd.xvar)
    x = pd.obj.array(10**(pd.xdata[0].d),code_units_x)
    x = x.in_units(plotting_x['units'])
    code_units_y = get_code_units(pd.yvar)
    y = pd.obj.array(10**(pd.ydata[0].d),code_units_y)
    y = y.in_units(plotting_y['units'])
    code_units_z = get_code_units(field)
    z = np.array(pd.zdata['hydro'][field_index][:,:,weight_index].d,order='F')
    z = pd.obj.array(z,code_units_z)
    sim_z = pd.obj.simulation.redshift
    if logscale:
        plot = ax.pcolormesh(x,y,
                            z[:,:,0].in_units(plotting_z['units']).T,
                            shading='auto',
                            cmap=plotting_z['cmap'],
                            norm=LogNorm(vmin=plotting_z['vmin'],
                            vmax=plotting_z['vmax']))
    else:
        plot = ax.pcolormesh(x,y,
                            z[:,:,0].in_units(plotting_z['units']).T,
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

    if gent:
        s1 = 4.4e+8
        s2 = 23.2e+8
        cv = 1.4e+8
        gamma = 5/3
        rho1_low = 1.673532784796145e-24 * (2/(np.exp(s1/cv))) ** (1/(gamma-1))
        rho1_high = 1.673532784796145e-24 * (1e+8/(np.exp(s1/cv))) ** (1/(gamma-1))
        ax.plot([rho1_low,rho1_high], [2,1e+8],color='k')
        rho2_low = 1.673532784796145e-24 * (2/(np.exp(s2/cv))) ** (1/(gamma-1))
        rho2_high = 1.673532784796145e-24 * (1e+8/(np.exp(s2/cv))) ** (1/(gamma-1))
        ax.plot([rho2_low,rho2_high], [2,1e+8],color='k')
        



    fig.subplots_adjust(top=0.97,bottom=0.1,left=0.1,right=0.99)
    fig.savefig(name+'.png',format='png',dpi=300)

def plot_compare_phase_diagram(pds,field,name,weightvar='cumulative',logscale=True,redshift=True,stats='none',extra_labels='none',gent=False,powell=False,doflows=False):

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
        if doflows:
            if field not in pd[0].zvars['hydro']:
                raise KeyError('The field %s is not included in this PhaseDiagram object. Aborting!'%field)
            else:
                for f in range(0,len(pd[0].zvars['hydro'])):
                    if pd[0].zvars['hydro'][f] == field:
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
    nrow = int(len(pds)/2)
    figsize = plt.figaspect(float((5.0 * nrow) / (5.0 * 2)))
    fig = plt.figure(figsize=2*figsize, facecolor='w',edgecolor='k')
    plot_grid = fig.add_gridspec(nrow, 2, wspace=0, hspace=0)#,  hspace=0,left=0,right=1, bottom=0, top=1)
    ax = []
    for i in range(0,nrow):
        ax.append([])
        for j in range(0,2):
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
            if doflows:
                pd = pds[ipd][0]
            else:
                pd = pds[ipd]
            plotting_x = plotting_dictionary[pd.xvar]
            plotting_y = plotting_dictionary[pd.yvar]
            plotting_z = plotting_dictionary[field]
            ax[i,j].set_xlabel(plotting_x['label'],fontsize=18)
            if ipd%2 == 0:
                ax[i,j].set_ylabel(plotting_y['label'],fontsize=18)
            else:
                ax[i,j].axes.yaxis.set_visible(False)

            ax[i,j].tick_params(labelsize=14,direction='in')
            ax[i,j].xaxis.set_ticks_position('both')
            ax[i,j].yaxis.set_ticks_position('both')
            ax[i,j].minorticks_on()
            ax[i,j].tick_params(which='major',axis="both",direction="in")
            
            code_units_x = get_code_units(pd.xvar)
            if logscale:
                ax[i,j].set_xscale('log')
                ax[i,j].set_yscale('log')
                x = pd.obj.array(10**(pd.xdata[0].d),code_units_x)
            else:
                x = pd.obj.array(pd.xdata[0].d,code_units_x)
            x = x.in_units(plotting_x['units'])
            code_units_y = get_code_units(pd.yvar)
            
            if logscale:
                y = pd.obj.array(10**(pd.ydata[0].d),code_units_y)
            else:
                y = pd.obj.array(pd.ydata[0].d,code_units_y)
            y = y.in_units(plotting_y['units'])
            code_units_z = get_code_units(field)
            z = np.array(pd.zdata['hydro'][field_indexes[i]][:,:,weight_indexes[i],0].d,order='F')
            z = pd.obj.array(z,code_units_z)
            sim_z = pd.obj.simulation.redshift
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
                ax[i,j].text(0.5, 0.9, extra_labels[ipd],
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
            
            if gent:
                ax[i,j].plot([1e-30,1e-23],[2e5,2e5],color='k',linewidth=1)
                ax[i,j].plot([1e-30,1e-23],[2e4,2e4],color='k',linewidth=1)

            if stats == 'mean':
                y_mean = np.zeros(len(x))
                z = z.T
                for k in range(0, len(x)):
                    a = np.nansum(y * z[:,k])
                    b = np.nansum(z[:,k])
                    y_mean[k] = a/b
                ax[i,j].plot(x,y_mean,color='k',linewidth=2, linestyle='--')
            
            if ipd==0:
                cbaxes = inset_axes(ax[i,j], width="200%", height="5%", loc='upper left',
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
            
        
        fig.subplots_adjust(top=0.92,bottom=0.05,left=0.1,right=0.95)
        fig.savefig(name+'.png',format='png',dpi=300)

def plot_compare_stacked_pd(pds,weights,field,name,weightvar='cumulative',logscale=True,redshift=True,stats='none',extra_labels='none',gent=False,powell=False,doflows=False):

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
    
    npds = len(pds)
    # Check that the required field is actually in the PDs provided
    field_indexes = []
    weight_indexes = []
    for pd in pds:
        if doflows:
            if field not in pd[0][0].zvars['hydro']:
                raise KeyError('The field %s is not included in this PhaseDiagram object. Aborting!'%field)
            else:
                for f in range(0,len(pd[0][0].zvars['hydro'])):
                    if pd[0][0].zvars['hydro'][f] == field:
                        field_indexes.append(f)
                        break
                if weightvar not in pd[0][0].weightvars['hydro']:
                    raise KeyError('The weight field %s is not included in this PhaseDiagram object. Aborting!'%weightvar)
                    exit
                else:
                    for w in range(0, len(pd[0][0].weightvars['hydro'])):
                        if pd[0][0].weightvars['hydro'][w] == weightvar:
                            weight_indexes.append(w)
                            break
        else:
            if field not in pd[0].zvars['hydro']:
                raise KeyError('The field %s is not included in this PhaseDiagram object. Aborting!'%field)
            else:
                for f in range(0,len(pd[0].zvars['hydro'])):
                    if pd[0].zvars['hydro'][f] == field:
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
    if len(field_indexes) != npds or len(weight_indexes) != npds:
        print('You should check the fields and weights available, not all your phase diagrams have them!')
        exit
    # With everything fine, we begin plotting…
    nrow = int(len(pds)/2)
    figsize = plt.figaspect(float((5.0 * nrow) / (5.0 * 2)))
    fig = plt.figure(figsize=2*figsize, facecolor='w',edgecolor='k')
    plot_grid = fig.add_gridspec(nrow, 2, wspace=0, hspace=0)#,  hspace=0,left=0,right=1, bottom=0, top=1)
    ax = []
    for i in range(0,nrow):
        ax.append([])
        for j in range(0,2):
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
            if doflows:
                plotting_x = plotting_dictionary[pd[0][0].xvar]
                plotting_y = plotting_dictionary[pd[0][0].yvar]
            else:
                plotting_x = plotting_dictionary[pd[0].xvar]
                plotting_y = plotting_dictionary[pd[0].yvar]
            plotting_z = plotting_dictionary[field]
            ax[i,j].set_xlabel(plotting_x['label'],fontsize=18)
            if ipd%2 == 0:
                ax[i,j].set_ylabel(plotting_y['label'],fontsize=18)
            else:
                ax[i,j].axes.yaxis.set_visible(False)
            ax[i,j].tick_params(labelsize=14,direction='in')
            ax[i,j].xaxis.set_ticks_position('both')
            ax[i,j].yaxis.set_ticks_position('both')
            ax[i,j].minorticks_on()
            ax[i,j].tick_params(which='major',axis="both",direction="in")

            ax[i,j].set_xscale('log')
            ax[i,j].set_yscale('log')
            if doflows:
                code_units_x = get_code_units(pd[0][0].xvar)
            else:
                code_units_x = get_code_units(pd[0].xvar)

            z = np.zeros((100,100))
            sim_z = np.zeros(len(pd_weight))
            tot_weight = 0
            for w in range(0, len(pd_weight)):
                if doflows:
                    temp_pd = pd[w][0]
                    temp_weight = pd_weight[w][0]
                else:
                    temp_pd = pd[w]
                    temp_weight = pd_weight[w]
                x = temp_pd.obj.array(10**(temp_pd.xdata[0].d),code_units_x)
                x = x.in_units(plotting_x['units'])
                code_units_y = get_code_units(temp_pd.yvar)
                y = temp_pd.obj.array(10**(temp_pd.ydata[0].d),code_units_y)
                y = y.in_units(plotting_y['units'])
                code_units_z = get_code_units(field)
                ztemp = np.array(temp_pd.zdata['hydro'][field_indexes[i]][:,:,weight_indexes[i],0].d,order='F')
                ztemp = temp_pd.obj.array(ztemp,code_units_z).in_units(plotting_z['units'])
                z = z + ztemp.d*temp_weight
                tot_weight = tot_weight + temp_weight
                sim_z[w] = temp_pd.obj.simulation.redshift
            z = z / tot_weight
            delta_z = sim_z.max() - sim_z.min()
            sim_z = np.mean(sim_z)
            if logscale:
                plot = ax[i,j].pcolormesh(x,y,
                                    z.T,
                                    shading='auto',
                                    cmap=plotting_z['cmap'],
                                    norm=LogNorm(vmin=plotting_z['vmin_galaxy'],
                                    vmax=plotting_z['vmax_galaxy']))
            else:
                plot = ax[i,j].pcolormesh(x,y,
                                    z.T,
                                    shading='auto',
                                    cmap=plotting_z['cmap'],
                                    vmin=plotting_z['vmin'],
                                    vmax=plotting_z['vmax'])
            if redshift:
                ax[i,j].text(0.05, 0.2, r'$\langle z \rangle = %s$'%str(round(sim_z, 2))+'\n'+
                                        r'$\delta z = %s$'%str(round(delta_z, 3)),
                            transform=ax[i,j].transAxes, fontsize=16,verticalalignment='top',
                            color='black')
            if isinstance(extra_labels,list):
                ax[i,j].text(0.5, 0.9, extra_labels[ipd],
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
            
            if gent:
                ax[i,j].plot([1e-30,1e-23],[2e5,2e5],color='k',linewidth=1)
                ax[i,j].plot([1e-30,1e-23],[2e4,2e4],color='k',linewidth=1)

            if stats == 'mean':
                y_mean = np.zeros(len(x))
                z = z.T
                for k in range(0, len(x)):
                    a = np.nansum(y * z[:,k])
                    b = np.nansum(z[:,k])
                    y_mean[k] = a/b
                ax[i,j].plot(x,y_mean,color='w',linewidth=2, linestyle='--')
            
            if ipd==0:
                cbaxes = inset_axes(ax[i,j], width="200%", height="5%", loc='upper left',
                                    bbox_to_anchor=(0.0, 0., 1.0, 1.05),
                                    bbox_transform=ax[i,j].transAxes,borderpad=0)
                cbar = fig.colorbar(plot, cax=cbaxes, orientation='horizontal')
                cbar.set_label(plotting_z['label'],fontsize=20)
                cbar.ax.tick_params(labelsize=10)
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
                    zoutflow_temp = np.array(outpd.zdata['hydro'][field_indexes[i]][:,:,weight_indexes[i],0].d,order='F')
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
                    zinflow_temp = np.array(inpd.zdata['hydro'][field_indexes[i]][:,:,weight_indexes[i],0].d,order='F')
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
                    zescape_temp = np.array(espd.zdata['hydro'][field_indexes[i]][:,:,weight_indexes[i],0].d,order='F')
                    zescape = zescape + espd.obj.array(zescape_temp,code_units_z).in_units(plotting_z['units']).d*weights[ipd][w][3]
                    tot_weight = tot_weight + weights[ipd][w][3]
                zescape = zescape / tot_weight
                ztot = np.sum(zescape)
                t = np.linspace(0, ztot, n)
                integral = ((zescape >= t[:, None, None]) * zescape).sum(axis=(1,2))
                f = interpolate.interp1d(integral, t)
                t_contours = f(np.array([0.8*ztot]))
                ax[i,j].contour(x, y, zescape.T, t_contours, colors='darkred', linewidths=2)
        
        fig.subplots_adjust(top=0.92,bottom=0.05,left=0.1,right=0.95)
        fig.savefig(name+'.png',format='png',dpi=300)
