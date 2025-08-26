import numpy as np
import h5py
import os
import ozy
from unyt import unyt_array,unyt_quantity
from ozy.dict_variables import common_variables,grid_variables,particle_variables,get_code_units
# TODO: Allow for parallel computation of flows.
from amr2 import amr_integrator,stats_utils,io_ramses
from ozy.utils import init_region,init_filter,structure_regions,get_code_bins,pdf_handler_to_stats
blacklist = [
    'data','weightvars'
]

class GalacticFlow(object):

    def __init__(self,group):
        self.group = group
        self.obj = group.obj
        self.type = None
        self.region = None
        self.filter = None
        self.data = {}
        self.weightvars = []
        self.rm_subs = False

    def _serialise(self,hdd):
        """This makes possible to save the group galactic flow attrs as dataset
            attributes of an HDF5 group.""" 
        
        for k,v in self.__dict__.items():
            if k in blacklist:
                continue
            if isinstance(v, (unyt_array, unyt_quantity)):
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
        self.region['r'] = self.obj.quantity(0.5*(reg.rmax+reg.rmin), 'code_length')
        self.region['dr'] = self.obj.quantity(reg.rmax-reg.rmin, 'code_length')
    
    def _get_python_filter(self,filt):
        """Save the Frotran derived type as a dictionary inside the GalacticFlow class."""
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

def get_flow_name(gf_group, gf, r, dr):
    """Create an individual galactic flow identifier name."""
    name = str(gf.type)
    name += '|'+str(gf.region['type'])
    name += '|'+str(r)
    name += '|'+str(dr)
    return name

def check_if_same_flow(hd,gf,r,dr):
    """This function looks at an OZY file for the flow data specified by gf for a particular object."""
    if gf.rm_subs:
        flow_key = '_data/flows_nosubs/'
    else:
        flow_key = '_data/flows/'
    if not str(gf.group.type)+flow_key+str(gf.group._index) in hd:
        return False, 'none'
    for g in hd[str(gf.group.type)+flow_key+str(gf.group._index)].keys():
        check_type = (g.split('|')[0] == gf.type)
        check_region = (g.split('|')[1] == gf.region['type'])
        check_r = (g.split('|')[2] == str(r))
        check_dr = (g.split('|')[3] == str(dr))
        if check_type and check_region and check_r and check_dr:
            # print(g, gf.type,gf.region['type'],str(r),str(dr))
            return True, g
    return False, 'none'

def write_flow(obj,ozy_file,gf,r,dr,verbose=False):
    """This function writes the resulting flow data to the original OZY file."""
    f = h5py.File(ozy_file, 'r+')
    if gf.rm_subs:
        flow_key = '_data/flows_nosubs/'
    else:
        flow_key = '_data/flows/'
    # Create group in HDF5 file
    try:
        flows = f.create_group(str(gf.group.type)+flow_key+str(gf.group._index))
    except:
        flows = f[str(gf.group.type)+flow_key+str(gf.group._index)]

    # Clean data and save to dataset
    gf_name = get_flow_name(flows,gf,r,dr)
    if verbose: print(f'Saving flow data {gf_name}...')
    hdgf = flows.create_group(gf_name)
    gf._serialise(hdgf)

    clean_data = hdgf.create_group('data')
    for var in gf.data:
        try:
            clean_data.create_dataset(var, data=gf.data[var].d)
        except:
            clean_data.create_dataset(var, data=gf.data[var])
        name = var.split('_'+r)[0]
        clean_data[var].attrs.create('units', get_code_units(name))
    
    f.close()
    return

def compute_flows(group,ozy_file,flow_type,rmin=(0.0,'rvir'), 
                  rmax=(1.0,'rvir'),recompute=False,save=False,
                  separate_phases=True,remove_subs=False,
                  verbose=False,pdf_bins=100):
    """Function which computes the analysis of galaxy wide flows, including outflows,
        inflows or AGN feedback."""

    use_neigh = True
    use_grav = True

    # Begin by setting up the necessary details of a given flow type
    if flow_type == 'outflow':
        reg = init_region(group, 'sphere', rmin=rmin, rmax=rmax)
        if separate_phases:
            all = init_filter(cond_strs='v_sphere_r/>/0/km*s**-1',name='all',group=group)
            # Powell et al. (2011)
            # clumpy = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','density/>/1e-23/g*cm**-3'],name='outflow_clumpy')
            # filament = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','density/>/1e-25/g*cm**-3','density/</1e-23/g*cm**-3','temperature/</2e4/K'],name='outflow_filament')
            # cold_diffuse = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','density/</1e-25/g*cm**-3','temperature/</2e4/K'],name='outflow_cold_diffuse')
            # cold_dense = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','density/>/1e-25/g*cm**-3','temperature/</2e4/K'],name='outflow_cold_dense')
            # warm = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','density/</1e-23/g*cm**-3','temperature/>/2e4/K','temperature/</2e5/K'],name='outflow_warm')
            # hot = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','density/</1e-23/g*cm**-3','temperature/>/2e5/K'],name='outflow_hot')

            # Sergio's suggestion
            hot = init_filter(cond_strs=['v_sphere_r/>/0/km*s**-1','temperature/>/1e5/K'],name='hot',group=group)
            warm_ionised = init_filter(cond_strs=['v_sphere_r/>/0/km*s**-1','temperature/</1e5/K','temperature/>/9e3/K'],name='warm_ionised',group=group)
            warm_neutral = init_filter(cond_strs=['v_sphere_r/>/0/km*s**-1','temperature/</9e3/K','temperature/>/1e3/K'],name='warm_neutral',group=group)
            cold = init_filter(cond_strs=['v_sphere_r/>/0/km*s**-1','temperature/</1e3/K'],name='cold',group=group)
            filt = [all,hot,warm_ionised,warm_neutral,cold]
        else:
            filt = init_filter(group,cond_strs='v_sphere_r/>/0/km*s**-1',name='outflow')
    elif flow_type == 'inflow':
        reg = init_region(group, 'sphere', rmin=rmin, rmax=rmax)
        if separate_phases:
            all = init_filter(cond_strs='v_sphere_r/<=/0/km*s**-1',name='all',group=group)
            # Powell et al. (2011)
            # clumpy = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/>/1e-23/g*cm**-3'],name='outflow_clumpy')
            # filament = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/>/1e-25/g*cm**-3','density/</1e-23/g*cm**-3','temperature/</2e4/K'],name='outflow_filament')
            # cold_diffuse = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/</1e-25/g*cm**-3','temperature/</2e4/K'],name='outflow_cold_diffuse')
            # cold_dense = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/>/1e-25/g*cm**-3','temperature/</2e4/K'],name='outflow_cold_dense')
            # warm = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/</1e-23/g*cm**-3','temperature/>/2e4/K','temperature/</2e5/K'],name='outflow_warm')
            # hot = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/</1e-23/g*cm**-3','temperature/>/2e5/K'],name='outflow_hot')

            hot = init_filter(cond_strs=['v_sphere_r/<=/0/km*s**-1','temperature/>/1e5/K'],name='hot',group=group)
            warm_ionised = init_filter(cond_strs=['v_sphere_r/<=/0/km*s**-1','temperature/</1e5/K','temperature/>/9e3/K'],name='warm_ionised',group=group)
            warm_neutral = init_filter(cond_strs=['v_sphere_r/<=/0/km*s**-1','temperature/</9e3/K','temperature/>/1e3/K'],name='warm_neutral',group=group)
            cold = init_filter(cond_strs=['v_sphere_r/<=/0/km*s**-1','temperature/</1e3/K'],name='cold',group=group)
            filt = [all,hot,warm_ionised,warm_neutral,cold]
        else:
            filt = init_filter(cond_strs='v_sphere_r/<=/0/km*s**-1',name='inflow',group=group)
    else:
        raise ValueError(f"This type of galactic flow is not supported: {flow_type}. Please check!")
    
    d_key = str(int(100*0.5*(rmax[0]+rmin[0])))
    shell_width = reg.rmax-reg.rmin
    width_key = str(round(100*(rmax[0]-rmin[0]),6))
    width_key = width_key.replace('.', '_')
    if verbose: print('width_key,shell_width: ',width_key,shell_width)
    shell_r = 0.5*(reg.rmax+reg.rmin)

    # Since the number of global quantities that can be computed
    # from the gas data depends on the specific configuration
    # of a simulation, this needs to be determined at the start
    # See: ozy/sim_attributes.py/assign_attributes

    output_path = group.obj.simulation.fullpath
    nvar = 0
    quantity_names = []
    do_binning = []
    weight_names = ['cumulative','volume','momentum_sphere_r']

    if group.obj.simulation.physics['hydro']:
        quantity_names += ['mass','density','temperature',
                            'momentum_sphere_r','v_sphere_r','v_tangential',
                            'thermal_energy','thermal_energy_specific',
                            'grav_therpfrsphere','grav_therpfrspherepos',
                            'grav_therpfrsphereneg']
        do_binning += [False,True,True,
                        False,True,True,
                        False,True,
                        True,True,
                        True]
    if group.obj.simulation.physics['metals']:
        quantity_names += ['metallicity']
        do_binning += [True]
    if group.obj.simulation.physics['magnetic']:
        quantity_names += ['magnetic_energy','magnetic_energy_specific',
                           'grav_magpfrsphere','grav_magpfrspherepos',
                           'grav_magpfrsphereneg']
        do_binning += [False,True,True,True,True]
        if not group.obj.simulation.physics['cr']:
            quantity_names += ['grav_totpfrsphere',
                                'grav_totpfrspherepos',
                                'grav_totpfrsphereneg']
            do_binning += [True,True,True]
    if group.obj.simulation.physics['cr']:
        quantity_names += ['cr_energy','cr_energy_specific',
                            'grav_crpfrsphere','grav_crpfrspherepos',
                            'grav_crpfrsphereneg','grav_totpfrsphere',
                            'grav_totpfrspherepos','grav_totpfrsphereneg']
        do_binning += [False, True, True,True, True,True,True,True]

    nvar = len(quantity_names)
    if group.obj.simulation.physics['rt']:
        quantity_names += ['xHII','xHeII','xHeIII']
        do_binning += [True, True, True]

    # Initialise GalacticFlow object
    gf = GalacticFlow(group)
    gf.type = flow_type
    gf.weightvars = weight_names
    gf.rm_subs = remove_subs

    # Save region and filter details to GalacticFlow object
    gf._get_python_region(reg)
    if separate_phases:
        gf._get_python_filter(filt[0])
    else:
        gf._get_python_filter(filt)

    # Check if flow data is already present and if it coincides with the new one
    f = h5py.File(ozy_file, 'r+')
    gf_present, gf_key = check_if_same_flow(f,gf,d_key,width_key)
    gfp_keys = []
    if gf_present and recompute:
        if remove_subs:
            del f[str(gf.group.type)+'_data/flows_nosubs/'+str(gf.group._index)+'/'+str(gf_key)]
        else:
            del f[str(gf.group.type)+'_data/flows/'+str(gf.group._index)+'/'+str(gf_key)]

        if verbose: print('Overwriting flow data in %s_data'%group.type)
    elif gf_present and not recompute:
        if verbose: print('Flow data with same details (%s) already present for galaxy %s. No overwritting!'%(gf_key,group._index))
        group._init_flows()
        if remove_subs:
            if verbose: print('Removing substructure!')
            for i,g in enumerate(group.flows_nosubs):
                gfp_keys.append(g.key)
                if g.key ==gf_key:
                    selected_gf = i
                    break
            try:
                # print(ozy_file,gf_key,gfp_keys,group.flows_nosubs)
                return group.flows_nosubs[selected_gf]
            except:
                print('Failed to recover old flow, so recomputing.')
                del f[str(gf.group.type)+'_data/flows_nosubs/'+str(gf.group._index)+'/'+str(gf_key)]
        else:
            for i,g in enumerate(group.flows):
                if g.key ==gf_key:
                    selected_gf = i
                    break
            return group.flows[selected_gf]
    else:
        if verbose: print('Writing flow data %s in %s_data'%(gf_key,group.type))
    f.close()

    # If substructre is removed, obtain regions
    if remove_subs:
        if verbose: print('Removing substructure!')
        subs = structure_regions(group, add_substructure=True, add_neighbours=False,
                            tidal_method='BT87_simple')
        nsubs = len(subs)
    else:
        nsubs = 0
    
    # Begin integration
    if verbose: print('Performing integration')
    if separate_phases:
        # Initialise Fortran derived type with attributes
        glob_attrs = amr_integrator.amr_region_attrs()
        glob_attrs.nvars = len(quantity_names)
        glob_attrs.nwvars = len(weight_names)
        glob_attrs.nfilter = len(filt)
        glob_attrs.nsubs = nsubs
        amr_integrator.allocate_amr_regions_attrs(glob_attrs)
        for i in range(0,len(quantity_names)):
            glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
        glob_attrs.result.nbins = pdf_bins
        glob_attrs.result.nvars = len(quantity_names)
        glob_attrs.result.nwvars = len(weight_names)
        glob_attrs.result.nfilter = len(filt)
        stats_utils.allocate_pdf(glob_attrs.result)
        for i in range(0, len(quantity_names)):
            mybins = get_code_bins(group.obj,'gas/'+quantity_names[i],pdf_bins)
            glob_attrs.result.varname.T.view('S128')[i] = quantity_names[i].ljust(128)
            glob_attrs.result.scaletype.T.view('S128')[i] = mybins[1].ljust(128)
            glob_attrs.result.bins[:,i] = mybins[0]
            glob_attrs.result.do_binning[i] = do_binning[i]
            glob_attrs.result.zero_index[i] = mybins[2]
            glob_attrs.result.linthresh[i] = mybins[3]
        for i in range(0, len(weight_names)):
            glob_attrs.result.wvarnames.T.view('S128')[i] = weight_names[i].ljust(128)

        for i in range(0, glob_attrs.nfilter):
            glob_attrs.filters[i] = filt[i]

        if remove_subs and nsubs>0:
            for i in range(0,nsubs):
                glob_attrs.subs[i] = subs[i]
        
        amr_integrator.integrate_region(output_path,reg,use_neigh,use_grav,glob_attrs)
        for i,f in enumerate(filt):
            phase_name = '_'+f.name.decode().split(' ')[0]
            # Assign results to galaxy object
            if group.obj.simulation.physics['hydro']:
                if verbose: print('Computing gas flow quantities')
                gf.data['mass_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.result.total[0,i,0,0], 'code_mass')
                gf.data['density_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,1,i)
                gf.data['temperature_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,2,i)
                # Powell et al. (2011) method
                # gf.data['massflow_rate_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[2,1,0]*4*np.pi*(shell_r**2), 'code_mass*code_velocity/code_length')
                # My method
                gf.data['massflow_rate_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.result.total[3,i,0,0]/shell_width, 'code_mass*code_velocity/code_length')
                gf.data['v_sphere_r_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,4,i)
                gf.data['v_tangential_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,5,i)
                gf.data['thermal_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.result.total[6,i,0,0], 'code_mass * code_velocity**2')
                gf.data['thermal_energy_specific_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,7,i)
                gf.data['grav_therpfrsphere_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,8,i)
                gf.data['grav_therpfrspherepos_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,9,i)
                gf.data['grav_therpfrsphereneg_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,10,i)
                if verbose: print('Massflow rate for %s phase  '%phase_name+str(gf.data['massflow_rate_'+d_key+'rvir'+phase_name].in_units('Msun/yr')))
                if verbose: print('Thermal support for %s phase  '%phase_name+str(gf.data['grav_therpfrsphere_'+d_key+'rvir'+phase_name]))
                if verbose: print('Radial velocity flow-weighted for %s phase  '%phase_name+str(gf.data['v_sphere_r_'+d_key+'rvir'+phase_name].to('km/s')))
                if group.obj.simulation.physics['metals']:
                    gf.data['metallicity_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,11,i)
            else:
                # TODO: Change shape for this empty array   
                gf.data['massflow_rate_'+d_key+'rvir'+phase_name] = group.obj.quantity(0.0, 'code_mass*code_velocity/code_length')
            
            if group.obj.simulation.physics['magnetic']:
                if verbose: print('Computing magnetic energies')
                if group.obj.simulation.physics['metals']:
                    gf.data['magnetic_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.result.total[12,i,0,0], 'code_mass * code_velocity**2')
                    gf.data['magnetic_energy_specific_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,13,i)
                    gf.data['grav_magpfrsphere_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,14,i)
                    gf.data['grav_magpfrspherepos_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,15,i)
                    gf.data['grav_magpfrsphereneg_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,16,i)
                    if verbose: print('Magnetic support for %s phase  '%phase_name+str(gf.data['grav_magpfrsphere_'+d_key+'rvir'+phase_name]))
                    if not group.obj.simulation.physics['cr']:
                        gf.data['grav_totpfrsphere_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,17,i)
                        gf.data['grav_totpfrspherepos_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,18,i)
                        gf.data['grav_totpfrsphereneg_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,19,i)
                        if verbose: print('Total support for %s phase  '%phase_name+str(gf.data['grav_totpfrsphere_'+d_key+'rvir'+phase_name]))
                else:
                    stat_array = pdf_handler_to_stats(group.obj,glob_attrs.result[12],i)
                    gf.data['magnetic_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.result.total[11,i,0,0], 'code_mass * code_velocity**2')
                    gf.data['magnetic_energy_specific_'+d_key+'rvir'+phase_name] = group.obj.array(stat_array, 'code_specific_energy')
                    gf.data['grav_magpfrsphere_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,13,i)
                    gf.data['grav_magpfrspherepos_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,14,i)
                    gf.data['grav_magpfrsphereneg_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,15,i)
                    if not group.obj.simulation.physics['cr']:
                        gf.data['grav_totpfrsphere_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,16,i)
                        gf.data['grav_totpfrspherepos_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,17,i)
                        gf.data['grav_totpfrsphereneg_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,18,i)

            if group.obj.simulation.physics['cr']:
                if verbose: print('Computing CR energies')
                if group.obj.simulation.physics['metals']:
                    gf.data['cr_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.result.total[17,i,0,0], 'code_mass * code_velocity**2')
                    gf.data['cr_energy_specific_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,18,i)
                    gf.data['grav_crpfrsphere_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,19,i)
                    gf.data['grav_crpfrspherepos_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,20,i)
                    gf.data['grav_crpfrsphereneg_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,21,i)
                    if verbose: print('CR support for %s phase  '%phase_name+str(gf.data['grav_crpfrsphere_'+d_key+'rvir'+phase_name]))
                    if verbose: print('CR support pos for %s phase  '%phase_name+str(gf.data['grav_crpfrspherepos_'+d_key+'rvir'+phase_name]))
                    if verbose: print('CR support neg for %s phase  '%phase_name+str(gf.data['grav_crpfrsphereneg_'+d_key+'rvir'+phase_name]))
                    gf.data['grav_totpfrsphere_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,22,i)
                    gf.data['grav_totpfrspherepos_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,23,i)
                    gf.data['grav_totpfrsphereneg_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,24,i)
                    if verbose: print('Total support for %s phase  '%phase_name+str(gf.data['grav_totpfrsphere_'+d_key+'rvir'+phase_name]))
                else:
                    gf.data['cr_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.result.total[16,i,0,0], 'code_mass * code_velocity**2')
                    gf.data['cr_energy_specific_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,17,i)
                    gf.data['grav_crpfrsphere_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,18,i)
                    gf.data['grav_crpfrspherepos_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,19,i)
                    gf.data['grav_crpfrsphereneg_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,20,i)
                    if verbose: print('CR support for %s phase  '%phase_name+str(gf.data['grav_crpfrsphere_'+d_key+'rvir'+phase_name]))
                    gf.data['grav_totpfrsphere_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,21,i)
                    gf.data['grav_totpfrsphere_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,22,i)
                    gf.data['grav_totpfrsphere_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,23,i)
                    if verbose: print('Total support for %s phase  '%phase_name+str(gf.data['grav_totpfrsphere_'+d_key+'rvir'+phase_name]))

            if group.obj.simulation.physics['rt']:
                if verbose: print('Computing ionisation fractions')
                gf.data['xHII_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,nvar,i)
                gf.data['xHeII_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,nvar+1,i)
                gf.data['xHeIII_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result,nvar+2,i)
    else:
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
        glob_attrs.nfilter = 1
        glob_attrs.nsubs = nsubs
        amr_integrator.allocate_amr_regions_attrs(glob_attrs)
        for i in range(0, len(quantity_names)):
            glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
            glob_attrs.result[i].nbins = pdf_bins
            glob_attrs.result[i].nfilter = 1
            glob_attrs.result[i].nwvars = len(weight_names)
            glob_attrs.result[i].varname = quantity_names[i]
            mybins = get_code_bins(group.obj,'gas/'+quantity_names[i],pdf_bins)
            glob_attrs.result[i].scaletype = mybins[1]
            stats_utils.allocate_pdf(glob_attrs.result[i])
            glob_attrs.result[i].bins = mybins[0]
            glob_attrs.result[i].do_binning = do_binning[i]
            for j in range(0, len(weight_names)):
                glob_attrs.result[i].wvarnames.T.view('S128')[j] = weight_names[j].ljust(128)

        if remove_subs and nsubs>0:
            for i in range(0,nsubs):
                glob_attrs.subs[i] = subs[i]

        amr_integrator.integrate_region(output_path,reg,filt,use_neigh,glob_attrs)

        # Assign results to galaxy object
        phase_name = ''
        i = 0
        if group.obj.simulation.physics['hydro']:
            if verbose: print('Computing gas flow quantities')
            stat_array = pdf_handler_to_stats(group.obj,glob_attrs.result[0],0)
            gf.data['density_'+d_key+'rvir'+phase_name] = group.obj.array(stat_array, 'code_density')
            stat_array = pdf_handler_to_stats(group.obj,glob_attrs.result[1],0)
            gf.data['temperature_'+d_key+'rvir'+phase_name] = group.obj.array(stat_array, 'code_temperature')
            # Powell et al. (2011) method
            # gf.data['massflow_rate_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[2,1,0]*4*np.pi*(shell_r**2), 'code_mass*code_velocity/code_length')
            # My method
            gf.data['massflow_rate_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.result[2].total[i,0,0]/shell_width, 'code_mass*code_velocity/code_length')
            stat_array = pdf_handler_to_stats(group.obj,glob_attrs.result[3],0)
            gf.data['v_sphere_r_'+d_key+'rvir'+phase_name] = group.obj.array(stat_array, 'code_velocity')
            gf.data['thermal_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.result[4].total[i,0,0], 'code_mass * code_velocity**2')
            stat_array = pdf_handler_to_stats(group.obj,glob_attrs.result[5],0)
            gf.data['thermal_energy_specific_'+d_key+'rvir'+phase_name] = group.obj.array(stat_array, 'code_specific_energy')
            stat_array = pdf_handler_to_stats(group.obj,glob_attrs.result[6],0)
            gf.data['grav_therpfrsphere_'+d_key+'rvir'+phase_name] = group.obj.array(stat_array, 'dimensionless')
            if verbose: print('Massflow rate for %s phase  '%phase_name+str(gf.data['massflow_rate_'+d_key+'rvir'+phase_name].in_units('Msun/yr')))
            if verbose: print('Thermal support for %s phase  '%phase_name+str(gf.data['grav_therpfrsphere_'+d_key+'rvir'+phase_name]))
            if verbose: print('Radial velocity flow-weighted for %s phase  '%phase_name+str(gf.data['v_sphere_r_'+d_key+'rvir'+phase_name].to('km/s')))
            if group.obj.simulation.physics['metals']:
                stat_array = pdf_handler_to_stats(group.obj,glob_attrs.result[7],)
                gf.data['metallicity_'+d_key+'rvir'+phase_name] = stat_array
        else:
            gf.data['massflow_rate_'+d_key+'rvir'+phase_name] = group.obj.quantity(0.0, 'code_mass*code_velocity/code_length')
        
        if group.obj.simulation.physics['magnetic']:
            if verbose: print('Computing magnetic energies')
            if group.obj.simulation.physics['metals']:
                stat_array = pdf_handler_to_stats(group.obj,glob_attrs.result[9],i)
                gf.data['magnetic_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.result[8].total[i,0,0], 'code_mass * code_velocity**2')
                gf.data['magnetic_energy_specific_'+d_key+'rvir'+phase_name] = group.obj.array(stat_array, 'code_specific_energy')
            else:
                stat_array = pdf_handler_to_stats(group.obj,glob_attrs.result[8],i)
                gf.data['magnetic_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.result[7].total[i,0,0], 'code_mass * code_velocity**2')
                gf.data['magnetic_energy_specific_'+d_key+'rvir'+phase_name] = group.obj.array(stat_array, 'code_specific_energy')

        if group.obj.simulation.physics['cr']:
            if verbose: print('Computing CR energies')
            if group.obj.simulation.physics['metals']:
                gf.data['cr_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.result[10].total[i,0,0], 'code_mass * code_velocity**2')
                stat_array = pdf_handler_to_stats(group.obj,glob_attrs.result[11],i)
                gf.data['cr_energy_specific_'+d_key+'rvir'+phase_name] = group.obj.array(stat_array, 'code_specific_energy')
                stat_array = pdf_handler_to_stats(group.obj,glob_attrs.result[12],i)
                gf.data['grav_crpfrsphere_'+d_key+'rvir'+phase_name] = group.obj.array(stat_array, 'dimensionless')
                if verbose: print('CR support for %s phase  '%phase_name+str(gf.data['grav_crpfrsphere_'+d_key+'rvir'+phase_name]))
            else:
                gf.data['cr_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.result[9].total[i,0,0], 'code_mass * code_velocity**2')
                stat_array = pdf_handler_to_stats(group.obj,glob_attrs.result[10],i)
                gf.data['cr_energy_specific_'+d_key+'rvir'+phase_name] = group.obj.array(stat_array, 'code_specific_energy')
                stat_array = pdf_handler_to_stats(group.obj,glob_attrs.result[11],i)
                gf.data['grav_crpfrsphere_'+d_key+'rvir'+phase_name] = group.obj.array(stat_array, 'dimensionless')
                if verbose: print('CR support for %s phase  '%phase_name+str(gf.data['grav_crpfrsphere_'+d_key+'rvir'+phase_name]))

        if group.obj.simulation.physics['rt']:
            if verbose: print('Computing ionisation fractions')
            gf.data['xHII_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result[nvar],i)
            gf.data['xHeII_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result[nvar+1],i)
            gf.data['xHeIII_'+d_key+'rvir'+phase_name] = pdf_handler_to_stats(group.obj,glob_attrs.result[nvar+2],i)

    if save:
        write_flow(group.obj, ozy_file, gf, d_key, width_key, verbose=verbose)
    
    return gf
