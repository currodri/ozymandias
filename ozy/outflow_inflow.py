import numpy as np
import h5py
import os
import ozy
from unyt import unyt_array,unyt_quantity
from ozy.dict_variables import common_variables,grid_variables,particle_variables,get_code_units
# TODO: Allow for parallel computation of flows.
from joblib import Parallel, delayed
import sys
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/part')
from amr2 import amr_integrator
from ozy.utils import init_region,init_filter
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

def get_flow_name(gf_group, gf, r):
    """Create an individual galactic flow identifier name."""
    name = str(gf.type)
    name += '|'+str(gf.region['type'])
    name += '|'+str(r)
    if name in gf_group:
        name += '|new'
    return name

def check_if_same_flow(hd,gf,r):
    """This function looks at a OZY file for the flow data specified by gf for a particular object."""
    if not str(gf.group.type)+'_data/flows/'+str(gf.group._index) in hd:
        return False, 'none'
    for g in hd[str(gf.group.type)+'_data/flows/'+str(gf.group._index)].keys():
        check_type = (g.split('|')[0] == gf.type)
        check_region = (g.split('|')[1] == gf.region['type'])
        check_r = (g.split('|')[2] == str(r))
        if check_type and check_region and check_r:
            return True, g
    return False, 'none'

def write_flow(obj,ozy_file,gf,r):
    """This function writes the resulting flow data to the original OZY file."""
    f = h5py.File(ozy_file, 'r+')

    # Create group in HDF5 file
    try:
        flows = f.create_group(str(gf.group.type)+'_data/flows/'+str(gf.group._index))
    except:
        flows = f[str(gf.group.type)+'_data/flows/'+str(gf.group._index)]

    # Clean data and save to dataset
    gf_name = get_flow_name(flows,gf,r)
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

def compute_flows(group,ozy_file,flow_type,rmin=(0.0,'rvir'), rmax=(1.0,'rvir'),recompute=False,save=False,separate_phases=True):
    """Function which computes the analysis of galaxy wide flows, including outflows,
        inflows or AGN feedback."""

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
        raise ValueError("This type of galactic flow is not supported. Please check!")
    
    d_key = str(int(100*0.5*(rmax[0]+rmin[0])))
    shell_width = reg.rmax-reg.rmin
    shell_r = 0.5*(reg.rmax+reg.rmin)

    # Since the number of global quantities that can be computed
    # from the gas data depends on the specific configuration
    # of a simulation, this needs to be determined at the start
    # See: ozy/sim_attributes.py/assign_attributes

    output_path = group.obj.simulation.fullpath
    nvar = 0
    quantity_names = []
    weight_names = ['cumulative','volume','massflux_rate_sphere_r']

    if group.obj.simulation.physics['hydro']:
        quantity_names += ['density','temperature',
                            'momentum_sphere_r','v_sphere_r',
                            'thermal_energy','thermal_energy_specific']
    if group.obj.simulation.physics['metals']:
        quantity_names += ['metallicity']
    if group.obj.simulation.physics['magnetic']:
        quantity_names += ['magnetic_energy','magnetic_energy_specific']
    if group.obj.simulation.physics['cr']:
        quantity_names += ['cr_energy','cr_energy_specific']

    nvar = len(quantity_names)
    if group.obj.simulation.physics['rt']:
        quantity_names += ['xHII','xHeII','xHeIII']

    # Initialise GalacticFlow object
    gf = GalacticFlow(group)
    gf.type = flow_type
    gf.weightvars = weight_names

    # Save region and filter details to GalacticFlow object
    gf._get_python_region(reg)
    if separate_phases:
        gf._get_python_filter(filt[0])
    else:
        gf._get_python_filter(filt)

    # Check if flow data is already present and if it coincides with the new one
    f = h5py.File(ozy_file, 'r+')
    gf_present, gf_key = check_if_same_flow(f,gf,d_key)
    if gf_present and recompute:
        del f[str(gf.group.type)+'_data/flows/'+str(gf.group._index)+'/'+str(gf_key)]
        print('Overwriting flow data in %s_data'%group.type)
    elif gf_present and not recompute:
        print('Flow data with same details already present for galaxy %s. No overwritting!'%group._index)
        group._init_flows()
        for i,g in enumerate(group.flows):
            if g.key ==gf_key:
                selected_gf = i
                break
        return group.flows[selected_gf]
    else:
        print('Writing flow data in %s_data'%group.type)
    f.close()

    # Begin integration
    print('Performing integration')
    if separate_phases:
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
        glob_attrs.nfilter = len(filt)
        amr_integrator.allocate_amr_regions_attrs(glob_attrs)
        for i in range(0, len(quantity_names)):
            glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
        for i in range(0, len(weight_names)):
            glob_attrs.wvarnames.T.view('S128')[i] = weight_names[i].ljust(128)
        for i in range(0, glob_attrs.nfilter):
            glob_attrs.filters[i] = filt[i]
        
        amr_integrator.integrate_region(output_path,reg,glob_attrs)
        for i,f in enumerate(filt):
            if f.name == 'all':
                phase_name = ''
            else:
                phase_name = '_'+f.name.decode().split(' ')[0]
            
            # Assign results to galaxy object
            if group.obj.simulation.physics['hydro']:
                print('Computing gas flow quantities')
                gf.data['density_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,0,2,0], 'code_density')
                gf.data['temperature_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,1,2,0], 'code_temperature')
                # Powell et al. (2011) method
                # gf.data['massflow_rate_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[2,1,0]*4*np.pi*(shell_r**2), 'code_mass*code_velocity/code_length')
                # My method
                gf.data['massflow_rate_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,2,0,0]/shell_width, 'code_mass*code_velocity/code_length')
                gf.data['v_sphere_r_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,3,2,0], 'code_velocity')
                gf.data['thermal_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,4,0,0], 'code_mass * code_velocity**2')
                gf.data['thermal_energy_specific_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,5,2,0], 'code_specific_energy')
                print('Massflow rate for %s phase'%phase_name+str(gf.data['massflow_rate_'+d_key+'rvir'+phase_name].in_units('Msun/yr')))
                if group.obj.simulation.physics['metals']:
                    gf.data['metallicity_'+d_key+'rvir'+phase_name] = glob_attrs.data[i,6,2,0]
            else:
                gf.data['massflow_rate_'+d_key+'rvir'+phase_name] = group.obj.quantity(0.0, 'code_mass*code_velocity/code_length')
            
            if group.obj.simulation.physics['magnetic']:
                print('Computing magnetic energies')
                if group.obj.simulation.physics['metals']:
                    gf.data['magnetic_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,7,0,0], 'code_mass * code_velocity**2')
                    gf.data['magnetic_energy_specific_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,8,2,0], 'code_specific_energy')
                else:
                    gf.data['magnetic_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,6,0,0], 'code_mass * code_velocity**2')
                    gf.data['magnetic_energy_specific_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,7,2,0], 'code_specific_energy')

            if group.obj.simulation.physics['cr']:
                print('Computing CR energies')
                if group.obj.simulation.physics['metals']:
                    gf.data['cr_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,9,0,0], 'code_mass * code_velocity**2')
                    gf.data['cr_energy_specific_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,10,2,0], 'code_specific_energy')
                else:
                    gf.data['cr_energy_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,8,0,0], 'code_mass * code_velocity**2')
                    gf.data['cr_energy_specific_'+d_key+'rvir'+phase_name] = group.obj.quantity(glob_attrs.data[i,9,2,0], 'code_specific_energy')

            if group.obj.simulation.physics['rt']:
                print('Computing ionisation fractions')
                gf.data['xHII_'+d_key+'rvir'+phase_name] = glob_attrs.data[i,nvar,2,0]
                gf.data['xHeII_'+d_key+'rvir'+phase_name] = glob_attrs.data[i,nvar+1,2,0]
                gf.data['xHeIII_'+d_key+'rvir'+phase_name] = glob_attrs.data[i,nvar+2,2,0]
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
        amr_integrator.allocate_amr_regions_attrs(glob_attrs)
        for i in range(0, len(quantity_names)):
            glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
        for i in range(0, len(weight_names)):
            glob_attrs.wvarnames.T.view('S128')[i] = weight_names[i].ljust(128)

        amr_integrator.integrate_region(output_path,reg,filt,glob_attrs)

        # Assign results to galaxy object
        if group.obj.simulation.physics['hydro']:
            print('Computing gas flow quantities')
            gf.data['density_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,0,2,0], 'code_density')
            gf.data['temperature_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,1,2,0], 'code_temperature')
            gf.data['massflow_rate_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,2,1,0]*4*np.pi*(shell_r**2), 'code_mass*code_velocity/code_length')
            gf.data['v_sphere_r_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,3,2,0], 'code_velocity')
            gf.data['thermal_energy_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,4,0,0], 'code_mass * code_velocity**2')
            gf.data['thermal_energy_specific_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,5,2,0], 'code_specific_energy')
            print('Massflow rate '+str(gf.data['massflow_rate_'+d_key+'rvir'].in_units('Msun/yr')))
            if group.obj.simulation.physics['metals']:
                gf.data['metallicity_'+d_key+'rvir'] = glob_attrs.data[0,6,2,0]
        else:
            gf.data['massflow_rate_'+d_key+'rvir'] = group.obj.quantity(0.0, 'code_mass*code_velocity/code_length')
        
        if group.obj.simulation.physics['magnetic']:
            print('Computing magnetic energies')
            if group.obj.simulation.physics['metals']:
                gf.data['magnetic_energy_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,7,0,0], 'code_mass * code_velocity**2')
                gf.data['magnetic_energy_specific_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,8,2,0], 'code_specific_energy')
            else:
                gf.data['magnetic_energy_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,6,0,0], 'code_mass * code_velocity**2')
                gf.data['magnetic_energy_specific_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,7,2,0], 'code_specific_energy')

        if group.obj.simulation.physics['cr']:
            print('Computing CR energies')
            if group.obj.simulation.physics['metals']:
                gf.data['cr_energy_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,9,0,0], 'code_mass * code_velocity**2')
                gf.data['cr_energy_specific_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,10,2,0], 'code_specific_energy')
            else:
                gf.data['cr_energy_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,8,0,0], 'code_mass * code_velocity**2')
                gf.data['cr_energy_specific_'+d_key+'rvir'] = group.obj.quantity(glob_attrs.data[0,9,2,0], 'code_specific_energy')

        if group.obj.simulation.physics['rt']:
            print('Computing ionisation fractions')
            gf.data['xHII_'+d_key+'rvir'] = glob_attrs.data[0,nvar,2,0]
            gf.data['xHeII_'+d_key+'rvir'] = glob_attrs.data[0,nvar+1,2,0]
            gf.data['xHeIII_'+d_key+'rvir'] = glob_attrs.data[0,nvar+2,2,0]

    if save:
        write_flow(group.obj, ozy_file, gf, d_key)
    
    return gf
