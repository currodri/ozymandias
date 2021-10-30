import numpy as np
import h5py
import os
from yt import YTArray,YTQuantity
import ozy
from ozy.dict_variables import common_variables,grid_variables,particle_variables,get_code_units
# TODO: Allow for parallel computation of flows.
from joblib import Parallel, delayed
import sys
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/part')
from amr2 import io_ramses
from amr2 import vectors
from amr2 import geometrical_regions as geo
from amr2 import filtering
from amr2 import amr_integrator
from ozy.saver import _write_attrib
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
        
        from yt import YTArray
        for k,v in self.__dict__.items():
            if k in blacklist:
                continue
            if isinstance(v, (YTArray, YTQuantity)):
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
        """Save the Fortran derived type as a dictionary inside the GalacticFlow class (only the necessary info)."""
        self.region = {}
        self.region['type'] = reg.name.decode().split(' ')[0]
        self.region['centre'] = YTArray([reg.centre.x, reg.centre.y, reg.centre.z], 'code_length', registry=self.obj.unit_registry)
        self.region['axis'] = YTArray([reg.axis.x, reg.axis.y, reg.axis.z], 'dimensionless', registry=self.obj.unit_registry)
        self.region['rmin'] = YTQuantity(reg.rmin, 'code_length', registry=self.obj.unit_registry)
        self.region['rmax'] = YTQuantity(reg.rmax, 'code_length', registry=self.obj.unit_registry)
        self.region['r'] = YTQuantity(0.5*(reg.rmax+reg.rmin), 'code_length', registry=self.obj.unit_registry)
    
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
                cond_value = YTQuantity(filt.cond_vals[i], str(cond_units), registry=self.obj.unit_registry)
                cond_str = cond_var+'/'+cond_op+'/'+str(cond_value.d)+'/'+cond_units
                self.filter['conditions'].append(cond_str)

def init_region(group, region_type, rmin=0.0, rmax = 1.0):
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
        try:
            velocity = group.velocity.in_units('code_velocity')
        except:
            velocity = YTArray(group.velocity,'km/s',registry=group.obj.unit_registry).in_units('code_velocity')
        bulk.x, bulk.y, bulk.z = velocity[0].d, velocity[1].d, velocity[2].d
        reg.bulk_velocity = bulk
        reg.rmin = rmin*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        reg.rmax = rmax*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
    else:
        raise KeyError('Region type not supported. Please check!')
    return reg
    
def init_filter(group,cond_strs, name):
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
    if not str(gf.group.obj_type)+'_data/flows/'+str(gf.group._index) in hd:
        return False, 'none'
    for g in hd[str(gf.group.obj_type)+'_data/flows/'+str(gf.group._index)].keys():
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
        flows = f.create_group(str(gf.group.obj_type)+'_data/flows/'+str(gf.group._index))
    except:
        flows = f[str(gf.group.obj_type)+'_data/flows/'+str(gf.group._index)]

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
        print(name,get_code_units(name))
        clean_data[var].attrs.create('units', get_code_units(name))
    
    f.close()
    return

def compute_flows(group,ozy_file,flow_type,rmin=0.0,rmax=1.0,recompute=False,save=False,separate_phases=True):
    """Function which computes the analysis of galaxy wide flows, including outflows,
        inflows or AGN feedback."""

    # Begin by setting up the necessary details of a given flow type
    if flow_type == 'outflow':
        reg = init_region(group, 'sphere', rmin=rmin, rmax=rmax)
        if separate_phases:
            all = init_filter(group,cond_strs='v_sphere_r/>/0/km*s**-1',name='all')
            # Powell et al. (2011)
            # clumpy = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','density/>/1e-23/g*cm**-3'],name='outflow_clumpy')
            # filament = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','density/>/1e-25/g*cm**-3','density/</1e-23/g*cm**-3','temperature/</2e4/K'],name='outflow_filament')
            # cold_diffuse = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','density/</1e-25/g*cm**-3','temperature/</2e4/K'],name='outflow_cold_diffuse')
            # cold_dense = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','density/>/1e-25/g*cm**-3','temperature/</2e4/K'],name='outflow_cold_dense')
            # warm = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','density/</1e-23/g*cm**-3','temperature/>/2e4/K','temperature/</2e5/K'],name='outflow_warm')
            # hot = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','density/</1e-23/g*cm**-3','temperature/>/2e5/K'],name='outflow_hot')

            # Sergio's suggestion
            hot = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','temperature/>/1e5/K'],name='hot')
            warm_ionised = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','temperature/</1e5/K','temperature/>/9e3/K'],name='warm_ionised')
            warm_neutral = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','temperature/</9e3/K','temperature/>/1e3/K'],name='warm_neutral')
            cold = init_filter(group,cond_strs=['v_sphere_r/>/0/km*s**-1','temperature/</1e3/K'],name='cold')
            filt = [all,hot,warm_ionised,warm_neutral,cold]
        else:
            filt = init_filter(group,cond_strs='v_sphere_r/>/0/km*s**-1',name='outflow')
    elif flow_type == 'inflow':
        reg = init_region(group, 'sphere', rmin=rmin, rmax=rmax)
        if separate_phases:
            all = init_filter(group,cond_strs='v_sphere_r/<=/0/km*s**-1',name='all')
            # Powell et al. (2011)
            # clumpy = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/>/1e-23/g*cm**-3'],name='outflow_clumpy')
            # filament = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/>/1e-25/g*cm**-3','density/</1e-23/g*cm**-3','temperature/</2e4/K'],name='outflow_filament')
            # cold_diffuse = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/</1e-25/g*cm**-3','temperature/</2e4/K'],name='outflow_cold_diffuse')
            # cold_dense = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/>/1e-25/g*cm**-3','temperature/</2e4/K'],name='outflow_cold_dense')
            # warm = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/</1e-23/g*cm**-3','temperature/>/2e4/K','temperature/</2e5/K'],name='outflow_warm')
            # hot = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/</1e-23/g*cm**-3','temperature/>/2e5/K'],name='outflow_hot')

            cold = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/>/1e-24/g*cm**-3','temperature/</500/K'],name='cold')
            warm = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/</1e-24/g*cm**-3','density/>/1e-26/g*cm**-3','temperature/>/500/K','temperature/</5e5/K'],name='warm')
            hot = init_filter(group,cond_strs=['v_sphere_r/<=/0/km*s**-1','density/</1e-26/g*cm**-3','temperature/>/5e5/K'],name='hot')
            filt = [all,cold,warm,hot]
        else:
            filt = init_filter(group,cond_strs='v_sphere_r/<=/0/km*s**-1',name='inflow')
    else:
        raise ValueError("This type of galactic flow is not supported. Please check!")
    
    d_key = str(int(100*0.5*(rmax+rmin)))
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
        del f[str(gf.group.obj_type)+'_data/flows/'+str(gf.group._index)+'/'+str(gf_key)]
        print('Overwriting flow data in %s_data'%group.obj_type)
    elif gf_present and not recompute:
        print('Flow data with same details already present for galaxy %s. No overwritting!'%group._index)
        group._init_flows()
        for i,g in enumerate(group.flows):
            if g.key ==gf_key:
                selected_gf = i
                break
        return group.flows[selected_gf]
    else:
        print('Writing flow data in %s_data'%group.obj_type)
    f.close()

    # Begin integration
    print('Performing integration')
    if separate_phases:
        for f in filt:
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
            amr_integrator.allocate_amr_regions_attrs(glob_attrs)
            for i in range(0, len(quantity_names)):
                glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
            for i in range(0, len(weight_names)):
                glob_attrs.wvarnames.T.view('S128')[i] = weight_names[i].ljust(128)
            
            amr_integrator.integrate_region(output_path,reg,f,glob_attrs)
            if f.name == 'all':
                phase_name = ''
            else:
                phase_name = f.name.decode().split(' ')[0]
            
            # Assign results to galaxy object
            if group.obj.simulation.physics['hydro']:
                print('Computing gas flow quantities')
                gf.data['density_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[0,2,0], 'code_density', registry=group.obj.unit_registry)
                gf.data['temperature_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[1,2,0], 'code_temperature', registry=group.obj.unit_registry)
                # Powell et al. (2011) method
                # gf.data['massflow_rate_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[2,1,0]*4*np.pi*(shell_r**2), 'code_mass*code_velocity/code_length', registry=group.obj.unit_registry)
                # My method
                gf.data['massflow_rate_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[2,0,0]/shell_width, 'code_mass*code_velocity/code_length', registry=group.obj.unit_registry)
                gf.data['v_sphere_r_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[3,2,0], 'code_velocity', registry=group.obj.unit_registry)
                gf.data['thermal_energy_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[4,0,0], 'code_mass * code_velocity**2', registry=group.obj.unit_registry)
                gf.data['thermal_energy_specific_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[5,2,0], 'code_specific_energy', registry=group.obj.unit_registry)
                print('Massflow rate '+str(gf.data['massflow_rate_'+d_key+'rvir_'+phase_name].in_units('Msun/yr')))
                if group.obj.simulation.physics['metals']:
                    gf.data['metallicity_'+d_key+'rvir_'+phase_name] = glob_attrs.data[6,2,0]
            else:
                gf.data['massflow_rate_'+d_key+'rvir_'+phase_name] = YTQuantity(0.0, 'code_mass*code_velocity/code_length', registry=group.obj.unit_registry)
            
            if group.obj.simulation.physics['magnetic']:
                print('Computing magnetic energies')
                if group.obj.simulation.physics['metals']:
                    gf.data['magnetic_energy_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[7,0,0], 'code_mass * code_velocity**2', registry=group.obj.unit_registry)
                    gf.data['magnetic_energy_specific_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[8,2,0], 'code_specific_energy', registry=group.obj.unit_registry)
                else:
                    gf.data['magnetic_energy_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[6,0,0], 'code_mass * code_velocity**2', registry=group.obj.unit_registry)
                    gf.data['magnetic_energy_specific_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[7,2,0], 'code_specific_energy', registry=group.obj.unit_registry)

            if group.obj.simulation.physics['cr']:
                print('Computing CR energies')
                if group.obj.simulation.physics['metals']:
                    gf.data['cr_energy_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[9,0,0], 'code_mass * code_velocity**2', registry=group.obj.unit_registry)
                    gf.data['cr_energy_specific_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[10,2,0], 'code_specific_energy', registry=group.obj.unit_registry)
                else:
                    gf.data['cr_energy_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[8,0,0], 'code_mass * code_velocity**2', registry=group.obj.unit_registry)
                    gf.data['cr_energy_specific_'+d_key+'rvir_'+phase_name] = YTQuantity(glob_attrs.data[9,2,0], 'code_specific_energy', registry=group.obj.unit_registry)

            if group.obj.simulation.physics['rt']:
                print('Computing ionisation fractions')
                gf.data['xHII_'+d_key+'rvir_'+phase_name] = glob_attrs.data[nvar,2,0]
                gf.data['xHeII_'+d_key+'rvir_'+phase_name] = glob_attrs.data[nvar+1,2,0]
                gf.data['xHeIII_'+d_key+'rvir_'+phase_name] = glob_attrs.data[nvar+2,2,0]
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
        amr_integrator.allocate_amr_regions_attrs(glob_attrs)
        for i in range(0, len(quantity_names)):
            glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
        for i in range(0, len(weight_names)):
            glob_attrs.wvarnames.T.view('S128')[i] = weight_names[i].ljust(128)

        amr_integrator.integrate_region(output_path,reg,filt,glob_attrs)

        # Assign results to galaxy object
        if group.obj.simulation.physics['hydro']:
            print('Computing gas flow quantities')
            gf.data['density_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[0,2,0], 'code_density', registry=group.obj.unit_registry)
            gf.data['temperature_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[1,2,0], 'code_temperature', registry=group.obj.unit_registry)
            gf.data['massflow_rate_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[2,1,0]*4*np.pi*(shell_r**2), 'code_mass*code_velocity/code_length', registry=group.obj.unit_registry)
            gf.data['v_sphere_r_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[3,2,0], 'code_velocity', registry=group.obj.unit_registry)
            gf.data['thermal_energy_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[4,0,0], 'code_mass * code_velocity**2', registry=group.obj.unit_registry)
            gf.data['thermal_energy_specific_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[5,2,0], 'code_specific_energy', registry=group.obj.unit_registry)
            print('Massflow rate '+str(gf.data['massflow_rate_'+d_key+'rvir'].in_units('Msun/yr')))
            if group.obj.simulation.physics['metals']:
                gf.data['metallicity_'+d_key+'rvir'] = glob_attrs.data[6,2,0]
        else:
            gf.data['massflow_rate_'+d_key+'rvir'] = YTQuantity(0.0, 'code_mass*code_velocity/code_length', registry=group.obj.unit_registry)
        
        if group.obj.simulation.physics['magnetic']:
            print('Computing magnetic energies')
            if group.obj.simulation.physics['metals']:
                gf.data['magnetic_energy_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[7,0,0], 'code_mass * code_velocity**2', registry=group.obj.unit_registry)
                gf.data['magnetic_energy_specific_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[8,2,0], 'code_specific_energy', registry=group.obj.unit_registry)
            else:
                gf.data['magnetic_energy_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[6,0,0], 'code_mass * code_velocity**2', registry=group.obj.unit_registry)
                gf.data['magnetic_energy_specific_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[7,2,0], 'code_specific_energy', registry=group.obj.unit_registry)

        if group.obj.simulation.physics['cr']:
            print('Computing CR energies')
            if group.obj.simulation.physics['metals']:
                gf.data['cr_energy_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[9,0,0], 'code_mass * code_velocity**2', registry=group.obj.unit_registry)
                gf.data['cr_energy_specific_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[10,2,0], 'code_specific_energy', registry=group.obj.unit_registry)
            else:
                gf.data['cr_energy_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[8,0,0], 'code_mass * code_velocity**2', registry=group.obj.unit_registry)
                gf.data['cr_energy_specific_'+d_key+'rvir'] = YTQuantity(glob_attrs.data[9,2,0], 'code_specific_energy', registry=group.obj.unit_registry)

        if group.obj.simulation.physics['rt']:
            print('Computing ionisation fractions')
            gf.data['xHII_'+d_key+'rvir'] = glob_attrs.data[nvar,2,0]
            gf.data['xHeII_'+d_key+'rvir'] = glob_attrs.data[nvar+1,2,0]
            gf.data['xHeIII_'+d_key+'rvir'] = glob_attrs.data[nvar+2,2,0]

    if save:
        write_flow(group.obj, ozy_file, gf, d_key)
    
    return gf
