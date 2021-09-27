import numpy as np
import h5py
import os
from yt import YTArray,YTQuantity
import ozy
from ozy.dict_variables import common_variables,grid_variables,particle_variables,get_code_units
# TODO: Allow for parallel computation of phase diagrams.
from joblib import Parallel, delayed
import sys
sys.path.append('/mnt/extraspace/currodri/Codes/ozymandias/ozy/amr')
sys.path.append('/mnt/extraspace/currodri/Codes/ozymandias/ozy/part')
from amr2 import vectors
from amr2 import geometrical_regions as geo
from amr2 import filtering
from amr2 import amr_integrator
from ozy.saver import _write_attrib
blacklist = []

class GalacticFlow(object):

    def __init__(self,group):
        self.group = group
        self.obj = group.obj
        self.type = None
        self.region = None
        self.filter = None
        self.data = {}
        self.weightvars = {}

    def _serialise(self,hdd):
        """This makes possible to save the group galactic flow attrs as dataset
            attributes of an HDF5 group.""" 
        
        from yt import YTArray
        for k,v in self.__dict__.items():
            if k in blacklist:
                continue
            if isinstance(v, YTArray, YTQuantity):
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

    if region_type == 'shell':
        reg.name = 'shell'
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
        velocity = YTArray(group.velocity,'km/s',registry=group.obj.unit_registry).in_units('code_velocity')
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

def compute_flows(group,ozy_file,flow_type,rmin=0.0,rmax=1.0,recompute=False,save=False):
    """Function which computes the analysis of galaxy wide flows, including outflows,
        inflows or AGN feedback."""

    # Begin by setting up the necessary details of a given flow type
    if flow_type == 'outflow':
        reg = init_region(self, 'shell', rmin=rmin, rmax=rmax)
        filt = init_filter(self,cond_strs='v_sphere_r/>/0/km*s**-1','outflow')
    elif flow_type == 'inflow':
        reg = init_region(self, 'shell', rmin=rmin, rmax=rmax)
        filt = init_filter(self,cond_strs='v_sphere_r/<=/0/km*s**-1','inflow')
    else:
        raise ValueError("This type of galactic flow is not supported. Please check!")

    # Since the number of global quantities that can be computed
    # from the gas data depends on the specific configuration
    # of a simulation, this needs to be determined at the start
    # See: ozy/sim_attributes.py/assign_attributes

    nvar = 0
    quantity_names = []
    weight_names = ['cumulative','massflux_rate_sphere_r']

    if self.obj.simulation.physics['hydro']:
        quantity_names += ['density','temperature',
                            'momentum_sphere_r','v_sphere_r',
                            'thermal_energy','thermal_energy_specific']
    if self.obj.simulation.physics['metals']:
        quantity_names += ['metallicity']
    if self.obj.simulation.physics['magnetic']:
        quantity_names += ['magnetic_energy','magnetic_energy_specific']
    if self.obj.simulation.physics['cr']:
        quantity_names += ['cr_energy','cr_energy_specific']

    nvar = len(quantity_names)
    if self.obj.simulation.physics['rt']:
        quantity_names += ['xHII','xHeII','xHeIII']

    # Initialise GalacticFlow object
    gf = GalacticFlow(group)
    gf.type = 