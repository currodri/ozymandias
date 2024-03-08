import functools
import os.path
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pprint import pprint

import h5py
import numpy as np
from unyt import UnitRegistry,unyt_array,unyt_quantity

from ozy.sim_attributes import SimulationAttributes
from ozy.utils import info_printer
from ozy.saver import _write_attrib, _write_dict

class LazyDataset:
    """A lazily-loaded HDF5 dataset.
    """
    def __init__(self, obj, dataset_path):
        self._obj = obj
        self._dataset_path = dataset_path
        self._data = None
    def __getitem__(self, index):
        if self._data is None:
            with h5py.File(self._obj.data_file, 'r') as hd:
                dataset = hd[self._dataset_path]
                if 'unit' in dataset.attrs:
                    self._data = unyt_array(dataset[:],
                                         dataset.attrs['unit'],
                                         registry=self._obj.unit_registry)
                else:
                    self._data = dataset[:]
        return self._data.__getitem__(index)
    
    def __setitem__(self, index, value):
        if self._data is not None:
            self._data.__setitem__(index, value)
        else:
            raise ValueError("Cannot perform item assignment before loading the dataset.")

class LazyList(Sequence):
    """This type should be indistinguishable from the built-in list.
    Any observable difference except the explicit type and performance
    is considered a bug.
    The implementation wraps a list which is initially filled with None,
    which is very fast to create at any size because None is a singleton.
    The initial elements are replaced by calling the passed-in callable
    as they are accessed.
    """
    def __init__(self, length, builder):
        self._inner = [None] * length
        self._builder = builder

    def __contains__(self, value):
        for i in range(len(self)):
            if self[i] == value:
                return True
        return False

    def __getitem__(self, index):
        trial_output = self._inner[index]

        # Handle uninitialized elements for integer indices
        if trial_output is None:
            self._inner[index] = self._builder(index)

        # And for all other kinds of indices
        if isinstance(trial_output, list) and None in trial_output:
            for i in range(len(self))[index]:
                if self._inner[i] is None:
                    self._inner[i] = self._builder(i)

        return self._inner[index]

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __reversed__(self):
        LazyList(len(self), lambda i: self._builder(len(self) - i))

    def count(self, value):
        return sum(i == value for i in self)

    def index(self, value, start=0, stop=None):
        if stop is None:
            stop = len(self)
        for i in range(start, stop):
            if self[i] == value:
                return i
        raise ValueError('{} is not in LazyList'.format(value))

    def __len__(self):
        return len(self._inner)

class LazyDict(Mapping):
    """This type should be indistinguishable from the built-in dict.
    Any observable difference except the explicit type and performance
    is considered a bug.
    The implementation wraps a dict which initially maps every key to None,
    and are replaced by calling the passed-in callable as they are accessed.
    """
    def __init__(self, keys, builder):
        self._inner = {k: None for k in keys}
        self._builder = builder

    def _init_all(self):
        """Internal use only, for operations that need all values"""
        for k in self._inner:
            self[k]

    def __contains__(self, key):
        return key in self._inner

    def __eq__(self, other):
        self._init_all()
        return self._inner == other

    def __getitem__(self, key):
        value = self._inner[key]
        if value is None:
            value = self._builder(key)
            self._inner[key] = value
        return value
    
    def __setitem__(self, key, value):
        self._inner[key] = value

    def get(self, key, default=None):
        if key in self._inner:
            return self[key]
        return default

    def items(self):
        for k in self._inner:
            yield (k, self[k])

    def keys(self):
        return self._inner.keys()

    def values(self):
        for k in self._inner:
            yield self[k]

    def __len__(self):
        return len(self._inner)

    def __iter__(self):
        return iter(self._inner)

    def __str__(self):
        return str(dict(self))

    def __repr__(self):
        return repr(dict(self))

    def _repr_pretty_(self, p, cycle):
        p.pretty(dict(self))

class Profile:
    def __init__(self, obj, index, group_type, key, hd, nosubs):
        self.obj = obj
        self._index = index
        self.group_type = group_type
        self.key = key
        self.nbins = 0
        self.xvar = None
        self.region = None
        self.filter = {}
        self.lmax = 0
        self.yvars = {}
        self.weightvars = {}
        self.xdata = None
        self.ydata = {}
        self.nosubs = nosubs

        self.blacklist = [
            'obj','_index','group_type','region',
            'filter','yvars','weightvars','xdata','ydata'
        ]
        self._unpack(hd)
    
    def _unpack(self, hd):
        if self.nosubs:
            path = str(self.group_type+'_data/profiles_nosubs/'+str(self._index)+'/'+self.key)
        else:
            path = str(self.group_type+'_data/profiles/'+str(self._index)+'/'+self.key)
        profile_gp = hd[path]
        for k,v in profile_gp.attrs.items():
            if k in self.blacklist:
                continue
            else:
                setattr(self, k, v)
        self.region = {}
        self.region['type'] = profile_gp.attrs['type']
        self.region['centre'] = profile_gp.attrs['centre']
        self.region['axis'] = profile_gp.attrs['axis']
        self.filter = {}
        self.filter['name'] = profile_gp.attrs['name']
        self.filter['conditions'] = profile_gp['conditions'][:]

        for k in profile_gp.keys():
            if k != 'xdata' and k != 'conditions':
                self.yvars[k] = []
                self.ydata[k] = []
                self.weightvars[k] = []
                for j in profile_gp[k].keys():
                    self.yvars[k].append(j)
                    data = profile_gp[k+'/'+j]
                    self.ydata[k].append(unyt_array(data[:], str(data.attrs['units']), registry=self.obj.unit_registry))
                    self.weightvars[k] = list(data.attrs['weightvars'])
            else:
                self.xdata = unyt_array(profile_gp['xdata'][:], profile_gp['xdata'].attrs['units'], registry=self.obj.unit_registry)

class PhaseDiagram:
    def __init__(self,obj,index,group_type,key,hd):
        self.obj = obj
        self._index = index
        self.group_type = group_type
        self.key = key
        self.nbins = [0,0]
        self.xvar = None
        self.yvar = None
        self.region = None
        self.filter = {}
        self.lmax = 0
        self.zvars = {}
        self.weightvars = {}
        self.xdata = None
        self.ydata = None
        self.zdata = {}

        self.blacklist = [
            'obj','_index','group_type','region',
            'filter','zvars','weightvars','xdata','ydata'
        ]
        self._unpack(hd)

    def _unpack(self, hd):
        path = str(self.group_type+'_data/phase_diagrams/'+str(self._index)+'/'+self.key)
        phase_diagrams_gp = hd[path]
        for k,v in phase_diagrams_gp.attrs.items():
            if k in self.blacklist:
                continue
            else:
                setattr(self, k, v)
        self.region = {}
        self.region['type'] = phase_diagrams_gp.attrs['type']
        self.region['centre'] = phase_diagrams_gp.attrs['centre']
        self.region['axis'] = phase_diagrams_gp.attrs['axis']
        self.filter = {}
        self.filter['name'] = phase_diagrams_gp.attrs['name']
        self.filter['conditions'] = phase_diagrams_gp['conditions'][:]
        for k in phase_diagrams_gp.keys():
            if k != 'xdata' and k != 'ydata' and k != 'conditions':
                self.zvars[k] = []
                self.zdata[k] = []
                for j in phase_diagrams_gp[k].keys():
                    if j == 'weightvars':
                        self.weightvars[k] = phase_diagrams_gp[k+'/'+j][:]
                        for w in range(0, len(self.weightvars[k])):
                            self.weightvars[k][w] = self.weightvars[k][w].decode('utf8')
                        continue
                    self.zvars[k].append(j)
                    data = phase_diagrams_gp[k+'/'+j]
                    self.zdata[k].append(unyt_array(data[:], str(data.attrs['units']), registry=self.obj.unit_registry))
            else:
                self.xdata = unyt_array(phase_diagrams_gp['xdata'][:], phase_diagrams_gp['xdata'].attrs['units'], registry=self.obj.unit_registry)
                self.ydata = unyt_array(phase_diagrams_gp['ydata'][:], phase_diagrams_gp['ydata'].attrs['units'], registry=self.obj.unit_registry)

class GalacticFlow:
    def __init__(self,obj,index,group_type,key,hd,nosubs):
        self.obj = obj
        self._index = index
        self.group_type = group_type
        self.key = key
        self.region = None
        self.type = 'none'
        self.filter = {}
        self.data = {}
        self.nosubs = nosubs

        self.blacklist = [
            'obj','_index','group_type','region',
            'filter'
        ]
        self._unpack(hd)

    def _unpack(self,hd):
        if self.nosubs:
            path = str(self.group_type+'_data/flows_nosubs/'+str(self._index)+'/'+self.key)
        else:
            path = str(self.group_type+'_data/flows/'+str(self._index)+'/'+self.key)
        flows_gp = hd[path]
        for k,v in flows_gp.attrs.items():
            if k not in self.blacklist:
                setattr(self, k, v)
        self.region = {}
        self.region['type'] = flows_gp.attrs['type']
        self.region['centre'] = flows_gp.attrs['centre']
        self.region['axis'] = flows_gp.attrs['axis']
        self.region['r'] = unyt_quantity(flows_gp.attrs['r'], 'code_length', registry=self.obj.unit_registry)
        self.region['rmin'] = unyt_quantity(flows_gp.attrs['rmin'], 'code_length', registry=self.obj.unit_registry)
        self.region['rmax'] = unyt_quantity(flows_gp.attrs['rmax'], 'code_length', registry=self.obj.unit_registry)
        self.region['dr'] = unyt_quantity(flows_gp.attrs['dr'], 'code_length', registry=self.obj.unit_registry)
        self.filter = {}
        self.type = flows_gp.attrs['name']
        self.filter['name'] = flows_gp.attrs['name']
        self.filter['conditions'] = flows_gp['conditions'][:]

        for k in flows_gp.keys():
            if k != 'conditions':
                data = flows_gp[k]
                for j in data.keys():
                    unit = data[j].attrs['units']
                    if data[j].shape == 1:
                        if unit =='code_energy':
                            unit = 'code_mass * code_velocity**2'
                        self.data[j] = unyt_quantity(data[j][()], unit, registry=self.obj.unit_registry)
                    else:
                        if unit =='code_energy':
                            unit = 'code_mass * code_velocity**2'
                        self.data[j] = unyt_array(data[j][()], unit, registry=self.obj.unit_registry)

class OZY:
    def __init__(self, filename, read_mode='r'):
        self._ds = None
        self.data_file = os.path.abspath(filename)
        self.snapname = self.data_file.split('/')[-1]
        self.snapID = self.snapname.split('.')[0].split('_')[-1]
        self._galaxy_slist = LazyDataset(self, 'galaxy_data/lists/slist')

        with h5py.File(filename, read_mode) as hd:
            
            # This should be the ozy_version with which the dataset was created.
            self.ozy = hd.attrs['ozy']
            self.unit_registry = UnitRegistry.from_json(
                hd.attrs['unit_registry_json'])
            
            # Load the simulation attributes.
            self.simulation = SimulationAttributes()
            self.simulation._unpack(self, hd)
            # self.simulation._update(self, hd)
            
            
            # Halo data is loaded ALWAYS.
            # TODO: Change this so that it's just a flag.
            self._galaxy_index_list = None
            if 'halo_data/lists/galaxy_index_list' in hd:
                self._galaxy_index_list = LazyDataset(
                    self, 'halo_data/lists/galaxy_index_list')
            
            self._halo_data = {}
            for k, v in hd['halo_data'].items():
                if type(v) is h5py.Dataset:
                    self._halo_data[k] = LazyDataset(self, 'halo_data/' + k)

            self._halo_dicts = defaultdict(dict)
            for k in hd['halo_data/dicts']:
                dictname, arrname = k.split('.')
                self._halo_dicts[dictname][arrname] = LazyDataset(
                    self, 'halo_data/dicts/' + k)
            
            self.nhalos = hd.attrs['nhalos']
            self.halos  = LazyList(self.nhalos, lambda i: Halo(self, i))
            
            # Not all snaps will have galaxies, so we need to load firstly 
            # default values for everything.
            self.have_flows = False
            self.have_flows_nosubs = False
            self.have_profiles = False
            self.have_profiles_nosubs = False
            self._galaxy_data  = {}
            self._galaxy_dicts = defaultdict(dict)
            self.ngalaxies     = 0
            self.galaxies      = LazyList(self.ngalaxies, lambda i: Galaxy(self, i))
            if 'galaxy_data' in hd:
                self._cloud_index_list = None
                if 'galaxy_data/lists/cloud_index_list' in hd:
                    self._cloud_index_list = LazyDataset(
                        self, 'galaxy_data/lists/cloud_index_list')

                if 'tree_data/progen_galaxy_star' in hd:
                    self._galaxy_data['progen_galaxy_star'] = self._progen_galaxy_star = LazyDataset(
                        self, 'tree_data/progen_galaxy_star')
                    
                if 'tree_data/descend_galaxy_star' in hd:
                    self._galaxy_data['descend_galaxy_star'] = self._descend_galaxy_star = LazyDataset(
                        self, 'tree_data/descend_galaxy_star')

                for k, v in hd['galaxy_data'].items():
                    if type(v) is h5py.Dataset:
                        self._galaxy_data[k] = LazyDataset(
                            self, 'galaxy_data/' + k)

                for k in hd['galaxy_data/dicts']:
                    dictname, arrname = k.split('.')
                    self._galaxy_dicts[dictname][arrname] = LazyDataset(
                        self, 'galaxy_data/dicts/' + k)

                self.ngalaxies = hd.attrs['ngalaxies']
                self.galaxies = LazyList(self.ngalaxies,
                                         lambda i: Galaxy(self, i))

                if 'galaxy_data/profiles' in hd:
                    self.have_profiles = True
                    prof_indices = []
                    prof_keys = []
                    for k in hd['galaxy_data/profiles'].keys():
                        for j in hd['galaxy_data/profiles/'+k].keys():
                            prof_indices.append(int(k))
                            prof_keys.append(j)
                    self._galaxy_profile_index_list = LazyList(
                        len(prof_indices), lambda i: int(prof_indices[i])
                    )
                    self._galaxy_profiles = [Profile(self,int(prof_indices[i]),'galaxy', prof_keys[i],hd,False) for i in range(0, len(prof_indices))]
                if 'galaxy_data/profiles_nosubs' in hd:
                    self.have_profiles_nosubs = True
                    prof_nosubs_indices = []
                    prof_nosubs_keys = []
                    for k in hd['galaxy_data/profiles_nosubs'].keys():
                        for j in hd['galaxy_data/profiles_nosubs/'+k].keys():
                            prof_nosubs_indices.append(int(k))
                            prof_nosubs_keys.append(j)
                    self._galaxy_profile_nosubs_index_list = LazyList(
                        len(prof_nosubs_indices), lambda i: int(prof_nosubs_indices[i])
                    )
                    self._galaxy_profiles_nosubs = [Profile(self,int(prof_nosubs_indices[i]),'galaxy', prof_nosubs_keys[i],hd,True) for i in range(0, len(prof_nosubs_indices))]
                if 'galaxy_data/phase_diagrams' in hd:
                    pd_indices = []
                    pd_keys = []
                    for k in hd['galaxy_data/phase_diagrams'].keys():
                        for j in hd['galaxy_data/phase_diagrams/'+k].keys():
                            pd_indices.append(int(k))
                            pd_keys.append(j)
                    self._galaxy_phasediag_index_list = LazyList(
                        len(pd_indices), lambda i: int(pd_indices[i])
                    )
                    self._galaxy_phasediag = [PhaseDiagram(self,int(pd_indices[i]),'galaxy', pd_keys[i],hd) for i in range(0, len(pd_indices))]
                if 'galaxy_data/flows' in hd:
                    self.have_flows = True
                    gf_indices = []
                    gf_keys = []
                    for k in hd['galaxy_data/flows'].keys():
                        for j in hd['galaxy_data/flows/'+k].keys():
                            gf_indices.append(int(k))
                            gf_keys.append(j)
                    self._galaxy_flows_index_list = LazyList(
                        len(gf_indices), lambda i: int(gf_indices[i])
                    )
                    self._galaxy_flows = [GalacticFlow(self,int(gf_indices[i]),'galaxy',gf_keys[i],hd,False) for i in range(0, len(gf_indices))]
                if 'galaxy_data/flows_nosubs' in hd:
                    self.have_flows_nosubs = True
                    gf_nosubs_indices = []
                    gf_nosubs_keys = []
                    for k in hd['galaxy_data/flows_nosubs'].keys():
                        for j in hd['galaxy_data/flows_nosubs/'+k].keys():
                            gf_nosubs_indices.append(int(k))
                            gf_nosubs_keys.append(j)
                    self._galaxy_flows_nosubs_index_list = LazyList(
                        len(gf_nosubs_indices), lambda i: int(gf_nosubs_indices[i])
                    )
                    self._galaxy_flows_nosubs = [GalacticFlow(self,int(gf_nosubs_indices[i]),'galaxy',gf_nosubs_keys[i],hd,True) for i in range(0, len(gf_nosubs_indices))]
    
    @property
    def central_galaxies(self):
        return [h.central_galaxy for h in self.halos]

    @property
    def satellite_galaxies(self):
        galaxies = []
        for h in self.halos:
            galaxies.extend(h.satellite_galaxies)
    
    def most_massive_galaxy(self,return_index=False):
        mstellar = [i.mass['stellar'] for i in self.galaxies]
        ibig = np.argmax(mstellar)
        if return_index:
            return ibig
        else:
            return self.galaxies[ibig]
        
    def most_massive_system(self,return_index=False):
        mvir = [i.virial_quantities['mass'] for i in self.galaxies]
        ibig = np.argmax(mvir)
        if return_index:
            return ibig
        else:
            return self.galaxies[ibig]
    def array(self, value, units):
        return unyt_array(value, units, registry=self.unit_registry)

    def quantity(self, value, units):
        return unyt_quantity(value, units, registry=self.unit_registry)

    def galinfo(self, top=10):
        info_printer(self, 'galaxy', top)

    def haloinfo(self, top=10):
        info_printer(self, 'halo', top)

    def cloudinfo(self, top=10):
        info_printer(self, 'cloud', top)

# TODO: All of these classes need to be adapted to be able to load particle lists.
        
class Group:
    
    def info(self):
        pdict = {}
        for k in getattr(self.obj, '_{}_data'.format(self.type)):
            pdict[k] = getattr(self, k)
        for k in getattr(self.obj, '_{}_dicts'.format(self.type)):
            pdict[k] = dict(getattr(self, k))
        pprint(pdict)
        
class Halo(Group):
    def __init__(self, obj, index):
        self.obj_type = 'halo'
        self.obj = obj
        self._index = index
        self._galaxies = None
        self._satellite_galaxies = None
        self._central_galaxy = None

    def __dir__(self):
        return dir(type(self)) + list(self.__dict__) + list(
            self.obj._halo_data) + list(self.obj._halo_dicts)
        
    @property
    def galaxy_index_list(self):
        return self.obj._galaxy_index_list[self.galaxy_index_list_start:self.
                                           galaxy_index_list_end]

    def _init_galaxies(self):
        self._galaxies = []
        self._satellite_galaxies = []
        for galaxy_index in self.galaxy_index_list:
            galaxy = self.obj.galaxies[galaxy_index]
            self._galaxies.append(galaxy)
            if galaxy.central:
                self._central_galaxy = galaxy
            else:
                self._satellite_galaxies.append(galaxy)

    @property
    def galaxies(self):
        if self._galaxies is None:
            self._init_galaxies()
        return self._galaxies

    @property
    def central_galaxy(self):
        if self._central_galaxy is None:
            self._init_galaxies()
        return self._central_galaxy

    @property
    def satellite_galaxies(self):
        if self._satellite_galaxies is None:
            self._init_galaxies()
        return self._satellite_galaxies

    @property
    def substructure_list(self):
        subs = []
        if self.nextsub == 0:
            print('This halo does not seem to have a substructure assigned!')
            return subs
        haloIDs = self.obj._halo_data['ID'][:]
        nexti = self.nextsub
        while nexti != -1:
            nextindex = np.where(nexti == haloIDs)[0][0]
            subs.append(self.obj.halos[nextindex])
            nexti = self.obj.halos[nextindex].nextsub
        return subs


    @functools.lru_cache(maxsize=None)
    def __getattr__(self, attr):
        if attr in self.obj._halo_data:
            return self.obj._halo_data[attr][self._index]
        if attr in self.obj._halo_dicts:
            return LazyDict(
                self.obj._halo_dicts[attr].keys(),
                lambda d: self.obj._halo_dicts[attr][d][self._index])
        raise AttributeError("'{}' object has no attribute '{}'".format(
            self.__class__.__name__, attr))

class Galaxy(Group):
    def __init__(self, obj, index):
        self.type = 'galaxy'
        self.obj = obj
        self._index = index
        self.halo = obj.halos[self.parent_halo_index]
        self.profiles = None
        self.phase_diagrams = None
        self.flows = None

    def __dir__(self):
        return dir(type(self)) + list(self.__dict__) + list(
            self.obj._galaxy_data) + list(self.obj._galaxy_dicts)

    def _init_profiles(self):
        self.profiles = []
        if self.obj.have_profiles:
            for p,profile_index in enumerate(self.obj._galaxy_profile_index_list):
                if profile_index == self._index:
                    profile = self.obj._galaxy_profiles[p]
                    self.profiles.append(profile)
                
        self.profiles_nosubs = []
        if self.obj.have_profiles_nosubs:
            for p,profile_index in enumerate(self.obj._galaxy_profile_nosubs_index_list):
                if profile_index == self._index:
                    profile = self.obj._galaxy_profiles_nosubs[p]
                    self.profiles_nosubs.append(profile)
    
    def _init_phase_diagrams(self):
        self.phase_diagrams = []
        for p,phasediag_index in enumerate(self.obj._galaxy_phasediag_index_list):
            if phasediag_index == self._index:
                phasediag = self.obj._galaxy_phasediag[p]
                self.phase_diagrams.append(phasediag)

    def _init_flows(self):
        self.flows = []
        if self.obj.have_flows:
            for g,flow_index in enumerate(self.obj._galaxy_flows_index_list):
                if flow_index == self._index:
                    flow = self.obj._galaxy_flows[g]
                    self.flows.append(flow)

        self.flows_nosubs = []
        if self.obj.have_flows_nosubs:
            for g,flow_index in enumerate(self.obj._galaxy_flows_nosubs_index_list):
                if flow_index == self._index:
                    flow = self.obj._galaxy_flows_nosubs[g]
                    self.flows_nosubs.append(flow)

    @property
    def slist(self):
        return self.obj._galaxy_slist[self.slist_start:self.slist_end]

    @property
    def satellites(self):
        if self.central:
            return self.halo.satellite_galaxies
        return []

    @functools.lru_cache(maxsize=None)
    def __getattr__(self, attr):
        if attr in self.obj._galaxy_data:
            return self.obj._galaxy_data[attr][self._index]
        if attr in self.obj._galaxy_dicts:
            return LazyDict(
                self.obj._galaxy_dicts[attr].keys(),
                lambda d: self.obj._galaxy_dicts[attr][d][self._index])
        raise AttributeError("'{}' object has no attribute '{}'".format(
            self.__class__.__name__, attr))
        
    def clear_cache(self):
        """Clear the cache for the __getattr__ method."""
        self.__getattr__.cache_clear()
        
    def update_attribute(self,attr,value):
        if attr in self.obj._galaxy_data:
            # A. The attribute is a single data point (no dict)
            # 1. We update the value in the currently loaded OZY object
            data = self.obj._galaxy_data[attr][:]
            data[self._index] = value
            self.clear_cache()
            
            # 2. Update the original HDF5 file for future reference
            with h5py.File(self.obj.data_file, 'a') as hd:
                _write_attrib(self.obj.galaxies, attr, value, hd['galaxy_data'])
        elif attr in self.obj._galaxy_dicts and isinstance(value, dict):
            # B. The attribute is a dictionary with keys and values
            # 1. Update the full dictionary in the currently loaded OZY object
            for kk, vv in value.items():
                data = self.obj._galaxy_dicts[attr][kk][:]
                data[self._index] = vv
            self.clear_cache()
            
            # 2. Update the original HDF5 file for future reference
            with h5py.File(self.obj.data_file, 'a') as hd:
                _write_dict(self.obj.galaxies, attr, value, hd['galaxy_data/dicts'])
        elif attr.split('.')[0] in self.obj._galaxy_dicts:
            # C. The attribute is a single key of a dictionary already present
            # 1. Update the single key of the dictionary in the currently loaded OZY object
            dictname = attr.split('.')[0]
            keyname = attr.split('.')[1]
            data = self.obj._galaxy_dicts[dictname][keyname][:]
            data[self._index] = value
            self.clear_cache()
            
            # 2. Update the original HDF5 file for future reference
            with h5py.File(self.obj.data_file, 'a') as hd:
                _write_dict(self.obj.galaxies, dictname, self.obj._galaxy_dicts[dictname], hd['galaxy_data/dicts'])
        else:
            # D. This is the case of an inexistent attribute/dict
            #    Choose what type of data is it
            if isinstance(value, dict):
                # This is the case of adding a new dictionary
                # 1. Add the new data to the global OZY object
                for kk, vv in value.items():
                    k = attr + '.' + kk
                    self.obj._galaxy_dicts[attr][kk] = LazyDataset(
                        self.obj, 'galaxy_data/dicts/' + k
                    )
                    if isinstance(vv, unyt_quantity):
                        empty_array = np.full(self.obj.ngalaxies, 0.0)
                        self.obj._galaxy_dicts[attr][kk]._data = self.obj.array(empty_array, vv.units)
                    elif isinstance(vv, unyt_array):
                        empty_array = np.full((self.obj.ngalaxies,)+np.shape(vv), 0.0)
                        self.obj._galaxy_dicts[attr][kk]._data = self.obj.array(empty_array, vv.units)
                    elif isinstance(vv, (np.ndarray,list)):
                        empty_array = np.full((self.obj.ngalaxies,)+np.shape(vv), 0.0)
                        self.obj._galaxy_dicts[attr][kk]._data = empty_array
                    else:
                        empty_array = np.zeros(self.obj.ngalaxies, dtype=type(value))
                        self.obj._galaxy_dicts[attr][kk]._data = empty_array
                    self.obj._galaxy_dicts[attr][kk][self._index] = vv
                    
                # 2. Clean cache to reload Galaxy details
                self.clear_cache()
                
                # 3. Update the original HDF5 file
                with h5py.File(self.obj.data_file, 'r+') as hd:
                    _write_dict(self.obj.galaxies, attr, value, hd['galaxy_data/dicts'])
            elif len(attr.split('.')) > 1:
                # This is the case that we want to add a key to an
                # already existent dictionary
                dictname = attr.split('.')[0]
                keyname = attr.split('.')[1]
                # 1. Add the new data to the global OZY object
                self.obj._galaxy_dicts[dictname][keyname] = LazyDataset(
                        self.obj, 'galaxy_data/dicts/' + attr
                    )
                if isinstance(vv, unyt_quantity):
                    empty_array = np.full(self.obj.ngalaxies, 0.0)
                    self.obj._galaxy_dicts[attr][kk]._data = self.obj.array(empty_array, vv.units)
                elif isinstance(vv, unyt_array):
                    empty_array = np.full((self.obj.ngalaxies,)+np.shape(vv), 0.0)
                    self.obj._galaxy_dicts[attr][kk]._data = self.obj.array(empty_array, vv.units)
                elif isinstance(vv, (np.ndarray,list)):
                    empty_array = np.full((self.obj.ngalaxies,)+np.shape(vv), 0.0)
                    self.obj._galaxy_dicts[attr][kk]._data = empty_array
                else:
                    empty_array = np.zeros(self.obj.ngalaxies, dtype=type(value))
                    self.obj._galaxy_dicts[attr][kk]._data = empty_array
                self.obj._galaxy_dicts[attr][kk][self._index] = value

                # 2. Clean cache to reload Galaxy details
                self.clear_cache()
                
                # 3. Update the original HDF5 file
                with h5py.File(self.obj.data_file, 'r+') as hd:
                    _write_dict(self.obj.galaxies, dictname, self.obj._galaxy_dicts[dictname], hd['galaxy_data/dicts'])
                
            else:
                # This is the final case in which we only want to add
                # a single attribute
                # 1. Add the new data to the global OZY object
                self.obj._galaxy_data[attr] = LazyDataset(
                        self.obj, 'galaxy_data/' + attr
                    )
                if isinstance(value, unyt_quantity):
                    empty_array = np.full(self.obj.ngalaxies, 0.0)
                    self.obj._galaxy_data[attr]._data = self.obj.array(empty_array, value.units)
                elif isinstance(value, unyt_array):
                    empty_array = np.full((self.obj.ngalaxies,)+np.shape(vv), 0.0)
                    self.obj._galaxy_data[attr]._data = self.obj.array(empty_array, value.units)
                elif isinstance(value, (np.ndarray,list)):
                    empty_array = np.full((self.obj.ngalaxies,)+np.shape(vv), 0.0)
                    self.obj._galaxy_data[attr]._data = empty_array
                else:
                    empty_array = np.zeros(self.obj.ngalaxies, dtype=type(value))
                    self.obj._galaxy_data[attr]._data = empty_array
                self.obj._galaxy_data[attr][self._index] = value

                # 2. Clean cache to reload Galaxy details
                self.clear_cache()
                
                with h5py.File(self.obj.data_file, 'r+') as hd:
                    _write_attrib(self.obj.galaxies, attr, value, hd['galaxy_data'])
            

class Cloud(Group):
    def __init__(self, obj, index):
        self.type = 'cloud'
        self.obj = obj
        self._index = index
        self.galaxy = obj.galaxies[self.parent_galaxy_index]
        self.halo = self.galaxy.halo

    def __dir__(self):
        return dir(type(self)) + list(self.__dict__) + list(
            self.obj._cloud_data) + list(self.obj._cloud_dicts)

    @functools.lru_cache(maxsize=None)
    def __getattr__(self, attr):
        if attr in self.obj._cloud_data:
            return self.obj._cloud_data[attr][self._index]
        if attr in self.obj._cloud_dicts:
            return LazyDict(
                self.obj._cloud_dicts[attr].keys(),
                lambda d: self.obj._cloud_dicts[attr][d][self._index])
        raise AttributeError("'{}' object has no attribute '{}'".format(
            self.__class__.__name__, attr))
        

# FINALLY, the function that we want!
def load(filename):
    return OZY(filename)
