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
    def __init__(self, obj, index, group_type, key,hd):
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

        self.blacklist = [
            'obj','_index','group_type','region',
            'filter','yvars','weightvars','xdata','ydata'
        ]
        self._unpack(hd)
    
    def _unpack(self, hd):
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
    def __init__(self,obj,index,group_type,key,hd):
        self.obj = obj
        self._index = index
        self.group_type = group_type
        self.key = key
        self.region = None
        self.type = 'none'
        self.filter = {}
        self.data = {}

        self.blacklist = [
            'obj','_index','group_type','region',
            'filter'
        ]
        self._unpack(hd)

    def _unpack(self,hd):
        path = str(self.group_type+'_data/flows/'+str(self._index)+'/'+self.key)
        flows_gp = hd[path]
        for k,v in flows_gp.attrs.items():
            if k not in self.blacklist:
                setattr(self, k, v)
        self.region = {}
        self.region['type'] = flows_gp.attrs['type']
        self.region['centre'] = flows_gp.attrs['centre']
        self.region['axis'] = flows_gp.attrs['axis']
        self.filter = {}
        self.type = flows_gp.attrs['name']
        self.filter['name'] = flows_gp.attrs['name']
        self.filter['conditions'] = flows_gp['conditions'][:]

        for k in flows_gp.keys():
            if k != 'conditions':
                data = flows_gp[k]
                for j in data.keys():
                    unit = data[j].attrs['units']
                    if unit =='code_energy':
                        unit = 'code_mass * code_velocity**2'
                    self.data[j] = unyt_quantity(data[j][()], unit, registry=self.obj.unit_registry)

class OZY:
    def __init__(self, filename, read_mode='r'):
        self._ds = None
        self.data_file = os.path.abspath(filename)
        
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
                    prof_indices = []
                    prof_keys = []
                    for k in hd['galaxy_data/profiles'].keys():
                        for j in hd['galaxy_data/profiles/'+k].keys():
                            prof_indices.append(int(k))
                            prof_keys.append(j)
                    self._galaxy_profile_index_list = LazyList(
                        len(prof_indices), lambda i: int(prof_indices[i])
                    )
                    self._galaxy_profiles = [Profile(self,int(prof_indices[i]),'galaxy', prof_keys[i],hd) for i in range(0, len(prof_indices))]
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
                    gf_indices = []
                    gf_keys = []
                    for k in hd['galaxy_data/flows'].keys():
                        for j in hd['galaxy_data/flows/'+k].keys():
                            gf_indices.append(int(k))
                            gf_keys.append(j)
                    self._galaxy_flows_index_list = LazyList(
                        len(gf_indices), lambda i: int(gf_indices[i])
                    )
                    self._galaxy_flows = [GalacticFlow(self,int(gf_indices[i]),'galaxy',gf_keys[i],hd) for i in range(0, len(gf_indices))]
    
    @property
    def central_galaxies(self):
        return [h.central_galaxy for h in self.halos]

    @property
    def satellite_galaxies(self):
        galaxies = []
        for h in self.halos:
            galaxies.extend(h.satellite_galaxies)
    
    @property
    def most_massive_galaxy(self):
        mstellar = [i.mass['stellar'] for i in self.galaxies]
        ibig = np.argmax(mstellar)
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
        for p,profile_index in enumerate(self.obj._galaxy_profile_index_list):
            if profile_index == self._index:
                profile = self.obj._galaxy_profiles[p]
                self.profiles.append(profile)
    
    def _init_phase_diagrams(self):
        self.phase_diagrams = []
        for p,phasediag_index in enumerate(self.obj._galaxy_phasediag_index_list):
            if phasediag_index == self._index:
                phasediag = self.obj._galaxy_phasediag[p]
                self.phase_diagrams.append(phasediag)

    def _init_flows(self):
        self.flows = []
        for g,flow_index in enumerate(self.obj._galaxy_flows_index_list):
            if flow_index == self._index:
                flow = self.obj._galaxy_flows[g]
                self.flows.append(flow)

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
