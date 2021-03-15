import os.path
import functools
from pprint import pprint
from collections import defaultdict
from collections.abc import Sequence, Mapping

import h5py
import numpy as np
from yt.units.yt_array import YTArray, UnitRegistry
from ozy.utils import info_printer
from ozy.sim_attributes import SimulationAttributes

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
                    self._data = YTArray(dataset[:],
                                         dataset.attrs['unit'],
                                         registry=self._obj.unit_registry)
                else:
                    self._data = dataset[:]
        return self._data.__getitem__(index)

class LazyList(Sequence):
    """This type should be indistinguishable from the built-in list.
    Any observable difference except the explicit type and performance
    is considered a bug.
    The implementation wraps a list which is intially filled with None,
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

class OZY:
    def __init__(self, filename):
        self._ds = None
        self.data_file = os.path.abspath(filename)
        
        self._galaxy_slist = LazyDataset(self, 'galaxy_data/lists/slist')

        with h5py.File(filename, 'r') as hd:
            
            # This should be the ozy_version with which the dataset was created.
            self.ozy = hd.attrs['ozy']
            self.unit_registry = UnitRegistry.from_json(
                hd.attrs['unit_registry_json'])
            
            # Load the simulation attributes.
            self.simulation = SimulationAttributes()
            self.simulation._unpack(self, hd)
            
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
                

    @property
    def yt_dataset(self):
        """The yt dataset to perform actions on."""
        if self._ds is None:
            raise Exception('No yt_dataset assigned!\nPlease assign '
                            'one via `obj.yt_dataset=<YT DATASET>` '
                            'if you want to do further analysis.')
        return self._ds
    
    @property
    def central_galaxies(self):
        return [h.central_galaxy for h in self.halos]

    @property
    def satellite_galaxies(self):
        galaxies = []
        for h in self.halos:
            galaxies.extend(h.satellite_galaxies)

    def galinfo(self, top=10):
        info_printer(self, 'galaxy', top)

    def haloinfo(self, top=10):
        info_printer(self, 'halo', top)

    def cloudinfo(self, top=10):
        info_printer(self, 'cloud', top)

# TODO: All of these classes need to be addapted to be able to load particle lists.
        
class Group:
    
    def info(self):
        pdict = {}
        for k in getattr(self.obj, '_{}_data'.format(self.obj_type)):
            pdict[k] = getattr(self, k)
        for k in getattr(self.obj, '_{}_dicts'.format(self.obj_type)):
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

    @functools.lru_cache(maxsize=None)
    def __getattr__(self, attr):
        if attr in self.obj._halo_data:
            return self.obj._halo_data[attr][self._index]
        if attr in self.obj._halo_dicts:
            return LazyDict(
                self.obj._halo_dicts[attr].keys(),
                lambda d: self.obj._halo_dicts[attr][d][self._index])
        raise AttributeError("'{}' object as no attribute '{}'".format(
            self.__class__.__name__, attr))

class Galaxy(Group):
    def __init__(self, obj, index):
        self.obj_type = 'galaxy'
        self.obj = obj
        self._index = index
        self.halo = obj.halos[self.parent_halo_index]

    def __dir__(self):
        return dir(type(self)) + list(self.__dict__) + list(
            self.obj._galaxy_data) + list(self.obj._galaxy_dicts)

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
        self.obj_type = 'cloud'
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