
from collections.abc import Mapping, Sequence
from unyt import UnitRegistry,unyt_array,unyt_quantity
import h5py

class LazyDataset:
    """A lazily-loaded HDF5 dataset.
    """
    def __init__(self, obj, dataset_path):
        self._obj = obj
        self._dataset_path = dataset_path
        self._data = None
        self.units = None
    def __getitem__(self, index):
        if self._data is None:
            with h5py.File(self._obj.data_file, 'r') as hd:
                dataset = hd[self._dataset_path]
                if 'unit' in dataset.attrs:
                    self._data = unyt_array(dataset[:],
                                         dataset.attrs['unit'],
                                         registry=self._obj.unit_registry)
                    self.units = dataset.attrs['unit']
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

