"""
Module part_integrator


Defined at integrator_module.fpp lines 24-299

"""
from __future__ import print_function, absolute_import, division
import _part2_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("part2_pkg.part_region_attrs")
class part_region_attrs(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=part_region_attrs)
    
    
    Defined at integrator_module.fpp lines 29-36
    
    """
    def __init__(self, handle=None):
        """
        self = Part_Region_Attrs()
        
        
        Defined at integrator_module.fpp lines 29-36
        
        
        Returns
        -------
        this : Part_Region_Attrs
        	Object to be constructed
        
        
        Automatically generated constructor for part_region_attrs
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _part2_pkg.f90wrap_part_region_attrs_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Part_Region_Attrs
        
        
        Defined at integrator_module.fpp lines 29-36
        
        Parameters
        ----------
        this : Part_Region_Attrs
        	Object to be destructed
        
        
        Automatically generated destructor for part_region_attrs
        """
        if self._alloc:
            _part2_pkg.f90wrap_part_region_attrs_finalise(this=self._handle)
    
    @property
    def nvars(self):
        """
        Element nvars ftype=integer  pytype=int
        
        
        Defined at integrator_module.fpp line 30
        
        """
        return _part2_pkg.f90wrap_part_region_attrs__get__nvars(self._handle)
    
    @nvars.setter
    def nvars(self, nvars):
        _part2_pkg.f90wrap_part_region_attrs__set__nvars(self._handle, nvars)
    
    @property
    def varnames(self):
        """
        Element varnames ftype=character(128) pytype=str
        
        
        Defined at integrator_module.fpp line 31
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _part2_pkg.f90wrap_part_region_attrs__array__varnames(self._handle)
        if array_handle in self._arrays:
            varnames = self._arrays[array_handle]
        else:
            varnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _part2_pkg.f90wrap_part_region_attrs__array__varnames)
            self._arrays[array_handle] = varnames
        return varnames
    
    @varnames.setter
    def varnames(self, varnames):
        self.varnames[...] = varnames
    
    @property
    def nwvars(self):
        """
        Element nwvars ftype=integer  pytype=int
        
        
        Defined at integrator_module.fpp line 32
        
        """
        return _part2_pkg.f90wrap_part_region_attrs__get__nwvars(self._handle)
    
    @nwvars.setter
    def nwvars(self, nwvars):
        _part2_pkg.f90wrap_part_region_attrs__set__nwvars(self._handle, nwvars)
    
    @property
    def wvarnames(self):
        """
        Element wvarnames ftype=character(128) pytype=str
        
        
        Defined at integrator_module.fpp line 33
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _part2_pkg.f90wrap_part_region_attrs__array__wvarnames(self._handle)
        if array_handle in self._arrays:
            wvarnames = self._arrays[array_handle]
        else:
            wvarnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _part2_pkg.f90wrap_part_region_attrs__array__wvarnames)
            self._arrays[array_handle] = wvarnames
        return wvarnames
    
    @wvarnames.setter
    def wvarnames(self, wvarnames):
        self.wvarnames[...] = wvarnames
    
    @property
    def data(self):
        """
        Element data ftype=real(dbl) pytype=float
        
        
        Defined at integrator_module.fpp line 34
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _part2_pkg.f90wrap_part_region_attrs__array__data(self._handle)
        if array_handle in self._arrays:
            data = self._arrays[array_handle]
        else:
            data = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _part2_pkg.f90wrap_part_region_attrs__array__data)
            self._arrays[array_handle] = data
        return data
    
    @data.setter
    def data(self, data):
        self.data[...] = data
    
    @property
    def ndm(self):
        """
        Element ndm ftype=integer(irg) pytype=int
        
        
        Defined at integrator_module.fpp line 35
        
        """
        return _part2_pkg.f90wrap_part_region_attrs__get__ndm(self._handle)
    
    @ndm.setter
    def ndm(self, ndm):
        _part2_pkg.f90wrap_part_region_attrs__set__ndm(self._handle, ndm)
    
    @property
    def nstar(self):
        """
        Element nstar ftype=integer(irg) pytype=int
        
        
        Defined at integrator_module.fpp line 35
        
        """
        return _part2_pkg.f90wrap_part_region_attrs__get__nstar(self._handle)
    
    @nstar.setter
    def nstar(self, nstar):
        _part2_pkg.f90wrap_part_region_attrs__set__nstar(self._handle, nstar)
    
    @property
    def nids(self):
        """
        Element nids ftype=integer(irg) pytype=int
        
        
        Defined at integrator_module.fpp line 35
        
        """
        return _part2_pkg.f90wrap_part_region_attrs__get__nids(self._handle)
    
    @nids.setter
    def nids(self, nids):
        _part2_pkg.f90wrap_part_region_attrs__set__nids(self._handle, nids)
    
    @property
    def ids(self):
        """
        Element ids ftype=integer(ilg) pytype=int
        
        
        Defined at integrator_module.fpp line 36
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _part2_pkg.f90wrap_part_region_attrs__array__ids(self._handle)
        if array_handle in self._arrays:
            ids = self._arrays[array_handle]
        else:
            ids = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _part2_pkg.f90wrap_part_region_attrs__array__ids)
            self._arrays[array_handle] = ids
        return ids
    
    @ids.setter
    def ids(self, ids):
        self.ids[...] = ids
    
    def __str__(self):
        ret = ['<part_region_attrs>{\n']
        ret.append('    nvars : ')
        ret.append(repr(self.nvars))
        ret.append(',\n    varnames : ')
        ret.append(repr(self.varnames))
        ret.append(',\n    nwvars : ')
        ret.append(repr(self.nwvars))
        ret.append(',\n    wvarnames : ')
        ret.append(repr(self.wvarnames))
        ret.append(',\n    data : ')
        ret.append(repr(self.data))
        ret.append(',\n    ndm : ')
        ret.append(repr(self.ndm))
        ret.append(',\n    nstar : ')
        ret.append(repr(self.nstar))
        ret.append(',\n    nids : ')
        ret.append(repr(self.nids))
        ret.append(',\n    ids : ')
        ret.append(repr(self.ids))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def allocate_part_regions_attrs(self):
    """
    allocate_part_regions_attrs(self)
    
    
    Defined at integrator_module.fpp lines 39-44
    
    Parameters
    ----------
    attrs : Part_Region_Attrs
    
    """
    _part2_pkg.f90wrap_allocate_part_regions_attrs(attrs=self._handle)

def extract_data(self, part, attrs):
    """
    extract_data(self, part, attrs)
    
    
    Defined at integrator_module.fpp lines 46-88
    
    Parameters
    ----------
    reg : Region
    part : Particle
    attrs : Part_Region_Attrs
    
    """
    _part2_pkg.f90wrap_extract_data(reg=self._handle, part=part._handle, \
        attrs=attrs._handle)

def renormalise(self):
    """
    renormalise(self)
    
    
    Defined at integrator_module.fpp lines 90-113
    
    Parameters
    ----------
    attrs : Part_Region_Attrs
    
    """
    _part2_pkg.f90wrap_renormalise(attrs=self._handle)

def integrate_region(repository, reg, filt, attrs, get_ids=None):
    """
    integrate_region(repository, reg, filt, attrs[, get_ids])
    
    
    Defined at integrator_module.fpp lines 115-299
    
    Parameters
    ----------
    repository : str
    reg : Region
    filt : Filter
    attrs : Part_Region_Attrs
    get_ids : bool
    
    """
    _part2_pkg.f90wrap_integrate_region(repository=repository, reg=reg._handle, \
        filt=filt._handle, attrs=attrs._handle, get_ids=get_ids)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "part_integrator".')

for func in _dt_array_initialisers:
    func()
