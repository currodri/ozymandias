"""
Module amr_integrator


Defined at integrator_module.fpp lines 24-310

"""
from __future__ import print_function, absolute_import, division
import _amr2_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("amr2_pkg.amr_region_attrs")
class amr_region_attrs(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=amr_region_attrs)
    
    
    Defined at integrator_module.fpp lines 28-33
    
    """
    def __init__(self, handle=None):
        """
        self = Amr_Region_Attrs()
        
        
        Defined at integrator_module.fpp lines 28-33
        
        
        Returns
        -------
        this : Amr_Region_Attrs
        	Object to be constructed
        
        
        Automatically generated constructor for amr_region_attrs
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _amr2_pkg.f90wrap_amr_region_attrs_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Amr_Region_Attrs
        
        
        Defined at integrator_module.fpp lines 28-33
        
        Parameters
        ----------
        this : Amr_Region_Attrs
        	Object to be destructed
        
        
        Automatically generated destructor for amr_region_attrs
        """
        if self._alloc:
            _amr2_pkg.f90wrap_amr_region_attrs_finalise(this=self._handle)
    
    @property
    def nvars(self):
        """
        Element nvars ftype=integer  pytype=int
        
        
        Defined at integrator_module.fpp line 29
        
        """
        return _amr2_pkg.f90wrap_amr_region_attrs__get__nvars(self._handle)
    
    @nvars.setter
    def nvars(self, nvars):
        _amr2_pkg.f90wrap_amr_region_attrs__set__nvars(self._handle, nvars)
    
    @property
    def varnames(self):
        """
        Element varnames ftype=character(128) pytype=str
        
        
        Defined at integrator_module.fpp line 30
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_amr_region_attrs__array__varnames(self._handle)
        if array_handle in self._arrays:
            varnames = self._arrays[array_handle]
        else:
            varnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_amr_region_attrs__array__varnames)
            self._arrays[array_handle] = varnames
        return varnames
    
    @varnames.setter
    def varnames(self, varnames):
        self.varnames[...] = varnames
    
    @property
    def nwvars(self):
        """
        Element nwvars ftype=integer  pytype=int
        
        
        Defined at integrator_module.fpp line 31
        
        """
        return _amr2_pkg.f90wrap_amr_region_attrs__get__nwvars(self._handle)
    
    @nwvars.setter
    def nwvars(self, nwvars):
        _amr2_pkg.f90wrap_amr_region_attrs__set__nwvars(self._handle, nwvars)
    
    @property
    def wvarnames(self):
        """
        Element wvarnames ftype=character(128) pytype=str
        
        
        Defined at integrator_module.fpp line 32
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_amr_region_attrs__array__wvarnames(self._handle)
        if array_handle in self._arrays:
            wvarnames = self._arrays[array_handle]
        else:
            wvarnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_amr_region_attrs__array__wvarnames)
            self._arrays[array_handle] = wvarnames
        return wvarnames
    
    @wvarnames.setter
    def wvarnames(self, wvarnames):
        self.wvarnames[...] = wvarnames
    
    @property
    def data(self):
        """
        Element data ftype=real(dbl) pytype=float
        
        
        Defined at integrator_module.fpp line 33
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_amr_region_attrs__array__data(self._handle)
        if array_handle in self._arrays:
            data = self._arrays[array_handle]
        else:
            data = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_amr_region_attrs__array__data)
            self._arrays[array_handle] = data
        return data
    
    @data.setter
    def data(self, data):
        self.data[...] = data
    
    def __str__(self):
        ret = ['<amr_region_attrs>{\n']
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
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def allocate_amr_regions_attrs(self):
    """
    allocate_amr_regions_attrs(self)
    
    
    Defined at integrator_module.fpp lines 36-41
    
    Parameters
    ----------
    attrs : Amr_Region_Attrs
    
    """
    _amr2_pkg.f90wrap_allocate_amr_regions_attrs(attrs=self._handle)

def extract_data(self, varids, pos, cellvars, cellsize, attrs):
    """
    extract_data(self, varids, pos, cellvars, cellsize, attrs)
    
    
    Defined at integrator_module.fpp lines 43-79
    
    Parameters
    ----------
    reg : Region
    varids : Hydroid
    pos : float array
    cellvars : float array
    cellsize : float
    attrs : Amr_Region_Attrs
    
    """
    _amr2_pkg.f90wrap_extract_data(reg=self._handle, varids=varids._handle, pos=pos, \
        cellvars=cellvars, cellsize=cellsize, attrs=attrs._handle)

def renormalise(self):
    """
    renormalise(self)
    
    
    Defined at integrator_module.fpp lines 81-93
    
    Parameters
    ----------
    attrs : Amr_Region_Attrs
    
    """
    _amr2_pkg.f90wrap_renormalise(attrs=self._handle)

def integrate_region(repository, reg, filt, attrs):
    """
    integrate_region(repository, reg, filt, attrs)
    
    
    Defined at integrator_module.fpp lines 95-310
    
    Parameters
    ----------
    repository : str
    reg : Region
    filt : Filter
    attrs : Amr_Region_Attrs
    
    """
    _amr2_pkg.f90wrap_integrate_region(repository=repository, reg=reg._handle, \
        filt=filt._handle, attrs=attrs._handle)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "amr_integrator".')

for func in _dt_array_initialisers:
    func()
