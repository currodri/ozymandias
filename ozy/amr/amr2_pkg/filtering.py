"""
Module filtering


Defined at read_amr_module.fpp lines 1345-1437

"""
from __future__ import print_function, absolute_import, division
import _amr2_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("amr2_pkg.filter")
class filter(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=filter)
    
    
    Defined at read_amr_module.fpp lines 1348-1353
    
    """
    def __init__(self, handle=None):
        """
        self = Filter()
        
        
        Defined at read_amr_module.fpp lines 1348-1353
        
        
        Returns
        -------
        this : Filter
        	Object to be constructed
        
        
        Automatically generated constructor for filter
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _amr2_pkg.f90wrap_filter_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Filter
        
        
        Defined at read_amr_module.fpp lines 1348-1353
        
        Parameters
        ----------
        this : Filter
        	Object to be destructed
        
        
        Automatically generated destructor for filter
        """
        if self._alloc:
            _amr2_pkg.f90wrap_filter_finalise(this=self._handle)
    
    @property
    def name(self):
        """
        Element name ftype=character(128) pytype=str
        
        
        Defined at read_amr_module.fpp line 1349
        
        """
        return _amr2_pkg.f90wrap_filter__get__name(self._handle)
    
    @name.setter
    def name(self, name):
        _amr2_pkg.f90wrap_filter__set__name(self._handle, name)
    
    @property
    def ncond(self):
        """
        Element ncond ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 1350
        
        """
        return _amr2_pkg.f90wrap_filter__get__ncond(self._handle)
    
    @ncond.setter
    def ncond(self, ncond):
        _amr2_pkg.f90wrap_filter__set__ncond(self._handle, ncond)
    
    @property
    def cond_vars(self):
        """
        Element cond_vars ftype=character(128) pytype=str
        
        
        Defined at read_amr_module.fpp line 1351
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_filter__array__cond_vars(self._handle)
        if array_handle in self._arrays:
            cond_vars = self._arrays[array_handle]
        else:
            cond_vars = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_filter__array__cond_vars)
            self._arrays[array_handle] = cond_vars
        return cond_vars
    
    @cond_vars.setter
    def cond_vars(self, cond_vars):
        self.cond_vars[...] = cond_vars
    
    @property
    def cond_ops(self):
        """
        Element cond_ops ftype=character(2) pytype=str
        
        
        Defined at read_amr_module.fpp line 1352
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_filter__array__cond_ops(self._handle)
        if array_handle in self._arrays:
            cond_ops = self._arrays[array_handle]
        else:
            cond_ops = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_filter__array__cond_ops)
            self._arrays[array_handle] = cond_ops
        return cond_ops
    
    @cond_ops.setter
    def cond_ops(self, cond_ops):
        self.cond_ops[...] = cond_ops
    
    @property
    def cond_vals(self):
        """
        Element cond_vals ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 1353
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_filter__array__cond_vals(self._handle)
        if array_handle in self._arrays:
            cond_vals = self._arrays[array_handle]
        else:
            cond_vals = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_filter__array__cond_vals)
            self._arrays[array_handle] = cond_vals
        return cond_vals
    
    @cond_vals.setter
    def cond_vals(self, cond_vals):
        self.cond_vals[...] = cond_vals
    
    def __str__(self):
        ret = ['<filter>{\n']
        ret.append('    name : ')
        ret.append(repr(self.name))
        ret.append(',\n    ncond : ')
        ret.append(repr(self.ncond))
        ret.append(',\n    cond_vars : ')
        ret.append(repr(self.cond_vars))
        ret.append(',\n    cond_ops : ')
        ret.append(repr(self.cond_ops))
        ret.append(',\n    cond_vals : ')
        ret.append(repr(self.cond_vals))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def allocate_filter(self):
    """
    allocate_filter(self)
    
    
    Defined at read_amr_module.fpp lines 1356-1361
    
    Parameters
    ----------
    filt : Filter
    
    """
    _amr2_pkg.f90wrap_allocate_filter(filt=self._handle)

def cond_string_to_filter(str, filt):
    """
    cond_string_to_filter(str, filt)
    
    
    Defined at read_amr_module.fpp lines 1363-1367
    
    Parameters
    ----------
    str : str
    filt : Filter
    
    """
    _amr2_pkg.f90wrap_cond_string_to_filter(str=str, filt=filt._handle)

def filter_cell(self, filt, cell_x, cell_dx, cell_var):
    """
    filter_cell = filter_cell(self, filt, cell_x, cell_dx, cell_var)
    
    
    Defined at read_amr_module.fpp lines 1369-1401
    
    Parameters
    ----------
    reg : Region
    filt : Filter
    cell_x : Vector
    cell_dx : float
    cell_var : float array
    
    Returns
    -------
    filter_cell : bool
    
    """
    filter_cell = _amr2_pkg.f90wrap_filter_cell(reg=self._handle, filt=filt._handle, \
        cell_x=cell_x._handle, cell_dx=cell_dx, cell_var=cell_var)
    return filter_cell

def filter_particle(self, filt, part, dx=None):
    """
    filter_particle = filter_particle(self, filt, part[, dx])
    
    
    Defined at read_amr_module.fpp lines 1403-1437
    
    Parameters
    ----------
    reg : Region
    filt : Filter
    part : Particle
    dx : Vector
    
    Returns
    -------
    filter_particle : bool
    
    """
    filter_particle = _amr2_pkg.f90wrap_filter_particle(reg=self._handle, \
        filt=filt._handle, part=part._handle, dx=None if dx is None else dx._handle)
    return filter_particle


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "filtering".')

for func in _dt_array_initialisers:
    func()
