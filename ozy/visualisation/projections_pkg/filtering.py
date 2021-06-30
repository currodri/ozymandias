"""
Module filtering


Defined at read_amr_module.fpp lines 882-972

"""
from __future__ import print_function, absolute_import, division
import _projections_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("projections_pkg.filter")
class filter(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=filter)
    
    
    Defined at read_amr_module.fpp lines 885-890
    
    """
    def __init__(self, handle=None):
        """
        self = Filter()
        
        
        Defined at read_amr_module.fpp lines 885-890
        
        
        Returns
        -------
        this : Filter
        	Object to be constructed
        
        
        Automatically generated constructor for filter
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _projections_pkg.f90wrap_filter_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Filter
        
        
        Defined at read_amr_module.fpp lines 885-890
        
        Parameters
        ----------
        this : Filter
        	Object to be destructed
        
        
        Automatically generated destructor for filter
        """
        if self._alloc:
            _projections_pkg.f90wrap_filter_finalise(this=self._handle)
    
    @property
    def name(self):
        """
        Element name ftype=character(128) pytype=str
        
        
        Defined at read_amr_module.fpp line 886
        
        """
        return _projections_pkg.f90wrap_filter__get__name(self._handle)
    
    @name.setter
    def name(self, name):
        _projections_pkg.f90wrap_filter__set__name(self._handle, name)
    
    @property
    def ncond(self):
        """
        Element ncond ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 887
        
        """
        return _projections_pkg.f90wrap_filter__get__ncond(self._handle)
    
    @ncond.setter
    def ncond(self, ncond):
        _projections_pkg.f90wrap_filter__set__ncond(self._handle, ncond)
    
    @property
    def cond_vars(self):
        """
        Element cond_vars ftype=character(128) pytype=str
        
        
        Defined at read_amr_module.fpp line 888
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _projections_pkg.f90wrap_filter__array__cond_vars(self._handle)
        if array_handle in self._arrays:
            cond_vars = self._arrays[array_handle]
        else:
            cond_vars = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _projections_pkg.f90wrap_filter__array__cond_vars)
            self._arrays[array_handle] = cond_vars
        return cond_vars
    
    @cond_vars.setter
    def cond_vars(self, cond_vars):
        self.cond_vars[...] = cond_vars
    
    @property
    def cond_ops(self):
        """
        Element cond_ops ftype=character(2) pytype=str
        
        
        Defined at read_amr_module.fpp line 889
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _projections_pkg.f90wrap_filter__array__cond_ops(self._handle)
        if array_handle in self._arrays:
            cond_ops = self._arrays[array_handle]
        else:
            cond_ops = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _projections_pkg.f90wrap_filter__array__cond_ops)
            self._arrays[array_handle] = cond_ops
        return cond_ops
    
    @cond_ops.setter
    def cond_ops(self, cond_ops):
        self.cond_ops[...] = cond_ops
    
    @property
    def cond_vals(self):
        """
        Element cond_vals ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 890
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _projections_pkg.f90wrap_filter__array__cond_vals(self._handle)
        if array_handle in self._arrays:
            cond_vals = self._arrays[array_handle]
        else:
            cond_vals = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _projections_pkg.f90wrap_filter__array__cond_vals)
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
    
    
    Defined at read_amr_module.fpp lines 893-898
    
    Parameters
    ----------
    filt : Filter
    
    """
    _projections_pkg.f90wrap_allocate_filter(filt=self._handle)

def cond_string_to_filter(str, filt):
    """
    cond_string_to_filter(str, filt)
    
    
    Defined at read_amr_module.fpp lines 900-904
    
    Parameters
    ----------
    str : str
    filt : Filter
    
    """
    _projections_pkg.f90wrap_cond_string_to_filter(str=str, filt=filt._handle)

def filter_cell(self, reg, filt, cell_x, cell_dx, cell_var):
    """
    filter_cell = filter_cell(self, reg, filt, cell_x, cell_dx, cell_var)
    
    
    Defined at read_amr_module.fpp lines 906-939
    
    Parameters
    ----------
    varids : Hydroid
    reg : Region
    filt : Filter
    cell_x : Vector
    cell_dx : float
    cell_var : float array
    
    Returns
    -------
    filter_cell : bool
    
    """
    filter_cell = _projections_pkg.f90wrap_filter_cell(varids=self._handle, \
        reg=reg._handle, filt=filt._handle, cell_x=cell_x._handle, cell_dx=cell_dx, \
        cell_var=cell_var)
    return filter_cell

def filter_particle(self, reg, filt, dx, part):
    """
    filter_particle = filter_particle(self, reg, filt, dx, part)
    
    
    Defined at read_amr_module.fpp lines 941-972
    
    Parameters
    ----------
    sim : Sim_Info
    reg : Region
    filt : Filter
    dx : Vector
    part : Particle
    
    Returns
    -------
    filter_particle : bool
    
    """
    filter_particle = _projections_pkg.f90wrap_filter_particle(sim=self._handle, \
        reg=reg._handle, filt=filt._handle, dx=dx._handle, part=part._handle)
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
