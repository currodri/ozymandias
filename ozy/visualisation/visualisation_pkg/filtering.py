"""
Module filtering


Defined at read_amr_module.fpp lines 2516-2648

"""
from __future__ import print_function, absolute_import, division
import _visualisation_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("visualisation_pkg.filter")
class filter(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=filter)
    
    
    Defined at read_amr_module.fpp lines 2519-2526
    
    """
    def __init__(self, handle=None):
        """
        self = Filter()
        
        
        Defined at read_amr_module.fpp lines 2519-2526
        
        
        Returns
        -------
        this : Filter
        	Object to be constructed
        
        
        Automatically generated constructor for filter
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _visualisation_pkg.f90wrap_filter_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Filter
        
        
        Defined at read_amr_module.fpp lines 2519-2526
        
        Parameters
        ----------
        this : Filter
        	Object to be destructed
        
        
        Automatically generated destructor for filter
        """
        if self._alloc:
            _visualisation_pkg.f90wrap_filter_finalise(this=self._handle)
    
    @property
    def name(self):
        """
        Element name ftype=character(128) pytype=str
        
        
        Defined at read_amr_module.fpp line 2520
        
        """
        return _visualisation_pkg.f90wrap_filter__get__name(self._handle)
    
    @name.setter
    def name(self, name):
        _visualisation_pkg.f90wrap_filter__set__name(self._handle, name)
    
    @property
    def ncond(self):
        """
        Element ncond ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 2521
        
        """
        return _visualisation_pkg.f90wrap_filter__get__ncond(self._handle)
    
    @ncond.setter
    def ncond(self, ncond):
        _visualisation_pkg.f90wrap_filter__set__ncond(self._handle, ncond)
    
    @property
    def cond_vars(self):
        """
        Element cond_vars ftype=character(128) pytype=str
        
        
        Defined at read_amr_module.fpp line 2522
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_filter__array__cond_vars(self._handle)
        if array_handle in self._arrays:
            cond_vars = self._arrays[array_handle]
        else:
            cond_vars = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_filter__array__cond_vars)
            self._arrays[array_handle] = cond_vars
        return cond_vars
    
    @cond_vars.setter
    def cond_vars(self, cond_vars):
        self.cond_vars[...] = cond_vars
    
    @property
    def cond_vars_comp(self):
        """
        Element cond_vars_comp ftype=character(128) pytype=str
        
        
        Defined at read_amr_module.fpp line 2523
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_filter__array__cond_vars_comp(self._handle)
        if array_handle in self._arrays:
            cond_vars_comp = self._arrays[array_handle]
        else:
            cond_vars_comp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_filter__array__cond_vars_comp)
            self._arrays[array_handle] = cond_vars_comp
        return cond_vars_comp
    
    @cond_vars_comp.setter
    def cond_vars_comp(self, cond_vars_comp):
        self.cond_vars_comp[...] = cond_vars_comp
    
    @property
    def cond_ops(self):
        """
        Element cond_ops ftype=character(2) pytype=str
        
        
        Defined at read_amr_module.fpp line 2524
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_filter__array__cond_ops(self._handle)
        if array_handle in self._arrays:
            cond_ops = self._arrays[array_handle]
        else:
            cond_ops = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_filter__array__cond_ops)
            self._arrays[array_handle] = cond_ops
        return cond_ops
    
    @cond_ops.setter
    def cond_ops(self, cond_ops):
        self.cond_ops[...] = cond_ops
    
    @property
    def cond_vals(self):
        """
        Element cond_vals ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 2525
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_filter__array__cond_vals(self._handle)
        if array_handle in self._arrays:
            cond_vals = self._arrays[array_handle]
        else:
            cond_vals = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_filter__array__cond_vals)
            self._arrays[array_handle] = cond_vals
        return cond_vals
    
    @cond_vals.setter
    def cond_vals(self, cond_vals):
        self.cond_vals[...] = cond_vals
    
    @property
    def use_var(self):
        """
        Element use_var ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 2526
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_filter__array__use_var(self._handle)
        if array_handle in self._arrays:
            use_var = self._arrays[array_handle]
        else:
            use_var = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_filter__array__use_var)
            self._arrays[array_handle] = use_var
        return use_var
    
    @use_var.setter
    def use_var(self, use_var):
        self.use_var[...] = use_var
    
    def __str__(self):
        ret = ['<filter>{\n']
        ret.append('    name : ')
        ret.append(repr(self.name))
        ret.append(',\n    ncond : ')
        ret.append(repr(self.ncond))
        ret.append(',\n    cond_vars : ')
        ret.append(repr(self.cond_vars))
        ret.append(',\n    cond_vars_comp : ')
        ret.append(repr(self.cond_vars_comp))
        ret.append(',\n    cond_ops : ')
        ret.append(repr(self.cond_ops))
        ret.append(',\n    cond_vals : ')
        ret.append(repr(self.cond_vals))
        ret.append(',\n    use_var : ')
        ret.append(repr(self.use_var))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def allocate_filter(self):
    """
    allocate_filter(self)
    
    
    Defined at read_amr_module.fpp lines 2529-2537
    
    Parameters
    ----------
    filt : Filter
    
    """
    _visualisation_pkg.f90wrap_allocate_filter(filt=self._handle)

def cond_string_to_filter(str, filt):
    """
    cond_string_to_filter(str, filt)
    
    
    Defined at read_amr_module.fpp lines 2539-2544
    
    Parameters
    ----------
    str : str
    filt : Filter
    
    """
    _visualisation_pkg.f90wrap_cond_string_to_filter(str=str, filt=filt._handle)

def filter_cell(self, filt, cell_x, cell_dx, cell_var, cell_son, trans_matrix, \
    grav_var=None):
    """
    filter_cell = filter_cell(self, filt, cell_x, cell_dx, cell_var, cell_son, \
        trans_matrix[, grav_var])
    
    
    Defined at read_amr_module.fpp lines 2546-2599
    
    Parameters
    ----------
    reg : Region
    filt : Filter
    cell_x : Vector
    cell_dx : float
    cell_var : float array
    cell_son : int array
    trans_matrix : float array
    grav_var : float array
    
    Returns
    -------
    filter_cell : bool
    
    """
    filter_cell = _visualisation_pkg.f90wrap_filter_cell(reg=self._handle, \
        filt=filt._handle, cell_x=cell_x._handle, cell_dx=cell_dx, \
        cell_var=cell_var, cell_son=cell_son, trans_matrix=trans_matrix, \
        grav_var=grav_var)
    return filter_cell

def filter_particle(self, filt, part, dx=None):
    """
    filter_particle = filter_particle(self, filt, part[, dx])
    
    
    Defined at read_amr_module.fpp lines 2601-2635
    
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
    filter_particle = _visualisation_pkg.f90wrap_filter_particle(reg=self._handle, \
        filt=filt._handle, part=part._handle, dx=None if dx is None else dx._handle)
    return filter_particle

def filter_sub(self, cell_x):
    """
    filter_sub = filter_sub(self, cell_x)
    
    
    Defined at read_amr_module.fpp lines 2637-2648
    
    Parameters
    ----------
    sub : Region
    cell_x : float array
    
    Returns
    -------
    filter_sub : bool
    
    """
    filter_sub = _visualisation_pkg.f90wrap_filter_sub(sub=self._handle, \
        cell_x=cell_x)
    return filter_sub


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "filtering".')

for func in _dt_array_initialisers:
    func()
