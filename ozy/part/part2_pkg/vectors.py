"""
Module vectors


Defined at linalg_module.fpp lines 23-159

"""
from __future__ import print_function, absolute_import, division
import _part2_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("part2_pkg.vector")
class vector(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=vector)
    
    
    Defined at linalg_module.fpp lines 26-29
    
    """
    def __init__(self, handle=None):
        """
        self = Vector()
        
        
        Defined at linalg_module.fpp lines 26-29
        
        
        Returns
        -------
        this : Vector
        	Object to be constructed
        
        
        Automatically generated constructor for vector
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _part2_pkg.f90wrap_vector_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Vector
        
        
        Defined at linalg_module.fpp lines 26-29
        
        Parameters
        ----------
        this : Vector
        	Object to be destructed
        
        
        Automatically generated destructor for vector
        """
        if self._alloc:
            _part2_pkg.f90wrap_vector_finalise(this=self._handle)
    
    @property
    def x(self):
        """
        Element x ftype=real(dbl) pytype=float
        
        
        Defined at linalg_module.fpp line 27
        
        """
        return _part2_pkg.f90wrap_vector__get__x(self._handle)
    
    @x.setter
    def x(self, x):
        _part2_pkg.f90wrap_vector__set__x(self._handle, x)
    
    @property
    def y(self):
        """
        Element y ftype=real(dbl) pytype=float
        
        
        Defined at linalg_module.fpp line 28
        
        """
        return _part2_pkg.f90wrap_vector__get__y(self._handle)
    
    @y.setter
    def y(self, y):
        _part2_pkg.f90wrap_vector__set__y(self._handle, y)
    
    @property
    def z(self):
        """
        Element z ftype=real(dbl) pytype=float
        
        
        Defined at linalg_module.fpp line 29
        
        """
        return _part2_pkg.f90wrap_vector__get__z(self._handle)
    
    @z.setter
    def z(self, z):
        _part2_pkg.f90wrap_vector__set__z(self._handle, z)
    
    def __str__(self):
        ret = ['<vector>{\n']
        ret.append('    x : ')
        ret.append(repr(self.x))
        ret.append(',\n    y : ')
        ret.append(repr(self.y))
        ret.append(',\n    z : ')
        ret.append(repr(self.z))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("part2_pkg.array_vectors")
class array_vectors(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=array_vectors)
    
    
    Defined at linalg_module.fpp lines 31-33
    
    """
    def __init__(self, handle=None):
        """
        self = Array_Vectors()
        
        
        Defined at linalg_module.fpp lines 31-33
        
        
        Returns
        -------
        this : Array_Vectors
        	Object to be constructed
        
        
        Automatically generated constructor for array_vectors
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _part2_pkg.f90wrap_array_vectors_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Array_Vectors
        
        
        Defined at linalg_module.fpp lines 31-33
        
        Parameters
        ----------
        this : Array_Vectors
        	Object to be destructed
        
        
        Automatically generated destructor for array_vectors
        """
        if self._alloc:
            _part2_pkg.f90wrap_array_vectors_finalise(this=self._handle)
    
    @property
    def n(self):
        """
        Element n ftype=integer  pytype=int
        
        
        Defined at linalg_module.fpp line 32
        
        """
        return _part2_pkg.f90wrap_array_vectors__get__n(self._handle)
    
    @n.setter
    def n(self, n):
        _part2_pkg.f90wrap_array_vectors__set__n(self._handle, n)
    
    def init_array_list(self):
        self.list = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _part2_pkg.f90wrap_array_vectors__array_getitem__list,
                                        _part2_pkg.f90wrap_array_vectors__array_setitem__list,
                                        _part2_pkg.f90wrap_array_vectors__array_len__list,
                                        """
        Element list ftype=type(vector) pytype=Vector
        
        
        Defined at linalg_module.fpp line 33
        
        """, vector)
        return self.list
    
    def __str__(self):
        ret = ['<array_vectors>{\n']
        ret.append('    n : ')
        ret.append(repr(self.n))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = [init_array_list]
    

def magnitude(self):
    """
    magnitude = magnitude(self)
    
    
    Defined at linalg_module.fpp lines 150-152
    
    Parameters
    ----------
    vec_1 : Vector
    
    Returns
    -------
    magnitude : float
    
    """
    magnitude = _part2_pkg.f90wrap_magnitude(vec_1=self._handle)
    return magnitude

def _array_to_vector(self, array):
    """
    _array_to_vector(self, array)
    
    
    Defined at linalg_module.fpp lines 64-69
    
    Parameters
    ----------
    vec_result : Vector
    array : float array
    
    """
    _part2_pkg.f90wrap_array_to_vector(vec_result=self._handle, array=array)

def _vector_to_array(array_result, vec_1):
    """
    _vector_to_array(array_result, vec_1)
    
    
    Defined at linalg_module.fpp lines 71-76
    
    Parameters
    ----------
    array_result : float array
    vec_1 : Vector
    
    """
    _part2_pkg.f90wrap_vector_to_array(array_result=array_result, \
        vec_1=vec_1._handle)

def assignment(*args, **kwargs):
    """
    assignment(*args, **kwargs)
    
    
    Defined at linalg_module.fpp lines 35-37
    
    Overloaded interface containing the following procedures:
      _array_to_vector
      _vector_to_array
    
    """
    for proc in [_array_to_vector, _vector_to_array]:
        try:
            return proc(*args, **kwargs)
        except TypeError:
            continue


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "vectors".')

for func in _dt_array_initialisers:
    func()
