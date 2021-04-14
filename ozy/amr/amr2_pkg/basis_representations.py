"""
Module basis_representations


Defined at linalg_module.fpp lines 257-290

"""
from __future__ import print_function, absolute_import, division
import _amr2_pkg
import f90wrap.runtime
import logging
from amr2_pkg.vectors import vector

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("amr2_pkg.basis")
class basis(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=basis)
    
    
    Defined at linalg_module.fpp lines 260-261
    
    """
    def __init__(self, handle=None):
        """
        self = Basis()
        
        
        Defined at linalg_module.fpp lines 260-261
        
        
        Returns
        -------
        this : Basis
        	Object to be constructed
        
        
        Automatically generated constructor for basis
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _amr2_pkg.f90wrap_basis_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Basis
        
        
        Defined at linalg_module.fpp lines 260-261
        
        Parameters
        ----------
        this : Basis
        	Object to be destructed
        
        
        Automatically generated destructor for basis
        """
        if self._alloc:
            _amr2_pkg.f90wrap_basis_finalise(this=self._handle)
    
    def init_array_u(self):
        self.u = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _amr2_pkg.f90wrap_basis__array_getitem__u,
                                        _amr2_pkg.f90wrap_basis__array_setitem__u,
                                        _amr2_pkg.f90wrap_basis__array_len__u,
                                        """
        Element u ftype=type(vector) pytype=Vector
        
        
        Defined at linalg_module.fpp line 261
        
        """, vector)
        return self.u
    
    _dt_array_initialisers = [init_array_u]
    

def initialise_basis(self):
    """
    initialise_basis(self)
    
    
    Defined at linalg_module.fpp lines 264-268
    
    Parameters
    ----------
    this : Basis
    
    """
    _amr2_pkg.f90wrap_initialise_basis(this=self._handle)

def mgramschmidt(self, e):
    """
    mgramschmidt(self, e)
    
    
    Defined at linalg_module.fpp lines 277-289
    
    Parameters
    ----------
    vecs : Basis
    e : Basis
    
    """
    _amr2_pkg.f90wrap_mgramschmidt(vecs=self._handle, e=e._handle)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "basis_representations".')

for func in _dt_array_initialisers:
    func()
