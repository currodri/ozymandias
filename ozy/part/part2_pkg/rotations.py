"""
Module rotations


Defined at linalg_module.fpp lines 176-249

"""
from __future__ import print_function, absolute_import, division
import _part2_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def euler_matrix(r, dim, angle, ap):
    """
    euler_matrix(r, dim, angle, ap)
    
    
    Defined at linalg_module.fpp lines 193-226
    
    Parameters
    ----------
    r : float array
    dim : int
    angle : float
    ap : str
    
    """
    _part2_pkg.f90wrap_euler_matrix(r=r, dim=dim, angle=angle, ap=ap)

def _rotate_vector_single(self, rotation_matrix):
    """
    _rotate_vector_single(self, rotation_matrix)
    
    
    Defined at linalg_module.fpp lines 237-240
    
    Parameters
    ----------
    vec : Vector
    rotation_matrix : float array
    
    """
    _part2_pkg.f90wrap_rotate_vector_single(vec=self._handle, \
        rotation_matrix=rotation_matrix)

def _rotate_vector_array(self, rotation_matrix):
    """
    _rotate_vector_array(self, rotation_matrix)
    
    
    Defined at linalg_module.fpp lines 242-248
    
    Parameters
    ----------
    vec : Array_Vectors
    rotation_matrix : float array
    
    """
    _part2_pkg.f90wrap_rotate_vector_array(vec=self._handle, \
        rotation_matrix=rotation_matrix)

def rotate_vector(*args, **kwargs):
    """
    rotate_vector(*args, **kwargs)
    
    
    Defined at linalg_module.fpp lines 180-182
    
    Overloaded interface containing the following procedures:
      _rotate_vector_single
      _rotate_vector_array
    
    """
    for proc in [_rotate_vector_single, _rotate_vector_array]:
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
    logging.debug('unallocated array(s) detected on import of module "rotations".')

for func in _dt_array_initialisers:
    func()
