"""
Module amr_map


Defined at amr2map.fpp lines 208-554

"""
from __future__ import print_function, absolute_import, division
import _projections_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def projection(repository, cam, bulk_velocity):
    """
    projection(repository, cam, bulk_velocity)
    
    
    Defined at amr2map.fpp lines 216-251
    
    Parameters
    ----------
    repository : str
    cam : Camera
    bulk_velocity : Vector
    
    """
    _projections_pkg.f90wrap_projection(repository=repository, cam=cam._handle, \
        bulk_velocity=bulk_velocity._handle)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "amr_map".')

for func in _dt_array_initialisers:
    func()
