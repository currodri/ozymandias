"""
Module coordinate_systems


Defined at coordinates_module.fpp lines 23-173

"""
from __future__ import print_function, absolute_import, division
import _projections_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def r_sphere(self):
    """
    r_sphere = r_sphere(self)
    
    
    Defined at coordinates_module.fpp lines 31-33
    
    Parameters
    ----------
    p : Vector
    
    Returns
    -------
    r_sphere : float
    
    """
    r_sphere = _projections_pkg.f90wrap_r_sphere(p=self._handle)
    return r_sphere

def theta_sphere(self):
    """
    theta_sphere = theta_sphere(self)
    
    
    Defined at coordinates_module.fpp lines 35-39
    
    Parameters
    ----------
    p : Vector
    
    Returns
    -------
    theta_sphere : float
    
    """
    theta_sphere = _projections_pkg.f90wrap_theta_sphere(p=self._handle)
    return theta_sphere

def phi_sphere(self):
    """
    phi_sphere = phi_sphere(self)
    
    
    Defined at coordinates_module.fpp lines 41-43
    
    Parameters
    ----------
    p : Vector
    
    Returns
    -------
    phi_sphere : float
    
    """
    phi_sphere = _projections_pkg.f90wrap_phi_sphere(p=self._handle)
    return phi_sphere

def r_cyl(self):
    """
    r_cyl = r_cyl(self)
    
    
    Defined at coordinates_module.fpp lines 45-47
    
    Parameters
    ----------
    p : Vector
    
    Returns
    -------
    r_cyl : float
    
    """
    r_cyl = _projections_pkg.f90wrap_r_cyl(p=self._handle)
    return r_cyl

def phi_cyl(self):
    """
    phi_cyl = phi_cyl(self)
    
    
    Defined at coordinates_module.fpp lines 49-61
    
    Parameters
    ----------
    p : Vector
    
    Returns
    -------
    phi_cyl : float
    
    """
    phi_cyl = _projections_pkg.f90wrap_phi_cyl(p=self._handle)
    return phi_cyl

def spherical_basis_from_cartesian(self, spher_basis):
    """
    spherical_basis_from_cartesian(self, spher_basis)
    
    
    Defined at coordinates_module.fpp lines 69-95
    
    Parameters
    ----------
    p : Vector
    spher_basis : Basis
    
    """
    _projections_pkg.f90wrap_spherical_basis_from_cartesian(p=self._handle, \
        spher_basis=spher_basis._handle)

def cylindrical_basis_from_cartesian(self, cyl_basis):
    """
    cylindrical_basis_from_cartesian(self, cyl_basis)
    
    
    Defined at coordinates_module.fpp lines 103-111
    
    Parameters
    ----------
    p : Vector
    cyl_basis : Basis
    
    """
    _projections_pkg.f90wrap_cylindrical_basis_from_cartesian(p=self._handle, \
        cyl_basis=cyl_basis._handle)

def new_z_coordinates(self, transformation_matrix, errormsg):
    """
    new_z_coordinates(self, transformation_matrix, errormsg)
    
    
    Defined at coordinates_module.fpp lines 123-172
    
    Parameters
    ----------
    axis : Vector
    transformation_matrix : float array
    errormsg : int
    
    """
    _projections_pkg.f90wrap_new_z_coordinates(axis=self._handle, \
        transformation_matrix=transformation_matrix, errormsg=errormsg)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "coordinate_systems".')

for func in _dt_array_initialisers:
    func()
