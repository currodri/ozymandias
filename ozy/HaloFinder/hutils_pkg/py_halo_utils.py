"""
Module py_halo_utils


Defined at py_halo_utils.fpp lines 5-251

"""
from __future__ import print_function, absolute_import, division
import _hutils_pkg
import f90wrap.runtime
import logging
import numpy

_arrays = {}
_objs = {}

def get_nb_halos(file):
    """
    nhalo, nsubhalo = get_nb_halos(file)
    
    
    Defined at py_halo_utils.fpp lines 9-22
    
    Parameters
    ----------
    file : str
    
    Returns
    -------
    nhalo : int
    nsubhalo : int
    
    """
    nhalo, nsubhalo = _hutils_pkg.f90wrap_py_halo_utils__get_nb_halos(file=file)
    return nhalo, nsubhalo

def read_all_halos_with_contam(file, halos, n):
    """
    read_all_halos_with_contam(file, halos, n)
    
    
    Defined at py_halo_utils.fpp lines 25-62
    
    Parameters
    ----------
    file : str
    halos : float array
    n : int
    
    """
    _hutils_pkg.f90wrap_py_halo_utils__read_all_halos_with_contam(file=file, \
        halos=halos, n=n)

def read_all_halos(file, halos, n):
    """
    read_all_halos(file, halos, n)
    
    
    Defined at py_halo_utils.fpp lines 65-101
    
    Parameters
    ----------
    file : str
    halos : float array
    n : int
    
    """
    _hutils_pkg.f90wrap_py_halo_utils__read_all_halos(file=file, halos=halos, n=n)

def read_all_galaxies(file, halos, n):
    """
    read_all_galaxies(file, halos, n)
    
    
    Defined at py_halo_utils.fpp lines 104-147
    
    Parameters
    ----------
    file : str
    halos : float array
    n : int
    
    """
    _hutils_pkg.f90wrap_py_halo_utils__read_all_galaxies(file=file, halos=halos, \
        n=n)

def get_sfrs(starfile, lookback_myr, ngals, sfrs):
    """
    get_sfrs(starfile, lookback_myr, ngals, sfrs)
    
    
    Defined at py_halo_utils.fpp lines 150-199
    
    Parameters
    ----------
    starfile : str
    lookback_myr : float
    ngals : int
    sfrs : float array
    
    -----------------------------------------------------------------------
    """
    _hutils_pkg.f90wrap_py_halo_utils__get_sfrs(starfile=starfile, \
        lookback_myr=lookback_myr, ngals=ngals, sfrs=sfrs)

def get_ages_metallicities(starfile, ngals, ages_myr, zs):
    """
    get_ages_metallicities(starfile, ngals, ages_myr, zs)
    
    
    Defined at py_halo_utils.fpp lines 202-251
    
    Parameters
    ----------
    starfile : str
    ngals : int
    ages_myr : float array
    zs : float array
    
    -----------------------------------------------------------------------
    """
    _hutils_pkg.f90wrap_py_halo_utils__get_ages_metallicities(starfile=starfile, \
        ngals=ngals, ages_myr=ages_myr, zs=zs)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "py_halo_utils".')

for func in _dt_array_initialisers:
    func()
