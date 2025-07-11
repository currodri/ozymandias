"""
Module py_tree_ios


Defined at py_tree_ios.fpp lines 5-109

"""
from __future__ import print_function, absolute_import, division
import _tutils_pkg
import f90wrap.runtime
import logging
import numpy

_arrays = {}
_objs = {}

def get_nsteps(treedir, treenum):
    """
    nsteps = get_nsteps(treedir, treenum)
    
    
    Defined at py_tree_ios.fpp lines 8-18
    
    Parameters
    ----------
    treedir : str
    treenum : int
    
    Returns
    -------
    nsteps : int
    
    """
    nsteps = _tutils_pkg.f90wrap_py_tree_ios__get_nsteps(treedir=treedir, \
        treenum=treenum)
    return nsteps

def get_nids(treedir, treenum):
    """
    nids = get_nids(treedir, treenum)
    
    
    Defined at py_tree_ios.fpp lines 20-31
    
    Parameters
    ----------
    treedir : str
    treenum : int
    
    Returns
    -------
    nids : int
    
    """
    nids = _tutils_pkg.f90wrap_py_tree_ios__get_nids(treedir=treedir, \
        treenum=treenum)
    return nids

def get_nprops(treedir, treenum):
    """
    nprops = get_nprops(treedir, treenum)
    
    
    Defined at py_tree_ios.fpp lines 33-44
    
    Parameters
    ----------
    treedir : str
    treenum : int
    
    Returns
    -------
    nprops : int
    
    """
    nprops = _tutils_pkg.f90wrap_py_tree_ios__get_nprops(treedir=treedir, \
        treenum=treenum)
    return nprops

def read_timestep_props(treedir, treenum, nsteps, nhalos, aexp, age_univ):
    """
    read_timestep_props(treedir, treenum, nsteps, nhalos, aexp, age_univ)
    
    
    Defined at py_tree_ios.fpp lines 46-61
    
    Parameters
    ----------
    treedir : str
    treenum : int
    nsteps : int
    nhalos : int array
    aexp : float array
    age_univ : float array
    
    """
    _tutils_pkg.f90wrap_py_tree_ios__read_timestep_props(treedir=treedir, \
        treenum=treenum, nsteps=nsteps, nhalos=nhalos, aexp=aexp, age_univ=age_univ)

def read_tree_struct(treedir, treenum, nsteps, nhalos, nh, nids, ids):
    """
    read_tree_struct(treedir, treenum, nsteps, nhalos, nh, nids, ids)
    
    
    Defined at py_tree_ios.fpp lines 63-85
    
    Parameters
    ----------
    treedir : str
    treenum : int
    nsteps : int
    nhalos : int array
    nh : int
    nids : int
    ids : int array
    
    """
    _tutils_pkg.f90wrap_py_tree_ios__read_tree_struct(treedir=treedir, \
        treenum=treenum, nsteps=nsteps, nhalos=nhalos, nh=nh, nids=nids, ids=ids)

def read_props(treedir, treenum, nsteps, nhalos, nh, nprops, props):
    """
    read_props(treedir, treenum, nsteps, nhalos, nh, nprops, props)
    
    
    Defined at py_tree_ios.fpp lines 87-108
    
    Parameters
    ----------
    treedir : str
    treenum : int
    nsteps : int
    nhalos : int array
    nh : int
    nprops : int
    props : float array
    
    """
    _tutils_pkg.f90wrap_py_tree_ios__read_props(treedir=treedir, treenum=treenum, \
        nsteps=nsteps, nhalos=nhalos, nh=nh, nprops=nprops, props=props)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "py_tree_ios".')

for func in _dt_array_initialisers:
    func()
