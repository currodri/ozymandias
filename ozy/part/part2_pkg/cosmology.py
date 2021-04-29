"""
Module cosmology


Defined at cosmology_module.fpp lines 5-138

"""
from __future__ import print_function, absolute_import, division
import _part2_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def cosmology_model(self):
    """
    cosmology_model(self)
    
    
    Defined at cosmology_module.fpp lines 16-42
    
    Parameters
    ----------
    sim : Sim_Info
    
    """
    _part2_pkg.f90wrap_cosmology_model(sim=self._handle)

def friedmann(o_mat_0, o_vac_0, o_k_0, alpha, axp_min, axp_out, hexp_out, \
    tau_out, t_out, ntable, age_tot):
    """
    friedmann(o_mat_0, o_vac_0, o_k_0, alpha, axp_min, axp_out, hexp_out, tau_out, \
        t_out, ntable, age_tot)
    
    
    Defined at cosmology_module.fpp lines 57-120
    
    Parameters
    ----------
    o_mat_0 : float
    o_vac_0 : float
    o_k_0 : float
    alpha : float
    axp_min : float
    axp_out : float array
    hexp_out : float array
    tau_out : float array
    t_out : float array
    ntable : int
    age_tot : float
    
    """
    _part2_pkg.f90wrap_friedmann(o_mat_0=o_mat_0, o_vac_0=o_vac_0, o_k_0=o_k_0, \
        alpha=alpha, axp_min=axp_min, axp_out=axp_out, hexp_out=hexp_out, \
        tau_out=tau_out, t_out=t_out, ntable=ntable, age_tot=age_tot)

def dadtau(axp_tau, o_mat_0, o_vac_0, o_k_0):
    """
    dadtau = dadtau(axp_tau, o_mat_0, o_vac_0, o_k_0)
    
    
    Defined at cosmology_module.fpp lines 122-129
    
    Parameters
    ----------
    axp_tau : float
    o_mat_0 : float
    o_vac_0 : float
    o_k_0 : float
    
    Returns
    -------
    dadtau : float
    
    """
    dadtau = _part2_pkg.f90wrap_dadtau(axp_tau=axp_tau, o_mat_0=o_mat_0, \
        o_vac_0=o_vac_0, o_k_0=o_k_0)
    return dadtau

def dadt(axp_t, o_mat_0, o_vac_0, o_k_0):
    """
    dadt = dadt(axp_t, o_mat_0, o_vac_0, o_k_0)
    
    
    Defined at cosmology_module.fpp lines 131-138
    
    Parameters
    ----------
    axp_t : float
    o_mat_0 : float
    o_vac_0 : float
    o_k_0 : float
    
    Returns
    -------
    dadt : float
    
    """
    dadt = _part2_pkg.f90wrap_dadt(axp_t=axp_t, o_mat_0=o_mat_0, o_vac_0=o_vac_0, \
        o_k_0=o_k_0)
    return dadt


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "cosmology".')

for func in _dt_array_initialisers:
    func()
