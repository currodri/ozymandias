from __future__ import print_function, absolute_import, division
import _hutils
import f90wrap.runtime
import logging
import numpy

class Py_Halo_Utils(f90wrap.runtime.FortranModule):
    """
    Module py_halo_utils
    
    
    Defined at py_halo_utils.fpp lines 5-251
    
    """
    @staticmethod
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
        nhalo, nsubhalo = _hutils.f90wrap_py_halo_utils__get_nb_halos(file=file)
        return nhalo, nsubhalo
    
    @staticmethod
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
        _hutils.f90wrap_py_halo_utils__read_all_halos_with_contam(file=file, \
            halos=halos, n=n)
    
    @staticmethod
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
        _hutils.f90wrap_py_halo_utils__read_all_halos(file=file, halos=halos, n=n)
    
    @staticmethod
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
        _hutils.f90wrap_py_halo_utils__read_all_galaxies(file=file, halos=halos, n=n)
    
    @staticmethod
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
        _hutils.f90wrap_py_halo_utils__get_sfrs(starfile=starfile, \
            lookback_myr=lookback_myr, ngals=ngals, sfrs=sfrs)
    
    @staticmethod
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
        _hutils.f90wrap_py_halo_utils__get_ages_metallicities(starfile=starfile, \
            ngals=ngals, ages_myr=ages_myr, zs=zs)
    
    _dt_array_initialisers = []
    

py_halo_utils = Py_Halo_Utils()

