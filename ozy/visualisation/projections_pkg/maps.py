"""
Module maps


Defined at ramses2map.fpp lines 194-1011

"""
from __future__ import print_function, absolute_import, division
import _projections_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("projections_pkg.projection_handler")
class projection_handler(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=projection_handler)
    
    
    Defined at ramses2map.fpp lines 201-206
    
    """
    def __init__(self, handle=None):
        """
        self = Projection_Handler()
        
        
        Defined at ramses2map.fpp lines 201-206
        
        
        Returns
        -------
        this : Projection_Handler
        	Object to be constructed
        
        
        Automatically generated constructor for projection_handler
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _projections_pkg.f90wrap_projection_handler_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Projection_Handler
        
        
        Defined at ramses2map.fpp lines 201-206
        
        Parameters
        ----------
        this : Projection_Handler
        	Object to be destructed
        
        
        Automatically generated destructor for projection_handler
        """
        if self._alloc:
            _projections_pkg.f90wrap_projection_handler_finalise(this=self._handle)
    
    @property
    def pov(self):
        """
        Element pov ftype=character(128) pytype=str
        
        
        Defined at ramses2map.fpp line 202
        
        """
        return _projections_pkg.f90wrap_projection_handler__get__pov(self._handle)
    
    @pov.setter
    def pov(self, pov):
        _projections_pkg.f90wrap_projection_handler__set__pov(self._handle, pov)
    
    @property
    def nvars(self):
        """
        Element nvars ftype=integer  pytype=int
        
        
        Defined at ramses2map.fpp line 203
        
        """
        return _projections_pkg.f90wrap_projection_handler__get__nvars(self._handle)
    
    @nvars.setter
    def nvars(self, nvars):
        _projections_pkg.f90wrap_projection_handler__set__nvars(self._handle, nvars)
    
    @property
    def varnames(self):
        """
        Element varnames ftype=character(128) pytype=str
        
        
        Defined at ramses2map.fpp line 204
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _projections_pkg.f90wrap_projection_handler__array__varnames(self._handle)
        if array_handle in self._arrays:
            varnames = self._arrays[array_handle]
        else:
            varnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _projections_pkg.f90wrap_projection_handler__array__varnames)
            self._arrays[array_handle] = varnames
        return varnames
    
    @varnames.setter
    def varnames(self, varnames):
        self.varnames[...] = varnames
    
    @property
    def weightvar(self):
        """
        Element weightvar ftype=character(128) pytype=str
        
        
        Defined at ramses2map.fpp line 205
        
        """
        return _projections_pkg.f90wrap_projection_handler__get__weightvar(self._handle)
    
    @weightvar.setter
    def weightvar(self, weightvar):
        _projections_pkg.f90wrap_projection_handler__set__weightvar(self._handle, \
            weightvar)
    
    @property
    def toto(self):
        """
        Element toto ftype=real(dbl) pytype=float
        
        
        Defined at ramses2map.fpp line 206
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _projections_pkg.f90wrap_projection_handler__array__toto(self._handle)
        if array_handle in self._arrays:
            toto = self._arrays[array_handle]
        else:
            toto = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _projections_pkg.f90wrap_projection_handler__array__toto)
            self._arrays[array_handle] = toto
        return toto
    
    @toto.setter
    def toto(self, toto):
        self.toto[...] = toto
    
    def __str__(self):
        ret = ['<projection_handler>{\n']
        ret.append('    pov : ')
        ret.append(repr(self.pov))
        ret.append(',\n    nvars : ')
        ret.append(repr(self.nvars))
        ret.append(',\n    varnames : ')
        ret.append(repr(self.varnames))
        ret.append(',\n    weightvar : ')
        ret.append(repr(self.weightvar))
        ret.append(',\n    toto : ')
        ret.append(repr(self.toto))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def allocate_projection_handler(self):
    """
    allocate_projection_handler(self)
    
    
    Defined at ramses2map.fpp lines 209-212
    
    Parameters
    ----------
    proj : Projection_Handler
    
    """
    _projections_pkg.f90wrap_allocate_projection_handler(proj=self._handle)

def projection_hydro(repository, cam, bulk_velocity, proj):
    """
    projection_hydro(repository, cam, bulk_velocity, proj)
    
    
    Defined at ramses2map.fpp lines 214-242
    
    Parameters
    ----------
    repository : str
    cam : Camera
    bulk_velocity : Vector
    proj : Projection_Handler
    
    """
    _projections_pkg.f90wrap_projection_hydro(repository=repository, \
        cam=cam._handle, bulk_velocity=bulk_velocity._handle, proj=proj._handle)

def project_cells(repository, amr, bbox, varids, cam, proj):
    """
    project_cells(repository, amr, bbox, varids, cam, proj)
    
    
    Defined at ramses2map.fpp lines 244-539
    
    Parameters
    ----------
    repository : str
    amr : Amr_Info
    bbox : Region
    varids : Hydroid
    cam : Camera
    proj : Projection_Handler
    
    """
    _projections_pkg.f90wrap_project_cells(repository=repository, amr=amr._handle, \
        bbox=bbox._handle, varids=varids._handle, cam=cam._handle, \
        proj=proj._handle)

def projection_parts(repository, cam, bulk_velocity, proj):
    """
    projection_parts(repository, cam, bulk_velocity, proj)
    
    
    Defined at ramses2map.fpp lines 541-561
    
    Parameters
    ----------
    repository : str
    cam : Camera
    bulk_velocity : Vector
    proj : Projection_Handler
    
    """
    _projections_pkg.f90wrap_projection_parts(repository=repository, \
        cam=cam._handle, bulk_velocity=bulk_velocity._handle, proj=proj._handle)

def project_particles(repository, amr, sim, bbox, cam, proj):
    """
    project_particles(repository, amr, sim, bbox, cam, proj)
    
    
    Defined at ramses2map.fpp lines 563-740
    
    Parameters
    ----------
    repository : str
    amr : Amr_Info
    sim : Sim_Info
    bbox : Region
    cam : Camera
    proj : Projection_Handler
    
    """
    _projections_pkg.f90wrap_project_particles(repository=repository, \
        amr=amr._handle, sim=sim._handle, bbox=bbox._handle, cam=cam._handle, \
        proj=proj._handle)

def healpix_hydro(repository, reg, nside, proj):
    """
    healpix_hydro(repository, reg, nside, proj)
    
    
    Defined at ramses2map.fpp lines 742-757
    
    Parameters
    ----------
    repository : str
    reg : Region
    nside : int
    proj : Projection_Handler
    
    """
    _projections_pkg.f90wrap_healpix_hydro(repository=repository, reg=reg._handle, \
        nside=nside, proj=proj._handle)

def project_cells_hpix(repository, amr, reg, varids, nside, proj):
    """
    project_cells_hpix(repository, amr, reg, varids, nside, proj)
    
    
    Defined at ramses2map.fpp lines 759-1011
    
    Parameters
    ----------
    repository : str
    amr : Amr_Info
    reg : Region
    varids : Hydroid
    nside : int
    proj : Projection_Handler
    
    """
    _projections_pkg.f90wrap_project_cells_hpix(repository=repository, \
        amr=amr._handle, reg=reg._handle, varids=varids._handle, nside=nside, \
        proj=proj._handle)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "maps".')

for func in _dt_array_initialisers:
    func()
