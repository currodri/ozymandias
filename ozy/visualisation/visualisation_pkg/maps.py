"""
Module maps


Defined at ramses2map.fpp lines 209-1733

"""
from __future__ import print_function, absolute_import, division
import _visualisation_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("visualisation_pkg.projection_handler")
class projection_handler(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=projection_handler)
    
    
    Defined at ramses2map.fpp lines 217-222
    
    """
    def __init__(self, handle=None):
        """
        self = Projection_Handler()
        
        
        Defined at ramses2map.fpp lines 217-222
        
        
        Returns
        -------
        this : Projection_Handler
        	Object to be constructed
        
        
        Automatically generated constructor for projection_handler
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _visualisation_pkg.f90wrap_projection_handler_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Projection_Handler
        
        
        Defined at ramses2map.fpp lines 217-222
        
        Parameters
        ----------
        this : Projection_Handler
        	Object to be destructed
        
        
        Automatically generated destructor for projection_handler
        """
        if self._alloc:
            _visualisation_pkg.f90wrap_projection_handler_finalise(this=self._handle)
    
    @property
    def pov(self):
        """
        Element pov ftype=character(128) pytype=str
        
        
        Defined at ramses2map.fpp line 218
        
        """
        return _visualisation_pkg.f90wrap_projection_handler__get__pov(self._handle)
    
    @pov.setter
    def pov(self, pov):
        _visualisation_pkg.f90wrap_projection_handler__set__pov(self._handle, pov)
    
    @property
    def nvars(self):
        """
        Element nvars ftype=integer  pytype=int
        
        
        Defined at ramses2map.fpp line 219
        
        """
        return _visualisation_pkg.f90wrap_projection_handler__get__nvars(self._handle)
    
    @nvars.setter
    def nvars(self, nvars):
        _visualisation_pkg.f90wrap_projection_handler__set__nvars(self._handle, nvars)
    
    @property
    def nfilter(self):
        """
        Element nfilter ftype=integer  pytype=int
        
        
        Defined at ramses2map.fpp line 219
        
        """
        return _visualisation_pkg.f90wrap_projection_handler__get__nfilter(self._handle)
    
    @nfilter.setter
    def nfilter(self, nfilter):
        _visualisation_pkg.f90wrap_projection_handler__set__nfilter(self._handle, \
            nfilter)
    
    @property
    def varnames(self):
        """
        Element varnames ftype=character(128) pytype=str
        
        
        Defined at ramses2map.fpp line 220
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_projection_handler__array__varnames(self._handle)
        if array_handle in self._arrays:
            varnames = self._arrays[array_handle]
        else:
            varnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_projection_handler__array__varnames)
            self._arrays[array_handle] = varnames
        return varnames
    
    @varnames.setter
    def varnames(self, varnames):
        self.varnames[...] = varnames
    
    @property
    def weightvar(self):
        """
        Element weightvar ftype=character(128) pytype=str
        
        
        Defined at ramses2map.fpp line 221
        
        """
        return \
            _visualisation_pkg.f90wrap_projection_handler__get__weightvar(self._handle)
    
    @weightvar.setter
    def weightvar(self, weightvar):
        _visualisation_pkg.f90wrap_projection_handler__set__weightvar(self._handle, \
            weightvar)
    
    @property
    def toto(self):
        """
        Element toto ftype=real(dbl) pytype=float
        
        
        Defined at ramses2map.fpp line 222
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_projection_handler__array__toto(self._handle)
        if array_handle in self._arrays:
            toto = self._arrays[array_handle]
        else:
            toto = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_projection_handler__array__toto)
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
        ret.append(',\n    nfilter : ')
        ret.append(repr(self.nfilter))
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
    
    
    Defined at ramses2map.fpp lines 225-228
    
    Parameters
    ----------
    proj : Projection_Handler
    
    """
    _visualisation_pkg.f90wrap_allocate_projection_handler(proj=self._handle)

def projection_hydro(repository, cam, bulk_velocity, use_neigh, proj, lmax=None, \
    lmin=None):
    """
    projection_hydro(repository, cam, bulk_velocity, use_neigh, proj[, lmax, lmin])
    
    
    Defined at ramses2map.fpp lines 230-1131
    
    Parameters
    ----------
    repository : str
    cam : Camera
    bulk_velocity : Vector
    use_neigh : bool
    proj : Projection_Handler
    lmax : int
    lmin : int
    
    """
    _visualisation_pkg.f90wrap_projection_hydro(repository=repository, \
        cam=cam._handle, bulk_velocity=bulk_velocity._handle, use_neigh=use_neigh, \
        proj=proj._handle, lmax=lmax, lmin=lmin)

def projection_parts(repository, cam, bulk_velocity, proj, tag_file=None, \
    inverse_tag=None):
    """
    projection_parts(repository, cam, bulk_velocity, proj[, tag_file, inverse_tag])
    
    
    Defined at ramses2map.fpp lines 1133-1163
    
    Parameters
    ----------
    repository : str
    cam : Camera
    bulk_velocity : Vector
    proj : Projection_Handler
    tag_file : str
    inverse_tag : bool
    
    """
    _visualisation_pkg.f90wrap_projection_parts(repository=repository, \
        cam=cam._handle, bulk_velocity=bulk_velocity._handle, proj=proj._handle, \
        tag_file=tag_file, inverse_tag=inverse_tag)

def project_particles(repository, bbox, cam, proj, tag_file=None, \
    inverse_tag=None):
    """
    project_particles(repository, bbox, cam, proj[, tag_file, inverse_tag])
    
    
    Defined at ramses2map.fpp lines 1165-1429
    
    Parameters
    ----------
    repository : str
    bbox : Region
    cam : Camera
    proj : Projection_Handler
    tag_file : str
    inverse_tag : bool
    
    """
    _visualisation_pkg.f90wrap_project_particles(repository=repository, \
        bbox=bbox._handle, cam=cam._handle, proj=proj._handle, tag_file=tag_file, \
        inverse_tag=inverse_tag)

def healpix_hydro(repository, reg, nside, proj):
    """
    healpix_hydro(repository, reg, nside, proj)
    
    
    Defined at ramses2map.fpp lines 1431-1443
    
    Parameters
    ----------
    repository : str
    reg : Region
    nside : int
    proj : Projection_Handler
    
    """
    _visualisation_pkg.f90wrap_healpix_hydro(repository=repository, reg=reg._handle, \
        nside=nside, proj=proj._handle)

def project_cells_hpix(repository, reg, nside, proj):
    """
    project_cells_hpix(repository, reg, nside, proj)
    
    
    Defined at ramses2map.fpp lines 1445-1733
    
    Parameters
    ----------
    repository : str
    reg : Region
    nside : int
    proj : Projection_Handler
    
    """
    _visualisation_pkg.f90wrap_project_cells_hpix(repository=repository, \
        reg=reg._handle, nside=nside, proj=proj._handle)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "maps".')

for func in _dt_array_initialisers:
    func()
