"""
Module amr_profiles


Defined at profiles_module.fpp lines 5-367

"""
from __future__ import print_function, absolute_import, division
import _amr2_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("amr2_pkg.profile_handler")
class profile_handler(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=profile_handler)
    
    
    Defined at profiles_module.fpp lines 9-18
    
    """
    def __init__(self, handle=None):
        """
        self = Profile_Handler()
        
        
        Defined at profiles_module.fpp lines 9-18
        
        
        Returns
        -------
        this : Profile_Handler
        	Object to be constructed
        
        
        Automatically generated constructor for profile_handler
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _amr2_pkg.f90wrap_profile_handler_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Profile_Handler
        
        
        Defined at profiles_module.fpp lines 9-18
        
        Parameters
        ----------
        this : Profile_Handler
        	Object to be destructed
        
        
        Automatically generated destructor for profile_handler
        """
        if self._alloc:
            _amr2_pkg.f90wrap_profile_handler_finalise(this=self._handle)
    
    @property
    def profdim(self):
        """
        Element profdim ftype=integer  pytype=int
        
        
        Defined at profiles_module.fpp line 10
        
        """
        return _amr2_pkg.f90wrap_profile_handler__get__profdim(self._handle)
    
    @profdim.setter
    def profdim(self, profdim):
        _amr2_pkg.f90wrap_profile_handler__set__profdim(self._handle, profdim)
    
    @property
    def xvarname(self):
        """
        Element xvarname ftype=character(128) pytype=str
        
        
        Defined at profiles_module.fpp line 11
        
        """
        return _amr2_pkg.f90wrap_profile_handler__get__xvarname(self._handle)
    
    @xvarname.setter
    def xvarname(self, xvarname):
        _amr2_pkg.f90wrap_profile_handler__set__xvarname(self._handle, xvarname)
    
    @property
    def nyvar(self):
        """
        Element nyvar ftype=integer  pytype=int
        
        
        Defined at profiles_module.fpp line 12
        
        """
        return _amr2_pkg.f90wrap_profile_handler__get__nyvar(self._handle)
    
    @nyvar.setter
    def nyvar(self, nyvar):
        _amr2_pkg.f90wrap_profile_handler__set__nyvar(self._handle, nyvar)
    
    @property
    def yvarnames(self):
        """
        Element yvarnames ftype=character(128) pytype=str
        
        
        Defined at profiles_module.fpp line 13
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_profile_handler__array__yvarnames(self._handle)
        if array_handle in self._arrays:
            yvarnames = self._arrays[array_handle]
        else:
            yvarnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_profile_handler__array__yvarnames)
            self._arrays[array_handle] = yvarnames
        return yvarnames
    
    @yvarnames.setter
    def yvarnames(self, yvarnames):
        self.yvarnames[...] = yvarnames
    
    @property
    def nbins(self):
        """
        Element nbins ftype=integer  pytype=int
        
        
        Defined at profiles_module.fpp line 14
        
        """
        return _amr2_pkg.f90wrap_profile_handler__get__nbins(self._handle)
    
    @nbins.setter
    def nbins(self, nbins):
        _amr2_pkg.f90wrap_profile_handler__set__nbins(self._handle, nbins)
    
    @property
    def nwvar(self):
        """
        Element nwvar ftype=integer  pytype=int
        
        
        Defined at profiles_module.fpp line 15
        
        """
        return _amr2_pkg.f90wrap_profile_handler__get__nwvar(self._handle)
    
    @nwvar.setter
    def nwvar(self, nwvar):
        _amr2_pkg.f90wrap_profile_handler__set__nwvar(self._handle, nwvar)
    
    @property
    def wvarnames(self):
        """
        Element wvarnames ftype=character(128) pytype=str
        
        
        Defined at profiles_module.fpp line 16
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_profile_handler__array__wvarnames(self._handle)
        if array_handle in self._arrays:
            wvarnames = self._arrays[array_handle]
        else:
            wvarnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_profile_handler__array__wvarnames)
            self._arrays[array_handle] = wvarnames
        return wvarnames
    
    @wvarnames.setter
    def wvarnames(self, wvarnames):
        self.wvarnames[...] = wvarnames
    
    @property
    def xdata(self):
        """
        Element xdata ftype=real(dbl) pytype=float
        
        
        Defined at profiles_module.fpp line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_profile_handler__array__xdata(self._handle)
        if array_handle in self._arrays:
            xdata = self._arrays[array_handle]
        else:
            xdata = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_profile_handler__array__xdata)
            self._arrays[array_handle] = xdata
        return xdata
    
    @xdata.setter
    def xdata(self, xdata):
        self.xdata[...] = xdata
    
    @property
    def ydata(self):
        """
        Element ydata ftype=real(dbl) pytype=float
        
        
        Defined at profiles_module.fpp line 18
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_profile_handler__array__ydata(self._handle)
        if array_handle in self._arrays:
            ydata = self._arrays[array_handle]
        else:
            ydata = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_profile_handler__array__ydata)
            self._arrays[array_handle] = ydata
        return ydata
    
    @ydata.setter
    def ydata(self, ydata):
        self.ydata[...] = ydata
    
    def __str__(self):
        ret = ['<profile_handler>{\n']
        ret.append('    profdim : ')
        ret.append(repr(self.profdim))
        ret.append(',\n    xvarname : ')
        ret.append(repr(self.xvarname))
        ret.append(',\n    nyvar : ')
        ret.append(repr(self.nyvar))
        ret.append(',\n    yvarnames : ')
        ret.append(repr(self.yvarnames))
        ret.append(',\n    nbins : ')
        ret.append(repr(self.nbins))
        ret.append(',\n    nwvar : ')
        ret.append(repr(self.nwvar))
        ret.append(',\n    wvarnames : ')
        ret.append(repr(self.wvarnames))
        ret.append(',\n    xdata : ')
        ret.append(repr(self.xdata))
        ret.append(',\n    ydata : ')
        ret.append(repr(self.ydata))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def allocate_profile_handler(self):
    """
    allocate_profile_handler(self)
    
    
    Defined at profiles_module.fpp lines 21-27
    
    Parameters
    ----------
    prof : Profile_Handler
    
    """
    _amr2_pkg.f90wrap_allocate_profile_handler(prof=self._handle)

def makebins(self, varname, nbins, bins):
    """
    makebins(self, varname, nbins, bins)
    
    
    Defined at profiles_module.fpp lines 29-43
    
    Parameters
    ----------
    reg : Region
    varname : str
    nbins : int
    bins : float array
    
    """
    _amr2_pkg.f90wrap_makebins(reg=self._handle, varname=varname, nbins=nbins, \
        bins=bins)

def findbinpos(self, varids, distance, pos, cellvars, cellsize, prof, ibin):
    """
    findbinpos(self, varids, distance, pos, cellvars, cellsize, prof, ibin)
    
    
    Defined at profiles_module.fpp lines 45-65
    
    Parameters
    ----------
    reg : Region
    varids : Hydroid
    distance : float
    pos : float array
    cellvars : float array
    cellsize : float
    prof : Profile_Handler
    ibin : int
    
    """
    _amr2_pkg.f90wrap_findbinpos(reg=self._handle, varids=varids._handle, \
        distance=distance, pos=pos, cellvars=cellvars, cellsize=cellsize, \
        prof=prof._handle, ibin=ibin)

def bindata(self, varids, pos, cellvars, cellsize, prof, ibin):
    """
    bindata(self, varids, pos, cellvars, cellsize, prof, ibin)
    
    
    Defined at profiles_module.fpp lines 67-114
    
    Parameters
    ----------
    reg : Region
    varids : Hydroid
    pos : float array
    cellvars : float array
    cellsize : float
    prof : Profile_Handler
    ibin : int
    
    """
    _amr2_pkg.f90wrap_bindata(reg=self._handle, varids=varids._handle, pos=pos, \
        cellvars=cellvars, cellsize=cellsize, prof=prof._handle, ibin=ibin)

def renormalise_bins(self):
    """
    renormalise_bins(self)
    
    
    Defined at profiles_module.fpp lines 116-141
    
    Parameters
    ----------
    prof_data : Profile_Handler
    
    """
    _amr2_pkg.f90wrap_renormalise_bins(prof_data=self._handle)

def get_cells_onedprofile(repository, amr, reg, filt, varids, prof_data):
    """
    get_cells_onedprofile(repository, amr, reg, filt, varids, prof_data)
    
    
    Defined at profiles_module.fpp lines 143-344
    
    Parameters
    ----------
    repository : str
    amr : Amr_Info
    reg : Region
    filt : Filter
    varids : Hydroid
    prof_data : Profile_Handler
    
    """
    _amr2_pkg.f90wrap_get_cells_onedprofile(repository=repository, amr=amr._handle, \
        reg=reg._handle, filt=filt._handle, varids=varids._handle, \
        prof_data=prof_data._handle)

def onedprofile(repository, reg, filt, prof_data, lmax):
    """
    onedprofile(repository, reg, filt, prof_data, lmax)
    
    
    Defined at profiles_module.fpp lines 346-367
    
    Parameters
    ----------
    repository : str
    reg : Region
    filt : Filter
    prof_data : Profile_Handler
    lmax : int
    
    """
    _amr2_pkg.f90wrap_onedprofile(repository=repository, reg=reg._handle, \
        filt=filt._handle, prof_data=prof_data._handle, lmax=lmax)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "amr_profiles".')

for func in _dt_array_initialisers:
    func()