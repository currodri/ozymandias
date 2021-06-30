"""
Module amr_profiles


Defined at profiles_module.fpp lines 25-770

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
    
    
    Defined at profiles_module.fpp lines 29-38
    
    """
    def __init__(self, handle=None):
        """
        self = Profile_Handler()
        
        
        Defined at profiles_module.fpp lines 29-38
        
        
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
        
        
        Defined at profiles_module.fpp lines 29-38
        
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
        
        
        Defined at profiles_module.fpp line 30
        
        """
        return _amr2_pkg.f90wrap_profile_handler__get__profdim(self._handle)
    
    @profdim.setter
    def profdim(self, profdim):
        _amr2_pkg.f90wrap_profile_handler__set__profdim(self._handle, profdim)
    
    @property
    def xvarname(self):
        """
        Element xvarname ftype=character(128) pytype=str
        
        
        Defined at profiles_module.fpp line 31
        
        """
        return _amr2_pkg.f90wrap_profile_handler__get__xvarname(self._handle)
    
    @xvarname.setter
    def xvarname(self, xvarname):
        _amr2_pkg.f90wrap_profile_handler__set__xvarname(self._handle, xvarname)
    
    @property
    def nyvar(self):
        """
        Element nyvar ftype=integer  pytype=int
        
        
        Defined at profiles_module.fpp line 32
        
        """
        return _amr2_pkg.f90wrap_profile_handler__get__nyvar(self._handle)
    
    @nyvar.setter
    def nyvar(self, nyvar):
        _amr2_pkg.f90wrap_profile_handler__set__nyvar(self._handle, nyvar)
    
    @property
    def yvarnames(self):
        """
        Element yvarnames ftype=character(128) pytype=str
        
        
        Defined at profiles_module.fpp line 33
        
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
        
        
        Defined at profiles_module.fpp line 34
        
        """
        return _amr2_pkg.f90wrap_profile_handler__get__nbins(self._handle)
    
    @nbins.setter
    def nbins(self, nbins):
        _amr2_pkg.f90wrap_profile_handler__set__nbins(self._handle, nbins)
    
    @property
    def nwvar(self):
        """
        Element nwvar ftype=integer  pytype=int
        
        
        Defined at profiles_module.fpp line 35
        
        """
        return _amr2_pkg.f90wrap_profile_handler__get__nwvar(self._handle)
    
    @nwvar.setter
    def nwvar(self, nwvar):
        _amr2_pkg.f90wrap_profile_handler__set__nwvar(self._handle, nwvar)
    
    @property
    def wvarnames(self):
        """
        Element wvarnames ftype=character(128) pytype=str
        
        
        Defined at profiles_module.fpp line 36
        
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
        
        
        Defined at profiles_module.fpp line 37
        
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
        
        
        Defined at profiles_module.fpp line 38
        
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
    

@f90wrap.runtime.register_class("amr2_pkg.profile_handler_twod")
class profile_handler_twod(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=profile_handler_twod)
    
    
    Defined at profiles_module.fpp lines 40-50
    
    """
    def __init__(self, handle=None):
        """
        self = Profile_Handler_Twod()
        
        
        Defined at profiles_module.fpp lines 40-50
        
        
        Returns
        -------
        this : Profile_Handler_Twod
        	Object to be constructed
        
        
        Automatically generated constructor for profile_handler_twod
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _amr2_pkg.f90wrap_profile_handler_twod_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Profile_Handler_Twod
        
        
        Defined at profiles_module.fpp lines 40-50
        
        Parameters
        ----------
        this : Profile_Handler_Twod
        	Object to be destructed
        
        
        Automatically generated destructor for profile_handler_twod
        """
        if self._alloc:
            _amr2_pkg.f90wrap_profile_handler_twod_finalise(this=self._handle)
    
    @property
    def profdim(self):
        """
        Element profdim ftype=integer  pytype=int
        
        
        Defined at profiles_module.fpp line 41
        
        """
        return _amr2_pkg.f90wrap_profile_handler_twod__get__profdim(self._handle)
    
    @profdim.setter
    def profdim(self, profdim):
        _amr2_pkg.f90wrap_profile_handler_twod__set__profdim(self._handle, profdim)
    
    @property
    def xvarname(self):
        """
        Element xvarname ftype=character(128) pytype=str
        
        
        Defined at profiles_module.fpp line 42
        
        """
        return _amr2_pkg.f90wrap_profile_handler_twod__get__xvarname(self._handle)
    
    @xvarname.setter
    def xvarname(self, xvarname):
        _amr2_pkg.f90wrap_profile_handler_twod__set__xvarname(self._handle, xvarname)
    
    @property
    def yvarname(self):
        """
        Element yvarname ftype=character(128) pytype=str
        
        
        Defined at profiles_module.fpp line 43
        
        """
        return _amr2_pkg.f90wrap_profile_handler_twod__get__yvarname(self._handle)
    
    @yvarname.setter
    def yvarname(self, yvarname):
        _amr2_pkg.f90wrap_profile_handler_twod__set__yvarname(self._handle, yvarname)
    
    @property
    def nzvar(self):
        """
        Element nzvar ftype=integer  pytype=int
        
        
        Defined at profiles_module.fpp line 44
        
        """
        return _amr2_pkg.f90wrap_profile_handler_twod__get__nzvar(self._handle)
    
    @nzvar.setter
    def nzvar(self, nzvar):
        _amr2_pkg.f90wrap_profile_handler_twod__set__nzvar(self._handle, nzvar)
    
    @property
    def zvarnames(self):
        """
        Element zvarnames ftype=character(128) pytype=str
        
        
        Defined at profiles_module.fpp line 45
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_profile_handler_twod__array__zvarnames(self._handle)
        if array_handle in self._arrays:
            zvarnames = self._arrays[array_handle]
        else:
            zvarnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_profile_handler_twod__array__zvarnames)
            self._arrays[array_handle] = zvarnames
        return zvarnames
    
    @zvarnames.setter
    def zvarnames(self, zvarnames):
        self.zvarnames[...] = zvarnames
    
    @property
    def nbins(self):
        """
        Element nbins ftype=integer pytype=int
        
        
        Defined at profiles_module.fpp line 46
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_profile_handler_twod__array__nbins(self._handle)
        if array_handle in self._arrays:
            nbins = self._arrays[array_handle]
        else:
            nbins = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_profile_handler_twod__array__nbins)
            self._arrays[array_handle] = nbins
        return nbins
    
    @nbins.setter
    def nbins(self, nbins):
        self.nbins[...] = nbins
    
    @property
    def nwvar(self):
        """
        Element nwvar ftype=integer  pytype=int
        
        
        Defined at profiles_module.fpp line 47
        
        """
        return _amr2_pkg.f90wrap_profile_handler_twod__get__nwvar(self._handle)
    
    @nwvar.setter
    def nwvar(self, nwvar):
        _amr2_pkg.f90wrap_profile_handler_twod__set__nwvar(self._handle, nwvar)
    
    @property
    def wvarnames(self):
        """
        Element wvarnames ftype=character(128) pytype=str
        
        
        Defined at profiles_module.fpp line 48
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_profile_handler_twod__array__wvarnames(self._handle)
        if array_handle in self._arrays:
            wvarnames = self._arrays[array_handle]
        else:
            wvarnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_profile_handler_twod__array__wvarnames)
            self._arrays[array_handle] = wvarnames
        return wvarnames
    
    @wvarnames.setter
    def wvarnames(self, wvarnames):
        self.wvarnames[...] = wvarnames
    
    @property
    def xdata(self):
        """
        Element xdata ftype=real(dbl) pytype=float
        
        
        Defined at profiles_module.fpp line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_profile_handler_twod__array__xdata(self._handle)
        if array_handle in self._arrays:
            xdata = self._arrays[array_handle]
        else:
            xdata = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_profile_handler_twod__array__xdata)
            self._arrays[array_handle] = xdata
        return xdata
    
    @xdata.setter
    def xdata(self, xdata):
        self.xdata[...] = xdata
    
    @property
    def ydata(self):
        """
        Element ydata ftype=real(dbl) pytype=float
        
        
        Defined at profiles_module.fpp line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_profile_handler_twod__array__ydata(self._handle)
        if array_handle in self._arrays:
            ydata = self._arrays[array_handle]
        else:
            ydata = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_profile_handler_twod__array__ydata)
            self._arrays[array_handle] = ydata
        return ydata
    
    @ydata.setter
    def ydata(self, ydata):
        self.ydata[...] = ydata
    
    @property
    def zdata(self):
        """
        Element zdata ftype=real(dbl) pytype=float
        
        
        Defined at profiles_module.fpp line 50
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_profile_handler_twod__array__zdata(self._handle)
        if array_handle in self._arrays:
            zdata = self._arrays[array_handle]
        else:
            zdata = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_profile_handler_twod__array__zdata)
            self._arrays[array_handle] = zdata
        return zdata
    
    @zdata.setter
    def zdata(self, zdata):
        self.zdata[...] = zdata
    
    def __str__(self):
        ret = ['<profile_handler_twod>{\n']
        ret.append('    profdim : ')
        ret.append(repr(self.profdim))
        ret.append(',\n    xvarname : ')
        ret.append(repr(self.xvarname))
        ret.append(',\n    yvarname : ')
        ret.append(repr(self.yvarname))
        ret.append(',\n    nzvar : ')
        ret.append(repr(self.nzvar))
        ret.append(',\n    zvarnames : ')
        ret.append(repr(self.zvarnames))
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
        ret.append(',\n    zdata : ')
        ret.append(repr(self.zdata))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def allocate_profile_handler(self):
    """
    allocate_profile_handler(self)
    
    
    Defined at profiles_module.fpp lines 53-59
    
    Parameters
    ----------
    prof : Profile_Handler
    
    """
    _amr2_pkg.f90wrap_allocate_profile_handler(prof=self._handle)

def allocate_profile_handler_twod(self):
    """
    allocate_profile_handler_twod(self)
    
    
    Defined at profiles_module.fpp lines 61-68
    
    Parameters
    ----------
    prof : Profile_Handler_Twod
    
    """
    _amr2_pkg.f90wrap_allocate_profile_handler_twod(prof=self._handle)

def makebins(self, sim, varname, nbins, bins, logscale):
    """
    makebins(self, sim, varname, nbins, bins, logscale)
    
    
    Defined at profiles_module.fpp lines 70-109
    
    Parameters
    ----------
    reg : Region
    sim : Sim_Info
    varname : str
    nbins : int
    bins : float array
    logscale : bool
    
    """
    _amr2_pkg.f90wrap_makebins(reg=self._handle, sim=sim._handle, varname=varname, \
        nbins=nbins, bins=bins, logscale=logscale)

def findbinpos(self, varids, distance, pos, cellvars, cellsize, prof, ibin):
    """
    findbinpos(self, varids, distance, pos, cellvars, cellsize, prof, ibin)
    
    
    Defined at profiles_module.fpp lines 111-133
    
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

def findbinpos_twod(self, varids, distance, pos, cellvars, cellsize, prof, \
    logscale, ibinx, ibiny):
    """
    findbinpos_twod(self, varids, distance, pos, cellvars, cellsize, prof, logscale, \
        ibinx, ibiny)
    
    
    Defined at profiles_module.fpp lines 135-174
    
    Parameters
    ----------
    reg : Region
    varids : Hydroid
    distance : float
    pos : float array
    cellvars : float array
    cellsize : float
    prof : Profile_Handler_Twod
    logscale : bool
    ibinx : int
    ibiny : int
    
    """
    _amr2_pkg.f90wrap_findbinpos_twod(reg=self._handle, varids=varids._handle, \
        distance=distance, pos=pos, cellvars=cellvars, cellsize=cellsize, \
        prof=prof._handle, logscale=logscale, ibinx=ibinx, ibiny=ibiny)

def bindata(self, varids, pos, cellvars, cellsize, prof, ibin):
    """
    bindata(self, varids, pos, cellvars, cellsize, prof, ibin)
    
    
    Defined at profiles_module.fpp lines 176-223
    
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

def bindata_twod(self, varids, pos, cellvars, cellsize, prof, ibinx, ibiny):
    """
    bindata_twod(self, varids, pos, cellvars, cellsize, prof, ibinx, ibiny)
    
    
    Defined at profiles_module.fpp lines 225-258
    
    Parameters
    ----------
    reg : Region
    varids : Hydroid
    pos : float array
    cellvars : float array
    cellsize : float
    prof : Profile_Handler_Twod
    ibinx : int
    ibiny : int
    
    """
    _amr2_pkg.f90wrap_bindata_twod(reg=self._handle, varids=varids._handle, pos=pos, \
        cellvars=cellvars, cellsize=cellsize, prof=prof._handle, ibinx=ibinx, \
        ibiny=ibiny)

def renormalise_bins(self):
    """
    renormalise_bins(self)
    
    
    Defined at profiles_module.fpp lines 260-285
    
    Parameters
    ----------
    prof_data : Profile_Handler
    
    """
    _amr2_pkg.f90wrap_renormalise_bins(prof_data=self._handle)

def renormalise_bins_twod(self):
    """
    renormalise_bins_twod(self)
    
    
    Defined at profiles_module.fpp lines 287-311
    
    Parameters
    ----------
    prof_data : Profile_Handler_Twod
    
    """
    _amr2_pkg.f90wrap_renormalise_bins_twod(prof_data=self._handle)

def get_cells_onedprofile(repository, amr, reg, filt, varids, prof_data):
    """
    get_cells_onedprofile(repository, amr, reg, filt, varids, prof_data)
    
    
    Defined at profiles_module.fpp lines 313-515
    
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

def onedprofile(repository, reg, filt, prof_data, lmax, logscale):
    """
    onedprofile(repository, reg, filt, prof_data, lmax, logscale)
    
    
    Defined at profiles_module.fpp lines 517-539
    
    Parameters
    ----------
    repository : str
    reg : Region
    filt : Filter
    prof_data : Profile_Handler
    lmax : int
    logscale : bool
    
    """
    _amr2_pkg.f90wrap_onedprofile(repository=repository, reg=reg._handle, \
        filt=filt._handle, prof_data=prof_data._handle, lmax=lmax, \
        logscale=logscale)

def twodprofile(repository, reg, filt, prof_data, lmax, logscale):
    """
    twodprofile(repository, reg, filt, prof_data, lmax, logscale)
    
    
    Defined at profiles_module.fpp lines 541-566
    
    Parameters
    ----------
    repository : str
    reg : Region
    filt : Filter
    prof_data : Profile_Handler_Twod
    lmax : int
    logscale : bool
    
    """
    _amr2_pkg.f90wrap_twodprofile(repository=repository, reg=reg._handle, \
        filt=filt._handle, prof_data=prof_data._handle, lmax=lmax, \
        logscale=logscale)

def get_cells_twodprofile(repository, amr, reg, filt, varids, prof_data, \
    logscale):
    """
    get_cells_twodprofile(repository, amr, reg, filt, varids, prof_data, logscale)
    
    
    Defined at profiles_module.fpp lines 568-770
    
    Parameters
    ----------
    repository : str
    amr : Amr_Info
    reg : Region
    filt : Filter
    varids : Hydroid
    prof_data : Profile_Handler_Twod
    logscale : bool
    
    """
    _amr2_pkg.f90wrap_get_cells_twodprofile(repository=repository, amr=amr._handle, \
        reg=reg._handle, filt=filt._handle, varids=varids._handle, \
        prof_data=prof_data._handle, logscale=logscale)


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
