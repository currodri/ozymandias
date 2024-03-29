"""
Module io_ramses


Defined at read_amr_module.fpp lines 24-2514

"""
from __future__ import print_function, absolute_import, division
import _visualisation_pkg
import f90wrap.runtime
import logging
from visualisation_pkg.vectors import vector

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("visualisation_pkg.hydroID")
class hydroID(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=hydroid)
    
    
    Defined at read_amr_module.fpp lines 29-36
    
    """
    def __init__(self, handle=None):
        """
        self = Hydroid()
        
        
        Defined at read_amr_module.fpp lines 29-36
        
        
        Returns
        -------
        this : Hydroid
        	Object to be constructed
        
        
        Automatically generated constructor for hydroid
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _visualisation_pkg.f90wrap_hydroid_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Hydroid
        
        
        Defined at read_amr_module.fpp lines 29-36
        
        Parameters
        ----------
        this : Hydroid
        	Object to be destructed
        
        
        Automatically generated destructor for hydroid
        """
        if self._alloc:
            _visualisation_pkg.f90wrap_hydroid_finalise(this=self._handle)
    
    @property
    def nvar(self):
        """
        Element nvar ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 30
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__nvar(self._handle)
    
    @nvar.setter
    def nvar(self, nvar):
        _visualisation_pkg.f90wrap_hydroid__set__nvar(self._handle, nvar)
    
    @property
    def density(self):
        """
        Element density ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 31
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__density(self._handle)
    
    @density.setter
    def density(self, density):
        _visualisation_pkg.f90wrap_hydroid__set__density(self._handle, density)
    
    @property
    def vx(self):
        """
        Element vx ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 31
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__vx(self._handle)
    
    @vx.setter
    def vx(self, vx):
        _visualisation_pkg.f90wrap_hydroid__set__vx(self._handle, vx)
    
    @property
    def vy(self):
        """
        Element vy ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 31
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__vy(self._handle)
    
    @vy.setter
    def vy(self, vy):
        _visualisation_pkg.f90wrap_hydroid__set__vy(self._handle, vy)
    
    @property
    def vz(self):
        """
        Element vz ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 31
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__vz(self._handle)
    
    @vz.setter
    def vz(self, vz):
        _visualisation_pkg.f90wrap_hydroid__set__vz(self._handle, vz)
    
    @property
    def thermal_pressure(self):
        """
        Element thermal_pressure ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 31
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__thermal_pressure(self._handle)
    
    @thermal_pressure.setter
    def thermal_pressure(self, thermal_pressure):
        _visualisation_pkg.f90wrap_hydroid__set__thermal_pressure(self._handle, \
            thermal_pressure)
    
    @property
    def metallicity(self):
        """
        Element metallicity ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 31
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__metallicity(self._handle)
    
    @metallicity.setter
    def metallicity(self, metallicity):
        _visualisation_pkg.f90wrap_hydroid__set__metallicity(self._handle, metallicity)
    
    @property
    def blx(self):
        """
        Element blx ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 32
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__blx(self._handle)
    
    @blx.setter
    def blx(self, blx):
        _visualisation_pkg.f90wrap_hydroid__set__blx(self._handle, blx)
    
    @property
    def bly(self):
        """
        Element bly ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 32
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__bly(self._handle)
    
    @bly.setter
    def bly(self, bly):
        _visualisation_pkg.f90wrap_hydroid__set__bly(self._handle, bly)
    
    @property
    def blz(self):
        """
        Element blz ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 32
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__blz(self._handle)
    
    @blz.setter
    def blz(self, blz):
        _visualisation_pkg.f90wrap_hydroid__set__blz(self._handle, blz)
    
    @property
    def brx(self):
        """
        Element brx ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 32
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__brx(self._handle)
    
    @brx.setter
    def brx(self, brx):
        _visualisation_pkg.f90wrap_hydroid__set__brx(self._handle, brx)
    
    @property
    def bry(self):
        """
        Element bry ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 32
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__bry(self._handle)
    
    @bry.setter
    def bry(self, bry):
        _visualisation_pkg.f90wrap_hydroid__set__bry(self._handle, bry)
    
    @property
    def brz(self):
        """
        Element brz ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 32
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__brz(self._handle)
    
    @brz.setter
    def brz(self, brz):
        _visualisation_pkg.f90wrap_hydroid__set__brz(self._handle, brz)
    
    @property
    def cr_pressure(self):
        """
        Element cr_pressure ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 33
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__cr_pressure(self._handle)
    
    @cr_pressure.setter
    def cr_pressure(self, cr_pressure):
        _visualisation_pkg.f90wrap_hydroid__set__cr_pressure(self._handle, cr_pressure)
    
    @property
    def xhii(self):
        """
        Element xhii ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 34
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__xhii(self._handle)
    
    @xhii.setter
    def xhii(self, xhii):
        _visualisation_pkg.f90wrap_hydroid__set__xhii(self._handle, xhii)
    
    @property
    def xheii(self):
        """
        Element xheii ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 34
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__xheii(self._handle)
    
    @xheii.setter
    def xheii(self, xheii):
        _visualisation_pkg.f90wrap_hydroid__set__xheii(self._handle, xheii)
    
    @property
    def xheiii(self):
        """
        Element xheiii ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 34
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__xheiii(self._handle)
    
    @xheiii.setter
    def xheiii(self, xheiii):
        _visualisation_pkg.f90wrap_hydroid__set__xheiii(self._handle, xheiii)
    
    @property
    def dust_density(self):
        """
        Element dust_density ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 35
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__dust_density(self._handle)
    
    @dust_density.setter
    def dust_density(self, dust_density):
        _visualisation_pkg.f90wrap_hydroid__set__dust_density(self._handle, \
            dust_density)
    
    @property
    def sigma2(self):
        """
        Element sigma2 ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 36
        
        """
        return _visualisation_pkg.f90wrap_hydroid__get__sigma2(self._handle)
    
    @sigma2.setter
    def sigma2(self, sigma2):
        _visualisation_pkg.f90wrap_hydroid__set__sigma2(self._handle, sigma2)
    
    def __str__(self):
        ret = ['<hydroid>{\n']
        ret.append('    nvar : ')
        ret.append(repr(self.nvar))
        ret.append(',\n    density : ')
        ret.append(repr(self.density))
        ret.append(',\n    vx : ')
        ret.append(repr(self.vx))
        ret.append(',\n    vy : ')
        ret.append(repr(self.vy))
        ret.append(',\n    vz : ')
        ret.append(repr(self.vz))
        ret.append(',\n    thermal_pressure : ')
        ret.append(repr(self.thermal_pressure))
        ret.append(',\n    metallicity : ')
        ret.append(repr(self.metallicity))
        ret.append(',\n    blx : ')
        ret.append(repr(self.blx))
        ret.append(',\n    bly : ')
        ret.append(repr(self.bly))
        ret.append(',\n    blz : ')
        ret.append(repr(self.blz))
        ret.append(',\n    brx : ')
        ret.append(repr(self.brx))
        ret.append(',\n    bry : ')
        ret.append(repr(self.bry))
        ret.append(',\n    brz : ')
        ret.append(repr(self.brz))
        ret.append(',\n    cr_pressure : ')
        ret.append(repr(self.cr_pressure))
        ret.append(',\n    xhii : ')
        ret.append(repr(self.xhii))
        ret.append(',\n    xheii : ')
        ret.append(repr(self.xheii))
        ret.append(',\n    xheiii : ')
        ret.append(repr(self.xheiii))
        ret.append(',\n    dust_density : ')
        ret.append(repr(self.dust_density))
        ret.append(',\n    sigma2 : ')
        ret.append(repr(self.sigma2))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("visualisation_pkg.amr_info")
class amr_info(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=amr_info)
    
    
    Defined at read_amr_module.fpp lines 38-47
    
    """
    def __init__(self, handle=None):
        """
        self = Amr_Info()
        
        
        Defined at read_amr_module.fpp lines 38-47
        
        
        Returns
        -------
        this : Amr_Info
        	Object to be constructed
        
        
        Automatically generated constructor for amr_info
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _visualisation_pkg.f90wrap_amr_info_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Amr_Info
        
        
        Defined at read_amr_module.fpp lines 38-47
        
        Parameters
        ----------
        this : Amr_Info
        	Object to be destructed
        
        
        Automatically generated destructor for amr_info
        """
        if self._alloc:
            _visualisation_pkg.f90wrap_amr_info_finalise(this=self._handle)
    
    @property
    def ncpu(self):
        """
        Element ncpu ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 39
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__ncpu(self._handle)
    
    @ncpu.setter
    def ncpu(self, ncpu):
        _visualisation_pkg.f90wrap_amr_info__set__ncpu(self._handle, ncpu)
    
    @property
    def ndim(self):
        """
        Element ndim ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 39
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__ndim(self._handle)
    
    @ndim.setter
    def ndim(self, ndim):
        _visualisation_pkg.f90wrap_amr_info__set__ndim(self._handle, ndim)
    
    @property
    def nlevelmax(self):
        """
        Element nlevelmax ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 39
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__nlevelmax(self._handle)
    
    @nlevelmax.setter
    def nlevelmax(self, nlevelmax):
        _visualisation_pkg.f90wrap_amr_info__set__nlevelmax(self._handle, nlevelmax)
    
    @property
    def nboundary(self):
        """
        Element nboundary ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 39
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__nboundary(self._handle)
    
    @nboundary.setter
    def nboundary(self, nboundary):
        _visualisation_pkg.f90wrap_amr_info__set__nboundary(self._handle, nboundary)
    
    @property
    def twotondim(self):
        """
        Element twotondim ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 39
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__twotondim(self._handle)
    
    @twotondim.setter
    def twotondim(self, twotondim):
        _visualisation_pkg.f90wrap_amr_info__set__twotondim(self._handle, twotondim)
    
    @property
    def ndom(self):
        """
        Element ndom ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 39
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__ndom(self._handle)
    
    @ndom.setter
    def ndom(self, ndom):
        _visualisation_pkg.f90wrap_amr_info__set__ndom(self._handle, ndom)
    
    @property
    def twondim(self):
        """
        Element twondim ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 40
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__twondim(self._handle)
    
    @twondim.setter
    def twondim(self, twondim):
        _visualisation_pkg.f90wrap_amr_info__set__twondim(self._handle, twondim)
    
    @property
    def ngridmax(self):
        """
        Element ngridmax ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 40
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__ngridmax(self._handle)
    
    @ngridmax.setter
    def ngridmax(self, ngridmax):
        _visualisation_pkg.f90wrap_amr_info__set__ngridmax(self._handle, ngridmax)
    
    @property
    def ncoarse(self):
        """
        Element ncoarse ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 40
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__ncoarse(self._handle)
    
    @ncoarse.setter
    def ncoarse(self, ncoarse):
        _visualisation_pkg.f90wrap_amr_info__set__ncoarse(self._handle, ncoarse)
    
    @property
    def levelmin(self):
        """
        Element levelmin ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 41
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__levelmin(self._handle)
    
    @levelmin.setter
    def levelmin(self, levelmin):
        _visualisation_pkg.f90wrap_amr_info__set__levelmin(self._handle, levelmin)
    
    @property
    def levelmax(self):
        """
        Element levelmax ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 41
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__levelmax(self._handle)
    
    @levelmax.setter
    def levelmax(self, levelmax):
        _visualisation_pkg.f90wrap_amr_info__set__levelmax(self._handle, levelmax)
    
    @property
    def lmax(self):
        """
        Element lmax ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 41
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__lmax(self._handle)
    
    @lmax.setter
    def lmax(self, lmax):
        _visualisation_pkg.f90wrap_amr_info__set__lmax(self._handle, lmax)
    
    @property
    def lmin(self):
        """
        Element lmin ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 41
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__lmin(self._handle)
    
    @lmin.setter
    def lmin(self, lmin):
        _visualisation_pkg.f90wrap_amr_info__set__lmin(self._handle, lmin)
    
    @property
    def active_lmax(self):
        """
        Element active_lmax ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 41
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__active_lmax(self._handle)
    
    @active_lmax.setter
    def active_lmax(self, active_lmax):
        _visualisation_pkg.f90wrap_amr_info__set__active_lmax(self._handle, active_lmax)
    
    @property
    def ncpu_read(self):
        """
        Element ncpu_read ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 42
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__ncpu_read(self._handle)
    
    @ncpu_read.setter
    def ncpu_read(self, ncpu_read):
        _visualisation_pkg.f90wrap_amr_info__set__ncpu_read(self._handle, ncpu_read)
    
    @property
    def ordering(self):
        """
        Element ordering ftype=character(80) pytype=str
        
        
        Defined at read_amr_module.fpp line 43
        
        """
        return _visualisation_pkg.f90wrap_amr_info__get__ordering(self._handle)
    
    @ordering.setter
    def ordering(self, ordering):
        _visualisation_pkg.f90wrap_amr_info__set__ordering(self._handle, ordering)
    
    @property
    def cpu_list(self):
        """
        Element cpu_list ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 44
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_amr_info__array__cpu_list(self._handle)
        if array_handle in self._arrays:
            cpu_list = self._arrays[array_handle]
        else:
            cpu_list = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_amr_info__array__cpu_list)
            self._arrays[array_handle] = cpu_list
        return cpu_list
    
    @cpu_list.setter
    def cpu_list(self, cpu_list):
        self.cpu_list[...] = cpu_list
    
    @property
    def bound_key(self):
        """
        Element bound_key ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 45
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_amr_info__array__bound_key(self._handle)
        if array_handle in self._arrays:
            bound_key = self._arrays[array_handle]
        else:
            bound_key = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_amr_info__array__bound_key)
            self._arrays[array_handle] = bound_key
        return bound_key
    
    @bound_key.setter
    def bound_key(self, bound_key):
        self.bound_key[...] = bound_key
    
    @property
    def cpu_read(self):
        """
        Element cpu_read ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 46
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_amr_info__array__cpu_read(self._handle)
        if array_handle in self._arrays:
            cpu_read = self._arrays[array_handle]
        else:
            cpu_read = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_amr_info__array__cpu_read)
            self._arrays[array_handle] = cpu_read
        return cpu_read
    
    @cpu_read.setter
    def cpu_read(self, cpu_read):
        self.cpu_read[...] = cpu_read
    
    @property
    def xbound(self):
        """
        Element xbound ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 47
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_amr_info__array__xbound(self._handle)
        if array_handle in self._arrays:
            xbound = self._arrays[array_handle]
        else:
            xbound = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_amr_info__array__xbound)
            self._arrays[array_handle] = xbound
        return xbound
    
    @xbound.setter
    def xbound(self, xbound):
        self.xbound[...] = xbound
    
    def __str__(self):
        ret = ['<amr_info>{\n']
        ret.append('    ncpu : ')
        ret.append(repr(self.ncpu))
        ret.append(',\n    ndim : ')
        ret.append(repr(self.ndim))
        ret.append(',\n    nlevelmax : ')
        ret.append(repr(self.nlevelmax))
        ret.append(',\n    nboundary : ')
        ret.append(repr(self.nboundary))
        ret.append(',\n    twotondim : ')
        ret.append(repr(self.twotondim))
        ret.append(',\n    ndom : ')
        ret.append(repr(self.ndom))
        ret.append(',\n    twondim : ')
        ret.append(repr(self.twondim))
        ret.append(',\n    ngridmax : ')
        ret.append(repr(self.ngridmax))
        ret.append(',\n    ncoarse : ')
        ret.append(repr(self.ncoarse))
        ret.append(',\n    levelmin : ')
        ret.append(repr(self.levelmin))
        ret.append(',\n    levelmax : ')
        ret.append(repr(self.levelmax))
        ret.append(',\n    lmax : ')
        ret.append(repr(self.lmax))
        ret.append(',\n    lmin : ')
        ret.append(repr(self.lmin))
        ret.append(',\n    active_lmax : ')
        ret.append(repr(self.active_lmax))
        ret.append(',\n    ncpu_read : ')
        ret.append(repr(self.ncpu_read))
        ret.append(',\n    ordering : ')
        ret.append(repr(self.ordering))
        ret.append(',\n    cpu_list : ')
        ret.append(repr(self.cpu_list))
        ret.append(',\n    bound_key : ')
        ret.append(repr(self.bound_key))
        ret.append(',\n    cpu_read : ')
        ret.append(repr(self.cpu_read))
        ret.append(',\n    xbound : ')
        ret.append(repr(self.xbound))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("visualisation_pkg.sim_info")
class sim_info(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=sim_info)
    
    
    Defined at read_amr_module.fpp lines 49-57
    
    """
    def __init__(self, handle=None):
        """
        self = Sim_Info()
        
        
        Defined at read_amr_module.fpp lines 49-57
        
        
        Returns
        -------
        this : Sim_Info
        	Object to be constructed
        
        
        Automatically generated constructor for sim_info
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _visualisation_pkg.f90wrap_sim_info_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Sim_Info
        
        
        Defined at read_amr_module.fpp lines 49-57
        
        Parameters
        ----------
        this : Sim_Info
        	Object to be destructed
        
        
        Automatically generated destructor for sim_info
        """
        if self._alloc:
            _visualisation_pkg.f90wrap_sim_info_finalise(this=self._handle)
    
    @property
    def cosmo(self):
        """
        Element cosmo ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 50
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__cosmo(self._handle)
    
    @cosmo.setter
    def cosmo(self, cosmo):
        _visualisation_pkg.f90wrap_sim_info__set__cosmo(self._handle, cosmo)
    
    @property
    def family(self):
        """
        Element family ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 50
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__family(self._handle)
    
    @family.setter
    def family(self, family):
        _visualisation_pkg.f90wrap_sim_info__set__family(self._handle, family)
    
    @property
    def dm(self):
        """
        Element dm ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 51
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__dm(self._handle)
    
    @dm.setter
    def dm(self, dm):
        _visualisation_pkg.f90wrap_sim_info__set__dm(self._handle, dm)
    
    @property
    def hydro(self):
        """
        Element hydro ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 51
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__hydro(self._handle)
    
    @hydro.setter
    def hydro(self, hydro):
        _visualisation_pkg.f90wrap_sim_info__set__hydro(self._handle, hydro)
    
    @property
    def mhd(self):
        """
        Element mhd ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 51
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__mhd(self._handle)
    
    @mhd.setter
    def mhd(self, mhd):
        _visualisation_pkg.f90wrap_sim_info__set__mhd(self._handle, mhd)
    
    @property
    def cr(self):
        """
        Element cr ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 51
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__cr(self._handle)
    
    @cr.setter
    def cr(self, cr):
        _visualisation_pkg.f90wrap_sim_info__set__cr(self._handle, cr)
    
    @property
    def rt(self):
        """
        Element rt ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 51
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__rt(self._handle)
    
    @rt.setter
    def rt(self, rt):
        _visualisation_pkg.f90wrap_sim_info__set__rt(self._handle, rt)
    
    @property
    def bh(self):
        """
        Element bh ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 51
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__bh(self._handle)
    
    @bh.setter
    def bh(self, bh):
        _visualisation_pkg.f90wrap_sim_info__set__bh(self._handle, bh)
    
    @property
    def cr_st(self):
        """
        Element cr_st ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 52
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__cr_st(self._handle)
    
    @cr_st.setter
    def cr_st(self, cr_st):
        _visualisation_pkg.f90wrap_sim_info__set__cr_st(self._handle, cr_st)
    
    @property
    def cr_heat(self):
        """
        Element cr_heat ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 52
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__cr_heat(self._handle)
    
    @cr_heat.setter
    def cr_heat(self, cr_heat):
        _visualisation_pkg.f90wrap_sim_info__set__cr_heat(self._handle, cr_heat)
    
    @property
    def dust(self):
        """
        Element dust ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 52
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__dust(self._handle)
    
    @dust.setter
    def dust(self, dust):
        _visualisation_pkg.f90wrap_sim_info__set__dust(self._handle, dust)
    
    @property
    def h0(self):
        """
        Element h0 ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 53
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__h0(self._handle)
    
    @h0.setter
    def h0(self, h0):
        _visualisation_pkg.f90wrap_sim_info__set__h0(self._handle, h0)
    
    @property
    def t(self):
        """
        Element t ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 53
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__t(self._handle)
    
    @t.setter
    def t(self, t):
        _visualisation_pkg.f90wrap_sim_info__set__t(self._handle, t)
    
    @property
    def aexp(self):
        """
        Element aexp ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 53
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__aexp(self._handle)
    
    @aexp.setter
    def aexp(self, aexp):
        _visualisation_pkg.f90wrap_sim_info__set__aexp(self._handle, aexp)
    
    @property
    def unit_l(self):
        """
        Element unit_l ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 53
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__unit_l(self._handle)
    
    @unit_l.setter
    def unit_l(self, unit_l):
        _visualisation_pkg.f90wrap_sim_info__set__unit_l(self._handle, unit_l)
    
    @property
    def unit_d(self):
        """
        Element unit_d ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 53
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__unit_d(self._handle)
    
    @unit_d.setter
    def unit_d(self, unit_d):
        _visualisation_pkg.f90wrap_sim_info__set__unit_d(self._handle, unit_d)
    
    @property
    def unit_t(self):
        """
        Element unit_t ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 53
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__unit_t(self._handle)
    
    @unit_t.setter
    def unit_t(self, unit_t):
        _visualisation_pkg.f90wrap_sim_info__set__unit_t(self._handle, unit_t)
    
    @property
    def unit_m(self):
        """
        Element unit_m ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 53
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__unit_m(self._handle)
    
    @unit_m.setter
    def unit_m(self, unit_m):
        _visualisation_pkg.f90wrap_sim_info__set__unit_m(self._handle, unit_m)
    
    @property
    def unit_v(self):
        """
        Element unit_v ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 53
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__unit_v(self._handle)
    
    @unit_v.setter
    def unit_v(self, unit_v):
        _visualisation_pkg.f90wrap_sim_info__set__unit_v(self._handle, unit_v)
    
    @property
    def boxlen(self):
        """
        Element boxlen ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 53
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__boxlen(self._handle)
    
    @boxlen.setter
    def boxlen(self, boxlen):
        _visualisation_pkg.f90wrap_sim_info__set__boxlen(self._handle, boxlen)
    
    @property
    def omega_m(self):
        """
        Element omega_m ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 53
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__omega_m(self._handle)
    
    @omega_m.setter
    def omega_m(self, omega_m):
        _visualisation_pkg.f90wrap_sim_info__set__omega_m(self._handle, omega_m)
    
    @property
    def omega_l(self):
        """
        Element omega_l ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 53
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__omega_l(self._handle)
    
    @omega_l.setter
    def omega_l(self, omega_l):
        _visualisation_pkg.f90wrap_sim_info__set__omega_l(self._handle, omega_l)
    
    @property
    def omega_k(self):
        """
        Element omega_k ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 53
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__omega_k(self._handle)
    
    @omega_k.setter
    def omega_k(self, omega_k):
        _visualisation_pkg.f90wrap_sim_info__set__omega_k(self._handle, omega_k)
    
    @property
    def omega_b(self):
        """
        Element omega_b ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 53
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__omega_b(self._handle)
    
    @omega_b.setter
    def omega_b(self, omega_b):
        _visualisation_pkg.f90wrap_sim_info__set__omega_b(self._handle, omega_b)
    
    @property
    def time_tot(self):
        """
        Element time_tot ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 54
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__time_tot(self._handle)
    
    @time_tot.setter
    def time_tot(self, time_tot):
        _visualisation_pkg.f90wrap_sim_info__set__time_tot(self._handle, time_tot)
    
    @property
    def time_simu(self):
        """
        Element time_simu ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 54
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__time_simu(self._handle)
    
    @time_simu.setter
    def time_simu(self, time_simu):
        _visualisation_pkg.f90wrap_sim_info__set__time_simu(self._handle, time_simu)
    
    @property
    def redshift(self):
        """
        Element redshift ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 54
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__redshift(self._handle)
    
    @redshift.setter
    def redshift(self, redshift):
        _visualisation_pkg.f90wrap_sim_info__set__redshift(self._handle, redshift)
    
    @property
    def t2(self):
        """
        Element t2 ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 54
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__t2(self._handle)
    
    @t2.setter
    def t2(self, t2):
        _visualisation_pkg.f90wrap_sim_info__set__t2(self._handle, t2)
    
    @property
    def nh(self):
        """
        Element nh ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 54
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__nh(self._handle)
    
    @nh.setter
    def nh(self, nh):
        _visualisation_pkg.f90wrap_sim_info__set__nh(self._handle, nh)
    
    @property
    def n_frw(self):
        """
        Element n_frw ftype=integer  pytype=int
        
        
        Defined at read_amr_module.fpp line 55
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__n_frw(self._handle)
    
    @n_frw.setter
    def n_frw(self, n_frw):
        _visualisation_pkg.f90wrap_sim_info__set__n_frw(self._handle, n_frw)
    
    @property
    def aexp_frw(self):
        """
        Element aexp_frw ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 56
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_sim_info__array__aexp_frw(self._handle)
        if array_handle in self._arrays:
            aexp_frw = self._arrays[array_handle]
        else:
            aexp_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_sim_info__array__aexp_frw)
            self._arrays[array_handle] = aexp_frw
        return aexp_frw
    
    @aexp_frw.setter
    def aexp_frw(self, aexp_frw):
        self.aexp_frw[...] = aexp_frw
    
    @property
    def hexp_frw(self):
        """
        Element hexp_frw ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 56
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_sim_info__array__hexp_frw(self._handle)
        if array_handle in self._arrays:
            hexp_frw = self._arrays[array_handle]
        else:
            hexp_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_sim_info__array__hexp_frw)
            self._arrays[array_handle] = hexp_frw
        return hexp_frw
    
    @hexp_frw.setter
    def hexp_frw(self, hexp_frw):
        self.hexp_frw[...] = hexp_frw
    
    @property
    def tau_frw(self):
        """
        Element tau_frw ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 56
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_sim_info__array__tau_frw(self._handle)
        if array_handle in self._arrays:
            tau_frw = self._arrays[array_handle]
        else:
            tau_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_sim_info__array__tau_frw)
            self._arrays[array_handle] = tau_frw
        return tau_frw
    
    @tau_frw.setter
    def tau_frw(self, tau_frw):
        self.tau_frw[...] = tau_frw
    
    @property
    def t_frw(self):
        """
        Element t_frw ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 56
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_sim_info__array__t_frw(self._handle)
        if array_handle in self._arrays:
            t_frw = self._arrays[array_handle]
        else:
            t_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_sim_info__array__t_frw)
            self._arrays[array_handle] = t_frw
        return t_frw
    
    @t_frw.setter
    def t_frw(self, t_frw):
        self.t_frw[...] = t_frw
    
    @property
    def eta_sn(self):
        """
        Element eta_sn ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 57
        
        """
        return _visualisation_pkg.f90wrap_sim_info__get__eta_sn(self._handle)
    
    @eta_sn.setter
    def eta_sn(self, eta_sn):
        _visualisation_pkg.f90wrap_sim_info__set__eta_sn(self._handle, eta_sn)
    
    def __str__(self):
        ret = ['<sim_info>{\n']
        ret.append('    cosmo : ')
        ret.append(repr(self.cosmo))
        ret.append(',\n    family : ')
        ret.append(repr(self.family))
        ret.append(',\n    dm : ')
        ret.append(repr(self.dm))
        ret.append(',\n    hydro : ')
        ret.append(repr(self.hydro))
        ret.append(',\n    mhd : ')
        ret.append(repr(self.mhd))
        ret.append(',\n    cr : ')
        ret.append(repr(self.cr))
        ret.append(',\n    rt : ')
        ret.append(repr(self.rt))
        ret.append(',\n    bh : ')
        ret.append(repr(self.bh))
        ret.append(',\n    cr_st : ')
        ret.append(repr(self.cr_st))
        ret.append(',\n    cr_heat : ')
        ret.append(repr(self.cr_heat))
        ret.append(',\n    dust : ')
        ret.append(repr(self.dust))
        ret.append(',\n    h0 : ')
        ret.append(repr(self.h0))
        ret.append(',\n    t : ')
        ret.append(repr(self.t))
        ret.append(',\n    aexp : ')
        ret.append(repr(self.aexp))
        ret.append(',\n    unit_l : ')
        ret.append(repr(self.unit_l))
        ret.append(',\n    unit_d : ')
        ret.append(repr(self.unit_d))
        ret.append(',\n    unit_t : ')
        ret.append(repr(self.unit_t))
        ret.append(',\n    unit_m : ')
        ret.append(repr(self.unit_m))
        ret.append(',\n    unit_v : ')
        ret.append(repr(self.unit_v))
        ret.append(',\n    boxlen : ')
        ret.append(repr(self.boxlen))
        ret.append(',\n    omega_m : ')
        ret.append(repr(self.omega_m))
        ret.append(',\n    omega_l : ')
        ret.append(repr(self.omega_l))
        ret.append(',\n    omega_k : ')
        ret.append(repr(self.omega_k))
        ret.append(',\n    omega_b : ')
        ret.append(repr(self.omega_b))
        ret.append(',\n    time_tot : ')
        ret.append(repr(self.time_tot))
        ret.append(',\n    time_simu : ')
        ret.append(repr(self.time_simu))
        ret.append(',\n    redshift : ')
        ret.append(repr(self.redshift))
        ret.append(',\n    t2 : ')
        ret.append(repr(self.t2))
        ret.append(',\n    nh : ')
        ret.append(repr(self.nh))
        ret.append(',\n    n_frw : ')
        ret.append(repr(self.n_frw))
        ret.append(',\n    aexp_frw : ')
        ret.append(repr(self.aexp_frw))
        ret.append(',\n    hexp_frw : ')
        ret.append(repr(self.hexp_frw))
        ret.append(',\n    tau_frw : ')
        ret.append(repr(self.tau_frw))
        ret.append(',\n    t_frw : ')
        ret.append(repr(self.t_frw))
        ret.append(',\n    eta_sn : ')
        ret.append(repr(self.eta_sn))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("visualisation_pkg.level")
class level(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=level)
    
    
    Defined at read_amr_module.fpp lines 59-73
    
    """
    def __init__(self, handle=None):
        """
        self = Level()
        
        
        Defined at read_amr_module.fpp lines 59-73
        
        
        Returns
        -------
        this : Level
        	Object to be constructed
        
        
        Automatically generated constructor for level
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _visualisation_pkg.f90wrap_level_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Level
        
        
        Defined at read_amr_module.fpp lines 59-73
        
        Parameters
        ----------
        this : Level
        	Object to be destructed
        
        
        Automatically generated destructor for level
        """
        if self._alloc:
            _visualisation_pkg.f90wrap_level_finalise(this=self._handle)
    
    @property
    def ilevel(self):
        """
        Element ilevel ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 60
        
        """
        return _visualisation_pkg.f90wrap_level__get__ilevel(self._handle)
    
    @ilevel.setter
    def ilevel(self, ilevel):
        _visualisation_pkg.f90wrap_level__set__ilevel(self._handle, ilevel)
    
    @property
    def ngrid(self):
        """
        Element ngrid ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 61
        
        """
        return _visualisation_pkg.f90wrap_level__get__ngrid(self._handle)
    
    @ngrid.setter
    def ngrid(self, ngrid):
        _visualisation_pkg.f90wrap_level__set__ngrid(self._handle, ngrid)
    
    @property
    def ind_grid(self):
        """
        Element ind_grid ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 62
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_level__array__ind_grid(self._handle)
        if array_handle in self._arrays:
            ind_grid = self._arrays[array_handle]
        else:
            ind_grid = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_level__array__ind_grid)
            self._arrays[array_handle] = ind_grid
        return ind_grid
    
    @ind_grid.setter
    def ind_grid(self, ind_grid):
        self.ind_grid[...] = ind_grid
    
    @property
    def real_ind(self):
        """
        Element real_ind ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 63
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_level__array__real_ind(self._handle)
        if array_handle in self._arrays:
            real_ind = self._arrays[array_handle]
        else:
            real_ind = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_level__array__real_ind)
            self._arrays[array_handle] = real_ind
        return real_ind
    
    @real_ind.setter
    def real_ind(self, real_ind):
        self.real_ind[...] = real_ind
    
    @property
    def xg(self):
        """
        Element xg ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 64
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_level__array__xg(self._handle)
        if array_handle in self._arrays:
            xg = self._arrays[array_handle]
        else:
            xg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_level__array__xg)
            self._arrays[array_handle] = xg
        return xg
    
    @xg.setter
    def xg(self, xg):
        self.xg[...] = xg
    
    @property
    def cube(self):
        """
        Element cube ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 65
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_level__array__cube(self._handle)
        if array_handle in self._arrays:
            cube = self._arrays[array_handle]
        else:
            cube = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_level__array__cube)
            self._arrays[array_handle] = cube
        return cube
    
    @cube.setter
    def cube(self, cube):
        self.cube[...] = cube
    
    @property
    def map(self):
        """
        Element map ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 66
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_level__array__map(self._handle)
        if array_handle in self._arrays:
            map = self._arrays[array_handle]
        else:
            map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_level__array__map)
            self._arrays[array_handle] = map
        return map
    
    @map.setter
    def map(self, map):
        self.map[...] = map
    
    @property
    def rho(self):
        """
        Element rho ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 67
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_level__array__rho(self._handle)
        if array_handle in self._arrays:
            rho = self._arrays[array_handle]
        else:
            rho = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_level__array__rho)
            self._arrays[array_handle] = rho
        return rho
    
    @rho.setter
    def rho(self, rho):
        self.rho[...] = rho
    
    @property
    def imin(self):
        """
        Element imin ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 68
        
        """
        return _visualisation_pkg.f90wrap_level__get__imin(self._handle)
    
    @imin.setter
    def imin(self, imin):
        _visualisation_pkg.f90wrap_level__set__imin(self._handle, imin)
    
    @property
    def imax(self):
        """
        Element imax ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 69
        
        """
        return _visualisation_pkg.f90wrap_level__get__imax(self._handle)
    
    @imax.setter
    def imax(self, imax):
        _visualisation_pkg.f90wrap_level__set__imax(self._handle, imax)
    
    @property
    def jmin(self):
        """
        Element jmin ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 70
        
        """
        return _visualisation_pkg.f90wrap_level__get__jmin(self._handle)
    
    @jmin.setter
    def jmin(self, jmin):
        _visualisation_pkg.f90wrap_level__set__jmin(self._handle, jmin)
    
    @property
    def jmax(self):
        """
        Element jmax ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 71
        
        """
        return _visualisation_pkg.f90wrap_level__get__jmax(self._handle)
    
    @jmax.setter
    def jmax(self, jmax):
        _visualisation_pkg.f90wrap_level__set__jmax(self._handle, jmax)
    
    @property
    def kmin(self):
        """
        Element kmin ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 72
        
        """
        return _visualisation_pkg.f90wrap_level__get__kmin(self._handle)
    
    @kmin.setter
    def kmin(self, kmin):
        _visualisation_pkg.f90wrap_level__set__kmin(self._handle, kmin)
    
    @property
    def kmax(self):
        """
        Element kmax ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 73
        
        """
        return _visualisation_pkg.f90wrap_level__get__kmax(self._handle)
    
    @kmax.setter
    def kmax(self, kmax):
        _visualisation_pkg.f90wrap_level__set__kmax(self._handle, kmax)
    
    def __str__(self):
        ret = ['<level>{\n']
        ret.append('    ilevel : ')
        ret.append(repr(self.ilevel))
        ret.append(',\n    ngrid : ')
        ret.append(repr(self.ngrid))
        ret.append(',\n    ind_grid : ')
        ret.append(repr(self.ind_grid))
        ret.append(',\n    real_ind : ')
        ret.append(repr(self.real_ind))
        ret.append(',\n    xg : ')
        ret.append(repr(self.xg))
        ret.append(',\n    cube : ')
        ret.append(repr(self.cube))
        ret.append(',\n    map : ')
        ret.append(repr(self.map))
        ret.append(',\n    rho : ')
        ret.append(repr(self.rho))
        ret.append(',\n    imin : ')
        ret.append(repr(self.imin))
        ret.append(',\n    imax : ')
        ret.append(repr(self.imax))
        ret.append(',\n    jmin : ')
        ret.append(repr(self.jmin))
        ret.append(',\n    jmax : ')
        ret.append(repr(self.jmax))
        ret.append(',\n    kmin : ')
        ret.append(repr(self.kmin))
        ret.append(',\n    kmax : ')
        ret.append(repr(self.kmax))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("visualisation_pkg.data_handler")
class data_handler(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=data_handler)
    
    
    Defined at read_amr_module.fpp lines 75-79
    
    """
    def __init__(self, handle=None):
        """
        self = Data_Handler()
        
        
        Defined at read_amr_module.fpp lines 75-79
        
        
        Returns
        -------
        this : Data_Handler
        	Object to be constructed
        
        
        Automatically generated constructor for data_handler
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _visualisation_pkg.f90wrap_data_handler_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Data_Handler
        
        
        Defined at read_amr_module.fpp lines 75-79
        
        Parameters
        ----------
        this : Data_Handler
        	Object to be destructed
        
        
        Automatically generated destructor for data_handler
        """
        if self._alloc:
            _visualisation_pkg.f90wrap_data_handler_finalise(this=self._handle)
    
    @property
    def name(self):
        """
        Element name ftype=character(80) pytype=str
        
        
        Defined at read_amr_module.fpp line 76
        
        """
        return _visualisation_pkg.f90wrap_data_handler__get__name(self._handle)
    
    @name.setter
    def name(self, name):
        _visualisation_pkg.f90wrap_data_handler__set__name(self._handle, name)
    
    @property
    def x_data(self):
        """
        Element x_data ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 77
        
        """
        return _visualisation_pkg.f90wrap_data_handler__get__x_data(self._handle)
    
    @x_data.setter
    def x_data(self, x_data):
        _visualisation_pkg.f90wrap_data_handler__set__x_data(self._handle, x_data)
    
    @property
    def y_data(self):
        """
        Element y_data ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 77
        
        """
        return _visualisation_pkg.f90wrap_data_handler__get__y_data(self._handle)
    
    @y_data.setter
    def y_data(self, y_data):
        _visualisation_pkg.f90wrap_data_handler__set__y_data(self._handle, y_data)
    
    @property
    def z_data(self):
        """
        Element z_data ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 77
        
        """
        return _visualisation_pkg.f90wrap_data_handler__get__z_data(self._handle)
    
    @z_data.setter
    def z_data(self, z_data):
        _visualisation_pkg.f90wrap_data_handler__set__z_data(self._handle, z_data)
    
    @property
    def nx(self):
        """
        Element nx ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 78
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_data_handler__array__nx(self._handle)
        if array_handle in self._arrays:
            nx = self._arrays[array_handle]
        else:
            nx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_data_handler__array__nx)
            self._arrays[array_handle] = nx
        return nx
    
    @nx.setter
    def nx(self, nx):
        self.nx[...] = nx
    
    @property
    def ny(self):
        """
        Element ny ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 78
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_data_handler__array__ny(self._handle)
        if array_handle in self._arrays:
            ny = self._arrays[array_handle]
        else:
            ny = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_data_handler__array__ny)
            self._arrays[array_handle] = ny
        return ny
    
    @ny.setter
    def ny(self, ny):
        self.ny[...] = ny
    
    @property
    def nz(self):
        """
        Element nz ftype=integer pytype=int
        
        
        Defined at read_amr_module.fpp line 78
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_data_handler__array__nz(self._handle)
        if array_handle in self._arrays:
            nz = self._arrays[array_handle]
        else:
            nz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_data_handler__array__nz)
            self._arrays[array_handle] = nz
        return nz
    
    @nz.setter
    def nz(self, nz):
        self.nz[...] = nz
    
    @property
    def x(self):
        """
        Element x ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 79
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_data_handler__array__x(self._handle)
        if array_handle in self._arrays:
            x = self._arrays[array_handle]
        else:
            x = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_data_handler__array__x)
            self._arrays[array_handle] = x
        return x
    
    @x.setter
    def x(self, x):
        self.x[...] = x
    
    @property
    def y(self):
        """
        Element y ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 79
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_data_handler__array__y(self._handle)
        if array_handle in self._arrays:
            y = self._arrays[array_handle]
        else:
            y = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_data_handler__array__y)
            self._arrays[array_handle] = y
        return y
    
    @y.setter
    def y(self, y):
        self.y[...] = y
    
    @property
    def z(self):
        """
        Element z ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 79
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_data_handler__array__z(self._handle)
        if array_handle in self._arrays:
            z = self._arrays[array_handle]
        else:
            z = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_data_handler__array__z)
            self._arrays[array_handle] = z
        return z
    
    @z.setter
    def z(self, z):
        self.z[...] = z
    
    def __str__(self):
        ret = ['<data_handler>{\n']
        ret.append('    name : ')
        ret.append(repr(self.name))
        ret.append(',\n    x_data : ')
        ret.append(repr(self.x_data))
        ret.append(',\n    y_data : ')
        ret.append(repr(self.y_data))
        ret.append(',\n    z_data : ')
        ret.append(repr(self.z_data))
        ret.append(',\n    nx : ')
        ret.append(repr(self.nx))
        ret.append(',\n    ny : ')
        ret.append(repr(self.ny))
        ret.append(',\n    nz : ')
        ret.append(repr(self.nz))
        ret.append(',\n    x : ')
        ret.append(repr(self.x))
        ret.append(',\n    y : ')
        ret.append(repr(self.y))
        ret.append(',\n    z : ')
        ret.append(repr(self.z))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("visualisation_pkg.particle")
class particle(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=particle)
    
    
    Defined at read_amr_module.fpp lines 81-84
    
    """
    def __init__(self, handle=None):
        """
        self = Particle()
        
        
        Defined at read_amr_module.fpp lines 81-84
        
        
        Returns
        -------
        this : Particle
        	Object to be constructed
        
        
        Automatically generated constructor for particle
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _visualisation_pkg.f90wrap_particle_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Particle
        
        
        Defined at read_amr_module.fpp lines 81-84
        
        Parameters
        ----------
        this : Particle
        	Object to be destructed
        
        
        Automatically generated destructor for particle
        """
        if self._alloc:
            _visualisation_pkg.f90wrap_particle_finalise(this=self._handle)
    
    @property
    def id(self):
        """
        Element id ftype=integer(irg) pytype=int
        
        
        Defined at read_amr_module.fpp line 82
        
        """
        return _visualisation_pkg.f90wrap_particle__get__id(self._handle)
    
    @id.setter
    def id(self, id):
        _visualisation_pkg.f90wrap_particle__set__id(self._handle, id)
    
    @property
    def x(self):
        """
        Element x ftype=type(vector) pytype=Vector
        
        
        Defined at read_amr_module.fpp line 83
        
        """
        x_handle = _visualisation_pkg.f90wrap_particle__get__x(self._handle)
        if tuple(x_handle) in self._objs:
            x = self._objs[tuple(x_handle)]
        else:
            x = vector.from_handle(x_handle)
            self._objs[tuple(x_handle)] = x
        return x
    
    @x.setter
    def x(self, x):
        x = x._handle
        _visualisation_pkg.f90wrap_particle__set__x(self._handle, x)
    
    @property
    def v(self):
        """
        Element v ftype=type(vector) pytype=Vector
        
        
        Defined at read_amr_module.fpp line 83
        
        """
        v_handle = _visualisation_pkg.f90wrap_particle__get__v(self._handle)
        if tuple(v_handle) in self._objs:
            v = self._objs[tuple(v_handle)]
        else:
            v = vector.from_handle(v_handle)
            self._objs[tuple(v_handle)] = v
        return v
    
    @v.setter
    def v(self, v):
        v = v._handle
        _visualisation_pkg.f90wrap_particle__set__v(self._handle, v)
    
    @property
    def m(self):
        """
        Element m ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 84
        
        """
        return _visualisation_pkg.f90wrap_particle__get__m(self._handle)
    
    @m.setter
    def m(self, m):
        _visualisation_pkg.f90wrap_particle__set__m(self._handle, m)
    
    @property
    def met(self):
        """
        Element met ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 84
        
        """
        return _visualisation_pkg.f90wrap_particle__get__met(self._handle)
    
    @met.setter
    def met(self, met):
        _visualisation_pkg.f90wrap_particle__set__met(self._handle, met)
    
    @property
    def imass(self):
        """
        Element imass ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 84
        
        """
        return _visualisation_pkg.f90wrap_particle__get__imass(self._handle)
    
    @imass.setter
    def imass(self, imass):
        _visualisation_pkg.f90wrap_particle__set__imass(self._handle, imass)
    
    @property
    def age(self):
        """
        Element age ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 84
        
        """
        return _visualisation_pkg.f90wrap_particle__get__age(self._handle)
    
    @age.setter
    def age(self, age):
        _visualisation_pkg.f90wrap_particle__set__age(self._handle, age)
    
    @property
    def tform(self):
        """
        Element tform ftype=real(dbl) pytype=float
        
        
        Defined at read_amr_module.fpp line 84
        
        """
        return _visualisation_pkg.f90wrap_particle__get__tform(self._handle)
    
    @tform.setter
    def tform(self, tform):
        _visualisation_pkg.f90wrap_particle__set__tform(self._handle, tform)
    
    def __str__(self):
        ret = ['<particle>{\n']
        ret.append('    id : ')
        ret.append(repr(self.id))
        ret.append(',\n    x : ')
        ret.append(repr(self.x))
        ret.append(',\n    v : ')
        ret.append(repr(self.v))
        ret.append(',\n    m : ')
        ret.append(repr(self.m))
        ret.append(',\n    met : ')
        ret.append(repr(self.met))
        ret.append(',\n    imass : ')
        ret.append(repr(self.imass))
        ret.append(',\n    age : ')
        ret.append(repr(self.age))
        ret.append(',\n    tform : ')
        ret.append(repr(self.tform))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def retrieve_vars(repository, myvars):
    """
    retrieve_vars(repository, myvars)
    
    
    Defined at read_amr_module.fpp lines 91-96
    
    Parameters
    ----------
    repository : str
    myvars : Hydroid
    
    """
    _visualisation_pkg.f90wrap_retrieve_vars(repository=repository, \
        myvars=myvars._handle)

def title(n, nchar):
    """
    title(n, nchar)
    
    
    Defined at read_amr_module.fpp lines 103-127
    
    Parameters
    ----------
    n : int
    nchar : str
    
    """
    _visualisation_pkg.f90wrap_title(n=n, nchar=nchar)

def hilbert3d(x, y, z, order, bit_length, npoint):
    """
    hilbert3d(x, y, z, order, bit_length, npoint)
    
    
    Defined at read_amr_module.fpp lines 134-206
    
    Parameters
    ----------
    x : int array
    y : int array
    z : int array
    order : float array
    bit_length : int
    npoint : int
    
    """
    _visualisation_pkg.f90wrap_hilbert3d(x=x, y=y, z=z, order=order, \
        bit_length=bit_length, npoint=npoint)

def check_lmax(ngridfile):
    """
    check_lmax(ngridfile)
    
    
    Defined at read_amr_module.fpp lines 214-226
    
    Parameters
    ----------
    ngridfile : int array
    
    """
    _visualisation_pkg.f90wrap_check_lmax(ngridfile=ngridfile)

def check_families(repository):
    """
    check_families(repository)
    
    
    Defined at read_amr_module.fpp lines 235-260
    
    Parameters
    ----------
    repository : str
    
    """
    _visualisation_pkg.f90wrap_check_families(repository=repository)

def read_hydrofile_descriptor(repository):
    """
    read_hydrofile_descriptor(repository)
    
    
    Defined at read_amr_module.fpp lines 269-303
    
    Parameters
    ----------
    repository : str
    
    """
    _visualisation_pkg.f90wrap_read_hydrofile_descriptor(repository=repository)

def read_hydrofile_descriptor_old(repository):
    """
    read_hydrofile_descriptor_old(repository)
    
    
    Defined at read_amr_module.fpp lines 313-369
    
    Parameters
    ----------
    repository : str
    
    """
    _visualisation_pkg.f90wrap_read_hydrofile_descriptor_old(repository=repository)

def select_from_descriptor_ids(newvar, newid):
    """
    select_from_descriptor_ids(newvar, newid)
    
    
    Defined at read_amr_module.fpp lines 371-472
    
    Parameters
    ----------
    newvar : str
    newid : int
    
    """
    _visualisation_pkg.f90wrap_select_from_descriptor_ids(newvar=newvar, \
        newid=newid)

def read_hydrofile_descriptor_new(repository):
    """
    read_hydrofile_descriptor_new(repository)
    
    
    Defined at read_amr_module.fpp lines 482-513
    
    Parameters
    ----------
    repository : str
    
    """
    _visualisation_pkg.f90wrap_read_hydrofile_descriptor_new(repository=repository)

def getvarvalue(self, dx, x, var, son, varname, value, trans_matrix=None, \
    grav_var=None):
    """
    getvarvalue(self, dx, x, var, son, varname, value[, trans_matrix, grav_var])
    
    
    Defined at read_amr_module.fpp lines 522-1166
    
    Parameters
    ----------
    reg : Region
    dx : float
    x : Vector
    var : float array
    son : int array
    varname : str
    value : float
    trans_matrix : float array
    grav_var : float array
    
    """
    _visualisation_pkg.f90wrap_getvarvalue(reg=self._handle, dx=dx, x=x._handle, \
        var=var, son=son, varname=varname, value=value, trans_matrix=trans_matrix, \
        grav_var=grav_var)

def get_eta_sn(repository):
    """
    get_eta_sn(repository)
    
    
    Defined at read_amr_module.fpp lines 1176-1206
    
    Parameters
    ----------
    repository : str
    
    """
    _visualisation_pkg.f90wrap_get_eta_sn(repository=repository)

def sf_eff(self, dx, x, var, star_maker):
    """
    sf_eff = sf_eff(self, dx, x, var, star_maker)
    
    
    Defined at read_amr_module.fpp lines 1403-1772
    
    Parameters
    ----------
    reg : Region
    dx : float
    x : Vector
    var : float array
    star_maker : str
    
    Returns
    -------
    sf_eff : float
    
    """
    sf_eff = _visualisation_pkg.f90wrap_sf_eff(reg=self._handle, dx=dx, x=x._handle, \
        var=var, star_maker=star_maker)
    return sf_eff

def init_amr_read(repository):
    """
    init_amr_read(repository)
    
    
    Defined at read_amr_module.fpp lines 1781-1883
    
    Parameters
    ----------
    repository : str
    
    """
    _visualisation_pkg.f90wrap_init_amr_read(repository=repository)

def get_cpu_map(self):
    """
    get_cpu_map(self)
    
    
    Defined at read_amr_module.fpp lines 1891-1979
    
    Parameters
    ----------
    reg : Region
    
    """
    _visualisation_pkg.f90wrap_get_cpu_map(reg=self._handle)

def getnborgrids(son, nbor, igrid, igridn, ngrid):
    """
    getnborgrids(son, nbor, igrid, igridn, ngrid)
    
    
    Defined at read_amr_module.fpp lines 1981-2004
    
    Parameters
    ----------
    son : int array
    nbor : int array
    igrid : int array
    igridn : int array
    ngrid : int
    
    ---------------------------------------------------------
     This routine computes the index of the 6 neighboring
     grids for grid igrid(:). The index for the central
     grid is stored in igridn(:,0). If for some reasons
     the neighboring grids don't exist, then igridn(:,j) = 0.
    ---------------------------------------------------------
    """
    _visualisation_pkg.f90wrap_getnborgrids(son=son, nbor=nbor, igrid=igrid, \
        igridn=igridn, ngrid=ngrid)

def getnborcells(igridn, ind, icelln, ng):
    """
    getnborcells(igridn, ind, icelln, ng)
    
    
    Defined at read_amr_module.fpp lines 2006-2037
    
    Parameters
    ----------
    igridn : int array
    ind : int
    icelln : int array
    ng : int
    
    --------------------------------------------------------------
     This routine computes the index of 6-neighboring cells
     The user must provide igridn = index of the 6 neighboring
     grids and the cell's grid(see routine getnborgrids).
     ind is the cell index in the grid.
    --------------------------------------------------------------
    """
    _visualisation_pkg.f90wrap_getnborcells(igridn=igridn, ind=ind, icelln=icelln, \
        ng=ng)

def getnbor(son, nbor, ind_cell, ind_father, ncell):
    """
    getnbor(son, nbor, ind_cell, ind_father, ncell)
    
    
    Defined at read_amr_module.fpp lines 2039-2101
    
    Parameters
    ----------
    son : int array
    nbor : int array
    ind_cell : int array
    ind_father : int array
    ncell : int
    
    -----------------------------------------------------------------
     This subroutine determines the 2*ndim neighboring cells
     cells of the input cell(ind_cell).
     If for some reasons they don't exist, the routine returns
     the input cell.
    -----------------------------------------------------------------
    """
    _visualisation_pkg.f90wrap_getnbor(son=son, nbor=nbor, ind_cell=ind_cell, \
        ind_father=ind_father, ncell=ncell)

def getparttype(self, ptype):
    """
    getparttype(self, ptype)
    
    
    Defined at read_amr_module.fpp lines 2103-2111
    
    Parameters
    ----------
    part : Particle
    ptype : str
    
    """
    _visualisation_pkg.f90wrap_getparttype(part=self._handle, ptype=ptype)

def getpartvalue(self, part, var, value, dx=None):
    """
    getpartvalue(self, part, var, value[, dx])
    
    
    Defined at read_amr_module.fpp lines 2113-2513
    
    Parameters
    ----------
    reg : Region
    part : Particle
    var : str
    value : float
    dx : Vector
    
    """
    _visualisation_pkg.f90wrap_getpartvalue(reg=self._handle, part=part._handle, \
        var=var, value=value, dx=None if dx is None else dx._handle)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "io_ramses".')

for func in _dt_array_initialisers:
    func()
