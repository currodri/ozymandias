"""
Module cooling_module


Defined at cooling_module.fpp lines 24-163

"""
from __future__ import print_function, absolute_import, division
import _visualisation_pkg
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("visualisation_pkg.cooling_table")
class cooling_table(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=cooling_table)
    
    
    Defined at cooling_module.fpp lines 27-43
    
    """
    def __init__(self, handle=None):
        """
        self = Cooling_Table()
        
        
        Defined at cooling_module.fpp lines 27-43
        
        
        Returns
        -------
        this : Cooling_Table
        	Object to be constructed
        
        
        Automatically generated constructor for cooling_table
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _visualisation_pkg.f90wrap_cooling_table_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Cooling_Table
        
        
        Defined at cooling_module.fpp lines 27-43
        
        Parameters
        ----------
        this : Cooling_Table
        	Object to be destructed
        
        
        Automatically generated destructor for cooling_table
        """
        if self._alloc:
            _visualisation_pkg.f90wrap_cooling_table_finalise(this=self._handle)
    
    @property
    def n1(self):
        """
        Element n1 ftype=integer pytype=int
        
        
        Defined at cooling_module.fpp line 28
        
        """
        return _visualisation_pkg.f90wrap_cooling_table__get__n1(self._handle)
    
    @n1.setter
    def n1(self, n1):
        _visualisation_pkg.f90wrap_cooling_table__set__n1(self._handle, n1)
    
    @property
    def n2(self):
        """
        Element n2 ftype=integer pytype=int
        
        
        Defined at cooling_module.fpp line 29
        
        """
        return _visualisation_pkg.f90wrap_cooling_table__get__n2(self._handle)
    
    @n2.setter
    def n2(self, n2):
        _visualisation_pkg.f90wrap_cooling_table__set__n2(self._handle, n2)
    
    @property
    def nh(self):
        """
        Element nh ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 30
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__nh(self._handle)
        if array_handle in self._arrays:
            nh = self._arrays[array_handle]
        else:
            nh = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__nh)
            self._arrays[array_handle] = nh
        return nh
    
    @nh.setter
    def nh(self, nh):
        self.nh[...] = nh
    
    @property
    def t2(self):
        """
        Element t2 ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 31
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__t2(self._handle)
        if array_handle in self._arrays:
            t2 = self._arrays[array_handle]
        else:
            t2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__t2)
            self._arrays[array_handle] = t2
        return t2
    
    @t2.setter
    def t2(self, t2):
        self.t2[...] = t2
    
    @property
    def cool(self):
        """
        Element cool ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 32
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__cool(self._handle)
        if array_handle in self._arrays:
            cool = self._arrays[array_handle]
        else:
            cool = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__cool)
            self._arrays[array_handle] = cool
        return cool
    
    @cool.setter
    def cool(self, cool):
        self.cool[...] = cool
    
    @property
    def heat(self):
        """
        Element heat ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 33
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__heat(self._handle)
        if array_handle in self._arrays:
            heat = self._arrays[array_handle]
        else:
            heat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__heat)
            self._arrays[array_handle] = heat
        return heat
    
    @heat.setter
    def heat(self, heat):
        self.heat[...] = heat
    
    @property
    def cool_com(self):
        """
        Element cool_com ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 34
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__cool_com(self._handle)
        if array_handle in self._arrays:
            cool_com = self._arrays[array_handle]
        else:
            cool_com = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__cool_com)
            self._arrays[array_handle] = cool_com
        return cool_com
    
    @cool_com.setter
    def cool_com(self, cool_com):
        self.cool_com[...] = cool_com
    
    @property
    def heat_com(self):
        """
        Element heat_com ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 35
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__heat_com(self._handle)
        if array_handle in self._arrays:
            heat_com = self._arrays[array_handle]
        else:
            heat_com = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__heat_com)
            self._arrays[array_handle] = heat_com
        return heat_com
    
    @heat_com.setter
    def heat_com(self, heat_com):
        self.heat_com[...] = heat_com
    
    @property
    def metal(self):
        """
        Element metal ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 36
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__metal(self._handle)
        if array_handle in self._arrays:
            metal = self._arrays[array_handle]
        else:
            metal = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__metal)
            self._arrays[array_handle] = metal
        return metal
    
    @metal.setter
    def metal(self, metal):
        self.metal[...] = metal
    
    @property
    def cool_prime(self):
        """
        Element cool_prime ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 37
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__cool_prime(self._handle)
        if array_handle in self._arrays:
            cool_prime = self._arrays[array_handle]
        else:
            cool_prime = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__cool_prime)
            self._arrays[array_handle] = cool_prime
        return cool_prime
    
    @cool_prime.setter
    def cool_prime(self, cool_prime):
        self.cool_prime[...] = cool_prime
    
    @property
    def heat_prime(self):
        """
        Element heat_prime ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 38
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__heat_prime(self._handle)
        if array_handle in self._arrays:
            heat_prime = self._arrays[array_handle]
        else:
            heat_prime = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__heat_prime)
            self._arrays[array_handle] = heat_prime
        return heat_prime
    
    @heat_prime.setter
    def heat_prime(self, heat_prime):
        self.heat_prime[...] = heat_prime
    
    @property
    def cool_com_prime(self):
        """
        Element cool_com_prime ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 39
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__cool_com_prime(self._handle)
        if array_handle in self._arrays:
            cool_com_prime = self._arrays[array_handle]
        else:
            cool_com_prime = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__cool_com_prime)
            self._arrays[array_handle] = cool_com_prime
        return cool_com_prime
    
    @cool_com_prime.setter
    def cool_com_prime(self, cool_com_prime):
        self.cool_com_prime[...] = cool_com_prime
    
    @property
    def heat_com_prime(self):
        """
        Element heat_com_prime ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 40
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__heat_com_prime(self._handle)
        if array_handle in self._arrays:
            heat_com_prime = self._arrays[array_handle]
        else:
            heat_com_prime = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__heat_com_prime)
            self._arrays[array_handle] = heat_com_prime
        return heat_com_prime
    
    @heat_com_prime.setter
    def heat_com_prime(self, heat_com_prime):
        self.heat_com_prime[...] = heat_com_prime
    
    @property
    def metal_prime(self):
        """
        Element metal_prime ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 41
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__metal_prime(self._handle)
        if array_handle in self._arrays:
            metal_prime = self._arrays[array_handle]
        else:
            metal_prime = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__metal_prime)
            self._arrays[array_handle] = metal_prime
        return metal_prime
    
    @metal_prime.setter
    def metal_prime(self, metal_prime):
        self.metal_prime[...] = metal_prime
    
    @property
    def mu(self):
        """
        Element mu ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 42
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__mu(self._handle)
        if array_handle in self._arrays:
            mu = self._arrays[array_handle]
        else:
            mu = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__mu)
            self._arrays[array_handle] = mu
        return mu
    
    @mu.setter
    def mu(self, mu):
        self.mu[...] = mu
    
    @property
    def n_spec(self):
        """
        Element n_spec ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 43
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_cooling_table__array__n_spec(self._handle)
        if array_handle in self._arrays:
            n_spec = self._arrays[array_handle]
        else:
            n_spec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_cooling_table__array__n_spec)
            self._arrays[array_handle] = n_spec
        return n_spec
    
    @n_spec.setter
    def n_spec(self, n_spec):
        self.n_spec[...] = n_spec
    
    def __str__(self):
        ret = ['<cooling_table>{\n']
        ret.append('    n1 : ')
        ret.append(repr(self.n1))
        ret.append(',\n    n2 : ')
        ret.append(repr(self.n2))
        ret.append(',\n    nh : ')
        ret.append(repr(self.nh))
        ret.append(',\n    t2 : ')
        ret.append(repr(self.t2))
        ret.append(',\n    cool : ')
        ret.append(repr(self.cool))
        ret.append(',\n    heat : ')
        ret.append(repr(self.heat))
        ret.append(',\n    cool_com : ')
        ret.append(repr(self.cool_com))
        ret.append(',\n    heat_com : ')
        ret.append(repr(self.heat_com))
        ret.append(',\n    metal : ')
        ret.append(repr(self.metal))
        ret.append(',\n    cool_prime : ')
        ret.append(repr(self.cool_prime))
        ret.append(',\n    heat_prime : ')
        ret.append(repr(self.heat_prime))
        ret.append(',\n    cool_com_prime : ')
        ret.append(repr(self.cool_com_prime))
        ret.append(',\n    heat_com_prime : ')
        ret.append(repr(self.heat_com_prime))
        ret.append(',\n    metal_prime : ')
        ret.append(repr(self.metal_prime))
        ret.append(',\n    mu : ')
        ret.append(repr(self.mu))
        ret.append(',\n    n_spec : ')
        ret.append(repr(self.n_spec))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def retrieve_table(repository, mytable):
    """
    retrieve_table(repository, mytable)
    
    
    Defined at cooling_module.fpp lines 58-63
    
    Parameters
    ----------
    repository : str
    mytable : Cooling_Table
    
    """
    _visualisation_pkg.f90wrap_retrieve_table(repository=repository, \
        mytable=mytable._handle)

def read_cool(filename):
    """
    read_cool(filename)
    
    
    Defined at cooling_module.fpp lines 65-111
    
    Parameters
    ----------
    filename : str
    
    """
    _visualisation_pkg.f90wrap_read_cool(filename=filename)

def solve_cooling(nh, t2, zsolar):
    """
    lambda_, lambda_prime = solve_cooling(nh, t2, zsolar)
    
    
    Defined at cooling_module.fpp lines 113-163
    
    Parameters
    ----------
    nh : float
    t2 : float
    zsolar : float
    
    Returns
    -------
    lambda_ : float
    lambda_prime : float
    
    """
    lambda_, lambda_prime = _visualisation_pkg.f90wrap_solve_cooling(nh=nh, t2=t2, \
        zsolar=zsolar)
    return lambda_, lambda_prime

def get_if_species_abundances():
    """
    Element if_species_abundances ftype=logical pytype=bool
    
    
    Defined at cooling_module.fpp line 47
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__if_species_abundances()

if_species_abundances = get_if_species_abundances()

def get_self_shielding():
    """
    Element self_shielding ftype=logical pytype=bool
    
    
    Defined at cooling_module.fpp line 48
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__self_shielding()

self_shielding = get_self_shielding()

def get_nbin_t_fix():
    """
    Element nbin_t_fix ftype=integer pytype=int
    
    
    Defined at cooling_module.fpp line 50
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__nbin_t_fix()

nbin_T_fix = get_nbin_t_fix()

def get_nbin_n_fix():
    """
    Element nbin_n_fix ftype=integer pytype=int
    
    
    Defined at cooling_module.fpp line 51
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__nbin_n_fix()

nbin_n_fix = get_nbin_n_fix()

def get_nh_min_fix():
    """
    Element nh_min_fix ftype=real(kind=8) pytype=float
    
    
    Defined at cooling_module.fpp line 52
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__nh_min_fix()

nH_min_fix = get_nh_min_fix()

def get_nh_max_fix():
    """
    Element nh_max_fix ftype=real(kind=8) pytype=float
    
    
    Defined at cooling_module.fpp line 53
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__nh_max_fix()

nH_max_fix = get_nh_max_fix()

def get_t2_min_fix():
    """
    Element t2_min_fix ftype=real(kind=8) pytype=float
    
    
    Defined at cooling_module.fpp line 54
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__t2_min_fix()

T2_min_fix = get_t2_min_fix()

def get_t2_max_fix():
    """
    Element t2_max_fix ftype=real(kind=8) pytype=float
    
    
    Defined at cooling_module.fpp line 55
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__t2_max_fix()

T2_max_fix = get_t2_max_fix()

def get_logt2max():
    """
    Element logt2max ftype=real(kind=8) pytype=float
    
    
    Defined at cooling_module.fpp line 56
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__logt2max()

def set_logt2max(logt2max):
    _visualisation_pkg.f90wrap_cooling_module__set__logt2max(logt2max)

def get_dlog_nh():
    """
    Element dlog_nh ftype=real(kind=8) pytype=float
    
    
    Defined at cooling_module.fpp line 56
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__dlog_nh()

def set_dlog_nh(dlog_nh):
    _visualisation_pkg.f90wrap_cooling_module__set__dlog_nh(dlog_nh)

def get_dlog_t2():
    """
    Element dlog_t2 ftype=real(kind=8) pytype=float
    
    
    Defined at cooling_module.fpp line 56
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__dlog_t2()

def set_dlog_t2(dlog_t2):
    _visualisation_pkg.f90wrap_cooling_module__set__dlog_t2(dlog_t2)

def get_h():
    """
    Element h ftype=real(kind=8) pytype=float
    
    
    Defined at cooling_module.fpp line 56
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__h()

def set_h(h):
    _visualisation_pkg.f90wrap_cooling_module__set__h(h)

def get_h2():
    """
    Element h2 ftype=real(kind=8) pytype=float
    
    
    Defined at cooling_module.fpp line 56
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__h2()

def set_h2(h2):
    _visualisation_pkg.f90wrap_cooling_module__set__h2(h2)

def get_h3():
    """
    Element h3 ftype=real(kind=8) pytype=float
    
    
    Defined at cooling_module.fpp line 56
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__h3()

def set_h3(h3):
    _visualisation_pkg.f90wrap_cooling_module__set__h3(h3)

def get_precoeff():
    """
    Element precoeff ftype=real(kind=8) pytype=float
    
    
    Defined at cooling_module.fpp line 56
    
    """
    return _visualisation_pkg.f90wrap_cooling_module__get__precoeff()

def set_precoeff(precoeff):
    _visualisation_pkg.f90wrap_cooling_module__set__precoeff(precoeff)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "cooling_module".')

for func in _dt_array_initialisers:
    func()
