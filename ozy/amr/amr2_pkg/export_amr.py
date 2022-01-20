"""
Module export_amr


Defined at export_module.fpp lines 22-699

"""
from __future__ import print_function, absolute_import, division
import _amr2_pkg
import f90wrap.runtime
import logging
from amr2_pkg.filtering import filter

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("amr2_pkg.chunk_handler")
class chunk_handler(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=chunk_handler)
    
    
    Defined at export_module.fpp lines 26-30
    
    """
    def __init__(self, handle=None):
        """
        self = Chunk_Handler()
        
        
        Defined at export_module.fpp lines 26-30
        
        
        Returns
        -------
        this : Chunk_Handler
        	Object to be constructed
        
        
        Automatically generated constructor for chunk_handler
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _amr2_pkg.f90wrap_chunk_handler_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Chunk_Handler
        
        
        Defined at export_module.fpp lines 26-30
        
        Parameters
        ----------
        this : Chunk_Handler
        	Object to be destructed
        
        
        Automatically generated destructor for chunk_handler
        """
        if self._alloc:
            _amr2_pkg.f90wrap_chunk_handler_finalise(this=self._handle)
    
    @property
    def nvars(self):
        """
        Element nvars ftype=integer  pytype=int
        
        
        Defined at export_module.fpp line 27
        
        """
        return _amr2_pkg.f90wrap_chunk_handler__get__nvars(self._handle)
    
    @nvars.setter
    def nvars(self, nvars):
        _amr2_pkg.f90wrap_chunk_handler__set__nvars(self._handle, nvars)
    
    @property
    def nx(self):
        """
        Element nx ftype=integer  pytype=int
        
        
        Defined at export_module.fpp line 27
        
        """
        return _amr2_pkg.f90wrap_chunk_handler__get__nx(self._handle)
    
    @nx.setter
    def nx(self, nx):
        _amr2_pkg.f90wrap_chunk_handler__set__nx(self._handle, nx)
    
    @property
    def ny(self):
        """
        Element ny ftype=integer  pytype=int
        
        
        Defined at export_module.fpp line 27
        
        """
        return _amr2_pkg.f90wrap_chunk_handler__get__ny(self._handle)
    
    @ny.setter
    def ny(self, ny):
        _amr2_pkg.f90wrap_chunk_handler__set__ny(self._handle, ny)
    
    @property
    def nz(self):
        """
        Element nz ftype=integer  pytype=int
        
        
        Defined at export_module.fpp line 27
        
        """
        return _amr2_pkg.f90wrap_chunk_handler__get__nz(self._handle)
    
    @nz.setter
    def nz(self, nz):
        _amr2_pkg.f90wrap_chunk_handler__set__nz(self._handle, nz)
    
    @property
    def varnames(self):
        """
        Element varnames ftype=character(128) pytype=str
        
        
        Defined at export_module.fpp line 28
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_chunk_handler__array__varnames(self._handle)
        if array_handle in self._arrays:
            varnames = self._arrays[array_handle]
        else:
            varnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_chunk_handler__array__varnames)
            self._arrays[array_handle] = varnames
        return varnames
    
    @varnames.setter
    def varnames(self, varnames):
        self.varnames[...] = varnames
    
    @property
    def filt(self):
        """
        Element filt ftype=type(filter) pytype=Filter
        
        
        Defined at export_module.fpp line 29
        
        """
        filt_handle = _amr2_pkg.f90wrap_chunk_handler__get__filt(self._handle)
        if tuple(filt_handle) in self._objs:
            filt = self._objs[tuple(filt_handle)]
        else:
            filt = filter.from_handle(filt_handle)
            self._objs[tuple(filt_handle)] = filt
        return filt
    
    @filt.setter
    def filt(self, filt):
        filt = filt._handle
        _amr2_pkg.f90wrap_chunk_handler__set__filt(self._handle, filt)
    
    @property
    def data(self):
        """
        Element data ftype=real(dbl) pytype=float
        
        
        Defined at export_module.fpp line 30
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _amr2_pkg.f90wrap_chunk_handler__array__data(self._handle)
        if array_handle in self._arrays:
            data = self._arrays[array_handle]
        else:
            data = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _amr2_pkg.f90wrap_chunk_handler__array__data)
            self._arrays[array_handle] = data
        return data
    
    @data.setter
    def data(self, data):
        self.data[...] = data
    
    def __str__(self):
        ret = ['<chunk_handler>{\n']
        ret.append('    nvars : ')
        ret.append(repr(self.nvars))
        ret.append(',\n    nx : ')
        ret.append(repr(self.nx))
        ret.append(',\n    ny : ')
        ret.append(repr(self.ny))
        ret.append(',\n    nz : ')
        ret.append(repr(self.nz))
        ret.append(',\n    varnames : ')
        ret.append(repr(self.varnames))
        ret.append(',\n    filt : ')
        ret.append(repr(self.filt))
        ret.append(',\n    data : ')
        ret.append(repr(self.data))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def allocate_chunk_handler(self):
    """
    allocate_chunk_handler(self)
    
    
    Defined at export_module.fpp lines 33-37
    
    Parameters
    ----------
    chunk : Chunk_Handler
    
    """
    _amr2_pkg.f90wrap_allocate_chunk_handler(chunk=self._handle)

def get_unigrid_old(repository, reg, chunk):
    """
    get_unigrid_old(repository, reg, chunk)
    
    
    Defined at export_module.fpp lines 39-357
    
    Parameters
    ----------
    repository : str
    reg : Region
    chunk : Chunk_Handler
    
    """
    _amr2_pkg.f90wrap_get_unigrid_old(repository=repository, reg=reg._handle, \
        chunk=chunk._handle)

def get_unigrid(repository, reg, lmax, symlog, chunk):
    """
    get_unigrid(repository, reg, lmax, symlog, chunk)
    
    
    Defined at export_module.fpp lines 359-699
    
    Parameters
    ----------
    repository : str
    reg : Region
    lmax : int
    symlog : bool
    chunk : Chunk_Handler
    
    """
    _amr2_pkg.f90wrap_get_unigrid(repository=repository, reg=reg._handle, lmax=lmax, \
        symlog=symlog, chunk=chunk._handle)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "export_amr".')

for func in _dt_array_initialisers:
    func()
