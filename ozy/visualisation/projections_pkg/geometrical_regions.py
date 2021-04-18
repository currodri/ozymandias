"""
Module geometrical_regions


Defined at coordinates_module.fpp lines 191-420

"""
from __future__ import print_function, absolute_import, division
import _projections_pkg
import f90wrap.runtime
import logging
from projections_pkg.vectors import vector

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("projections_pkg.region")
class region(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=region)
    
    
    Defined at coordinates_module.fpp lines 195-201
    
    """
    def __init__(self, handle=None):
        """
        self = Region()
        
        
        Defined at coordinates_module.fpp lines 195-201
        
        
        Returns
        -------
        this : Region
        	Object to be constructed
        
        
        Automatically generated constructor for region
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _projections_pkg.f90wrap_region_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Region
        
        
        Defined at coordinates_module.fpp lines 195-201
        
        Parameters
        ----------
        this : Region
        	Object to be destructed
        
        
        Automatically generated destructor for region
        """
        if self._alloc:
            _projections_pkg.f90wrap_region_finalise(this=self._handle)
    
    @property
    def name(self):
        """
        Element name ftype=character(128) pytype=str
        
        
        Defined at coordinates_module.fpp line 196
        
        """
        return _projections_pkg.f90wrap_region__get__name(self._handle)
    
    @name.setter
    def name(self, name):
        _projections_pkg.f90wrap_region__set__name(self._handle, name)
    
    @property
    def centre(self):
        """
        Element centre ftype=type(vector) pytype=Vector
        
        
        Defined at coordinates_module.fpp line 197
        
        """
        centre_handle = _projections_pkg.f90wrap_region__get__centre(self._handle)
        if tuple(centre_handle) in self._objs:
            centre = self._objs[tuple(centre_handle)]
        else:
            centre = vector.from_handle(centre_handle)
            self._objs[tuple(centre_handle)] = centre
        return centre
    
    @centre.setter
    def centre(self, centre):
        centre = centre._handle
        _projections_pkg.f90wrap_region__set__centre(self._handle, centre)
    
    @property
    def axis(self):
        """
        Element axis ftype=type(vector) pytype=Vector
        
        
        Defined at coordinates_module.fpp line 197
        
        """
        axis_handle = _projections_pkg.f90wrap_region__get__axis(self._handle)
        if tuple(axis_handle) in self._objs:
            axis = self._objs[tuple(axis_handle)]
        else:
            axis = vector.from_handle(axis_handle)
            self._objs[tuple(axis_handle)] = axis
        return axis
    
    @axis.setter
    def axis(self, axis):
        axis = axis._handle
        _projections_pkg.f90wrap_region__set__axis(self._handle, axis)
    
    @property
    def bulk_velocity(self):
        """
        Element bulk_velocity ftype=type(vector) pytype=Vector
        
        
        Defined at coordinates_module.fpp line 197
        
        """
        bulk_velocity_handle = \
            _projections_pkg.f90wrap_region__get__bulk_velocity(self._handle)
        if tuple(bulk_velocity_handle) in self._objs:
            bulk_velocity = self._objs[tuple(bulk_velocity_handle)]
        else:
            bulk_velocity = vector.from_handle(bulk_velocity_handle)
            self._objs[tuple(bulk_velocity_handle)] = bulk_velocity
        return bulk_velocity
    
    @bulk_velocity.setter
    def bulk_velocity(self, bulk_velocity):
        bulk_velocity = bulk_velocity._handle
        _projections_pkg.f90wrap_region__set__bulk_velocity(self._handle, bulk_velocity)
    
    @property
    def xmin(self):
        """
        Element xmin ftype=real(dbl) pytype=float
        
        
        Defined at coordinates_module.fpp line 198
        
        """
        return _projections_pkg.f90wrap_region__get__xmin(self._handle)
    
    @xmin.setter
    def xmin(self, xmin):
        _projections_pkg.f90wrap_region__set__xmin(self._handle, xmin)
    
    @property
    def xmax(self):
        """
        Element xmax ftype=real(dbl) pytype=float
        
        
        Defined at coordinates_module.fpp line 198
        
        """
        return _projections_pkg.f90wrap_region__get__xmax(self._handle)
    
    @xmax.setter
    def xmax(self, xmax):
        _projections_pkg.f90wrap_region__set__xmax(self._handle, xmax)
    
    @property
    def ymin(self):
        """
        Element ymin ftype=real(dbl) pytype=float
        
        
        Defined at coordinates_module.fpp line 198
        
        """
        return _projections_pkg.f90wrap_region__get__ymin(self._handle)
    
    @ymin.setter
    def ymin(self, ymin):
        _projections_pkg.f90wrap_region__set__ymin(self._handle, ymin)
    
    @property
    def ymax(self):
        """
        Element ymax ftype=real(dbl) pytype=float
        
        
        Defined at coordinates_module.fpp line 198
        
        """
        return _projections_pkg.f90wrap_region__get__ymax(self._handle)
    
    @ymax.setter
    def ymax(self, ymax):
        _projections_pkg.f90wrap_region__set__ymax(self._handle, ymax)
    
    @property
    def zmin(self):
        """
        Element zmin ftype=real(dbl) pytype=float
        
        
        Defined at coordinates_module.fpp line 198
        
        """
        return _projections_pkg.f90wrap_region__get__zmin(self._handle)
    
    @zmin.setter
    def zmin(self, zmin):
        _projections_pkg.f90wrap_region__set__zmin(self._handle, zmin)
    
    @property
    def zmax(self):
        """
        Element zmax ftype=real(dbl) pytype=float
        
        
        Defined at coordinates_module.fpp line 198
        
        """
        return _projections_pkg.f90wrap_region__get__zmax(self._handle)
    
    @zmax.setter
    def zmax(self, zmax):
        _projections_pkg.f90wrap_region__set__zmax(self._handle, zmax)
    
    @property
    def rmin(self):
        """
        Element rmin ftype=real(dbl) pytype=float
        
        
        Defined at coordinates_module.fpp line 199
        
        """
        return _projections_pkg.f90wrap_region__get__rmin(self._handle)
    
    @rmin.setter
    def rmin(self, rmin):
        _projections_pkg.f90wrap_region__set__rmin(self._handle, rmin)
    
    @property
    def rmax(self):
        """
        Element rmax ftype=real(dbl) pytype=float
        
        
        Defined at coordinates_module.fpp line 199
        
        """
        return _projections_pkg.f90wrap_region__get__rmax(self._handle)
    
    @rmax.setter
    def rmax(self, rmax):
        _projections_pkg.f90wrap_region__set__rmax(self._handle, rmax)
    
    @property
    def angle(self):
        """
        Element angle ftype=real(dbl) pytype=float
        
        
        Defined at coordinates_module.fpp line 200
        
        """
        return _projections_pkg.f90wrap_region__get__angle(self._handle)
    
    @angle.setter
    def angle(self, angle):
        _projections_pkg.f90wrap_region__set__angle(self._handle, angle)
    
    @property
    def criteria_name(self):
        """
        Element criteria_name ftype=character(128) pytype=str
        
        
        Defined at coordinates_module.fpp line 201
        
        """
        return _projections_pkg.f90wrap_region__get__criteria_name(self._handle)
    
    @criteria_name.setter
    def criteria_name(self, criteria_name):
        _projections_pkg.f90wrap_region__set__criteria_name(self._handle, criteria_name)
    
    def __str__(self):
        ret = ['<region>{\n']
        ret.append('    name : ')
        ret.append(repr(self.name))
        ret.append(',\n    centre : ')
        ret.append(repr(self.centre))
        ret.append(',\n    axis : ')
        ret.append(repr(self.axis))
        ret.append(',\n    bulk_velocity : ')
        ret.append(repr(self.bulk_velocity))
        ret.append(',\n    xmin : ')
        ret.append(repr(self.xmin))
        ret.append(',\n    xmax : ')
        ret.append(repr(self.xmax))
        ret.append(',\n    ymin : ')
        ret.append(repr(self.ymin))
        ret.append(',\n    ymax : ')
        ret.append(repr(self.ymax))
        ret.append(',\n    zmin : ')
        ret.append(repr(self.zmin))
        ret.append(',\n    zmax : ')
        ret.append(repr(self.zmax))
        ret.append(',\n    rmin : ')
        ret.append(repr(self.rmin))
        ret.append(',\n    rmax : ')
        ret.append(repr(self.rmax))
        ret.append(',\n    angle : ')
        ret.append(repr(self.angle))
        ret.append(',\n    criteria_name : ')
        ret.append(repr(self.criteria_name))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def limits(self, lim):
    """
    limits(self, lim)
    
    
    Defined at coordinates_module.fpp lines 212-342
    
    Parameters
    ----------
    reg : Region
    lim : float array
    
    """
    _projections_pkg.f90wrap_limits(reg=self._handle, lim=lim)

def checkifinside(pos, reg, ok, distance):
    """
    checkifinside(pos, reg, ok, distance)
    
    
    Defined at coordinates_module.fpp lines 354-374
    
    Parameters
    ----------
    pos : float array
    reg : Region
    ok : bool
    distance : float
    
    """
    _projections_pkg.f90wrap_checkifinside(pos=pos, reg=reg._handle, ok=ok, \
        distance=distance)

def cube(self, reg, ok, distance):
    """
    cube(self, reg, ok, distance)
    
    
    Defined at coordinates_module.fpp lines 376-389
    
    Parameters
    ----------
    p : Vector
    reg : Region
    ok : bool
    distance : float
    
    """
    _projections_pkg.f90wrap_cube(p=self._handle, reg=reg._handle, ok=ok, \
        distance=distance)

def sphere(self, reg, ok, distance):
    """
    sphere(self, reg, ok, distance)
    
    
    Defined at coordinates_module.fpp lines 391-398
    
    Parameters
    ----------
    p : Vector
    reg : Region
    ok : bool
    distance : float
    
    """
    _projections_pkg.f90wrap_sphere(p=self._handle, reg=reg._handle, ok=ok, \
        distance=distance)

def cylinder(self, reg, ok, distance):
    """
    cylinder(self, reg, ok, distance)
    
    
    Defined at coordinates_module.fpp lines 400-408
    
    Parameters
    ----------
    p : Vector
    reg : Region
    ok : bool
    distance : float
    
    """
    _projections_pkg.f90wrap_cylinder(p=self._handle, reg=reg._handle, ok=ok, \
        distance=distance)

def cone(self, reg, ok, distance):
    """
    cone(self, reg, ok, distance)
    
    
    Defined at coordinates_module.fpp lines 410-420
    
    Parameters
    ----------
    p : Vector
    reg : Region
    ok : bool
    distance : float
    
    """
    _projections_pkg.f90wrap_cone(p=self._handle, reg=reg._handle, ok=ok, \
        distance=distance)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "geometrical_regions".')

for func in _dt_array_initialisers:
    func()
