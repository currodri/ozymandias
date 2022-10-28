"""
Module obs_instruments


Defined at ramses2map.fpp lines 5-207

"""
from __future__ import print_function, absolute_import, division
import _visualisation_pkg
import f90wrap.runtime
import logging
from visualisation_pkg.filtering import filter
from visualisation_pkg.geometrical_regions import region
from visualisation_pkg.vectors import vector

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("visualisation_pkg.camera")
class camera(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=camera)
    
    
    Defined at ramses2map.fpp lines 11-19
    
    """
    def __init__(self, handle=None):
        """
        self = Camera()
        
        
        Defined at ramses2map.fpp lines 11-19
        
        
        Returns
        -------
        this : Camera
        	Object to be constructed
        
        
        Automatically generated constructor for camera
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _visualisation_pkg.f90wrap_camera_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Camera
        
        
        Defined at ramses2map.fpp lines 11-19
        
        Parameters
        ----------
        this : Camera
        	Object to be destructed
        
        
        Automatically generated destructor for camera
        """
        if self._alloc:
            _visualisation_pkg.f90wrap_camera_finalise(this=self._handle)
    
    @property
    def centre(self):
        """
        Element centre ftype=type(vector) pytype=Vector
        
        
        Defined at ramses2map.fpp line 12
        
        """
        centre_handle = _visualisation_pkg.f90wrap_camera__get__centre(self._handle)
        if tuple(centre_handle) in self._objs:
            centre = self._objs[tuple(centre_handle)]
        else:
            centre = vector.from_handle(centre_handle)
            self._objs[tuple(centre_handle)] = centre
        return centre
    
    @centre.setter
    def centre(self, centre):
        centre = centre._handle
        _visualisation_pkg.f90wrap_camera__set__centre(self._handle, centre)
    
    @property
    def los_axis(self):
        """
        Element los_axis ftype=type(vector) pytype=Vector
        
        
        Defined at ramses2map.fpp line 12
        
        """
        los_axis_handle = _visualisation_pkg.f90wrap_camera__get__los_axis(self._handle)
        if tuple(los_axis_handle) in self._objs:
            los_axis = self._objs[tuple(los_axis_handle)]
        else:
            los_axis = vector.from_handle(los_axis_handle)
            self._objs[tuple(los_axis_handle)] = los_axis
        return los_axis
    
    @los_axis.setter
    def los_axis(self, los_axis):
        los_axis = los_axis._handle
        _visualisation_pkg.f90wrap_camera__set__los_axis(self._handle, los_axis)
    
    @property
    def up_vector(self):
        """
        Element up_vector ftype=type(vector) pytype=Vector
        
        
        Defined at ramses2map.fpp line 12
        
        """
        up_vector_handle = \
            _visualisation_pkg.f90wrap_camera__get__up_vector(self._handle)
        if tuple(up_vector_handle) in self._objs:
            up_vector = self._objs[tuple(up_vector_handle)]
        else:
            up_vector = vector.from_handle(up_vector_handle)
            self._objs[tuple(up_vector_handle)] = up_vector
        return up_vector
    
    @up_vector.setter
    def up_vector(self, up_vector):
        up_vector = up_vector._handle
        _visualisation_pkg.f90wrap_camera__set__up_vector(self._handle, up_vector)
    
    @property
    def region_axis(self):
        """
        Element region_axis ftype=type(vector) pytype=Vector
        
        
        Defined at ramses2map.fpp line 12
        
        """
        region_axis_handle = \
            _visualisation_pkg.f90wrap_camera__get__region_axis(self._handle)
        if tuple(region_axis_handle) in self._objs:
            region_axis = self._objs[tuple(region_axis_handle)]
        else:
            region_axis = vector.from_handle(region_axis_handle)
            self._objs[tuple(region_axis_handle)] = region_axis
        return region_axis
    
    @region_axis.setter
    def region_axis(self, region_axis):
        region_axis = region_axis._handle
        _visualisation_pkg.f90wrap_camera__set__region_axis(self._handle, region_axis)
    
    @property
    def region_size(self):
        """
        Element region_size ftype=real(dbl) pytype=float
        
        
        Defined at ramses2map.fpp line 13
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _visualisation_pkg.f90wrap_camera__array__region_size(self._handle)
        if array_handle in self._arrays:
            region_size = self._arrays[array_handle]
        else:
            region_size = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _visualisation_pkg.f90wrap_camera__array__region_size)
            self._arrays[array_handle] = region_size
        return region_size
    
    @region_size.setter
    def region_size(self, region_size):
        self.region_size[...] = region_size
    
    @property
    def distance(self):
        """
        Element distance ftype=real(dbl) pytype=float
        
        
        Defined at ramses2map.fpp line 14
        
        """
        return _visualisation_pkg.f90wrap_camera__get__distance(self._handle)
    
    @distance.setter
    def distance(self, distance):
        _visualisation_pkg.f90wrap_camera__set__distance(self._handle, distance)
    
    @property
    def far_cut_depth(self):
        """
        Element far_cut_depth ftype=real(dbl) pytype=float
        
        
        Defined at ramses2map.fpp line 14
        
        """
        return _visualisation_pkg.f90wrap_camera__get__far_cut_depth(self._handle)
    
    @far_cut_depth.setter
    def far_cut_depth(self, far_cut_depth):
        _visualisation_pkg.f90wrap_camera__set__far_cut_depth(self._handle, \
            far_cut_depth)
    
    @property
    def map_max_size(self):
        """
        Element map_max_size ftype=integer  pytype=int
        
        
        Defined at ramses2map.fpp line 15
        
        """
        return _visualisation_pkg.f90wrap_camera__get__map_max_size(self._handle)
    
    @map_max_size.setter
    def map_max_size(self, map_max_size):
        _visualisation_pkg.f90wrap_camera__set__map_max_size(self._handle, map_max_size)
    
    @property
    def nfilter(self):
        """
        Element nfilter ftype=integer  pytype=int
        
        
        Defined at ramses2map.fpp line 16
        
        """
        return _visualisation_pkg.f90wrap_camera__get__nfilter(self._handle)
    
    @nfilter.setter
    def nfilter(self, nfilter):
        _visualisation_pkg.f90wrap_camera__set__nfilter(self._handle, nfilter)
    
    @property
    def nsubs(self):
        """
        Element nsubs ftype=integer  pytype=int
        
        
        Defined at ramses2map.fpp line 16
        
        """
        return _visualisation_pkg.f90wrap_camera__get__nsubs(self._handle)
    
    @nsubs.setter
    def nsubs(self, nsubs):
        _visualisation_pkg.f90wrap_camera__set__nsubs(self._handle, nsubs)
    
    @property
    def lmin(self):
        """
        Element lmin ftype=integer  pytype=int
        
        
        Defined at ramses2map.fpp line 17
        
        """
        return _visualisation_pkg.f90wrap_camera__get__lmin(self._handle)
    
    @lmin.setter
    def lmin(self, lmin):
        _visualisation_pkg.f90wrap_camera__set__lmin(self._handle, lmin)
    
    @property
    def lmax(self):
        """
        Element lmax ftype=integer  pytype=int
        
        
        Defined at ramses2map.fpp line 17
        
        """
        return _visualisation_pkg.f90wrap_camera__get__lmax(self._handle)
    
    @lmax.setter
    def lmax(self, lmax):
        _visualisation_pkg.f90wrap_camera__set__lmax(self._handle, lmax)
    
    def init_array_filters(self):
        self.filters = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _visualisation_pkg.f90wrap_camera__array_getitem__filters,
                                        _visualisation_pkg.f90wrap_camera__array_setitem__filters,
                                        _visualisation_pkg.f90wrap_camera__array_len__filters,
                                        """
        Element filters ftype=type(filter) pytype=Filter
        
        
        Defined at ramses2map.fpp line 18
        
        """, filter)
        return self.filters
    
    def init_array_subs(self):
        self.subs = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _visualisation_pkg.f90wrap_camera__array_getitem__subs,
                                        _visualisation_pkg.f90wrap_camera__array_setitem__subs,
                                        _visualisation_pkg.f90wrap_camera__array_len__subs,
                                        """
        Element subs ftype=type(region) pytype=Region
        
        
        Defined at ramses2map.fpp line 19
        
        """, region)
        return self.subs
    
    def __str__(self):
        ret = ['<camera>{\n']
        ret.append('    centre : ')
        ret.append(repr(self.centre))
        ret.append(',\n    los_axis : ')
        ret.append(repr(self.los_axis))
        ret.append(',\n    up_vector : ')
        ret.append(repr(self.up_vector))
        ret.append(',\n    region_axis : ')
        ret.append(repr(self.region_axis))
        ret.append(',\n    region_size : ')
        ret.append(repr(self.region_size))
        ret.append(',\n    distance : ')
        ret.append(repr(self.distance))
        ret.append(',\n    far_cut_depth : ')
        ret.append(repr(self.far_cut_depth))
        ret.append(',\n    map_max_size : ')
        ret.append(repr(self.map_max_size))
        ret.append(',\n    nfilter : ')
        ret.append(repr(self.nfilter))
        ret.append(',\n    nsubs : ')
        ret.append(repr(self.nsubs))
        ret.append(',\n    lmin : ')
        ret.append(repr(self.lmin))
        ret.append(',\n    lmax : ')
        ret.append(repr(self.lmax))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = [init_array_filters, init_array_subs]
    

def log2(x):
    """
    log2 = log2(x)
    
    
    Defined at ramses2map.fpp lines 23-27
    
    Parameters
    ----------
    x : float
    
    Returns
    -------
    log2 : float
    
    """
    log2 = _visualisation_pkg.f90wrap_log2(x=x)
    return log2

def init_camera(self, los_axis, up_vector, region_size, region_axis, distance, \
    far_cut_depth, map_max_size, nfilter, nsubs):
    """
    init_camera = init_camera(self, los_axis, up_vector, region_size, region_axis, \
        distance, far_cut_depth, map_max_size, nfilter, nsubs)
    
    
    Defined at ramses2map.fpp lines 29-48
    
    Parameters
    ----------
    centre : Vector
    los_axis : Vector
    up_vector : Vector
    region_size : float array
    region_axis : Vector
    distance : float
    far_cut_depth : float
    map_max_size : int
    nfilter : int
    nsubs : int
    
    Returns
    -------
    init_camera : Camera
    
    """
    init_camera = _visualisation_pkg.f90wrap_init_camera(centre=self._handle, \
        los_axis=los_axis._handle, up_vector=up_vector._handle, \
        region_size=region_size, region_axis=region_axis._handle, distance=distance, \
        far_cut_depth=far_cut_depth, map_max_size=map_max_size, nfilter=nfilter, \
        nsubs=nsubs)
    init_camera = \
        f90wrap.runtime.lookup_class("visualisation_pkg.camera").from_handle(init_camera, \
        alloc=True)
    return init_camera

def get_required_resolution(self):
    """
    get_required_resolution = get_required_resolution(self)
    
    
    Defined at ramses2map.fpp lines 50-53
    
    Parameters
    ----------
    cam : Camera
    
    Returns
    -------
    get_required_resolution : int
    
    """
    get_required_resolution = \
        _visualisation_pkg.f90wrap_get_required_resolution(cam=self._handle)
    return get_required_resolution

def get_map_size(self, n_map):
    """
    get_map_size(self, n_map)
    
    
    Defined at ramses2map.fpp lines 55-67
    
    Parameters
    ----------
    cam : Camera
    n_map : int array
    
    """
    _visualisation_pkg.f90wrap_get_map_size(cam=self._handle, n_map=n_map)

def get_map_box(self, box):
    """
    get_map_box(self, box)
    
    
    Defined at ramses2map.fpp lines 69-79
    
    Parameters
    ----------
    cam : Camera
    box : Region
    
    """
    _visualisation_pkg.f90wrap_get_map_box(cam=self._handle, box=box._handle)

def get_camera_basis(self, cam_basis):
    """
    get_camera_basis(self, cam_basis)
    
    
    Defined at ramses2map.fpp lines 101-109
    
    Parameters
    ----------
    cam : Camera
    cam_basis : Basis
    
    """
    _visualisation_pkg.f90wrap_get_camera_basis(cam=self._handle, \
        cam_basis=cam_basis._handle)

def los_transformation(self, trans_matrix):
    """
    los_transformation(self, trans_matrix)
    
    
    Defined at ramses2map.fpp lines 111-122
    
    Parameters
    ----------
    cam : Camera
    trans_matrix : float array
    
    """
    _visualisation_pkg.f90wrap_los_transformation(cam=self._handle, \
        trans_matrix=trans_matrix)

def get_bounding_box(self, bbox):
    """
    get_bounding_box(self, bbox)
    
    
    Defined at ramses2map.fpp lines 124-168
    
    Parameters
    ----------
    cam : Camera
    bbox : Region
    
    """
    _visualisation_pkg.f90wrap_get_bounding_box(cam=self._handle, bbox=bbox._handle)

def deproject_points(self, npoints, points):
    """
    deproject_points(self, npoints, points)
    
    
    Defined at ramses2map.fpp lines 170-187
    
    Parameters
    ----------
    cam : Camera
    npoints : int
    points : float array
    
    """
    _visualisation_pkg.f90wrap_deproject_points(cam=self._handle, npoints=npoints, \
        points=points)

def project_points(self, npoints, points):
    """
    project_points(self, npoints, points)
    
    
    Defined at ramses2map.fpp lines 189-206
    
    Parameters
    ----------
    cam : Camera
    npoints : int
    points : float array
    
    """
    _visualisation_pkg.f90wrap_project_points(cam=self._handle, npoints=npoints, \
        points=points)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "obs_instruments".')

for func in _dt_array_initialisers:
    func()
