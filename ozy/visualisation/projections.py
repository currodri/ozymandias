from __future__ import print_function, absolute_import, division
import _projections
import f90wrap.runtime
import logging

class Vectors(f90wrap.runtime.FortranModule):
    """
    Module vectors
    
    
    Defined at linalg_module.fpp lines 23-159
    
    """
    @f90wrap.runtime.register_class("projections.vector")
    class vector(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=vector)
        
        
        Defined at linalg_module.fpp lines 26-29
        
        """
        def __init__(self, handle=None):
            """
            self = Vector()
            
            
            Defined at linalg_module.fpp lines 26-29
            
            
            Returns
            -------
            this : Vector
            	Object to be constructed
            
            
            Automatically generated constructor for vector
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _projections.f90wrap_vector_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Vector
            
            
            Defined at linalg_module.fpp lines 26-29
            
            Parameters
            ----------
            this : Vector
            	Object to be destructed
            
            
            Automatically generated destructor for vector
            """
            if self._alloc:
                _projections.f90wrap_vector_finalise(this=self._handle)
        
        @property
        def x(self):
            """
            Element x ftype=real(dbl) pytype=float
            
            
            Defined at linalg_module.fpp line 27
            
            """
            return _projections.f90wrap_vector__get__x(self._handle)
        
        @x.setter
        def x(self, x):
            _projections.f90wrap_vector__set__x(self._handle, x)
        
        @property
        def y(self):
            """
            Element y ftype=real(dbl) pytype=float
            
            
            Defined at linalg_module.fpp line 28
            
            """
            return _projections.f90wrap_vector__get__y(self._handle)
        
        @y.setter
        def y(self, y):
            _projections.f90wrap_vector__set__y(self._handle, y)
        
        @property
        def z(self):
            """
            Element z ftype=real(dbl) pytype=float
            
            
            Defined at linalg_module.fpp line 29
            
            """
            return _projections.f90wrap_vector__get__z(self._handle)
        
        @z.setter
        def z(self, z):
            _projections.f90wrap_vector__set__z(self._handle, z)
        
        def __str__(self):
            ret = ['<vector>{\n']
            ret.append('    x : ')
            ret.append(repr(self.x))
            ret.append(',\n    y : ')
            ret.append(repr(self.y))
            ret.append(',\n    z : ')
            ret.append(repr(self.z))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("projections.array_vectors")
    class array_vectors(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=array_vectors)
        
        
        Defined at linalg_module.fpp lines 31-33
        
        """
        def __init__(self, handle=None):
            """
            self = Array_Vectors()
            
            
            Defined at linalg_module.fpp lines 31-33
            
            
            Returns
            -------
            this : Array_Vectors
            	Object to be constructed
            
            
            Automatically generated constructor for array_vectors
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _projections.f90wrap_array_vectors_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Array_Vectors
            
            
            Defined at linalg_module.fpp lines 31-33
            
            Parameters
            ----------
            this : Array_Vectors
            	Object to be destructed
            
            
            Automatically generated destructor for array_vectors
            """
            if self._alloc:
                _projections.f90wrap_array_vectors_finalise(this=self._handle)
        
        @property
        def n(self):
            """
            Element n ftype=integer  pytype=int
            
            
            Defined at linalg_module.fpp line 32
            
            """
            return _projections.f90wrap_array_vectors__get__n(self._handle)
        
        @n.setter
        def n(self, n):
            _projections.f90wrap_array_vectors__set__n(self._handle, n)
        
        def init_array_list(self):
            self.list = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _projections.f90wrap_array_vectors__array_getitem__list,
                                            _projections.f90wrap_array_vectors__array_setitem__list,
                                            _projections.f90wrap_array_vectors__array_len__list,
                                            """
            Element list ftype=type(vector) pytype=Vector
            
            
            Defined at linalg_module.fpp line 33
            
            """, Vectors.vector)
            return self.list
        
        def __str__(self):
            ret = ['<array_vectors>{\n']
            ret.append('    n : ')
            ret.append(repr(self.n))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = [init_array_list]
        
    
    @staticmethod
    def magnitude(self):
        """
        magnitude = magnitude(self)
        
        
        Defined at linalg_module.fpp lines 150-152
        
        Parameters
        ----------
        vec_1 : Vector
        
        Returns
        -------
        magnitude : float
        
        """
        magnitude = _projections.f90wrap_magnitude(vec_1=self._handle)
        return magnitude
    
    @staticmethod
    def _array_to_vector(self, array):
        """
        _array_to_vector(self, array)
        
        
        Defined at linalg_module.fpp lines 64-69
        
        Parameters
        ----------
        vec_result : Vector
        array : float array
        
        """
        _projections.f90wrap_array_to_vector(vec_result=self._handle, array=array)
    
    @staticmethod
    def _vector_to_array(array_result, vec_1):
        """
        _vector_to_array(array_result, vec_1)
        
        
        Defined at linalg_module.fpp lines 71-76
        
        Parameters
        ----------
        array_result : float array
        vec_1 : Vector
        
        """
        _projections.f90wrap_vector_to_array(array_result=array_result, \
            vec_1=vec_1._handle)
    
    @staticmethod
    def assignment(*args, **kwargs):
        """
        assignment(*args, **kwargs)
        
        
        Defined at linalg_module.fpp lines 35-37
        
        Overloaded interface containing the following procedures:
          _array_to_vector
          _vector_to_array
        
        """
        for proc in [Vectors._array_to_vector, Vectors._vector_to_array]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    _dt_array_initialisers = []
    

vectors = Vectors()

class Rotations(f90wrap.runtime.FortranModule):
    """
    Module rotations
    
    
    Defined at linalg_module.fpp lines 176-249
    
    """
    @staticmethod
    def euler_matrix(r, dim, angle, ap):
        """
        euler_matrix(r, dim, angle, ap)
        
        
        Defined at linalg_module.fpp lines 193-226
        
        Parameters
        ----------
        r : float array
        dim : int
        angle : float
        ap : str
        
        """
        _projections.f90wrap_euler_matrix(r=r, dim=dim, angle=angle, ap=ap)
    
    @staticmethod
    def _rotate_vector_single(self, rotation_matrix):
        """
        _rotate_vector_single(self, rotation_matrix)
        
        
        Defined at linalg_module.fpp lines 237-240
        
        Parameters
        ----------
        vec : Vector
        rotation_matrix : float array
        
        """
        _projections.f90wrap_rotate_vector_single(vec=self._handle, \
            rotation_matrix=rotation_matrix)
    
    @staticmethod
    def _rotate_vector_array(self, rotation_matrix):
        """
        _rotate_vector_array(self, rotation_matrix)
        
        
        Defined at linalg_module.fpp lines 242-248
        
        Parameters
        ----------
        vec : Array_Vectors
        rotation_matrix : float array
        
        """
        _projections.f90wrap_rotate_vector_array(vec=self._handle, \
            rotation_matrix=rotation_matrix)
    
    @staticmethod
    def rotate_vector(*args, **kwargs):
        """
        rotate_vector(*args, **kwargs)
        
        
        Defined at linalg_module.fpp lines 180-182
        
        Overloaded interface containing the following procedures:
          _rotate_vector_single
          _rotate_vector_array
        
        """
        for proc in [Rotations._rotate_vector_single, Rotations._rotate_vector_array]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
    
    _dt_array_initialisers = []
    

rotations = Rotations()

class Basis_Representations(f90wrap.runtime.FortranModule):
    """
    Module basis_representations
    
    
    Defined at linalg_module.fpp lines 266-299
    
    """
    @f90wrap.runtime.register_class("projections.basis")
    class basis(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=basis)
        
        
        Defined at linalg_module.fpp lines 269-270
        
        """
        def __init__(self, handle=None):
            """
            self = Basis()
            
            
            Defined at linalg_module.fpp lines 269-270
            
            
            Returns
            -------
            this : Basis
            	Object to be constructed
            
            
            Automatically generated constructor for basis
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _projections.f90wrap_basis_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Basis
            
            
            Defined at linalg_module.fpp lines 269-270
            
            Parameters
            ----------
            this : Basis
            	Object to be destructed
            
            
            Automatically generated destructor for basis
            """
            if self._alloc:
                _projections.f90wrap_basis_finalise(this=self._handle)
        
        def init_array_u(self):
            self.u = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _projections.f90wrap_basis__array_getitem__u,
                                            _projections.f90wrap_basis__array_setitem__u,
                                            _projections.f90wrap_basis__array_len__u,
                                            """
            Element u ftype=type(vector) pytype=Vector
            
            
            Defined at linalg_module.fpp line 270
            
            """, Vectors.vector)
            return self.u
        
        _dt_array_initialisers = [init_array_u]
        
    
    @staticmethod
    def initialise_basis(self):
        """
        initialise_basis(self)
        
        
        Defined at linalg_module.fpp lines 273-277
        
        Parameters
        ----------
        this : Basis
        
        """
        _projections.f90wrap_initialise_basis(this=self._handle)
    
    @staticmethod
    def mgramschmidt(self, e):
        """
        mgramschmidt(self, e)
        
        
        Defined at linalg_module.fpp lines 286-298
        
        Parameters
        ----------
        vecs : Basis
        e : Basis
        
        """
        _projections.f90wrap_mgramschmidt(vecs=self._handle, e=e._handle)
    
    _dt_array_initialisers = []
    

basis_representations = Basis_Representations()

class Coordinate_Systems(f90wrap.runtime.FortranModule):
    """
    Module coordinate_systems
    
    
    Defined at coordinates_module.fpp lines 23-173
    
    """
    @staticmethod
    def r_sphere(self):
        """
        r_sphere = r_sphere(self)
        
        
        Defined at coordinates_module.fpp lines 31-33
        
        Parameters
        ----------
        p : Vector
        
        Returns
        -------
        r_sphere : float
        
        """
        r_sphere = _projections.f90wrap_r_sphere(p=self._handle)
        return r_sphere
    
    @staticmethod
    def theta_sphere(self):
        """
        theta_sphere = theta_sphere(self)
        
        
        Defined at coordinates_module.fpp lines 35-39
        
        Parameters
        ----------
        p : Vector
        
        Returns
        -------
        theta_sphere : float
        
        """
        theta_sphere = _projections.f90wrap_theta_sphere(p=self._handle)
        return theta_sphere
    
    @staticmethod
    def phi_sphere(self):
        """
        phi_sphere = phi_sphere(self)
        
        
        Defined at coordinates_module.fpp lines 41-43
        
        Parameters
        ----------
        p : Vector
        
        Returns
        -------
        phi_sphere : float
        
        """
        phi_sphere = _projections.f90wrap_phi_sphere(p=self._handle)
        return phi_sphere
    
    @staticmethod
    def r_cyl(self):
        """
        r_cyl = r_cyl(self)
        
        
        Defined at coordinates_module.fpp lines 45-47
        
        Parameters
        ----------
        p : Vector
        
        Returns
        -------
        r_cyl : float
        
        """
        r_cyl = _projections.f90wrap_r_cyl(p=self._handle)
        return r_cyl
    
    @staticmethod
    def phi_cyl(self):
        """
        phi_cyl = phi_cyl(self)
        
        
        Defined at coordinates_module.fpp lines 49-61
        
        Parameters
        ----------
        p : Vector
        
        Returns
        -------
        phi_cyl : float
        
        """
        phi_cyl = _projections.f90wrap_phi_cyl(p=self._handle)
        return phi_cyl
    
    @staticmethod
    def spherical_basis_from_cartesian(self, spher_basis):
        """
        spherical_basis_from_cartesian(self, spher_basis)
        
        
        Defined at coordinates_module.fpp lines 69-95
        
        Parameters
        ----------
        p : Vector
        spher_basis : Basis
        
        """
        _projections.f90wrap_spherical_basis_from_cartesian(p=self._handle, \
            spher_basis=spher_basis._handle)
    
    @staticmethod
    def cylindrical_basis_from_cartesian(self, cyl_basis):
        """
        cylindrical_basis_from_cartesian(self, cyl_basis)
        
        
        Defined at coordinates_module.fpp lines 103-111
        
        Parameters
        ----------
        p : Vector
        cyl_basis : Basis
        
        """
        _projections.f90wrap_cylindrical_basis_from_cartesian(p=self._handle, \
            cyl_basis=cyl_basis._handle)
    
    @staticmethod
    def new_z_coordinates(self, transformation_matrix, errormsg):
        """
        new_z_coordinates(self, transformation_matrix, errormsg)
        
        
        Defined at coordinates_module.fpp lines 123-172
        
        Parameters
        ----------
        axis : Vector
        transformation_matrix : float array
        errormsg : int
        
        """
        _projections.f90wrap_new_z_coordinates(axis=self._handle, \
            transformation_matrix=transformation_matrix, errormsg=errormsg)
    
    _dt_array_initialisers = []
    

coordinate_systems = Coordinate_Systems()

class Geometrical_Regions(f90wrap.runtime.FortranModule):
    """
    Module geometrical_regions
    
    
    Defined at coordinates_module.fpp lines 191-420
    
    """
    @f90wrap.runtime.register_class("projections.region")
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
            result = _projections.f90wrap_region_initialise()
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
                _projections.f90wrap_region_finalise(this=self._handle)
        
        @property
        def name(self):
            """
            Element name ftype=character(128) pytype=str
            
            
            Defined at coordinates_module.fpp line 196
            
            """
            return _projections.f90wrap_region__get__name(self._handle)
        
        @name.setter
        def name(self, name):
            _projections.f90wrap_region__set__name(self._handle, name)
        
        @property
        def centre(self):
            """
            Element centre ftype=type(vector) pytype=Vector
            
            
            Defined at coordinates_module.fpp line 197
            
            """
            centre_handle = _projections.f90wrap_region__get__centre(self._handle)
            if tuple(centre_handle) in self._objs:
                centre = self._objs[tuple(centre_handle)]
            else:
                centre = vectors.vector.from_handle(centre_handle)
                self._objs[tuple(centre_handle)] = centre
            return centre
        
        @centre.setter
        def centre(self, centre):
            centre = centre._handle
            _projections.f90wrap_region__set__centre(self._handle, centre)
        
        @property
        def axis(self):
            """
            Element axis ftype=type(vector) pytype=Vector
            
            
            Defined at coordinates_module.fpp line 197
            
            """
            axis_handle = _projections.f90wrap_region__get__axis(self._handle)
            if tuple(axis_handle) in self._objs:
                axis = self._objs[tuple(axis_handle)]
            else:
                axis = vectors.vector.from_handle(axis_handle)
                self._objs[tuple(axis_handle)] = axis
            return axis
        
        @axis.setter
        def axis(self, axis):
            axis = axis._handle
            _projections.f90wrap_region__set__axis(self._handle, axis)
        
        @property
        def bulk_velocity(self):
            """
            Element bulk_velocity ftype=type(vector) pytype=Vector
            
            
            Defined at coordinates_module.fpp line 197
            
            """
            bulk_velocity_handle = \
                _projections.f90wrap_region__get__bulk_velocity(self._handle)
            if tuple(bulk_velocity_handle) in self._objs:
                bulk_velocity = self._objs[tuple(bulk_velocity_handle)]
            else:
                bulk_velocity = vectors.vector.from_handle(bulk_velocity_handle)
                self._objs[tuple(bulk_velocity_handle)] = bulk_velocity
            return bulk_velocity
        
        @bulk_velocity.setter
        def bulk_velocity(self, bulk_velocity):
            bulk_velocity = bulk_velocity._handle
            _projections.f90wrap_region__set__bulk_velocity(self._handle, bulk_velocity)
        
        @property
        def xmin(self):
            """
            Element xmin ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _projections.f90wrap_region__get__xmin(self._handle)
        
        @xmin.setter
        def xmin(self, xmin):
            _projections.f90wrap_region__set__xmin(self._handle, xmin)
        
        @property
        def xmax(self):
            """
            Element xmax ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _projections.f90wrap_region__get__xmax(self._handle)
        
        @xmax.setter
        def xmax(self, xmax):
            _projections.f90wrap_region__set__xmax(self._handle, xmax)
        
        @property
        def ymin(self):
            """
            Element ymin ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _projections.f90wrap_region__get__ymin(self._handle)
        
        @ymin.setter
        def ymin(self, ymin):
            _projections.f90wrap_region__set__ymin(self._handle, ymin)
        
        @property
        def ymax(self):
            """
            Element ymax ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _projections.f90wrap_region__get__ymax(self._handle)
        
        @ymax.setter
        def ymax(self, ymax):
            _projections.f90wrap_region__set__ymax(self._handle, ymax)
        
        @property
        def zmin(self):
            """
            Element zmin ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _projections.f90wrap_region__get__zmin(self._handle)
        
        @zmin.setter
        def zmin(self, zmin):
            _projections.f90wrap_region__set__zmin(self._handle, zmin)
        
        @property
        def zmax(self):
            """
            Element zmax ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _projections.f90wrap_region__get__zmax(self._handle)
        
        @zmax.setter
        def zmax(self, zmax):
            _projections.f90wrap_region__set__zmax(self._handle, zmax)
        
        @property
        def rmin(self):
            """
            Element rmin ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 199
            
            """
            return _projections.f90wrap_region__get__rmin(self._handle)
        
        @rmin.setter
        def rmin(self, rmin):
            _projections.f90wrap_region__set__rmin(self._handle, rmin)
        
        @property
        def rmax(self):
            """
            Element rmax ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 199
            
            """
            return _projections.f90wrap_region__get__rmax(self._handle)
        
        @rmax.setter
        def rmax(self, rmax):
            _projections.f90wrap_region__set__rmax(self._handle, rmax)
        
        @property
        def angle(self):
            """
            Element angle ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 200
            
            """
            return _projections.f90wrap_region__get__angle(self._handle)
        
        @angle.setter
        def angle(self, angle):
            _projections.f90wrap_region__set__angle(self._handle, angle)
        
        @property
        def criteria_name(self):
            """
            Element criteria_name ftype=character(128) pytype=str
            
            
            Defined at coordinates_module.fpp line 201
            
            """
            return _projections.f90wrap_region__get__criteria_name(self._handle)
        
        @criteria_name.setter
        def criteria_name(self, criteria_name):
            _projections.f90wrap_region__set__criteria_name(self._handle, criteria_name)
        
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
        
    
    @staticmethod
    def limits(self, lim):
        """
        limits(self, lim)
        
        
        Defined at coordinates_module.fpp lines 212-342
        
        Parameters
        ----------
        reg : Region
        lim : float array
        
        """
        _projections.f90wrap_limits(reg=self._handle, lim=lim)
    
    @staticmethod
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
        _projections.f90wrap_checkifinside(pos=pos, reg=reg._handle, ok=ok, \
            distance=distance)
    
    @staticmethod
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
        _projections.f90wrap_cube(p=self._handle, reg=reg._handle, ok=ok, \
            distance=distance)
    
    @staticmethod
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
        _projections.f90wrap_sphere(p=self._handle, reg=reg._handle, ok=ok, \
            distance=distance)
    
    @staticmethod
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
        _projections.f90wrap_cylinder(p=self._handle, reg=reg._handle, ok=ok, \
            distance=distance)
    
    @staticmethod
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
        _projections.f90wrap_cone(p=self._handle, reg=reg._handle, ok=ok, \
            distance=distance)
    
    _dt_array_initialisers = []
    

geometrical_regions = Geometrical_Regions()

class Io_Ramses(f90wrap.runtime.FortranModule):
    """
    Module io_ramses
    
    
    Defined at read_amr_module.fpp lines 24-1291
    
    """
    @f90wrap.runtime.register_class("projections.hydroID")
    class hydroID(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=hydroid)
        
        
        Defined at read_amr_module.fpp lines 27-32
        
        """
        def __init__(self, handle=None):
            """
            self = Hydroid()
            
            
            Defined at read_amr_module.fpp lines 27-32
            
            
            Returns
            -------
            this : Hydroid
            	Object to be constructed
            
            
            Automatically generated constructor for hydroid
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _projections.f90wrap_hydroid_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Hydroid
            
            
            Defined at read_amr_module.fpp lines 27-32
            
            Parameters
            ----------
            this : Hydroid
            	Object to be destructed
            
            
            Automatically generated destructor for hydroid
            """
            if self._alloc:
                _projections.f90wrap_hydroid_finalise(this=self._handle)
        
        @property
        def nvar(self):
            """
            Element nvar ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 28
            
            """
            return _projections.f90wrap_hydroid__get__nvar(self._handle)
        
        @nvar.setter
        def nvar(self, nvar):
            _projections.f90wrap_hydroid__set__nvar(self._handle, nvar)
        
        @property
        def density(self):
            """
            Element density ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 29
            
            """
            return _projections.f90wrap_hydroid__get__density(self._handle)
        
        @density.setter
        def density(self, density):
            _projections.f90wrap_hydroid__set__density(self._handle, density)
        
        @property
        def vx(self):
            """
            Element vx ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 29
            
            """
            return _projections.f90wrap_hydroid__get__vx(self._handle)
        
        @vx.setter
        def vx(self, vx):
            _projections.f90wrap_hydroid__set__vx(self._handle, vx)
        
        @property
        def vy(self):
            """
            Element vy ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 29
            
            """
            return _projections.f90wrap_hydroid__get__vy(self._handle)
        
        @vy.setter
        def vy(self, vy):
            _projections.f90wrap_hydroid__set__vy(self._handle, vy)
        
        @property
        def vz(self):
            """
            Element vz ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 29
            
            """
            return _projections.f90wrap_hydroid__get__vz(self._handle)
        
        @vz.setter
        def vz(self, vz):
            _projections.f90wrap_hydroid__set__vz(self._handle, vz)
        
        @property
        def thermal_pressure(self):
            """
            Element thermal_pressure ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 29
            
            """
            return _projections.f90wrap_hydroid__get__thermal_pressure(self._handle)
        
        @thermal_pressure.setter
        def thermal_pressure(self, thermal_pressure):
            _projections.f90wrap_hydroid__set__thermal_pressure(self._handle, \
                thermal_pressure)
        
        @property
        def metallicity(self):
            """
            Element metallicity ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 29
            
            """
            return _projections.f90wrap_hydroid__get__metallicity(self._handle)
        
        @metallicity.setter
        def metallicity(self, metallicity):
            _projections.f90wrap_hydroid__set__metallicity(self._handle, metallicity)
        
        @property
        def blx(self):
            """
            Element blx ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 30
            
            """
            return _projections.f90wrap_hydroid__get__blx(self._handle)
        
        @blx.setter
        def blx(self, blx):
            _projections.f90wrap_hydroid__set__blx(self._handle, blx)
        
        @property
        def bly(self):
            """
            Element bly ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 30
            
            """
            return _projections.f90wrap_hydroid__get__bly(self._handle)
        
        @bly.setter
        def bly(self, bly):
            _projections.f90wrap_hydroid__set__bly(self._handle, bly)
        
        @property
        def blz(self):
            """
            Element blz ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 30
            
            """
            return _projections.f90wrap_hydroid__get__blz(self._handle)
        
        @blz.setter
        def blz(self, blz):
            _projections.f90wrap_hydroid__set__blz(self._handle, blz)
        
        @property
        def brx(self):
            """
            Element brx ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 30
            
            """
            return _projections.f90wrap_hydroid__get__brx(self._handle)
        
        @brx.setter
        def brx(self, brx):
            _projections.f90wrap_hydroid__set__brx(self._handle, brx)
        
        @property
        def bry(self):
            """
            Element bry ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 30
            
            """
            return _projections.f90wrap_hydroid__get__bry(self._handle)
        
        @bry.setter
        def bry(self, bry):
            _projections.f90wrap_hydroid__set__bry(self._handle, bry)
        
        @property
        def brz(self):
            """
            Element brz ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 30
            
            """
            return _projections.f90wrap_hydroid__get__brz(self._handle)
        
        @brz.setter
        def brz(self, brz):
            _projections.f90wrap_hydroid__set__brz(self._handle, brz)
        
        @property
        def ecr(self):
            """
            Element ecr ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 31
            
            """
            return _projections.f90wrap_hydroid__get__ecr(self._handle)
        
        @ecr.setter
        def ecr(self, ecr):
            _projections.f90wrap_hydroid__set__ecr(self._handle, ecr)
        
        @property
        def xhii(self):
            """
            Element xhii ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 32
            
            """
            return _projections.f90wrap_hydroid__get__xhii(self._handle)
        
        @xhii.setter
        def xhii(self, xhii):
            _projections.f90wrap_hydroid__set__xhii(self._handle, xhii)
        
        @property
        def xheii(self):
            """
            Element xheii ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 32
            
            """
            return _projections.f90wrap_hydroid__get__xheii(self._handle)
        
        @xheii.setter
        def xheii(self, xheii):
            _projections.f90wrap_hydroid__set__xheii(self._handle, xheii)
        
        @property
        def xheiii(self):
            """
            Element xheiii ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 32
            
            """
            return _projections.f90wrap_hydroid__get__xheiii(self._handle)
        
        @xheiii.setter
        def xheiii(self, xheiii):
            _projections.f90wrap_hydroid__set__xheiii(self._handle, xheiii)
        
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
            ret.append(',\n    ecr : ')
            ret.append(repr(self.ecr))
            ret.append(',\n    xhii : ')
            ret.append(repr(self.xhii))
            ret.append(',\n    xheii : ')
            ret.append(repr(self.xheii))
            ret.append(',\n    xheiii : ')
            ret.append(repr(self.xheiii))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("projections.amr_info")
    class amr_info(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=amr_info)
        
        
        Defined at read_amr_module.fpp lines 34-42
        
        """
        def __init__(self, handle=None):
            """
            self = Amr_Info()
            
            
            Defined at read_amr_module.fpp lines 34-42
            
            
            Returns
            -------
            this : Amr_Info
            	Object to be constructed
            
            
            Automatically generated constructor for amr_info
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _projections.f90wrap_amr_info_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Amr_Info
            
            
            Defined at read_amr_module.fpp lines 34-42
            
            Parameters
            ----------
            this : Amr_Info
            	Object to be destructed
            
            
            Automatically generated destructor for amr_info
            """
            if self._alloc:
                _projections.f90wrap_amr_info_finalise(this=self._handle)
        
        @property
        def ncpu(self):
            """
            Element ncpu ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 35
            
            """
            return _projections.f90wrap_amr_info__get__ncpu(self._handle)
        
        @ncpu.setter
        def ncpu(self, ncpu):
            _projections.f90wrap_amr_info__set__ncpu(self._handle, ncpu)
        
        @property
        def ndim(self):
            """
            Element ndim ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 35
            
            """
            return _projections.f90wrap_amr_info__get__ndim(self._handle)
        
        @ndim.setter
        def ndim(self, ndim):
            _projections.f90wrap_amr_info__set__ndim(self._handle, ndim)
        
        @property
        def nlevelmax(self):
            """
            Element nlevelmax ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 35
            
            """
            return _projections.f90wrap_amr_info__get__nlevelmax(self._handle)
        
        @nlevelmax.setter
        def nlevelmax(self, nlevelmax):
            _projections.f90wrap_amr_info__set__nlevelmax(self._handle, nlevelmax)
        
        @property
        def nboundary(self):
            """
            Element nboundary ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 35
            
            """
            return _projections.f90wrap_amr_info__get__nboundary(self._handle)
        
        @nboundary.setter
        def nboundary(self, nboundary):
            _projections.f90wrap_amr_info__set__nboundary(self._handle, nboundary)
        
        @property
        def twotondim(self):
            """
            Element twotondim ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 35
            
            """
            return _projections.f90wrap_amr_info__get__twotondim(self._handle)
        
        @twotondim.setter
        def twotondim(self, twotondim):
            _projections.f90wrap_amr_info__set__twotondim(self._handle, twotondim)
        
        @property
        def ndom(self):
            """
            Element ndom ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 35
            
            """
            return _projections.f90wrap_amr_info__get__ndom(self._handle)
        
        @ndom.setter
        def ndom(self, ndom):
            _projections.f90wrap_amr_info__set__ndom(self._handle, ndom)
        
        @property
        def levelmin(self):
            """
            Element levelmin ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 36
            
            """
            return _projections.f90wrap_amr_info__get__levelmin(self._handle)
        
        @levelmin.setter
        def levelmin(self, levelmin):
            _projections.f90wrap_amr_info__set__levelmin(self._handle, levelmin)
        
        @property
        def levelmax(self):
            """
            Element levelmax ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 36
            
            """
            return _projections.f90wrap_amr_info__get__levelmax(self._handle)
        
        @levelmax.setter
        def levelmax(self, levelmax):
            _projections.f90wrap_amr_info__set__levelmax(self._handle, levelmax)
        
        @property
        def lmax(self):
            """
            Element lmax ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 36
            
            """
            return _projections.f90wrap_amr_info__get__lmax(self._handle)
        
        @lmax.setter
        def lmax(self, lmax):
            _projections.f90wrap_amr_info__set__lmax(self._handle, lmax)
        
        @property
        def lmin(self):
            """
            Element lmin ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 36
            
            """
            return _projections.f90wrap_amr_info__get__lmin(self._handle)
        
        @lmin.setter
        def lmin(self, lmin):
            _projections.f90wrap_amr_info__set__lmin(self._handle, lmin)
        
        @property
        def ncpu_read(self):
            """
            Element ncpu_read ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 37
            
            """
            return _projections.f90wrap_amr_info__get__ncpu_read(self._handle)
        
        @ncpu_read.setter
        def ncpu_read(self, ncpu_read):
            _projections.f90wrap_amr_info__set__ncpu_read(self._handle, ncpu_read)
        
        @property
        def ordering(self):
            """
            Element ordering ftype=character(80) pytype=str
            
            
            Defined at read_amr_module.fpp line 38
            
            """
            return _projections.f90wrap_amr_info__get__ordering(self._handle)
        
        @ordering.setter
        def ordering(self, ordering):
            _projections.f90wrap_amr_info__set__ordering(self._handle, ordering)
        
        @property
        def cpu_list(self):
            """
            Element cpu_list ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 39
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_amr_info__array__cpu_list(self._handle)
            if array_handle in self._arrays:
                cpu_list = self._arrays[array_handle]
            else:
                cpu_list = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_amr_info__array__cpu_list)
                self._arrays[array_handle] = cpu_list
            return cpu_list
        
        @cpu_list.setter
        def cpu_list(self, cpu_list):
            self.cpu_list[...] = cpu_list
        
        @property
        def bound_key(self):
            """
            Element bound_key ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 40
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_amr_info__array__bound_key(self._handle)
            if array_handle in self._arrays:
                bound_key = self._arrays[array_handle]
            else:
                bound_key = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_amr_info__array__bound_key)
                self._arrays[array_handle] = bound_key
            return bound_key
        
        @bound_key.setter
        def bound_key(self, bound_key):
            self.bound_key[...] = bound_key
        
        @property
        def cpu_read(self):
            """
            Element cpu_read ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 41
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_amr_info__array__cpu_read(self._handle)
            if array_handle in self._arrays:
                cpu_read = self._arrays[array_handle]
            else:
                cpu_read = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_amr_info__array__cpu_read)
                self._arrays[array_handle] = cpu_read
            return cpu_read
        
        @cpu_read.setter
        def cpu_read(self, cpu_read):
            self.cpu_read[...] = cpu_read
        
        @property
        def xbound(self):
            """
            Element xbound ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 42
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_amr_info__array__xbound(self._handle)
            if array_handle in self._arrays:
                xbound = self._arrays[array_handle]
            else:
                xbound = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_amr_info__array__xbound)
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
            ret.append(',\n    levelmin : ')
            ret.append(repr(self.levelmin))
            ret.append(',\n    levelmax : ')
            ret.append(repr(self.levelmax))
            ret.append(',\n    lmax : ')
            ret.append(repr(self.lmax))
            ret.append(',\n    lmin : ')
            ret.append(repr(self.lmin))
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
        
    
    @f90wrap.runtime.register_class("projections.sim_info")
    class sim_info(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=sim_info)
        
        
        Defined at read_amr_module.fpp lines 44-49
        
        """
        def __init__(self, handle=None):
            """
            self = Sim_Info()
            
            
            Defined at read_amr_module.fpp lines 44-49
            
            
            Returns
            -------
            this : Sim_Info
            	Object to be constructed
            
            
            Automatically generated constructor for sim_info
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _projections.f90wrap_sim_info_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Sim_Info
            
            
            Defined at read_amr_module.fpp lines 44-49
            
            Parameters
            ----------
            this : Sim_Info
            	Object to be destructed
            
            
            Automatically generated destructor for sim_info
            """
            if self._alloc:
                _projections.f90wrap_sim_info_finalise(this=self._handle)
        
        @property
        def cosmo(self):
            """
            Element cosmo ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 45
            
            """
            return _projections.f90wrap_sim_info__get__cosmo(self._handle)
        
        @cosmo.setter
        def cosmo(self, cosmo):
            _projections.f90wrap_sim_info__set__cosmo(self._handle, cosmo)
        
        @property
        def family(self):
            """
            Element family ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 45
            
            """
            return _projections.f90wrap_sim_info__get__family(self._handle)
        
        @family.setter
        def family(self, family):
            _projections.f90wrap_sim_info__set__family(self._handle, family)
        
        @property
        def h0(self):
            """
            Element h0 ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _projections.f90wrap_sim_info__get__h0(self._handle)
        
        @h0.setter
        def h0(self, h0):
            _projections.f90wrap_sim_info__set__h0(self._handle, h0)
        
        @property
        def t(self):
            """
            Element t ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _projections.f90wrap_sim_info__get__t(self._handle)
        
        @t.setter
        def t(self, t):
            _projections.f90wrap_sim_info__set__t(self._handle, t)
        
        @property
        def aexp(self):
            """
            Element aexp ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _projections.f90wrap_sim_info__get__aexp(self._handle)
        
        @aexp.setter
        def aexp(self, aexp):
            _projections.f90wrap_sim_info__set__aexp(self._handle, aexp)
        
        @property
        def unit_l(self):
            """
            Element unit_l ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _projections.f90wrap_sim_info__get__unit_l(self._handle)
        
        @unit_l.setter
        def unit_l(self, unit_l):
            _projections.f90wrap_sim_info__set__unit_l(self._handle, unit_l)
        
        @property
        def unit_d(self):
            """
            Element unit_d ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _projections.f90wrap_sim_info__get__unit_d(self._handle)
        
        @unit_d.setter
        def unit_d(self, unit_d):
            _projections.f90wrap_sim_info__set__unit_d(self._handle, unit_d)
        
        @property
        def unit_t(self):
            """
            Element unit_t ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _projections.f90wrap_sim_info__get__unit_t(self._handle)
        
        @unit_t.setter
        def unit_t(self, unit_t):
            _projections.f90wrap_sim_info__set__unit_t(self._handle, unit_t)
        
        @property
        def unit_m(self):
            """
            Element unit_m ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _projections.f90wrap_sim_info__get__unit_m(self._handle)
        
        @unit_m.setter
        def unit_m(self, unit_m):
            _projections.f90wrap_sim_info__set__unit_m(self._handle, unit_m)
        
        @property
        def boxlen(self):
            """
            Element boxlen ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _projections.f90wrap_sim_info__get__boxlen(self._handle)
        
        @boxlen.setter
        def boxlen(self, boxlen):
            _projections.f90wrap_sim_info__set__boxlen(self._handle, boxlen)
        
        @property
        def omega_m(self):
            """
            Element omega_m ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _projections.f90wrap_sim_info__get__omega_m(self._handle)
        
        @omega_m.setter
        def omega_m(self, omega_m):
            _projections.f90wrap_sim_info__set__omega_m(self._handle, omega_m)
        
        @property
        def omega_l(self):
            """
            Element omega_l ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _projections.f90wrap_sim_info__get__omega_l(self._handle)
        
        @omega_l.setter
        def omega_l(self, omega_l):
            _projections.f90wrap_sim_info__set__omega_l(self._handle, omega_l)
        
        @property
        def omega_k(self):
            """
            Element omega_k ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _projections.f90wrap_sim_info__get__omega_k(self._handle)
        
        @omega_k.setter
        def omega_k(self, omega_k):
            _projections.f90wrap_sim_info__set__omega_k(self._handle, omega_k)
        
        @property
        def omega_b(self):
            """
            Element omega_b ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _projections.f90wrap_sim_info__get__omega_b(self._handle)
        
        @omega_b.setter
        def omega_b(self, omega_b):
            _projections.f90wrap_sim_info__set__omega_b(self._handle, omega_b)
        
        @property
        def time_tot(self):
            """
            Element time_tot ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 47
            
            """
            return _projections.f90wrap_sim_info__get__time_tot(self._handle)
        
        @time_tot.setter
        def time_tot(self, time_tot):
            _projections.f90wrap_sim_info__set__time_tot(self._handle, time_tot)
        
        @property
        def time_simu(self):
            """
            Element time_simu ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 47
            
            """
            return _projections.f90wrap_sim_info__get__time_simu(self._handle)
        
        @time_simu.setter
        def time_simu(self, time_simu):
            _projections.f90wrap_sim_info__set__time_simu(self._handle, time_simu)
        
        @property
        def n_frw(self):
            """
            Element n_frw ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 48
            
            """
            return _projections.f90wrap_sim_info__get__n_frw(self._handle)
        
        @n_frw.setter
        def n_frw(self, n_frw):
            _projections.f90wrap_sim_info__set__n_frw(self._handle, n_frw)
        
        @property
        def aexp_frw(self):
            """
            Element aexp_frw ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 49
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_sim_info__array__aexp_frw(self._handle)
            if array_handle in self._arrays:
                aexp_frw = self._arrays[array_handle]
            else:
                aexp_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_sim_info__array__aexp_frw)
                self._arrays[array_handle] = aexp_frw
            return aexp_frw
        
        @aexp_frw.setter
        def aexp_frw(self, aexp_frw):
            self.aexp_frw[...] = aexp_frw
        
        @property
        def hexp_frw(self):
            """
            Element hexp_frw ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 49
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_sim_info__array__hexp_frw(self._handle)
            if array_handle in self._arrays:
                hexp_frw = self._arrays[array_handle]
            else:
                hexp_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_sim_info__array__hexp_frw)
                self._arrays[array_handle] = hexp_frw
            return hexp_frw
        
        @hexp_frw.setter
        def hexp_frw(self, hexp_frw):
            self.hexp_frw[...] = hexp_frw
        
        @property
        def tau_frw(self):
            """
            Element tau_frw ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 49
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_sim_info__array__tau_frw(self._handle)
            if array_handle in self._arrays:
                tau_frw = self._arrays[array_handle]
            else:
                tau_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_sim_info__array__tau_frw)
                self._arrays[array_handle] = tau_frw
            return tau_frw
        
        @tau_frw.setter
        def tau_frw(self, tau_frw):
            self.tau_frw[...] = tau_frw
        
        @property
        def t_frw(self):
            """
            Element t_frw ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 49
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_sim_info__array__t_frw(self._handle)
            if array_handle in self._arrays:
                t_frw = self._arrays[array_handle]
            else:
                t_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_sim_info__array__t_frw)
                self._arrays[array_handle] = t_frw
            return t_frw
        
        @t_frw.setter
        def t_frw(self, t_frw):
            self.t_frw[...] = t_frw
        
        def __str__(self):
            ret = ['<sim_info>{\n']
            ret.append('    cosmo : ')
            ret.append(repr(self.cosmo))
            ret.append(',\n    family : ')
            ret.append(repr(self.family))
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
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("projections.level")
    class level(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=level)
        
        
        Defined at read_amr_module.fpp lines 51-62
        
        """
        def __init__(self, handle=None):
            """
            self = Level()
            
            
            Defined at read_amr_module.fpp lines 51-62
            
            
            Returns
            -------
            this : Level
            	Object to be constructed
            
            
            Automatically generated constructor for level
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _projections.f90wrap_level_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Level
            
            
            Defined at read_amr_module.fpp lines 51-62
            
            Parameters
            ----------
            this : Level
            	Object to be destructed
            
            
            Automatically generated destructor for level
            """
            if self._alloc:
                _projections.f90wrap_level_finalise(this=self._handle)
        
        @property
        def ilevel(self):
            """
            Element ilevel ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 52
            
            """
            return _projections.f90wrap_level__get__ilevel(self._handle)
        
        @ilevel.setter
        def ilevel(self, ilevel):
            _projections.f90wrap_level__set__ilevel(self._handle, ilevel)
        
        @property
        def ngrid(self):
            """
            Element ngrid ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _projections.f90wrap_level__get__ngrid(self._handle)
        
        @ngrid.setter
        def ngrid(self, ngrid):
            _projections.f90wrap_level__set__ngrid(self._handle, ngrid)
        
        @property
        def cube(self):
            """
            Element cube ftype=real(sgl) pytype=float
            
            
            Defined at read_amr_module.fpp line 54
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_level__array__cube(self._handle)
            if array_handle in self._arrays:
                cube = self._arrays[array_handle]
            else:
                cube = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_level__array__cube)
                self._arrays[array_handle] = cube
            return cube
        
        @cube.setter
        def cube(self, cube):
            self.cube[...] = cube
        
        @property
        def map(self):
            """
            Element map ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 55
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_level__array__map(self._handle)
            if array_handle in self._arrays:
                map = self._arrays[array_handle]
            else:
                map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_level__array__map)
                self._arrays[array_handle] = map
            return map
        
        @map.setter
        def map(self, map):
            self.map[...] = map
        
        @property
        def rho(self):
            """
            Element rho ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 56
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_level__array__rho(self._handle)
            if array_handle in self._arrays:
                rho = self._arrays[array_handle]
            else:
                rho = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_level__array__rho)
                self._arrays[array_handle] = rho
            return rho
        
        @rho.setter
        def rho(self, rho):
            self.rho[...] = rho
        
        @property
        def imin(self):
            """
            Element imin ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 57
            
            """
            return _projections.f90wrap_level__get__imin(self._handle)
        
        @imin.setter
        def imin(self, imin):
            _projections.f90wrap_level__set__imin(self._handle, imin)
        
        @property
        def imax(self):
            """
            Element imax ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 58
            
            """
            return _projections.f90wrap_level__get__imax(self._handle)
        
        @imax.setter
        def imax(self, imax):
            _projections.f90wrap_level__set__imax(self._handle, imax)
        
        @property
        def jmin(self):
            """
            Element jmin ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 59
            
            """
            return _projections.f90wrap_level__get__jmin(self._handle)
        
        @jmin.setter
        def jmin(self, jmin):
            _projections.f90wrap_level__set__jmin(self._handle, jmin)
        
        @property
        def jmax(self):
            """
            Element jmax ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 60
            
            """
            return _projections.f90wrap_level__get__jmax(self._handle)
        
        @jmax.setter
        def jmax(self, jmax):
            _projections.f90wrap_level__set__jmax(self._handle, jmax)
        
        @property
        def kmin(self):
            """
            Element kmin ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 61
            
            """
            return _projections.f90wrap_level__get__kmin(self._handle)
        
        @kmin.setter
        def kmin(self, kmin):
            _projections.f90wrap_level__set__kmin(self._handle, kmin)
        
        @property
        def kmax(self):
            """
            Element kmax ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 62
            
            """
            return _projections.f90wrap_level__get__kmax(self._handle)
        
        @kmax.setter
        def kmax(self, kmax):
            _projections.f90wrap_level__set__kmax(self._handle, kmax)
        
        def __str__(self):
            ret = ['<level>{\n']
            ret.append('    ilevel : ')
            ret.append(repr(self.ilevel))
            ret.append(',\n    ngrid : ')
            ret.append(repr(self.ngrid))
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
        
    
    @f90wrap.runtime.register_class("projections.data_handler")
    class data_handler(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=data_handler)
        
        
        Defined at read_amr_module.fpp lines 64-68
        
        """
        def __init__(self, handle=None):
            """
            self = Data_Handler()
            
            
            Defined at read_amr_module.fpp lines 64-68
            
            
            Returns
            -------
            this : Data_Handler
            	Object to be constructed
            
            
            Automatically generated constructor for data_handler
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _projections.f90wrap_data_handler_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Data_Handler
            
            
            Defined at read_amr_module.fpp lines 64-68
            
            Parameters
            ----------
            this : Data_Handler
            	Object to be destructed
            
            
            Automatically generated destructor for data_handler
            """
            if self._alloc:
                _projections.f90wrap_data_handler_finalise(this=self._handle)
        
        @property
        def name(self):
            """
            Element name ftype=character(80) pytype=str
            
            
            Defined at read_amr_module.fpp line 65
            
            """
            return _projections.f90wrap_data_handler__get__name(self._handle)
        
        @name.setter
        def name(self, name):
            _projections.f90wrap_data_handler__set__name(self._handle, name)
        
        @property
        def x_data(self):
            """
            Element x_data ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 66
            
            """
            return _projections.f90wrap_data_handler__get__x_data(self._handle)
        
        @x_data.setter
        def x_data(self, x_data):
            _projections.f90wrap_data_handler__set__x_data(self._handle, x_data)
        
        @property
        def y_data(self):
            """
            Element y_data ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 66
            
            """
            return _projections.f90wrap_data_handler__get__y_data(self._handle)
        
        @y_data.setter
        def y_data(self, y_data):
            _projections.f90wrap_data_handler__set__y_data(self._handle, y_data)
        
        @property
        def z_data(self):
            """
            Element z_data ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 66
            
            """
            return _projections.f90wrap_data_handler__get__z_data(self._handle)
        
        @z_data.setter
        def z_data(self, z_data):
            _projections.f90wrap_data_handler__set__z_data(self._handle, z_data)
        
        @property
        def nx(self):
            """
            Element nx ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 67
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_data_handler__array__nx(self._handle)
            if array_handle in self._arrays:
                nx = self._arrays[array_handle]
            else:
                nx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_data_handler__array__nx)
                self._arrays[array_handle] = nx
            return nx
        
        @nx.setter
        def nx(self, nx):
            self.nx[...] = nx
        
        @property
        def ny(self):
            """
            Element ny ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 67
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_data_handler__array__ny(self._handle)
            if array_handle in self._arrays:
                ny = self._arrays[array_handle]
            else:
                ny = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_data_handler__array__ny)
                self._arrays[array_handle] = ny
            return ny
        
        @ny.setter
        def ny(self, ny):
            self.ny[...] = ny
        
        @property
        def nz(self):
            """
            Element nz ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 67
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_data_handler__array__nz(self._handle)
            if array_handle in self._arrays:
                nz = self._arrays[array_handle]
            else:
                nz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_data_handler__array__nz)
                self._arrays[array_handle] = nz
            return nz
        
        @nz.setter
        def nz(self, nz):
            self.nz[...] = nz
        
        @property
        def x(self):
            """
            Element x ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 68
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_data_handler__array__x(self._handle)
            if array_handle in self._arrays:
                x = self._arrays[array_handle]
            else:
                x = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_data_handler__array__x)
                self._arrays[array_handle] = x
            return x
        
        @x.setter
        def x(self, x):
            self.x[...] = x
        
        @property
        def y(self):
            """
            Element y ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 68
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_data_handler__array__y(self._handle)
            if array_handle in self._arrays:
                y = self._arrays[array_handle]
            else:
                y = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_data_handler__array__y)
                self._arrays[array_handle] = y
            return y
        
        @y.setter
        def y(self, y):
            self.y[...] = y
        
        @property
        def z(self):
            """
            Element z ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 68
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_data_handler__array__z(self._handle)
            if array_handle in self._arrays:
                z = self._arrays[array_handle]
            else:
                z = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_data_handler__array__z)
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
        
    
    @f90wrap.runtime.register_class("projections.particle")
    class particle(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=particle)
        
        
        Defined at read_amr_module.fpp lines 70-73
        
        """
        def __init__(self, handle=None):
            """
            self = Particle()
            
            
            Defined at read_amr_module.fpp lines 70-73
            
            
            Returns
            -------
            this : Particle
            	Object to be constructed
            
            
            Automatically generated constructor for particle
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _projections.f90wrap_particle_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Particle
            
            
            Defined at read_amr_module.fpp lines 70-73
            
            Parameters
            ----------
            this : Particle
            	Object to be destructed
            
            
            Automatically generated destructor for particle
            """
            if self._alloc:
                _projections.f90wrap_particle_finalise(this=self._handle)
        
        @property
        def id(self):
            """
            Element id ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 71
            
            """
            return _projections.f90wrap_particle__get__id(self._handle)
        
        @id.setter
        def id(self, id):
            _projections.f90wrap_particle__set__id(self._handle, id)
        
        @property
        def x(self):
            """
            Element x ftype=type(vector) pytype=Vector
            
            
            Defined at read_amr_module.fpp line 72
            
            """
            x_handle = _projections.f90wrap_particle__get__x(self._handle)
            if tuple(x_handle) in self._objs:
                x = self._objs[tuple(x_handle)]
            else:
                x = vectors.vector.from_handle(x_handle)
                self._objs[tuple(x_handle)] = x
            return x
        
        @x.setter
        def x(self, x):
            x = x._handle
            _projections.f90wrap_particle__set__x(self._handle, x)
        
        @property
        def v(self):
            """
            Element v ftype=type(vector) pytype=Vector
            
            
            Defined at read_amr_module.fpp line 72
            
            """
            v_handle = _projections.f90wrap_particle__get__v(self._handle)
            if tuple(v_handle) in self._objs:
                v = self._objs[tuple(v_handle)]
            else:
                v = vectors.vector.from_handle(v_handle)
                self._objs[tuple(v_handle)] = v
            return v
        
        @v.setter
        def v(self, v):
            v = v._handle
            _projections.f90wrap_particle__set__v(self._handle, v)
        
        @property
        def m(self):
            """
            Element m ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 73
            
            """
            return _projections.f90wrap_particle__get__m(self._handle)
        
        @m.setter
        def m(self, m):
            _projections.f90wrap_particle__set__m(self._handle, m)
        
        @property
        def met(self):
            """
            Element met ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 73
            
            """
            return _projections.f90wrap_particle__get__met(self._handle)
        
        @met.setter
        def met(self, met):
            _projections.f90wrap_particle__set__met(self._handle, met)
        
        @property
        def imass(self):
            """
            Element imass ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 73
            
            """
            return _projections.f90wrap_particle__get__imass(self._handle)
        
        @imass.setter
        def imass(self, imass):
            _projections.f90wrap_particle__set__imass(self._handle, imass)
        
        @property
        def age(self):
            """
            Element age ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 73
            
            """
            return _projections.f90wrap_particle__get__age(self._handle)
        
        @age.setter
        def age(self, age):
            _projections.f90wrap_particle__set__age(self._handle, age)
        
        @property
        def tform(self):
            """
            Element tform ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 73
            
            """
            return _projections.f90wrap_particle__get__tform(self._handle)
        
        @tform.setter
        def tform(self, tform):
            _projections.f90wrap_particle__set__tform(self._handle, tform)
        
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
        
    
    @staticmethod
    def title(n, nchar):
        """
        title(n, nchar)
        
        
        Defined at read_amr_module.fpp lines 81-105
        
        Parameters
        ----------
        n : int
        nchar : str
        
        """
        _projections.f90wrap_title(n=n, nchar=nchar)
    
    @staticmethod
    def hilbert3d(x, y, z, order, bit_length, npoint):
        """
        hilbert3d(x, y, z, order, bit_length, npoint)
        
        
        Defined at read_amr_module.fpp lines 112-184
        
        Parameters
        ----------
        x : int array
        y : int array
        z : int array
        order : float array
        bit_length : int
        npoint : int
        
        """
        _projections.f90wrap_hilbert3d(x=x, y=y, z=z, order=order, \
            bit_length=bit_length, npoint=npoint)
    
    @staticmethod
    def check_lmax(ngridfile, amr):
        """
        check_lmax(ngridfile, amr)
        
        
        Defined at read_amr_module.fpp lines 192-205
        
        Parameters
        ----------
        ngridfile : int array
        amr : Amr_Info
        
        """
        _projections.f90wrap_check_lmax(ngridfile=ngridfile, amr=amr._handle)
    
    @staticmethod
    def check_families(repository, sim):
        """
        check_families(repository, sim)
        
        
        Defined at read_amr_module.fpp lines 214-240
        
        Parameters
        ----------
        repository : str
        sim : Sim_Info
        
        """
        _projections.f90wrap_check_families(repository=repository, sim=sim._handle)
    
    @staticmethod
    def read_hydrofile_descriptor(repository, varids):
        """
        read_hydrofile_descriptor(repository, varids)
        
        
        Defined at read_amr_module.fpp lines 249-284
        
        Parameters
        ----------
        repository : str
        varids : Hydroid
        
        """
        _projections.f90wrap_read_hydrofile_descriptor(repository=repository, \
            varids=varids._handle)
    
    @staticmethod
    def read_hydrofile_descriptor_old(repository, varids):
        """
        read_hydrofile_descriptor_old(repository, varids)
        
        
        Defined at read_amr_module.fpp lines 294-351
        
        Parameters
        ----------
        repository : str
        varids : Hydroid
        
        """
        _projections.f90wrap_read_hydrofile_descriptor_old(repository=repository, \
            varids=varids._handle)
    
    @staticmethod
    def select_from_descriptor_ids(self, newvar, newid):
        """
        select_from_descriptor_ids(self, newvar, newid)
        
        
        Defined at read_amr_module.fpp lines 353-430
        
        Parameters
        ----------
        varids : Hydroid
        newvar : str
        newid : int
        
        """
        _projections.f90wrap_select_from_descriptor_ids(varids=self._handle, \
            newvar=newvar, newid=newid)
    
    @staticmethod
    def read_hydrofile_descriptor_new(repository, varids):
        """
        read_hydrofile_descriptor_new(repository, varids)
        
        
        Defined at read_amr_module.fpp lines 440-472
        
        Parameters
        ----------
        repository : str
        varids : Hydroid
        
        """
        _projections.f90wrap_read_hydrofile_descriptor_new(repository=repository, \
            varids=varids._handle)
    
    @staticmethod
    def getvarvalue(self, reg, dx, x, var, varname, value):
        """
        getvarvalue(self, reg, dx, x, var, varname, value)
        
        
        Defined at read_amr_module.fpp lines 481-704
        
        Parameters
        ----------
        varids : Hydroid
        reg : Region
        dx : float
        x : Vector
        var : float array
        varname : str
        value : float
        
        """
        _projections.f90wrap_getvarvalue(varids=self._handle, reg=reg._handle, dx=dx, \
            x=x._handle, var=var, varname=varname, value=value)
    
    @staticmethod
    def init_amr_read(repository, amr, sim):
        """
        init_amr_read(repository, amr, sim)
        
        
        Defined at read_amr_module.fpp lines 713-789
        
        Parameters
        ----------
        repository : str
        amr : Amr_Info
        sim : Sim_Info
        
        """
        _projections.f90wrap_init_amr_read(repository=repository, amr=amr._handle, \
            sim=sim._handle)
    
    @staticmethod
    def get_cpu_map(self, amr):
        """
        get_cpu_map(self, amr)
        
        
        Defined at read_amr_module.fpp lines 797-886
        
        Parameters
        ----------
        reg : Region
        amr : Amr_Info
        
        """
        _projections.f90wrap_get_cpu_map(reg=self._handle, amr=amr._handle)
    
    @staticmethod
    def getparttype(self, ptype):
        """
        getparttype(self, ptype)
        
        
        Defined at read_amr_module.fpp lines 888-896
        
        Parameters
        ----------
        part : Particle
        ptype : str
        
        """
        _projections.f90wrap_getparttype(part=self._handle, ptype=ptype)
    
    @staticmethod
    def getpartvalue(self, reg, part, var, value, dx=None):
        """
        getpartvalue(self, reg, part, var, value[, dx])
        
        
        Defined at read_amr_module.fpp lines 898-1290
        
        Parameters
        ----------
        sim : Sim_Info
        reg : Region
        part : Particle
        var : str
        value : float
        dx : Vector
        
        """
        _projections.f90wrap_getpartvalue(sim=self._handle, reg=reg._handle, \
            part=part._handle, var=var, value=value, dx=None if dx is None else \
            dx._handle)
    
    _dt_array_initialisers = []
    

io_ramses = Io_Ramses()

class Filtering(f90wrap.runtime.FortranModule):
    """
    Module filtering
    
    
    Defined at read_amr_module.fpp lines 1293-1387
    
    """
    @f90wrap.runtime.register_class("projections.filter")
    class filter(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=filter)
        
        
        Defined at read_amr_module.fpp lines 1296-1301
        
        """
        def __init__(self, handle=None):
            """
            self = Filter()
            
            
            Defined at read_amr_module.fpp lines 1296-1301
            
            
            Returns
            -------
            this : Filter
            	Object to be constructed
            
            
            Automatically generated constructor for filter
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _projections.f90wrap_filter_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Filter
            
            
            Defined at read_amr_module.fpp lines 1296-1301
            
            Parameters
            ----------
            this : Filter
            	Object to be destructed
            
            
            Automatically generated destructor for filter
            """
            if self._alloc:
                _projections.f90wrap_filter_finalise(this=self._handle)
        
        @property
        def name(self):
            """
            Element name ftype=character(128) pytype=str
            
            
            Defined at read_amr_module.fpp line 1297
            
            """
            return _projections.f90wrap_filter__get__name(self._handle)
        
        @name.setter
        def name(self, name):
            _projections.f90wrap_filter__set__name(self._handle, name)
        
        @property
        def ncond(self):
            """
            Element ncond ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 1298
            
            """
            return _projections.f90wrap_filter__get__ncond(self._handle)
        
        @ncond.setter
        def ncond(self, ncond):
            _projections.f90wrap_filter__set__ncond(self._handle, ncond)
        
        @property
        def cond_vars(self):
            """
            Element cond_vars ftype=character(128) pytype=str
            
            
            Defined at read_amr_module.fpp line 1299
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_filter__array__cond_vars(self._handle)
            if array_handle in self._arrays:
                cond_vars = self._arrays[array_handle]
            else:
                cond_vars = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_filter__array__cond_vars)
                self._arrays[array_handle] = cond_vars
            return cond_vars
        
        @cond_vars.setter
        def cond_vars(self, cond_vars):
            self.cond_vars[...] = cond_vars
        
        @property
        def cond_ops(self):
            """
            Element cond_ops ftype=character(2) pytype=str
            
            
            Defined at read_amr_module.fpp line 1300
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_filter__array__cond_ops(self._handle)
            if array_handle in self._arrays:
                cond_ops = self._arrays[array_handle]
            else:
                cond_ops = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_filter__array__cond_ops)
                self._arrays[array_handle] = cond_ops
            return cond_ops
        
        @cond_ops.setter
        def cond_ops(self, cond_ops):
            self.cond_ops[...] = cond_ops
        
        @property
        def cond_vals(self):
            """
            Element cond_vals ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 1301
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_filter__array__cond_vals(self._handle)
            if array_handle in self._arrays:
                cond_vals = self._arrays[array_handle]
            else:
                cond_vals = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_filter__array__cond_vals)
                self._arrays[array_handle] = cond_vals
            return cond_vals
        
        @cond_vals.setter
        def cond_vals(self, cond_vals):
            self.cond_vals[...] = cond_vals
        
        def __str__(self):
            ret = ['<filter>{\n']
            ret.append('    name : ')
            ret.append(repr(self.name))
            ret.append(',\n    ncond : ')
            ret.append(repr(self.ncond))
            ret.append(',\n    cond_vars : ')
            ret.append(repr(self.cond_vars))
            ret.append(',\n    cond_ops : ')
            ret.append(repr(self.cond_ops))
            ret.append(',\n    cond_vals : ')
            ret.append(repr(self.cond_vals))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def allocate_filter(self):
        """
        allocate_filter(self)
        
        
        Defined at read_amr_module.fpp lines 1304-1309
        
        Parameters
        ----------
        filt : Filter
        
        """
        _projections.f90wrap_allocate_filter(filt=self._handle)
    
    @staticmethod
    def cond_string_to_filter(str, filt):
        """
        cond_string_to_filter(str, filt)
        
        
        Defined at read_amr_module.fpp lines 1311-1315
        
        Parameters
        ----------
        str : str
        filt : Filter
        
        """
        _projections.f90wrap_cond_string_to_filter(str=str, filt=filt._handle)
    
    @staticmethod
    def filter_cell(self, reg, filt, cell_x, cell_dx, cell_var):
        """
        filter_cell = filter_cell(self, reg, filt, cell_x, cell_dx, cell_var)
        
        
        Defined at read_amr_module.fpp lines 1317-1350
        
        Parameters
        ----------
        varids : Hydroid
        reg : Region
        filt : Filter
        cell_x : Vector
        cell_dx : float
        cell_var : float array
        
        Returns
        -------
        filter_cell : bool
        
        """
        filter_cell = _projections.f90wrap_filter_cell(varids=self._handle, \
            reg=reg._handle, filt=filt._handle, cell_x=cell_x._handle, cell_dx=cell_dx, \
            cell_var=cell_var)
        return filter_cell
    
    @staticmethod
    def filter_particle(self, reg, filt, part, dx=None):
        """
        filter_particle = filter_particle(self, reg, filt, part[, dx])
        
        
        Defined at read_amr_module.fpp lines 1352-1387
        
        Parameters
        ----------
        sim : Sim_Info
        reg : Region
        filt : Filter
        part : Particle
        dx : Vector
        
        Returns
        -------
        filter_particle : bool
        
        """
        filter_particle = _projections.f90wrap_filter_particle(sim=self._handle, \
            reg=reg._handle, filt=filt._handle, part=part._handle, dx=None if dx is None \
            else dx._handle)
        return filter_particle
    
    _dt_array_initialisers = []
    

filtering = Filtering()

class Obs_Instruments(f90wrap.runtime.FortranModule):
    """
    Module obs_instruments
    
    
    Defined at ramses2map.fpp lines 5-192
    
    """
    @f90wrap.runtime.register_class("projections.camera")
    class camera(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=camera)
        
        
        Defined at ramses2map.fpp lines 8-12
        
        """
        def __init__(self, handle=None):
            """
            self = Camera()
            
            
            Defined at ramses2map.fpp lines 8-12
            
            
            Returns
            -------
            this : Camera
            	Object to be constructed
            
            
            Automatically generated constructor for camera
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _projections.f90wrap_camera_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Camera
            
            
            Defined at ramses2map.fpp lines 8-12
            
            Parameters
            ----------
            this : Camera
            	Object to be destructed
            
            
            Automatically generated destructor for camera
            """
            if self._alloc:
                _projections.f90wrap_camera_finalise(this=self._handle)
        
        @property
        def centre(self):
            """
            Element centre ftype=type(vector) pytype=Vector
            
            
            Defined at ramses2map.fpp line 9
            
            """
            centre_handle = _projections.f90wrap_camera__get__centre(self._handle)
            if tuple(centre_handle) in self._objs:
                centre = self._objs[tuple(centre_handle)]
            else:
                centre = vectors.vector.from_handle(centre_handle)
                self._objs[tuple(centre_handle)] = centre
            return centre
        
        @centre.setter
        def centre(self, centre):
            centre = centre._handle
            _projections.f90wrap_camera__set__centre(self._handle, centre)
        
        @property
        def los_axis(self):
            """
            Element los_axis ftype=type(vector) pytype=Vector
            
            
            Defined at ramses2map.fpp line 9
            
            """
            los_axis_handle = _projections.f90wrap_camera__get__los_axis(self._handle)
            if tuple(los_axis_handle) in self._objs:
                los_axis = self._objs[tuple(los_axis_handle)]
            else:
                los_axis = vectors.vector.from_handle(los_axis_handle)
                self._objs[tuple(los_axis_handle)] = los_axis
            return los_axis
        
        @los_axis.setter
        def los_axis(self, los_axis):
            los_axis = los_axis._handle
            _projections.f90wrap_camera__set__los_axis(self._handle, los_axis)
        
        @property
        def up_vector(self):
            """
            Element up_vector ftype=type(vector) pytype=Vector
            
            
            Defined at ramses2map.fpp line 9
            
            """
            up_vector_handle = _projections.f90wrap_camera__get__up_vector(self._handle)
            if tuple(up_vector_handle) in self._objs:
                up_vector = self._objs[tuple(up_vector_handle)]
            else:
                up_vector = vectors.vector.from_handle(up_vector_handle)
                self._objs[tuple(up_vector_handle)] = up_vector
            return up_vector
        
        @up_vector.setter
        def up_vector(self, up_vector):
            up_vector = up_vector._handle
            _projections.f90wrap_camera__set__up_vector(self._handle, up_vector)
        
        @property
        def region_size(self):
            """
            Element region_size ftype=real(dbl) pytype=float
            
            
            Defined at ramses2map.fpp line 10
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_camera__array__region_size(self._handle)
            if array_handle in self._arrays:
                region_size = self._arrays[array_handle]
            else:
                region_size = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_camera__array__region_size)
                self._arrays[array_handle] = region_size
            return region_size
        
        @region_size.setter
        def region_size(self, region_size):
            self.region_size[...] = region_size
        
        @property
        def distance(self):
            """
            Element distance ftype=real(dbl) pytype=float
            
            
            Defined at ramses2map.fpp line 11
            
            """
            return _projections.f90wrap_camera__get__distance(self._handle)
        
        @distance.setter
        def distance(self, distance):
            _projections.f90wrap_camera__set__distance(self._handle, distance)
        
        @property
        def far_cut_depth(self):
            """
            Element far_cut_depth ftype=real(dbl) pytype=float
            
            
            Defined at ramses2map.fpp line 11
            
            """
            return _projections.f90wrap_camera__get__far_cut_depth(self._handle)
        
        @far_cut_depth.setter
        def far_cut_depth(self, far_cut_depth):
            _projections.f90wrap_camera__set__far_cut_depth(self._handle, far_cut_depth)
        
        @property
        def map_max_size(self):
            """
            Element map_max_size ftype=integer  pytype=int
            
            
            Defined at ramses2map.fpp line 12
            
            """
            return _projections.f90wrap_camera__get__map_max_size(self._handle)
        
        @map_max_size.setter
        def map_max_size(self, map_max_size):
            _projections.f90wrap_camera__set__map_max_size(self._handle, map_max_size)
        
        def __str__(self):
            ret = ['<camera>{\n']
            ret.append('    centre : ')
            ret.append(repr(self.centre))
            ret.append(',\n    los_axis : ')
            ret.append(repr(self.los_axis))
            ret.append(',\n    up_vector : ')
            ret.append(repr(self.up_vector))
            ret.append(',\n    region_size : ')
            ret.append(repr(self.region_size))
            ret.append(',\n    distance : ')
            ret.append(repr(self.distance))
            ret.append(',\n    far_cut_depth : ')
            ret.append(repr(self.far_cut_depth))
            ret.append(',\n    map_max_size : ')
            ret.append(repr(self.map_max_size))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def log2(x):
        """
        log2 = log2(x)
        
        
        Defined at ramses2map.fpp lines 16-19
        
        Parameters
        ----------
        x : float
        
        Returns
        -------
        log2 : float
        
        """
        log2 = _projections.f90wrap_log2(x=x)
        return log2
    
    @staticmethod
    def init_camera(self, los_axis, up_vector, region_size, distance, far_cut_depth, \
        map_max_size):
        """
        init_camera = init_camera(self, los_axis, up_vector, region_size, distance, \
            far_cut_depth, map_max_size)
        
        
        Defined at ramses2map.fpp lines 21-33
        
        Parameters
        ----------
        centre : Vector
        los_axis : Vector
        up_vector : Vector
        region_size : float array
        distance : float
        far_cut_depth : float
        map_max_size : int
        
        Returns
        -------
        init_camera : Camera
        
        """
        init_camera = _projections.f90wrap_init_camera(centre=self._handle, \
            los_axis=los_axis._handle, up_vector=up_vector._handle, \
            region_size=region_size, distance=distance, far_cut_depth=far_cut_depth, \
            map_max_size=map_max_size)
        init_camera = \
            f90wrap.runtime.lookup_class("projections.camera").from_handle(init_camera, \
            alloc=True)
        return init_camera
    
    @staticmethod
    def get_required_resolution(self):
        """
        get_required_resolution = get_required_resolution(self)
        
        
        Defined at ramses2map.fpp lines 35-38
        
        Parameters
        ----------
        cam : Camera
        
        Returns
        -------
        get_required_resolution : int
        
        """
        get_required_resolution = \
            _projections.f90wrap_get_required_resolution(cam=self._handle)
        return get_required_resolution
    
    @staticmethod
    def get_map_size(self, n_map):
        """
        get_map_size(self, n_map)
        
        
        Defined at ramses2map.fpp lines 40-52
        
        Parameters
        ----------
        cam : Camera
        n_map : int array
        
        """
        _projections.f90wrap_get_map_size(cam=self._handle, n_map=n_map)
    
    @staticmethod
    def get_map_box(self, box):
        """
        get_map_box(self, box)
        
        
        Defined at ramses2map.fpp lines 54-64
        
        Parameters
        ----------
        cam : Camera
        box : Region
        
        """
        _projections.f90wrap_get_map_box(cam=self._handle, box=box._handle)
    
    @staticmethod
    def get_camera_basis(self, cam_basis):
        """
        get_camera_basis(self, cam_basis)
        
        
        Defined at ramses2map.fpp lines 86-94
        
        Parameters
        ----------
        cam : Camera
        cam_basis : Basis
        
        """
        _projections.f90wrap_get_camera_basis(cam=self._handle, \
            cam_basis=cam_basis._handle)
    
    @staticmethod
    def los_transformation(self, trans_matrix):
        """
        los_transformation(self, trans_matrix)
        
        
        Defined at ramses2map.fpp lines 96-107
        
        Parameters
        ----------
        cam : Camera
        trans_matrix : float array
        
        """
        _projections.f90wrap_los_transformation(cam=self._handle, \
            trans_matrix=trans_matrix)
    
    @staticmethod
    def get_bounding_box(self, bbox):
        """
        get_bounding_box(self, bbox)
        
        
        Defined at ramses2map.fpp lines 109-153
        
        Parameters
        ----------
        cam : Camera
        bbox : Region
        
        """
        _projections.f90wrap_get_bounding_box(cam=self._handle, bbox=bbox._handle)
    
    @staticmethod
    def deproject_points(self, npoints, points):
        """
        deproject_points(self, npoints, points)
        
        
        Defined at ramses2map.fpp lines 155-172
        
        Parameters
        ----------
        cam : Camera
        npoints : int
        points : float array
        
        """
        _projections.f90wrap_deproject_points(cam=self._handle, npoints=npoints, \
            points=points)
    
    @staticmethod
    def project_points(self, npoints, points):
        """
        project_points(self, npoints, points)
        
        
        Defined at ramses2map.fpp lines 174-191
        
        Parameters
        ----------
        cam : Camera
        npoints : int
        points : float array
        
        """
        _projections.f90wrap_project_points(cam=self._handle, npoints=npoints, \
            points=points)
    
    _dt_array_initialisers = []
    

obs_instruments = Obs_Instruments()

class Maps(f90wrap.runtime.FortranModule):
    """
    Module maps
    
    
    Defined at ramses2map.fpp lines 194-1016
    
    """
    @f90wrap.runtime.register_class("projections.projection_handler")
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
            result = _projections.f90wrap_projection_handler_initialise()
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
                _projections.f90wrap_projection_handler_finalise(this=self._handle)
        
        @property
        def pov(self):
            """
            Element pov ftype=character(128) pytype=str
            
            
            Defined at ramses2map.fpp line 202
            
            """
            return _projections.f90wrap_projection_handler__get__pov(self._handle)
        
        @pov.setter
        def pov(self, pov):
            _projections.f90wrap_projection_handler__set__pov(self._handle, pov)
        
        @property
        def nvars(self):
            """
            Element nvars ftype=integer  pytype=int
            
            
            Defined at ramses2map.fpp line 203
            
            """
            return _projections.f90wrap_projection_handler__get__nvars(self._handle)
        
        @nvars.setter
        def nvars(self, nvars):
            _projections.f90wrap_projection_handler__set__nvars(self._handle, nvars)
        
        @property
        def varnames(self):
            """
            Element varnames ftype=character(128) pytype=str
            
            
            Defined at ramses2map.fpp line 204
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_projection_handler__array__varnames(self._handle)
            if array_handle in self._arrays:
                varnames = self._arrays[array_handle]
            else:
                varnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_projection_handler__array__varnames)
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
            return _projections.f90wrap_projection_handler__get__weightvar(self._handle)
        
        @weightvar.setter
        def weightvar(self, weightvar):
            _projections.f90wrap_projection_handler__set__weightvar(self._handle, weightvar)
        
        @property
        def toto(self):
            """
            Element toto ftype=real(dbl) pytype=float
            
            
            Defined at ramses2map.fpp line 206
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _projections.f90wrap_projection_handler__array__toto(self._handle)
            if array_handle in self._arrays:
                toto = self._arrays[array_handle]
            else:
                toto = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _projections.f90wrap_projection_handler__array__toto)
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
        
    
    @staticmethod
    def allocate_projection_handler(self):
        """
        allocate_projection_handler(self)
        
        
        Defined at ramses2map.fpp lines 209-212
        
        Parameters
        ----------
        proj : Projection_Handler
        
        """
        _projections.f90wrap_allocate_projection_handler(proj=self._handle)
    
    @staticmethod
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
        _projections.f90wrap_projection_hydro(repository=repository, cam=cam._handle, \
            bulk_velocity=bulk_velocity._handle, proj=proj._handle)
    
    @staticmethod
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
        _projections.f90wrap_project_cells(repository=repository, amr=amr._handle, \
            bbox=bbox._handle, varids=varids._handle, cam=cam._handle, \
            proj=proj._handle)
    
    @staticmethod
    def projection_parts(repository, cam, bulk_velocity, proj):
        """
        projection_parts(repository, cam, bulk_velocity, proj)
        
        
        Defined at ramses2map.fpp lines 541-562
        
        Parameters
        ----------
        repository : str
        cam : Camera
        bulk_velocity : Vector
        proj : Projection_Handler
        
        """
        _projections.f90wrap_projection_parts(repository=repository, cam=cam._handle, \
            bulk_velocity=bulk_velocity._handle, proj=proj._handle)
    
    @staticmethod
    def project_particles(repository, amr, sim, bbox, cam, proj):
        """
        project_particles(repository, amr, sim, bbox, cam, proj)
        
        
        Defined at ramses2map.fpp lines 564-745
        
        Parameters
        ----------
        repository : str
        amr : Amr_Info
        sim : Sim_Info
        bbox : Region
        cam : Camera
        proj : Projection_Handler
        
        """
        _projections.f90wrap_project_particles(repository=repository, amr=amr._handle, \
            sim=sim._handle, bbox=bbox._handle, cam=cam._handle, proj=proj._handle)
    
    @staticmethod
    def healpix_hydro(repository, reg, nside, proj):
        """
        healpix_hydro(repository, reg, nside, proj)
        
        
        Defined at ramses2map.fpp lines 747-762
        
        Parameters
        ----------
        repository : str
        reg : Region
        nside : int
        proj : Projection_Handler
        
        """
        _projections.f90wrap_healpix_hydro(repository=repository, reg=reg._handle, \
            nside=nside, proj=proj._handle)
    
    @staticmethod
    def project_cells_hpix(repository, amr, reg, varids, nside, proj):
        """
        project_cells_hpix(repository, amr, reg, varids, nside, proj)
        
        
        Defined at ramses2map.fpp lines 764-1016
        
        Parameters
        ----------
        repository : str
        amr : Amr_Info
        reg : Region
        varids : Hydroid
        nside : int
        proj : Projection_Handler
        
        """
        _projections.f90wrap_project_cells_hpix(repository=repository, amr=amr._handle, \
            reg=reg._handle, varids=varids._handle, nside=nside, proj=proj._handle)
    
    _dt_array_initialisers = []
    

maps = Maps()

