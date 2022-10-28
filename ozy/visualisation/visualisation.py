from __future__ import print_function, absolute_import, division
import _visualisation
import f90wrap.runtime
import logging

class Vectors(f90wrap.runtime.FortranModule):
    """
    Module vectors
    
    
    Defined at linalg_module.fpp lines 23-159
    
    """
    @f90wrap.runtime.register_class("visualisation.vector")
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
            result = _visualisation.f90wrap_vector_initialise()
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
                _visualisation.f90wrap_vector_finalise(this=self._handle)
        
        @property
        def x(self):
            """
            Element x ftype=real(dbl) pytype=float
            
            
            Defined at linalg_module.fpp line 27
            
            """
            return _visualisation.f90wrap_vector__get__x(self._handle)
        
        @x.setter
        def x(self, x):
            _visualisation.f90wrap_vector__set__x(self._handle, x)
        
        @property
        def y(self):
            """
            Element y ftype=real(dbl) pytype=float
            
            
            Defined at linalg_module.fpp line 28
            
            """
            return _visualisation.f90wrap_vector__get__y(self._handle)
        
        @y.setter
        def y(self, y):
            _visualisation.f90wrap_vector__set__y(self._handle, y)
        
        @property
        def z(self):
            """
            Element z ftype=real(dbl) pytype=float
            
            
            Defined at linalg_module.fpp line 29
            
            """
            return _visualisation.f90wrap_vector__get__z(self._handle)
        
        @z.setter
        def z(self, z):
            _visualisation.f90wrap_vector__set__z(self._handle, z)
        
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
        
    
    @f90wrap.runtime.register_class("visualisation.array_vectors")
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
            result = _visualisation.f90wrap_array_vectors_initialise()
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
                _visualisation.f90wrap_array_vectors_finalise(this=self._handle)
        
        @property
        def n(self):
            """
            Element n ftype=integer  pytype=int
            
            
            Defined at linalg_module.fpp line 32
            
            """
            return _visualisation.f90wrap_array_vectors__get__n(self._handle)
        
        @n.setter
        def n(self, n):
            _visualisation.f90wrap_array_vectors__set__n(self._handle, n)
        
        def init_array_list(self):
            self.list = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _visualisation.f90wrap_array_vectors__array_getitem__list,
                                            _visualisation.f90wrap_array_vectors__array_setitem__list,
                                            _visualisation.f90wrap_array_vectors__array_len__list,
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
        magnitude = _visualisation.f90wrap_magnitude(vec_1=self._handle)
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
        _visualisation.f90wrap_array_to_vector(vec_result=self._handle, array=array)
    
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
        _visualisation.f90wrap_vector_to_array(array_result=array_result, \
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
        _visualisation.f90wrap_euler_matrix(r=r, dim=dim, angle=angle, ap=ap)
    
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
        _visualisation.f90wrap_rotate_vector_single(vec=self._handle, \
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
        _visualisation.f90wrap_rotate_vector_array(vec=self._handle, \
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
    @f90wrap.runtime.register_class("visualisation.basis")
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
            result = _visualisation.f90wrap_basis_initialise()
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
                _visualisation.f90wrap_basis_finalise(this=self._handle)
        
        def init_array_u(self):
            self.u = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _visualisation.f90wrap_basis__array_getitem__u,
                                            _visualisation.f90wrap_basis__array_setitem__u,
                                            _visualisation.f90wrap_basis__array_len__u,
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
        _visualisation.f90wrap_initialise_basis(this=self._handle)
    
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
        _visualisation.f90wrap_mgramschmidt(vecs=self._handle, e=e._handle)
    
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
        r_sphere = _visualisation.f90wrap_r_sphere(p=self._handle)
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
        theta_sphere = _visualisation.f90wrap_theta_sphere(p=self._handle)
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
        phi_sphere = _visualisation.f90wrap_phi_sphere(p=self._handle)
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
        r_cyl = _visualisation.f90wrap_r_cyl(p=self._handle)
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
        phi_cyl = _visualisation.f90wrap_phi_cyl(p=self._handle)
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
        _visualisation.f90wrap_spherical_basis_from_cartesian(p=self._handle, \
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
        _visualisation.f90wrap_cylindrical_basis_from_cartesian(p=self._handle, \
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
        _visualisation.f90wrap_new_z_coordinates(axis=self._handle, \
            transformation_matrix=transformation_matrix, errormsg=errormsg)
    
    _dt_array_initialisers = []
    

coordinate_systems = Coordinate_Systems()

class Geometrical_Regions(f90wrap.runtime.FortranModule):
    """
    Module geometrical_regions
    
    
    Defined at coordinates_module.fpp lines 191-420
    
    """
    @f90wrap.runtime.register_class("visualisation.region")
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
            result = _visualisation.f90wrap_region_initialise()
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
                _visualisation.f90wrap_region_finalise(this=self._handle)
        
        @property
        def name(self):
            """
            Element name ftype=character(128) pytype=str
            
            
            Defined at coordinates_module.fpp line 196
            
            """
            return _visualisation.f90wrap_region__get__name(self._handle)
        
        @name.setter
        def name(self, name):
            _visualisation.f90wrap_region__set__name(self._handle, name)
        
        @property
        def centre(self):
            """
            Element centre ftype=type(vector) pytype=Vector
            
            
            Defined at coordinates_module.fpp line 197
            
            """
            centre_handle = _visualisation.f90wrap_region__get__centre(self._handle)
            if tuple(centre_handle) in self._objs:
                centre = self._objs[tuple(centre_handle)]
            else:
                centre = vectors.vector.from_handle(centre_handle)
                self._objs[tuple(centre_handle)] = centre
            return centre
        
        @centre.setter
        def centre(self, centre):
            centre = centre._handle
            _visualisation.f90wrap_region__set__centre(self._handle, centre)
        
        @property
        def axis(self):
            """
            Element axis ftype=type(vector) pytype=Vector
            
            
            Defined at coordinates_module.fpp line 197
            
            """
            axis_handle = _visualisation.f90wrap_region__get__axis(self._handle)
            if tuple(axis_handle) in self._objs:
                axis = self._objs[tuple(axis_handle)]
            else:
                axis = vectors.vector.from_handle(axis_handle)
                self._objs[tuple(axis_handle)] = axis
            return axis
        
        @axis.setter
        def axis(self, axis):
            axis = axis._handle
            _visualisation.f90wrap_region__set__axis(self._handle, axis)
        
        @property
        def bulk_velocity(self):
            """
            Element bulk_velocity ftype=type(vector) pytype=Vector
            
            
            Defined at coordinates_module.fpp line 197
            
            """
            bulk_velocity_handle = \
                _visualisation.f90wrap_region__get__bulk_velocity(self._handle)
            if tuple(bulk_velocity_handle) in self._objs:
                bulk_velocity = self._objs[tuple(bulk_velocity_handle)]
            else:
                bulk_velocity = vectors.vector.from_handle(bulk_velocity_handle)
                self._objs[tuple(bulk_velocity_handle)] = bulk_velocity
            return bulk_velocity
        
        @bulk_velocity.setter
        def bulk_velocity(self, bulk_velocity):
            bulk_velocity = bulk_velocity._handle
            _visualisation.f90wrap_region__set__bulk_velocity(self._handle, bulk_velocity)
        
        @property
        def xmin(self):
            """
            Element xmin ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _visualisation.f90wrap_region__get__xmin(self._handle)
        
        @xmin.setter
        def xmin(self, xmin):
            _visualisation.f90wrap_region__set__xmin(self._handle, xmin)
        
        @property
        def xmax(self):
            """
            Element xmax ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _visualisation.f90wrap_region__get__xmax(self._handle)
        
        @xmax.setter
        def xmax(self, xmax):
            _visualisation.f90wrap_region__set__xmax(self._handle, xmax)
        
        @property
        def ymin(self):
            """
            Element ymin ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _visualisation.f90wrap_region__get__ymin(self._handle)
        
        @ymin.setter
        def ymin(self, ymin):
            _visualisation.f90wrap_region__set__ymin(self._handle, ymin)
        
        @property
        def ymax(self):
            """
            Element ymax ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _visualisation.f90wrap_region__get__ymax(self._handle)
        
        @ymax.setter
        def ymax(self, ymax):
            _visualisation.f90wrap_region__set__ymax(self._handle, ymax)
        
        @property
        def zmin(self):
            """
            Element zmin ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _visualisation.f90wrap_region__get__zmin(self._handle)
        
        @zmin.setter
        def zmin(self, zmin):
            _visualisation.f90wrap_region__set__zmin(self._handle, zmin)
        
        @property
        def zmax(self):
            """
            Element zmax ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _visualisation.f90wrap_region__get__zmax(self._handle)
        
        @zmax.setter
        def zmax(self, zmax):
            _visualisation.f90wrap_region__set__zmax(self._handle, zmax)
        
        @property
        def rmin(self):
            """
            Element rmin ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 199
            
            """
            return _visualisation.f90wrap_region__get__rmin(self._handle)
        
        @rmin.setter
        def rmin(self, rmin):
            _visualisation.f90wrap_region__set__rmin(self._handle, rmin)
        
        @property
        def rmax(self):
            """
            Element rmax ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 199
            
            """
            return _visualisation.f90wrap_region__get__rmax(self._handle)
        
        @rmax.setter
        def rmax(self, rmax):
            _visualisation.f90wrap_region__set__rmax(self._handle, rmax)
        
        @property
        def angle(self):
            """
            Element angle ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 200
            
            """
            return _visualisation.f90wrap_region__get__angle(self._handle)
        
        @angle.setter
        def angle(self, angle):
            _visualisation.f90wrap_region__set__angle(self._handle, angle)
        
        @property
        def criteria_name(self):
            """
            Element criteria_name ftype=character(128) pytype=str
            
            
            Defined at coordinates_module.fpp line 201
            
            """
            return _visualisation.f90wrap_region__get__criteria_name(self._handle)
        
        @criteria_name.setter
        def criteria_name(self, criteria_name):
            _visualisation.f90wrap_region__set__criteria_name(self._handle, criteria_name)
        
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
        _visualisation.f90wrap_limits(reg=self._handle, lim=lim)
    
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
        _visualisation.f90wrap_checkifinside(pos=pos, reg=reg._handle, ok=ok, \
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
        _visualisation.f90wrap_cube(p=self._handle, reg=reg._handle, ok=ok, \
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
        _visualisation.f90wrap_sphere(p=self._handle, reg=reg._handle, ok=ok, \
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
        _visualisation.f90wrap_cylinder(p=self._handle, reg=reg._handle, ok=ok, \
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
        _visualisation.f90wrap_cone(p=self._handle, reg=reg._handle, ok=ok, \
            distance=distance)
    
    _dt_array_initialisers = []
    

geometrical_regions = Geometrical_Regions()

class Cooling_Module(f90wrap.runtime.FortranModule):
    """
    Module cooling_module
    
    
    Defined at cooling_module.fpp lines 24-163
    
    """
    @f90wrap.runtime.register_class("visualisation.cooling_table")
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
            result = _visualisation.f90wrap_cooling_table_initialise()
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
                _visualisation.f90wrap_cooling_table_finalise(this=self._handle)
        
        @property
        def n1(self):
            """
            Element n1 ftype=integer pytype=int
            
            
            Defined at cooling_module.fpp line 28
            
            """
            return _visualisation.f90wrap_cooling_table__get__n1(self._handle)
        
        @n1.setter
        def n1(self, n1):
            _visualisation.f90wrap_cooling_table__set__n1(self._handle, n1)
        
        @property
        def n2(self):
            """
            Element n2 ftype=integer pytype=int
            
            
            Defined at cooling_module.fpp line 29
            
            """
            return _visualisation.f90wrap_cooling_table__get__n2(self._handle)
        
        @n2.setter
        def n2(self, n2):
            _visualisation.f90wrap_cooling_table__set__n2(self._handle, n2)
        
        @property
        def nh(self):
            """
            Element nh ftype=real(kind=8) pytype=float
            
            
            Defined at cooling_module.fpp line 30
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _visualisation.f90wrap_cooling_table__array__nh(self._handle)
            if array_handle in self._arrays:
                nh = self._arrays[array_handle]
            else:
                nh = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__nh)
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
                _visualisation.f90wrap_cooling_table__array__t2(self._handle)
            if array_handle in self._arrays:
                t2 = self._arrays[array_handle]
            else:
                t2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__t2)
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
                _visualisation.f90wrap_cooling_table__array__cool(self._handle)
            if array_handle in self._arrays:
                cool = self._arrays[array_handle]
            else:
                cool = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__cool)
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
                _visualisation.f90wrap_cooling_table__array__heat(self._handle)
            if array_handle in self._arrays:
                heat = self._arrays[array_handle]
            else:
                heat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__heat)
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
                _visualisation.f90wrap_cooling_table__array__cool_com(self._handle)
            if array_handle in self._arrays:
                cool_com = self._arrays[array_handle]
            else:
                cool_com = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__cool_com)
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
                _visualisation.f90wrap_cooling_table__array__heat_com(self._handle)
            if array_handle in self._arrays:
                heat_com = self._arrays[array_handle]
            else:
                heat_com = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__heat_com)
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
                _visualisation.f90wrap_cooling_table__array__metal(self._handle)
            if array_handle in self._arrays:
                metal = self._arrays[array_handle]
            else:
                metal = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__metal)
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
                _visualisation.f90wrap_cooling_table__array__cool_prime(self._handle)
            if array_handle in self._arrays:
                cool_prime = self._arrays[array_handle]
            else:
                cool_prime = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__cool_prime)
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
                _visualisation.f90wrap_cooling_table__array__heat_prime(self._handle)
            if array_handle in self._arrays:
                heat_prime = self._arrays[array_handle]
            else:
                heat_prime = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__heat_prime)
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
                _visualisation.f90wrap_cooling_table__array__cool_com_prime(self._handle)
            if array_handle in self._arrays:
                cool_com_prime = self._arrays[array_handle]
            else:
                cool_com_prime = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__cool_com_prime)
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
                _visualisation.f90wrap_cooling_table__array__heat_com_prime(self._handle)
            if array_handle in self._arrays:
                heat_com_prime = self._arrays[array_handle]
            else:
                heat_com_prime = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__heat_com_prime)
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
                _visualisation.f90wrap_cooling_table__array__metal_prime(self._handle)
            if array_handle in self._arrays:
                metal_prime = self._arrays[array_handle]
            else:
                metal_prime = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__metal_prime)
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
                _visualisation.f90wrap_cooling_table__array__mu(self._handle)
            if array_handle in self._arrays:
                mu = self._arrays[array_handle]
            else:
                mu = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__mu)
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
                _visualisation.f90wrap_cooling_table__array__n_spec(self._handle)
            if array_handle in self._arrays:
                n_spec = self._arrays[array_handle]
            else:
                n_spec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_cooling_table__array__n_spec)
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
        
    
    @staticmethod
    def retrieve_table(repository, mytable):
        """
        retrieve_table(repository, mytable)
        
        
        Defined at cooling_module.fpp lines 58-63
        
        Parameters
        ----------
        repository : str
        mytable : Cooling_Table
        
        """
        _visualisation.f90wrap_retrieve_table(repository=repository, \
            mytable=mytable._handle)
    
    @staticmethod
    def read_cool(filename):
        """
        read_cool(filename)
        
        
        Defined at cooling_module.fpp lines 65-111
        
        Parameters
        ----------
        filename : str
        
        """
        _visualisation.f90wrap_read_cool(filename=filename)
    
    @staticmethod
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
        lambda_, lambda_prime = _visualisation.f90wrap_solve_cooling(nh=nh, t2=t2, \
            zsolar=zsolar)
        return lambda_, lambda_prime
    
    @property
    def if_species_abundances(self):
        """
        Element if_species_abundances ftype=logical pytype=bool
        
        
        Defined at cooling_module.fpp line 47
        
        """
        return _visualisation.f90wrap_cooling_module__get__if_species_abundances()
    
    @property
    def self_shielding(self):
        """
        Element self_shielding ftype=logical pytype=bool
        
        
        Defined at cooling_module.fpp line 48
        
        """
        return _visualisation.f90wrap_cooling_module__get__self_shielding()
    
    @property
    def nbin_t_fix(self):
        """
        Element nbin_t_fix ftype=integer pytype=int
        
        
        Defined at cooling_module.fpp line 50
        
        """
        return _visualisation.f90wrap_cooling_module__get__nbin_t_fix()
    
    @property
    def nbin_n_fix(self):
        """
        Element nbin_n_fix ftype=integer pytype=int
        
        
        Defined at cooling_module.fpp line 51
        
        """
        return _visualisation.f90wrap_cooling_module__get__nbin_n_fix()
    
    @property
    def nh_min_fix(self):
        """
        Element nh_min_fix ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 52
        
        """
        return _visualisation.f90wrap_cooling_module__get__nh_min_fix()
    
    @property
    def nh_max_fix(self):
        """
        Element nh_max_fix ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 53
        
        """
        return _visualisation.f90wrap_cooling_module__get__nh_max_fix()
    
    @property
    def t2_min_fix(self):
        """
        Element t2_min_fix ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 54
        
        """
        return _visualisation.f90wrap_cooling_module__get__t2_min_fix()
    
    @property
    def t2_max_fix(self):
        """
        Element t2_max_fix ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 55
        
        """
        return _visualisation.f90wrap_cooling_module__get__t2_max_fix()
    
    @property
    def logt2max(self):
        """
        Element logt2max ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 56
        
        """
        return _visualisation.f90wrap_cooling_module__get__logt2max()
    
    @logt2max.setter
    def logt2max(self, logt2max):
        _visualisation.f90wrap_cooling_module__set__logt2max(logt2max)
    
    @property
    def dlog_nh(self):
        """
        Element dlog_nh ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 56
        
        """
        return _visualisation.f90wrap_cooling_module__get__dlog_nh()
    
    @dlog_nh.setter
    def dlog_nh(self, dlog_nh):
        _visualisation.f90wrap_cooling_module__set__dlog_nh(dlog_nh)
    
    @property
    def dlog_t2(self):
        """
        Element dlog_t2 ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 56
        
        """
        return _visualisation.f90wrap_cooling_module__get__dlog_t2()
    
    @dlog_t2.setter
    def dlog_t2(self, dlog_t2):
        _visualisation.f90wrap_cooling_module__set__dlog_t2(dlog_t2)
    
    @property
    def h(self):
        """
        Element h ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 56
        
        """
        return _visualisation.f90wrap_cooling_module__get__h()
    
    @h.setter
    def h(self, h):
        _visualisation.f90wrap_cooling_module__set__h(h)
    
    @property
    def h2(self):
        """
        Element h2 ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 56
        
        """
        return _visualisation.f90wrap_cooling_module__get__h2()
    
    @h2.setter
    def h2(self, h2):
        _visualisation.f90wrap_cooling_module__set__h2(h2)
    
    @property
    def h3(self):
        """
        Element h3 ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 56
        
        """
        return _visualisation.f90wrap_cooling_module__get__h3()
    
    @h3.setter
    def h3(self, h3):
        _visualisation.f90wrap_cooling_module__set__h3(h3)
    
    @property
    def precoeff(self):
        """
        Element precoeff ftype=real(kind=8) pytype=float
        
        
        Defined at cooling_module.fpp line 56
        
        """
        return _visualisation.f90wrap_cooling_module__get__precoeff()
    
    @precoeff.setter
    def precoeff(self, precoeff):
        _visualisation.f90wrap_cooling_module__set__precoeff(precoeff)
    
    def __str__(self):
        ret = ['<cooling_module>{\n']
        ret.append('    if_species_abundances : ')
        ret.append(repr(self.if_species_abundances))
        ret.append(',\n    self_shielding : ')
        ret.append(repr(self.self_shielding))
        ret.append(',\n    nbin_t_fix : ')
        ret.append(repr(self.nbin_t_fix))
        ret.append(',\n    nbin_n_fix : ')
        ret.append(repr(self.nbin_n_fix))
        ret.append(',\n    nh_min_fix : ')
        ret.append(repr(self.nh_min_fix))
        ret.append(',\n    nh_max_fix : ')
        ret.append(repr(self.nh_max_fix))
        ret.append(',\n    t2_min_fix : ')
        ret.append(repr(self.t2_min_fix))
        ret.append(',\n    t2_max_fix : ')
        ret.append(repr(self.t2_max_fix))
        ret.append(',\n    logt2max : ')
        ret.append(repr(self.logt2max))
        ret.append(',\n    dlog_nh : ')
        ret.append(repr(self.dlog_nh))
        ret.append(',\n    dlog_t2 : ')
        ret.append(repr(self.dlog_t2))
        ret.append(',\n    h : ')
        ret.append(repr(self.h))
        ret.append(',\n    h2 : ')
        ret.append(repr(self.h2))
        ret.append(',\n    h3 : ')
        ret.append(repr(self.h3))
        ret.append(',\n    precoeff : ')
        ret.append(repr(self.precoeff))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

cooling_module = Cooling_Module()

class Io_Ramses(f90wrap.runtime.FortranModule):
    """
    Module io_ramses
    
    
    Defined at read_amr_module.fpp lines 24-2514
    
    """
    @f90wrap.runtime.register_class("visualisation.hydroID")
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
            result = _visualisation.f90wrap_hydroid_initialise()
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
                _visualisation.f90wrap_hydroid_finalise(this=self._handle)
        
        @property
        def nvar(self):
            """
            Element nvar ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 30
            
            """
            return _visualisation.f90wrap_hydroid__get__nvar(self._handle)
        
        @nvar.setter
        def nvar(self, nvar):
            _visualisation.f90wrap_hydroid__set__nvar(self._handle, nvar)
        
        @property
        def density(self):
            """
            Element density ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 31
            
            """
            return _visualisation.f90wrap_hydroid__get__density(self._handle)
        
        @density.setter
        def density(self, density):
            _visualisation.f90wrap_hydroid__set__density(self._handle, density)
        
        @property
        def vx(self):
            """
            Element vx ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 31
            
            """
            return _visualisation.f90wrap_hydroid__get__vx(self._handle)
        
        @vx.setter
        def vx(self, vx):
            _visualisation.f90wrap_hydroid__set__vx(self._handle, vx)
        
        @property
        def vy(self):
            """
            Element vy ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 31
            
            """
            return _visualisation.f90wrap_hydroid__get__vy(self._handle)
        
        @vy.setter
        def vy(self, vy):
            _visualisation.f90wrap_hydroid__set__vy(self._handle, vy)
        
        @property
        def vz(self):
            """
            Element vz ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 31
            
            """
            return _visualisation.f90wrap_hydroid__get__vz(self._handle)
        
        @vz.setter
        def vz(self, vz):
            _visualisation.f90wrap_hydroid__set__vz(self._handle, vz)
        
        @property
        def thermal_pressure(self):
            """
            Element thermal_pressure ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 31
            
            """
            return _visualisation.f90wrap_hydroid__get__thermal_pressure(self._handle)
        
        @thermal_pressure.setter
        def thermal_pressure(self, thermal_pressure):
            _visualisation.f90wrap_hydroid__set__thermal_pressure(self._handle, \
                thermal_pressure)
        
        @property
        def metallicity(self):
            """
            Element metallicity ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 31
            
            """
            return _visualisation.f90wrap_hydroid__get__metallicity(self._handle)
        
        @metallicity.setter
        def metallicity(self, metallicity):
            _visualisation.f90wrap_hydroid__set__metallicity(self._handle, metallicity)
        
        @property
        def blx(self):
            """
            Element blx ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 32
            
            """
            return _visualisation.f90wrap_hydroid__get__blx(self._handle)
        
        @blx.setter
        def blx(self, blx):
            _visualisation.f90wrap_hydroid__set__blx(self._handle, blx)
        
        @property
        def bly(self):
            """
            Element bly ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 32
            
            """
            return _visualisation.f90wrap_hydroid__get__bly(self._handle)
        
        @bly.setter
        def bly(self, bly):
            _visualisation.f90wrap_hydroid__set__bly(self._handle, bly)
        
        @property
        def blz(self):
            """
            Element blz ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 32
            
            """
            return _visualisation.f90wrap_hydroid__get__blz(self._handle)
        
        @blz.setter
        def blz(self, blz):
            _visualisation.f90wrap_hydroid__set__blz(self._handle, blz)
        
        @property
        def brx(self):
            """
            Element brx ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 32
            
            """
            return _visualisation.f90wrap_hydroid__get__brx(self._handle)
        
        @brx.setter
        def brx(self, brx):
            _visualisation.f90wrap_hydroid__set__brx(self._handle, brx)
        
        @property
        def bry(self):
            """
            Element bry ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 32
            
            """
            return _visualisation.f90wrap_hydroid__get__bry(self._handle)
        
        @bry.setter
        def bry(self, bry):
            _visualisation.f90wrap_hydroid__set__bry(self._handle, bry)
        
        @property
        def brz(self):
            """
            Element brz ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 32
            
            """
            return _visualisation.f90wrap_hydroid__get__brz(self._handle)
        
        @brz.setter
        def brz(self, brz):
            _visualisation.f90wrap_hydroid__set__brz(self._handle, brz)
        
        @property
        def cr_pressure(self):
            """
            Element cr_pressure ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 33
            
            """
            return _visualisation.f90wrap_hydroid__get__cr_pressure(self._handle)
        
        @cr_pressure.setter
        def cr_pressure(self, cr_pressure):
            _visualisation.f90wrap_hydroid__set__cr_pressure(self._handle, cr_pressure)
        
        @property
        def xhii(self):
            """
            Element xhii ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 34
            
            """
            return _visualisation.f90wrap_hydroid__get__xhii(self._handle)
        
        @xhii.setter
        def xhii(self, xhii):
            _visualisation.f90wrap_hydroid__set__xhii(self._handle, xhii)
        
        @property
        def xheii(self):
            """
            Element xheii ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 34
            
            """
            return _visualisation.f90wrap_hydroid__get__xheii(self._handle)
        
        @xheii.setter
        def xheii(self, xheii):
            _visualisation.f90wrap_hydroid__set__xheii(self._handle, xheii)
        
        @property
        def xheiii(self):
            """
            Element xheiii ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 34
            
            """
            return _visualisation.f90wrap_hydroid__get__xheiii(self._handle)
        
        @xheiii.setter
        def xheiii(self, xheiii):
            _visualisation.f90wrap_hydroid__set__xheiii(self._handle, xheiii)
        
        @property
        def dust_density(self):
            """
            Element dust_density ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 35
            
            """
            return _visualisation.f90wrap_hydroid__get__dust_density(self._handle)
        
        @dust_density.setter
        def dust_density(self, dust_density):
            _visualisation.f90wrap_hydroid__set__dust_density(self._handle, dust_density)
        
        @property
        def sigma2(self):
            """
            Element sigma2 ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 36
            
            """
            return _visualisation.f90wrap_hydroid__get__sigma2(self._handle)
        
        @sigma2.setter
        def sigma2(self, sigma2):
            _visualisation.f90wrap_hydroid__set__sigma2(self._handle, sigma2)
        
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
        
    
    @f90wrap.runtime.register_class("visualisation.amr_info")
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
            result = _visualisation.f90wrap_amr_info_initialise()
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
                _visualisation.f90wrap_amr_info_finalise(this=self._handle)
        
        @property
        def ncpu(self):
            """
            Element ncpu ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 39
            
            """
            return _visualisation.f90wrap_amr_info__get__ncpu(self._handle)
        
        @ncpu.setter
        def ncpu(self, ncpu):
            _visualisation.f90wrap_amr_info__set__ncpu(self._handle, ncpu)
        
        @property
        def ndim(self):
            """
            Element ndim ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 39
            
            """
            return _visualisation.f90wrap_amr_info__get__ndim(self._handle)
        
        @ndim.setter
        def ndim(self, ndim):
            _visualisation.f90wrap_amr_info__set__ndim(self._handle, ndim)
        
        @property
        def nlevelmax(self):
            """
            Element nlevelmax ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 39
            
            """
            return _visualisation.f90wrap_amr_info__get__nlevelmax(self._handle)
        
        @nlevelmax.setter
        def nlevelmax(self, nlevelmax):
            _visualisation.f90wrap_amr_info__set__nlevelmax(self._handle, nlevelmax)
        
        @property
        def nboundary(self):
            """
            Element nboundary ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 39
            
            """
            return _visualisation.f90wrap_amr_info__get__nboundary(self._handle)
        
        @nboundary.setter
        def nboundary(self, nboundary):
            _visualisation.f90wrap_amr_info__set__nboundary(self._handle, nboundary)
        
        @property
        def twotondim(self):
            """
            Element twotondim ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 39
            
            """
            return _visualisation.f90wrap_amr_info__get__twotondim(self._handle)
        
        @twotondim.setter
        def twotondim(self, twotondim):
            _visualisation.f90wrap_amr_info__set__twotondim(self._handle, twotondim)
        
        @property
        def ndom(self):
            """
            Element ndom ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 39
            
            """
            return _visualisation.f90wrap_amr_info__get__ndom(self._handle)
        
        @ndom.setter
        def ndom(self, ndom):
            _visualisation.f90wrap_amr_info__set__ndom(self._handle, ndom)
        
        @property
        def twondim(self):
            """
            Element twondim ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 40
            
            """
            return _visualisation.f90wrap_amr_info__get__twondim(self._handle)
        
        @twondim.setter
        def twondim(self, twondim):
            _visualisation.f90wrap_amr_info__set__twondim(self._handle, twondim)
        
        @property
        def ngridmax(self):
            """
            Element ngridmax ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 40
            
            """
            return _visualisation.f90wrap_amr_info__get__ngridmax(self._handle)
        
        @ngridmax.setter
        def ngridmax(self, ngridmax):
            _visualisation.f90wrap_amr_info__set__ngridmax(self._handle, ngridmax)
        
        @property
        def ncoarse(self):
            """
            Element ncoarse ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 40
            
            """
            return _visualisation.f90wrap_amr_info__get__ncoarse(self._handle)
        
        @ncoarse.setter
        def ncoarse(self, ncoarse):
            _visualisation.f90wrap_amr_info__set__ncoarse(self._handle, ncoarse)
        
        @property
        def levelmin(self):
            """
            Element levelmin ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 41
            
            """
            return _visualisation.f90wrap_amr_info__get__levelmin(self._handle)
        
        @levelmin.setter
        def levelmin(self, levelmin):
            _visualisation.f90wrap_amr_info__set__levelmin(self._handle, levelmin)
        
        @property
        def levelmax(self):
            """
            Element levelmax ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 41
            
            """
            return _visualisation.f90wrap_amr_info__get__levelmax(self._handle)
        
        @levelmax.setter
        def levelmax(self, levelmax):
            _visualisation.f90wrap_amr_info__set__levelmax(self._handle, levelmax)
        
        @property
        def lmax(self):
            """
            Element lmax ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 41
            
            """
            return _visualisation.f90wrap_amr_info__get__lmax(self._handle)
        
        @lmax.setter
        def lmax(self, lmax):
            _visualisation.f90wrap_amr_info__set__lmax(self._handle, lmax)
        
        @property
        def lmin(self):
            """
            Element lmin ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 41
            
            """
            return _visualisation.f90wrap_amr_info__get__lmin(self._handle)
        
        @lmin.setter
        def lmin(self, lmin):
            _visualisation.f90wrap_amr_info__set__lmin(self._handle, lmin)
        
        @property
        def active_lmax(self):
            """
            Element active_lmax ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 41
            
            """
            return _visualisation.f90wrap_amr_info__get__active_lmax(self._handle)
        
        @active_lmax.setter
        def active_lmax(self, active_lmax):
            _visualisation.f90wrap_amr_info__set__active_lmax(self._handle, active_lmax)
        
        @property
        def ncpu_read(self):
            """
            Element ncpu_read ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 42
            
            """
            return _visualisation.f90wrap_amr_info__get__ncpu_read(self._handle)
        
        @ncpu_read.setter
        def ncpu_read(self, ncpu_read):
            _visualisation.f90wrap_amr_info__set__ncpu_read(self._handle, ncpu_read)
        
        @property
        def ordering(self):
            """
            Element ordering ftype=character(80) pytype=str
            
            
            Defined at read_amr_module.fpp line 43
            
            """
            return _visualisation.f90wrap_amr_info__get__ordering(self._handle)
        
        @ordering.setter
        def ordering(self, ordering):
            _visualisation.f90wrap_amr_info__set__ordering(self._handle, ordering)
        
        @property
        def cpu_list(self):
            """
            Element cpu_list ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 44
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _visualisation.f90wrap_amr_info__array__cpu_list(self._handle)
            if array_handle in self._arrays:
                cpu_list = self._arrays[array_handle]
            else:
                cpu_list = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_amr_info__array__cpu_list)
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
                _visualisation.f90wrap_amr_info__array__bound_key(self._handle)
            if array_handle in self._arrays:
                bound_key = self._arrays[array_handle]
            else:
                bound_key = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_amr_info__array__bound_key)
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
                _visualisation.f90wrap_amr_info__array__cpu_read(self._handle)
            if array_handle in self._arrays:
                cpu_read = self._arrays[array_handle]
            else:
                cpu_read = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_amr_info__array__cpu_read)
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
                _visualisation.f90wrap_amr_info__array__xbound(self._handle)
            if array_handle in self._arrays:
                xbound = self._arrays[array_handle]
            else:
                xbound = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_amr_info__array__xbound)
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
        
    
    @f90wrap.runtime.register_class("visualisation.sim_info")
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
            result = _visualisation.f90wrap_sim_info_initialise()
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
                _visualisation.f90wrap_sim_info_finalise(this=self._handle)
        
        @property
        def cosmo(self):
            """
            Element cosmo ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 50
            
            """
            return _visualisation.f90wrap_sim_info__get__cosmo(self._handle)
        
        @cosmo.setter
        def cosmo(self, cosmo):
            _visualisation.f90wrap_sim_info__set__cosmo(self._handle, cosmo)
        
        @property
        def family(self):
            """
            Element family ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 50
            
            """
            return _visualisation.f90wrap_sim_info__get__family(self._handle)
        
        @family.setter
        def family(self, family):
            _visualisation.f90wrap_sim_info__set__family(self._handle, family)
        
        @property
        def dm(self):
            """
            Element dm ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 51
            
            """
            return _visualisation.f90wrap_sim_info__get__dm(self._handle)
        
        @dm.setter
        def dm(self, dm):
            _visualisation.f90wrap_sim_info__set__dm(self._handle, dm)
        
        @property
        def hydro(self):
            """
            Element hydro ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 51
            
            """
            return _visualisation.f90wrap_sim_info__get__hydro(self._handle)
        
        @hydro.setter
        def hydro(self, hydro):
            _visualisation.f90wrap_sim_info__set__hydro(self._handle, hydro)
        
        @property
        def mhd(self):
            """
            Element mhd ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 51
            
            """
            return _visualisation.f90wrap_sim_info__get__mhd(self._handle)
        
        @mhd.setter
        def mhd(self, mhd):
            _visualisation.f90wrap_sim_info__set__mhd(self._handle, mhd)
        
        @property
        def cr(self):
            """
            Element cr ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 51
            
            """
            return _visualisation.f90wrap_sim_info__get__cr(self._handle)
        
        @cr.setter
        def cr(self, cr):
            _visualisation.f90wrap_sim_info__set__cr(self._handle, cr)
        
        @property
        def rt(self):
            """
            Element rt ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 51
            
            """
            return _visualisation.f90wrap_sim_info__get__rt(self._handle)
        
        @rt.setter
        def rt(self, rt):
            _visualisation.f90wrap_sim_info__set__rt(self._handle, rt)
        
        @property
        def bh(self):
            """
            Element bh ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 51
            
            """
            return _visualisation.f90wrap_sim_info__get__bh(self._handle)
        
        @bh.setter
        def bh(self, bh):
            _visualisation.f90wrap_sim_info__set__bh(self._handle, bh)
        
        @property
        def cr_st(self):
            """
            Element cr_st ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 52
            
            """
            return _visualisation.f90wrap_sim_info__get__cr_st(self._handle)
        
        @cr_st.setter
        def cr_st(self, cr_st):
            _visualisation.f90wrap_sim_info__set__cr_st(self._handle, cr_st)
        
        @property
        def cr_heat(self):
            """
            Element cr_heat ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 52
            
            """
            return _visualisation.f90wrap_sim_info__get__cr_heat(self._handle)
        
        @cr_heat.setter
        def cr_heat(self, cr_heat):
            _visualisation.f90wrap_sim_info__set__cr_heat(self._handle, cr_heat)
        
        @property
        def dust(self):
            """
            Element dust ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 52
            
            """
            return _visualisation.f90wrap_sim_info__get__dust(self._handle)
        
        @dust.setter
        def dust(self, dust):
            _visualisation.f90wrap_sim_info__set__dust(self._handle, dust)
        
        @property
        def h0(self):
            """
            Element h0 ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _visualisation.f90wrap_sim_info__get__h0(self._handle)
        
        @h0.setter
        def h0(self, h0):
            _visualisation.f90wrap_sim_info__set__h0(self._handle, h0)
        
        @property
        def t(self):
            """
            Element t ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _visualisation.f90wrap_sim_info__get__t(self._handle)
        
        @t.setter
        def t(self, t):
            _visualisation.f90wrap_sim_info__set__t(self._handle, t)
        
        @property
        def aexp(self):
            """
            Element aexp ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _visualisation.f90wrap_sim_info__get__aexp(self._handle)
        
        @aexp.setter
        def aexp(self, aexp):
            _visualisation.f90wrap_sim_info__set__aexp(self._handle, aexp)
        
        @property
        def unit_l(self):
            """
            Element unit_l ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _visualisation.f90wrap_sim_info__get__unit_l(self._handle)
        
        @unit_l.setter
        def unit_l(self, unit_l):
            _visualisation.f90wrap_sim_info__set__unit_l(self._handle, unit_l)
        
        @property
        def unit_d(self):
            """
            Element unit_d ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _visualisation.f90wrap_sim_info__get__unit_d(self._handle)
        
        @unit_d.setter
        def unit_d(self, unit_d):
            _visualisation.f90wrap_sim_info__set__unit_d(self._handle, unit_d)
        
        @property
        def unit_t(self):
            """
            Element unit_t ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _visualisation.f90wrap_sim_info__get__unit_t(self._handle)
        
        @unit_t.setter
        def unit_t(self, unit_t):
            _visualisation.f90wrap_sim_info__set__unit_t(self._handle, unit_t)
        
        @property
        def unit_m(self):
            """
            Element unit_m ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _visualisation.f90wrap_sim_info__get__unit_m(self._handle)
        
        @unit_m.setter
        def unit_m(self, unit_m):
            _visualisation.f90wrap_sim_info__set__unit_m(self._handle, unit_m)
        
        @property
        def unit_v(self):
            """
            Element unit_v ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _visualisation.f90wrap_sim_info__get__unit_v(self._handle)
        
        @unit_v.setter
        def unit_v(self, unit_v):
            _visualisation.f90wrap_sim_info__set__unit_v(self._handle, unit_v)
        
        @property
        def boxlen(self):
            """
            Element boxlen ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _visualisation.f90wrap_sim_info__get__boxlen(self._handle)
        
        @boxlen.setter
        def boxlen(self, boxlen):
            _visualisation.f90wrap_sim_info__set__boxlen(self._handle, boxlen)
        
        @property
        def omega_m(self):
            """
            Element omega_m ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _visualisation.f90wrap_sim_info__get__omega_m(self._handle)
        
        @omega_m.setter
        def omega_m(self, omega_m):
            _visualisation.f90wrap_sim_info__set__omega_m(self._handle, omega_m)
        
        @property
        def omega_l(self):
            """
            Element omega_l ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _visualisation.f90wrap_sim_info__get__omega_l(self._handle)
        
        @omega_l.setter
        def omega_l(self, omega_l):
            _visualisation.f90wrap_sim_info__set__omega_l(self._handle, omega_l)
        
        @property
        def omega_k(self):
            """
            Element omega_k ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _visualisation.f90wrap_sim_info__get__omega_k(self._handle)
        
        @omega_k.setter
        def omega_k(self, omega_k):
            _visualisation.f90wrap_sim_info__set__omega_k(self._handle, omega_k)
        
        @property
        def omega_b(self):
            """
            Element omega_b ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 53
            
            """
            return _visualisation.f90wrap_sim_info__get__omega_b(self._handle)
        
        @omega_b.setter
        def omega_b(self, omega_b):
            _visualisation.f90wrap_sim_info__set__omega_b(self._handle, omega_b)
        
        @property
        def time_tot(self):
            """
            Element time_tot ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 54
            
            """
            return _visualisation.f90wrap_sim_info__get__time_tot(self._handle)
        
        @time_tot.setter
        def time_tot(self, time_tot):
            _visualisation.f90wrap_sim_info__set__time_tot(self._handle, time_tot)
        
        @property
        def time_simu(self):
            """
            Element time_simu ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 54
            
            """
            return _visualisation.f90wrap_sim_info__get__time_simu(self._handle)
        
        @time_simu.setter
        def time_simu(self, time_simu):
            _visualisation.f90wrap_sim_info__set__time_simu(self._handle, time_simu)
        
        @property
        def redshift(self):
            """
            Element redshift ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 54
            
            """
            return _visualisation.f90wrap_sim_info__get__redshift(self._handle)
        
        @redshift.setter
        def redshift(self, redshift):
            _visualisation.f90wrap_sim_info__set__redshift(self._handle, redshift)
        
        @property
        def t2(self):
            """
            Element t2 ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 54
            
            """
            return _visualisation.f90wrap_sim_info__get__t2(self._handle)
        
        @t2.setter
        def t2(self, t2):
            _visualisation.f90wrap_sim_info__set__t2(self._handle, t2)
        
        @property
        def nh(self):
            """
            Element nh ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 54
            
            """
            return _visualisation.f90wrap_sim_info__get__nh(self._handle)
        
        @nh.setter
        def nh(self, nh):
            _visualisation.f90wrap_sim_info__set__nh(self._handle, nh)
        
        @property
        def n_frw(self):
            """
            Element n_frw ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 55
            
            """
            return _visualisation.f90wrap_sim_info__get__n_frw(self._handle)
        
        @n_frw.setter
        def n_frw(self, n_frw):
            _visualisation.f90wrap_sim_info__set__n_frw(self._handle, n_frw)
        
        @property
        def aexp_frw(self):
            """
            Element aexp_frw ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 56
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _visualisation.f90wrap_sim_info__array__aexp_frw(self._handle)
            if array_handle in self._arrays:
                aexp_frw = self._arrays[array_handle]
            else:
                aexp_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_sim_info__array__aexp_frw)
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
                _visualisation.f90wrap_sim_info__array__hexp_frw(self._handle)
            if array_handle in self._arrays:
                hexp_frw = self._arrays[array_handle]
            else:
                hexp_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_sim_info__array__hexp_frw)
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
                _visualisation.f90wrap_sim_info__array__tau_frw(self._handle)
            if array_handle in self._arrays:
                tau_frw = self._arrays[array_handle]
            else:
                tau_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_sim_info__array__tau_frw)
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
                _visualisation.f90wrap_sim_info__array__t_frw(self._handle)
            if array_handle in self._arrays:
                t_frw = self._arrays[array_handle]
            else:
                t_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_sim_info__array__t_frw)
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
            return _visualisation.f90wrap_sim_info__get__eta_sn(self._handle)
        
        @eta_sn.setter
        def eta_sn(self, eta_sn):
            _visualisation.f90wrap_sim_info__set__eta_sn(self._handle, eta_sn)
        
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
        
    
    @f90wrap.runtime.register_class("visualisation.level")
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
            result = _visualisation.f90wrap_level_initialise()
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
                _visualisation.f90wrap_level_finalise(this=self._handle)
        
        @property
        def ilevel(self):
            """
            Element ilevel ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 60
            
            """
            return _visualisation.f90wrap_level__get__ilevel(self._handle)
        
        @ilevel.setter
        def ilevel(self, ilevel):
            _visualisation.f90wrap_level__set__ilevel(self._handle, ilevel)
        
        @property
        def ngrid(self):
            """
            Element ngrid ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 61
            
            """
            return _visualisation.f90wrap_level__get__ngrid(self._handle)
        
        @ngrid.setter
        def ngrid(self, ngrid):
            _visualisation.f90wrap_level__set__ngrid(self._handle, ngrid)
        
        @property
        def ind_grid(self):
            """
            Element ind_grid ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 62
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _visualisation.f90wrap_level__array__ind_grid(self._handle)
            if array_handle in self._arrays:
                ind_grid = self._arrays[array_handle]
            else:
                ind_grid = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_level__array__ind_grid)
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
                _visualisation.f90wrap_level__array__real_ind(self._handle)
            if array_handle in self._arrays:
                real_ind = self._arrays[array_handle]
            else:
                real_ind = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_level__array__real_ind)
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
                _visualisation.f90wrap_level__array__xg(self._handle)
            if array_handle in self._arrays:
                xg = self._arrays[array_handle]
            else:
                xg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_level__array__xg)
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
                _visualisation.f90wrap_level__array__cube(self._handle)
            if array_handle in self._arrays:
                cube = self._arrays[array_handle]
            else:
                cube = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_level__array__cube)
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
                _visualisation.f90wrap_level__array__map(self._handle)
            if array_handle in self._arrays:
                map = self._arrays[array_handle]
            else:
                map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_level__array__map)
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
                _visualisation.f90wrap_level__array__rho(self._handle)
            if array_handle in self._arrays:
                rho = self._arrays[array_handle]
            else:
                rho = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_level__array__rho)
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
            return _visualisation.f90wrap_level__get__imin(self._handle)
        
        @imin.setter
        def imin(self, imin):
            _visualisation.f90wrap_level__set__imin(self._handle, imin)
        
        @property
        def imax(self):
            """
            Element imax ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 69
            
            """
            return _visualisation.f90wrap_level__get__imax(self._handle)
        
        @imax.setter
        def imax(self, imax):
            _visualisation.f90wrap_level__set__imax(self._handle, imax)
        
        @property
        def jmin(self):
            """
            Element jmin ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 70
            
            """
            return _visualisation.f90wrap_level__get__jmin(self._handle)
        
        @jmin.setter
        def jmin(self, jmin):
            _visualisation.f90wrap_level__set__jmin(self._handle, jmin)
        
        @property
        def jmax(self):
            """
            Element jmax ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 71
            
            """
            return _visualisation.f90wrap_level__get__jmax(self._handle)
        
        @jmax.setter
        def jmax(self, jmax):
            _visualisation.f90wrap_level__set__jmax(self._handle, jmax)
        
        @property
        def kmin(self):
            """
            Element kmin ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 72
            
            """
            return _visualisation.f90wrap_level__get__kmin(self._handle)
        
        @kmin.setter
        def kmin(self, kmin):
            _visualisation.f90wrap_level__set__kmin(self._handle, kmin)
        
        @property
        def kmax(self):
            """
            Element kmax ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 73
            
            """
            return _visualisation.f90wrap_level__get__kmax(self._handle)
        
        @kmax.setter
        def kmax(self, kmax):
            _visualisation.f90wrap_level__set__kmax(self._handle, kmax)
        
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
        
    
    @f90wrap.runtime.register_class("visualisation.data_handler")
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
            result = _visualisation.f90wrap_data_handler_initialise()
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
                _visualisation.f90wrap_data_handler_finalise(this=self._handle)
        
        @property
        def name(self):
            """
            Element name ftype=character(80) pytype=str
            
            
            Defined at read_amr_module.fpp line 76
            
            """
            return _visualisation.f90wrap_data_handler__get__name(self._handle)
        
        @name.setter
        def name(self, name):
            _visualisation.f90wrap_data_handler__set__name(self._handle, name)
        
        @property
        def x_data(self):
            """
            Element x_data ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 77
            
            """
            return _visualisation.f90wrap_data_handler__get__x_data(self._handle)
        
        @x_data.setter
        def x_data(self, x_data):
            _visualisation.f90wrap_data_handler__set__x_data(self._handle, x_data)
        
        @property
        def y_data(self):
            """
            Element y_data ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 77
            
            """
            return _visualisation.f90wrap_data_handler__get__y_data(self._handle)
        
        @y_data.setter
        def y_data(self, y_data):
            _visualisation.f90wrap_data_handler__set__y_data(self._handle, y_data)
        
        @property
        def z_data(self):
            """
            Element z_data ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 77
            
            """
            return _visualisation.f90wrap_data_handler__get__z_data(self._handle)
        
        @z_data.setter
        def z_data(self, z_data):
            _visualisation.f90wrap_data_handler__set__z_data(self._handle, z_data)
        
        @property
        def nx(self):
            """
            Element nx ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 78
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _visualisation.f90wrap_data_handler__array__nx(self._handle)
            if array_handle in self._arrays:
                nx = self._arrays[array_handle]
            else:
                nx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_data_handler__array__nx)
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
                _visualisation.f90wrap_data_handler__array__ny(self._handle)
            if array_handle in self._arrays:
                ny = self._arrays[array_handle]
            else:
                ny = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_data_handler__array__ny)
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
                _visualisation.f90wrap_data_handler__array__nz(self._handle)
            if array_handle in self._arrays:
                nz = self._arrays[array_handle]
            else:
                nz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_data_handler__array__nz)
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
                _visualisation.f90wrap_data_handler__array__x(self._handle)
            if array_handle in self._arrays:
                x = self._arrays[array_handle]
            else:
                x = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_data_handler__array__x)
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
                _visualisation.f90wrap_data_handler__array__y(self._handle)
            if array_handle in self._arrays:
                y = self._arrays[array_handle]
            else:
                y = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_data_handler__array__y)
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
                _visualisation.f90wrap_data_handler__array__z(self._handle)
            if array_handle in self._arrays:
                z = self._arrays[array_handle]
            else:
                z = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_data_handler__array__z)
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
        
    
    @f90wrap.runtime.register_class("visualisation.particle")
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
            result = _visualisation.f90wrap_particle_initialise()
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
                _visualisation.f90wrap_particle_finalise(this=self._handle)
        
        @property
        def id(self):
            """
            Element id ftype=integer(irg) pytype=int
            
            
            Defined at read_amr_module.fpp line 82
            
            """
            return _visualisation.f90wrap_particle__get__id(self._handle)
        
        @id.setter
        def id(self, id):
            _visualisation.f90wrap_particle__set__id(self._handle, id)
        
        @property
        def x(self):
            """
            Element x ftype=type(vector) pytype=Vector
            
            
            Defined at read_amr_module.fpp line 83
            
            """
            x_handle = _visualisation.f90wrap_particle__get__x(self._handle)
            if tuple(x_handle) in self._objs:
                x = self._objs[tuple(x_handle)]
            else:
                x = vectors.vector.from_handle(x_handle)
                self._objs[tuple(x_handle)] = x
            return x
        
        @x.setter
        def x(self, x):
            x = x._handle
            _visualisation.f90wrap_particle__set__x(self._handle, x)
        
        @property
        def v(self):
            """
            Element v ftype=type(vector) pytype=Vector
            
            
            Defined at read_amr_module.fpp line 83
            
            """
            v_handle = _visualisation.f90wrap_particle__get__v(self._handle)
            if tuple(v_handle) in self._objs:
                v = self._objs[tuple(v_handle)]
            else:
                v = vectors.vector.from_handle(v_handle)
                self._objs[tuple(v_handle)] = v
            return v
        
        @v.setter
        def v(self, v):
            v = v._handle
            _visualisation.f90wrap_particle__set__v(self._handle, v)
        
        @property
        def m(self):
            """
            Element m ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 84
            
            """
            return _visualisation.f90wrap_particle__get__m(self._handle)
        
        @m.setter
        def m(self, m):
            _visualisation.f90wrap_particle__set__m(self._handle, m)
        
        @property
        def met(self):
            """
            Element met ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 84
            
            """
            return _visualisation.f90wrap_particle__get__met(self._handle)
        
        @met.setter
        def met(self, met):
            _visualisation.f90wrap_particle__set__met(self._handle, met)
        
        @property
        def imass(self):
            """
            Element imass ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 84
            
            """
            return _visualisation.f90wrap_particle__get__imass(self._handle)
        
        @imass.setter
        def imass(self, imass):
            _visualisation.f90wrap_particle__set__imass(self._handle, imass)
        
        @property
        def age(self):
            """
            Element age ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 84
            
            """
            return _visualisation.f90wrap_particle__get__age(self._handle)
        
        @age.setter
        def age(self, age):
            _visualisation.f90wrap_particle__set__age(self._handle, age)
        
        @property
        def tform(self):
            """
            Element tform ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 84
            
            """
            return _visualisation.f90wrap_particle__get__tform(self._handle)
        
        @tform.setter
        def tform(self, tform):
            _visualisation.f90wrap_particle__set__tform(self._handle, tform)
        
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
    def retrieve_vars(repository, myvars):
        """
        retrieve_vars(repository, myvars)
        
        
        Defined at read_amr_module.fpp lines 91-96
        
        Parameters
        ----------
        repository : str
        myvars : Hydroid
        
        """
        _visualisation.f90wrap_retrieve_vars(repository=repository, \
            myvars=myvars._handle)
    
    @staticmethod
    def title(n, nchar):
        """
        title(n, nchar)
        
        
        Defined at read_amr_module.fpp lines 103-127
        
        Parameters
        ----------
        n : int
        nchar : str
        
        """
        _visualisation.f90wrap_title(n=n, nchar=nchar)
    
    @staticmethod
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
        _visualisation.f90wrap_hilbert3d(x=x, y=y, z=z, order=order, \
            bit_length=bit_length, npoint=npoint)
    
    @staticmethod
    def check_lmax(ngridfile):
        """
        check_lmax(ngridfile)
        
        
        Defined at read_amr_module.fpp lines 214-226
        
        Parameters
        ----------
        ngridfile : int array
        
        """
        _visualisation.f90wrap_check_lmax(ngridfile=ngridfile)
    
    @staticmethod
    def check_families(repository):
        """
        check_families(repository)
        
        
        Defined at read_amr_module.fpp lines 235-260
        
        Parameters
        ----------
        repository : str
        
        """
        _visualisation.f90wrap_check_families(repository=repository)
    
    @staticmethod
    def read_hydrofile_descriptor(repository):
        """
        read_hydrofile_descriptor(repository)
        
        
        Defined at read_amr_module.fpp lines 269-303
        
        Parameters
        ----------
        repository : str
        
        """
        _visualisation.f90wrap_read_hydrofile_descriptor(repository=repository)
    
    @staticmethod
    def read_hydrofile_descriptor_old(repository):
        """
        read_hydrofile_descriptor_old(repository)
        
        
        Defined at read_amr_module.fpp lines 313-369
        
        Parameters
        ----------
        repository : str
        
        """
        _visualisation.f90wrap_read_hydrofile_descriptor_old(repository=repository)
    
    @staticmethod
    def select_from_descriptor_ids(newvar, newid):
        """
        select_from_descriptor_ids(newvar, newid)
        
        
        Defined at read_amr_module.fpp lines 371-472
        
        Parameters
        ----------
        newvar : str
        newid : int
        
        """
        _visualisation.f90wrap_select_from_descriptor_ids(newvar=newvar, newid=newid)
    
    @staticmethod
    def read_hydrofile_descriptor_new(repository):
        """
        read_hydrofile_descriptor_new(repository)
        
        
        Defined at read_amr_module.fpp lines 482-513
        
        Parameters
        ----------
        repository : str
        
        """
        _visualisation.f90wrap_read_hydrofile_descriptor_new(repository=repository)
    
    @staticmethod
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
        _visualisation.f90wrap_getvarvalue(reg=self._handle, dx=dx, x=x._handle, \
            var=var, son=son, varname=varname, value=value, trans_matrix=trans_matrix, \
            grav_var=grav_var)
    
    @staticmethod
    def get_eta_sn(repository):
        """
        get_eta_sn(repository)
        
        
        Defined at read_amr_module.fpp lines 1176-1206
        
        Parameters
        ----------
        repository : str
        
        """
        _visualisation.f90wrap_get_eta_sn(repository=repository)
    
    @staticmethod
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
        sf_eff = _visualisation.f90wrap_sf_eff(reg=self._handle, dx=dx, x=x._handle, \
            var=var, star_maker=star_maker)
        return sf_eff
    
    @staticmethod
    def init_amr_read(repository):
        """
        init_amr_read(repository)
        
        
        Defined at read_amr_module.fpp lines 1781-1883
        
        Parameters
        ----------
        repository : str
        
        """
        _visualisation.f90wrap_init_amr_read(repository=repository)
    
    @staticmethod
    def get_cpu_map(self):
        """
        get_cpu_map(self)
        
        
        Defined at read_amr_module.fpp lines 1891-1979
        
        Parameters
        ----------
        reg : Region
        
        """
        _visualisation.f90wrap_get_cpu_map(reg=self._handle)
    
    @staticmethod
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
        _visualisation.f90wrap_getnborgrids(son=son, nbor=nbor, igrid=igrid, \
            igridn=igridn, ngrid=ngrid)
    
    @staticmethod
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
        _visualisation.f90wrap_getnborcells(igridn=igridn, ind=ind, icelln=icelln, \
            ng=ng)
    
    @staticmethod
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
        _visualisation.f90wrap_getnbor(son=son, nbor=nbor, ind_cell=ind_cell, \
            ind_father=ind_father, ncell=ncell)
    
    @staticmethod
    def getparttype(self, ptype):
        """
        getparttype(self, ptype)
        
        
        Defined at read_amr_module.fpp lines 2103-2111
        
        Parameters
        ----------
        part : Particle
        ptype : str
        
        """
        _visualisation.f90wrap_getparttype(part=self._handle, ptype=ptype)
    
    @staticmethod
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
        _visualisation.f90wrap_getpartvalue(reg=self._handle, part=part._handle, \
            var=var, value=value, dx=None if dx is None else dx._handle)
    
    _dt_array_initialisers = []
    

io_ramses = Io_Ramses()

class Filtering(f90wrap.runtime.FortranModule):
    """
    Module filtering
    
    
    Defined at read_amr_module.fpp lines 2516-2648
    
    """
    @f90wrap.runtime.register_class("visualisation.filter")
    class filter(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=filter)
        
        
        Defined at read_amr_module.fpp lines 2519-2526
        
        """
        def __init__(self, handle=None):
            """
            self = Filter()
            
            
            Defined at read_amr_module.fpp lines 2519-2526
            
            
            Returns
            -------
            this : Filter
            	Object to be constructed
            
            
            Automatically generated constructor for filter
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _visualisation.f90wrap_filter_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Filter
            
            
            Defined at read_amr_module.fpp lines 2519-2526
            
            Parameters
            ----------
            this : Filter
            	Object to be destructed
            
            
            Automatically generated destructor for filter
            """
            if self._alloc:
                _visualisation.f90wrap_filter_finalise(this=self._handle)
        
        @property
        def name(self):
            """
            Element name ftype=character(128) pytype=str
            
            
            Defined at read_amr_module.fpp line 2520
            
            """
            return _visualisation.f90wrap_filter__get__name(self._handle)
        
        @name.setter
        def name(self, name):
            _visualisation.f90wrap_filter__set__name(self._handle, name)
        
        @property
        def ncond(self):
            """
            Element ncond ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 2521
            
            """
            return _visualisation.f90wrap_filter__get__ncond(self._handle)
        
        @ncond.setter
        def ncond(self, ncond):
            _visualisation.f90wrap_filter__set__ncond(self._handle, ncond)
        
        @property
        def cond_vars(self):
            """
            Element cond_vars ftype=character(128) pytype=str
            
            
            Defined at read_amr_module.fpp line 2522
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _visualisation.f90wrap_filter__array__cond_vars(self._handle)
            if array_handle in self._arrays:
                cond_vars = self._arrays[array_handle]
            else:
                cond_vars = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_filter__array__cond_vars)
                self._arrays[array_handle] = cond_vars
            return cond_vars
        
        @cond_vars.setter
        def cond_vars(self, cond_vars):
            self.cond_vars[...] = cond_vars
        
        @property
        def cond_vars_comp(self):
            """
            Element cond_vars_comp ftype=character(128) pytype=str
            
            
            Defined at read_amr_module.fpp line 2523
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _visualisation.f90wrap_filter__array__cond_vars_comp(self._handle)
            if array_handle in self._arrays:
                cond_vars_comp = self._arrays[array_handle]
            else:
                cond_vars_comp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_filter__array__cond_vars_comp)
                self._arrays[array_handle] = cond_vars_comp
            return cond_vars_comp
        
        @cond_vars_comp.setter
        def cond_vars_comp(self, cond_vars_comp):
            self.cond_vars_comp[...] = cond_vars_comp
        
        @property
        def cond_ops(self):
            """
            Element cond_ops ftype=character(2) pytype=str
            
            
            Defined at read_amr_module.fpp line 2524
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _visualisation.f90wrap_filter__array__cond_ops(self._handle)
            if array_handle in self._arrays:
                cond_ops = self._arrays[array_handle]
            else:
                cond_ops = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_filter__array__cond_ops)
                self._arrays[array_handle] = cond_ops
            return cond_ops
        
        @cond_ops.setter
        def cond_ops(self, cond_ops):
            self.cond_ops[...] = cond_ops
        
        @property
        def cond_vals(self):
            """
            Element cond_vals ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 2525
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _visualisation.f90wrap_filter__array__cond_vals(self._handle)
            if array_handle in self._arrays:
                cond_vals = self._arrays[array_handle]
            else:
                cond_vals = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_filter__array__cond_vals)
                self._arrays[array_handle] = cond_vals
            return cond_vals
        
        @cond_vals.setter
        def cond_vals(self, cond_vals):
            self.cond_vals[...] = cond_vals
        
        @property
        def use_var(self):
            """
            Element use_var ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 2526
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _visualisation.f90wrap_filter__array__use_var(self._handle)
            if array_handle in self._arrays:
                use_var = self._arrays[array_handle]
            else:
                use_var = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_filter__array__use_var)
                self._arrays[array_handle] = use_var
            return use_var
        
        @use_var.setter
        def use_var(self, use_var):
            self.use_var[...] = use_var
        
        def __str__(self):
            ret = ['<filter>{\n']
            ret.append('    name : ')
            ret.append(repr(self.name))
            ret.append(',\n    ncond : ')
            ret.append(repr(self.ncond))
            ret.append(',\n    cond_vars : ')
            ret.append(repr(self.cond_vars))
            ret.append(',\n    cond_vars_comp : ')
            ret.append(repr(self.cond_vars_comp))
            ret.append(',\n    cond_ops : ')
            ret.append(repr(self.cond_ops))
            ret.append(',\n    cond_vals : ')
            ret.append(repr(self.cond_vals))
            ret.append(',\n    use_var : ')
            ret.append(repr(self.use_var))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def allocate_filter(self):
        """
        allocate_filter(self)
        
        
        Defined at read_amr_module.fpp lines 2529-2537
        
        Parameters
        ----------
        filt : Filter
        
        """
        _visualisation.f90wrap_allocate_filter(filt=self._handle)
    
    @staticmethod
    def cond_string_to_filter(str, filt):
        """
        cond_string_to_filter(str, filt)
        
        
        Defined at read_amr_module.fpp lines 2539-2544
        
        Parameters
        ----------
        str : str
        filt : Filter
        
        """
        _visualisation.f90wrap_cond_string_to_filter(str=str, filt=filt._handle)
    
    @staticmethod
    def filter_cell(self, filt, cell_x, cell_dx, cell_var, cell_son, trans_matrix, \
        grav_var=None):
        """
        filter_cell = filter_cell(self, filt, cell_x, cell_dx, cell_var, cell_son, \
            trans_matrix[, grav_var])
        
        
        Defined at read_amr_module.fpp lines 2546-2599
        
        Parameters
        ----------
        reg : Region
        filt : Filter
        cell_x : Vector
        cell_dx : float
        cell_var : float array
        cell_son : int array
        trans_matrix : float array
        grav_var : float array
        
        Returns
        -------
        filter_cell : bool
        
        """
        filter_cell = _visualisation.f90wrap_filter_cell(reg=self._handle, \
            filt=filt._handle, cell_x=cell_x._handle, cell_dx=cell_dx, \
            cell_var=cell_var, cell_son=cell_son, trans_matrix=trans_matrix, \
            grav_var=grav_var)
        return filter_cell
    
    @staticmethod
    def filter_particle(self, filt, part, dx=None):
        """
        filter_particle = filter_particle(self, filt, part[, dx])
        
        
        Defined at read_amr_module.fpp lines 2601-2635
        
        Parameters
        ----------
        reg : Region
        filt : Filter
        part : Particle
        dx : Vector
        
        Returns
        -------
        filter_particle : bool
        
        """
        filter_particle = _visualisation.f90wrap_filter_particle(reg=self._handle, \
            filt=filt._handle, part=part._handle, dx=None if dx is None else dx._handle)
        return filter_particle
    
    @staticmethod
    def filter_sub(self, cell_x):
        """
        filter_sub = filter_sub(self, cell_x)
        
        
        Defined at read_amr_module.fpp lines 2637-2648
        
        Parameters
        ----------
        sub : Region
        cell_x : float array
        
        Returns
        -------
        filter_sub : bool
        
        """
        filter_sub = _visualisation.f90wrap_filter_sub(sub=self._handle, cell_x=cell_x)
        return filter_sub
    
    _dt_array_initialisers = []
    

filtering = Filtering()

class Obs_Instruments(f90wrap.runtime.FortranModule):
    """
    Module obs_instruments
    
    
    Defined at ramses2map.fpp lines 5-207
    
    """
    @f90wrap.runtime.register_class("visualisation.camera")
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
            result = _visualisation.f90wrap_camera_initialise()
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
                _visualisation.f90wrap_camera_finalise(this=self._handle)
        
        @property
        def centre(self):
            """
            Element centre ftype=type(vector) pytype=Vector
            
            
            Defined at ramses2map.fpp line 12
            
            """
            centre_handle = _visualisation.f90wrap_camera__get__centre(self._handle)
            if tuple(centre_handle) in self._objs:
                centre = self._objs[tuple(centre_handle)]
            else:
                centre = vectors.vector.from_handle(centre_handle)
                self._objs[tuple(centre_handle)] = centre
            return centre
        
        @centre.setter
        def centre(self, centre):
            centre = centre._handle
            _visualisation.f90wrap_camera__set__centre(self._handle, centre)
        
        @property
        def los_axis(self):
            """
            Element los_axis ftype=type(vector) pytype=Vector
            
            
            Defined at ramses2map.fpp line 12
            
            """
            los_axis_handle = _visualisation.f90wrap_camera__get__los_axis(self._handle)
            if tuple(los_axis_handle) in self._objs:
                los_axis = self._objs[tuple(los_axis_handle)]
            else:
                los_axis = vectors.vector.from_handle(los_axis_handle)
                self._objs[tuple(los_axis_handle)] = los_axis
            return los_axis
        
        @los_axis.setter
        def los_axis(self, los_axis):
            los_axis = los_axis._handle
            _visualisation.f90wrap_camera__set__los_axis(self._handle, los_axis)
        
        @property
        def up_vector(self):
            """
            Element up_vector ftype=type(vector) pytype=Vector
            
            
            Defined at ramses2map.fpp line 12
            
            """
            up_vector_handle = _visualisation.f90wrap_camera__get__up_vector(self._handle)
            if tuple(up_vector_handle) in self._objs:
                up_vector = self._objs[tuple(up_vector_handle)]
            else:
                up_vector = vectors.vector.from_handle(up_vector_handle)
                self._objs[tuple(up_vector_handle)] = up_vector
            return up_vector
        
        @up_vector.setter
        def up_vector(self, up_vector):
            up_vector = up_vector._handle
            _visualisation.f90wrap_camera__set__up_vector(self._handle, up_vector)
        
        @property
        def region_axis(self):
            """
            Element region_axis ftype=type(vector) pytype=Vector
            
            
            Defined at ramses2map.fpp line 12
            
            """
            region_axis_handle = \
                _visualisation.f90wrap_camera__get__region_axis(self._handle)
            if tuple(region_axis_handle) in self._objs:
                region_axis = self._objs[tuple(region_axis_handle)]
            else:
                region_axis = vectors.vector.from_handle(region_axis_handle)
                self._objs[tuple(region_axis_handle)] = region_axis
            return region_axis
        
        @region_axis.setter
        def region_axis(self, region_axis):
            region_axis = region_axis._handle
            _visualisation.f90wrap_camera__set__region_axis(self._handle, region_axis)
        
        @property
        def region_size(self):
            """
            Element region_size ftype=real(dbl) pytype=float
            
            
            Defined at ramses2map.fpp line 13
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _visualisation.f90wrap_camera__array__region_size(self._handle)
            if array_handle in self._arrays:
                region_size = self._arrays[array_handle]
            else:
                region_size = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_camera__array__region_size)
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
            return _visualisation.f90wrap_camera__get__distance(self._handle)
        
        @distance.setter
        def distance(self, distance):
            _visualisation.f90wrap_camera__set__distance(self._handle, distance)
        
        @property
        def far_cut_depth(self):
            """
            Element far_cut_depth ftype=real(dbl) pytype=float
            
            
            Defined at ramses2map.fpp line 14
            
            """
            return _visualisation.f90wrap_camera__get__far_cut_depth(self._handle)
        
        @far_cut_depth.setter
        def far_cut_depth(self, far_cut_depth):
            _visualisation.f90wrap_camera__set__far_cut_depth(self._handle, far_cut_depth)
        
        @property
        def map_max_size(self):
            """
            Element map_max_size ftype=integer  pytype=int
            
            
            Defined at ramses2map.fpp line 15
            
            """
            return _visualisation.f90wrap_camera__get__map_max_size(self._handle)
        
        @map_max_size.setter
        def map_max_size(self, map_max_size):
            _visualisation.f90wrap_camera__set__map_max_size(self._handle, map_max_size)
        
        @property
        def nfilter(self):
            """
            Element nfilter ftype=integer  pytype=int
            
            
            Defined at ramses2map.fpp line 16
            
            """
            return _visualisation.f90wrap_camera__get__nfilter(self._handle)
        
        @nfilter.setter
        def nfilter(self, nfilter):
            _visualisation.f90wrap_camera__set__nfilter(self._handle, nfilter)
        
        @property
        def nsubs(self):
            """
            Element nsubs ftype=integer  pytype=int
            
            
            Defined at ramses2map.fpp line 16
            
            """
            return _visualisation.f90wrap_camera__get__nsubs(self._handle)
        
        @nsubs.setter
        def nsubs(self, nsubs):
            _visualisation.f90wrap_camera__set__nsubs(self._handle, nsubs)
        
        @property
        def lmin(self):
            """
            Element lmin ftype=integer  pytype=int
            
            
            Defined at ramses2map.fpp line 17
            
            """
            return _visualisation.f90wrap_camera__get__lmin(self._handle)
        
        @lmin.setter
        def lmin(self, lmin):
            _visualisation.f90wrap_camera__set__lmin(self._handle, lmin)
        
        @property
        def lmax(self):
            """
            Element lmax ftype=integer  pytype=int
            
            
            Defined at ramses2map.fpp line 17
            
            """
            return _visualisation.f90wrap_camera__get__lmax(self._handle)
        
        @lmax.setter
        def lmax(self, lmax):
            _visualisation.f90wrap_camera__set__lmax(self._handle, lmax)
        
        def init_array_filters(self):
            self.filters = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _visualisation.f90wrap_camera__array_getitem__filters,
                                            _visualisation.f90wrap_camera__array_setitem__filters,
                                            _visualisation.f90wrap_camera__array_len__filters,
                                            """
            Element filters ftype=type(filter) pytype=Filter
            
            
            Defined at ramses2map.fpp line 18
            
            """, Filtering.filter)
            return self.filters
        
        def init_array_subs(self):
            self.subs = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _visualisation.f90wrap_camera__array_getitem__subs,
                                            _visualisation.f90wrap_camera__array_setitem__subs,
                                            _visualisation.f90wrap_camera__array_len__subs,
                                            """
            Element subs ftype=type(region) pytype=Region
            
            
            Defined at ramses2map.fpp line 19
            
            """, Geometrical_Regions.region)
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
        
    
    @staticmethod
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
        log2 = _visualisation.f90wrap_log2(x=x)
        return log2
    
    @staticmethod
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
        init_camera = _visualisation.f90wrap_init_camera(centre=self._handle, \
            los_axis=los_axis._handle, up_vector=up_vector._handle, \
            region_size=region_size, region_axis=region_axis._handle, distance=distance, \
            far_cut_depth=far_cut_depth, map_max_size=map_max_size, nfilter=nfilter, \
            nsubs=nsubs)
        init_camera = \
            f90wrap.runtime.lookup_class("visualisation.camera").from_handle(init_camera, \
            alloc=True)
        return init_camera
    
    @staticmethod
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
            _visualisation.f90wrap_get_required_resolution(cam=self._handle)
        return get_required_resolution
    
    @staticmethod
    def get_map_size(self, n_map):
        """
        get_map_size(self, n_map)
        
        
        Defined at ramses2map.fpp lines 55-67
        
        Parameters
        ----------
        cam : Camera
        n_map : int array
        
        """
        _visualisation.f90wrap_get_map_size(cam=self._handle, n_map=n_map)
    
    @staticmethod
    def get_map_box(self, box):
        """
        get_map_box(self, box)
        
        
        Defined at ramses2map.fpp lines 69-79
        
        Parameters
        ----------
        cam : Camera
        box : Region
        
        """
        _visualisation.f90wrap_get_map_box(cam=self._handle, box=box._handle)
    
    @staticmethod
    def get_camera_basis(self, cam_basis):
        """
        get_camera_basis(self, cam_basis)
        
        
        Defined at ramses2map.fpp lines 101-109
        
        Parameters
        ----------
        cam : Camera
        cam_basis : Basis
        
        """
        _visualisation.f90wrap_get_camera_basis(cam=self._handle, \
            cam_basis=cam_basis._handle)
    
    @staticmethod
    def los_transformation(self, trans_matrix):
        """
        los_transformation(self, trans_matrix)
        
        
        Defined at ramses2map.fpp lines 111-122
        
        Parameters
        ----------
        cam : Camera
        trans_matrix : float array
        
        """
        _visualisation.f90wrap_los_transformation(cam=self._handle, \
            trans_matrix=trans_matrix)
    
    @staticmethod
    def get_bounding_box(self, bbox):
        """
        get_bounding_box(self, bbox)
        
        
        Defined at ramses2map.fpp lines 124-168
        
        Parameters
        ----------
        cam : Camera
        bbox : Region
        
        """
        _visualisation.f90wrap_get_bounding_box(cam=self._handle, bbox=bbox._handle)
    
    @staticmethod
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
        _visualisation.f90wrap_deproject_points(cam=self._handle, npoints=npoints, \
            points=points)
    
    @staticmethod
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
        _visualisation.f90wrap_project_points(cam=self._handle, npoints=npoints, \
            points=points)
    
    _dt_array_initialisers = []
    

obs_instruments = Obs_Instruments()

class Maps(f90wrap.runtime.FortranModule):
    """
    Module maps
    
    
    Defined at ramses2map.fpp lines 209-1733
    
    """
    @f90wrap.runtime.register_class("visualisation.projection_handler")
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
            result = _visualisation.f90wrap_projection_handler_initialise()
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
                _visualisation.f90wrap_projection_handler_finalise(this=self._handle)
        
        @property
        def pov(self):
            """
            Element pov ftype=character(128) pytype=str
            
            
            Defined at ramses2map.fpp line 218
            
            """
            return _visualisation.f90wrap_projection_handler__get__pov(self._handle)
        
        @pov.setter
        def pov(self, pov):
            _visualisation.f90wrap_projection_handler__set__pov(self._handle, pov)
        
        @property
        def nvars(self):
            """
            Element nvars ftype=integer  pytype=int
            
            
            Defined at ramses2map.fpp line 219
            
            """
            return _visualisation.f90wrap_projection_handler__get__nvars(self._handle)
        
        @nvars.setter
        def nvars(self, nvars):
            _visualisation.f90wrap_projection_handler__set__nvars(self._handle, nvars)
        
        @property
        def nfilter(self):
            """
            Element nfilter ftype=integer  pytype=int
            
            
            Defined at ramses2map.fpp line 219
            
            """
            return _visualisation.f90wrap_projection_handler__get__nfilter(self._handle)
        
        @nfilter.setter
        def nfilter(self, nfilter):
            _visualisation.f90wrap_projection_handler__set__nfilter(self._handle, nfilter)
        
        @property
        def varnames(self):
            """
            Element varnames ftype=character(128) pytype=str
            
            
            Defined at ramses2map.fpp line 220
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _visualisation.f90wrap_projection_handler__array__varnames(self._handle)
            if array_handle in self._arrays:
                varnames = self._arrays[array_handle]
            else:
                varnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_projection_handler__array__varnames)
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
            return _visualisation.f90wrap_projection_handler__get__weightvar(self._handle)
        
        @weightvar.setter
        def weightvar(self, weightvar):
            _visualisation.f90wrap_projection_handler__set__weightvar(self._handle, \
                weightvar)
        
        @property
        def toto(self):
            """
            Element toto ftype=real(dbl) pytype=float
            
            
            Defined at ramses2map.fpp line 222
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _visualisation.f90wrap_projection_handler__array__toto(self._handle)
            if array_handle in self._arrays:
                toto = self._arrays[array_handle]
            else:
                toto = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _visualisation.f90wrap_projection_handler__array__toto)
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
        
    
    @staticmethod
    def allocate_projection_handler(self):
        """
        allocate_projection_handler(self)
        
        
        Defined at ramses2map.fpp lines 225-228
        
        Parameters
        ----------
        proj : Projection_Handler
        
        """
        _visualisation.f90wrap_allocate_projection_handler(proj=self._handle)
    
    @staticmethod
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
        _visualisation.f90wrap_projection_hydro(repository=repository, cam=cam._handle, \
            bulk_velocity=bulk_velocity._handle, use_neigh=use_neigh, proj=proj._handle, \
            lmax=lmax, lmin=lmin)
    
    @staticmethod
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
        _visualisation.f90wrap_projection_parts(repository=repository, cam=cam._handle, \
            bulk_velocity=bulk_velocity._handle, proj=proj._handle, tag_file=tag_file, \
            inverse_tag=inverse_tag)
    
    @staticmethod
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
        _visualisation.f90wrap_project_particles(repository=repository, \
            bbox=bbox._handle, cam=cam._handle, proj=proj._handle, tag_file=tag_file, \
            inverse_tag=inverse_tag)
    
    @staticmethod
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
        _visualisation.f90wrap_healpix_hydro(repository=repository, reg=reg._handle, \
            nside=nside, proj=proj._handle)
    
    @staticmethod
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
        _visualisation.f90wrap_project_cells_hpix(repository=repository, \
            reg=reg._handle, nside=nside, proj=proj._handle)
    
    _dt_array_initialisers = []
    

maps = Maps()

