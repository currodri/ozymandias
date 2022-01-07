from __future__ import print_function, absolute_import, division
import _amr2
import f90wrap.runtime
import logging

class Vectors(f90wrap.runtime.FortranModule):
    """
    Module vectors
    
    
    Defined at linalg_module.fpp lines 23-150
    
    """
    @f90wrap.runtime.register_class("amr2.vector")
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
            result = _amr2.f90wrap_vector_initialise()
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
                _amr2.f90wrap_vector_finalise(this=self._handle)
        
        @property
        def x(self):
            """
            Element x ftype=real(dbl) pytype=float
            
            
            Defined at linalg_module.fpp line 27
            
            """
            return _amr2.f90wrap_vector__get__x(self._handle)
        
        @x.setter
        def x(self, x):
            _amr2.f90wrap_vector__set__x(self._handle, x)
        
        @property
        def y(self):
            """
            Element y ftype=real(dbl) pytype=float
            
            
            Defined at linalg_module.fpp line 28
            
            """
            return _amr2.f90wrap_vector__get__y(self._handle)
        
        @y.setter
        def y(self, y):
            _amr2.f90wrap_vector__set__y(self._handle, y)
        
        @property
        def z(self):
            """
            Element z ftype=real(dbl) pytype=float
            
            
            Defined at linalg_module.fpp line 29
            
            """
            return _amr2.f90wrap_vector__get__z(self._handle)
        
        @z.setter
        def z(self, z):
            _amr2.f90wrap_vector__set__z(self._handle, z)
        
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
        
    
    @f90wrap.runtime.register_class("amr2.array_vectors")
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
            result = _amr2.f90wrap_array_vectors_initialise()
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
                _amr2.f90wrap_array_vectors_finalise(this=self._handle)
        
        @property
        def n(self):
            """
            Element n ftype=integer  pytype=int
            
            
            Defined at linalg_module.fpp line 32
            
            """
            return _amr2.f90wrap_array_vectors__get__n(self._handle)
        
        @n.setter
        def n(self, n):
            _amr2.f90wrap_array_vectors__set__n(self._handle, n)
        
        def init_array_list(self):
            self.list = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _amr2.f90wrap_array_vectors__array_getitem__list,
                                            _amr2.f90wrap_array_vectors__array_setitem__list,
                                            _amr2.f90wrap_array_vectors__array_len__list,
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
        
        
        Defined at linalg_module.fpp lines 147-149
        
        Parameters
        ----------
        vec_1 : Vector
        
        Returns
        -------
        magnitude : float
        
        """
        magnitude = _amr2.f90wrap_magnitude(vec_1=self._handle)
        return magnitude
    
    @staticmethod
    def _array_to_vector(self, array):
        """
        _array_to_vector(self, array)
        
        
        Defined at linalg_module.fpp lines 61-66
        
        Parameters
        ----------
        vec_result : Vector
        array : float array
        
        """
        _amr2.f90wrap_array_to_vector(vec_result=self._handle, array=array)
    
    @staticmethod
    def _vector_to_array(array_result, vec_1):
        """
        _vector_to_array(array_result, vec_1)
        
        
        Defined at linalg_module.fpp lines 68-73
        
        Parameters
        ----------
        array_result : float array
        vec_1 : Vector
        
        """
        _amr2.f90wrap_vector_to_array(array_result=array_result, vec_1=vec_1._handle)
    
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
    
    
    Defined at linalg_module.fpp lines 167-240
    
    """
    @staticmethod
    def euler_matrix(r, dim, angle, ap):
        """
        euler_matrix(r, dim, angle, ap)
        
        
        Defined at linalg_module.fpp lines 184-217
        
        Parameters
        ----------
        r : float array
        dim : int
        angle : float
        ap : str
        
        """
        _amr2.f90wrap_euler_matrix(r=r, dim=dim, angle=angle, ap=ap)
    
    @staticmethod
    def _rotate_vector_single(self, rotation_matrix):
        """
        _rotate_vector_single(self, rotation_matrix)
        
        
        Defined at linalg_module.fpp lines 228-231
        
        Parameters
        ----------
        vec : Vector
        rotation_matrix : float array
        
        """
        _amr2.f90wrap_rotate_vector_single(vec=self._handle, \
            rotation_matrix=rotation_matrix)
    
    @staticmethod
    def _rotate_vector_array(self, rotation_matrix):
        """
        _rotate_vector_array(self, rotation_matrix)
        
        
        Defined at linalg_module.fpp lines 233-239
        
        Parameters
        ----------
        vec : Array_Vectors
        rotation_matrix : float array
        
        """
        _amr2.f90wrap_rotate_vector_array(vec=self._handle, \
            rotation_matrix=rotation_matrix)
    
    @staticmethod
    def rotate_vector(*args, **kwargs):
        """
        rotate_vector(*args, **kwargs)
        
        
        Defined at linalg_module.fpp lines 171-173
        
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
    
    
    Defined at linalg_module.fpp lines 257-290
    
    """
    @f90wrap.runtime.register_class("amr2.basis")
    class basis(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=basis)
        
        
        Defined at linalg_module.fpp lines 260-261
        
        """
        def __init__(self, handle=None):
            """
            self = Basis()
            
            
            Defined at linalg_module.fpp lines 260-261
            
            
            Returns
            -------
            this : Basis
            	Object to be constructed
            
            
            Automatically generated constructor for basis
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _amr2.f90wrap_basis_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Basis
            
            
            Defined at linalg_module.fpp lines 260-261
            
            Parameters
            ----------
            this : Basis
            	Object to be destructed
            
            
            Automatically generated destructor for basis
            """
            if self._alloc:
                _amr2.f90wrap_basis_finalise(this=self._handle)
        
        def init_array_u(self):
            self.u = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _amr2.f90wrap_basis__array_getitem__u,
                                            _amr2.f90wrap_basis__array_setitem__u,
                                            _amr2.f90wrap_basis__array_len__u,
                                            """
            Element u ftype=type(vector) pytype=Vector
            
            
            Defined at linalg_module.fpp line 261
            
            """, Vectors.vector)
            return self.u
        
        _dt_array_initialisers = [init_array_u]
        
    
    @staticmethod
    def initialise_basis(self):
        """
        initialise_basis(self)
        
        
        Defined at linalg_module.fpp lines 264-268
        
        Parameters
        ----------
        this : Basis
        
        """
        _amr2.f90wrap_initialise_basis(this=self._handle)
    
    @staticmethod
    def mgramschmidt(self, e):
        """
        mgramschmidt(self, e)
        
        
        Defined at linalg_module.fpp lines 277-289
        
        Parameters
        ----------
        vecs : Basis
        e : Basis
        
        """
        _amr2.f90wrap_mgramschmidt(vecs=self._handle, e=e._handle)
    
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
        r_sphere = _amr2.f90wrap_r_sphere(p=self._handle)
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
        theta_sphere = _amr2.f90wrap_theta_sphere(p=self._handle)
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
        phi_sphere = _amr2.f90wrap_phi_sphere(p=self._handle)
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
        r_cyl = _amr2.f90wrap_r_cyl(p=self._handle)
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
        phi_cyl = _amr2.f90wrap_phi_cyl(p=self._handle)
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
        _amr2.f90wrap_spherical_basis_from_cartesian(p=self._handle, \
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
        _amr2.f90wrap_cylindrical_basis_from_cartesian(p=self._handle, \
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
        _amr2.f90wrap_new_z_coordinates(axis=self._handle, \
            transformation_matrix=transformation_matrix, errormsg=errormsg)
    
    _dt_array_initialisers = []
    

coordinate_systems = Coordinate_Systems()

class Geometrical_Regions(f90wrap.runtime.FortranModule):
    """
    Module geometrical_regions
    
    
    Defined at coordinates_module.fpp lines 191-363
    
    """
    @f90wrap.runtime.register_class("amr2.region")
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
            result = _amr2.f90wrap_region_initialise()
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
                _amr2.f90wrap_region_finalise(this=self._handle)
        
        @property
        def name(self):
            """
            Element name ftype=character(128) pytype=str
            
            
            Defined at coordinates_module.fpp line 196
            
            """
            return _amr2.f90wrap_region__get__name(self._handle)
        
        @name.setter
        def name(self, name):
            _amr2.f90wrap_region__set__name(self._handle, name)
        
        @property
        def centre(self):
            """
            Element centre ftype=type(vector) pytype=Vector
            
            
            Defined at coordinates_module.fpp line 197
            
            """
            centre_handle = _amr2.f90wrap_region__get__centre(self._handle)
            if tuple(centre_handle) in self._objs:
                centre = self._objs[tuple(centre_handle)]
            else:
                centre = vectors.vector.from_handle(centre_handle)
                self._objs[tuple(centre_handle)] = centre
            return centre
        
        @centre.setter
        def centre(self, centre):
            centre = centre._handle
            _amr2.f90wrap_region__set__centre(self._handle, centre)
        
        @property
        def axis(self):
            """
            Element axis ftype=type(vector) pytype=Vector
            
            
            Defined at coordinates_module.fpp line 197
            
            """
            axis_handle = _amr2.f90wrap_region__get__axis(self._handle)
            if tuple(axis_handle) in self._objs:
                axis = self._objs[tuple(axis_handle)]
            else:
                axis = vectors.vector.from_handle(axis_handle)
                self._objs[tuple(axis_handle)] = axis
            return axis
        
        @axis.setter
        def axis(self, axis):
            axis = axis._handle
            _amr2.f90wrap_region__set__axis(self._handle, axis)
        
        @property
        def bulk_velocity(self):
            """
            Element bulk_velocity ftype=type(vector) pytype=Vector
            
            
            Defined at coordinates_module.fpp line 197
            
            """
            bulk_velocity_handle = _amr2.f90wrap_region__get__bulk_velocity(self._handle)
            if tuple(bulk_velocity_handle) in self._objs:
                bulk_velocity = self._objs[tuple(bulk_velocity_handle)]
            else:
                bulk_velocity = vectors.vector.from_handle(bulk_velocity_handle)
                self._objs[tuple(bulk_velocity_handle)] = bulk_velocity
            return bulk_velocity
        
        @bulk_velocity.setter
        def bulk_velocity(self, bulk_velocity):
            bulk_velocity = bulk_velocity._handle
            _amr2.f90wrap_region__set__bulk_velocity(self._handle, bulk_velocity)
        
        @property
        def xmin(self):
            """
            Element xmin ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _amr2.f90wrap_region__get__xmin(self._handle)
        
        @xmin.setter
        def xmin(self, xmin):
            _amr2.f90wrap_region__set__xmin(self._handle, xmin)
        
        @property
        def xmax(self):
            """
            Element xmax ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _amr2.f90wrap_region__get__xmax(self._handle)
        
        @xmax.setter
        def xmax(self, xmax):
            _amr2.f90wrap_region__set__xmax(self._handle, xmax)
        
        @property
        def ymin(self):
            """
            Element ymin ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _amr2.f90wrap_region__get__ymin(self._handle)
        
        @ymin.setter
        def ymin(self, ymin):
            _amr2.f90wrap_region__set__ymin(self._handle, ymin)
        
        @property
        def ymax(self):
            """
            Element ymax ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _amr2.f90wrap_region__get__ymax(self._handle)
        
        @ymax.setter
        def ymax(self, ymax):
            _amr2.f90wrap_region__set__ymax(self._handle, ymax)
        
        @property
        def zmin(self):
            """
            Element zmin ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _amr2.f90wrap_region__get__zmin(self._handle)
        
        @zmin.setter
        def zmin(self, zmin):
            _amr2.f90wrap_region__set__zmin(self._handle, zmin)
        
        @property
        def zmax(self):
            """
            Element zmax ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 198
            
            """
            return _amr2.f90wrap_region__get__zmax(self._handle)
        
        @zmax.setter
        def zmax(self, zmax):
            _amr2.f90wrap_region__set__zmax(self._handle, zmax)
        
        @property
        def rmin(self):
            """
            Element rmin ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 199
            
            """
            return _amr2.f90wrap_region__get__rmin(self._handle)
        
        @rmin.setter
        def rmin(self, rmin):
            _amr2.f90wrap_region__set__rmin(self._handle, rmin)
        
        @property
        def rmax(self):
            """
            Element rmax ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 199
            
            """
            return _amr2.f90wrap_region__get__rmax(self._handle)
        
        @rmax.setter
        def rmax(self, rmax):
            _amr2.f90wrap_region__set__rmax(self._handle, rmax)
        
        @property
        def angle(self):
            """
            Element angle ftype=real(dbl) pytype=float
            
            
            Defined at coordinates_module.fpp line 200
            
            """
            return _amr2.f90wrap_region__get__angle(self._handle)
        
        @angle.setter
        def angle(self, angle):
            _amr2.f90wrap_region__set__angle(self._handle, angle)
        
        @property
        def criteria_name(self):
            """
            Element criteria_name ftype=character(128) pytype=str
            
            
            Defined at coordinates_module.fpp line 201
            
            """
            return _amr2.f90wrap_region__get__criteria_name(self._handle)
        
        @criteria_name.setter
        def criteria_name(self, criteria_name):
            _amr2.f90wrap_region__set__criteria_name(self._handle, criteria_name)
        
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
        
        
        Defined at coordinates_module.fpp lines 212-285
        
        Parameters
        ----------
        reg : Region
        lim : float array
        
        """
        _amr2.f90wrap_limits(reg=self._handle, lim=lim)
    
    @staticmethod
    def checkifinside(pos, reg, ok, distance):
        """
        checkifinside(pos, reg, ok, distance)
        
        
        Defined at coordinates_module.fpp lines 297-317
        
        Parameters
        ----------
        pos : float array
        reg : Region
        ok : bool
        distance : float
        
        """
        _amr2.f90wrap_checkifinside(pos=pos, reg=reg._handle, ok=ok, distance=distance)
    
    @staticmethod
    def cube(self, reg, ok, distance):
        """
        cube(self, reg, ok, distance)
        
        
        Defined at coordinates_module.fpp lines 319-332
        
        Parameters
        ----------
        p : Vector
        reg : Region
        ok : bool
        distance : float
        
        """
        _amr2.f90wrap_cube(p=self._handle, reg=reg._handle, ok=ok, distance=distance)
    
    @staticmethod
    def sphere(self, reg, ok, distance):
        """
        sphere(self, reg, ok, distance)
        
        
        Defined at coordinates_module.fpp lines 334-341
        
        Parameters
        ----------
        p : Vector
        reg : Region
        ok : bool
        distance : float
        
        """
        _amr2.f90wrap_sphere(p=self._handle, reg=reg._handle, ok=ok, distance=distance)
    
    @staticmethod
    def cylinder(self, reg, ok, distance):
        """
        cylinder(self, reg, ok, distance)
        
        
        Defined at coordinates_module.fpp lines 343-351
        
        Parameters
        ----------
        p : Vector
        reg : Region
        ok : bool
        distance : float
        
        """
        _amr2.f90wrap_cylinder(p=self._handle, reg=reg._handle, ok=ok, \
            distance=distance)
    
    @staticmethod
    def cone(self, reg, ok, distance):
        """
        cone(self, reg, ok, distance)
        
        
        Defined at coordinates_module.fpp lines 353-363
        
        Parameters
        ----------
        p : Vector
        reg : Region
        ok : bool
        distance : float
        
        """
        _amr2.f90wrap_cone(p=self._handle, reg=reg._handle, ok=ok, distance=distance)
    
    _dt_array_initialisers = []
    

geometrical_regions = Geometrical_Regions()

class Io_Ramses(f90wrap.runtime.FortranModule):
    """
    Module io_ramses
    
    
    Defined at read_amr_module.fpp lines 24-1343
    
    """
    @f90wrap.runtime.register_class("amr2.hydroID")
    class hydroID(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=hydroid)
        
        
        Defined at read_amr_module.fpp lines 28-33
        
        """
        def __init__(self, handle=None):
            """
            self = Hydroid()
            
            
            Defined at read_amr_module.fpp lines 28-33
            
            
            Returns
            -------
            this : Hydroid
            	Object to be constructed
            
            
            Automatically generated constructor for hydroid
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _amr2.f90wrap_hydroid_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Hydroid
            
            
            Defined at read_amr_module.fpp lines 28-33
            
            Parameters
            ----------
            this : Hydroid
            	Object to be destructed
            
            
            Automatically generated destructor for hydroid
            """
            if self._alloc:
                _amr2.f90wrap_hydroid_finalise(this=self._handle)
        
        @property
        def nvar(self):
            """
            Element nvar ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 29
            
            """
            return _amr2.f90wrap_hydroid__get__nvar(self._handle)
        
        @nvar.setter
        def nvar(self, nvar):
            _amr2.f90wrap_hydroid__set__nvar(self._handle, nvar)
        
        @property
        def density(self):
            """
            Element density ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 30
            
            """
            return _amr2.f90wrap_hydroid__get__density(self._handle)
        
        @density.setter
        def density(self, density):
            _amr2.f90wrap_hydroid__set__density(self._handle, density)
        
        @property
        def vx(self):
            """
            Element vx ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 30
            
            """
            return _amr2.f90wrap_hydroid__get__vx(self._handle)
        
        @vx.setter
        def vx(self, vx):
            _amr2.f90wrap_hydroid__set__vx(self._handle, vx)
        
        @property
        def vy(self):
            """
            Element vy ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 30
            
            """
            return _amr2.f90wrap_hydroid__get__vy(self._handle)
        
        @vy.setter
        def vy(self, vy):
            _amr2.f90wrap_hydroid__set__vy(self._handle, vy)
        
        @property
        def vz(self):
            """
            Element vz ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 30
            
            """
            return _amr2.f90wrap_hydroid__get__vz(self._handle)
        
        @vz.setter
        def vz(self, vz):
            _amr2.f90wrap_hydroid__set__vz(self._handle, vz)
        
        @property
        def thermal_pressure(self):
            """
            Element thermal_pressure ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 30
            
            """
            return _amr2.f90wrap_hydroid__get__thermal_pressure(self._handle)
        
        @thermal_pressure.setter
        def thermal_pressure(self, thermal_pressure):
            _amr2.f90wrap_hydroid__set__thermal_pressure(self._handle, thermal_pressure)
        
        @property
        def metallicity(self):
            """
            Element metallicity ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 30
            
            """
            return _amr2.f90wrap_hydroid__get__metallicity(self._handle)
        
        @metallicity.setter
        def metallicity(self, metallicity):
            _amr2.f90wrap_hydroid__set__metallicity(self._handle, metallicity)
        
        @property
        def blx(self):
            """
            Element blx ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 31
            
            """
            return _amr2.f90wrap_hydroid__get__blx(self._handle)
        
        @blx.setter
        def blx(self, blx):
            _amr2.f90wrap_hydroid__set__blx(self._handle, blx)
        
        @property
        def bly(self):
            """
            Element bly ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 31
            
            """
            return _amr2.f90wrap_hydroid__get__bly(self._handle)
        
        @bly.setter
        def bly(self, bly):
            _amr2.f90wrap_hydroid__set__bly(self._handle, bly)
        
        @property
        def blz(self):
            """
            Element blz ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 31
            
            """
            return _amr2.f90wrap_hydroid__get__blz(self._handle)
        
        @blz.setter
        def blz(self, blz):
            _amr2.f90wrap_hydroid__set__blz(self._handle, blz)
        
        @property
        def brx(self):
            """
            Element brx ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 31
            
            """
            return _amr2.f90wrap_hydroid__get__brx(self._handle)
        
        @brx.setter
        def brx(self, brx):
            _amr2.f90wrap_hydroid__set__brx(self._handle, brx)
        
        @property
        def bry(self):
            """
            Element bry ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 31
            
            """
            return _amr2.f90wrap_hydroid__get__bry(self._handle)
        
        @bry.setter
        def bry(self, bry):
            _amr2.f90wrap_hydroid__set__bry(self._handle, bry)
        
        @property
        def brz(self):
            """
            Element brz ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 31
            
            """
            return _amr2.f90wrap_hydroid__get__brz(self._handle)
        
        @brz.setter
        def brz(self, brz):
            _amr2.f90wrap_hydroid__set__brz(self._handle, brz)
        
        @property
        def cr_pressure(self):
            """
            Element cr_pressure ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 32
            
            """
            return _amr2.f90wrap_hydroid__get__cr_pressure(self._handle)
        
        @cr_pressure.setter
        def cr_pressure(self, cr_pressure):
            _amr2.f90wrap_hydroid__set__cr_pressure(self._handle, cr_pressure)
        
        @property
        def xhii(self):
            """
            Element xhii ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 33
            
            """
            return _amr2.f90wrap_hydroid__get__xhii(self._handle)
        
        @xhii.setter
        def xhii(self, xhii):
            _amr2.f90wrap_hydroid__set__xhii(self._handle, xhii)
        
        @property
        def xheii(self):
            """
            Element xheii ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 33
            
            """
            return _amr2.f90wrap_hydroid__get__xheii(self._handle)
        
        @xheii.setter
        def xheii(self, xheii):
            _amr2.f90wrap_hydroid__set__xheii(self._handle, xheii)
        
        @property
        def xheiii(self):
            """
            Element xheiii ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 33
            
            """
            return _amr2.f90wrap_hydroid__get__xheiii(self._handle)
        
        @xheiii.setter
        def xheiii(self, xheiii):
            _amr2.f90wrap_hydroid__set__xheiii(self._handle, xheiii)
        
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
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("amr2.amr_info")
    class amr_info(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=amr_info)
        
        
        Defined at read_amr_module.fpp lines 35-43
        
        """
        def __init__(self, handle=None):
            """
            self = Amr_Info()
            
            
            Defined at read_amr_module.fpp lines 35-43
            
            
            Returns
            -------
            this : Amr_Info
            	Object to be constructed
            
            
            Automatically generated constructor for amr_info
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _amr2.f90wrap_amr_info_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Amr_Info
            
            
            Defined at read_amr_module.fpp lines 35-43
            
            Parameters
            ----------
            this : Amr_Info
            	Object to be destructed
            
            
            Automatically generated destructor for amr_info
            """
            if self._alloc:
                _amr2.f90wrap_amr_info_finalise(this=self._handle)
        
        @property
        def ncpu(self):
            """
            Element ncpu ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 36
            
            """
            return _amr2.f90wrap_amr_info__get__ncpu(self._handle)
        
        @ncpu.setter
        def ncpu(self, ncpu):
            _amr2.f90wrap_amr_info__set__ncpu(self._handle, ncpu)
        
        @property
        def ndim(self):
            """
            Element ndim ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 36
            
            """
            return _amr2.f90wrap_amr_info__get__ndim(self._handle)
        
        @ndim.setter
        def ndim(self, ndim):
            _amr2.f90wrap_amr_info__set__ndim(self._handle, ndim)
        
        @property
        def nlevelmax(self):
            """
            Element nlevelmax ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 36
            
            """
            return _amr2.f90wrap_amr_info__get__nlevelmax(self._handle)
        
        @nlevelmax.setter
        def nlevelmax(self, nlevelmax):
            _amr2.f90wrap_amr_info__set__nlevelmax(self._handle, nlevelmax)
        
        @property
        def nboundary(self):
            """
            Element nboundary ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 36
            
            """
            return _amr2.f90wrap_amr_info__get__nboundary(self._handle)
        
        @nboundary.setter
        def nboundary(self, nboundary):
            _amr2.f90wrap_amr_info__set__nboundary(self._handle, nboundary)
        
        @property
        def twotondim(self):
            """
            Element twotondim ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 36
            
            """
            return _amr2.f90wrap_amr_info__get__twotondim(self._handle)
        
        @twotondim.setter
        def twotondim(self, twotondim):
            _amr2.f90wrap_amr_info__set__twotondim(self._handle, twotondim)
        
        @property
        def ndom(self):
            """
            Element ndom ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 36
            
            """
            return _amr2.f90wrap_amr_info__get__ndom(self._handle)
        
        @ndom.setter
        def ndom(self, ndom):
            _amr2.f90wrap_amr_info__set__ndom(self._handle, ndom)
        
        @property
        def levelmin(self):
            """
            Element levelmin ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 37
            
            """
            return _amr2.f90wrap_amr_info__get__levelmin(self._handle)
        
        @levelmin.setter
        def levelmin(self, levelmin):
            _amr2.f90wrap_amr_info__set__levelmin(self._handle, levelmin)
        
        @property
        def levelmax(self):
            """
            Element levelmax ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 37
            
            """
            return _amr2.f90wrap_amr_info__get__levelmax(self._handle)
        
        @levelmax.setter
        def levelmax(self, levelmax):
            _amr2.f90wrap_amr_info__set__levelmax(self._handle, levelmax)
        
        @property
        def lmax(self):
            """
            Element lmax ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 37
            
            """
            return _amr2.f90wrap_amr_info__get__lmax(self._handle)
        
        @lmax.setter
        def lmax(self, lmax):
            _amr2.f90wrap_amr_info__set__lmax(self._handle, lmax)
        
        @property
        def lmin(self):
            """
            Element lmin ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 37
            
            """
            return _amr2.f90wrap_amr_info__get__lmin(self._handle)
        
        @lmin.setter
        def lmin(self, lmin):
            _amr2.f90wrap_amr_info__set__lmin(self._handle, lmin)
        
        @property
        def ncpu_read(self):
            """
            Element ncpu_read ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 38
            
            """
            return _amr2.f90wrap_amr_info__get__ncpu_read(self._handle)
        
        @ncpu_read.setter
        def ncpu_read(self, ncpu_read):
            _amr2.f90wrap_amr_info__set__ncpu_read(self._handle, ncpu_read)
        
        @property
        def ordering(self):
            """
            Element ordering ftype=character(80) pytype=str
            
            
            Defined at read_amr_module.fpp line 39
            
            """
            return _amr2.f90wrap_amr_info__get__ordering(self._handle)
        
        @ordering.setter
        def ordering(self, ordering):
            _amr2.f90wrap_amr_info__set__ordering(self._handle, ordering)
        
        @property
        def cpu_list(self):
            """
            Element cpu_list ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 40
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_amr_info__array__cpu_list(self._handle)
            if array_handle in self._arrays:
                cpu_list = self._arrays[array_handle]
            else:
                cpu_list = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_amr_info__array__cpu_list)
                self._arrays[array_handle] = cpu_list
            return cpu_list
        
        @cpu_list.setter
        def cpu_list(self, cpu_list):
            self.cpu_list[...] = cpu_list
        
        @property
        def bound_key(self):
            """
            Element bound_key ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 41
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_amr_info__array__bound_key(self._handle)
            if array_handle in self._arrays:
                bound_key = self._arrays[array_handle]
            else:
                bound_key = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_amr_info__array__bound_key)
                self._arrays[array_handle] = bound_key
            return bound_key
        
        @bound_key.setter
        def bound_key(self, bound_key):
            self.bound_key[...] = bound_key
        
        @property
        def cpu_read(self):
            """
            Element cpu_read ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 42
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_amr_info__array__cpu_read(self._handle)
            if array_handle in self._arrays:
                cpu_read = self._arrays[array_handle]
            else:
                cpu_read = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_amr_info__array__cpu_read)
                self._arrays[array_handle] = cpu_read
            return cpu_read
        
        @cpu_read.setter
        def cpu_read(self, cpu_read):
            self.cpu_read[...] = cpu_read
        
        @property
        def xbound(self):
            """
            Element xbound ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 43
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_amr_info__array__xbound(self._handle)
            if array_handle in self._arrays:
                xbound = self._arrays[array_handle]
            else:
                xbound = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_amr_info__array__xbound)
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
        
    
    @f90wrap.runtime.register_class("amr2.sim_info")
    class sim_info(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=sim_info)
        
        
        Defined at read_amr_module.fpp lines 45-51
        
        """
        def __init__(self, handle=None):
            """
            self = Sim_Info()
            
            
            Defined at read_amr_module.fpp lines 45-51
            
            
            Returns
            -------
            this : Sim_Info
            	Object to be constructed
            
            
            Automatically generated constructor for sim_info
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _amr2.f90wrap_sim_info_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Sim_Info
            
            
            Defined at read_amr_module.fpp lines 45-51
            
            Parameters
            ----------
            this : Sim_Info
            	Object to be destructed
            
            
            Automatically generated destructor for sim_info
            """
            if self._alloc:
                _amr2.f90wrap_sim_info_finalise(this=self._handle)
        
        @property
        def cosmo(self):
            """
            Element cosmo ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _amr2.f90wrap_sim_info__get__cosmo(self._handle)
        
        @cosmo.setter
        def cosmo(self, cosmo):
            _amr2.f90wrap_sim_info__set__cosmo(self._handle, cosmo)
        
        @property
        def family(self):
            """
            Element family ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 46
            
            """
            return _amr2.f90wrap_sim_info__get__family(self._handle)
        
        @family.setter
        def family(self, family):
            _amr2.f90wrap_sim_info__set__family(self._handle, family)
        
        @property
        def dm(self):
            """
            Element dm ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 47
            
            """
            return _amr2.f90wrap_sim_info__get__dm(self._handle)
        
        @dm.setter
        def dm(self, dm):
            _amr2.f90wrap_sim_info__set__dm(self._handle, dm)
        
        @property
        def hydro(self):
            """
            Element hydro ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 47
            
            """
            return _amr2.f90wrap_sim_info__get__hydro(self._handle)
        
        @hydro.setter
        def hydro(self, hydro):
            _amr2.f90wrap_sim_info__set__hydro(self._handle, hydro)
        
        @property
        def mhd(self):
            """
            Element mhd ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 47
            
            """
            return _amr2.f90wrap_sim_info__get__mhd(self._handle)
        
        @mhd.setter
        def mhd(self, mhd):
            _amr2.f90wrap_sim_info__set__mhd(self._handle, mhd)
        
        @property
        def cr(self):
            """
            Element cr ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 47
            
            """
            return _amr2.f90wrap_sim_info__get__cr(self._handle)
        
        @cr.setter
        def cr(self, cr):
            _amr2.f90wrap_sim_info__set__cr(self._handle, cr)
        
        @property
        def rt(self):
            """
            Element rt ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 47
            
            """
            return _amr2.f90wrap_sim_info__get__rt(self._handle)
        
        @rt.setter
        def rt(self, rt):
            _amr2.f90wrap_sim_info__set__rt(self._handle, rt)
        
        @property
        def bh(self):
            """
            Element bh ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 47
            
            """
            return _amr2.f90wrap_sim_info__get__bh(self._handle)
        
        @bh.setter
        def bh(self, bh):
            _amr2.f90wrap_sim_info__set__bh(self._handle, bh)
        
        @property
        def h0(self):
            """
            Element h0 ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 48
            
            """
            return _amr2.f90wrap_sim_info__get__h0(self._handle)
        
        @h0.setter
        def h0(self, h0):
            _amr2.f90wrap_sim_info__set__h0(self._handle, h0)
        
        @property
        def t(self):
            """
            Element t ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 48
            
            """
            return _amr2.f90wrap_sim_info__get__t(self._handle)
        
        @t.setter
        def t(self, t):
            _amr2.f90wrap_sim_info__set__t(self._handle, t)
        
        @property
        def aexp(self):
            """
            Element aexp ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 48
            
            """
            return _amr2.f90wrap_sim_info__get__aexp(self._handle)
        
        @aexp.setter
        def aexp(self, aexp):
            _amr2.f90wrap_sim_info__set__aexp(self._handle, aexp)
        
        @property
        def unit_l(self):
            """
            Element unit_l ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 48
            
            """
            return _amr2.f90wrap_sim_info__get__unit_l(self._handle)
        
        @unit_l.setter
        def unit_l(self, unit_l):
            _amr2.f90wrap_sim_info__set__unit_l(self._handle, unit_l)
        
        @property
        def unit_d(self):
            """
            Element unit_d ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 48
            
            """
            return _amr2.f90wrap_sim_info__get__unit_d(self._handle)
        
        @unit_d.setter
        def unit_d(self, unit_d):
            _amr2.f90wrap_sim_info__set__unit_d(self._handle, unit_d)
        
        @property
        def unit_t(self):
            """
            Element unit_t ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 48
            
            """
            return _amr2.f90wrap_sim_info__get__unit_t(self._handle)
        
        @unit_t.setter
        def unit_t(self, unit_t):
            _amr2.f90wrap_sim_info__set__unit_t(self._handle, unit_t)
        
        @property
        def unit_m(self):
            """
            Element unit_m ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 48
            
            """
            return _amr2.f90wrap_sim_info__get__unit_m(self._handle)
        
        @unit_m.setter
        def unit_m(self, unit_m):
            _amr2.f90wrap_sim_info__set__unit_m(self._handle, unit_m)
        
        @property
        def boxlen(self):
            """
            Element boxlen ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 48
            
            """
            return _amr2.f90wrap_sim_info__get__boxlen(self._handle)
        
        @boxlen.setter
        def boxlen(self, boxlen):
            _amr2.f90wrap_sim_info__set__boxlen(self._handle, boxlen)
        
        @property
        def omega_m(self):
            """
            Element omega_m ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 48
            
            """
            return _amr2.f90wrap_sim_info__get__omega_m(self._handle)
        
        @omega_m.setter
        def omega_m(self, omega_m):
            _amr2.f90wrap_sim_info__set__omega_m(self._handle, omega_m)
        
        @property
        def omega_l(self):
            """
            Element omega_l ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 48
            
            """
            return _amr2.f90wrap_sim_info__get__omega_l(self._handle)
        
        @omega_l.setter
        def omega_l(self, omega_l):
            _amr2.f90wrap_sim_info__set__omega_l(self._handle, omega_l)
        
        @property
        def omega_k(self):
            """
            Element omega_k ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 48
            
            """
            return _amr2.f90wrap_sim_info__get__omega_k(self._handle)
        
        @omega_k.setter
        def omega_k(self, omega_k):
            _amr2.f90wrap_sim_info__set__omega_k(self._handle, omega_k)
        
        @property
        def omega_b(self):
            """
            Element omega_b ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 48
            
            """
            return _amr2.f90wrap_sim_info__get__omega_b(self._handle)
        
        @omega_b.setter
        def omega_b(self, omega_b):
            _amr2.f90wrap_sim_info__set__omega_b(self._handle, omega_b)
        
        @property
        def time_tot(self):
            """
            Element time_tot ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 49
            
            """
            return _amr2.f90wrap_sim_info__get__time_tot(self._handle)
        
        @time_tot.setter
        def time_tot(self, time_tot):
            _amr2.f90wrap_sim_info__set__time_tot(self._handle, time_tot)
        
        @property
        def time_simu(self):
            """
            Element time_simu ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 49
            
            """
            return _amr2.f90wrap_sim_info__get__time_simu(self._handle)
        
        @time_simu.setter
        def time_simu(self, time_simu):
            _amr2.f90wrap_sim_info__set__time_simu(self._handle, time_simu)
        
        @property
        def n_frw(self):
            """
            Element n_frw ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 50
            
            """
            return _amr2.f90wrap_sim_info__get__n_frw(self._handle)
        
        @n_frw.setter
        def n_frw(self, n_frw):
            _amr2.f90wrap_sim_info__set__n_frw(self._handle, n_frw)
        
        @property
        def aexp_frw(self):
            """
            Element aexp_frw ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 51
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_sim_info__array__aexp_frw(self._handle)
            if array_handle in self._arrays:
                aexp_frw = self._arrays[array_handle]
            else:
                aexp_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_sim_info__array__aexp_frw)
                self._arrays[array_handle] = aexp_frw
            return aexp_frw
        
        @aexp_frw.setter
        def aexp_frw(self, aexp_frw):
            self.aexp_frw[...] = aexp_frw
        
        @property
        def hexp_frw(self):
            """
            Element hexp_frw ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 51
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_sim_info__array__hexp_frw(self._handle)
            if array_handle in self._arrays:
                hexp_frw = self._arrays[array_handle]
            else:
                hexp_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_sim_info__array__hexp_frw)
                self._arrays[array_handle] = hexp_frw
            return hexp_frw
        
        @hexp_frw.setter
        def hexp_frw(self, hexp_frw):
            self.hexp_frw[...] = hexp_frw
        
        @property
        def tau_frw(self):
            """
            Element tau_frw ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 51
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_sim_info__array__tau_frw(self._handle)
            if array_handle in self._arrays:
                tau_frw = self._arrays[array_handle]
            else:
                tau_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_sim_info__array__tau_frw)
                self._arrays[array_handle] = tau_frw
            return tau_frw
        
        @tau_frw.setter
        def tau_frw(self, tau_frw):
            self.tau_frw[...] = tau_frw
        
        @property
        def t_frw(self):
            """
            Element t_frw ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 51
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_sim_info__array__t_frw(self._handle)
            if array_handle in self._arrays:
                t_frw = self._arrays[array_handle]
            else:
                t_frw = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_sim_info__array__t_frw)
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
        
    
    @f90wrap.runtime.register_class("amr2.level")
    class level(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=level)
        
        
        Defined at read_amr_module.fpp lines 53-64
        
        """
        def __init__(self, handle=None):
            """
            self = Level()
            
            
            Defined at read_amr_module.fpp lines 53-64
            
            
            Returns
            -------
            this : Level
            	Object to be constructed
            
            
            Automatically generated constructor for level
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _amr2.f90wrap_level_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Level
            
            
            Defined at read_amr_module.fpp lines 53-64
            
            Parameters
            ----------
            this : Level
            	Object to be destructed
            
            
            Automatically generated destructor for level
            """
            if self._alloc:
                _amr2.f90wrap_level_finalise(this=self._handle)
        
        @property
        def ilevel(self):
            """
            Element ilevel ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 54
            
            """
            return _amr2.f90wrap_level__get__ilevel(self._handle)
        
        @ilevel.setter
        def ilevel(self, ilevel):
            _amr2.f90wrap_level__set__ilevel(self._handle, ilevel)
        
        @property
        def ngrid(self):
            """
            Element ngrid ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 55
            
            """
            return _amr2.f90wrap_level__get__ngrid(self._handle)
        
        @ngrid.setter
        def ngrid(self, ngrid):
            _amr2.f90wrap_level__set__ngrid(self._handle, ngrid)
        
        @property
        def cube(self):
            """
            Element cube ftype=real(sgl) pytype=float
            
            
            Defined at read_amr_module.fpp line 56
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_level__array__cube(self._handle)
            if array_handle in self._arrays:
                cube = self._arrays[array_handle]
            else:
                cube = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_level__array__cube)
                self._arrays[array_handle] = cube
            return cube
        
        @cube.setter
        def cube(self, cube):
            self.cube[...] = cube
        
        @property
        def map(self):
            """
            Element map ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 57
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_level__array__map(self._handle)
            if array_handle in self._arrays:
                map = self._arrays[array_handle]
            else:
                map = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_level__array__map)
                self._arrays[array_handle] = map
            return map
        
        @map.setter
        def map(self, map):
            self.map[...] = map
        
        @property
        def rho(self):
            """
            Element rho ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 58
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_level__array__rho(self._handle)
            if array_handle in self._arrays:
                rho = self._arrays[array_handle]
            else:
                rho = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_level__array__rho)
                self._arrays[array_handle] = rho
            return rho
        
        @rho.setter
        def rho(self, rho):
            self.rho[...] = rho
        
        @property
        def imin(self):
            """
            Element imin ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 59
            
            """
            return _amr2.f90wrap_level__get__imin(self._handle)
        
        @imin.setter
        def imin(self, imin):
            _amr2.f90wrap_level__set__imin(self._handle, imin)
        
        @property
        def imax(self):
            """
            Element imax ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 60
            
            """
            return _amr2.f90wrap_level__get__imax(self._handle)
        
        @imax.setter
        def imax(self, imax):
            _amr2.f90wrap_level__set__imax(self._handle, imax)
        
        @property
        def jmin(self):
            """
            Element jmin ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 61
            
            """
            return _amr2.f90wrap_level__get__jmin(self._handle)
        
        @jmin.setter
        def jmin(self, jmin):
            _amr2.f90wrap_level__set__jmin(self._handle, jmin)
        
        @property
        def jmax(self):
            """
            Element jmax ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 62
            
            """
            return _amr2.f90wrap_level__get__jmax(self._handle)
        
        @jmax.setter
        def jmax(self, jmax):
            _amr2.f90wrap_level__set__jmax(self._handle, jmax)
        
        @property
        def kmin(self):
            """
            Element kmin ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 63
            
            """
            return _amr2.f90wrap_level__get__kmin(self._handle)
        
        @kmin.setter
        def kmin(self, kmin):
            _amr2.f90wrap_level__set__kmin(self._handle, kmin)
        
        @property
        def kmax(self):
            """
            Element kmax ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 64
            
            """
            return _amr2.f90wrap_level__get__kmax(self._handle)
        
        @kmax.setter
        def kmax(self, kmax):
            _amr2.f90wrap_level__set__kmax(self._handle, kmax)
        
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
        
    
    @f90wrap.runtime.register_class("amr2.data_handler")
    class data_handler(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=data_handler)
        
        
        Defined at read_amr_module.fpp lines 66-70
        
        """
        def __init__(self, handle=None):
            """
            self = Data_Handler()
            
            
            Defined at read_amr_module.fpp lines 66-70
            
            
            Returns
            -------
            this : Data_Handler
            	Object to be constructed
            
            
            Automatically generated constructor for data_handler
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _amr2.f90wrap_data_handler_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Data_Handler
            
            
            Defined at read_amr_module.fpp lines 66-70
            
            Parameters
            ----------
            this : Data_Handler
            	Object to be destructed
            
            
            Automatically generated destructor for data_handler
            """
            if self._alloc:
                _amr2.f90wrap_data_handler_finalise(this=self._handle)
        
        @property
        def name(self):
            """
            Element name ftype=character(80) pytype=str
            
            
            Defined at read_amr_module.fpp line 67
            
            """
            return _amr2.f90wrap_data_handler__get__name(self._handle)
        
        @name.setter
        def name(self, name):
            _amr2.f90wrap_data_handler__set__name(self._handle, name)
        
        @property
        def x_data(self):
            """
            Element x_data ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 68
            
            """
            return _amr2.f90wrap_data_handler__get__x_data(self._handle)
        
        @x_data.setter
        def x_data(self, x_data):
            _amr2.f90wrap_data_handler__set__x_data(self._handle, x_data)
        
        @property
        def y_data(self):
            """
            Element y_data ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 68
            
            """
            return _amr2.f90wrap_data_handler__get__y_data(self._handle)
        
        @y_data.setter
        def y_data(self, y_data):
            _amr2.f90wrap_data_handler__set__y_data(self._handle, y_data)
        
        @property
        def z_data(self):
            """
            Element z_data ftype=logical pytype=bool
            
            
            Defined at read_amr_module.fpp line 68
            
            """
            return _amr2.f90wrap_data_handler__get__z_data(self._handle)
        
        @z_data.setter
        def z_data(self, z_data):
            _amr2.f90wrap_data_handler__set__z_data(self._handle, z_data)
        
        @property
        def nx(self):
            """
            Element nx ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 69
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_data_handler__array__nx(self._handle)
            if array_handle in self._arrays:
                nx = self._arrays[array_handle]
            else:
                nx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_data_handler__array__nx)
                self._arrays[array_handle] = nx
            return nx
        
        @nx.setter
        def nx(self, nx):
            self.nx[...] = nx
        
        @property
        def ny(self):
            """
            Element ny ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 69
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_data_handler__array__ny(self._handle)
            if array_handle in self._arrays:
                ny = self._arrays[array_handle]
            else:
                ny = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_data_handler__array__ny)
                self._arrays[array_handle] = ny
            return ny
        
        @ny.setter
        def ny(self, ny):
            self.ny[...] = ny
        
        @property
        def nz(self):
            """
            Element nz ftype=integer pytype=int
            
            
            Defined at read_amr_module.fpp line 69
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_data_handler__array__nz(self._handle)
            if array_handle in self._arrays:
                nz = self._arrays[array_handle]
            else:
                nz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_data_handler__array__nz)
                self._arrays[array_handle] = nz
            return nz
        
        @nz.setter
        def nz(self, nz):
            self.nz[...] = nz
        
        @property
        def x(self):
            """
            Element x ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 70
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_data_handler__array__x(self._handle)
            if array_handle in self._arrays:
                x = self._arrays[array_handle]
            else:
                x = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_data_handler__array__x)
                self._arrays[array_handle] = x
            return x
        
        @x.setter
        def x(self, x):
            self.x[...] = x
        
        @property
        def y(self):
            """
            Element y ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 70
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_data_handler__array__y(self._handle)
            if array_handle in self._arrays:
                y = self._arrays[array_handle]
            else:
                y = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_data_handler__array__y)
                self._arrays[array_handle] = y
            return y
        
        @y.setter
        def y(self, y):
            self.y[...] = y
        
        @property
        def z(self):
            """
            Element z ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 70
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_data_handler__array__z(self._handle)
            if array_handle in self._arrays:
                z = self._arrays[array_handle]
            else:
                z = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_data_handler__array__z)
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
        
    
    @f90wrap.runtime.register_class("amr2.particle")
    class particle(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=particle)
        
        
        Defined at read_amr_module.fpp lines 72-75
        
        """
        def __init__(self, handle=None):
            """
            self = Particle()
            
            
            Defined at read_amr_module.fpp lines 72-75
            
            
            Returns
            -------
            this : Particle
            	Object to be constructed
            
            
            Automatically generated constructor for particle
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _amr2.f90wrap_particle_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Particle
            
            
            Defined at read_amr_module.fpp lines 72-75
            
            Parameters
            ----------
            this : Particle
            	Object to be destructed
            
            
            Automatically generated destructor for particle
            """
            if self._alloc:
                _amr2.f90wrap_particle_finalise(this=self._handle)
        
        @property
        def id(self):
            """
            Element id ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 73
            
            """
            return _amr2.f90wrap_particle__get__id(self._handle)
        
        @id.setter
        def id(self, id):
            _amr2.f90wrap_particle__set__id(self._handle, id)
        
        @property
        def x(self):
            """
            Element x ftype=type(vector) pytype=Vector
            
            
            Defined at read_amr_module.fpp line 74
            
            """
            x_handle = _amr2.f90wrap_particle__get__x(self._handle)
            if tuple(x_handle) in self._objs:
                x = self._objs[tuple(x_handle)]
            else:
                x = vectors.vector.from_handle(x_handle)
                self._objs[tuple(x_handle)] = x
            return x
        
        @x.setter
        def x(self, x):
            x = x._handle
            _amr2.f90wrap_particle__set__x(self._handle, x)
        
        @property
        def v(self):
            """
            Element v ftype=type(vector) pytype=Vector
            
            
            Defined at read_amr_module.fpp line 74
            
            """
            v_handle = _amr2.f90wrap_particle__get__v(self._handle)
            if tuple(v_handle) in self._objs:
                v = self._objs[tuple(v_handle)]
            else:
                v = vectors.vector.from_handle(v_handle)
                self._objs[tuple(v_handle)] = v
            return v
        
        @v.setter
        def v(self, v):
            v = v._handle
            _amr2.f90wrap_particle__set__v(self._handle, v)
        
        @property
        def m(self):
            """
            Element m ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 75
            
            """
            return _amr2.f90wrap_particle__get__m(self._handle)
        
        @m.setter
        def m(self, m):
            _amr2.f90wrap_particle__set__m(self._handle, m)
        
        @property
        def met(self):
            """
            Element met ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 75
            
            """
            return _amr2.f90wrap_particle__get__met(self._handle)
        
        @met.setter
        def met(self, met):
            _amr2.f90wrap_particle__set__met(self._handle, met)
        
        @property
        def imass(self):
            """
            Element imass ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 75
            
            """
            return _amr2.f90wrap_particle__get__imass(self._handle)
        
        @imass.setter
        def imass(self, imass):
            _amr2.f90wrap_particle__set__imass(self._handle, imass)
        
        @property
        def age(self):
            """
            Element age ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 75
            
            """
            return _amr2.f90wrap_particle__get__age(self._handle)
        
        @age.setter
        def age(self, age):
            _amr2.f90wrap_particle__set__age(self._handle, age)
        
        @property
        def tform(self):
            """
            Element tform ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 75
            
            """
            return _amr2.f90wrap_particle__get__tform(self._handle)
        
        @tform.setter
        def tform(self, tform):
            _amr2.f90wrap_particle__set__tform(self._handle, tform)
        
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
        
        
        Defined at read_amr_module.fpp lines 89-113
        
        Parameters
        ----------
        n : int
        nchar : str
        
        """
        _amr2.f90wrap_title(n=n, nchar=nchar)
    
    @staticmethod
    def hilbert3d(x, y, z, order, bit_length, npoint):
        """
        hilbert3d(x, y, z, order, bit_length, npoint)
        
        
        Defined at read_amr_module.fpp lines 120-192
        
        Parameters
        ----------
        x : int array
        y : int array
        z : int array
        order : float array
        bit_length : int
        npoint : int
        
        """
        _amr2.f90wrap_hilbert3d(x=x, y=y, z=z, order=order, bit_length=bit_length, \
            npoint=npoint)
    
    @staticmethod
    def check_lmax(ngridfile):
        """
        check_lmax(ngridfile)
        
        
        Defined at read_amr_module.fpp lines 200-212
        
        Parameters
        ----------
        ngridfile : int array
        
        """
        _amr2.f90wrap_check_lmax(ngridfile=ngridfile)
    
    @staticmethod
    def check_families(repository):
        """
        check_families(repository)
        
        
        Defined at read_amr_module.fpp lines 221-246
        
        Parameters
        ----------
        repository : str
        
        """
        _amr2.f90wrap_check_families(repository=repository)
    
    @staticmethod
    def read_hydrofile_descriptor(repository):
        """
        read_hydrofile_descriptor(repository)
        
        
        Defined at read_amr_module.fpp lines 255-289
        
        Parameters
        ----------
        repository : str
        
        """
        _amr2.f90wrap_read_hydrofile_descriptor(repository=repository)
    
    @staticmethod
    def read_hydrofile_descriptor_old(repository):
        """
        read_hydrofile_descriptor_old(repository)
        
        
        Defined at read_amr_module.fpp lines 299-355
        
        Parameters
        ----------
        repository : str
        
        """
        _amr2.f90wrap_read_hydrofile_descriptor_old(repository=repository)
    
    @staticmethod
    def select_from_descriptor_ids(newvar, newid):
        """
        select_from_descriptor_ids(newvar, newid)
        
        
        Defined at read_amr_module.fpp lines 357-433
        
        Parameters
        ----------
        newvar : str
        newid : int
        
        """
        _amr2.f90wrap_select_from_descriptor_ids(newvar=newvar, newid=newid)
    
    @staticmethod
    def read_hydrofile_descriptor_new(repository):
        """
        read_hydrofile_descriptor_new(repository)
        
        
        Defined at read_amr_module.fpp lines 443-474
        
        Parameters
        ----------
        repository : str
        
        """
        _amr2.f90wrap_read_hydrofile_descriptor_new(repository=repository)
    
    @staticmethod
    def getvarvalue(self, dx, x, var, varname, value):
        """
        getvarvalue(self, dx, x, var, varname, value)
        
        
        Defined at read_amr_module.fpp lines 483-742
        
        Parameters
        ----------
        reg : Region
        dx : float
        x : Vector
        var : float array
        varname : str
        value : float
        
        """
        _amr2.f90wrap_getvarvalue(reg=self._handle, dx=dx, x=x._handle, var=var, \
            varname=varname, value=value)
    
    @staticmethod
    def init_amr_read(repository):
        """
        init_amr_read(repository)
        
        
        Defined at read_amr_module.fpp lines 751-839
        
        Parameters
        ----------
        repository : str
        
        """
        _amr2.f90wrap_init_amr_read(repository=repository)
    
    @staticmethod
    def get_cpu_map(self):
        """
        get_cpu_map(self)
        
        
        Defined at read_amr_module.fpp lines 847-935
        
        Parameters
        ----------
        reg : Region
        
        """
        _amr2.f90wrap_get_cpu_map(reg=self._handle)
    
    @staticmethod
    def getparttype(self, ptype):
        """
        getparttype(self, ptype)
        
        
        Defined at read_amr_module.fpp lines 937-945
        
        Parameters
        ----------
        part : Particle
        ptype : str
        
        """
        _amr2.f90wrap_getparttype(part=self._handle, ptype=ptype)
    
    @staticmethod
    def getpartvalue(self, part, var, value, dx=None):
        """
        getpartvalue(self, part, var, value[, dx])
        
        
        Defined at read_amr_module.fpp lines 947-1342
        
        Parameters
        ----------
        reg : Region
        part : Particle
        var : str
        value : float
        dx : Vector
        
        """
        _amr2.f90wrap_getpartvalue(reg=self._handle, part=part._handle, var=var, \
            value=value, dx=None if dx is None else dx._handle)
    
    @property
    def longint(self):
        """
        Element longint ftype=logical pytype=bool
        
        
        Defined at read_amr_module.fpp line 82
        
        """
        return _amr2.f90wrap_io_ramses__get__longint()
    
    @longint.setter
    def longint(self, longint):
        _amr2.f90wrap_io_ramses__set__longint(longint)
    
    def __str__(self):
        ret = ['<io_ramses>{\n']
        ret.append('    longint : ')
        ret.append(repr(self.longint))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

io_ramses = Io_Ramses()

class Filtering(f90wrap.runtime.FortranModule):
    """
    Module filtering
    
    
    Defined at read_amr_module.fpp lines 1345-1437
    
    """
    @f90wrap.runtime.register_class("amr2.filter")
    class filter(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=filter)
        
        
        Defined at read_amr_module.fpp lines 1348-1353
        
        """
        def __init__(self, handle=None):
            """
            self = Filter()
            
            
            Defined at read_amr_module.fpp lines 1348-1353
            
            
            Returns
            -------
            this : Filter
            	Object to be constructed
            
            
            Automatically generated constructor for filter
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _amr2.f90wrap_filter_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Filter
            
            
            Defined at read_amr_module.fpp lines 1348-1353
            
            Parameters
            ----------
            this : Filter
            	Object to be destructed
            
            
            Automatically generated destructor for filter
            """
            if self._alloc:
                _amr2.f90wrap_filter_finalise(this=self._handle)
        
        @property
        def name(self):
            """
            Element name ftype=character(128) pytype=str
            
            
            Defined at read_amr_module.fpp line 1349
            
            """
            return _amr2.f90wrap_filter__get__name(self._handle)
        
        @name.setter
        def name(self, name):
            _amr2.f90wrap_filter__set__name(self._handle, name)
        
        @property
        def ncond(self):
            """
            Element ncond ftype=integer  pytype=int
            
            
            Defined at read_amr_module.fpp line 1350
            
            """
            return _amr2.f90wrap_filter__get__ncond(self._handle)
        
        @ncond.setter
        def ncond(self, ncond):
            _amr2.f90wrap_filter__set__ncond(self._handle, ncond)
        
        @property
        def cond_vars(self):
            """
            Element cond_vars ftype=character(128) pytype=str
            
            
            Defined at read_amr_module.fpp line 1351
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_filter__array__cond_vars(self._handle)
            if array_handle in self._arrays:
                cond_vars = self._arrays[array_handle]
            else:
                cond_vars = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_filter__array__cond_vars)
                self._arrays[array_handle] = cond_vars
            return cond_vars
        
        @cond_vars.setter
        def cond_vars(self, cond_vars):
            self.cond_vars[...] = cond_vars
        
        @property
        def cond_ops(self):
            """
            Element cond_ops ftype=character(2) pytype=str
            
            
            Defined at read_amr_module.fpp line 1352
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_filter__array__cond_ops(self._handle)
            if array_handle in self._arrays:
                cond_ops = self._arrays[array_handle]
            else:
                cond_ops = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_filter__array__cond_ops)
                self._arrays[array_handle] = cond_ops
            return cond_ops
        
        @cond_ops.setter
        def cond_ops(self, cond_ops):
            self.cond_ops[...] = cond_ops
        
        @property
        def cond_vals(self):
            """
            Element cond_vals ftype=real(dbl) pytype=float
            
            
            Defined at read_amr_module.fpp line 1353
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_filter__array__cond_vals(self._handle)
            if array_handle in self._arrays:
                cond_vals = self._arrays[array_handle]
            else:
                cond_vals = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_filter__array__cond_vals)
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
        
        
        Defined at read_amr_module.fpp lines 1356-1361
        
        Parameters
        ----------
        filt : Filter
        
        """
        _amr2.f90wrap_allocate_filter(filt=self._handle)
    
    @staticmethod
    def cond_string_to_filter(str, filt):
        """
        cond_string_to_filter(str, filt)
        
        
        Defined at read_amr_module.fpp lines 1363-1367
        
        Parameters
        ----------
        str : str
        filt : Filter
        
        """
        _amr2.f90wrap_cond_string_to_filter(str=str, filt=filt._handle)
    
    @staticmethod
    def filter_cell(self, filt, cell_x, cell_dx, cell_var):
        """
        filter_cell = filter_cell(self, filt, cell_x, cell_dx, cell_var)
        
        
        Defined at read_amr_module.fpp lines 1369-1401
        
        Parameters
        ----------
        reg : Region
        filt : Filter
        cell_x : Vector
        cell_dx : float
        cell_var : float array
        
        Returns
        -------
        filter_cell : bool
        
        """
        filter_cell = _amr2.f90wrap_filter_cell(reg=self._handle, filt=filt._handle, \
            cell_x=cell_x._handle, cell_dx=cell_dx, cell_var=cell_var)
        return filter_cell
    
    @staticmethod
    def filter_particle(self, filt, part, dx=None):
        """
        filter_particle = filter_particle(self, filt, part[, dx])
        
        
        Defined at read_amr_module.fpp lines 1403-1437
        
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
        filter_particle = _amr2.f90wrap_filter_particle(reg=self._handle, \
            filt=filt._handle, part=part._handle, dx=None if dx is None else dx._handle)
        return filter_particle
    
    _dt_array_initialisers = []
    

filtering = Filtering()

class Amr_Profiles(f90wrap.runtime.FortranModule):
    """
    Module amr_profiles
    
    
    Defined at profiles_module.fpp lines 25-756
    
    """
    @f90wrap.runtime.register_class("amr2.profile_handler")
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
            result = _amr2.f90wrap_profile_handler_initialise()
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
                _amr2.f90wrap_profile_handler_finalise(this=self._handle)
        
        @property
        def profdim(self):
            """
            Element profdim ftype=integer  pytype=int
            
            
            Defined at profiles_module.fpp line 30
            
            """
            return _amr2.f90wrap_profile_handler__get__profdim(self._handle)
        
        @profdim.setter
        def profdim(self, profdim):
            _amr2.f90wrap_profile_handler__set__profdim(self._handle, profdim)
        
        @property
        def xvarname(self):
            """
            Element xvarname ftype=character(128) pytype=str
            
            
            Defined at profiles_module.fpp line 31
            
            """
            return _amr2.f90wrap_profile_handler__get__xvarname(self._handle)
        
        @xvarname.setter
        def xvarname(self, xvarname):
            _amr2.f90wrap_profile_handler__set__xvarname(self._handle, xvarname)
        
        @property
        def nyvar(self):
            """
            Element nyvar ftype=integer  pytype=int
            
            
            Defined at profiles_module.fpp line 32
            
            """
            return _amr2.f90wrap_profile_handler__get__nyvar(self._handle)
        
        @nyvar.setter
        def nyvar(self, nyvar):
            _amr2.f90wrap_profile_handler__set__nyvar(self._handle, nyvar)
        
        @property
        def yvarnames(self):
            """
            Element yvarnames ftype=character(128) pytype=str
            
            
            Defined at profiles_module.fpp line 33
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_profile_handler__array__yvarnames(self._handle)
            if array_handle in self._arrays:
                yvarnames = self._arrays[array_handle]
            else:
                yvarnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_profile_handler__array__yvarnames)
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
            return _amr2.f90wrap_profile_handler__get__nbins(self._handle)
        
        @nbins.setter
        def nbins(self, nbins):
            _amr2.f90wrap_profile_handler__set__nbins(self._handle, nbins)
        
        @property
        def nwvar(self):
            """
            Element nwvar ftype=integer  pytype=int
            
            
            Defined at profiles_module.fpp line 35
            
            """
            return _amr2.f90wrap_profile_handler__get__nwvar(self._handle)
        
        @nwvar.setter
        def nwvar(self, nwvar):
            _amr2.f90wrap_profile_handler__set__nwvar(self._handle, nwvar)
        
        @property
        def wvarnames(self):
            """
            Element wvarnames ftype=character(128) pytype=str
            
            
            Defined at profiles_module.fpp line 36
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_profile_handler__array__wvarnames(self._handle)
            if array_handle in self._arrays:
                wvarnames = self._arrays[array_handle]
            else:
                wvarnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_profile_handler__array__wvarnames)
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
                _amr2.f90wrap_profile_handler__array__xdata(self._handle)
            if array_handle in self._arrays:
                xdata = self._arrays[array_handle]
            else:
                xdata = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_profile_handler__array__xdata)
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
                _amr2.f90wrap_profile_handler__array__ydata(self._handle)
            if array_handle in self._arrays:
                ydata = self._arrays[array_handle]
            else:
                ydata = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_profile_handler__array__ydata)
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
        
    
    @f90wrap.runtime.register_class("amr2.profile_handler_twod")
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
            result = _amr2.f90wrap_profile_handler_twod_initialise()
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
                _amr2.f90wrap_profile_handler_twod_finalise(this=self._handle)
        
        @property
        def profdim(self):
            """
            Element profdim ftype=integer  pytype=int
            
            
            Defined at profiles_module.fpp line 41
            
            """
            return _amr2.f90wrap_profile_handler_twod__get__profdim(self._handle)
        
        @profdim.setter
        def profdim(self, profdim):
            _amr2.f90wrap_profile_handler_twod__set__profdim(self._handle, profdim)
        
        @property
        def xvarname(self):
            """
            Element xvarname ftype=character(128) pytype=str
            
            
            Defined at profiles_module.fpp line 42
            
            """
            return _amr2.f90wrap_profile_handler_twod__get__xvarname(self._handle)
        
        @xvarname.setter
        def xvarname(self, xvarname):
            _amr2.f90wrap_profile_handler_twod__set__xvarname(self._handle, xvarname)
        
        @property
        def yvarname(self):
            """
            Element yvarname ftype=character(128) pytype=str
            
            
            Defined at profiles_module.fpp line 43
            
            """
            return _amr2.f90wrap_profile_handler_twod__get__yvarname(self._handle)
        
        @yvarname.setter
        def yvarname(self, yvarname):
            _amr2.f90wrap_profile_handler_twod__set__yvarname(self._handle, yvarname)
        
        @property
        def nzvar(self):
            """
            Element nzvar ftype=integer  pytype=int
            
            
            Defined at profiles_module.fpp line 44
            
            """
            return _amr2.f90wrap_profile_handler_twod__get__nzvar(self._handle)
        
        @nzvar.setter
        def nzvar(self, nzvar):
            _amr2.f90wrap_profile_handler_twod__set__nzvar(self._handle, nzvar)
        
        @property
        def zvarnames(self):
            """
            Element zvarnames ftype=character(128) pytype=str
            
            
            Defined at profiles_module.fpp line 45
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_profile_handler_twod__array__zvarnames(self._handle)
            if array_handle in self._arrays:
                zvarnames = self._arrays[array_handle]
            else:
                zvarnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_profile_handler_twod__array__zvarnames)
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
                _amr2.f90wrap_profile_handler_twod__array__nbins(self._handle)
            if array_handle in self._arrays:
                nbins = self._arrays[array_handle]
            else:
                nbins = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_profile_handler_twod__array__nbins)
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
            return _amr2.f90wrap_profile_handler_twod__get__nwvar(self._handle)
        
        @nwvar.setter
        def nwvar(self, nwvar):
            _amr2.f90wrap_profile_handler_twod__set__nwvar(self._handle, nwvar)
        
        @property
        def wvarnames(self):
            """
            Element wvarnames ftype=character(128) pytype=str
            
            
            Defined at profiles_module.fpp line 48
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_profile_handler_twod__array__wvarnames(self._handle)
            if array_handle in self._arrays:
                wvarnames = self._arrays[array_handle]
            else:
                wvarnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_profile_handler_twod__array__wvarnames)
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
                _amr2.f90wrap_profile_handler_twod__array__xdata(self._handle)
            if array_handle in self._arrays:
                xdata = self._arrays[array_handle]
            else:
                xdata = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_profile_handler_twod__array__xdata)
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
                _amr2.f90wrap_profile_handler_twod__array__ydata(self._handle)
            if array_handle in self._arrays:
                ydata = self._arrays[array_handle]
            else:
                ydata = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_profile_handler_twod__array__ydata)
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
                _amr2.f90wrap_profile_handler_twod__array__zdata(self._handle)
            if array_handle in self._arrays:
                zdata = self._arrays[array_handle]
            else:
                zdata = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_profile_handler_twod__array__zdata)
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
        
    
    @staticmethod
    def allocate_profile_handler(self):
        """
        allocate_profile_handler(self)
        
        
        Defined at profiles_module.fpp lines 53-59
        
        Parameters
        ----------
        prof : Profile_Handler
        
        """
        _amr2.f90wrap_allocate_profile_handler(prof=self._handle)
    
    @staticmethod
    def allocate_profile_handler_twod(self):
        """
        allocate_profile_handler_twod(self)
        
        
        Defined at profiles_module.fpp lines 61-68
        
        Parameters
        ----------
        prof : Profile_Handler_Twod
        
        """
        _amr2.f90wrap_allocate_profile_handler_twod(prof=self._handle)
    
    @staticmethod
    def makebins(self, varname, nbins, bins, logscale):
        """
        makebins(self, varname, nbins, bins, logscale)
        
        
        Defined at profiles_module.fpp lines 70-108
        
        Parameters
        ----------
        reg : Region
        varname : str
        nbins : int
        bins : float array
        logscale : bool
        
        """
        _amr2.f90wrap_makebins(reg=self._handle, varname=varname, nbins=nbins, \
            bins=bins, logscale=logscale)
    
    @staticmethod
    def findbinpos(self, distance, pos, cellvars, cellsize, prof, ibin):
        """
        findbinpos(self, distance, pos, cellvars, cellsize, prof, ibin)
        
        
        Defined at profiles_module.fpp lines 110-131
        
        Parameters
        ----------
        reg : Region
        distance : float
        pos : float array
        cellvars : float array
        cellsize : float
        prof : Profile_Handler
        ibin : int
        
        """
        _amr2.f90wrap_findbinpos(reg=self._handle, distance=distance, pos=pos, \
            cellvars=cellvars, cellsize=cellsize, prof=prof._handle, ibin=ibin)
    
    @staticmethod
    def findbinpos_twod(self, distance, pos, cellvars, cellsize, prof, logscale, \
        ibinx, ibiny):
        """
        findbinpos_twod(self, distance, pos, cellvars, cellsize, prof, logscale, ibinx, \
            ibiny)
        
        
        Defined at profiles_module.fpp lines 133-171
        
        Parameters
        ----------
        reg : Region
        distance : float
        pos : float array
        cellvars : float array
        cellsize : float
        prof : Profile_Handler_Twod
        logscale : bool
        ibinx : int
        ibiny : int
        
        """
        _amr2.f90wrap_findbinpos_twod(reg=self._handle, distance=distance, pos=pos, \
            cellvars=cellvars, cellsize=cellsize, prof=prof._handle, logscale=logscale, \
            ibinx=ibinx, ibiny=ibiny)
    
    @staticmethod
    def bindata(self, pos, cellvars, cellsize, prof, ibin):
        """
        bindata(self, pos, cellvars, cellsize, prof, ibin)
        
        
        Defined at profiles_module.fpp lines 173-219
        
        Parameters
        ----------
        reg : Region
        pos : float array
        cellvars : float array
        cellsize : float
        prof : Profile_Handler
        ibin : int
        
        """
        _amr2.f90wrap_bindata(reg=self._handle, pos=pos, cellvars=cellvars, \
            cellsize=cellsize, prof=prof._handle, ibin=ibin)
    
    @staticmethod
    def bindata_twod(self, pos, cellvars, cellsize, prof, ibinx, ibiny):
        """
        bindata_twod(self, pos, cellvars, cellsize, prof, ibinx, ibiny)
        
        
        Defined at profiles_module.fpp lines 221-253
        
        Parameters
        ----------
        reg : Region
        pos : float array
        cellvars : float array
        cellsize : float
        prof : Profile_Handler_Twod
        ibinx : int
        ibiny : int
        
        """
        _amr2.f90wrap_bindata_twod(reg=self._handle, pos=pos, cellvars=cellvars, \
            cellsize=cellsize, prof=prof._handle, ibinx=ibinx, ibiny=ibiny)
    
    @staticmethod
    def renormalise_bins(self):
        """
        renormalise_bins(self)
        
        
        Defined at profiles_module.fpp lines 255-280
        
        Parameters
        ----------
        prof_data : Profile_Handler
        
        """
        _amr2.f90wrap_renormalise_bins(prof_data=self._handle)
    
    @staticmethod
    def renormalise_bins_twod(self):
        """
        renormalise_bins_twod(self)
        
        
        Defined at profiles_module.fpp lines 282-306
        
        Parameters
        ----------
        prof_data : Profile_Handler_Twod
        
        """
        _amr2.f90wrap_renormalise_bins_twod(prof_data=self._handle)
    
    @staticmethod
    def get_cells_onedprofile(repository, reg, filt, prof_data):
        """
        get_cells_onedprofile(repository, reg, filt, prof_data)
        
        
        Defined at profiles_module.fpp lines 308-508
        
        Parameters
        ----------
        repository : str
        reg : Region
        filt : Filter
        prof_data : Profile_Handler
        
        """
        _amr2.f90wrap_get_cells_onedprofile(repository=repository, reg=reg._handle, \
            filt=filt._handle, prof_data=prof_data._handle)
    
    @staticmethod
    def onedprofile(repository, reg, filt, prof_data, lmax, logscale):
        """
        onedprofile(repository, reg, filt, prof_data, lmax, logscale)
        
        
        Defined at profiles_module.fpp lines 510-529
        
        Parameters
        ----------
        repository : str
        reg : Region
        filt : Filter
        prof_data : Profile_Handler
        lmax : int
        logscale : bool
        
        """
        _amr2.f90wrap_onedprofile(repository=repository, reg=reg._handle, \
            filt=filt._handle, prof_data=prof_data._handle, lmax=lmax, \
            logscale=logscale)
    
    @staticmethod
    def twodprofile(repository, reg, filt, prof_data, lmax, logscale):
        """
        twodprofile(repository, reg, filt, prof_data, lmax, logscale)
        
        
        Defined at profiles_module.fpp lines 531-554
        
        Parameters
        ----------
        repository : str
        reg : Region
        filt : Filter
        prof_data : Profile_Handler_Twod
        lmax : int
        logscale : bool
        
        """
        _amr2.f90wrap_twodprofile(repository=repository, reg=reg._handle, \
            filt=filt._handle, prof_data=prof_data._handle, lmax=lmax, \
            logscale=logscale)
    
    @staticmethod
    def get_cells_twodprofile(repository, reg, filt, prof_data, logscale):
        """
        get_cells_twodprofile(repository, reg, filt, prof_data, logscale)
        
        
        Defined at profiles_module.fpp lines 556-756
        
        Parameters
        ----------
        repository : str
        reg : Region
        filt : Filter
        prof_data : Profile_Handler_Twod
        logscale : bool
        
        """
        _amr2.f90wrap_get_cells_twodprofile(repository=repository, reg=reg._handle, \
            filt=filt._handle, prof_data=prof_data._handle, logscale=logscale)
    
    _dt_array_initialisers = []
    

amr_profiles = Amr_Profiles()

class Amr_Integrator(f90wrap.runtime.FortranModule):
    """
    Module amr_integrator
    
    
    Defined at integrator_module.fpp lines 24-311
    
    """
    @f90wrap.runtime.register_class("amr2.amr_region_attrs")
    class amr_region_attrs(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=amr_region_attrs)
        
        
        Defined at integrator_module.fpp lines 28-34
        
        """
        def __init__(self, handle=None):
            """
            self = Amr_Region_Attrs()
            
            
            Defined at integrator_module.fpp lines 28-34
            
            
            Returns
            -------
            this : Amr_Region_Attrs
            	Object to be constructed
            
            
            Automatically generated constructor for amr_region_attrs
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _amr2.f90wrap_amr_region_attrs_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Amr_Region_Attrs
            
            
            Defined at integrator_module.fpp lines 28-34
            
            Parameters
            ----------
            this : Amr_Region_Attrs
            	Object to be destructed
            
            
            Automatically generated destructor for amr_region_attrs
            """
            if self._alloc:
                _amr2.f90wrap_amr_region_attrs_finalise(this=self._handle)
        
        @property
        def nvars(self):
            """
            Element nvars ftype=integer  pytype=int
            
            
            Defined at integrator_module.fpp line 29
            
            """
            return _amr2.f90wrap_amr_region_attrs__get__nvars(self._handle)
        
        @nvars.setter
        def nvars(self, nvars):
            _amr2.f90wrap_amr_region_attrs__set__nvars(self._handle, nvars)
        
        @property
        def nfilter(self):
            """
            Element nfilter ftype=integer  pytype=int
            
            
            Defined at integrator_module.fpp line 29
            
            """
            return _amr2.f90wrap_amr_region_attrs__get__nfilter(self._handle)
        
        @nfilter.setter
        def nfilter(self, nfilter):
            _amr2.f90wrap_amr_region_attrs__set__nfilter(self._handle, nfilter)
        
        @property
        def varnames(self):
            """
            Element varnames ftype=character(128) pytype=str
            
            
            Defined at integrator_module.fpp line 30
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_amr_region_attrs__array__varnames(self._handle)
            if array_handle in self._arrays:
                varnames = self._arrays[array_handle]
            else:
                varnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_amr_region_attrs__array__varnames)
                self._arrays[array_handle] = varnames
            return varnames
        
        @varnames.setter
        def varnames(self, varnames):
            self.varnames[...] = varnames
        
        @property
        def nwvars(self):
            """
            Element nwvars ftype=integer  pytype=int
            
            
            Defined at integrator_module.fpp line 31
            
            """
            return _amr2.f90wrap_amr_region_attrs__get__nwvars(self._handle)
        
        @nwvars.setter
        def nwvars(self, nwvars):
            _amr2.f90wrap_amr_region_attrs__set__nwvars(self._handle, nwvars)
        
        def init_array_filters(self):
            self.filters = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _amr2.f90wrap_amr_region_attrs__array_getitem__filters,
                                            _amr2.f90wrap_amr_region_attrs__array_setitem__filters,
                                            _amr2.f90wrap_amr_region_attrs__array_len__filters,
                                            """
            Element filters ftype=type(filter) pytype=Filter
            
            
            Defined at integrator_module.fpp line 32
            
            """, Filtering.filter)
            return self.filters
        
        @property
        def wvarnames(self):
            """
            Element wvarnames ftype=character(128) pytype=str
            
            
            Defined at integrator_module.fpp line 33
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_amr_region_attrs__array__wvarnames(self._handle)
            if array_handle in self._arrays:
                wvarnames = self._arrays[array_handle]
            else:
                wvarnames = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_amr_region_attrs__array__wvarnames)
                self._arrays[array_handle] = wvarnames
            return wvarnames
        
        @wvarnames.setter
        def wvarnames(self, wvarnames):
            self.wvarnames[...] = wvarnames
        
        @property
        def data(self):
            """
            Element data ftype=real(dbl) pytype=float
            
            
            Defined at integrator_module.fpp line 34
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _amr2.f90wrap_amr_region_attrs__array__data(self._handle)
            if array_handle in self._arrays:
                data = self._arrays[array_handle]
            else:
                data = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _amr2.f90wrap_amr_region_attrs__array__data)
                self._arrays[array_handle] = data
            return data
        
        @data.setter
        def data(self, data):
            self.data[...] = data
        
        def __str__(self):
            ret = ['<amr_region_attrs>{\n']
            ret.append('    nvars : ')
            ret.append(repr(self.nvars))
            ret.append(',\n    nfilter : ')
            ret.append(repr(self.nfilter))
            ret.append(',\n    varnames : ')
            ret.append(repr(self.varnames))
            ret.append(',\n    nwvars : ')
            ret.append(repr(self.nwvars))
            ret.append(',\n    wvarnames : ')
            ret.append(repr(self.wvarnames))
            ret.append(',\n    data : ')
            ret.append(repr(self.data))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = [init_array_filters]
        
    
    @staticmethod
    def allocate_amr_regions_attrs(self):
        """
        allocate_amr_regions_attrs(self)
        
        
        Defined at integrator_module.fpp lines 37-43
        
        Parameters
        ----------
        attrs : Amr_Region_Attrs
        
        """
        _amr2.f90wrap_allocate_amr_regions_attrs(attrs=self._handle)
    
    @staticmethod
    def extract_data(self, pos, cellvars, cellsize, attrs, ifilt):
        """
        extract_data(self, pos, cellvars, cellsize, attrs, ifilt)
        
        
        Defined at integrator_module.fpp lines 45-81
        
        Parameters
        ----------
        reg : Region
        pos : float array
        cellvars : float array
        cellsize : float
        attrs : Amr_Region_Attrs
        ifilt : int
        
        """
        _amr2.f90wrap_extract_data(reg=self._handle, pos=pos, cellvars=cellvars, \
            cellsize=cellsize, attrs=attrs._handle, ifilt=ifilt)
    
    @staticmethod
    def renormalise(self):
        """
        renormalise(self)
        
        
        Defined at integrator_module.fpp lines 83-97
        
        Parameters
        ----------
        attrs : Amr_Region_Attrs
        
        """
        _amr2.f90wrap_renormalise(attrs=self._handle)
    
    @staticmethod
    def integrate_region(repository, reg, attrs):
        """
        integrate_region(repository, reg, attrs)
        
        
        Defined at integrator_module.fpp lines 99-311
        
        Parameters
        ----------
        repository : str
        reg : Region
        attrs : Amr_Region_Attrs
        
        """
        _amr2.f90wrap_integrate_region(repository=repository, reg=reg._handle, \
            attrs=attrs._handle)
    
    _dt_array_initialisers = []
    

amr_integrator = Amr_Integrator()

