! Module coordinate_systems defined in file coordinates_module.fpp

subroutine f90wrap_r_sphere(ret_r_sphere, p)
    use vectors, only: vector
    use coordinate_systems, only: r_sphere
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    real(8), intent(out) :: ret_r_sphere
    type(vector_ptr_type) :: p_ptr
    integer, intent(in), dimension(2) :: p
    p_ptr = transfer(p, p_ptr)
    ret_r_sphere = r_sphere(p=p_ptr%p)
end subroutine f90wrap_r_sphere

subroutine f90wrap_theta_sphere(ret_theta_sphere, p)
    use vectors, only: vector
    use coordinate_systems, only: theta_sphere
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    real(8), intent(out) :: ret_theta_sphere
    type(vector_ptr_type) :: p_ptr
    integer, intent(in), dimension(2) :: p
    p_ptr = transfer(p, p_ptr)
    ret_theta_sphere = theta_sphere(p=p_ptr%p)
end subroutine f90wrap_theta_sphere

subroutine f90wrap_phi_sphere(ret_phi_sphere, p)
    use vectors, only: vector
    use coordinate_systems, only: phi_sphere
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    real(8), intent(out) :: ret_phi_sphere
    type(vector_ptr_type) :: p_ptr
    integer, intent(in), dimension(2) :: p
    p_ptr = transfer(p, p_ptr)
    ret_phi_sphere = phi_sphere(p=p_ptr%p)
end subroutine f90wrap_phi_sphere

subroutine f90wrap_r_cyl(ret_r_cyl, p)
    use vectors, only: vector
    use coordinate_systems, only: r_cyl
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    real(8), intent(out) :: ret_r_cyl
    type(vector_ptr_type) :: p_ptr
    integer, intent(in), dimension(2) :: p
    p_ptr = transfer(p, p_ptr)
    ret_r_cyl = r_cyl(p=p_ptr%p)
end subroutine f90wrap_r_cyl

subroutine f90wrap_phi_cyl(ret_phi_cyl, p)
    use vectors, only: vector
    use coordinate_systems, only: phi_cyl
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    real(8), intent(out) :: ret_phi_cyl
    type(vector_ptr_type) :: p_ptr
    integer, intent(in), dimension(2) :: p
    p_ptr = transfer(p, p_ptr)
    ret_phi_cyl = phi_cyl(p=p_ptr%p)
end subroutine f90wrap_phi_cyl

subroutine f90wrap_spherical_basis_from_cartesian(p, spher_basis)
    use vectors, only: vector
    use coordinate_systems, only: spherical_basis_from_cartesian
    use basis_representations, only: basis
    implicit none
    
    type basis_ptr_type
        type(basis), pointer :: p => NULL()
    end type basis_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type(vector_ptr_type) :: p_ptr
    integer, intent(in), dimension(2) :: p
    type(basis_ptr_type) :: spher_basis_ptr
    integer, intent(in), dimension(2) :: spher_basis
    p_ptr = transfer(p, p_ptr)
    spher_basis_ptr = transfer(spher_basis, spher_basis_ptr)
    call spherical_basis_from_cartesian(p=p_ptr%p, spher_basis=spher_basis_ptr%p)
end subroutine f90wrap_spherical_basis_from_cartesian

subroutine f90wrap_cylindrical_basis_from_cartesian(p, cyl_basis)
    use vectors, only: vector
    use basis_representations, only: basis
    use coordinate_systems, only: cylindrical_basis_from_cartesian
    implicit none
    
    type basis_ptr_type
        type(basis), pointer :: p => NULL()
    end type basis_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type(vector_ptr_type) :: p_ptr
    integer, intent(in), dimension(2) :: p
    type(basis_ptr_type) :: cyl_basis_ptr
    integer, intent(in), dimension(2) :: cyl_basis
    p_ptr = transfer(p, p_ptr)
    cyl_basis_ptr = transfer(cyl_basis, cyl_basis_ptr)
    call cylindrical_basis_from_cartesian(p=p_ptr%p, cyl_basis=cyl_basis_ptr%p)
end subroutine f90wrap_cylindrical_basis_from_cartesian

subroutine f90wrap_new_z_coordinates(axis, transformation_matrix, errormsg)
    use vectors, only: vector
    use coordinate_systems, only: new_z_coordinates
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type(vector_ptr_type) :: axis_ptr
    integer, intent(in), dimension(2) :: axis
    real(8), dimension(3,3), intent(inout) :: transformation_matrix
    integer, intent(inout) :: errormsg
    axis_ptr = transfer(axis, axis_ptr)
    call new_z_coordinates(axis=axis_ptr%p, transformation_matrix=transformation_matrix, errormsg=errormsg)
end subroutine f90wrap_new_z_coordinates

! End of module coordinate_systems defined in file coordinates_module.fpp

! Module geometrical_regions defined in file coordinates_module.fpp

subroutine f90wrap_region__get__name(this, f90wrap_name)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    character(128), intent(out) :: f90wrap_name
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_name = this_ptr%p%name
end subroutine f90wrap_region__get__name

subroutine f90wrap_region__set__name(this, f90wrap_name)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    character(128), intent(in) :: f90wrap_name
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%name = f90wrap_name
end subroutine f90wrap_region__set__name

subroutine f90wrap_region__get__centre(this, f90wrap_centre)
    use geometrical_regions, only: region
    use vectors, only: vector
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_centre(2)
    type(vector_ptr_type) :: centre_ptr
    
    this_ptr = transfer(this, this_ptr)
    centre_ptr%p => this_ptr%p%centre
    f90wrap_centre = transfer(centre_ptr,f90wrap_centre)
end subroutine f90wrap_region__get__centre

subroutine f90wrap_region__set__centre(this, f90wrap_centre)
    use geometrical_regions, only: region
    use vectors, only: vector
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_centre(2)
    type(vector_ptr_type) :: centre_ptr
    
    this_ptr = transfer(this, this_ptr)
    centre_ptr = transfer(f90wrap_centre,centre_ptr)
    this_ptr%p%centre = centre_ptr%p
end subroutine f90wrap_region__set__centre

subroutine f90wrap_region__get__axis(this, f90wrap_axis)
    use geometrical_regions, only: region
    use vectors, only: vector
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_axis(2)
    type(vector_ptr_type) :: axis_ptr
    
    this_ptr = transfer(this, this_ptr)
    axis_ptr%p => this_ptr%p%axis
    f90wrap_axis = transfer(axis_ptr,f90wrap_axis)
end subroutine f90wrap_region__get__axis

subroutine f90wrap_region__set__axis(this, f90wrap_axis)
    use geometrical_regions, only: region
    use vectors, only: vector
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_axis(2)
    type(vector_ptr_type) :: axis_ptr
    
    this_ptr = transfer(this, this_ptr)
    axis_ptr = transfer(f90wrap_axis,axis_ptr)
    this_ptr%p%axis = axis_ptr%p
end subroutine f90wrap_region__set__axis

subroutine f90wrap_region__get__bulk_velocity(this, f90wrap_bulk_velocity)
    use geometrical_regions, only: region
    use vectors, only: vector
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_bulk_velocity(2)
    type(vector_ptr_type) :: bulk_velocity_ptr
    
    this_ptr = transfer(this, this_ptr)
    bulk_velocity_ptr%p => this_ptr%p%bulk_velocity
    f90wrap_bulk_velocity = transfer(bulk_velocity_ptr,f90wrap_bulk_velocity)
end subroutine f90wrap_region__get__bulk_velocity

subroutine f90wrap_region__set__bulk_velocity(this, f90wrap_bulk_velocity)
    use geometrical_regions, only: region
    use vectors, only: vector
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_bulk_velocity(2)
    type(vector_ptr_type) :: bulk_velocity_ptr
    
    this_ptr = transfer(this, this_ptr)
    bulk_velocity_ptr = transfer(f90wrap_bulk_velocity,bulk_velocity_ptr)
    this_ptr%p%bulk_velocity = bulk_velocity_ptr%p
end subroutine f90wrap_region__set__bulk_velocity

subroutine f90wrap_region__get__xmin(this, f90wrap_xmin)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_xmin
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_xmin = this_ptr%p%xmin
end subroutine f90wrap_region__get__xmin

subroutine f90wrap_region__set__xmin(this, f90wrap_xmin)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_xmin
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%xmin = f90wrap_xmin
end subroutine f90wrap_region__set__xmin

subroutine f90wrap_region__get__xmax(this, f90wrap_xmax)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_xmax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_xmax = this_ptr%p%xmax
end subroutine f90wrap_region__get__xmax

subroutine f90wrap_region__set__xmax(this, f90wrap_xmax)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_xmax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%xmax = f90wrap_xmax
end subroutine f90wrap_region__set__xmax

subroutine f90wrap_region__get__ymin(this, f90wrap_ymin)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_ymin
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ymin = this_ptr%p%ymin
end subroutine f90wrap_region__get__ymin

subroutine f90wrap_region__set__ymin(this, f90wrap_ymin)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_ymin
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ymin = f90wrap_ymin
end subroutine f90wrap_region__set__ymin

subroutine f90wrap_region__get__ymax(this, f90wrap_ymax)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_ymax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ymax = this_ptr%p%ymax
end subroutine f90wrap_region__get__ymax

subroutine f90wrap_region__set__ymax(this, f90wrap_ymax)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_ymax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ymax = f90wrap_ymax
end subroutine f90wrap_region__set__ymax

subroutine f90wrap_region__get__zmin(this, f90wrap_zmin)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_zmin
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_zmin = this_ptr%p%zmin
end subroutine f90wrap_region__get__zmin

subroutine f90wrap_region__set__zmin(this, f90wrap_zmin)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_zmin
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%zmin = f90wrap_zmin
end subroutine f90wrap_region__set__zmin

subroutine f90wrap_region__get__zmax(this, f90wrap_zmax)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_zmax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_zmax = this_ptr%p%zmax
end subroutine f90wrap_region__get__zmax

subroutine f90wrap_region__set__zmax(this, f90wrap_zmax)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_zmax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%zmax = f90wrap_zmax
end subroutine f90wrap_region__set__zmax

subroutine f90wrap_region__get__rmin(this, f90wrap_rmin)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_rmin
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_rmin = this_ptr%p%rmin
end subroutine f90wrap_region__get__rmin

subroutine f90wrap_region__set__rmin(this, f90wrap_rmin)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_rmin
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%rmin = f90wrap_rmin
end subroutine f90wrap_region__set__rmin

subroutine f90wrap_region__get__rmax(this, f90wrap_rmax)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_rmax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_rmax = this_ptr%p%rmax
end subroutine f90wrap_region__get__rmax

subroutine f90wrap_region__set__rmax(this, f90wrap_rmax)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_rmax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%rmax = f90wrap_rmax
end subroutine f90wrap_region__set__rmax

subroutine f90wrap_region__get__angle(this, f90wrap_angle)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_angle
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_angle = this_ptr%p%angle
end subroutine f90wrap_region__get__angle

subroutine f90wrap_region__set__angle(this, f90wrap_angle)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_angle
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%angle = f90wrap_angle
end subroutine f90wrap_region__set__angle

subroutine f90wrap_region__get__criteria_name(this, f90wrap_criteria_name)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    character(128), intent(out) :: f90wrap_criteria_name
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_criteria_name = this_ptr%p%criteria_name
end subroutine f90wrap_region__get__criteria_name

subroutine f90wrap_region__set__criteria_name(this, f90wrap_criteria_name)
    use geometrical_regions, only: region
    implicit none
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    integer, intent(in)   :: this(2)
    type(region_ptr_type) :: this_ptr
    character(128), intent(in) :: f90wrap_criteria_name
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%criteria_name = f90wrap_criteria_name
end subroutine f90wrap_region__set__criteria_name

subroutine f90wrap_region_initialise(this)
    use geometrical_regions, only: region
    implicit none
    
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type(region_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_region_initialise

subroutine f90wrap_region_finalise(this)
    use geometrical_regions, only: region
    implicit none
    
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type(region_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_region_finalise

subroutine f90wrap_limits(reg, lim, n0, n1)
    use geometrical_regions, only: region, limits
    implicit none
    
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    real(8), intent(inout), dimension(n0,n1) :: lim
    integer :: n0
    !f2py intent(hide), depend(lim) :: n0 = shape(lim,0)
    integer :: n1
    !f2py intent(hide), depend(lim) :: n1 = shape(lim,1)
    reg_ptr = transfer(reg, reg_ptr)
    call limits(reg=reg_ptr%p, lim=lim)
end subroutine f90wrap_limits

subroutine f90wrap_checkifinside(pos, reg, ok, distance, n0)
    use geometrical_regions, only: checkifinside, region
    implicit none
    
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    real(8), intent(in), dimension(n0) :: pos
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    logical, intent(inout) :: ok
    real(8), intent(inout) :: distance
    integer :: n0
    !f2py intent(hide), depend(pos) :: n0 = shape(pos,0)
    reg_ptr = transfer(reg, reg_ptr)
    call checkifinside(pos=pos, reg=reg_ptr%p, ok=ok, distance=distance)
end subroutine f90wrap_checkifinside

subroutine f90wrap_cube(p, reg, ok, distance)
    use geometrical_regions, only: region, cube
    use vectors, only: vector
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type(vector_ptr_type) :: p_ptr
    integer, intent(in), dimension(2) :: p
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    logical, intent(inout) :: ok
    real(8), intent(inout) :: distance
    p_ptr = transfer(p, p_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    call cube(p=p_ptr%p, reg=reg_ptr%p, ok=ok, distance=distance)
end subroutine f90wrap_cube

subroutine f90wrap_sphere(p, reg, ok, distance)
    use geometrical_regions, only: region, sphere
    use vectors, only: vector
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type(vector_ptr_type) :: p_ptr
    integer, intent(in), dimension(2) :: p
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    logical, intent(inout) :: ok
    real(8), intent(inout) :: distance
    p_ptr = transfer(p, p_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    call sphere(p=p_ptr%p, reg=reg_ptr%p, ok=ok, distance=distance)
end subroutine f90wrap_sphere

subroutine f90wrap_cylinder(p, reg, ok, distance)
    use geometrical_regions, only: region, cylinder
    use vectors, only: vector
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type(vector_ptr_type) :: p_ptr
    integer, intent(in), dimension(2) :: p
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    logical, intent(inout) :: ok
    real(8), intent(inout) :: distance
    p_ptr = transfer(p, p_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    call cylinder(p=p_ptr%p, reg=reg_ptr%p, ok=ok, distance=distance)
end subroutine f90wrap_cylinder

subroutine f90wrap_cone(p, reg, ok, distance)
    use geometrical_regions, only: cone, region
    use vectors, only: vector
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type(vector_ptr_type) :: p_ptr
    integer, intent(in), dimension(2) :: p
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    logical, intent(inout) :: ok
    real(8), intent(inout) :: distance
    p_ptr = transfer(p, p_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    call cone(p=p_ptr%p, reg=reg_ptr%p, ok=ok, distance=distance)
end subroutine f90wrap_cone

! End of module geometrical_regions defined in file coordinates_module.fpp

