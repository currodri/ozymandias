! Module obs_instruments defined in file amr2map.fpp

subroutine f90wrap_camera__get__centre(this, f90wrap_centre)
    use obs_instruments, only: camera
    use vectors, only: vector
    implicit none
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(camera_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_centre(2)
    type(vector_ptr_type) :: centre_ptr
    
    this_ptr = transfer(this, this_ptr)
    centre_ptr%p => this_ptr%p%centre
    f90wrap_centre = transfer(centre_ptr,f90wrap_centre)
end subroutine f90wrap_camera__get__centre

subroutine f90wrap_camera__set__centre(this, f90wrap_centre)
    use obs_instruments, only: camera
    use vectors, only: vector
    implicit none
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(camera_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_centre(2)
    type(vector_ptr_type) :: centre_ptr
    
    this_ptr = transfer(this, this_ptr)
    centre_ptr = transfer(f90wrap_centre,centre_ptr)
    this_ptr%p%centre = centre_ptr%p
end subroutine f90wrap_camera__set__centre

subroutine f90wrap_camera__get__los_axis(this, f90wrap_los_axis)
    use obs_instruments, only: camera
    use vectors, only: vector
    implicit none
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(camera_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_los_axis(2)
    type(vector_ptr_type) :: los_axis_ptr
    
    this_ptr = transfer(this, this_ptr)
    los_axis_ptr%p => this_ptr%p%los_axis
    f90wrap_los_axis = transfer(los_axis_ptr,f90wrap_los_axis)
end subroutine f90wrap_camera__get__los_axis

subroutine f90wrap_camera__set__los_axis(this, f90wrap_los_axis)
    use obs_instruments, only: camera
    use vectors, only: vector
    implicit none
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(camera_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_los_axis(2)
    type(vector_ptr_type) :: los_axis_ptr
    
    this_ptr = transfer(this, this_ptr)
    los_axis_ptr = transfer(f90wrap_los_axis,los_axis_ptr)
    this_ptr%p%los_axis = los_axis_ptr%p
end subroutine f90wrap_camera__set__los_axis

subroutine f90wrap_camera__get__up_vector(this, f90wrap_up_vector)
    use obs_instruments, only: camera
    use vectors, only: vector
    implicit none
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(camera_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_up_vector(2)
    type(vector_ptr_type) :: up_vector_ptr
    
    this_ptr = transfer(this, this_ptr)
    up_vector_ptr%p => this_ptr%p%up_vector
    f90wrap_up_vector = transfer(up_vector_ptr,f90wrap_up_vector)
end subroutine f90wrap_camera__get__up_vector

subroutine f90wrap_camera__set__up_vector(this, f90wrap_up_vector)
    use obs_instruments, only: camera
    use vectors, only: vector
    implicit none
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(camera_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_up_vector(2)
    type(vector_ptr_type) :: up_vector_ptr
    
    this_ptr = transfer(this, this_ptr)
    up_vector_ptr = transfer(f90wrap_up_vector,up_vector_ptr)
    this_ptr%p%up_vector = up_vector_ptr%p
end subroutine f90wrap_camera__set__up_vector

subroutine f90wrap_camera__array__region_size(this, nd, dtype, dshape, dloc)
    use obs_instruments, only: camera
    implicit none
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    integer, intent(in) :: this(2)
    type(camera_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%region_size)
    dloc = loc(this_ptr%p%region_size)
end subroutine f90wrap_camera__array__region_size

subroutine f90wrap_camera__get__distance(this, f90wrap_distance)
    use obs_instruments, only: camera
    implicit none
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    integer, intent(in)   :: this(2)
    type(camera_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_distance
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_distance = this_ptr%p%distance
end subroutine f90wrap_camera__get__distance

subroutine f90wrap_camera__set__distance(this, f90wrap_distance)
    use obs_instruments, only: camera
    implicit none
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    integer, intent(in)   :: this(2)
    type(camera_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_distance
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%distance = f90wrap_distance
end subroutine f90wrap_camera__set__distance

subroutine f90wrap_camera__get__far_cut_depth(this, f90wrap_far_cut_depth)
    use obs_instruments, only: camera
    implicit none
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    integer, intent(in)   :: this(2)
    type(camera_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_far_cut_depth
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_far_cut_depth = this_ptr%p%far_cut_depth
end subroutine f90wrap_camera__get__far_cut_depth

subroutine f90wrap_camera__set__far_cut_depth(this, f90wrap_far_cut_depth)
    use obs_instruments, only: camera
    implicit none
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    integer, intent(in)   :: this(2)
    type(camera_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_far_cut_depth
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%far_cut_depth = f90wrap_far_cut_depth
end subroutine f90wrap_camera__set__far_cut_depth

subroutine f90wrap_camera__get__map_max_size(this, f90wrap_map_max_size)
    use obs_instruments, only: camera
    implicit none
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    integer, intent(in)   :: this(2)
    type(camera_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_map_max_size
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_map_max_size = this_ptr%p%map_max_size
end subroutine f90wrap_camera__get__map_max_size

subroutine f90wrap_camera__set__map_max_size(this, f90wrap_map_max_size)
    use obs_instruments, only: camera
    implicit none
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    integer, intent(in)   :: this(2)
    type(camera_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_map_max_size
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%map_max_size = f90wrap_map_max_size
end subroutine f90wrap_camera__set__map_max_size

subroutine f90wrap_camera_initialise(this)
    use obs_instruments, only: camera
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type(camera_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_camera_initialise

subroutine f90wrap_camera_finalise(this)
    use obs_instruments, only: camera
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type(camera_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_camera_finalise

subroutine f90wrap_log2(ret_log2, x)
    use obs_instruments, only: log2
    implicit none
    
    real, intent(out) :: ret_log2
    real(8), intent(in) :: x
    ret_log2 = log2(x=x)
end subroutine f90wrap_log2

subroutine f90wrap_init_camera(centre, los_axis, up_vector, region_size, distance, far_cut_depth, ret_init_camera, &
    map_max_size, n0)
    use vectors, only: vector
    use obs_instruments, only: camera, init_camera
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type(vector_ptr_type) :: centre_ptr
    integer, intent(in), dimension(2) :: centre
    type(vector_ptr_type) :: los_axis_ptr
    integer, intent(in), dimension(2) :: los_axis
    type(vector_ptr_type) :: up_vector_ptr
    integer, intent(in), dimension(2) :: up_vector
    real(8), intent(in), dimension(n0) :: region_size
    real(8), intent(in) :: distance
    real(8), intent(in) :: far_cut_depth
    type(camera_ptr_type) :: ret_init_camera_ptr
    integer, intent(out), dimension(2) :: ret_init_camera
    integer, intent(in) :: map_max_size
    integer :: n0
    !f2py intent(hide), depend(region_size) :: n0 = shape(region_size,0)
    centre_ptr = transfer(centre, centre_ptr)
    los_axis_ptr = transfer(los_axis, los_axis_ptr)
    up_vector_ptr = transfer(up_vector, up_vector_ptr)
    allocate(ret_init_camera_ptr%p)
    ret_init_camera_ptr%p = init_camera(centre=centre_ptr%p, los_axis=los_axis_ptr%p, up_vector=up_vector_ptr%p, &
        region_size=region_size, distance=distance, far_cut_depth=far_cut_depth, map_max_size=map_max_size)
    ret_init_camera = transfer(ret_init_camera_ptr, ret_init_camera)
end subroutine f90wrap_init_camera

subroutine f90wrap_get_required_resolution(ret_get_required_resolution, cam)
    use obs_instruments, only: camera, get_required_resolution
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    integer, intent(out) :: ret_get_required_resolution
    type(camera_ptr_type) :: cam_ptr
    integer, intent(in), dimension(2) :: cam
    cam_ptr = transfer(cam, cam_ptr)
    ret_get_required_resolution = get_required_resolution(cam=cam_ptr%p)
end subroutine f90wrap_get_required_resolution

subroutine f90wrap_get_map_size(cam, n_map, n0)
    use obs_instruments, only: camera, get_map_size
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type(camera_ptr_type) :: cam_ptr
    integer, intent(in), dimension(2) :: cam
    integer, intent(inout), dimension(n0) :: n_map
    integer :: n0
    !f2py intent(hide), depend(n_map) :: n0 = shape(n_map,0)
    cam_ptr = transfer(cam, cam_ptr)
    call get_map_size(cam=cam_ptr%p, n_map=n_map)
end subroutine f90wrap_get_map_size

subroutine f90wrap_get_map_box(cam, box)
    use obs_instruments, only: camera, get_map_box
    use geometrical_regions, only: region
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type(camera_ptr_type) :: cam_ptr
    integer, intent(in), dimension(2) :: cam
    type(region_ptr_type) :: box_ptr
    integer, intent(in), dimension(2) :: box
    cam_ptr = transfer(cam, cam_ptr)
    box_ptr = transfer(box, box_ptr)
    call get_map_box(cam=cam_ptr%p, box=box_ptr%p)
end subroutine f90wrap_get_map_box

subroutine f90wrap_get_camera_basis(cam, cam_basis)
    use obs_instruments, only: camera, get_camera_basis
    use basis_representations, only: basis
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type basis_ptr_type
        type(basis), pointer :: p => NULL()
    end type basis_ptr_type
    type(camera_ptr_type) :: cam_ptr
    integer, intent(in), dimension(2) :: cam
    type(basis_ptr_type) :: cam_basis_ptr
    integer, intent(in), dimension(2) :: cam_basis
    cam_ptr = transfer(cam, cam_ptr)
    cam_basis_ptr = transfer(cam_basis, cam_basis_ptr)
    call get_camera_basis(cam=cam_ptr%p, cam_basis=cam_basis_ptr%p)
end subroutine f90wrap_get_camera_basis

subroutine f90wrap_los_transformation(cam, trans_matrix, n0, n1)
    use obs_instruments, only: camera, los_transformation
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type(camera_ptr_type) :: cam_ptr
    integer, intent(in), dimension(2) :: cam
    real(8), intent(inout), dimension(n0,n1) :: trans_matrix
    integer :: n0
    !f2py intent(hide), depend(trans_matrix) :: n0 = shape(trans_matrix,0)
    integer :: n1
    !f2py intent(hide), depend(trans_matrix) :: n1 = shape(trans_matrix,1)
    cam_ptr = transfer(cam, cam_ptr)
    call los_transformation(cam=cam_ptr%p, trans_matrix=trans_matrix)
end subroutine f90wrap_los_transformation

subroutine f90wrap_get_bounding_box(cam, bbox)
    use obs_instruments, only: camera, get_bounding_box
    use geometrical_regions, only: region
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type(camera_ptr_type) :: cam_ptr
    integer, intent(in), dimension(2) :: cam
    type(region_ptr_type) :: bbox_ptr
    integer, intent(in), dimension(2) :: bbox
    cam_ptr = transfer(cam, cam_ptr)
    bbox_ptr = transfer(bbox, bbox_ptr)
    call get_bounding_box(cam=cam_ptr%p, bbox=bbox_ptr%p)
end subroutine f90wrap_get_bounding_box

subroutine f90wrap_deproject_points(cam, npoints, points, n0, n1)
    use obs_instruments, only: deproject_points, camera
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type(camera_ptr_type) :: cam_ptr
    integer, intent(in), dimension(2) :: cam
    integer, intent(in) :: npoints
    real(8), intent(inout), dimension(n0,n1) :: points
    integer :: n0
    !f2py intent(hide), depend(points) :: n0 = shape(points,0)
    integer :: n1
    !f2py intent(hide), depend(points) :: n1 = shape(points,1)
    cam_ptr = transfer(cam, cam_ptr)
    call deproject_points(cam=cam_ptr%p, npoints=npoints, points=points)
end subroutine f90wrap_deproject_points

subroutine f90wrap_project_points(cam, npoints, points, n0, n1)
    use obs_instruments, only: camera, project_points
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type(camera_ptr_type) :: cam_ptr
    integer, intent(in), dimension(2) :: cam
    integer, intent(in) :: npoints
    real(8), intent(inout), dimension(n0,n1) :: points
    integer :: n0
    !f2py intent(hide), depend(points) :: n0 = shape(points,0)
    integer :: n1
    !f2py intent(hide), depend(points) :: n1 = shape(points,1)
    cam_ptr = transfer(cam, cam_ptr)
    call project_points(cam=cam_ptr%p, npoints=npoints, points=points)
end subroutine f90wrap_project_points

! End of module obs_instruments defined in file amr2map.fpp

! Module amr_map defined in file amr2map.fpp

subroutine f90wrap_projection(repository, cam, bulk_velocity)
    use amr_map, only: projection
    use obs_instruments, only: camera
    use vectors, only: vector
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    character(128), intent(in) :: repository
    type(camera_ptr_type) :: cam_ptr
    integer, intent(in), dimension(2) :: cam
    type(vector_ptr_type) :: bulk_velocity_ptr
    integer, intent(in), dimension(2) :: bulk_velocity
    cam_ptr = transfer(cam, cam_ptr)
    bulk_velocity_ptr = transfer(bulk_velocity, bulk_velocity_ptr)
    call projection(repository=repository, cam=cam_ptr%p, bulk_velocity=bulk_velocity_ptr%p)
end subroutine f90wrap_projection

! End of module amr_map defined in file amr2map.fpp

