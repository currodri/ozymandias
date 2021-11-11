! Module obs_instruments defined in file ramses2map.fpp

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
    use obs_instruments, only: camera, init_camera
    use vectors, only: vector
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
    use obs_instruments, only: get_bounding_box, camera
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
    use obs_instruments, only: camera, deproject_points
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

! End of module obs_instruments defined in file ramses2map.fpp

! Module maps defined in file ramses2map.fpp

subroutine f90wrap_projection_handler__get__pov(this, f90wrap_pov)
    use maps, only: projection_handler
    implicit none
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(projection_handler_ptr_type) :: this_ptr
    character(128), intent(out) :: f90wrap_pov
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_pov = this_ptr%p%pov
end subroutine f90wrap_projection_handler__get__pov

subroutine f90wrap_projection_handler__set__pov(this, f90wrap_pov)
    use maps, only: projection_handler
    implicit none
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(projection_handler_ptr_type) :: this_ptr
    character(128), intent(in) :: f90wrap_pov
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%pov = f90wrap_pov
end subroutine f90wrap_projection_handler__set__pov

subroutine f90wrap_projection_handler__get__nvars(this, f90wrap_nvars)
    use maps, only: projection_handler
    implicit none
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(projection_handler_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nvars
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nvars = this_ptr%p%nvars
end subroutine f90wrap_projection_handler__get__nvars

subroutine f90wrap_projection_handler__set__nvars(this, f90wrap_nvars)
    use maps, only: projection_handler
    implicit none
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(projection_handler_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nvars
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nvars = f90wrap_nvars
end subroutine f90wrap_projection_handler__set__nvars

subroutine f90wrap_projection_handler__array__varnames(this, nd, dtype, dshape, dloc)
    use maps, only: projection_handler
    implicit none
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    integer, intent(in) :: this(2)
    type(projection_handler_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%varnames)) then
        dshape(1:2) = (/len(this_ptr%p%varnames(1)), shape(this_ptr%p%varnames)/)
        dloc = loc(this_ptr%p%varnames)
    else
        dloc = 0
    end if
end subroutine f90wrap_projection_handler__array__varnames

subroutine f90wrap_projection_handler__get__weightvar(this, f90wrap_weightvar)
    use maps, only: projection_handler
    implicit none
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(projection_handler_ptr_type) :: this_ptr
    character(128), intent(out) :: f90wrap_weightvar
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_weightvar = this_ptr%p%weightvar
end subroutine f90wrap_projection_handler__get__weightvar

subroutine f90wrap_projection_handler__set__weightvar(this, f90wrap_weightvar)
    use maps, only: projection_handler
    implicit none
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(projection_handler_ptr_type) :: this_ptr
    character(128), intent(in) :: f90wrap_weightvar
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%weightvar = f90wrap_weightvar
end subroutine f90wrap_projection_handler__set__weightvar

subroutine f90wrap_projection_handler__array__toto(this, nd, dtype, dshape, dloc)
    use maps, only: projection_handler
    implicit none
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    integer, intent(in) :: this(2)
    type(projection_handler_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%toto)) then
        dshape(1:3) = shape(this_ptr%p%toto)
        dloc = loc(this_ptr%p%toto)
    else
        dloc = 0
    end if
end subroutine f90wrap_projection_handler__array__toto

subroutine f90wrap_projection_handler_initialise(this)
    use maps, only: projection_handler
    implicit none
    
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    type(projection_handler_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_projection_handler_initialise

subroutine f90wrap_projection_handler_finalise(this)
    use maps, only: projection_handler
    implicit none
    
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    type(projection_handler_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_projection_handler_finalise

subroutine f90wrap_allocate_projection_handler(proj)
    use maps, only: allocate_projection_handler, projection_handler
    implicit none
    
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    type(projection_handler_ptr_type) :: proj_ptr
    integer, intent(in), dimension(2) :: proj
    proj_ptr = transfer(proj, proj_ptr)
    call allocate_projection_handler(proj=proj_ptr%p)
end subroutine f90wrap_allocate_projection_handler

subroutine f90wrap_projection_hydro(repository, cam, bulk_velocity, proj)
    use maps, only: projection_hydro, projection_handler
    use obs_instruments, only: camera
    use vectors, only: vector
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    character(128), intent(in) :: repository
    type(camera_ptr_type) :: cam_ptr
    integer, intent(in), dimension(2) :: cam
    type(vector_ptr_type) :: bulk_velocity_ptr
    integer, intent(in), dimension(2) :: bulk_velocity
    type(projection_handler_ptr_type) :: proj_ptr
    integer, intent(in), dimension(2) :: proj
    cam_ptr = transfer(cam, cam_ptr)
    bulk_velocity_ptr = transfer(bulk_velocity, bulk_velocity_ptr)
    proj_ptr = transfer(proj, proj_ptr)
    call projection_hydro(repository=repository, cam=cam_ptr%p, bulk_velocity=bulk_velocity_ptr%p, proj=proj_ptr%p)
end subroutine f90wrap_projection_hydro

subroutine f90wrap_project_cells(repository, amr, bbox, varids, cam, proj)
    use obs_instruments, only: camera
    use io_ramses, only: hydroid, amr_info
    use maps, only: project_cells, projection_handler
    use geometrical_regions, only: region
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    character(128), intent(in) :: repository
    type(amr_info_ptr_type) :: amr_ptr
    integer, intent(in), dimension(2) :: amr
    type(region_ptr_type) :: bbox_ptr
    integer, intent(in), dimension(2) :: bbox
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    type(camera_ptr_type) :: cam_ptr
    integer, intent(in), dimension(2) :: cam
    type(projection_handler_ptr_type) :: proj_ptr
    integer, intent(in), dimension(2) :: proj
    amr_ptr = transfer(amr, amr_ptr)
    bbox_ptr = transfer(bbox, bbox_ptr)
    varids_ptr = transfer(varids, varids_ptr)
    cam_ptr = transfer(cam, cam_ptr)
    proj_ptr = transfer(proj, proj_ptr)
    call project_cells(repository=repository, amr=amr_ptr%p, bbox=bbox_ptr%p, varIDs=varids_ptr%p, cam=cam_ptr%p, &
        proj=proj_ptr%p)
end subroutine f90wrap_project_cells

subroutine f90wrap_projection_parts(repository, cam, bulk_velocity, proj, tag_file)
    use maps, only: projection_parts, projection_handler
    use obs_instruments, only: camera
    use vectors, only: vector
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    character(128), intent(in) :: repository
    type(camera_ptr_type) :: cam_ptr
    integer, intent(in), dimension(2) :: cam
    type(vector_ptr_type) :: bulk_velocity_ptr
    integer, intent(in), dimension(2) :: bulk_velocity
    type(projection_handler_ptr_type) :: proj_ptr
    integer, intent(in), dimension(2) :: proj
    character(128), intent(in), optional :: tag_file
    cam_ptr = transfer(cam, cam_ptr)
    bulk_velocity_ptr = transfer(bulk_velocity, bulk_velocity_ptr)
    proj_ptr = transfer(proj, proj_ptr)
    call projection_parts(repository=repository, cam=cam_ptr%p, bulk_velocity=bulk_velocity_ptr%p, proj=proj_ptr%p, &
        tag_file=tag_file)
end subroutine f90wrap_projection_parts

subroutine f90wrap_project_particles(repository, amr, sim, bbox, cam, proj, tag_file)
    use maps, only: projection_handler, project_particles
    use obs_instruments, only: camera
    use io_ramses, only: amr_info, sim_info
    use geometrical_regions, only: region
    implicit none
    
    type camera_ptr_type
        type(camera), pointer :: p => NULL()
    end type camera_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    character(128), intent(in) :: repository
    type(amr_info_ptr_type) :: amr_ptr
    integer, intent(in), dimension(2) :: amr
    type(sim_info_ptr_type) :: sim_ptr
    integer, intent(in), dimension(2) :: sim
    type(region_ptr_type) :: bbox_ptr
    integer, intent(in), dimension(2) :: bbox
    type(camera_ptr_type) :: cam_ptr
    integer, intent(in), dimension(2) :: cam
    type(projection_handler_ptr_type) :: proj_ptr
    integer, intent(in), dimension(2) :: proj
    character(128), intent(in), optional :: tag_file
    amr_ptr = transfer(amr, amr_ptr)
    sim_ptr = transfer(sim, sim_ptr)
    bbox_ptr = transfer(bbox, bbox_ptr)
    cam_ptr = transfer(cam, cam_ptr)
    proj_ptr = transfer(proj, proj_ptr)
    call project_particles(repository=repository, amr=amr_ptr%p, sim=sim_ptr%p, bbox=bbox_ptr%p, cam=cam_ptr%p, &
        proj=proj_ptr%p, tag_file=tag_file)
end subroutine f90wrap_project_particles

subroutine f90wrap_healpix_hydro(repository, reg, nside, proj)
    use maps, only: projection_handler, healpix_hydro
    use geometrical_regions, only: region
    implicit none
    
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    character(128), intent(in) :: repository
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    integer, intent(in) :: nside
    type(projection_handler_ptr_type) :: proj_ptr
    integer, intent(in), dimension(2) :: proj
    reg_ptr = transfer(reg, reg_ptr)
    proj_ptr = transfer(proj, proj_ptr)
    call healpix_hydro(repository=repository, reg=reg_ptr%p, nside=nside, proj=proj_ptr%p)
end subroutine f90wrap_healpix_hydro

subroutine f90wrap_project_cells_hpix(repository, amr, reg, varids, nside, proj)
    use io_ramses, only: hydroid, amr_info
    use maps, only: project_cells_hpix, projection_handler
    use geometrical_regions, only: region
    implicit none
    
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type projection_handler_ptr_type
        type(projection_handler), pointer :: p => NULL()
    end type projection_handler_ptr_type
    character(128), intent(in) :: repository
    type(amr_info_ptr_type) :: amr_ptr
    integer, intent(in), dimension(2) :: amr
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    integer, intent(in) :: nside
    type(projection_handler_ptr_type) :: proj_ptr
    integer, intent(in), dimension(2) :: proj
    amr_ptr = transfer(amr, amr_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    varids_ptr = transfer(varids, varids_ptr)
    proj_ptr = transfer(proj, proj_ptr)
    call project_cells_hpix(repository=repository, amr=amr_ptr%p, reg=reg_ptr%p, varIDs=varids_ptr%p, nside=nside, &
        proj=proj_ptr%p)
end subroutine f90wrap_project_cells_hpix

! End of module maps defined in file ramses2map.fpp

