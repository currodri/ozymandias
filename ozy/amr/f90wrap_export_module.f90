! Module export_amr defined in file export_module.fpp

subroutine f90wrap_chunk_handler__get__nvars(this, f90wrap_nvars)
    use export_amr, only: chunk_handler
    implicit none
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(chunk_handler_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nvars
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nvars = this_ptr%p%nvars
end subroutine f90wrap_chunk_handler__get__nvars

subroutine f90wrap_chunk_handler__set__nvars(this, f90wrap_nvars)
    use export_amr, only: chunk_handler
    implicit none
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(chunk_handler_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nvars
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nvars = f90wrap_nvars
end subroutine f90wrap_chunk_handler__set__nvars

subroutine f90wrap_chunk_handler__get__nx(this, f90wrap_nx)
    use export_amr, only: chunk_handler
    implicit none
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(chunk_handler_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nx = this_ptr%p%nx
end subroutine f90wrap_chunk_handler__get__nx

subroutine f90wrap_chunk_handler__set__nx(this, f90wrap_nx)
    use export_amr, only: chunk_handler
    implicit none
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(chunk_handler_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nx = f90wrap_nx
end subroutine f90wrap_chunk_handler__set__nx

subroutine f90wrap_chunk_handler__get__ny(this, f90wrap_ny)
    use export_amr, only: chunk_handler
    implicit none
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(chunk_handler_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_ny
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ny = this_ptr%p%ny
end subroutine f90wrap_chunk_handler__get__ny

subroutine f90wrap_chunk_handler__set__ny(this, f90wrap_ny)
    use export_amr, only: chunk_handler
    implicit none
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(chunk_handler_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_ny
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ny = f90wrap_ny
end subroutine f90wrap_chunk_handler__set__ny

subroutine f90wrap_chunk_handler__get__nz(this, f90wrap_nz)
    use export_amr, only: chunk_handler
    implicit none
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(chunk_handler_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nz
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nz = this_ptr%p%nz
end subroutine f90wrap_chunk_handler__get__nz

subroutine f90wrap_chunk_handler__set__nz(this, f90wrap_nz)
    use export_amr, only: chunk_handler
    implicit none
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(chunk_handler_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nz
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nz = f90wrap_nz
end subroutine f90wrap_chunk_handler__set__nz

subroutine f90wrap_chunk_handler__array__varnames(this, nd, dtype, dshape, dloc)
    use export_amr, only: chunk_handler
    implicit none
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    integer, intent(in) :: this(2)
    type(chunk_handler_ptr_type) :: this_ptr
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
end subroutine f90wrap_chunk_handler__array__varnames

subroutine f90wrap_chunk_handler__get__filt(this, f90wrap_filt)
    use export_amr, only: chunk_handler
    use filtering, only: filter
    implicit none
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    integer, intent(in)   :: this(2)
    type(chunk_handler_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_filt(2)
    type(filter_ptr_type) :: filt_ptr
    
    this_ptr = transfer(this, this_ptr)
    filt_ptr%p => this_ptr%p%filt
    f90wrap_filt = transfer(filt_ptr,f90wrap_filt)
end subroutine f90wrap_chunk_handler__get__filt

subroutine f90wrap_chunk_handler__set__filt(this, f90wrap_filt)
    use export_amr, only: chunk_handler
    use filtering, only: filter
    implicit none
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    integer, intent(in)   :: this(2)
    type(chunk_handler_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_filt(2)
    type(filter_ptr_type) :: filt_ptr
    
    this_ptr = transfer(this, this_ptr)
    filt_ptr = transfer(f90wrap_filt,filt_ptr)
    this_ptr%p%filt = filt_ptr%p
end subroutine f90wrap_chunk_handler__set__filt

subroutine f90wrap_chunk_handler__array__data(this, nd, dtype, dshape, dloc)
    use export_amr, only: chunk_handler
    implicit none
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    integer, intent(in) :: this(2)
    type(chunk_handler_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%data)) then
        dshape(1:4) = shape(this_ptr%p%data)
        dloc = loc(this_ptr%p%data)
    else
        dloc = 0
    end if
end subroutine f90wrap_chunk_handler__array__data

subroutine f90wrap_chunk_handler_initialise(this)
    use export_amr, only: chunk_handler
    implicit none
    
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    type(chunk_handler_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_chunk_handler_initialise

subroutine f90wrap_chunk_handler_finalise(this)
    use export_amr, only: chunk_handler
    implicit none
    
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    type(chunk_handler_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_chunk_handler_finalise

subroutine f90wrap_allocate_chunk_handler(chunk)
    use export_amr, only: allocate_chunk_handler, chunk_handler
    implicit none
    
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    type(chunk_handler_ptr_type) :: chunk_ptr
    integer, intent(in), dimension(2) :: chunk
    chunk_ptr = transfer(chunk, chunk_ptr)
    call allocate_chunk_handler(chunk=chunk_ptr%p)
end subroutine f90wrap_allocate_chunk_handler

subroutine f90wrap_get_unigrid_old(repository, reg, chunk)
    use export_amr, only: get_unigrid_old, chunk_handler
    use geometrical_regions, only: region
    implicit none
    
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    character(128), intent(in) :: repository
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(chunk_handler_ptr_type) :: chunk_ptr
    integer, intent(in), dimension(2) :: chunk
    reg_ptr = transfer(reg, reg_ptr)
    chunk_ptr = transfer(chunk, chunk_ptr)
    call get_unigrid_old(repository=repository, reg=reg_ptr%p, chunk=chunk_ptr%p)
end subroutine f90wrap_get_unigrid_old

subroutine f90wrap_get_unigrid(repository, reg, lmax, symlog, chunk)
    use export_amr, only: chunk_handler, get_unigrid
    use geometrical_regions, only: region
    implicit none
    
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type chunk_handler_ptr_type
        type(chunk_handler), pointer :: p => NULL()
    end type chunk_handler_ptr_type
    character(128), intent(in) :: repository
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    integer, intent(in) :: lmax
    logical, intent(in) :: symlog
    type(chunk_handler_ptr_type) :: chunk_ptr
    integer, intent(in), dimension(2) :: chunk
    reg_ptr = transfer(reg, reg_ptr)
    chunk_ptr = transfer(chunk, chunk_ptr)
    call get_unigrid(repository=repository, reg=reg_ptr%p, lmax=lmax, symlog=symlog, chunk=chunk_ptr%p)
end subroutine f90wrap_get_unigrid

! End of module export_amr defined in file export_module.fpp

