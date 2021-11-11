! Module part_integrator defined in file integrator_module.fpp

subroutine f90wrap_part_region_attrs__get__nvars(this, f90wrap_nvars)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in)   :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nvars
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nvars = this_ptr%p%nvars
end subroutine f90wrap_part_region_attrs__get__nvars

subroutine f90wrap_part_region_attrs__set__nvars(this, f90wrap_nvars)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in)   :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nvars
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nvars = f90wrap_nvars
end subroutine f90wrap_part_region_attrs__set__nvars

subroutine f90wrap_part_region_attrs__array__varnames(this, nd, dtype, dshape, dloc)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in) :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
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
end subroutine f90wrap_part_region_attrs__array__varnames

subroutine f90wrap_part_region_attrs__get__nwvars(this, f90wrap_nwvars)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in)   :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nwvars
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nwvars = this_ptr%p%nwvars
end subroutine f90wrap_part_region_attrs__get__nwvars

subroutine f90wrap_part_region_attrs__set__nwvars(this, f90wrap_nwvars)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in)   :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nwvars
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nwvars = f90wrap_nwvars
end subroutine f90wrap_part_region_attrs__set__nwvars

subroutine f90wrap_part_region_attrs__array__wvarnames(this, nd, dtype, dshape, dloc)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in) :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%wvarnames)) then
        dshape(1:2) = (/len(this_ptr%p%wvarnames(1)), shape(this_ptr%p%wvarnames)/)
        dloc = loc(this_ptr%p%wvarnames)
    else
        dloc = 0
    end if
end subroutine f90wrap_part_region_attrs__array__wvarnames

subroutine f90wrap_part_region_attrs__get__ndm(this, f90wrap_ndm)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in)   :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_ndm
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ndm = this_ptr%p%ndm
end subroutine f90wrap_part_region_attrs__get__ndm

subroutine f90wrap_part_region_attrs__set__ndm(this, f90wrap_ndm)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in)   :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_ndm
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ndm = f90wrap_ndm
end subroutine f90wrap_part_region_attrs__set__ndm

subroutine f90wrap_part_region_attrs__get__nstar(this, f90wrap_nstar)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in)   :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nstar
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nstar = this_ptr%p%nstar
end subroutine f90wrap_part_region_attrs__get__nstar

subroutine f90wrap_part_region_attrs__set__nstar(this, f90wrap_nstar)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in)   :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nstar
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nstar = f90wrap_nstar
end subroutine f90wrap_part_region_attrs__set__nstar

subroutine f90wrap_part_region_attrs__get__nids(this, f90wrap_nids)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in)   :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nids
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nids = this_ptr%p%nids
end subroutine f90wrap_part_region_attrs__get__nids

subroutine f90wrap_part_region_attrs__set__nids(this, f90wrap_nids)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in)   :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nids
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nids = f90wrap_nids
end subroutine f90wrap_part_region_attrs__set__nids

subroutine f90wrap_part_region_attrs__array__data(this, nd, dtype, dshape, dloc)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in) :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%data)) then
        dshape(1:3) = shape(this_ptr%p%data)
        dloc = loc(this_ptr%p%data)
    else
        dloc = 0
    end if
end subroutine f90wrap_part_region_attrs__array__data

subroutine f90wrap_part_region_attrs__array__ids(this, nd, dtype, dshape, dloc)
    use part_integrator, only: part_region_attrs
    implicit none
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    integer, intent(in) :: this(2)
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ids)) then
        dshape(1:1) = shape(this_ptr%p%ids)
        dloc = loc(this_ptr%p%ids)
    else
        dloc = 0
    end if
end subroutine f90wrap_part_region_attrs__array__ids

subroutine f90wrap_part_region_attrs_initialise(this)
    use part_integrator, only: part_region_attrs
    implicit none
    
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_part_region_attrs_initialise

subroutine f90wrap_part_region_attrs_finalise(this)
    use part_integrator, only: part_region_attrs
    implicit none
    
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    type(part_region_attrs_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_part_region_attrs_finalise

subroutine f90wrap_allocate_part_regions_attrs(attrs)
    use part_integrator, only: allocate_part_regions_attrs, part_region_attrs
    implicit none
    
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    type(part_region_attrs_ptr_type) :: attrs_ptr
    integer, intent(in), dimension(2) :: attrs
    attrs_ptr = transfer(attrs, attrs_ptr)
    call allocate_part_regions_attrs(attrs=attrs_ptr%p)
end subroutine f90wrap_allocate_part_regions_attrs

subroutine f90wrap_extract_data(sim, reg, part, attrs)
    use io_ramses, only: particle, sim_info
    use geometrical_regions, only: region
    use part_integrator, only: part_region_attrs, extract_data
    implicit none
    
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    type(sim_info_ptr_type) :: sim_ptr
    integer, intent(in), dimension(2) :: sim
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(particle_ptr_type) :: part_ptr
    integer, intent(in), dimension(2) :: part
    type(part_region_attrs_ptr_type) :: attrs_ptr
    integer, intent(in), dimension(2) :: attrs
    sim_ptr = transfer(sim, sim_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    part_ptr = transfer(part, part_ptr)
    attrs_ptr = transfer(attrs, attrs_ptr)
    call extract_data(sim=sim_ptr%p, reg=reg_ptr%p, part=part_ptr%p, attrs=attrs_ptr%p)
end subroutine f90wrap_extract_data

subroutine f90wrap_renormalise(sim, attrs)
    use part_integrator, only: part_region_attrs, renormalise
    use io_ramses, only: sim_info
    implicit none
    
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    type(sim_info_ptr_type) :: sim_ptr
    integer, intent(in), dimension(2) :: sim
    type(part_region_attrs_ptr_type) :: attrs_ptr
    integer, intent(in), dimension(2) :: attrs
    sim_ptr = transfer(sim, sim_ptr)
    attrs_ptr = transfer(attrs, attrs_ptr)
    call renormalise(sim=sim_ptr%p, attrs=attrs_ptr%p)
end subroutine f90wrap_renormalise

subroutine f90wrap_integrate_region(repository, reg, filt, attrs, get_ids)
    use part_integrator, only: integrate_region, part_region_attrs
    use filtering, only: filter
    use geometrical_regions, only: region
    implicit none
    
    type part_region_attrs_ptr_type
        type(part_region_attrs), pointer :: p => NULL()
    end type part_region_attrs_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    character(128), intent(in) :: repository
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(filter_ptr_type) :: filt_ptr
    integer, intent(in), dimension(2) :: filt
    type(part_region_attrs_ptr_type) :: attrs_ptr
    integer, intent(in), dimension(2) :: attrs
    logical, intent(in), optional :: get_ids
    reg_ptr = transfer(reg, reg_ptr)
    filt_ptr = transfer(filt, filt_ptr)
    attrs_ptr = transfer(attrs, attrs_ptr)
    call integrate_region(repository=repository, reg=reg_ptr%p, filt=filt_ptr%p, attrs=attrs_ptr%p, get_ids=get_ids)
end subroutine f90wrap_integrate_region

! End of module part_integrator defined in file integrator_module.fpp

