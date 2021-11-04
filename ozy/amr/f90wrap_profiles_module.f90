! Module amr_profiles defined in file profiles_module.fpp

subroutine f90wrap_profile_handler__get__profdim(this, f90wrap_profdim)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_profdim
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_profdim = this_ptr%p%profdim
end subroutine f90wrap_profile_handler__get__profdim

subroutine f90wrap_profile_handler__set__profdim(this, f90wrap_profdim)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_profdim
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%profdim = f90wrap_profdim
end subroutine f90wrap_profile_handler__set__profdim

subroutine f90wrap_profile_handler__get__xvarname(this, f90wrap_xvarname)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
    character(128), intent(out) :: f90wrap_xvarname
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_xvarname = this_ptr%p%xvarname
end subroutine f90wrap_profile_handler__get__xvarname

subroutine f90wrap_profile_handler__set__xvarname(this, f90wrap_xvarname)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
    character(128), intent(in) :: f90wrap_xvarname
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%xvarname = f90wrap_xvarname
end subroutine f90wrap_profile_handler__set__xvarname

subroutine f90wrap_profile_handler__get__nyvar(this, f90wrap_nyvar)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nyvar
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nyvar = this_ptr%p%nyvar
end subroutine f90wrap_profile_handler__get__nyvar

subroutine f90wrap_profile_handler__set__nyvar(this, f90wrap_nyvar)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nyvar
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nyvar = f90wrap_nyvar
end subroutine f90wrap_profile_handler__set__nyvar

subroutine f90wrap_profile_handler__array__yvarnames(this, nd, dtype, dshape, dloc)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in) :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%yvarnames)) then
        dshape(1:2) = (/len(this_ptr%p%yvarnames(1)), shape(this_ptr%p%yvarnames)/)
        dloc = loc(this_ptr%p%yvarnames)
    else
        dloc = 0
    end if
end subroutine f90wrap_profile_handler__array__yvarnames

subroutine f90wrap_profile_handler__get__nbins(this, f90wrap_nbins)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nbins
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nbins = this_ptr%p%nbins
end subroutine f90wrap_profile_handler__get__nbins

subroutine f90wrap_profile_handler__set__nbins(this, f90wrap_nbins)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nbins
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nbins = f90wrap_nbins
end subroutine f90wrap_profile_handler__set__nbins

subroutine f90wrap_profile_handler__get__nwvar(this, f90wrap_nwvar)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nwvar
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nwvar = this_ptr%p%nwvar
end subroutine f90wrap_profile_handler__get__nwvar

subroutine f90wrap_profile_handler__set__nwvar(this, f90wrap_nwvar)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nwvar
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nwvar = f90wrap_nwvar
end subroutine f90wrap_profile_handler__set__nwvar

subroutine f90wrap_profile_handler__array__wvarnames(this, nd, dtype, dshape, dloc)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in) :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
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
end subroutine f90wrap_profile_handler__array__wvarnames

subroutine f90wrap_profile_handler__array__xdata(this, nd, dtype, dshape, dloc)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in) :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%xdata)) then
        dshape(1:1) = shape(this_ptr%p%xdata)
        dloc = loc(this_ptr%p%xdata)
    else
        dloc = 0
    end if
end subroutine f90wrap_profile_handler__array__xdata

subroutine f90wrap_profile_handler__array__ydata(this, nd, dtype, dshape, dloc)
    use amr_profiles, only: profile_handler
    implicit none
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    integer, intent(in) :: this(2)
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ydata)) then
        dshape(1:4) = shape(this_ptr%p%ydata)
        dloc = loc(this_ptr%p%ydata)
    else
        dloc = 0
    end if
end subroutine f90wrap_profile_handler__array__ydata

subroutine f90wrap_profile_handler_initialise(this)
    use amr_profiles, only: profile_handler
    implicit none
    
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_profile_handler_initialise

subroutine f90wrap_profile_handler_finalise(this)
    use amr_profiles, only: profile_handler
    implicit none
    
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_profile_handler_finalise

subroutine f90wrap_profile_handler_twod__get__profdim(this, f90wrap_profdim)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_profdim
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_profdim = this_ptr%p%profdim
end subroutine f90wrap_profile_handler_twod__get__profdim

subroutine f90wrap_profile_handler_twod__set__profdim(this, f90wrap_profdim)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_profdim
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%profdim = f90wrap_profdim
end subroutine f90wrap_profile_handler_twod__set__profdim

subroutine f90wrap_profile_handler_twod__get__xvarname(this, f90wrap_xvarname)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    character(128), intent(out) :: f90wrap_xvarname
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_xvarname = this_ptr%p%xvarname
end subroutine f90wrap_profile_handler_twod__get__xvarname

subroutine f90wrap_profile_handler_twod__set__xvarname(this, f90wrap_xvarname)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    character(128), intent(in) :: f90wrap_xvarname
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%xvarname = f90wrap_xvarname
end subroutine f90wrap_profile_handler_twod__set__xvarname

subroutine f90wrap_profile_handler_twod__get__yvarname(this, f90wrap_yvarname)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    character(128), intent(out) :: f90wrap_yvarname
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_yvarname = this_ptr%p%yvarname
end subroutine f90wrap_profile_handler_twod__get__yvarname

subroutine f90wrap_profile_handler_twod__set__yvarname(this, f90wrap_yvarname)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    character(128), intent(in) :: f90wrap_yvarname
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%yvarname = f90wrap_yvarname
end subroutine f90wrap_profile_handler_twod__set__yvarname

subroutine f90wrap_profile_handler_twod__get__nzvar(this, f90wrap_nzvar)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nzvar
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nzvar = this_ptr%p%nzvar
end subroutine f90wrap_profile_handler_twod__get__nzvar

subroutine f90wrap_profile_handler_twod__set__nzvar(this, f90wrap_nzvar)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nzvar
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nzvar = f90wrap_nzvar
end subroutine f90wrap_profile_handler_twod__set__nzvar

subroutine f90wrap_profile_handler_twod__array__zvarnames(this, nd, dtype, dshape, dloc)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in) :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%zvarnames)) then
        dshape(1:2) = (/len(this_ptr%p%zvarnames(1)), shape(this_ptr%p%zvarnames)/)
        dloc = loc(this_ptr%p%zvarnames)
    else
        dloc = 0
    end if
end subroutine f90wrap_profile_handler_twod__array__zvarnames

subroutine f90wrap_profile_handler_twod__array__nbins(this, nd, dtype, dshape, dloc)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in) :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%nbins)
    dloc = loc(this_ptr%p%nbins)
end subroutine f90wrap_profile_handler_twod__array__nbins

subroutine f90wrap_profile_handler_twod__get__nwvar(this, f90wrap_nwvar)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nwvar
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nwvar = this_ptr%p%nwvar
end subroutine f90wrap_profile_handler_twod__get__nwvar

subroutine f90wrap_profile_handler_twod__set__nwvar(this, f90wrap_nwvar)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in)   :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nwvar
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nwvar = f90wrap_nwvar
end subroutine f90wrap_profile_handler_twod__set__nwvar

subroutine f90wrap_profile_handler_twod__array__wvarnames(this, nd, dtype, dshape, dloc)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in) :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
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
end subroutine f90wrap_profile_handler_twod__array__wvarnames

subroutine f90wrap_profile_handler_twod__array__xdata(this, nd, dtype, dshape, dloc)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in) :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%xdata)) then
        dshape(1:1) = shape(this_ptr%p%xdata)
        dloc = loc(this_ptr%p%xdata)
    else
        dloc = 0
    end if
end subroutine f90wrap_profile_handler_twod__array__xdata

subroutine f90wrap_profile_handler_twod__array__ydata(this, nd, dtype, dshape, dloc)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in) :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ydata)) then
        dshape(1:1) = shape(this_ptr%p%ydata)
        dloc = loc(this_ptr%p%ydata)
    else
        dloc = 0
    end if
end subroutine f90wrap_profile_handler_twod__array__ydata

subroutine f90wrap_profile_handler_twod__array__zdata(this, nd, dtype, dshape, dloc)
    use amr_profiles, only: profile_handler_twod
    implicit none
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    integer, intent(in) :: this(2)
    type(profile_handler_twod_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 5
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%zdata)) then
        dshape(1:5) = shape(this_ptr%p%zdata)
        dloc = loc(this_ptr%p%zdata)
    else
        dloc = 0
    end if
end subroutine f90wrap_profile_handler_twod__array__zdata

subroutine f90wrap_profile_handler_twod_initialise(this)
    use amr_profiles, only: profile_handler_twod
    implicit none
    
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    type(profile_handler_twod_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_profile_handler_twod_initialise

subroutine f90wrap_profile_handler_twod_finalise(this)
    use amr_profiles, only: profile_handler_twod
    implicit none
    
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    type(profile_handler_twod_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_profile_handler_twod_finalise

subroutine f90wrap_allocate_profile_handler(prof)
    use amr_profiles, only: allocate_profile_handler, profile_handler
    implicit none
    
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type(profile_handler_ptr_type) :: prof_ptr
    integer, intent(in), dimension(2) :: prof
    prof_ptr = transfer(prof, prof_ptr)
    call allocate_profile_handler(prof=prof_ptr%p)
end subroutine f90wrap_allocate_profile_handler

subroutine f90wrap_allocate_profile_handler_twod(prof)
    use amr_profiles, only: profile_handler_twod, allocate_profile_handler_twod
    implicit none
    
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    type(profile_handler_twod_ptr_type) :: prof_ptr
    integer, intent(in), dimension(2) :: prof
    prof_ptr = transfer(prof, prof_ptr)
    call allocate_profile_handler_twod(prof=prof_ptr%p)
end subroutine f90wrap_allocate_profile_handler_twod

subroutine f90wrap_makebins(reg, sim, varname, nbins, bins, logscale, n0)
    use geometrical_regions, only: region
    use io_ramses, only: sim_info
    use amr_profiles, only: makebins
    implicit none
    
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(sim_info_ptr_type) :: sim_ptr
    integer, intent(in), dimension(2) :: sim
    character(128), intent(in) :: varname
    integer, intent(in) :: nbins
    real(8), intent(inout), dimension(n0) :: bins
    logical, intent(in) :: logscale
    integer :: n0
    !f2py intent(hide), depend(bins) :: n0 = shape(bins,0)
    reg_ptr = transfer(reg, reg_ptr)
    sim_ptr = transfer(sim, sim_ptr)
    call makebins(reg=reg_ptr%p, sim=sim_ptr%p, varname=varname, nbins=nbins, bins=bins, logscale=logscale)
end subroutine f90wrap_makebins

subroutine f90wrap_findbinpos(reg, varids, distance, pos, cellvars, cellsize, prof, ibin, n0, n1)
    use io_ramses, only: hydroid
    use amr_profiles, only: profile_handler, findbinpos
    use geometrical_regions, only: region
    implicit none
    
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    real(8), intent(in) :: distance
    real(8), intent(in), dimension(n0) :: pos
    real(8), intent(in), dimension(n1) :: cellvars
    real(8), intent(in) :: cellsize
    type(profile_handler_ptr_type) :: prof_ptr
    integer, intent(in), dimension(2) :: prof
    integer, intent(inout) :: ibin
    integer :: n0
    !f2py intent(hide), depend(pos) :: n0 = shape(pos,0)
    integer :: n1
    !f2py intent(hide), depend(cellvars) :: n1 = shape(cellvars,0)
    reg_ptr = transfer(reg, reg_ptr)
    varids_ptr = transfer(varids, varids_ptr)
    prof_ptr = transfer(prof, prof_ptr)
    call findbinpos(reg=reg_ptr%p, varIDs=varids_ptr%p, distance=distance, pos=pos, cellvars=cellvars, cellsize=cellsize, &
        prof=prof_ptr%p, ibin=ibin)
end subroutine f90wrap_findbinpos

subroutine f90wrap_findbinpos_twod(reg, varids, distance, pos, cellvars, cellsize, prof, logscale, ibinx, ibiny, n0, n1)
    use io_ramses, only: hydroid
    use geometrical_regions, only: region
    use amr_profiles, only: findbinpos_twod, profile_handler_twod
    implicit none
    
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    real(8), intent(in) :: distance
    real(8), intent(in), dimension(n0) :: pos
    real(8), intent(in), dimension(n1) :: cellvars
    real(8), intent(in) :: cellsize
    type(profile_handler_twod_ptr_type) :: prof_ptr
    integer, intent(in), dimension(2) :: prof
    logical, intent(in) :: logscale
    integer, intent(inout) :: ibinx
    integer, intent(inout) :: ibiny
    integer :: n0
    !f2py intent(hide), depend(pos) :: n0 = shape(pos,0)
    integer :: n1
    !f2py intent(hide), depend(cellvars) :: n1 = shape(cellvars,0)
    reg_ptr = transfer(reg, reg_ptr)
    varids_ptr = transfer(varids, varids_ptr)
    prof_ptr = transfer(prof, prof_ptr)
    call findbinpos_twod(reg=reg_ptr%p, varIDs=varids_ptr%p, distance=distance, pos=pos, cellvars=cellvars, &
        cellsize=cellsize, prof=prof_ptr%p, logscale=logscale, ibinx=ibinx, ibiny=ibiny)
end subroutine f90wrap_findbinpos_twod

subroutine f90wrap_bindata(reg, varids, pos, cellvars, cellsize, prof, ibin, n0, n1)
    use io_ramses, only: hydroid
    use amr_profiles, only: bindata, profile_handler
    use geometrical_regions, only: region
    implicit none
    
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    real(8), intent(in), dimension(n0) :: pos
    real(8), intent(in), dimension(n1) :: cellvars
    real(8), intent(in) :: cellsize
    type(profile_handler_ptr_type) :: prof_ptr
    integer, intent(in), dimension(2) :: prof
    integer, intent(in) :: ibin
    integer :: n0
    !f2py intent(hide), depend(pos) :: n0 = shape(pos,0)
    integer :: n1
    !f2py intent(hide), depend(cellvars) :: n1 = shape(cellvars,0)
    reg_ptr = transfer(reg, reg_ptr)
    varids_ptr = transfer(varids, varids_ptr)
    prof_ptr = transfer(prof, prof_ptr)
    call bindata(reg=reg_ptr%p, varIDs=varids_ptr%p, pos=pos, cellvars=cellvars, cellsize=cellsize, prof=prof_ptr%p, &
        ibin=ibin)
end subroutine f90wrap_bindata

subroutine f90wrap_bindata_twod(reg, varids, pos, cellvars, cellsize, prof, ibinx, ibiny, n0, n1)
    use io_ramses, only: hydroid
    use geometrical_regions, only: region
    use amr_profiles, only: profile_handler_twod, bindata_twod
    implicit none
    
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    real(8), intent(in), dimension(n0) :: pos
    real(8), intent(in), dimension(n1) :: cellvars
    real(8), intent(in) :: cellsize
    type(profile_handler_twod_ptr_type) :: prof_ptr
    integer, intent(in), dimension(2) :: prof
    integer, intent(in) :: ibinx
    integer, intent(in) :: ibiny
    integer :: n0
    !f2py intent(hide), depend(pos) :: n0 = shape(pos,0)
    integer :: n1
    !f2py intent(hide), depend(cellvars) :: n1 = shape(cellvars,0)
    reg_ptr = transfer(reg, reg_ptr)
    varids_ptr = transfer(varids, varids_ptr)
    prof_ptr = transfer(prof, prof_ptr)
    call bindata_twod(reg=reg_ptr%p, varIDs=varids_ptr%p, pos=pos, cellvars=cellvars, cellsize=cellsize, prof=prof_ptr%p, &
        ibinx=ibinx, ibiny=ibiny)
end subroutine f90wrap_bindata_twod

subroutine f90wrap_renormalise_bins(prof_data)
    use amr_profiles, only: renormalise_bins, profile_handler
    implicit none
    
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type(profile_handler_ptr_type) :: prof_data_ptr
    integer, intent(in), dimension(2) :: prof_data
    prof_data_ptr = transfer(prof_data, prof_data_ptr)
    call renormalise_bins(prof_data=prof_data_ptr%p)
end subroutine f90wrap_renormalise_bins

subroutine f90wrap_renormalise_bins_twod(prof_data)
    use amr_profiles, only: renormalise_bins_twod, profile_handler_twod
    implicit none
    
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    type(profile_handler_twod_ptr_type) :: prof_data_ptr
    integer, intent(in), dimension(2) :: prof_data
    prof_data_ptr = transfer(prof_data, prof_data_ptr)
    call renormalise_bins_twod(prof_data=prof_data_ptr%p)
end subroutine f90wrap_renormalise_bins_twod

subroutine f90wrap_get_cells_onedprofile(repository, amr, reg, filt, varids, prof_data)
    use amr_profiles, only: get_cells_onedprofile, profile_handler
    use io_ramses, only: hydroid, amr_info
    use filtering, only: filter
    use geometrical_regions, only: region
    implicit none
    
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    character(128), intent(in) :: repository
    type(amr_info_ptr_type) :: amr_ptr
    integer, intent(in), dimension(2) :: amr
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(filter_ptr_type) :: filt_ptr
    integer, intent(in), dimension(2) :: filt
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    type(profile_handler_ptr_type) :: prof_data_ptr
    integer, intent(in), dimension(2) :: prof_data
    amr_ptr = transfer(amr, amr_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    filt_ptr = transfer(filt, filt_ptr)
    varids_ptr = transfer(varids, varids_ptr)
    prof_data_ptr = transfer(prof_data, prof_data_ptr)
    call get_cells_onedprofile(repository=repository, amr=amr_ptr%p, reg=reg_ptr%p, filt=filt_ptr%p, varIDs=varids_ptr%p, &
        prof_data=prof_data_ptr%p)
end subroutine f90wrap_get_cells_onedprofile

subroutine f90wrap_onedprofile(repository, reg, filt, prof_data, lmax, logscale)
    use filtering, only: filter
    use amr_profiles, only: profile_handler, onedprofile
    use geometrical_regions, only: region
    implicit none
    
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    character(128), intent(in) :: repository
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(filter_ptr_type) :: filt_ptr
    integer, intent(in), dimension(2) :: filt
    type(profile_handler_ptr_type) :: prof_data_ptr
    integer, intent(in), dimension(2) :: prof_data
    integer, intent(in) :: lmax
    logical, intent(in) :: logscale
    reg_ptr = transfer(reg, reg_ptr)
    filt_ptr = transfer(filt, filt_ptr)
    prof_data_ptr = transfer(prof_data, prof_data_ptr)
    call onedprofile(repository=repository, reg=reg_ptr%p, filt=filt_ptr%p, prof_data=prof_data_ptr%p, lmax=lmax, &
        logscale=logscale)
end subroutine f90wrap_onedprofile

subroutine f90wrap_twodprofile(repository, reg, filt, prof_data, lmax, logscale)
    use filtering, only: filter
    use geometrical_regions, only: region
    use amr_profiles, only: twodprofile, profile_handler_twod
    implicit none
    
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    character(128), intent(in) :: repository
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(filter_ptr_type) :: filt_ptr
    integer, intent(in), dimension(2) :: filt
    type(profile_handler_twod_ptr_type) :: prof_data_ptr
    integer, intent(in), dimension(2) :: prof_data
    integer, intent(in) :: lmax
    logical, intent(in) :: logscale
    reg_ptr = transfer(reg, reg_ptr)
    filt_ptr = transfer(filt, filt_ptr)
    prof_data_ptr = transfer(prof_data, prof_data_ptr)
    call twodprofile(repository=repository, reg=reg_ptr%p, filt=filt_ptr%p, prof_data=prof_data_ptr%p, lmax=lmax, &
        logscale=logscale)
end subroutine f90wrap_twodprofile

subroutine f90wrap_get_cells_twodprofile(repository, amr, reg, filt, varids, prof_data, logscale)
    use amr_profiles, only: profile_handler_twod, get_cells_twodprofile
    use io_ramses, only: hydroid, amr_info
    use filtering, only: filter
    use geometrical_regions, only: region
    implicit none
    
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type profile_handler_twod_ptr_type
        type(profile_handler_twod), pointer :: p => NULL()
    end type profile_handler_twod_ptr_type
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    character(128), intent(in) :: repository
    type(amr_info_ptr_type) :: amr_ptr
    integer, intent(in), dimension(2) :: amr
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(filter_ptr_type) :: filt_ptr
    integer, intent(in), dimension(2) :: filt
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    type(profile_handler_twod_ptr_type) :: prof_data_ptr
    integer, intent(in), dimension(2) :: prof_data
    logical, intent(in) :: logscale
    amr_ptr = transfer(amr, amr_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    filt_ptr = transfer(filt, filt_ptr)
    varids_ptr = transfer(varids, varids_ptr)
    prof_data_ptr = transfer(prof_data, prof_data_ptr)
    call get_cells_twodprofile(repository=repository, amr=amr_ptr%p, reg=reg_ptr%p, filt=filt_ptr%p, varIDs=varids_ptr%p, &
        prof_data=prof_data_ptr%p, logscale=logscale)
end subroutine f90wrap_get_cells_twodprofile

! End of module amr_profiles defined in file profiles_module.fpp

