! Module part_profiles defined in file profiles_module.fpp

subroutine f90wrap_profile_handler__get__profdim(this, f90wrap_profdim)
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
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
    use part_profiles, only: profile_handler
    implicit none
    
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type(profile_handler_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_profile_handler_finalise

subroutine f90wrap_allocate_profile_handler(prof)
    use part_profiles, only: allocate_profile_handler, profile_handler
    implicit none
    
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type(profile_handler_ptr_type) :: prof_ptr
    integer, intent(in), dimension(2) :: prof
    prof_ptr = transfer(prof, prof_ptr)
    call allocate_profile_handler(prof=prof_ptr%p)
end subroutine f90wrap_allocate_profile_handler

subroutine f90wrap_makebins(reg, varname, nbins, bins, n0)
    use geometrical_regions, only: region
    use part_profiles, only: makebins
    implicit none
    
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    character(128), intent(in) :: varname
    integer, intent(in) :: nbins
    real(8), dimension(n0) :: bins
    integer :: n0
    !f2py intent(hide), depend(bins) :: n0 = shape(bins,0)
    reg_ptr = transfer(reg, reg_ptr)
    call makebins(reg=reg_ptr%p, varname=varname, nbins=nbins, bins=bins)
end subroutine f90wrap_makebins

subroutine f90wrap_findbinpos(sim, reg, distance, part, prof, ibin)
    use part_profiles, only: findbinpos, profile_handler
    use io_ramses, only: particle, sim_info
    use geometrical_regions, only: region
    implicit none
    
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type(sim_info_ptr_type) :: sim_ptr
    integer, intent(in), dimension(2) :: sim
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    real(8), intent(in) :: distance
    type(particle_ptr_type) :: part_ptr
    integer, intent(in), dimension(2) :: part
    type(profile_handler_ptr_type) :: prof_ptr
    integer, intent(in), dimension(2) :: prof
    integer, intent(inout) :: ibin
    sim_ptr = transfer(sim, sim_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    part_ptr = transfer(part, part_ptr)
    prof_ptr = transfer(prof, prof_ptr)
    call findbinpos(sim=sim_ptr%p, reg=reg_ptr%p, distance=distance, part=part_ptr%p, prof=prof_ptr%p, ibin=ibin)
end subroutine f90wrap_findbinpos

subroutine f90wrap_bindata(sim, reg, part, prof, ibin)
    use io_ramses, only: particle, sim_info
    use part_profiles, only: bindata, profile_handler
    use geometrical_regions, only: region
    implicit none
    
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type(sim_info_ptr_type) :: sim_ptr
    integer, intent(in), dimension(2) :: sim
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(particle_ptr_type) :: part_ptr
    integer, intent(in), dimension(2) :: part
    type(profile_handler_ptr_type) :: prof_ptr
    integer, intent(in), dimension(2) :: prof
    integer, intent(in) :: ibin
    sim_ptr = transfer(sim, sim_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    part_ptr = transfer(part, part_ptr)
    prof_ptr = transfer(prof, prof_ptr)
    call bindata(sim=sim_ptr%p, reg=reg_ptr%p, part=part_ptr%p, prof=prof_ptr%p, ibin=ibin)
end subroutine f90wrap_bindata

subroutine f90wrap_renormalise_bins(prof_data)
    use part_profiles, only: renormalise_bins, profile_handler
    implicit none
    
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type(profile_handler_ptr_type) :: prof_data_ptr
    integer, intent(in), dimension(2) :: prof_data
    prof_data_ptr = transfer(prof_data, prof_data_ptr)
    call renormalise_bins(prof_data=prof_data_ptr%p)
end subroutine f90wrap_renormalise_bins

subroutine f90wrap_get_parts_onedprofile(repository, amr, sim, reg, filt, prof_data)
    use filtering, only: filter
    use io_ramses, only: amr_info, sim_info
    use part_profiles, only: get_parts_onedprofile, profile_handler
    use geometrical_regions, only: region
    implicit none
    
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    character(128), intent(in) :: repository
    type(amr_info_ptr_type) :: amr_ptr
    integer, intent(in), dimension(2) :: amr
    type(sim_info_ptr_type) :: sim_ptr
    integer, intent(in), dimension(2) :: sim
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(filter_ptr_type) :: filt_ptr
    integer, intent(in), dimension(2) :: filt
    type(profile_handler_ptr_type) :: prof_data_ptr
    integer, intent(in), dimension(2) :: prof_data
    amr_ptr = transfer(amr, amr_ptr)
    sim_ptr = transfer(sim, sim_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    filt_ptr = transfer(filt, filt_ptr)
    prof_data_ptr = transfer(prof_data, prof_data_ptr)
    call get_parts_onedprofile(repository=repository, amr=amr_ptr%p, sim=sim_ptr%p, reg=reg_ptr%p, filt=filt_ptr%p, &
        prof_data=prof_data_ptr%p)
end subroutine f90wrap_get_parts_onedprofile

subroutine f90wrap_onedprofile(repository, reg, filt, prof_data, lmax)
    use geometrical_regions, only: region
    use part_profiles, only: onedprofile, profile_handler
    use filtering, only: filter
    implicit none
    
    type profile_handler_ptr_type
        type(profile_handler), pointer :: p => NULL()
    end type profile_handler_ptr_type
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
    type(profile_handler_ptr_type) :: prof_data_ptr
    integer, intent(in), dimension(2) :: prof_data
    integer, intent(in) :: lmax
    reg_ptr = transfer(reg, reg_ptr)
    filt_ptr = transfer(filt, filt_ptr)
    prof_data_ptr = transfer(prof_data, prof_data_ptr)
    call onedprofile(repository=repository, reg=reg_ptr%p, filt=filt_ptr%p, prof_data=prof_data_ptr%p, lmax=lmax)
end subroutine f90wrap_onedprofile

! End of module part_profiles defined in file profiles_module.fpp

