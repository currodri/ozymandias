! Module io_ramses defined in file read_amr_module.fpp

subroutine f90wrap_hydroid__get__nvar(this, f90wrap_nvar)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nvar
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nvar = this_ptr%p%nvar
end subroutine f90wrap_hydroid__get__nvar

subroutine f90wrap_hydroid__set__nvar(this, f90wrap_nvar)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nvar
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nvar = f90wrap_nvar
end subroutine f90wrap_hydroid__set__nvar

subroutine f90wrap_hydroid__get__density(this, f90wrap_density)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_density
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_density = this_ptr%p%density
end subroutine f90wrap_hydroid__get__density

subroutine f90wrap_hydroid__set__density(this, f90wrap_density)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_density
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%density = f90wrap_density
end subroutine f90wrap_hydroid__set__density

subroutine f90wrap_hydroid__get__vx(this, f90wrap_vx)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_vx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_vx = this_ptr%p%vx
end subroutine f90wrap_hydroid__get__vx

subroutine f90wrap_hydroid__set__vx(this, f90wrap_vx)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_vx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%vx = f90wrap_vx
end subroutine f90wrap_hydroid__set__vx

subroutine f90wrap_hydroid__get__vy(this, f90wrap_vy)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_vy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_vy = this_ptr%p%vy
end subroutine f90wrap_hydroid__get__vy

subroutine f90wrap_hydroid__set__vy(this, f90wrap_vy)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_vy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%vy = f90wrap_vy
end subroutine f90wrap_hydroid__set__vy

subroutine f90wrap_hydroid__get__vz(this, f90wrap_vz)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_vz
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_vz = this_ptr%p%vz
end subroutine f90wrap_hydroid__get__vz

subroutine f90wrap_hydroid__set__vz(this, f90wrap_vz)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_vz
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%vz = f90wrap_vz
end subroutine f90wrap_hydroid__set__vz

subroutine f90wrap_hydroid__get__thermal_pressure(this, f90wrap_thermal_pressure)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_thermal_pressure
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_thermal_pressure = this_ptr%p%thermal_pressure
end subroutine f90wrap_hydroid__get__thermal_pressure

subroutine f90wrap_hydroid__set__thermal_pressure(this, f90wrap_thermal_pressure)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_thermal_pressure
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%thermal_pressure = f90wrap_thermal_pressure
end subroutine f90wrap_hydroid__set__thermal_pressure

subroutine f90wrap_hydroid__get__metallicity(this, f90wrap_metallicity)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_metallicity
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_metallicity = this_ptr%p%metallicity
end subroutine f90wrap_hydroid__get__metallicity

subroutine f90wrap_hydroid__set__metallicity(this, f90wrap_metallicity)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_metallicity
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%metallicity = f90wrap_metallicity
end subroutine f90wrap_hydroid__set__metallicity

subroutine f90wrap_hydroid__get__Blx(this, f90wrap_Blx)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_Blx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Blx = this_ptr%p%Blx
end subroutine f90wrap_hydroid__get__Blx

subroutine f90wrap_hydroid__set__Blx(this, f90wrap_Blx)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_Blx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Blx = f90wrap_Blx
end subroutine f90wrap_hydroid__set__Blx

subroutine f90wrap_hydroid__get__Bly(this, f90wrap_Bly)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_Bly
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Bly = this_ptr%p%Bly
end subroutine f90wrap_hydroid__get__Bly

subroutine f90wrap_hydroid__set__Bly(this, f90wrap_Bly)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_Bly
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Bly = f90wrap_Bly
end subroutine f90wrap_hydroid__set__Bly

subroutine f90wrap_hydroid__get__Blz(this, f90wrap_Blz)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_Blz
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Blz = this_ptr%p%Blz
end subroutine f90wrap_hydroid__get__Blz

subroutine f90wrap_hydroid__set__Blz(this, f90wrap_Blz)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_Blz
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Blz = f90wrap_Blz
end subroutine f90wrap_hydroid__set__Blz

subroutine f90wrap_hydroid__get__Brx(this, f90wrap_Brx)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_Brx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Brx = this_ptr%p%Brx
end subroutine f90wrap_hydroid__get__Brx

subroutine f90wrap_hydroid__set__Brx(this, f90wrap_Brx)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_Brx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Brx = f90wrap_Brx
end subroutine f90wrap_hydroid__set__Brx

subroutine f90wrap_hydroid__get__Bry(this, f90wrap_Bry)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_Bry
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Bry = this_ptr%p%Bry
end subroutine f90wrap_hydroid__get__Bry

subroutine f90wrap_hydroid__set__Bry(this, f90wrap_Bry)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_Bry
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Bry = f90wrap_Bry
end subroutine f90wrap_hydroid__set__Bry

subroutine f90wrap_hydroid__get__Brz(this, f90wrap_Brz)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_Brz
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Brz = this_ptr%p%Brz
end subroutine f90wrap_hydroid__get__Brz

subroutine f90wrap_hydroid__set__Brz(this, f90wrap_Brz)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_Brz
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Brz = f90wrap_Brz
end subroutine f90wrap_hydroid__set__Brz

subroutine f90wrap_hydroid__get__eCR(this, f90wrap_eCR)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_eCR
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_eCR = this_ptr%p%eCR
end subroutine f90wrap_hydroid__get__eCR

subroutine f90wrap_hydroid__set__eCR(this, f90wrap_eCR)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_eCR
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eCR = f90wrap_eCR
end subroutine f90wrap_hydroid__set__eCR

subroutine f90wrap_hydroid__get__xHII(this, f90wrap_xHII)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_xHII
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_xHII = this_ptr%p%xHII
end subroutine f90wrap_hydroid__get__xHII

subroutine f90wrap_hydroid__set__xHII(this, f90wrap_xHII)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_xHII
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%xHII = f90wrap_xHII
end subroutine f90wrap_hydroid__set__xHII

subroutine f90wrap_hydroid__get__xHeII(this, f90wrap_xHeII)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_xHeII
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_xHeII = this_ptr%p%xHeII
end subroutine f90wrap_hydroid__get__xHeII

subroutine f90wrap_hydroid__set__xHeII(this, f90wrap_xHeII)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_xHeII
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%xHeII = f90wrap_xHeII
end subroutine f90wrap_hydroid__set__xHeII

subroutine f90wrap_hydroid__get__xHeIII(this, f90wrap_xHeIII)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_xHeIII
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_xHeIII = this_ptr%p%xHeIII
end subroutine f90wrap_hydroid__get__xHeIII

subroutine f90wrap_hydroid__set__xHeIII(this, f90wrap_xHeIII)
    use io_ramses, only: hydroid
    implicit none
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    integer, intent(in)   :: this(2)
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_xHeIII
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%xHeIII = f90wrap_xHeIII
end subroutine f90wrap_hydroid__set__xHeIII

subroutine f90wrap_hydroid_initialise(this)
    use io_ramses, only: hydroid
    implicit none
    
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_hydroid_initialise

subroutine f90wrap_hydroid_finalise(this)
    use io_ramses, only: hydroid
    implicit none
    
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    type(hydroid_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_hydroid_finalise

subroutine f90wrap_amr_info__get__ncpu(this, f90wrap_ncpu)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_ncpu
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ncpu = this_ptr%p%ncpu
end subroutine f90wrap_amr_info__get__ncpu

subroutine f90wrap_amr_info__set__ncpu(this, f90wrap_ncpu)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_ncpu
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ncpu = f90wrap_ncpu
end subroutine f90wrap_amr_info__set__ncpu

subroutine f90wrap_amr_info__get__ndim(this, f90wrap_ndim)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_ndim
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ndim = this_ptr%p%ndim
end subroutine f90wrap_amr_info__get__ndim

subroutine f90wrap_amr_info__set__ndim(this, f90wrap_ndim)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_ndim
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ndim = f90wrap_ndim
end subroutine f90wrap_amr_info__set__ndim

subroutine f90wrap_amr_info__get__nlevelmax(this, f90wrap_nlevelmax)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nlevelmax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nlevelmax = this_ptr%p%nlevelmax
end subroutine f90wrap_amr_info__get__nlevelmax

subroutine f90wrap_amr_info__set__nlevelmax(this, f90wrap_nlevelmax)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nlevelmax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nlevelmax = f90wrap_nlevelmax
end subroutine f90wrap_amr_info__set__nlevelmax

subroutine f90wrap_amr_info__get__nboundary(this, f90wrap_nboundary)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nboundary
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nboundary = this_ptr%p%nboundary
end subroutine f90wrap_amr_info__get__nboundary

subroutine f90wrap_amr_info__set__nboundary(this, f90wrap_nboundary)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nboundary
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nboundary = f90wrap_nboundary
end subroutine f90wrap_amr_info__set__nboundary

subroutine f90wrap_amr_info__get__twotondim(this, f90wrap_twotondim)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_twotondim
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_twotondim = this_ptr%p%twotondim
end subroutine f90wrap_amr_info__get__twotondim

subroutine f90wrap_amr_info__set__twotondim(this, f90wrap_twotondim)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_twotondim
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%twotondim = f90wrap_twotondim
end subroutine f90wrap_amr_info__set__twotondim

subroutine f90wrap_amr_info__get__ndom(this, f90wrap_ndom)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_ndom
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ndom = this_ptr%p%ndom
end subroutine f90wrap_amr_info__get__ndom

subroutine f90wrap_amr_info__set__ndom(this, f90wrap_ndom)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_ndom
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ndom = f90wrap_ndom
end subroutine f90wrap_amr_info__set__ndom

subroutine f90wrap_amr_info__get__levelmin(this, f90wrap_levelmin)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_levelmin
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_levelmin = this_ptr%p%levelmin
end subroutine f90wrap_amr_info__get__levelmin

subroutine f90wrap_amr_info__set__levelmin(this, f90wrap_levelmin)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_levelmin
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%levelmin = f90wrap_levelmin
end subroutine f90wrap_amr_info__set__levelmin

subroutine f90wrap_amr_info__get__levelmax(this, f90wrap_levelmax)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_levelmax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_levelmax = this_ptr%p%levelmax
end subroutine f90wrap_amr_info__get__levelmax

subroutine f90wrap_amr_info__set__levelmax(this, f90wrap_levelmax)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_levelmax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%levelmax = f90wrap_levelmax
end subroutine f90wrap_amr_info__set__levelmax

subroutine f90wrap_amr_info__get__lmax(this, f90wrap_lmax)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_lmax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_lmax = this_ptr%p%lmax
end subroutine f90wrap_amr_info__get__lmax

subroutine f90wrap_amr_info__set__lmax(this, f90wrap_lmax)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_lmax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lmax = f90wrap_lmax
end subroutine f90wrap_amr_info__set__lmax

subroutine f90wrap_amr_info__get__lmin(this, f90wrap_lmin)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_lmin
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_lmin = this_ptr%p%lmin
end subroutine f90wrap_amr_info__get__lmin

subroutine f90wrap_amr_info__set__lmin(this, f90wrap_lmin)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_lmin
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lmin = f90wrap_lmin
end subroutine f90wrap_amr_info__set__lmin

subroutine f90wrap_amr_info__get__ncpu_read(this, f90wrap_ncpu_read)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_ncpu_read
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ncpu_read = this_ptr%p%ncpu_read
end subroutine f90wrap_amr_info__get__ncpu_read

subroutine f90wrap_amr_info__set__ncpu_read(this, f90wrap_ncpu_read)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_ncpu_read
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ncpu_read = f90wrap_ncpu_read
end subroutine f90wrap_amr_info__set__ncpu_read

subroutine f90wrap_amr_info__get__ordering(this, f90wrap_ordering)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    character(80), intent(out) :: f90wrap_ordering
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ordering = this_ptr%p%ordering
end subroutine f90wrap_amr_info__get__ordering

subroutine f90wrap_amr_info__set__ordering(this, f90wrap_ordering)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in)   :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    character(80), intent(in) :: f90wrap_ordering
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ordering = f90wrap_ordering
end subroutine f90wrap_amr_info__set__ordering

subroutine f90wrap_amr_info__array__cpu_list(this, nd, dtype, dshape, dloc)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in) :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%cpu_list)) then
        dshape(1:1) = shape(this_ptr%p%cpu_list)
        dloc = loc(this_ptr%p%cpu_list)
    else
        dloc = 0
    end if
end subroutine f90wrap_amr_info__array__cpu_list

subroutine f90wrap_amr_info__array__bound_key(this, nd, dtype, dshape, dloc)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in) :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%bound_key)) then
        dshape(1:1) = shape(this_ptr%p%bound_key)
        dloc = loc(this_ptr%p%bound_key)
    else
        dloc = 0
    end if
end subroutine f90wrap_amr_info__array__bound_key

subroutine f90wrap_amr_info__array__cpu_read(this, nd, dtype, dshape, dloc)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in) :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%cpu_read)) then
        dshape(1:1) = shape(this_ptr%p%cpu_read)
        dloc = loc(this_ptr%p%cpu_read)
    else
        dloc = 0
    end if
end subroutine f90wrap_amr_info__array__cpu_read

subroutine f90wrap_amr_info__array__xbound(this, nd, dtype, dshape, dloc)
    use io_ramses, only: amr_info
    implicit none
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in) :: this(2)
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%xbound)
    dloc = loc(this_ptr%p%xbound)
end subroutine f90wrap_amr_info__array__xbound

subroutine f90wrap_amr_info_initialise(this)
    use io_ramses, only: amr_info
    implicit none
    
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_amr_info_initialise

subroutine f90wrap_amr_info_finalise(this)
    use io_ramses, only: amr_info
    implicit none
    
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    type(amr_info_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_amr_info_finalise

subroutine f90wrap_sim_info__get__cosmo(this, f90wrap_cosmo)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_cosmo
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_cosmo = this_ptr%p%cosmo
end subroutine f90wrap_sim_info__get__cosmo

subroutine f90wrap_sim_info__set__cosmo(this, f90wrap_cosmo)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_cosmo
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%cosmo = f90wrap_cosmo
end subroutine f90wrap_sim_info__set__cosmo

subroutine f90wrap_sim_info__get__family(this, f90wrap_family)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_family
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_family = this_ptr%p%family
end subroutine f90wrap_sim_info__get__family

subroutine f90wrap_sim_info__set__family(this, f90wrap_family)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_family
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%family = f90wrap_family
end subroutine f90wrap_sim_info__set__family

subroutine f90wrap_sim_info__get__h0(this, f90wrap_h0)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_h0
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_h0 = this_ptr%p%h0
end subroutine f90wrap_sim_info__get__h0

subroutine f90wrap_sim_info__set__h0(this, f90wrap_h0)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_h0
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%h0 = f90wrap_h0
end subroutine f90wrap_sim_info__set__h0

subroutine f90wrap_sim_info__get__t(this, f90wrap_t)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_t
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_t = this_ptr%p%t
end subroutine f90wrap_sim_info__get__t

subroutine f90wrap_sim_info__set__t(this, f90wrap_t)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_t
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%t = f90wrap_t
end subroutine f90wrap_sim_info__set__t

subroutine f90wrap_sim_info__get__aexp(this, f90wrap_aexp)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_aexp
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_aexp = this_ptr%p%aexp
end subroutine f90wrap_sim_info__get__aexp

subroutine f90wrap_sim_info__set__aexp(this, f90wrap_aexp)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_aexp
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%aexp = f90wrap_aexp
end subroutine f90wrap_sim_info__set__aexp

subroutine f90wrap_sim_info__get__unit_l(this, f90wrap_unit_l)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_unit_l
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_unit_l = this_ptr%p%unit_l
end subroutine f90wrap_sim_info__get__unit_l

subroutine f90wrap_sim_info__set__unit_l(this, f90wrap_unit_l)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_unit_l
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%unit_l = f90wrap_unit_l
end subroutine f90wrap_sim_info__set__unit_l

subroutine f90wrap_sim_info__get__unit_d(this, f90wrap_unit_d)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_unit_d
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_unit_d = this_ptr%p%unit_d
end subroutine f90wrap_sim_info__get__unit_d

subroutine f90wrap_sim_info__set__unit_d(this, f90wrap_unit_d)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_unit_d
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%unit_d = f90wrap_unit_d
end subroutine f90wrap_sim_info__set__unit_d

subroutine f90wrap_sim_info__get__unit_t(this, f90wrap_unit_t)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_unit_t
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_unit_t = this_ptr%p%unit_t
end subroutine f90wrap_sim_info__get__unit_t

subroutine f90wrap_sim_info__set__unit_t(this, f90wrap_unit_t)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_unit_t
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%unit_t = f90wrap_unit_t
end subroutine f90wrap_sim_info__set__unit_t

subroutine f90wrap_sim_info__get__unit_m(this, f90wrap_unit_m)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_unit_m
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_unit_m = this_ptr%p%unit_m
end subroutine f90wrap_sim_info__get__unit_m

subroutine f90wrap_sim_info__set__unit_m(this, f90wrap_unit_m)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_unit_m
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%unit_m = f90wrap_unit_m
end subroutine f90wrap_sim_info__set__unit_m

subroutine f90wrap_sim_info__get__boxlen(this, f90wrap_boxlen)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_boxlen
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_boxlen = this_ptr%p%boxlen
end subroutine f90wrap_sim_info__get__boxlen

subroutine f90wrap_sim_info__set__boxlen(this, f90wrap_boxlen)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_boxlen
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%boxlen = f90wrap_boxlen
end subroutine f90wrap_sim_info__set__boxlen

subroutine f90wrap_sim_info__get__omega_m(this, f90wrap_omega_m)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_omega_m
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_omega_m = this_ptr%p%omega_m
end subroutine f90wrap_sim_info__get__omega_m

subroutine f90wrap_sim_info__set__omega_m(this, f90wrap_omega_m)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_omega_m
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%omega_m = f90wrap_omega_m
end subroutine f90wrap_sim_info__set__omega_m

subroutine f90wrap_sim_info__get__omega_l(this, f90wrap_omega_l)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_omega_l
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_omega_l = this_ptr%p%omega_l
end subroutine f90wrap_sim_info__get__omega_l

subroutine f90wrap_sim_info__set__omega_l(this, f90wrap_omega_l)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_omega_l
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%omega_l = f90wrap_omega_l
end subroutine f90wrap_sim_info__set__omega_l

subroutine f90wrap_sim_info__get__omega_k(this, f90wrap_omega_k)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_omega_k
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_omega_k = this_ptr%p%omega_k
end subroutine f90wrap_sim_info__get__omega_k

subroutine f90wrap_sim_info__set__omega_k(this, f90wrap_omega_k)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_omega_k
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%omega_k = f90wrap_omega_k
end subroutine f90wrap_sim_info__set__omega_k

subroutine f90wrap_sim_info__get__omega_b(this, f90wrap_omega_b)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_omega_b
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_omega_b = this_ptr%p%omega_b
end subroutine f90wrap_sim_info__get__omega_b

subroutine f90wrap_sim_info__set__omega_b(this, f90wrap_omega_b)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_omega_b
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%omega_b = f90wrap_omega_b
end subroutine f90wrap_sim_info__set__omega_b

subroutine f90wrap_sim_info__get__time_tot(this, f90wrap_time_tot)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_time_tot
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_time_tot = this_ptr%p%time_tot
end subroutine f90wrap_sim_info__get__time_tot

subroutine f90wrap_sim_info__set__time_tot(this, f90wrap_time_tot)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_time_tot
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%time_tot = f90wrap_time_tot
end subroutine f90wrap_sim_info__set__time_tot

subroutine f90wrap_sim_info__get__time_simu(this, f90wrap_time_simu)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_time_simu
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_time_simu = this_ptr%p%time_simu
end subroutine f90wrap_sim_info__get__time_simu

subroutine f90wrap_sim_info__set__time_simu(this, f90wrap_time_simu)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_time_simu
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%time_simu = f90wrap_time_simu
end subroutine f90wrap_sim_info__set__time_simu

subroutine f90wrap_sim_info__get__n_frw(this, f90wrap_n_frw)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_n_frw
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_n_frw = this_ptr%p%n_frw
end subroutine f90wrap_sim_info__get__n_frw

subroutine f90wrap_sim_info__set__n_frw(this, f90wrap_n_frw)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in)   :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_n_frw
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%n_frw = f90wrap_n_frw
end subroutine f90wrap_sim_info__set__n_frw

subroutine f90wrap_sim_info__array__aexp_frw(this, nd, dtype, dshape, dloc)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in) :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%aexp_frw)) then
        dshape(1:1) = shape(this_ptr%p%aexp_frw)
        dloc = loc(this_ptr%p%aexp_frw)
    else
        dloc = 0
    end if
end subroutine f90wrap_sim_info__array__aexp_frw

subroutine f90wrap_sim_info__array__hexp_frw(this, nd, dtype, dshape, dloc)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in) :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%hexp_frw)) then
        dshape(1:1) = shape(this_ptr%p%hexp_frw)
        dloc = loc(this_ptr%p%hexp_frw)
    else
        dloc = 0
    end if
end subroutine f90wrap_sim_info__array__hexp_frw

subroutine f90wrap_sim_info__array__tau_frw(this, nd, dtype, dshape, dloc)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in) :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%tau_frw)) then
        dshape(1:1) = shape(this_ptr%p%tau_frw)
        dloc = loc(this_ptr%p%tau_frw)
    else
        dloc = 0
    end if
end subroutine f90wrap_sim_info__array__tau_frw

subroutine f90wrap_sim_info__array__t_frw(this, nd, dtype, dshape, dloc)
    use io_ramses, only: sim_info
    implicit none
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    integer, intent(in) :: this(2)
    type(sim_info_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%t_frw)) then
        dshape(1:1) = shape(this_ptr%p%t_frw)
        dloc = loc(this_ptr%p%t_frw)
    else
        dloc = 0
    end if
end subroutine f90wrap_sim_info__array__t_frw

subroutine f90wrap_sim_info_initialise(this)
    use io_ramses, only: sim_info
    implicit none
    
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    type(sim_info_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_sim_info_initialise

subroutine f90wrap_sim_info_finalise(this)
    use io_ramses, only: sim_info
    implicit none
    
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    type(sim_info_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_sim_info_finalise

subroutine f90wrap_level__get__ilevel(this, f90wrap_ilevel)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_ilevel
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ilevel = this_ptr%p%ilevel
end subroutine f90wrap_level__get__ilevel

subroutine f90wrap_level__set__ilevel(this, f90wrap_ilevel)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_ilevel
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ilevel = f90wrap_ilevel
end subroutine f90wrap_level__set__ilevel

subroutine f90wrap_level__get__ngrid(this, f90wrap_ngrid)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_ngrid
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ngrid = this_ptr%p%ngrid
end subroutine f90wrap_level__get__ngrid

subroutine f90wrap_level__set__ngrid(this, f90wrap_ngrid)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_ngrid
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ngrid = f90wrap_ngrid
end subroutine f90wrap_level__set__ngrid

subroutine f90wrap_level__array__cube(this, nd, dtype, dshape, dloc)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in) :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:3) = shape(this_ptr%p%cube)
    dloc = loc(this_ptr%p%cube)
end subroutine f90wrap_level__array__cube

subroutine f90wrap_level__array__map(this, nd, dtype, dshape, dloc)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in) :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:3) = shape(this_ptr%p%map)
    dloc = loc(this_ptr%p%map)
end subroutine f90wrap_level__array__map

subroutine f90wrap_level__array__rho(this, nd, dtype, dshape, dloc)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in) :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%rho)
    dloc = loc(this_ptr%p%rho)
end subroutine f90wrap_level__array__rho

subroutine f90wrap_level__get__imin(this, f90wrap_imin)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_imin
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_imin = this_ptr%p%imin
end subroutine f90wrap_level__get__imin

subroutine f90wrap_level__set__imin(this, f90wrap_imin)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_imin
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%imin = f90wrap_imin
end subroutine f90wrap_level__set__imin

subroutine f90wrap_level__get__imax(this, f90wrap_imax)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_imax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_imax = this_ptr%p%imax
end subroutine f90wrap_level__get__imax

subroutine f90wrap_level__set__imax(this, f90wrap_imax)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_imax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%imax = f90wrap_imax
end subroutine f90wrap_level__set__imax

subroutine f90wrap_level__get__jmin(this, f90wrap_jmin)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_jmin
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_jmin = this_ptr%p%jmin
end subroutine f90wrap_level__get__jmin

subroutine f90wrap_level__set__jmin(this, f90wrap_jmin)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_jmin
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%jmin = f90wrap_jmin
end subroutine f90wrap_level__set__jmin

subroutine f90wrap_level__get__jmax(this, f90wrap_jmax)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_jmax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_jmax = this_ptr%p%jmax
end subroutine f90wrap_level__get__jmax

subroutine f90wrap_level__set__jmax(this, f90wrap_jmax)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_jmax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%jmax = f90wrap_jmax
end subroutine f90wrap_level__set__jmax

subroutine f90wrap_level__get__kmin(this, f90wrap_kmin)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_kmin
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_kmin = this_ptr%p%kmin
end subroutine f90wrap_level__get__kmin

subroutine f90wrap_level__set__kmin(this, f90wrap_kmin)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_kmin
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%kmin = f90wrap_kmin
end subroutine f90wrap_level__set__kmin

subroutine f90wrap_level__get__kmax(this, f90wrap_kmax)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_kmax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_kmax = this_ptr%p%kmax
end subroutine f90wrap_level__get__kmax

subroutine f90wrap_level__set__kmax(this, f90wrap_kmax)
    use io_ramses, only: level
    implicit none
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    integer, intent(in)   :: this(2)
    type(level_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_kmax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%kmax = f90wrap_kmax
end subroutine f90wrap_level__set__kmax

subroutine f90wrap_level_initialise(this)
    use io_ramses, only: level
    implicit none
    
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    type(level_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_level_initialise

subroutine f90wrap_level_finalise(this)
    use io_ramses, only: level
    implicit none
    
    type level_ptr_type
        type(level), pointer :: p => NULL()
    end type level_ptr_type
    type(level_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_level_finalise

subroutine f90wrap_data_handler__get__name(this, f90wrap_name)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    character(80), intent(out) :: f90wrap_name
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_name = this_ptr%p%name
end subroutine f90wrap_data_handler__get__name

subroutine f90wrap_data_handler__set__name(this, f90wrap_name)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    character(80), intent(in) :: f90wrap_name
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%name = f90wrap_name
end subroutine f90wrap_data_handler__set__name

subroutine f90wrap_data_handler__get__x_data(this, f90wrap_x_data)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_x_data
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_x_data = this_ptr%p%x_data
end subroutine f90wrap_data_handler__get__x_data

subroutine f90wrap_data_handler__set__x_data(this, f90wrap_x_data)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_x_data
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%x_data = f90wrap_x_data
end subroutine f90wrap_data_handler__set__x_data

subroutine f90wrap_data_handler__get__y_data(this, f90wrap_y_data)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_y_data
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_y_data = this_ptr%p%y_data
end subroutine f90wrap_data_handler__get__y_data

subroutine f90wrap_data_handler__set__y_data(this, f90wrap_y_data)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_y_data
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%y_data = f90wrap_y_data
end subroutine f90wrap_data_handler__set__y_data

subroutine f90wrap_data_handler__get__z_data(this, f90wrap_z_data)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_z_data
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_z_data = this_ptr%p%z_data
end subroutine f90wrap_data_handler__get__z_data

subroutine f90wrap_data_handler__set__z_data(this, f90wrap_z_data)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in)   :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_z_data
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%z_data = f90wrap_z_data
end subroutine f90wrap_data_handler__set__z_data

subroutine f90wrap_data_handler__array__nx(this, nd, dtype, dshape, dloc)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in) :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%nx)
    dloc = loc(this_ptr%p%nx)
end subroutine f90wrap_data_handler__array__nx

subroutine f90wrap_data_handler__array__ny(this, nd, dtype, dshape, dloc)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in) :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%ny)
    dloc = loc(this_ptr%p%ny)
end subroutine f90wrap_data_handler__array__ny

subroutine f90wrap_data_handler__array__nz(this, nd, dtype, dshape, dloc)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in) :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%nz)
    dloc = loc(this_ptr%p%nz)
end subroutine f90wrap_data_handler__array__nz

subroutine f90wrap_data_handler__array__x(this, nd, dtype, dshape, dloc)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in) :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%x)) then
        dshape(1:3) = shape(this_ptr%p%x)
        dloc = loc(this_ptr%p%x)
    else
        dloc = 0
    end if
end subroutine f90wrap_data_handler__array__x

subroutine f90wrap_data_handler__array__y(this, nd, dtype, dshape, dloc)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in) :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%y)) then
        dshape(1:3) = shape(this_ptr%p%y)
        dloc = loc(this_ptr%p%y)
    else
        dloc = 0
    end if
end subroutine f90wrap_data_handler__array__y

subroutine f90wrap_data_handler__array__z(this, nd, dtype, dshape, dloc)
    use io_ramses, only: data_handler
    implicit none
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    integer, intent(in) :: this(2)
    type(data_handler_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%z)) then
        dshape(1:3) = shape(this_ptr%p%z)
        dloc = loc(this_ptr%p%z)
    else
        dloc = 0
    end if
end subroutine f90wrap_data_handler__array__z

subroutine f90wrap_data_handler_initialise(this)
    use io_ramses, only: data_handler
    implicit none
    
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    type(data_handler_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_data_handler_initialise

subroutine f90wrap_data_handler_finalise(this)
    use io_ramses, only: data_handler
    implicit none
    
    type data_handler_ptr_type
        type(data_handler), pointer :: p => NULL()
    end type data_handler_ptr_type
    type(data_handler_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_data_handler_finalise

subroutine f90wrap_particle__get__id(this, f90wrap_id)
    use io_ramses, only: particle
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_id
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_id = this_ptr%p%id
end subroutine f90wrap_particle__get__id

subroutine f90wrap_particle__set__id(this, f90wrap_id)
    use io_ramses, only: particle
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_id
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%id = f90wrap_id
end subroutine f90wrap_particle__set__id

subroutine f90wrap_particle__get__x(this, f90wrap_x)
    use io_ramses, only: particle
    use vectors, only: vector
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_x(2)
    type(vector_ptr_type) :: x_ptr
    
    this_ptr = transfer(this, this_ptr)
    x_ptr%p => this_ptr%p%x
    f90wrap_x = transfer(x_ptr,f90wrap_x)
end subroutine f90wrap_particle__get__x

subroutine f90wrap_particle__set__x(this, f90wrap_x)
    use io_ramses, only: particle
    use vectors, only: vector
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_x(2)
    type(vector_ptr_type) :: x_ptr
    
    this_ptr = transfer(this, this_ptr)
    x_ptr = transfer(f90wrap_x,x_ptr)
    this_ptr%p%x = x_ptr%p
end subroutine f90wrap_particle__set__x

subroutine f90wrap_particle__get__v(this, f90wrap_v)
    use io_ramses, only: particle
    use vectors, only: vector
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_v(2)
    type(vector_ptr_type) :: v_ptr
    
    this_ptr = transfer(this, this_ptr)
    v_ptr%p => this_ptr%p%v
    f90wrap_v = transfer(v_ptr,f90wrap_v)
end subroutine f90wrap_particle__get__v

subroutine f90wrap_particle__set__v(this, f90wrap_v)
    use io_ramses, only: particle
    use vectors, only: vector
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_v(2)
    type(vector_ptr_type) :: v_ptr
    
    this_ptr = transfer(this, this_ptr)
    v_ptr = transfer(f90wrap_v,v_ptr)
    this_ptr%p%v = v_ptr%p
end subroutine f90wrap_particle__set__v

subroutine f90wrap_particle__get__m(this, f90wrap_m)
    use io_ramses, only: particle
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_m
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_m = this_ptr%p%m
end subroutine f90wrap_particle__get__m

subroutine f90wrap_particle__set__m(this, f90wrap_m)
    use io_ramses, only: particle
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_m
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%m = f90wrap_m
end subroutine f90wrap_particle__set__m

subroutine f90wrap_particle__get__met(this, f90wrap_met)
    use io_ramses, only: particle
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_met
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_met = this_ptr%p%met
end subroutine f90wrap_particle__get__met

subroutine f90wrap_particle__set__met(this, f90wrap_met)
    use io_ramses, only: particle
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_met
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%met = f90wrap_met
end subroutine f90wrap_particle__set__met

subroutine f90wrap_particle__get__imass(this, f90wrap_imass)
    use io_ramses, only: particle
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_imass
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_imass = this_ptr%p%imass
end subroutine f90wrap_particle__get__imass

subroutine f90wrap_particle__set__imass(this, f90wrap_imass)
    use io_ramses, only: particle
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_imass
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%imass = f90wrap_imass
end subroutine f90wrap_particle__set__imass

subroutine f90wrap_particle__get__age(this, f90wrap_age)
    use io_ramses, only: particle
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_age
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_age = this_ptr%p%age
end subroutine f90wrap_particle__get__age

subroutine f90wrap_particle__set__age(this, f90wrap_age)
    use io_ramses, only: particle
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_age
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%age = f90wrap_age
end subroutine f90wrap_particle__set__age

subroutine f90wrap_particle__get__tform(this, f90wrap_tform)
    use io_ramses, only: particle
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_tform
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_tform = this_ptr%p%tform
end subroutine f90wrap_particle__get__tform

subroutine f90wrap_particle__set__tform(this, f90wrap_tform)
    use io_ramses, only: particle
    implicit none
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    integer, intent(in)   :: this(2)
    type(particle_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_tform
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%tform = f90wrap_tform
end subroutine f90wrap_particle__set__tform

subroutine f90wrap_particle_initialise(this)
    use io_ramses, only: particle
    implicit none
    
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    type(particle_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_particle_initialise

subroutine f90wrap_particle_finalise(this)
    use io_ramses, only: particle
    implicit none
    
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    type(particle_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_particle_finalise

subroutine f90wrap_title(n, nchar)
    use io_ramses, only: title
    implicit none
    
    integer, intent(in) :: n
    character(5), intent(inout) :: nchar
    call title(n=n, nchar=nchar)
end subroutine f90wrap_title

subroutine f90wrap_hilbert3d(x, y, z, order, bit_length, npoint, n0, n1, n2, n3)
    use io_ramses, only: hilbert3d
    implicit none
    
    integer, intent(in), dimension(n0) :: x
    integer, intent(in), dimension(n1) :: y
    integer, intent(in), dimension(n2) :: z
    real(8), intent(inout), dimension(n3) :: order
    integer, intent(in) :: bit_length
    integer, intent(in) :: npoint
    integer :: n0
    !f2py intent(hide), depend(x) :: n0 = shape(x,0)
    integer :: n1
    !f2py intent(hide), depend(y) :: n1 = shape(y,0)
    integer :: n2
    !f2py intent(hide), depend(z) :: n2 = shape(z,0)
    integer :: n3
    !f2py intent(hide), depend(order) :: n3 = shape(order,0)
    call hilbert3d(x=x, y=y, z=z, order=order, bit_length=bit_length, npoint=npoint)
end subroutine f90wrap_hilbert3d

subroutine f90wrap_check_lmax(ngridfile, amr, n0, n1)
    use io_ramses, only: amr_info, check_lmax
    implicit none
    
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    integer, intent(in), dimension(n0,n1) :: ngridfile
    type(amr_info_ptr_type) :: amr_ptr
    integer, intent(in), dimension(2) :: amr
    integer :: n0
    !f2py intent(hide), depend(ngridfile) :: n0 = shape(ngridfile,0)
    integer :: n1
    !f2py intent(hide), depend(ngridfile) :: n1 = shape(ngridfile,1)
    amr_ptr = transfer(amr, amr_ptr)
    call check_lmax(ngridfile=ngridfile, amr=amr_ptr%p)
end subroutine f90wrap_check_lmax

subroutine f90wrap_check_families(repository, sim)
    use io_ramses, only: check_families, sim_info
    implicit none
    
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    character(128), intent(in) :: repository
    type(sim_info_ptr_type) :: sim_ptr
    integer, intent(in), dimension(2) :: sim
    sim_ptr = transfer(sim, sim_ptr)
    call check_families(repository=repository, sim=sim_ptr%p)
end subroutine f90wrap_check_families

subroutine f90wrap_read_hydrofile_descriptor(repository, varids)
    use io_ramses, only: read_hydrofile_descriptor, hydroid
    implicit none
    
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    character(128), intent(in) :: repository
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    varids_ptr = transfer(varids, varids_ptr)
    call read_hydrofile_descriptor(repository=repository, varIDs=varids_ptr%p)
end subroutine f90wrap_read_hydrofile_descriptor

subroutine f90wrap_read_hydrofile_descriptor_old(repository, varids)
    use io_ramses, only: read_hydrofile_descriptor_old, hydroid
    implicit none
    
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    character(128), intent(in) :: repository
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    varids_ptr = transfer(varids, varids_ptr)
    call read_hydrofile_descriptor_old(repository=repository, varIDs=varids_ptr%p)
end subroutine f90wrap_read_hydrofile_descriptor_old

subroutine f90wrap_select_from_descriptor_ids(varids, newvar, newid)
    use io_ramses, only: select_from_descriptor_ids, hydroid
    implicit none
    
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    character(25), intent(in) :: newvar
    integer, intent(in) :: newid
    varids_ptr = transfer(varids, varids_ptr)
    call select_from_descriptor_ids(varIDs=varids_ptr%p, newvar=newvar, newID=newid)
end subroutine f90wrap_select_from_descriptor_ids

subroutine f90wrap_read_hydrofile_descriptor_new(repository, varids)
    use io_ramses, only: read_hydrofile_descriptor_new, hydroid
    implicit none
    
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    character(128), intent(in) :: repository
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    varids_ptr = transfer(varids, varids_ptr)
    call read_hydrofile_descriptor_new(repository=repository, varIDs=varids_ptr%p)
end subroutine f90wrap_read_hydrofile_descriptor_new

subroutine f90wrap_getvarvalue(varids, reg, dx, x, var, varname, value, n0)
    use io_ramses, only: getvarvalue, hydroid
    use vectors, only: vector
    use geometrical_regions, only: region
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    real(8), intent(in) :: dx
    type(vector_ptr_type) :: x_ptr
    integer, intent(in), dimension(2) :: x
    real(8), intent(in), dimension(n0) :: var
    character(128), intent(in) :: varname
    real(8), intent(inout) :: value
    integer :: n0
    !f2py intent(hide), depend(var) :: n0 = shape(var,0)
    varids_ptr = transfer(varids, varids_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    x_ptr = transfer(x, x_ptr)
    call getvarvalue(varIDs=varids_ptr%p, reg=reg_ptr%p, dx=dx, x=x_ptr%p, var=var, varname=varname, value=value)
end subroutine f90wrap_getvarvalue

subroutine f90wrap_init_amr_read(repository, amr, sim)
    use io_ramses, only: init_amr_read, amr_info, sim_info
    implicit none
    
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    character(128), intent(in) :: repository
    type(amr_info_ptr_type) :: amr_ptr
    integer, intent(in), dimension(2) :: amr
    type(sim_info_ptr_type) :: sim_ptr
    integer, intent(in), dimension(2) :: sim
    amr_ptr = transfer(amr, amr_ptr)
    sim_ptr = transfer(sim, sim_ptr)
    call init_amr_read(repository=repository, amr=amr_ptr%p, sim=sim_ptr%p)
end subroutine f90wrap_init_amr_read

subroutine f90wrap_get_cpu_map(reg, amr)
    use io_ramses, only: get_cpu_map, amr_info
    use geometrical_regions, only: region
    implicit none
    
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type amr_info_ptr_type
        type(amr_info), pointer :: p => NULL()
    end type amr_info_ptr_type
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(amr_info_ptr_type) :: amr_ptr
    integer, intent(in), dimension(2) :: amr
    reg_ptr = transfer(reg, reg_ptr)
    amr_ptr = transfer(amr, amr_ptr)
    call get_cpu_map(reg=reg_ptr%p, amr=amr_ptr%p)
end subroutine f90wrap_get_cpu_map

subroutine f90wrap_getparttype(part, ptype)
    use io_ramses, only: particle, getparttype
    implicit none
    
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    type(particle_ptr_type) :: part_ptr
    integer, intent(in), dimension(2) :: part
    character(6), intent(inout) :: ptype
    part_ptr = transfer(part, part_ptr)
    call getparttype(part=part_ptr%p, ptype=ptype)
end subroutine f90wrap_getparttype

subroutine f90wrap_getpartvalue(sim, reg, part, var, value, dx)
    use vectors, only: vector
    use io_ramses, only: particle, getpartvalue, sim_info
    use geometrical_regions, only: region
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    type(sim_info_ptr_type) :: sim_ptr
    integer, intent(in), dimension(2) :: sim
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(particle_ptr_type) :: part_ptr
    integer, intent(in), dimension(2) :: part
    character(128), intent(in) :: var
    real(8), intent(inout) :: value
    type(vector_ptr_type) :: dx_ptr
    integer, optional, intent(in), dimension(2) :: dx
    sim_ptr = transfer(sim, sim_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    part_ptr = transfer(part, part_ptr)
    if (present(dx)) then
        dx_ptr = transfer(dx, dx_ptr)
    else
        dx_ptr%p => null()
    end if
    call getpartvalue(sim=sim_ptr%p, reg=reg_ptr%p, part=part_ptr%p, var=var, value=value, dx=dx_ptr%p)
end subroutine f90wrap_getpartvalue

! End of module io_ramses defined in file read_amr_module.fpp

! Module filtering defined in file read_amr_module.fpp

subroutine f90wrap_filter__get__name(this, f90wrap_name)
    use filtering, only: filter
    implicit none
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    integer, intent(in)   :: this(2)
    type(filter_ptr_type) :: this_ptr
    character(128), intent(out) :: f90wrap_name
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_name = this_ptr%p%name
end subroutine f90wrap_filter__get__name

subroutine f90wrap_filter__set__name(this, f90wrap_name)
    use filtering, only: filter
    implicit none
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    integer, intent(in)   :: this(2)
    type(filter_ptr_type) :: this_ptr
    character(128), intent(in) :: f90wrap_name
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%name = f90wrap_name
end subroutine f90wrap_filter__set__name

subroutine f90wrap_filter__get__ncond(this, f90wrap_ncond)
    use filtering, only: filter
    implicit none
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    integer, intent(in)   :: this(2)
    type(filter_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_ncond
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ncond = this_ptr%p%ncond
end subroutine f90wrap_filter__get__ncond

subroutine f90wrap_filter__set__ncond(this, f90wrap_ncond)
    use filtering, only: filter
    implicit none
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    integer, intent(in)   :: this(2)
    type(filter_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_ncond
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ncond = f90wrap_ncond
end subroutine f90wrap_filter__set__ncond

subroutine f90wrap_filter__array__cond_vars(this, nd, dtype, dshape, dloc)
    use filtering, only: filter
    implicit none
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    integer, intent(in) :: this(2)
    type(filter_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%cond_vars)) then
        dshape(1:2) = (/len(this_ptr%p%cond_vars(1)), shape(this_ptr%p%cond_vars)/)
        dloc = loc(this_ptr%p%cond_vars)
    else
        dloc = 0
    end if
end subroutine f90wrap_filter__array__cond_vars

subroutine f90wrap_filter__array__cond_ops(this, nd, dtype, dshape, dloc)
    use filtering, only: filter
    implicit none
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    integer, intent(in) :: this(2)
    type(filter_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%cond_ops)) then
        dshape(1:2) = (/len(this_ptr%p%cond_ops(1)), shape(this_ptr%p%cond_ops)/)
        dloc = loc(this_ptr%p%cond_ops)
    else
        dloc = 0
    end if
end subroutine f90wrap_filter__array__cond_ops

subroutine f90wrap_filter__array__cond_vals(this, nd, dtype, dshape, dloc)
    use filtering, only: filter
    implicit none
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    integer, intent(in) :: this(2)
    type(filter_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%cond_vals)) then
        dshape(1:1) = shape(this_ptr%p%cond_vals)
        dloc = loc(this_ptr%p%cond_vals)
    else
        dloc = 0
    end if
end subroutine f90wrap_filter__array__cond_vals

subroutine f90wrap_filter_initialise(this)
    use filtering, only: filter
    implicit none
    
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    type(filter_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_filter_initialise

subroutine f90wrap_filter_finalise(this)
    use filtering, only: filter
    implicit none
    
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    type(filter_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_filter_finalise

subroutine f90wrap_allocate_filter(filt)
    use filtering, only: allocate_filter, filter
    implicit none
    
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    type(filter_ptr_type) :: filt_ptr
    integer, intent(in), dimension(2) :: filt
    filt_ptr = transfer(filt, filt_ptr)
    call allocate_filter(filt=filt_ptr%p)
end subroutine f90wrap_allocate_filter

subroutine f90wrap_cond_string_to_filter(str, filt)
    use filtering, only: filter, cond_string_to_filter
    implicit none
    
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    character(128), intent(in) :: str
    type(filter_ptr_type) :: filt_ptr
    integer, intent(in), dimension(2) :: filt
    filt_ptr = transfer(filt, filt_ptr)
    call cond_string_to_filter(str=str, filt=filt_ptr%p)
end subroutine f90wrap_cond_string_to_filter

subroutine f90wrap_filter_cell(varids, reg, filt, cell_x, cell_dx, ret_filter_cell, cell_var, n0)
    use vectors, only: vector
    use filtering, only: filter_cell, filter
    use geometrical_regions, only: region
    use io_ramses, only: hydroid
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    type hydroid_ptr_type
        type(hydroid), pointer :: p => NULL()
    end type hydroid_ptr_type
    type(hydroid_ptr_type) :: varids_ptr
    integer, intent(in), dimension(2) :: varids
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(filter_ptr_type) :: filt_ptr
    integer, intent(in), dimension(2) :: filt
    type(vector_ptr_type) :: cell_x_ptr
    integer, intent(in), dimension(2) :: cell_x
    real(8), intent(in) :: cell_dx
    logical, intent(out) :: ret_filter_cell
    real(8), intent(in), dimension(n0) :: cell_var
    integer :: n0
    !f2py intent(hide), depend(cell_var) :: n0 = shape(cell_var,0)
    varids_ptr = transfer(varids, varids_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    filt_ptr = transfer(filt, filt_ptr)
    cell_x_ptr = transfer(cell_x, cell_x_ptr)
    ret_filter_cell = filter_cell(varIDs=varids_ptr%p, reg=reg_ptr%p, filt=filt_ptr%p, cell_x=cell_x_ptr%p, cell_dx=cell_dx, &
        cell_var=cell_var)
end subroutine f90wrap_filter_cell

subroutine f90wrap_filter_particle(sim, reg, filt, part, ret_filter_particle, dx)
    use vectors, only: vector
    use filtering, only: filter, filter_particle
    use io_ramses, only: particle, sim_info
    use geometrical_regions, only: region
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    type region_ptr_type
        type(region), pointer :: p => NULL()
    end type region_ptr_type
    type filter_ptr_type
        type(filter), pointer :: p => NULL()
    end type filter_ptr_type
    type particle_ptr_type
        type(particle), pointer :: p => NULL()
    end type particle_ptr_type
    type(sim_info_ptr_type) :: sim_ptr
    integer, intent(in), dimension(2) :: sim
    type(region_ptr_type) :: reg_ptr
    integer, intent(in), dimension(2) :: reg
    type(filter_ptr_type) :: filt_ptr
    integer, intent(in), dimension(2) :: filt
    type(particle_ptr_type) :: part_ptr
    integer, intent(in), dimension(2) :: part
    logical, intent(out) :: ret_filter_particle
    type(vector_ptr_type) :: dx_ptr
    integer, optional, intent(in), dimension(2) :: dx
    sim_ptr = transfer(sim, sim_ptr)
    reg_ptr = transfer(reg, reg_ptr)
    filt_ptr = transfer(filt, filt_ptr)
    part_ptr = transfer(part, part_ptr)
    if (present(dx)) then
        dx_ptr = transfer(dx, dx_ptr)
    else
        dx_ptr%p => null()
    end if
    ret_filter_particle = filter_particle(sim=sim_ptr%p, reg=reg_ptr%p, filt=filt_ptr%p, part=part_ptr%p, dx=dx_ptr%p)
end subroutine f90wrap_filter_particle

! End of module filtering defined in file read_amr_module.fpp

