!--------------------------------------------------------------------------
! ozymandias:statistics.f90
!--------------------------------------------------------------------------
!
! MODULE: stats_utils
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> types and routines used when computing statistical quantities from
!> raw RAMSES data.
!
!> @details  
!> 
! 
!
!> @date 31/10/2022   0.2 basic developer version ready
!--------------------------------------------------------------------------

module stats_utils
    use local
    use constants
    use io_ramses

    type pdf_handler
        integer :: nbins,nwvars,nfilter
        logical :: do_binning
        character(128) :: varname,scaletype
        character(128),dimension(:),allocatable :: wvarnames
        real(dbl),dimension(:),allocatable :: maxv, minv
        real(dbl),dimension(:),allocatable :: bins
        real(dbl),dimension(:,:,:),allocatable :: heights
        real(dbl),dimension(:,:),allocatable :: totweights
        real(dbl),dimension(:,:,:),allocatable :: total
    end type pdf_handler


    contains

    subroutine allocate_pdf(mypdf)
        implicit none
        type(pdf_handler),intent(inout) :: mypdf

        if (.not.allocated(mypdf%maxv)) allocate(mypdf%maxv(1:mypdf%nfilter))
        if (.not.allocated(mypdf%minv)) allocate(mypdf%minv(1:mypdf%nfilter))
        if (.not.allocated(mypdf%bins)) allocate(mypdf%bins(0:mypdf%nbins))
        if (.not.allocated(mypdf%heights)) allocate(mypdf%heights(1:mypdf%nfilter,1:mypdf%nwvars,1:mypdf%nbins))
        if (.not.allocated(mypdf%wvarnames)) allocate(mypdf%wvarnames(mypdf%nwvars))
        if (.not.allocated(mypdf%totweights)) allocate(mypdf%totweights(1:mypdf%nfilter,1:mypdf%nwvars))
        if (.not.allocated(mypdf%total)) allocate(mypdf%total(1:mypdf%nfilter,1:mypdf%nwvars,2))

        ! Just make sure that initial values are zero
        mypdf%do_binning = .true.
        mypdf%maxv(:) = 0D0
        mypdf%minv(:) = 0D0
        mypdf%bins(:) = 0D0
        mypdf%heights(:,:,:) = 0D0
        mypdf%totweights(:,:) = 0D0
        mypdf%total(:,:,:) = 0D0
    end subroutine allocate_pdf

    function makebins(reg,varname,nbins,scaletype)
        use geometrical_regions

        implicit none
        type(region),intent(in) :: reg
        character(128),intent(in) :: varname
        integer,intent(in) :: nbins
        character(128),intent(in) :: scaletype

        real(dbl),dimension(0:nbins) :: makebins

        integer :: n
        real(dbl) :: rmin
        logical :: logscale

        logscale = .false.
        if (trim(scaletype).eq.'log_even') logscale = .true.

        select case (TRIM(varname))
        case('r_sphere','r_cyl')
            rmin = max(1D0/(2D0**(amr%nlevelmax-1)),1D-3*reg%rmax)
            do n=0,nbins
                if (logscale) then
                    if (reg%rmin.eq.0D0) then
                        makebins(n) = dble(n)*(log10(reg%rmax)-log10(rmin))/dble(nbins) + log10(rmin)
                    else
                        makebins(n) = dble(n)*(log10(reg%rmax)-log10(reg%rmin))/dble(nbins) + log10(reg%rmin)
                    end if
                else
                    makebins(n) = dble(n)*(reg%rmax-reg%rmin)/dble(nbins) + reg%rmin
                endif
            end do
        case('z')
            rmin = max(1D0/(2D0**(amr%nlevelmax-1)),1D-3*reg%zmax)
            do n=0,nbins
                if (logscale) then
                    if (reg%zmin.eq.0D0) then
                        makebins(n) = dble(n)*(log10(reg%zmax)-log10(rmin))/dble(nbins) + log10(rmin)
                    else
                        makebins(n) = dble(n)*(log10(reg%zmax)-log10(reg%zmin))/dble(nbins) + log10(reg%zmin)
                    end if
                else
                    makebins(n) = dble(n)*(reg%zmax-reg%zmin)/dble(nbins) + reg%zmin
                endif
            end do
        case('density')
            do n=0,nbins
                if (logscale) then
                    makebins(n) = dble(n)*(log10(1D-20/sim%unit_d)-log10(1D-30/sim%unit_d))/dble(nbins) + log10(1D-30/sim%unit_d)
                else
                    makebins(n) = dble(n)*(1D-20/sim%unit_d - 1D-30/sim%unit_d)/dble(nbins) + 1D-30/sim%unit_d
                endif
            end do
        case('temperature')
            do n=0,nbins
                if (logscale) then
                    makebins(n) = dble(n)*(log10(1D8/sim%T2)-log10(1D0/sim%T2))/dble(nbins) + log10(1D0/sim%T2)
                else
                    makebins(n) = dble(n)*(1D8/sim%T2 - 1D0/sim%T2)/dble(nbins) + 1D0/sim%T2
                endif
            end do
        case('total_coolingtime')
            do n=0,nbins
                if (logscale) then
                    makebins(n) = dble(n)*(log10(3.15D18/sim%unit_t)-log10(3.15D13/sim%unit_t))/dble(nbins) + log10(3.15D13/sim%unit_t)
                else
                    makebins(n) = dble(n)*(3.15D18/sim%unit_t - 3.15D13/sim%unit_t)/dble(nbins) + 3.15D13/sim%unit_t
                endif
            end do
        case('thermal_pressure')
            do n=0,nbins
                if (logscale) then
                    makebins(n) = dble(n)*(log10(1D-9/sim%unit_p)-log10(1D-16/sim%unit_p))/dble(nbins) + log10(1D-16/sim%unit_p)
                else
                    makebins(n) = dble(n)*(1D-9/sim%unit_p - 1D-16/sim%unit_p)/dble(nbins) + 1D-16/sim%unit_p
                endif
            end do
        case('v_sphere_r')
            do n=0,nbins
                if (logscale) then
                    write(*,*)'You cannot use logscale for a velocity profile. Stopping!'
                    stop
                else
                    makebins(n) = dble(n)*(4D7/sim%unit_v + 4D7/sim%unit_v)/dble(nbins) - 4D7/sim%unit_v
                endif
            end do
        case('theta_sphere')
            do n=0,nbins
                if (logscale) then
                    makebins(n) = dble(n)*(log10(pi))/dble(nbins)
                else
                    makebins(n) = dble(n)*(pi + 0D0)/dble(nbins) - 0D0
                endif
            end do
        !TODO: Add more cases
        case default
            write(*,*)'Variable not supported for makebins: ',TRIM(varname)
            write(*,*)'Aborting!'
            stop
        end select
    end function makebins

    subroutine findbinpos(reg,x,cvars,csons,csize,&
                            & ibin,value,trans_matrix,&
                            & scaletype,nbins,bins,vname,&
                            & gvars)
        use vectors
        use geometrical_regions
        implicit none
        type(region),intent(in) :: reg
        type(vector),intent(in) :: x
        real(dbl),dimension(0:amr%twondim,1:varIDs%nvar),intent(in) :: cvars
        integer,dimension(0:amr%twondim),intent(in) :: csons
        real(dbl),intent(in) :: csize
        integer,intent(inout) :: ibin
        real(dbl),intent(inout) :: value
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix
        character(128),intent(in) :: scaletype
        integer,intent(in) :: nbins
        real(dbl) :: origvalue
        real(dbl),dimension(0:nbins) :: bins
        character(128),intent(in) :: vname
        real(dbl),dimension(0:amr%twondim,1:4),optional,intent(in) :: gvars

        ! Get variable value
        if (present(gvars)) then
            call getvarvalue(reg,csize,x,cvars,csons,vname,value,trans_matrix,gvars)
        else
            call getvarvalue(reg,csize,x,cvars,csons,vname,value,trans_matrix)
        end if
        origvalue = value

        ! Transform value depending how the bins are provided
        if (trim(scaletype).eq.'log_even') then
            value = log10(value)
            ibin = int(dble(nbins)*(value-bins(0))/(bins(nbins)-bins(0))) + 1
        else if (trim(scaletype).eq.'linear_even') then
            ibin = int(dble(nbins)*(value-bins(0))/(bins(nbins)-bins(0))) + 1
        else
            ibin = 1
            do while (ibin .lt. nbins)
                if (value .le. bins(ibin)) exit
                ibin = ibin + 1
            end do
        end if
        
        ! Make sure we are not out of the outer boundaries
        if (value .eq. bins(nbins)) then
            ibin = nbins
        else if (value<bins(0).or.value>bins(nbins)) then
            ibin = 0
        else if (ibin < 0) then
            ibin = 0
        end if
        value = origvalue
    end subroutine findbinpos
end module stats_utils