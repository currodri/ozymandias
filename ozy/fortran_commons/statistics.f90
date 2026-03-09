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
        integer :: nbins,nvars,nwvars,nfilter
        integer,dimension(:,:),allocatable :: nvalues,nout
        integer,dimension(:),allocatable :: zero_index
        logical,dimension(:),allocatable :: do_binning
        character(128),dimension(:),allocatable :: varname,scaletype
        character(128),dimension(:),allocatable :: wvarnames
        real(dbl),dimension(:),allocatable :: linthresh
        real(dbl),dimension(:,:),allocatable :: maxv, minv
        real(dbl),dimension(:,:),allocatable :: bins
        real(dbl),dimension(:,:,:,:),allocatable :: heights
        real(dbl),dimension(:,:,:),allocatable :: totweights
        real(dbl),dimension(:,:,:,:),allocatable :: total
    end type pdf_handler


    contains

    subroutine allocate_pdf(mypdf)
        implicit none
        type(pdf_handler),intent(inout) :: mypdf

        if (.not.allocated(mypdf%do_binning)) allocate(mypdf%do_binning(1:mypdf%nvars))
        if (.not.allocated(mypdf%maxv)) allocate(mypdf%maxv(1:mypdf%nvars,1:mypdf%nfilter))
        if (.not.allocated(mypdf%minv)) allocate(mypdf%minv(1:mypdf%nvars,1:mypdf%nfilter))
        if (.not.allocated(mypdf%bins)) allocate(mypdf%bins(0:mypdf%nbins,1:mypdf%nvars))
        if (.not.allocated(mypdf%zero_index)) allocate(mypdf%zero_index(1:mypdf%nvars))
        if (.not.allocated(mypdf%linthresh)) allocate(mypdf%linthresh(1:mypdf%nvars))
        if (.not.allocated(mypdf%heights)) allocate(mypdf%heights(1:mypdf%nvars,1:mypdf%nfilter,1:mypdf%nwvars,1:mypdf%nbins))
        if (.not.allocated(mypdf%wvarnames)) allocate(mypdf%wvarnames(mypdf%nwvars))
        if (.not.allocated(mypdf%varname)) allocate(mypdf%varname(mypdf%nvars))
        if (.not.allocated(mypdf%scaletype)) allocate(mypdf%scaletype(mypdf%nvars))
        if (.not.allocated(mypdf%totweights)) allocate(mypdf%totweights(1:mypdf%nvars,1:mypdf%nfilter,1:mypdf%nwvars))
        if (.not.allocated(mypdf%total)) allocate(mypdf%total(1:mypdf%nvars,1:mypdf%nfilter,1:mypdf%nwvars,2))
        if (.not.allocated(mypdf%nvalues)) allocate(mypdf%nvalues(1:mypdf%nvars,1:mypdf%nfilter))
        if (.not.allocated(mypdf%nout)) allocate(mypdf%nout(1:mypdf%nvars,1:mypdf%nfilter))
        ! Just make sure that initial values are zero
        mypdf%do_binning(:) = .true.
        mypdf%maxv(:,:) = 0D0
        mypdf%minv(:,:) = 0D0
        mypdf%bins(:,:) = 0D0
        mypdf%linthresh(:) = 0D0
        mypdf%zero_index(:) = 0
        mypdf%nvalues(:,:) = 0
        mypdf%nout(:,:) = 0
        mypdf%heights(:,:,:,:) = 0D0
        mypdf%totweights(:,:,:) = 0D0
        mypdf%total(:,:,:,:) = 0D0
    end subroutine allocate_pdf
end module stats_utils