!--------------------------------------------------------------------------
! ozymandias:profiles_module.f90
!--------------------------------------------------------------------------
!
! MODULE: amr_profiles
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> types and routines useful for the computation and storage of 1D and 2D
!> profiles using the AMR data of RAMSES.
!
!> @details  
!> defines profile_handler for 1D and profile_handler_twod, routines for
!> for the binning of cells and the computation of profiles
! 
!
!> @date 10/4/2021   0.1 fortran routines for ozymandias
!> @update 20/4/2021   0.1.1 adding 2D profiles
!--------------------------------------------------------------------------

module amr_profiles
    use local
    use io_ramses
    use filtering

    type profile_handler
        logical :: logscale
        logical :: cr_st=.false.,cr_heat=.false.
        integer :: profdim
        character(128) :: xvarname
        integer :: nyvar
        character(128),dimension(:),allocatable :: yvarnames
        integer :: nbins
        integer :: nwvar
        character(128),dimension(:),allocatable :: wvarnames
        real(dbl),dimension(:),allocatable :: xdata
        real(dbl),dimension(:,:,:,:),allocatable :: ydata
    end type profile_handler

    type profile_handler_twod
        logical :: cr_st=.false.,cr_heat=.false.
        integer :: profdim
        character(128) :: xvarname
        character(128) :: yvarname
        integer :: nzvar
        character(128),dimension(:),allocatable :: zvarnames
        integer,dimension(1:2) :: nbins
        integer :: nwvar
        character(128),dimension(:),allocatable :: wvarnames
        real(dbl),dimension(:),allocatable :: xdata,ydata
        real(dbl),dimension(:,:,:,:,:),allocatable :: zdata ! dimension(nx,ny,nz,nw,4)
    end type profile_handler_twod

    contains

    subroutine allocate_profile_handler(prof)
        implicit none
        type(profile_handler),intent(inout) :: prof

        if (.not.allocated(prof%yvarnames)) allocate(prof%yvarnames(prof%nyvar))
        if (.not.allocated(prof%wvarnames)) allocate(prof%wvarnames(prof%nwvar))
        if (.not.allocated(prof%xdata)) allocate(prof%xdata(0:prof%nbins))
        if (.not.allocated(prof%ydata)) allocate(prof%ydata(prof%nbins,prof%nyvar,prof%nwvar,4))

        if (prof%cr_st) sim%cr_st = .true.
        if (prof%cr_heat)sim%cr_heat = .true.
    end subroutine allocate_profile_handler

    subroutine allocate_profile_handler_twod(prof)
        implicit none
        type(profile_handler_twod),intent(inout) :: prof

        if (.not.allocated(prof%zvarnames)) allocate(prof%zvarnames(prof%nzvar))
        if (.not.allocated(prof%wvarnames)) allocate(prof%wvarnames(prof%nwvar))
        if (.not.allocated(prof%xdata)) allocate(prof%xdata(0:prof%nbins(1)))
        if (.not.allocated(prof%ydata)) allocate(prof%ydata(0:prof%nbins(2)))
        if (.not.allocated(prof%zdata)) allocate(prof%zdata(prof%nbins(1),prof%nbins(2),prof%nzvar,prof%nwvar,4))
    
        if (prof%cr_st) sim%cr_st = .true.
        if (prof%cr_heat)sim%cr_heat = .true.
    end subroutine allocate_profile_handler_twod

    subroutine makebins(reg,varname,nbins,bins,logscale)
        use geometrical_regions

        implicit none
        type(region),intent(in) :: reg
        character(128),intent(in) :: varname
        integer,intent(in) :: nbins
        real(dbl),dimension(0:nbins),intent(inout) :: bins
        logical,intent(in) :: logscale
        integer :: n
        real(dbl) :: temp_convfac,pressure_convfac,velocity_convfac
        real(dbl) :: rmin

        temp_convfac = ((sim%unit_l/sim%unit_t)**2)/1.38d-16*1.66d-24
        pressure_convfac = sim%unit_d*((sim%unit_l/sim%unit_t)**2)
        velocity_convfac = sim%unit_l/sim%unit_t

        select case (TRIM(varname))
        case('r_sphere','r_cyl')
            rmin = max(1D0/(2D0**(amr%nlevelmax-1)),1D-3*reg%rmax)
            do n=0,nbins
                if (logscale) then
                    if (reg%rmin.eq.0D0) then
                        bins(n) = dble(n)*(log10(reg%rmax)-log10(rmin))/dble(nbins) + log10(rmin)
                    else
                        bins(n) = dble(n)*(log10(reg%rmax)-log10(reg%rmin))/dble(nbins) + log10(reg%rmin)
                    end if
                else
                    bins(n) = dble(n)*(reg%rmax-reg%rmin)/dble(nbins) + reg%rmin
                endif
            end do
        case('z')
            do n=0,nbins
                if (logscale) then
                    bins(n) = dble(n)*(log10(reg%zmax)-log10(reg%zmin))/dble(nbins)
                    if (reg%zmin > 0D0) bins(n) = bins(n) + log10(reg%zmin)
                else
                    bins(n) = dble(n)*(reg%zmax-reg%zmin)/dble(nbins) + reg%zmin
                endif
            end do
        case('density')
            do n=0,nbins
                if (logscale) then
                    bins(n) = dble(n)*(log10(1D-20/sim%unit_d)-log10(1D-30/sim%unit_d))/dble(nbins) + log10(1D-30/sim%unit_d)
                else
                    bins(n) = dble(n)*(1D-20/sim%unit_d - 1D-30/sim%unit_d)/dble(nbins) + 1D-30/sim%unit_d
                endif
            end do
        case('temperature')
            do n=0,nbins
                if (logscale) then
                    bins(n) = dble(n)*(log10(1D8/temp_convfac)-log10(1D0/temp_convfac))/dble(nbins) + log10(1D0/temp_convfac)
                else
                    bins(n) = dble(n)*(1D8/temp_convfac - 1D0/temp_convfac)/dble(nbins) + 1D0/temp_convfac
                endif
            end do
        case('total_coolingtime')
            do n=0,nbins
                if (logscale) then
                    bins(n) = dble(n)*(log10(3.15D18/sim%unit_t)-log10(3.15D13/sim%unit_t))/dble(nbins) + log10(3.15D13/sim%unit_t)
                else
                    bins(n) = dble(n)*(3.15D18/sim%unit_t - 3.15D13/sim%unit_t)/dble(nbins) + 3.15D13/sim%unit_t
                endif
            end do
        case('thermal_pressure')
            do n=0,nbins
                if (logscale) then
                    bins(n) = dble(n)*(log10(1D-9/pressure_convfac)-log10(1D-16/pressure_convfac))/dble(nbins) + log10(1D-16/pressure_convfac)
                else
                    bins(n) = dble(n)*(1D-9/pressure_convfac - 1D-16/pressure_convfac)/dble(nbins) + 1D-16/pressure_convfac
                endif
            end do
        case('v_sphere_r')
            do n=0,nbins
                if (logscale) then
                    write(*,*)'You cannot use logscale for a velocity profile. Stopping!'
                    stop
                else
                    bins(n) = dble(n)*(4D7/velocity_convfac + 3D7/velocity_convfac)/dble(nbins) - 3D7/velocity_convfac
                endif
            end do
        case('theta_sphere')
            do n=0,nbins
                if (logscale) then
                    bins(n) = dble(n)*(log10(pi))/dble(nbins)
                else
                    bins(n) = dble(n)*(pi + 0D0)/dble(nbins) - 0D0
                endif
            end do
        !TODO: Add more cases
        end select
    end subroutine makebins

    subroutine findbinpos(reg,distance,pos,cellvars,cellsons,cellsize,prof,ibin,trans_matrix,grav_var)
        use vectors
        use geometrical_regions
        implicit none
        type(region),intent(in) :: reg
        real(dbl),dimension(1:3),intent(in) :: pos
        real(dbl),intent(in) :: distance
        real(dbl),dimension(0:amr%twondim,1:varIDs%nvar),intent(in) :: cellvars
        integer,dimension(0:amr%twondim),intent(in) :: cellsons
        real(dbl),intent(in) :: cellsize
        type(profile_handler),intent(in) :: prof
        integer,intent(inout) :: ibin
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix
        real(dbl),dimension(0:amr%twondim,1:4),optional,intent(in) :: grav_var
        real(dbl) :: value
        type(vector) :: x

        x = pos
        if (prof%xvarname.eq.reg%criteria_name) then
            value = distance
        else
            if (present(grav_var)) then
                call getvarvalue(reg,cellsize,x,cellvars,cellsons,prof%xvarname,value,trans_matrix,grav_var)
            else
                call getvarvalue(reg,cellsize,x,cellvars,cellsons,prof%xvarname,value,trans_matrix)
            end if
        endif
        if (prof%logscale) value = log10(value)

        ibin = int(dble(prof%nbins)*(value-prof%xdata(0))/(prof%xdata(prof%nbins)-prof%xdata(0))) + 1
        if (value .eq. prof%xdata(prof%nbins)) then
            ibin = prof%nbins
        else if (value<prof%xdata(0).or.value>prof%xdata(prof%nbins)) then
            ibin = 0
        end if
    end subroutine findbinpos

    subroutine findbinpos_twod(reg,distance,pos,cellvars,cellsons,cellsize,prof,logscale,ibinx,ibiny,trans_matrix,grav_var)
        use vectors
        use geometrical_regions
        implicit none
        type(region),intent(in) :: reg
        real(dbl),dimension(1:3),intent(in) :: pos
        real(dbl),intent(in) :: distance
        real(dbl),dimension(0:amr%twondim,1:varIDs%nvar),intent(in) :: cellvars
        integer,dimension(0:amr%twondim),intent(in) :: cellsons
        real(dbl),intent(in) :: cellsize
        type(profile_handler_twod),intent(in) :: prof
        logical,intent(in) :: logscale
        integer,intent(inout) :: ibinx,ibiny
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix
        real(dbl),dimension(0:amr%twondim,1:4),optional,intent(in) :: grav_var
        real(dbl) :: value
        type(vector) :: x

        x = pos
        if (prof%xvarname.eq.reg%criteria_name) then
            value = distance
        else
            if (present(grav_var)) then
                call getvarvalue(reg,cellsize,x,cellvars,cellsons,prof%xvarname,value,trans_matrix,grav_var)
            else
                call getvarvalue(reg,cellsize,x,cellvars,cellsons,prof%xvarname,value,trans_matrix)
            end if
        endif
        ! print*,value
        if (logscale) then
            if (value.le.0D0) then
                ibinx = 0
            else
                ibinx = int(dble(prof%nbins(1))*(log10(value)-prof%xdata(0))/(prof%xdata(prof%nbins(1))-prof%xdata(0))) + 1
                if (log10(value)<prof%xdata(0).or.log10(value)>prof%xdata(prof%nbins(1))) ibinx = 0
            end if
        else
            ibinx = int(dble(prof%nbins(1))*(value-prof%xdata(0))/(prof%xdata(prof%nbins(1))-prof%xdata(0))) + 1
            if (value<prof%xdata(0).or.value>prof%xdata(prof%nbins(1))) ibinx = 0
        endif

        if (prof%yvarname.eq.reg%criteria_name) then
            value = distance
        else
            if (present(grav_var)) then
                call getvarvalue(reg,cellsize,x,cellvars,cellsons,prof%yvarname,value,trans_matrix,grav_var)
            else
                call getvarvalue(reg,cellsize,x,cellvars,cellsons,prof%yvarname,value,trans_matrix)
            end if
        endif
        ! print*,value
        if (logscale) then
            if (value.le.0D0) then
                    ibiny = 0
            else
                ibiny = int(dble(prof%nbins(2))*(log10(value)-prof%ydata(0))/(prof%ydata(prof%nbins(2))-prof%ydata(0))) + 1
                if (log10(value)<prof%ydata(0).or.log10(value)>prof%ydata(prof%nbins(2))) ibiny = 0
            end if
        else            
            ibiny = int(dble(prof%nbins(2))*(value-prof%ydata(0))/(prof%ydata(prof%nbins(2))-prof%ydata(0))) + 1
            if (value<prof%ydata(0).or.value>prof%ydata(prof%nbins(2))) ibiny = 0
        endif
        ! print*,ibiny,ibinx,value
    end subroutine findbinpos_twod

    subroutine bindata(reg,pos,cellvars,cellsons,cellsize,prof,ibin,trans_matrix,grav_var)
        use vectors
        use geometrical_regions
        implicit none
        type(region),intent(in) :: reg
        real(dbl),dimension(1:3),intent(in) :: pos
        real(dbl),dimension(0:amr%twondim,1:varIDs%nvar),intent(in) :: cellvars
        integer,dimension(0:amr%twondim),intent(in) :: cellsons
        real(dbl),intent(in) :: cellsize
        type(profile_handler),intent(inout) :: prof
        integer,intent(in) :: ibin
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix
        real(dbl),dimension(0:amr%twondim,1:4),optional,intent(in) :: grav_var
        integer :: i,j
        real(dbl) :: ytemp,wtemp,bigwtemp,bigatemp
        type(vector) :: x
        x = pos
        yvarloop: do i=1,prof%nyvar
            if (present(grav_var)) then
                call getvarvalue(reg,cellsize,x,cellvars,cellsons,prof%yvarnames(i),ytemp,trans_matrix,grav_var)
            else
                call getvarvalue(reg,cellsize,x,cellvars,cellsons,prof%yvarnames(i),ytemp,trans_matrix)
            end if                
            wvarloop: do j=1,prof%nwvar
                if (prof%wvarnames(j)=='counts'.or.prof%wvarnames(j)=='cumulative') then
                    wtemp = 1D0
                else
                    call getvarvalue(reg,cellsize,x,cellvars,cellsons,prof%wvarnames(j),wtemp)
                endif
                ! Unbiased STD method. See: https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
                ! Q_k
                prof%ydata(ibin,i,j,1) = prof%ydata(ibin,i,j,1) + ytemp*wtemp

                ! W_k
                prof%ydata(ibin,i,j,2) = prof%ydata(ibin,i,j,2) + wtemp

                ! V_k
                prof%ydata(ibin,i,j,3) = prof%ydata(ibin,i,j,3) + wtemp**2

                !A_k
                prof%ydata(ibin,i,j,4) = prof%ydata(ibin,i,j,4) + wtemp*(ytemp**2)

            end do wvarloop
        end do yvarloop
    end subroutine bindata

    subroutine bindata_twod(reg,pos,cellvars,cellsons,cellsize,prof,ibinx,ibiny,trans_matrix,grav_var)
        use vectors
        use geometrical_regions
        implicit none
        type(region),intent(in) :: reg
        real(dbl),dimension(1:3),intent(in) :: pos
        real(dbl),dimension(0:amr%twondim,1:varIDs%nvar),intent(in) :: cellvars
        integer,dimension(0:amr%twondim),intent(in) :: cellsons
        real(dbl),intent(in) :: cellsize
        type(profile_handler_twod),intent(inout) :: prof
        integer,intent(in) :: ibinx,ibiny
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix
        real(dbl),dimension(0:amr%twondim,1:4),optional,intent(in) :: grav_var
        integer :: i,j
        real(dbl) :: ytemp,wtemp,bigwtemp,bigatemp
        type(vector) :: x
        x = pos
        zvarloop: do i=1,prof%nzvar
            if (present(grav_var)) then
                call getvarvalue(reg,cellsize,x,cellvars,cellsons,prof%zvarnames(i),ytemp,trans_matrix,grav_var)
            else
                call getvarvalue(reg,cellsize,x,cellvars,cellsons,prof%zvarnames(i),ytemp,trans_matrix)
            end if 
            wvarloop: do j=1,prof%nwvar
                if (prof%wvarnames(j)=='counts'.or.prof%wvarnames(j)=='cumulative') then
                    wtemp = 1D0
                else
                    call getvarvalue(reg,cellsize,x,cellvars,cellsons,prof%wvarnames(j),wtemp)
                endif
                ! Unbiased STD method. See: https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
                ! Q_k
                prof%zdata(ibinx,ibiny,i,j,1) = prof%zdata(ibinx,ibiny,i,j,1) + ytemp*wtemp

                ! W_k
                prof%zdata(ibinx,ibiny,i,j,2) = prof%zdata(ibinx,ibiny,i,j,2) + wtemp

                ! V_k
                prof%zdata(ibinx,ibiny,i,j,3) = prof%zdata(ibinx,ibiny,i,j,3) + wtemp**2

                !A_k
                prof%zdata(ibinx,ibiny,i,j,4) = prof%zdata(ibinx,ibiny,i,j,4) + wtemp*(ytemp**2)
            end do wvarloop
        end do zvarloop
    end subroutine bindata_twod

    subroutine renormalise_bins(prof_data)
        implicit none
        type(profile_handler),intent(inout) :: prof_data
        integer :: ibin,iy,iw
        real(dbl) :: Q_k,W_k,V_k,A_k

        binloop: do ibin=1,prof_data%nbins
            yloop: do iy=1,prof_data%nyvar
                wloop: do iw=1,prof_data%nwvar
                    Q_k = prof_data%ydata(ibin,iy,iw,1)
                    W_k = prof_data%ydata(ibin,iy,iw,2)
                    V_k = prof_data%ydata(ibin,iy,iw,3)
                    A_k = prof_data%ydata(ibin,iy,iw,4)
                    if (prof_data%wvarnames(iw) /= 'cumulative') then
                        ! Mean value or mean weighted value
                        prof_data%ydata(ibin,iy,iw,1) = Q_k / W_k
                        ! Standard deviation or weighted standard deviation
                        prof_data%ydata(ibin,iy,iw,2) = (A_k*W_k - Q_k**2) &
                                                        &/ (W_k**2 - V_k)
                        prof_data%ydata(ibin,iy,iw,2) = sqrt(prof_data%ydata(ibin,iy,iw,2))
                    endif
                end do wloop
            end do yloop
        end do binloop
    end subroutine renormalise_bins

    subroutine renormalise_bins_twod(prof_data)
        implicit none
        type(profile_handler_twod),intent(inout) :: prof_data
        integer :: ixbin,iybin,iz,iw
        real(dbl) :: Q_k,W_k,V_k,A_k
        xbinloop: do ixbin=1,prof_data%nbins(1)
                ybinloop: do iybin=1,prof_data%nbins(2)
                    zloop: do iz=1,prof_data%nzvar
                        wloop: do iw=1,prof_data%nwvar
                            Q_k = prof_data%zdata(ixbin,iybin,iz,iw,1)
                            W_k = prof_data%zdata(ixbin,iybin,iz,iw,2)
                            V_k = prof_data%zdata(ixbin,iybin,iz,iw,3)
                            A_k = prof_data%zdata(ixbin,iybin,iz,iw,4)
                            if (prof_data%wvarnames(iw) /= 'cumulative') then
                                ! Mean value or mean weighted value
                                prof_data%zdata(ixbin,iybin,iz,iw,1) = Q_k / W_k
                                ! Standard deviation or weighted standard deviation
                                prof_data%zdata(ixbin,iybin,iz,iw,2) = (A_k*W_k - Q_k**2) &
                                                                &/ (W_k**2 - V_k)
                                prof_data%zdata(ixbin,iybin,iz,iw,2) = sqrt(prof_data%zdata(ixbin,iybin,iz,iw,2))
                            endif
                        end do wloop
                    end do zloop
                end do ybinloop
            end do xbinloop
    end subroutine renormalise_bins_twod

    subroutine get_cells_onedprofile(repository,reg,filt,prof_data)
        use vectors
        use coordinate_systems
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(region), intent(in)  :: reg
        type(filter),intent(in) :: filt
        type(profile_handler),intent(inout) :: prof_data
        integer :: binpos
        logical :: ok_cell,ok_filter,read_gravity
        ! logical :: filter_cell
        integer :: i,j,k
        integer :: ipos,icpu,ilevel,ind,idim,ivar,iskip,inbor
        integer :: ix,iy,iz,ngrida,nx_full,ny_full,nz_full,cumngrida
        integer :: nvarh
        integer :: roterr
        character(5) :: nchar,ncharcpu
        character(128) :: nomfich
        real(dbl) :: distance,dx
        type(vector) :: xtemp,vtemp
        integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
        real(dbl),dimension(1:8,1:3) :: xc
        real(dbl),dimension(3,3) :: trans_matrix
        real(dbl),dimension(:,:),allocatable :: xg,x
        real(dbl),dimension(:,:),allocatable :: var
        real(dbl),dimension(:,:),allocatable :: grav_var
        real(dbl),dimension(:,:),allocatable :: tempvar
        real(dbl),dimension(:,:),allocatable :: tempgrav_var
        integer,dimension(:,:),allocatable :: nbor
        integer,dimension(:),allocatable :: son,tempson
        integer,dimension(:),allocatable :: ind_grid,ind_cell,ind_cell2
        integer ,dimension(:),allocatable :: ind_nbor
        logical,dimension(:),allocatable :: ref

        ! Check whether we need to read the gravity files
        read_gravity = .false.
        do ivar=1,prof_data%nyvar
            if (prof_data%yvarnames(ivar)(1:4) .eq. 'grav') then
                read_gravity = .true.
                write(*,*)'Reading gravity files...'
                exit
            endif
        end do

        allocate(ngridfile(1:amr%ncpu+amr%nboundary,1:amr%nlevelmax))
        allocate(ngridlevel(1:amr%ncpu,1:amr%nlevelmax))
        if(amr%nboundary>0)allocate(ngridbound(1:amr%nboundary,1:amr%nlevelmax))

        trans_matrix = 0D0
        call new_z_coordinates(reg%axis,trans_matrix,roterr)
        if (roterr.eq.1) then
            write(*,*) 'Incorrect CS transformation!'
            stop
        endif
        ipos=INDEX(repository,'output_')
        nchar=repository(ipos+7:ipos+13)
        ! Loop over processor files
        cpuloop: do k=1,amr%ncpu_read
            icpu = amr%cpu_list(k)
            call title(icpu,ncharcpu)

            ! Open AMR file and skip header
            nomfich = TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
            open(unit=10,file=nomfich,status='old',form='unformatted')
            do i=1,21
                read(10) ! Skip header
            end do
            ! Read grid numbers
            read(10)ngridlevel
            ngridfile(1:amr%ncpu,1:amr%nlevelmax) = ngridlevel
            read(10) ! Skip
            if(amr%nboundary>0) then
                do i=1,2
                    read(10)
                end do
                read(10)ngridbound
                ngridfile(amr%ncpu+1:amr%ncpu+amr%nboundary,1:amr%nlevelmax) = ngridbound
            endif
            read(10) ! Skip
            ! R. Teyssier: comment the single following line for old stuff
            read(10)
            if(TRIM(amr%ordering).eq.'bisection')then
                do i=1,5
                    read(10)
                end do
            else
                read(10)
            endif
            read(10)
            read(10)
            read(10)

            ! Make sure that we are not trying to access too far in the refinement map…
            ! call check_lmax(ngridfile)
            ! write(*,*)'active_lmax ', amr%active_lmax

            allocate(nbor(1:amr%ngridmax,1:amr%twondim))
            allocate(son(1:amr%ncoarse+amr%twotondim*amr%ngridmax))
            son = 0; nbor = 0
            cumngrida = 0

            ! Open HYDRO file and skip header
            nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
            open(unit=11,file=nomfich,status='old',form='unformatted')
            read(11)
            read(11)nvarh
            read(11)
            read(11)
            read(11)
            read(11)

            allocate(var(1:amr%ncoarse+amr%twotondim*amr%ngridmax,1:nvarh))
            
            if (read_gravity) then
                ! Open GRAV file and skip header
                nomfich=TRIM(repository)//'/grav_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
                open(unit=12,file=nomfich,status='old',form='unformatted')
                read(12) !ncpu
                read(12) !ndim
                read(12) !nlevelmax
                read(12) !nboundary 
                allocate(grav_var(1:amr%ncoarse+amr%twotondim*amr%ngridmax,1:4))
            endif

            ! Loop over levels
            levelloop: do ilevel=1,amr%lmax
                ! Geometry
                dx = 0.5**ilevel
                nx_full = 2**ilevel
                ny_full = 2**ilevel
                nz_full = 2**ilevel
                do ind=1,amr%twotondim
                    iz=(ind-1)/4
                    iy=(ind-1-4*iz)/2
                    ix=(ind-1-2*iy-4*iz)
                    xc(ind,1)=(dble(ix)-0.5D0)*dx
                    xc(ind,2)=(dble(iy)-0.5D0)*dx
                    xc(ind,3)=(dble(iz)-0.5D0)*dx
                end do
                
                ! Allocate work arrays
                ngrida = ngridfile(icpu,ilevel)
                if(ngrida>0) then
                    allocate(ind_grid(1:ngrida))
                    allocate(ind_cell(1:ngrida))
                    allocate(xg (1:ngrida,1:amr%ndim))
                    allocate(x  (1:ngrida,1:amr%ndim))
                    allocate(ref(1:ngrida))
                endif

                ! Loop over domains
                domloop: do j=1,amr%nboundary+amr%ncpu
                    
                    ! Read AMR data
                    if (ngridfile(j,ilevel)>0) then
                        if(j.eq.icpu)then
                            read(10) ind_grid
                        else
                            read(10)
                        end if
                        read(10) ! Skip next index
                        read(10) ! Skip prev index

                        ! Read grid center
                        do idim=1,amr%ndim
                            if(j.eq.icpu)then
                                read(10)xg(:,idim)
                            else
                                read(10)
                            endif
                        end do

                        read(10) ! Skip father index
                        do ind=1,amr%twondim
                            if(j.eq.icpu)then
                                read(10)nbor(ind_grid,ind)
                            else
                                read(10)
                            end if
                        end do

                        ! Read son index
                        do ind=1,amr%twotondim
                            iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                            if(j.eq.icpu)then
                                read(10)son(ind_grid+iskip)
                            else
                                read(10)
                            end if
                        end do

                        ! Skip cpu map
                        do ind=1,amr%twotondim
                            read(10)
                        end do

                        ! Skip refinement map
                        do ind=1,amr%twotondim
                            read(10)
                        end do
                    endif

                    ! Read HYDRO data
                    read(11)
                    read(11)
                    if(ngridfile(j,ilevel)>0)then
                        ! Read hydro variables
                        tndimloop: do ind=1,amr%twotondim
                            iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                            varloop: do ivar=1,nvarh
                                if (j.eq.icpu) then
                                    read(11)var(ind_grid+iskip,ivar)
                                else
                                    read(11)
                                endif
                            end do varloop
                        end do tndimloop
                    endif

                    if (read_gravity) then
                        ! Read GRAV data
                        read(12)
                        read(12)
                        if(ngridfile(j,ilevel)>0)then
                            do ind=1,amr%twotondim
                                iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                                if (j.eq.icpu) then
                                    read(12)grav_var(ind_grid+iskip,1)
                                    do ivar=1,amr%ndim
                                        read(12)grav_var(ind_grid+iskip,ivar+1)
                                    end do
                                else
                                    read(12)
                                    do ivar=1,amr%ndim
                                        read(12)
                                    end do
                                end if
                            end do
                        end if
                    end if
                end do domloop

                !Compute map
                if (ngrida>0) then
                    ! Loop over cells
                    cellloop: do ind=1,amr%twotondim

                        ! Compute cell center
                        do i=1,ngrida
                            x(i,1)=(xg(i,1)+xc(ind,1)-amr%xbound(1))
                            x(i,2)=(xg(i,2)+xc(ind,2)-amr%xbound(2))
                            x(i,3)=(xg(i,3)+xc(ind,3)-amr%xbound(3))
                        end do

                        ! Check if cell is refined
                        iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                        do i=1,ngrida
                            ref(i) = son(ind_grid(i)+iskip)>0.and.ilevel<amr%nlevelmax
                        end do
                        
                        ! Get cell indexes
                        do i=1,ngrida
                            ind_cell(i) = iskip+ind_grid(i)
                        end do

                        ngridaloop: do i=1,ngrida
                            ! Check if cell is inside the desired region
                            distance = 0D0
                            xtemp = x(i,:)
                            xtemp = xtemp - reg%centre
                            call rotate_vector(xtemp,trans_matrix)
                            x(i,:) = xtemp
                            call checkifinside(x(i,:),reg,ok_cell,distance)
                            
                            ! Velocity transformed --> ONLY FOR CENTRAL CELL
                            vtemp = var(ind_cell(i),varIDs%vx:varIDs%vz)
                            vtemp = vtemp - reg%bulk_velocity
                            call rotate_vector(vtemp,trans_matrix)
                            var(ind_cell(i),varIDs%vx:varIDs%vz) = vtemp

                            ! Gravitational acc --> ONLY FOR CENTRAL CELL
                            if (read_gravity) then
                                vtemp = grav_var(ind_cell(i),2:4)
                                call rotate_vector(vtemp,trans_matrix)
                                grav_var(ind_cell(i),2:4) = vtemp
                            endif

                            ! Get neighbours
                            allocate(ind_cell2(1))
                            ind_cell2(1) = ind_cell(i)
                            allocate(ind_nbor(0:amr%twondim))
                            call getnbor(son,nbor,ind_cell2,ind_nbor,1)
                            deallocate(ind_cell2)
                            allocate(tempvar(0:amr%twondim,nvarh))
                            allocate(tempson(0:amr%twondim))
                            if (read_gravity) allocate(tempgrav_var(0:amr%twondim,1:4))
                            do inbor=0,amr%twondim
                                tempvar(inbor,:) = var(ind_nbor(inbor),:)
                                tempson(inbor)       = son(ind_nbor(inbor))
                                if (read_gravity) tempgrav_var(inbor,:) = grav_var(ind_nbor(inbor),:)
                            end do
                            deallocate(ind_nbor)

                            ok_filter = filter_cell(reg,filt,xtemp,dx,tempvar,tempson)
                            ok_cell= ok_cell.and..not.ref(i).and.ok_filter
                            if (ok_cell) then
                                binpos = 0
                                if (read_gravity) then
                                    call findbinpos(reg,distance,x(i,:),tempvar,tempson,dx,prof_data,binpos,trans_matrix,tempgrav_var)
                                    if (binpos.ne.0) call bindata(reg,x(i,:),tempvar,tempson,dx,prof_data,binpos,trans_matrix,tempgrav_var)
                                else
                                    call findbinpos(reg,distance,x(i,:),tempvar,tempson,dx,prof_data,binpos,trans_matrix)
                                    if (binpos.ne.0) call bindata(reg,x(i,:),tempvar,tempson,dx,prof_data,binpos,trans_matrix)
                                end if
                            endif
                            deallocate(tempvar,tempson)
                            if (read_gravity) deallocate(tempgrav_var)
                        end do ngridaloop
                    end do cellloop
                    deallocate(xg,ref,x,ind_grid,ind_cell)
                endif
                cumngrida = cumngrida + ngrida
            end do levelloop
            deallocate(nbor,son,var)
            close(10)
            close(11)
            if (read_gravity) then
                close(12)
                deallocate(grav_var)
            end if
        end do cpuloop
    end subroutine get_cells_onedprofile

    subroutine onedprofile(repository,reg,filt,prof_data,lmax,logscale)
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        type(filter),intent(in) :: filt
        type(profile_handler),intent(inout) :: prof_data
        integer,intent(in) :: lmax
        logical,intent(in) :: logscale

        call read_hydrofile_descriptor(repository)

        call init_amr_read(repository)
        amr%lmax = lmax
        if (lmax.eq.0) amr%lmax = amr%nlevelmax
        prof_data%xdata = 0D0
        prof_data%ydata = 0D0
        prof_data%logscale = logscale
        call makebins(reg,prof_data%xvarname,prof_data%nbins,prof_data%xdata,logscale)
        
        call get_cpu_map(reg)
        write(*,*)'ncpu_read:',amr%ncpu_read
        call get_cells_onedprofile(repository,reg,filt,prof_data)

        call renormalise_bins(prof_data)

        if (logscale) prof_data%xdata = 10.**(prof_data%xdata)

    end subroutine onedprofile

    subroutine twodprofile(repository,reg,filt,prof_data,lmax,logscale)
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        type(filter),intent(in) :: filt
        type(profile_handler_twod),intent(inout) :: prof_data
        integer,intent(in) :: lmax
        logical,intent(in) :: logscale

        call read_hydrofile_descriptor(repository)
        call init_amr_read(repository)
        amr%lmax = lmax
        if (lmax.eq.0) amr%lmax = amr%nlevelmax
        
        prof_data%xdata = 0D0
        prof_data%ydata = 0D0
        prof_data%zdata = 0D0

        call makebins(reg,prof_data%xvarname,prof_data%nbins(1),prof_data%xdata,logscale)
        call makebins(reg,prof_data%yvarname,prof_data%nbins(2),prof_data%ydata,logscale)
        write(*,*)'lmax: ',amr%lmax
        call get_cpu_map(reg)
        write(*,*)'ncpu_read:',amr%ncpu_read
        call get_cells_twodprofile(repository,reg,filt,prof_data,logscale)
        call renormalise_bins_twod(prof_data)

    end subroutine twodprofile

    subroutine get_cells_twodprofile(repository,reg,filt,prof_data,logscale)
        use vectors
        use coordinate_systems
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(region), intent(in)  :: reg
        type(filter),intent(in) :: filt
        type(profile_handler_twod),intent(inout) :: prof_data
        logical,intent(in) :: logscale
        integer :: xbinpos,ybinpos
        logical :: ok_cell,ok_filter,read_gravity
        integer :: i,j,k
        integer :: ipos,icpu,ilevel,ind,idim,ivar,iskip,inbor
        integer :: ix,iy,iz,ngrida,nx_full,ny_full,nz_full,total_ncell
        integer :: nvarh
        integer :: roterr
        character(5) :: nchar,ncharcpu
        character(128) :: nomfich
        real(dbl) :: distance,dx
        type(vector) :: xtemp,vtemp
        integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
        real(dbl),dimension(1:8,1:3) :: xc
        real(dbl),dimension(3,3) :: trans_matrix
        real(dbl),dimension(:,:),allocatable :: xg,x
        real(dbl),dimension(:,:),allocatable :: var
        real(dbl),dimension(:,:),allocatable :: grav_var
        real(dbl),dimension(:,:),allocatable :: tempvar
        real(dbl),dimension(:,:),allocatable :: tempgrav_var
        integer,dimension(:,:),allocatable :: nbor
        integer,dimension(:),allocatable :: son,tempson
        integer,dimension(:),allocatable :: ind_grid,ind_cell,ind_cell2
        integer ,dimension(:),allocatable :: ind_nbor
        logical,dimension(:),allocatable :: ref

        total_ncell = 0
        
        ! Check whether we need to read the gravity files
        read_gravity = .false.
        do ivar=1,prof_data%nzvar
            if (prof_data%zvarnames(ivar)(1:4) .eq. 'grav') then
                read_gravity = .true.
                write(*,*)'Reading gravity files...'
                exit
            endif
        end do

        allocate(ngridfile(1:amr%ncpu+amr%nboundary,1:amr%nlevelmax))
        allocate(ngridlevel(1:amr%ncpu,1:amr%nlevelmax))
        if(amr%nboundary>0)allocate(ngridbound(1:amr%nboundary,1:amr%nlevelmax))

        trans_matrix = 0D0
        call new_z_coordinates(reg%axis,trans_matrix,roterr)
        if (roterr.eq.1) then
            write(*,*) 'Incorrect CS transformation!'
            stop
        endif

        ipos=INDEX(repository,'output_')
        nchar=repository(ipos+7:ipos+13)
        ! Loop over processor files
        cpuloop: do k=1,amr%ncpu_read
            icpu = amr%cpu_list(k)
            call title(icpu,ncharcpu)

            ! Open AMR file and skip header
            nomfich = TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
            open(unit=10,file=nomfich,status='old',form='unformatted')
            ! write(*,*)'Processing file '//TRIM(nomfich)
            do i=1,21
                read(10) ! Skip header
            end do
            ! Read grid numbers
            read(10)ngridlevel
            ngridfile(1:amr%ncpu,1:amr%nlevelmax) = ngridlevel
            read(10) ! Skip
            if(amr%nboundary>0) then
                do i=1,2
                    read(10)
                end do
                read(10)ngridbound
                ngridfile(amr%ncpu+1:amr%ncpu+amr%nboundary,1:amr%nlevelmax) = ngridbound
            endif
            read(10) ! Skip
            ! R. Teyssier: comment the single following line for old stuff
            read(10)
            if(TRIM(amr%ordering).eq.'bisection')then
                do i=1,5
                    read(10)
                end do
            else
                read(10)
            endif
            read(10)
            read(10)
            read(10)

            ! Make sure that we are not trying to access to far in the refinement map…
            ! call check_lmax(ngridfile)

            allocate(nbor(1:amr%ngridmax,1:amr%twondim))
            allocate(son(1:amr%ncoarse+amr%twotondim*amr%ngridmax))
            son = 0; nbor = 0

            ! Open HYDRO file and skip header
            nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
            open(unit=11,file=nomfich,status='old',form='unformatted')
            read(11)
            read(11)nvarh
            read(11)
            read(11)
            read(11)
            read(11)

            allocate(var(1:amr%ncoarse+amr%twotondim*amr%ngridmax,1:nvarh))

            if (read_gravity) then
                ! Open GRAV file and skip header
                nomfich=TRIM(repository)//'/grav_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
                open(unit=12,file=nomfich,status='old',form='unformatted')
                read(12) !ncpu
                read(12) !ndim
                read(12) !nlevelmax
                read(12) !nboundary 
                allocate(grav_var(1:amr%ncoarse+amr%twotondim*amr%ngridmax,1:4))
            endif

            ! Loop over levels
            levelloop: do ilevel=1,amr%lmax
                ! Geometry
                dx = 0.5**ilevel
                nx_full = 2**ilevel
                ny_full = 2**ilevel
                nz_full = 2**ilevel
                do ind=1,amr%twotondim
                    iz=(ind-1)/4
                    iy=(ind-1-4*iz)/2
                    ix=(ind-1-2*iy-4*iz)
                    xc(ind,1)=(dble(ix)-0.5D0)*dx
                    xc(ind,2)=(dble(iy)-0.5D0)*dx
                    xc(ind,3)=(dble(iz)-0.5D0)*dx
                end do
                ! Allocate work arrays
                ngrida = ngridfile(icpu,ilevel)
                if(ngrida>0) then
                    allocate(ind_grid(1:ngrida))
                    allocate(ind_cell(1:ngrida))
                    allocate(xg (1:ngrida,1:amr%ndim))
                    allocate(x  (1:ngrida,1:amr%ndim))
                    allocate(ref(1:ngrida))
                endif
                ! Loop over domains
                domloop: do j=1,amr%nboundary+amr%ncpu
                    ! Read AMR data
                    if (ngridfile(j,ilevel)>0) then
                        if(j.eq.icpu)then
                            read(10) ind_grid
                        else
                            read(10)
                        end if
                        read(10) ! Skip next index
                        read(10) ! Skip prev index

                        ! Read grid center
                        do idim=1,amr%ndim
                            if(j.eq.icpu)then
                                read(10)xg(:,idim)
                            else
                                read(10)
                            endif
                        end do

                        read(10) ! Skip father index
                        do ind=1,amr%twondim
                            if(j.eq.icpu)then
                                read(10)nbor(ind_grid,ind)
                            else
                                read(10)
                            end if
                        end do

                        ! Read son index
                        do ind=1,amr%twotondim
                            iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                            if(j.eq.icpu)then
                                read(10)son(ind_grid+iskip)
                            else
                                read(10)
                            end if
                        end do

                        ! Skip cpu map
                        do ind=1,amr%twotondim
                            read(10)
                        end do

                        ! Skip refinement map
                        do ind=1,amr%twotondim
                            read(10)
                        end do
                    endif
                    ! Read HYDRO data
                    read(11)
                    read(11)
                    if(ngridfile(j,ilevel)>0)then
                        ! Read hydro variables
                        tndimloop: do ind=1,amr%twotondim
                            iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                            varloop: do ivar=1,nvarh
                                if (j.eq.icpu) then
                                    read(11)var(ind_grid+iskip,ivar)
                                else
                                    read(11)
                                endif
                            end do varloop
                        end do tndimloop
                    endif

                    if (read_gravity) then
                        ! Read GRAV data
                        read(12)
                        read(12)
                        if(ngridfile(j,ilevel)>0)then
                            do ind=1,amr%twotondim
                                iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                                if (j.eq.icpu) then
                                    read(12)grav_var(ind_grid+iskip,1)
                                    do ivar=1,amr%ndim
                                        read(12)grav_var(ind_grid+iskip,ivar+1)
                                    end do
                                else
                                    read(12)
                                    do ivar=1,amr%ndim
                                        read(12)
                                    end do
                                end if
                            end do
                        end if
                    end if
                end do domloop
                !Compute map
                if (ngrida>0) then
                    ! Loop over cells
                    cellloop: do ind=1,amr%twotondim
                        ! Compute cell center
                        do i=1,ngrida
                            x(i,1)=(xg(i,1)+xc(ind,1)-amr%xbound(1))
                            x(i,2)=(xg(i,2)+xc(ind,2)-amr%xbound(2))
                            x(i,3)=(xg(i,3)+xc(ind,3)-amr%xbound(3))
                        end do

                        ! Check if cell is refined
                        iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                        do i=1,ngrida
                            ref(i) = son(ind_grid(i)+iskip)>0.and.ilevel<amr%lmax
                        end do
                        
                        ! Get cell indexes
                        do i=1,ngrida
                            ind_cell(i) = iskip+ind_grid(i)
                        end do
                        ngridaloop: do i=1,ngrida
                            ! Check if cell is inside the desired region
                            distance = 0D0
                            xtemp = x(i,:)
                            xtemp = xtemp - reg%centre
                            call rotate_vector(xtemp,trans_matrix)
                            x(i,:) = xtemp
                            call checkifinside(x(i,:),reg,ok_cell,distance)

                            ! Velocity transformed --> ONLY FOR CENTRAL CELL
                            vtemp = var(ind_cell(i),varIDs%vx:varIDs%vz)
                            vtemp = vtemp - reg%bulk_velocity
                            call rotate_vector(vtemp,trans_matrix)
                            var(ind_cell(i),varIDs%vx:varIDs%vz) = vtemp

                            ! Gravitational acc --> ONLY FOR CENTRAL CELL
                            if (read_gravity) then
                                vtemp = grav_var(ind_cell(i),2:4)
                                call rotate_vector(vtemp,trans_matrix)
                                grav_var(ind_cell(i),2:4) = vtemp
                            endif

                            ! Get neighbours
                            allocate(ind_cell2(1))
                            ind_cell2(1) = ind_cell(i)
                            allocate(ind_nbor(0:amr%twondim))
                            call getnbor(son,nbor,ind_cell2,ind_nbor,1)
                            deallocate(ind_cell2)
                            allocate(tempvar(0:amr%twondim,nvarh))
                            allocate(tempson(0:amr%twondim))
                            if (read_gravity) allocate(tempgrav_var(0:amr%twondim,1:4))
                            do inbor=0,amr%twondim
                                tempvar(inbor,:) = var(ind_nbor(inbor),:)
                                tempson(inbor)       = son(ind_nbor(inbor))
                                if (read_gravity) tempgrav_var(inbor,:) = grav_var(ind_nbor(inbor),:)
                            end do
                            deallocate(ind_nbor)

                            ok_filter = filter_cell(reg,filt,xtemp,dx,tempvar,tempson)
                            ok_cell= ok_cell.and..not.ref(i).and.ok_filter
                            if (ok_cell) then
                                xbinpos = 0; ybinpos=0
                                total_ncell = total_ncell + 1
                                if (read_gravity) then
                                    call findbinpos_twod(reg,distance,x(i,:),tempvar,tempson,&
                                                        &dx,prof_data,logscale,xbinpos,ybinpos,&
                                                        &trans_matrix,tempgrav_var)
                                    if (xbinpos.ne.0.and.ybinpos.ne.0) call bindata_twod(reg,x(i,:),&
                                                                            &tempvar,tempson,dx,prof_data,&
                                                                            &xbinpos,ybinpos,&
                                                                            &trans_matrix,tempgrav_var)
                                else
                                    call findbinpos_twod(reg,distance,x(i,:),tempvar,tempson,&
                                                    &dx,prof_data,logscale,xbinpos,ybinpos,trans_matrix)
                                    if (xbinpos.ne.0.and.ybinpos.ne.0) call bindata_twod(reg,x(i,:),&
                                                                            &tempvar,tempson,dx,prof_data,&
                                                                            &xbinpos,ybinpos,trans_matrix)
                                end if
                            endif
                            deallocate(tempvar,tempson)
                            if (read_gravity) deallocate(tempgrav_var)
                        end do ngridaloop
                    end do cellloop
                    deallocate(xg,ref,x,ind_grid,ind_cell)
                endif
            end do levelloop
            deallocate(nbor,son,var)
            close(10)
            close(11)
            if (read_gravity) then
                close(12)
                deallocate(grav_var)
            end if
        end do cpuloop
        write(*,*)'Total number of cells used: ', total_ncell
    end subroutine get_cells_twodprofile
end module amr_profiles