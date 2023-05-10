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
    use geometrical_regions
    use stats_utils

    type profile_handler
        character(128) :: scaletype
        logical :: cr_st=.false.,cr_heat=.false.
        integer :: profdim
        character(128) :: xvarname
        integer :: nyvar
        character(128),dimension(:),allocatable :: yvarnames
        integer :: nbins
        integer :: nwvar
        integer :: nsubs=0
        character(128),dimension(:),allocatable :: wvarnames
        real(dbl) :: Dcr = 3D28
        real(dbl),dimension(:),allocatable :: xdata
        real(dbl),dimension(:,:,:,:),allocatable :: ydata
        type(region),dimension(:),allocatable :: subs
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
        real(dbl) :: Dcr = 3D28
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
        if (.not.allocated(prof%subs).and.(prof%nsubs>0)) allocate(prof%subs(1:prof%nsubs))
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

    subroutine findbinpos_twod(reg,distance,pos,cellvars,cellsons,cellsize,prof,scaletype,ibinx,ibiny,trans_matrix,grav_var)
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
        character(128),intent(in) :: scaletype
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
        if (trim(scaletype).eq.'log_even') then
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
        if (trim(scaletype).eq.'log_even') then
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

    subroutine onedprofile(repository,reg,filt,prof_data,lmax,scaletype,use_neigh)
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        type(filter),intent(in) :: filt
        type(profile_handler),intent(inout) :: prof_data
        integer,intent(in) :: lmax
        character(128),intent(in) :: scaletype
        logical, intent(in) :: use_neigh

        call read_hydrofile_descriptor(repository)

        call init_amr_read(repository)
        amr%lmax = lmax
        if (lmax.eq.0) amr%lmax = amr%nlevelmax
        prof_data%xdata = 0D0
        prof_data%ydata = 0D0
        prof_data%scaletype = scaletype
        prof_data%xdata = makebins(reg,prof_data%xvarname,prof_data%nbins,scaletype)
        
        call get_cpu_map(reg)
        write(*,*)'ncpu_read:',amr%ncpu_read
        ! Choose type of onedprofile
        if (use_neigh) then
            write(*,*)'Loading neighbours...'
            call get_cells_onedprofile_neigh
        else
            call get_cells_onedprofile_fast
        end if

        call renormalise_bins(prof_data)

        if (trim(scaletype).eq.'log_even') prof_data%xdata = 10.**(prof_data%xdata)

        contains

        subroutine get_cells_onedprofile_neigh
            use vectors
            use coordinate_systems
            implicit none
            integer :: binpos
            logical :: ok_cell,ok_filter,read_gravity,ok_sub
            ! logical :: filter_cell
            integer :: i,j,k
            integer :: ipos,icpu,ilevel,ind,idim,ivar,iskip,inbor,ison,isub
            integer :: ix,iy,iz,ngrida,nx_full,ny_full,nz_full
            integer :: total_ncell
            integer :: nvarh
            integer :: roterr
            character(5) :: nchar,ncharcpu
            character(128) :: nomfich
            real(dbl) :: distance,dx,ytemp
            type(vector) :: xtemp,vtemp,gtemp
            integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
            real(dbl),dimension(:),allocatable :: xxg,son_dens
            real(dbl),dimension(1:8,1:3) :: xc
            real(dbl),dimension(1:3,1:3) :: trans_matrix
            real(dbl),dimension(:,:),allocatable :: x,xorig
            real(dbl),dimension(:,:),allocatable :: var
            real(dbl),dimension(:,:),allocatable :: grav_var
            real(dbl),dimension(:,:),allocatable :: tempvar
            real(dbl),dimension(:,:),allocatable :: tempgrav_var
            real(dbl),dimension(:,:),allocatable :: cellpos
            integer,dimension(:,:),allocatable :: nbor
            integer,dimension(:),allocatable :: son,tempson,iig
            integer,dimension(:),allocatable :: ind_cell,ind_cell2
            integer ,dimension(1,0:amr%twondim) :: ind_nbor
            logical,dimension(:),allocatable :: ref
            type(level),dimension(1:100) :: grid

            total_ncell = 0
            ! Check whether we need to read the gravity files
            read_gravity = .false.
            do ivar=1,prof_data%nyvar
                if (prof_data%yvarnames(ivar)(1:4) .eq. 'grav' .or.&
                & trim(prof_data%yvarnames(ivar)) .eq. 'neighbour_accuracy') then
                    read_gravity = .true.
                    write(*,*)'Reading gravity files...'
                    exit
                endif
            end do

            allocate(ngridfile(1:amr%ncpu+amr%nboundary,1:amr%nlevelmax))
            allocate(ngridlevel(1:amr%ncpu,1:amr%nlevelmax))
            if(amr%nboundary>0)allocate(ngridbound(1:amr%nboundary,1:amr%nlevelmax))
            ! Compute hierarchy
            do ilevel=1,amr%lmax
                grid(ilevel)%ngrid = 0
            end do

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

                allocate(nbor(1:amr%ngridmax,1:amr%twondim))
                allocate(son(1:amr%ncoarse+amr%twotondim*amr%ngridmax))
                nbor = 0
                son = 0

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
                read(10)son(1:amr%ncoarse)
                read(10)
                read(10)

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
                allocate(cellpos(1:amr%ncoarse+amr%twotondim*amr%ngridmax,1:3))
                cellpos = 0d0
                var = 0d0
                
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
                levelloop1: do ilevel=1,amr%lmax
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
                    grid(ilevel)%ngrid = 0
                    ! Loop over domains
                    domloop: do j=1,amr%nboundary+amr%ncpu
                        ! Allocate work arrays
                        ngrida = ngridfile(j,ilevel)
                        if(ngrida>0)then
                            if (allocated(grid(ilevel)%ind_grid)) deallocate(grid(ilevel)%ind_grid)
                            if (allocated(grid(ilevel)%xg)) deallocate(grid(ilevel)%xg)
                            allocate(grid(ilevel)%ind_grid(1:ngrida))
                            allocate(grid(ilevel)%xg (1:ngrida,1:amr%ndim))
                            allocate(iig(1:ngrida))
                            allocate(xxg(1:ngrida))
                            
                            
                            ! Read AMR data
                            read(10) grid(ilevel)%ind_grid
                            read(10) ! Skip next index
                            read(10) ! Skip prev index
                            if(j.eq.icpu) then
                                if (allocated(grid(ilevel)%real_ind)) deallocate(grid(ilevel)%real_ind)
                                allocate(grid(ilevel)%real_ind(1:ngrida))
                                grid(ilevel)%real_ind = grid(ilevel)%ind_grid
                                grid(ilevel)%ngrid = ngridfile(j,ilevel)
                            end if

                            
                            ! Read grid center
                            do idim=1,amr%ndim
                                read(10)xxg
                                grid(ilevel)%xg(:,idim) = xxg(:)
                            end do
                            
                            read(10) ! Skip father index
                            ! Read nbor index
                            do ind=1,amr%twondim
                                read(10)iig
                                nbor(grid(ilevel)%ind_grid(:),ind) = iig(:)
                            end do
                            ! Read son index
                            do ind=1,amr%twotondim
                                iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                                read(10)iig
                                son(grid(ilevel)%ind_grid(:)+iskip) = iig(:)
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
                        if(ngrida>0)then
                            ! Read hydro variables
                            tndimloop: do ind=1,amr%twotondim
                                iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                                varloop: do ivar=1,nvarh
                                    read(11)xxg
                                    var(grid(ilevel)%ind_grid(:)+iskip,ivar) = xxg(:)
                                end do varloop
                            end do tndimloop
                        endif

                        if (read_gravity) then
                            ! Read GRAV data
                            read(12)
                            read(12)
                            if(ngrida>0)then
                                do ind=1,amr%twotondim
                                    iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                                    read(12)xxg
                                    grav_var(grid(ilevel)%ind_grid(:)+iskip,1) = xxg(:)
                                    do ivar=1,amr%ndim
                                        read(12)xxg
                                        grav_var(grid(ilevel)%ind_grid(:)+iskip,ivar+1) = xxg(:)
                                    end do
                                end do
                            end if
                        end if

                        !Compute positions
                        if(ngrida>0)then
                            do ind=1,amr%twotondim
                                iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                                do i=1,ngrida
                                    do ivar=1,amr%ndim
                                        cellpos(grid(ilevel)%ind_grid(i)+iskip,ivar)=(grid(ilevel)%xg(i,ivar)+xc(ind,ivar)-amr%xbound(ivar))
                                    end do
                                end do
                            end do
                        end if
                        if (ngrida>0) deallocate(iig,xxg)

                    end do domloop
                end do levelloop1
                close(10)
                close(11)
                if (read_gravity) then
                    close(12)
                end if
                ! Loop over levels again now with arrays fully filled
                levelloop2: do ilevel=1,amr%lmax
                    ! Geometry
                    dx = 0.5**ilevel
                    nx_full = 2**ilevel
                    ny_full = 2**ilevel

                    ! Allocate work arrays
                    ngrida = grid(ilevel)%ngrid
                    if(ngrida>0)then
                        allocate(ind_cell(1:ngrida))
                        allocate(x  (1:ngrida,1:amr%ndim))
                        allocate(xorig(1:ngrida,1:amr%ndim))
                        allocate(ref(1:ngrida))
                    endif

                    !Compute map
                    if (ngrida>0) then
                        ! Loop over cells
                        cellloop: do ind=1,amr%twotondim

                            ! Get cell indexes
                            iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                            do i=1,ngrida
                                ind_cell(i) = iskip+grid(ilevel)%real_ind(i)
                            end do
                            ! Compute cell center
                            do i=1,ngrida
                                x(i,:)=cellpos(ind_cell(i),:)
                            end do

                            ! Check if cell is refined
                            do i=1,ngrida
                                ref(i) = son(ind_cell(i))>0.and.ilevel<amr%lmax
                            end do          
                            xorig = x
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

                                ! Gravitational acc --> ONLY FOR CENTRAL CELL
                                if (read_gravity) then
                                    gtemp = grav_var(ind_cell(i),2:4)
                                    call rotate_vector(gtemp,trans_matrix)
                                endif

                                ! Get neighbours
                                allocate(ind_cell2(1))
                                ind_cell2(1) = ind_cell(i)
                                call getnbor(son,nbor,ind_cell2,ind_nbor,1)
                                deallocate(ind_cell2)
                                allocate(tempvar(0:amr%twondim,nvarh))
                                allocate(tempson(0:amr%twondim))
                                if (read_gravity) allocate(tempgrav_var(0:amr%twondim,1:4))
                                ! Just correct central cell vectors for the region
                                tempvar(0,:) = var(ind_nbor(1,0),:)
                                tempson(0)       = son(ind_nbor(1,0))
                                if (read_gravity) tempgrav_var(0,:) = grav_var(ind_nbor(1,0),:)
                                tempvar(0,varIDs%vx:varIDs%vz) = vtemp
                                if (read_gravity) tempgrav_var(0,2:4) = gtemp

                                do inbor=1,amr%twondim
                                    tempvar(inbor,:) = var(ind_nbor(1,inbor),:)
                                    tempson(inbor)       = son(ind_nbor(1,inbor))
                                    if (read_gravity) tempgrav_var(inbor,:) = grav_var(ind_nbor(1,inbor),:)
                                end do
                                if (read_gravity) then
                                    ok_filter = filter_cell(reg,filt,xtemp,dx,tempvar,tempson,&
                                                            &trans_matrix,tempgrav_var)
                                else
                                    ok_filter = filter_cell(reg,filt,xtemp,dx,tempvar,tempson,&
                                                            &trans_matrix)
                                end if
                                ok_cell= ok_cell.and..not.ref(i).and.ok_filter

                                ! If we are avoiding substructure, check whether we are safe
                                if (prof_data%nsubs>0) then
                                    ok_sub = .true.
                                    do isub=1,prof_data%nsubs
                                        ok_sub = ok_sub .and. filter_sub(prof_data%subs(isub),xorig(i,:))
                                    end do
                                    ok_cell = ok_cell .and. ok_sub
                                end if
                                if (ok_cell) then
                                    binpos = 0
                                    if (read_gravity) then
                                        call findbinpos(reg,xtemp,tempvar,tempson,&
                                                        & dx,binpos,ytemp,trans_matrix,&
                                                        & prof_data%scaletype,prof_data%nbins,&
                                                        & prof_data%xdata,prof_data%xvarname,&
                                                        & tempgrav_var)
                                        if (binpos.ne.0) call bindata(reg,x(i,:),tempvar,tempson,dx,prof_data,binpos,trans_matrix,tempgrav_var)
                                    else
                                        call findbinpos(reg,xtemp,tempvar,tempson,&
                                                        & dx,binpos,ytemp,trans_matrix,&
                                                        & prof_data%scaletype,prof_data%nbins,&
                                                        & prof_data%xdata,prof_data%xvarname)
                                        if (binpos.ne.0) call bindata(reg,x(i,:),tempvar,tempson,dx,prof_data,binpos,trans_matrix)
                                    end if
                                    if (binpos.ne.0)total_ncell = total_ncell + 1
                                endif
                                deallocate(tempvar,tempson)
                                if (read_gravity) deallocate(tempgrav_var)
                            end do ngridaloop
                        end do cellloop
                        deallocate(ref,x,ind_cell,xorig)
                    endif
                end do levelloop2
                deallocate(nbor,son,var,cellpos)

                if (read_gravity) then
                    deallocate(grav_var)
                end if
            end do cpuloop
            write(*,*)'Total number of cells used: ', total_ncell
        end subroutine get_cells_onedprofile_neigh

        subroutine get_cells_onedprofile_fast
            use vectors
            use coordinate_systems
            implicit none

            ! Specific variables for this subroutine
            integer :: i,j,k,binpos
            integer :: ipos,icpu,ilevel,ind,idim,ivar,ifilt,isub
            integer :: ix,iy,iz,ngrida,nx_full,ny_full,nz_full
            integer :: tot_pos,tot_ref,total_ncell,tot_insubs
            integer :: tot_sel
            integer :: nvarh
            integer :: roterr
            character(5) :: nchar,ncharcpu
            character(128) :: nomfich
            real(dbl) :: distance,dx,ytemp
            type(vector) :: xtemp,vtemp,gtemp
            logical :: ok_cell,ok_filter,ok_cell_each,ok_sub,read_gravity
            integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
            real(dbl),dimension(1:8,1:3) :: xc
            real(dbl),dimension(3,3) :: trans_matrix
            real(dbl),dimension(:,:),allocatable :: xg,x,xorig
            real(dbl),dimension(:,:,:),allocatable :: var,grav_var
            real(dbl),dimension(:,:),allocatable :: tempvar
            real(dbl),dimension(:,:),allocatable :: tempgrav_var
            integer,dimension(:,:),allocatable :: son
            integer,dimension(:),allocatable :: tempson
            logical,dimension(:),allocatable :: ref

            total_ncell = 0
            tot_pos = 0
            tot_ref = 0
            tot_insubs = 0
            tot_sel = 0
            ! Check whether we need to read the gravity files
            read_gravity = .false.
            do ivar=1,prof_data%nyvar
                if (prof_data%yvarnames(ivar)(1:4) .eq. 'grav' .or.&
                & trim(prof_data%yvarnames(ivar)) .eq. 'neighbour_accuracy') then
                    read_gravity = .true.
                    write(*,*)'Reading gravity files...'
                    exit
                endif
            end do

            ! Allocate grids
            allocate(ngridfile(1:amr%ncpu+amr%nboundary,1:amr%nlevelmax))
            allocate(ngridlevel(1:amr%ncpu,1:amr%nlevelmax))
            if(amr%nboundary>0)allocate(ngridbound(1:amr%nboundary,1:amr%nlevelmax))

            ! Compute linear transformation
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

                ! Make sure that we are not trying to access to far in the refinement map…
                call check_lmax(ngridfile)
                ! Open HYDRO file and skip header
                nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
                open(unit=11,file=nomfich,status='old',form='unformatted')
                read(11)
                read(11)nvarh
                read(11)
                read(11)
                read(11)
                read(11)

                if (read_gravity) then
                    ! Open GRAV file and skip header
                    nomfich=TRIM(repository)//'/grav_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
                    open(unit=12,file=nomfich,status='old',form='unformatted')
                    read(12) !ncpu
                    read(12) !ndim
                    read(12) !nlevelmax
                    read(12) !nboundary 
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
                        allocate(xg (1:ngrida,1:amr%ndim))
                        allocate(son(1:ngrida,1:amr%twotondim))
                        allocate(var(1:ngrida,1:amr%twotondim,1:nvarh))
                        allocate(x  (1:ngrida,1:amr%ndim))
                        allocate(xorig(1:ngrida,1:amr%ndim))
                        allocate(ref(1:ngrida))
                        if (read_gravity)allocate(grav_var(1:ngrida,1:amr%twotondim,1:4))
                    endif

                    ! Loop over domains
                    domloop: do j=1,amr%nboundary+amr%ncpu
                        
                        ! Read AMR data
                        if (ngridfile(j,ilevel)>0) then
                            read(10) ! Skip grid index
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
                            do ind=1,2*amr%ndim
                                read(10) ! Skip nbor index
                            end do

                            ! Read son index
                            do ind=1,amr%twotondim
                                if(j.eq.icpu)then
                                    read(10)son(:,ind)
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
                                varloop: do ivar=1,nvarh
                                    if (j.eq.icpu) then
                                        read(11)var(:,ind,ivar)
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
                                    if (j.eq.icpu) then
                                        read(12)grav_var(:,ind,1)
                                    else
                                        read(12)
                                    end if
                                    do ivar=1,amr%ndim
                                        if (j.eq.icpu) then
                                            read(12)grav_var(:,ind,ivar+1)
                                        else
                                            read(12)
                                        end if
                                    end do
                                end do
                            end if
                        end if
                    end do domloop

                    ! Finally, get to every cell
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
                            do i=1,ngrida
                                ref(i) = son(i,ind)>0.and.ilevel<amr%lmax
                            end do
                            xorig  = x
                            ngridaloop: do i=1,ngrida
                                ! Check if cell is inside the desired region
                                distance = 0D0
                                xtemp = x(i,:)
                                xtemp = xtemp - reg%centre
                                call rotate_vector(xtemp,trans_matrix)
                                x(i,:) = xtemp
                                call checkifinside(x(i,:),reg,ok_cell,distance)
                                if(ok_cell) tot_pos = tot_pos + 1
                                if(.not.ref(i)) tot_ref = tot_ref + 1

                                ! If we are avoiding substructure, check whether we are safe
                                if (prof_data%nsubs>0) then
                                    ok_sub = .true.
                                    do isub=1,prof_data%nsubs
                                        ok_sub = ok_sub .and. filter_sub(prof_data%subs(isub),xorig(i,:))
                                    end do
                                    if (.not.ok_sub) tot_insubs = tot_insubs + 1
                                    ok_cell = ok_cell .and. ok_sub
                                end if
                                ok_cell = ok_cell.and.(.not.ref(i))
                                if (ok_cell) then
                                    ! Transform position to galaxy frame
                                    xtemp = xorig(i,:)
                                    xtemp = xtemp - reg%centre
                                    call rotate_vector(xtemp,trans_matrix)
                                    ! Velocity transformed
                                    vtemp = var(i,ind,varIDs%vx:varIDs%vz)
                                    vtemp = vtemp - reg%bulk_velocity
                                    call rotate_vector(vtemp,trans_matrix)

                                    ! Gravitational acc
                                    if (read_gravity) then
                                        gtemp = grav_var(i,ind,2:4)
                                        call rotate_vector(gtemp,trans_matrix)
                                    endif
                                    allocate(tempvar(0:amr%twondim,nvarh))
                                    allocate(tempson(0:amr%twondim))
                                    if (read_gravity) allocate(tempgrav_var(0:amr%twondim,1:4))
                                    ! Just add central cell as we do not want neighbours
                                    tempvar(0,:) = var(i,ind,:)
                                    tempson(0)       = son(i,ind)
                                    if (read_gravity) tempgrav_var(0,:) = grav_var(i,ind,:)
                                    tempvar(0,varIDs%vx:varIDs%vz) = vtemp
                                    if (read_gravity) tempgrav_var(0,2:4) = gtemp
                                    if (read_gravity) then
                                        ok_filter = filter_cell(reg,filt,xtemp,dx,tempvar,tempson,&
                                                                &trans_matrix,tempgrav_var)
                                    else
                                        ok_filter = filter_cell(reg,filt,xtemp,dx,tempvar,tempson,&
                                                                &trans_matrix)
                                    end if
                                    ok_cell= ok_cell.and.ok_filter
                                    tot_sel = tot_sel + 1
                                    if (ok_cell) then
                                        binpos = 0
                                        if (read_gravity) then
                                            call findbinpos(reg,xtemp,tempvar,tempson,&
                                                            & dx,binpos,ytemp,trans_matrix,&
                                                            & prof_data%scaletype,prof_data%nbins,&
                                                            & prof_data%xdata,prof_data%xvarname,&
                                                            & tempgrav_var)
                                            if (binpos.ne.0) call bindata(reg,x(i,:),tempvar,tempson,dx,prof_data,binpos,trans_matrix,tempgrav_var)
                                        else
                                            call findbinpos(reg,xtemp,tempvar,tempson,&
                                                            & dx,binpos,ytemp,trans_matrix,&
                                                            & prof_data%scaletype,prof_data%nbins,&
                                                            & prof_data%xdata,prof_data%xvarname)
                                            if (binpos.ne.0) call bindata(reg,x(i,:),tempvar,tempson,dx,prof_data,binpos,trans_matrix)
                                        end if
                                        total_ncell = total_ncell + 1
                                    endif
                                    deallocate(tempvar,tempson)
                                    if (read_gravity) deallocate(tempgrav_var)
                                end if
                            end do ngridaloop
                        end do cellloop
                        deallocate(xg,son,var,ref,x,xorig)
                        if (read_gravity) then
                            deallocate(grav_var)
                        end if
                    endif
                end do levelloop
                close(10)
                close(11)
            end do cpuloop
        write(*,*)'Total number of cells used: ', total_ncell
        write(*,*)'Total number of cells in region and refined: ', tot_sel
        write(*,*)'Total number of cells refined: ', tot_ref
        write(*,*)'Total number of cells in region: ', tot_pos
        write(*,*)'Total number of cells in substructures: ', tot_insubs
        end subroutine get_cells_onedprofile_fast
    end subroutine onedprofile

    subroutine twodprofile(repository,reg,filt,prof_data,lmax,scaletype)
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        type(filter),intent(in) :: filt
        type(profile_handler_twod),intent(inout) :: prof_data
        integer,intent(in) :: lmax
        character(128),intent(in) :: scaletype

        call read_hydrofile_descriptor(repository)
        call init_amr_read(repository)
        amr%lmax = lmax
        if (lmax.eq.0) amr%lmax = amr%nlevelmax
        
        prof_data%xdata = 0D0
        prof_data%ydata = 0D0
        prof_data%zdata = 0D0

        prof_data%xdata = makebins(reg,prof_data%xvarname,prof_data%nbins(1),scaletype)
        prof_data%ydata = makebins(reg,prof_data%yvarname,prof_data%nbins(2),scaletype)
        write(*,*)'lmax: ',amr%lmax
        call get_cpu_map(reg)
        write(*,*)'ncpu_read:',amr%ncpu_read
        call get_cells_twodprofile(repository,reg,filt,prof_data,scaletype)
        call renormalise_bins_twod(prof_data)

    end subroutine twodprofile

    subroutine get_cells_twodprofile(repository,reg,filt,prof_data,scaletype)
        use vectors
        use coordinate_systems
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(region), intent(in)  :: reg
        type(filter),intent(in) :: filt
        type(profile_handler_twod),intent(inout) :: prof_data
        character(128),intent(in) :: scaletype
        integer :: xbinpos,ybinpos
        logical :: ok_cell,ok_filter,read_gravity
        integer :: i,j,k
        integer :: ipos,icpu,ilevel,ind,idim,ivar,iskip,inbor,ison
        integer :: ix,iy,iz,ngrida,nx_full,ny_full,nz_full,total_ncell
        integer :: nvarh
        integer :: roterr
        character(5) :: nchar,ncharcpu
        character(128) :: nomfich
        real(dbl) :: distance,dx
        type(vector) :: xtemp,vtemp,gtemp
        integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
        real(dbl),dimension(:),allocatable :: xxg,son_dens
        real(dbl),dimension(1:8,1:3) :: xc
        real(dbl),dimension(1:3,1:3) :: trans_matrix
        real(dbl),dimension(:,:),allocatable :: x
        real(dbl),dimension(:,:),allocatable :: var
        real(dbl),dimension(:,:),allocatable :: grav_var
        real(dbl),dimension(:,:),allocatable :: tempvar
        real(dbl),dimension(:,:),allocatable :: tempgrav_var
        real(dbl),dimension(:,:),allocatable :: cellpos
        integer,dimension(:,:),allocatable :: nbor
        integer,dimension(:),allocatable :: son,tempson,iig
        integer,dimension(:),allocatable :: ind_cell,ind_cell2
        integer ,dimension(1,0:amr%twondim) :: ind_nbor
        logical,dimension(:),allocatable :: ref
        type(level),dimension(1:100) :: grid

        total_ncell = 0
        
        ! Check whether we need to read the gravity files
        read_gravity = .false.
        do ivar=1,prof_data%nzvar
            if (prof_data%zvarnames(ivar)(1:4) .eq. 'grav' .or.&
            & trim(prof_data%zvarnames(ivar)) .eq. 'neighbour_accuracy') then
                read_gravity = .true.
                write(*,*)'Reading gravity files...'
                exit
            endif
        end do

        allocate(ngridfile(1:amr%ncpu+amr%nboundary,1:amr%nlevelmax))
        allocate(ngridlevel(1:amr%ncpu,1:amr%nlevelmax))
        if(amr%nboundary>0)allocate(ngridbound(1:amr%nboundary,1:amr%nlevelmax))
        ! Compute hierarchy
        do ilevel=1,amr%lmax
            grid(ilevel)%ngrid = 0
        end do
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

            allocate(nbor(1:amr%ngridmax,1:amr%twondim))
            allocate(son(1:amr%ncoarse+amr%twotondim*amr%ngridmax))
            nbor = 0
            son = 0
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
            read(10)son(1:amr%ncoarse)
            read(10)
            read(10)

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
            allocate(cellpos(1:amr%ncoarse+amr%twotondim*amr%ngridmax,1:3))
            cellpos = 0d0
            var = 0d0
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
            levelloop1: do ilevel=1,amr%lmax
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
                grid(ilevel)%ngrid = 0
                ! Loop over domains
                domloop: do j=1,amr%nboundary+amr%ncpu
                    ! Allocate work arrays
                    ngrida = ngridfile(j,ilevel)
                    if(ngrida>0)then
                        if (allocated(grid(ilevel)%ind_grid)) deallocate(grid(ilevel)%ind_grid)
                        if (allocated(grid(ilevel)%xg)) deallocate(grid(ilevel)%xg)
                        allocate(grid(ilevel)%ind_grid(1:ngrida))
                        allocate(grid(ilevel)%xg (1:ngrida,1:amr%ndim))
                        allocate(iig(1:ngrida))
                        allocate(xxg(1:ngrida))
                        
                        ! Read AMR data
                        read(10) grid(ilevel)%ind_grid
                        read(10) ! Skip next index
                        read(10) ! Skip prev index
                        if(j.eq.icpu) then
                            if (allocated(grid(ilevel)%real_ind)) deallocate(grid(ilevel)%real_ind)
                            allocate(grid(ilevel)%real_ind(1:ngrida))
                            grid(ilevel)%real_ind = grid(ilevel)%ind_grid
                            grid(ilevel)%ngrid = ngridfile(j,ilevel)
                        end if
                        ! Read grid center
                        do idim=1,amr%ndim
                            read(10)xxg
                            grid(ilevel)%xg(:,idim) = xxg(:)
                        end do
                        
                        read(10) ! Skip father index
                        ! Read nbor index
                        do ind=1,amr%twondim
                            read(10)iig
                            nbor(grid(ilevel)%ind_grid(:),ind) = iig(:)
                        end do
                        ! Read son index
                        do ind=1,amr%twotondim
                            iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                            read(10)iig
                            son(grid(ilevel)%ind_grid(:)+iskip) = iig(:)
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
                    if(ngrida>0)then
                        ! Read hydro variables
                        tndimloop: do ind=1,amr%twotondim
                            iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                            varloop: do ivar=1,nvarh
                                read(11)xxg
                                var(grid(ilevel)%ind_grid(:)+iskip,ivar) = xxg(:)
                            end do varloop
                        end do tndimloop
                    endif

                    if (read_gravity) then
                        ! Read GRAV data
                        read(12)
                        read(12)
                        if(ngrida>0)then
                            do ind=1,amr%twotondim
                                iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                                read(12)xxg
                                grav_var(grid(ilevel)%ind_grid(:)+iskip,1) = xxg(:)
                                do ivar=1,amr%ndim
                                    read(12)xxg
                                    grav_var(grid(ilevel)%ind_grid(:)+iskip,ivar+1) = xxg(:)
                                end do
                            end do
                        end if
                    end if
                    !Compute positions
                    if(ngrida>0)then
                        do ind=1,amr%twotondim
                            iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                            do i=1,ngrida
                                do ivar=1,amr%ndim
                                    cellpos(grid(ilevel)%ind_grid(i)+iskip,ivar)=(grid(ilevel)%xg(i,ivar)+xc(ind,ivar)-amr%xbound(ivar))
                                end do
                            end do
                        end do
                    end if
                    if (ngrida>0) deallocate(iig,xxg)
                end do domloop
            end do levelloop1
            close(10)
            close(11)
            if (read_gravity) then
                close(12)
            end if
            ! Loop over levels again now with arrays fully filled
            levelloop2: do ilevel=1,amr%lmax
                ! Geometry
                dx = 0.5**ilevel
                nx_full = 2**ilevel
                ny_full = 2**ilevel

                ! Allocate work arrays
                ngrida = grid(ilevel)%ngrid
                if(ngrida>0)then
                    allocate(ind_cell(1:ngrida))
                    allocate(x  (1:ngrida,1:amr%ndim))
                    allocate(ref(1:ngrida))
                endif
                !Compute map
                if (ngrida>0) then
                    ! Loop over cells
                    cellloop: do ind=1,amr%twotondim
                        ! Get cell indexes
                        iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                        do i=1,ngrida
                            ind_cell(i) = iskip+grid(ilevel)%real_ind(i)
                        end do
                        ! Compute cell center
                        do i=1,ngrida
                            x(i,:)=cellpos(ind_cell(i),:)
                        end do
                        ! Check if cell is refined
                        do i=1,ngrida
                            ref(i) = son(ind_cell(i))>0.and.ilevel<amr%lmax
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

                            ! Gravitational acc --> ONLY FOR CENTRAL CELL
                            if (read_gravity) then
                                gtemp = grav_var(ind_cell(i),2:4)
                                call rotate_vector(gtemp,trans_matrix)
                            endif
                            ! Get neighbours
                            allocate(ind_cell2(1))
                            ind_cell2(1) = ind_cell(i)
                            call getnbor(son,nbor,ind_cell2,ind_nbor,1)
                            deallocate(ind_cell2)
                            allocate(tempvar(0:amr%twondim,nvarh))
                            allocate(tempson(0:amr%twondim))
                            if (read_gravity) allocate(tempgrav_var(0:amr%twondim,1:4))
                            ! Just correct central cell vectors for the region
                            tempvar(0,:) = var(ind_nbor(1,0),:)
                            tempson(0)       = son(ind_nbor(1,0))
                            if (read_gravity) tempgrav_var(0,:) = grav_var(ind_nbor(1,0),:)
                            tempvar(0,varIDs%vx:varIDs%vz) = vtemp
                            if (read_gravity) tempgrav_var(0,2:4) = gtemp
                            do inbor=1,amr%twondim
                                tempvar(inbor,:) = var(ind_nbor(1,inbor),:)
                                tempson(inbor)       = son(ind_nbor(1,inbor))
                                if (read_gravity) tempgrav_var(inbor,:) = grav_var(ind_nbor(1,inbor),:)
                            end do
                            if (read_gravity) then
                                ok_filter = filter_cell(reg,filt,xtemp,dx,tempvar,tempson,&
                                                        &trans_matrix,tempgrav_var)
                            else
                                ok_filter = filter_cell(reg,filt,xtemp,dx,tempvar,tempson,&
                                                        &trans_matrix)
                            end if
                            ok_cell= ok_cell.and..not.ref(i).and.ok_filter
                            if (ok_cell) then
                                xbinpos = 0; ybinpos=0
                                total_ncell = total_ncell + 1
                                if (read_gravity) then
                                    call findbinpos_twod(reg,distance,x(i,:),tempvar,tempson,&
                                                        &dx,prof_data,scaletype,xbinpos,ybinpos,&
                                                        &trans_matrix,tempgrav_var)
                                    if (xbinpos.ne.0.and.ybinpos.ne.0) call bindata_twod(reg,x(i,:),&
                                                                            &tempvar,tempson,dx,prof_data,&
                                                                            &xbinpos,ybinpos,&
                                                                            &trans_matrix,tempgrav_var)
                                else
                                    call findbinpos_twod(reg,distance,x(i,:),tempvar,tempson,&
                                                    &dx,prof_data,scaletype,xbinpos,ybinpos,trans_matrix)
                                    if (xbinpos.ne.0.and.ybinpos.ne.0) call bindata_twod(reg,x(i,:),&
                                                                            &tempvar,tempson,dx,prof_data,&
                                                                            &xbinpos,ybinpos,trans_matrix)
                                end if
                            endif
                            deallocate(tempvar,tempson)
                            if (read_gravity) deallocate(tempgrav_var)
                        end do ngridaloop
                    end do cellloop
                    deallocate(ref,x,ind_cell)
                endif
            end do levelloop2
            deallocate(nbor,son,var,cellpos)
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