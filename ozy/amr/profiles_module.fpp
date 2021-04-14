# 1 "profiles_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "profiles_module.f90"
module amr_profiles
    use local
    use io_ramses
    use filtering

    type profile_handler
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

    contains

    subroutine allocate_profile_handler(prof)
        implicit none
        type(profile_handler),intent(inout) :: prof

        if (.not.allocated(prof%yvarnames)) allocate(prof%yvarnames(prof%nyvar))
        if (.not.allocated(prof%wvarnames)) allocate(prof%wvarnames(prof%nwvar))
        if (.not.allocated(prof%xdata)) allocate(prof%xdata(prof%nbins))
        if (.not.allocated(prof%ydata)) allocate(prof%ydata(prof%nbins,prof%nyvar,prof%nwvar,4))
    end subroutine allocate_profile_handler

    subroutine makebins(reg,varname,nbins,bins)
        use geometrical_regions

        implicit none
        type(region),intent(in) :: reg
        character(128),intent(in) :: varname
        integer,intent(in) :: nbins
        real(dbl),dimension(1:nbins) :: bins
        integer :: n

        select case (TRIM(varname))
        case('r_sphere','r_cyl')
            do n=1,nbins
                bins(n) = dble(n)*(reg%rmax-reg%rmin)/dble(nbins)
            end do
        !TODO: Add more cases
        end select
    end subroutine makebins

    subroutine findbinpos(reg,varIDs,distance,pos,cellvars,cellsize,prof,ibin)
        use vectors
        use geometrical_regions
        implicit none
        type(region),intent(in) :: reg
        real(dbl),dimension(1:3),intent(in) :: pos
        real(dbl),intent(in) :: distance
        type(hydroID),intent(in) :: varIDs
        real(dbl),dimension(1:varIDs%nvar),intent(in) :: cellvars
        real(dbl),intent(in) :: cellsize
        type(profile_handler),intent(in) :: prof
        integer,intent(inout) :: ibin
        real(dbl) :: value
        type(vector) :: x

        x = pos
        if (prof%xvarname.eq.reg%criteria_name) then
            value = distance
        else
            call getvarvalue(varIDs,reg,cellsize,x,cellvars,prof%xvarname,value)
        endif
        ibin = int(dble(prof%nbins)*value/prof%xdata(prof%nbins)) + 1
    end subroutine findbinpos

    subroutine bindata(reg,varIDs,pos,cellvars,cellsize,prof,ibin)
        use vectors
        use geometrical_regions
        implicit none
        type(region),intent(in) :: reg
        real(dbl),dimension(1:3),intent(in) :: pos
        type(hydroID),intent(in) :: varIDs
        real(dbl),dimension(1:varIDs%nvar),intent(in) :: cellvars
        real(dbl),intent(in) :: cellsize
        type(profile_handler),intent(inout) :: prof
        integer,intent(in) :: ibin
        integer :: i,j
        real(dbl) :: ytemp,wtemp,bigwtemp,bigatemp
        type(vector) :: x
        x = pos
        yvarloop: do i=1,prof%nyvar
            call getvarvalue(varIDs,reg,cellsize,x,cellvars,prof%yvarnames(i),ytemp)
            wvarloop: do j=1,prof%nwvar
                if (prof%wvarnames(j)=='counts'.or.prof%wvarnames(j)=='cumulative') then
                    wtemp = 1D0
                else
                    call getvarvalue(varIDs,reg,cellsize,x,cellvars,prof%wvarnames(j),wtemp)
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

                ! OLD METHOD (deprecated)
                ! prof%ydata(ibin,i,j,1) = prof%ydata(ibin,i,j,1) + ytemp*wtemp
                ! ! W_k-1     
                ! bigwtemp = prof%ydata(ibin,i,j,2)
                ! ! Sum weight
                ! prof%ydata(ibin,i,j,2) = prof%ydata(ibin,i,j,2) + wtemp
                ! ! A_k-1
                ! bigatemp = prof%ydata(ibin,i,j,3)
                ! ! A_k = A_k-1 + wtemp/W_k * (ytemp - A_k-1)
                ! prof%ydata(ibin,i,j,3) = prof%ydata(ibin,i,j,3) + (wtemp / prof%ydata(ibin,i,j,2))&
                !                         &* (ytemp - prof%ydata(ibin,i,j,3))
                ! ! Q_k = Q_k-1 + wtemp
                ! prof%ydata(ibin,i,j,4) = prof%ydata(ibin,i,j,4) + (wtemp * bigwtemp &
                !                         &/ prof%ydata(ibin,i,j,2)) * (ytemp - bigatemp)**2
            end do wvarloop
        end do yvarloop
    end subroutine bindata

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
                        ! OLD METHOD (deprecated)
                        ! prof_data%ydata(ibin,iy,iw,2) =  prof_data%ydata(ibin,iy,iw,4) &
                        !                                 &/ (prof_data%ydata(ibin,iy,iw,2) - 1D0)
                    endif
                end do wloop
            end do yloop
        end do binloop
    end subroutine renormalise_bins

    subroutine get_cells_onedprofile(repository,amr,reg,filt,varIDs,prof_data)
        use vectors
        use coordinate_systems
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(amr_info),intent(inout) :: amr
        type(region), intent(in)  :: reg
        type(filter),intent(in) :: filt
        type(hydroID),intent(in) :: varIDs
        type(profile_handler),intent(inout) :: prof_data
        integer :: binpos
        logical :: ok_cell,ok_filter
        ! logical :: filter_cell
        integer :: i,j,k
        integer :: ipos,icpu,ilevel,ind,idim,ivar
        integer :: ix,iy,iz,ngrida,nx_full,ny_full,nz_full
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
        real(dbl),dimension(:,:,:),allocatable :: var
        integer,dimension(:,:),allocatable :: son
        logical,dimension(:),allocatable :: ref

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

            ! Make sure that we are not trying to access to far in the refinement map…
            call check_lmax(ngridfile,amr)
            ! Open HYDRO file and skip header
            nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
            open(unit=11,file=nomfich,status='old',form='unformatted')
            read(11)
            read(11)nvarh
            read(11)
            read(11)
            read(11)
            read(11)

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
                    allocate(ref(1:ngrida))
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
                        do i=1,ngrida
                            ref(i) = son(i,ind)>0.and.ilevel<amr%lmax
                        end do
                        ngridaloop: do i=1,ngrida
                            ! Check if cell is inside the desired region
                            distance = 0D0
                            xtemp = x(i,:)
                            xtemp = xtemp - reg%centre
                            call rotate_vector(xtemp,trans_matrix)
                            x(i,:) = xtemp
                            call checkifinside(x(i,:),reg,ok_cell,distance)
                            ok_filter = filter_cell(varIDs,reg,filt,xtemp,dx,var(i,ind,:))
                            ok_cell= ok_cell.and..not.ref(i).and.ok_filter
                            if (ok_cell) then
                                vtemp = var(i,ind,varIDs%vx:varIDs%vz)
                                call rotate_vector(vtemp,trans_matrix)
                                var(i,ind,varIDs%vx:varIDs%vz) = vtemp
                                binpos = 0
                                call findbinpos(reg,varIDs,distance,x(i,:),var(i,ind,:),dx,prof_data,binpos)
                                if (binpos.ne.0) call bindata(reg,varIDs,x(i,:),var(i,ind,:),dx,prof_data,binpos)
                            endif
                        end do ngridaloop
                    end do cellloop

                    deallocate(xg,son,var,ref,x)
                endif
            end do levelloop
            close(10)
            close(11)
        end do cpuloop
    end subroutine get_cells_onedprofile

    subroutine onedprofile(repository,reg,filt,prof_data,lmax)
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        type(filter),intent(in) :: filt
        type(profile_handler),intent(inout) :: prof_data
        integer,intent(in) :: lmax

        type(hydroID) :: varIDs
        type(amr_info) :: amr
        type(sim_info) :: sim

        call read_hydrofile_descriptor(repository,varIDs)

        call init_amr_read(repository,amr,sim)
        amr%lmax = lmax
        if (lmax.eq.0) amr%lmax = amr%nlevelmax
        prof_data%xdata = 0D0
        prof_data%ydata = 0D0
        call makebins(reg,prof_data%xvarname,prof_data%nbins,prof_data%xdata)
        
        call get_cpu_map(reg,amr)
        write(*,*)'ncpu_read:',amr%ncpu_read
        call get_cells_onedprofile(repository,amr,reg,filt,varIDs,prof_data)

        call renormalise_bins(prof_data)

    end subroutine onedprofile
end module amr_profiles
