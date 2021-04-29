module part_profiles
    use local
    use io_ramses
    use filtering
    use cosmology

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

    subroutine findbinpos(sim,reg,distance,part,prof,ibin)
        use vectors
        use geometrical_regions
        implicit none
        type(sim_info),intent(in) :: sim
        type(region),intent(in) :: reg
        real(dbl),intent(in) :: distance
        type(particle),intent(in) :: part
        type(profile_handler),intent(in) :: prof
        integer,intent(inout) :: ibin
        real(dbl) :: value

        if (prof%xvarname.eq.reg%criteria_name) then
            value = distance
        else
            call getpartvalue(sim,reg,part,prof%xvarname,value)
        endif
        ibin = int(dble(prof%nbins)*value/prof%xdata(prof%nbins)) + 1
    end subroutine findbinpos


    subroutine bindata(sim,reg,part,prof,ibin)
        use vectors
        use geometrical_regions
        implicit none
        type(sim_info),intent(in) :: sim
        type(region),intent(in) :: reg
        type(particle),intent(in) :: part
        type(profile_handler),intent(inout) :: prof
        integer,intent(in) :: ibin
        integer :: i,j
        real(dbl) :: ytemp,wtemp,bigwtemp,bigatemp
        yvarloop: do i=1,prof%nyvar
            call getpartvalue(sim,reg,part,prof%yvarnames(i),ytemp)
            wvarloop: do j=1,prof%nwvar
                if (prof%wvarnames(j)=='counts'.or.prof%wvarnames(j)=='cumulative') then
                    wtemp = 1D0
                else
                    call getpartvalue(sim,reg,part,prof%wvarnames(j),wtemp)
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

    subroutine get_parts_onedprofile(repository,amr,sim,reg,filt,prof_data)
        use vectors
        use coordinate_systems
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(amr_info),intent(inout) :: amr
        type(sim_info),intent(inout) :: sim
        type(region), intent(in)  :: reg
        type(filter),intent(in) :: filt
        type(profile_handler),intent(inout) :: prof_data

        logical :: ok_part,ok_filter
        integer :: roterr
        integer :: i,j,k
        integer :: ipos,icpu,binpos
        integer :: npart,npart2,nstar,ncpu2,ndim2
        real(dbl) :: distance
        real(dbl),dimension(1:3,1:3) :: trans_matrix
        character(5) :: nchar,ncharcpu
        character(128) :: nomfich
        type(vector) :: xtemp,vtemp
        type(particle) :: part
        integer,dimension(:),allocatable :: id
        real(dbl),dimension(:),allocatable :: m,age,met,imass
        real(dbl),dimension(:,:),allocatable :: x,v

        ! Compute rotation matrix
        trans_matrix = 0D0
        call new_z_coordinates(reg%axis,trans_matrix,roterr)
        if (roterr.eq.1) then
            write(*,*) 'Incorrect CS transformation!'
            stop
        endif

        ! Cosmological model
        if (sim%aexp.eq.1.and.sim%h0.eq.1)sim%cosmo=.false.
        if (sim%cosmo) then
            call cosmology_model(sim)
        else
            sim%time_simu = sim%t
            write(*,*)'Age simu=',sim%time_simu*sim%unit_t/(365.*24.*3600.*1d9)
        endif

        ! Check number of particles in selected CPUs
        ipos = INDEX(repository,'output_')
        nchar = repository(ipos+7:ipos+13)
        npart = 0
        do k=1,amr%ncpu_read
            icpu = amr%cpu_list(k)
            call title(icpu,ncharcpu)
            nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
            write(*,*)'Processing file '//TRIM(nomfich)
            open(unit=1,file=nomfich,status='old',form='unformatted')
            read(1)ncpu2
            read(1)ndim2
            read(1)npart2
            read(1)
            read(1)nstar
            close(1)
            npart=npart+npart2
        end do
        write(*,*)'Found ',npart,' particles.'
        if(nstar>0)then
            write(*,*)'Found ',nstar,' star particles.'
        endif

        ! Compute binned variables
        cpuloop: do k=1,amr%ncpu_read
            icpu = amr%cpu_list(k)
            call title(icpu,ncharcpu)
            nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
            open(unit=1,file=nomfich,status='old',form='unformatted')
            read(1)ncpu2
            read(1)ndim2
            read(1)npart2
            read(1)
            read(1)
            read(1)
            read(1)
            read(1)
            allocate(m(1:npart2))
            if(nstar>0)then
                allocate(age(1:npart2))
                allocate(id(1:npart2))
                allocate(met(1:npart2))
                allocate(imass(1:npart2))
            endif
            allocate(x(1:npart2,1:ndim2))
            allocate(v(1:npart2,1:ndim2))

            ! Read position
            do i=1,amr%ndim
                read(1)m
                x(1:npart2,i) = m/sim%boxlen
            end do

            ! Read velocity
            do i=1,amr%ndim
                read(1)m
                v(1:npart2,i) = m
            end do

            ! Read mass
            read(1)m
            if (nstar>0) then
                read(1)id
                read(1) ! Skip level
                read(1)age
                read(1)met
                read(1)imass
            endif
            close(1)

            ! Project variables into map for particles
            ! of interest in the region
            partloop: do i=1,npart2
                distance = 0D0
                part%x = x(i,:)
                part%v = v(i,:)
                part%m = m(i)
                if (nstar>0) then
                    part%id = id(i)
                    part%age = age(i)
                    part%met = met(i)
                    part%imass = imass(i)
                else
                    part%id = 0
                    part%age = 0D0
                    part%met = 0D0
                    part%imass = 0D0
                endif
                ! Check if particle is inside the desired region
                call rotate_vector(part%x,trans_matrix)
                x(i,:) = part%x
                call checkifinside(x(i,:),reg,ok_part,distance)
                ok_filter = filter_particle(sim,reg,filt,part)
                ok_part = ok_part.and.ok_filter
                if (ok_part) then
                    call rotate_vector(part%v,trans_matrix)
                    binpos = 0
                    call findbinpos(sim,reg,distance,part,prof_data,binpos)
                    if (binpos.ne.0)  call bindata(sim,reg,part,prof_data,binpos)
                endif
            end do partloop
            deallocate(m,x,v)
            if (nstar>0)deallocate(id,age,met,imass)
        end do cpuloop
    end subroutine get_parts_onedprofile

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
        call get_parts_onedprofile(repository,amr,sim,reg,filt,prof_data)

        call renormalise_bins(prof_data)
    end subroutine onedprofile
end module part_profiles