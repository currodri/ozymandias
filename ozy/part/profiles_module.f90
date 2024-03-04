module part_profiles
    use local
    use io_ramses
    use filtering
    use cosmology

    type profile_handler
        logical :: logscale
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
        if (.not.allocated(prof%xdata)) allocate(prof%xdata(0:prof%nbins))
        if (.not.allocated(prof%ydata)) allocate(prof%ydata(prof%nbins,prof%nyvar,prof%nwvar,4))
    end subroutine allocate_profile_handler

    subroutine makebins(reg,varname,nbins,bins,logscale)
        use geometrical_regions

        implicit none
        type(region),intent(in) :: reg
        character(128),intent(in) :: varname
        integer,intent(in) :: nbins
        real(dbl),dimension(0:nbins) :: bins
        logical,intent(in)::logscale
        integer :: n
        character(128) :: tempvar,vartype,var
        integer :: index
        real(dbl) :: rmin

        tempvar = TRIM(varname)
        index = scan(tempvar,'/')
        vartype = tempvar(1:index-1)
        var = tempvar(index+1:)
        select case (TRIM(var))
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
        !TODO: Add more cases
        end select
    end subroutine makebins

    subroutine findbinpos(reg,distance,part,prof,ibin)
        use vectors
        use geometrical_regions
        implicit none
        type(region),intent(in) :: reg
        real(dbl),intent(in) :: distance
        type(particle),intent(in) :: part
        type(profile_handler),intent(in) :: prof
        integer,intent(inout) :: ibin
        real(dbl) :: value

        if (prof%xvarname.eq.reg%criteria_name) then
            value = distance
        else
            call getpartvalue(reg,part,prof%xvarname,value)
        endif
        if (prof%logscale) value = log10(value)

        ibin = int(dble(prof%nbins)*(value-prof%xdata(0))/(prof%xdata(prof%nbins)-prof%xdata(0))) + 1
        if (value .eq. prof%xdata(prof%nbins)) then
            ibin = prof%nbins
        else if (value<prof%xdata(0).or.value>prof%xdata(prof%nbins)) then
            ibin = 0
        end if
    end subroutine findbinpos


    subroutine bindata(reg,part,prof,ibin)
        use vectors
        use geometrical_regions
        implicit none
        type(region),intent(in) :: reg
        type(particle),intent(in) :: part
        type(profile_handler),intent(inout) :: prof
        integer,intent(in) :: ibin
        integer :: i,j,index
        real(dbl) :: ytemp,wtemp,bigwtemp,bigatemp
        character(128) :: tempvar,vartype,varname
        yvarloop: do i=1,prof%nyvar
            call getpartvalue(reg,part,prof%yvarnames(i),ytemp)
            if (ytemp/=0D0) then
                wvarloop: do j=1,prof%nwvar
                    tempvar = TRIM(prof%wvarnames(j))
                    index = scan(tempvar,'/')
                    vartype = tempvar(1:index-1)
                    varname = tempvar(index+1:)
                    wtemp = 0D0
                    if (varname=='counts'.or.varname=='cumulative') then
                        wtemp = 1D0
                    else
                        call getpartvalue(reg,part,prof%wvarnames(j),wtemp)
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
            endif
        end do yvarloop
    end subroutine bindata

    subroutine renormalise_bins(prof_data)
        implicit none
        type(profile_handler),intent(inout) :: prof_data
        integer :: ibin,iy,iw,index
        real(dbl) :: Q_k,W_k,V_k,A_k
        character(128) :: tempvar,vartype,varname
        binloop: do ibin=1,prof_data%nbins
            yloop: do iy=1,prof_data%nyvar
                wloop: do iw=1,prof_data%nwvar
                    Q_k = prof_data%ydata(ibin,iy,iw,1)
                    W_k = prof_data%ydata(ibin,iy,iw,2)
                    V_k = prof_data%ydata(ibin,iy,iw,3)
                    A_k = prof_data%ydata(ibin,iy,iw,4)
                    tempvar = TRIM(prof_data%wvarnames(iw))
                    index = scan(tempvar,'/')
                    vartype = tempvar(1:index-1)
                    varname = tempvar(index+1:)
                    if (varname /= 'cumulative') then
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

    subroutine get_parts_onedprofile(repository,reg,filt,prof_data,tag_file,inverse_tag)
#ifndef LONGINT
        use utils, only:quick_sort_irg,binarysearch_irg
#else
        use utils, only:quick_sort_ilg,binarysearch_ilg
#endif
        use vectors
        use coordinate_systems
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(region), intent(in)  :: reg
        type(filter),intent(in) :: filt
        type(profile_handler),intent(inout) :: prof_data
        character(128),intent(in),optional :: tag_file
        logical,intent(in),optional :: inverse_tag

        logical :: ok_part,ok_filter,ok_tag
        integer :: roterr
        integer :: i,j,k,itag
        integer :: ipos,icpu,binpos
        integer :: npart,npart2,nstar,ntag
        integer :: ncpu2,ndim2
        real(dbl) :: distance
        real(dbl),dimension(1:3,1:3) :: trans_matrix
        character(5) :: nchar,ncharcpu
        character(128) :: nomfich
        type(vector) :: xtemp,vtemp
        type(particle) :: part
        character(6) :: ptype
        integer,dimension(:),allocatable :: order
#ifndef LONGINT
        integer(irg),dimension(:),allocatable :: id,tag_id
#else
        integer(ilg),dimension(:),allocatable :: id,tag_id
#endif
#ifndef IMASS
        integer(1),dimension(:), allocatable :: part_tags
#endif
        real(dbl),dimension(:),allocatable :: m,age,met,imass
        real(dbl),dimension(:,:),allocatable :: x,v

#ifndef IMASS
        if (sim%eta_sn .eq. -1D0) then
            write(*,*)': eta_sn=-1 and not IMASS --> should set this up!'
            stop
        end if
#endif
        ! If tagged particles file exists, read and allocate array
        if (present(tag_file)) then
            open(unit=58,file=TRIM(tag_file),status='old',form='formatted')
            write(*,*)'Reading particle tags file '//TRIM(tag_file)
            read(58,'(I11)')ntag
            write(*,*)'Number of tagged particles in file: ',ntag
            if (allocated(tag_id)) then
                deallocate(tag_id)
                allocate(tag_id(1:ntag))
            else
                allocate(tag_id(1:ntag))
            endif
            do itag=1,ntag
                read(58,'(I11)')tag_id(itag)
            end do
            allocate(order(1:ntag))
            write(*,*)'Sorting list of particle ids for binary search...'
#ifndef LONGINT
            call quick_sort_irg(tag_id,order,ntag)
#else
            call quick_sort_ilg(tag_id,order,ntag)
#endif
            deallocate(order)
            close(58)
        endif

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
            call cosmology_model
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
#ifndef IMASS
                allocate(part_tags(1:npart2))
#endif
            endif
            if (present(tag_file) .and. (.not. allocated(id))) allocate(id(1:npart2))
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
                if (sim%family) then
                    read(1) ! Skip family
#ifndef IMASS
                    read(1)part_tags
#else
                    read(1) ! Skip tags
#endif
                endif
                read(1)age
                read(1)met
#ifdef IMASS
                read(1)imass
#endif
            elseif (present(tag_file) .and. nstar .eq. 0) then
                read(1)id
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
#ifdef IMASS
                    part%imass = imass(i)
#else
                    part%imass = 0D0
                    if (part_tags(i)==1) then
                        part%imass = m(i)
                    elseif (part_tags(i)==0.or.part_tags(i)==-1) then
                        part%imass = m(i) / (1D0 - sim%eta_sn)
                    end if
#endif
                elseif (present(tag_file)) then
                    part%id = id(i)
                    part%age = 0D0
                    part%met = 0D0
                    part%imass = 0D0
                else
                    part%id = 0
                    part%age = 0D0
                    part%met = 0D0
                    part%imass = 0D0
                endif
                ! Check if particle is inside the desired region
                part%x = part%x - reg%centre
                call rotate_vector(part%x,trans_matrix)
                x(i,:) = part%x
                call checkifinside(x(i,:),reg,ok_part,distance)
                ok_filter = filter_particle(reg,filt,part)
                ok_part = ok_part.and.ok_filter
                ! Check if tags are present for particles
                if (present(tag_file) .and. ok_part) then
                    ok_tag = .false.      
#ifndef LONGINT
                    call binarysearch_irg(ntag,tag_id,part%id,ok_tag)
#else
                    call binarysearch_ilg(ntag,tag_id,part%id,ok_tag)
#endif
                    if (present(inverse_tag)) then
                        if (inverse_tag .and. ok_tag) ok_tag = .false.
                    end if
                    ok_part = ok_tag .and. ok_part
                endif
                if (ok_part) then
                    if (part%m < 2.842170943040401D-014) write(*,*)part%m
                    part%v = part%v -reg%bulk_velocity
                    call rotate_vector(part%v,trans_matrix)
                    binpos = 0
                    call findbinpos(reg,distance,part,prof_data,binpos)
                    if (binpos.ne.0)  call bindata(reg,part,prof_data,binpos)
                endif
            end do partloop
            deallocate(m,x,v)
            if (allocated(id))deallocate(id)
            if (nstar>0)deallocate(age,met,imass)
#ifndef IMASS
            if (nstar>0)deallocate(part_tags)
#endif
        end do cpuloop
    end subroutine get_parts_onedprofile

    subroutine onedprofile(repository,reg,filt,prof_data,lmax,logscale,tag_file,inverse_tag)
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        type(filter),intent(in) :: filt
        type(profile_handler),intent(inout) :: prof_data
        integer,intent(in) :: lmax
        logical,intent(in) :: logscale
        character(128),intent(in),optional :: tag_file
        logical,intent(in),optional :: inverse_tag

        call read_hydrofile_descriptor(repository)

        call init_amr_read(repository)
        amr%lmax = lmax
        if (lmax.eq.0) amr%lmax = amr%nlevelmax
        ! Check if particle data uses family
        if (sim%dm .and. sim%hydro) call check_families(repository)
        prof_data%xdata = 0D0
        prof_data%ydata = 0D0
        prof_data%logscale = logscale
        call makebins(reg,prof_data%xvarname,prof_data%nbins,prof_data%xdata,logscale)

        call get_cpu_map(reg)
        write(*,*)'ncpu_read:',amr%ncpu_read
        if (present(tag_file)) then
            if (present(inverse_tag)) then
                call get_parts_onedprofile(repository,reg,filt,prof_data,tag_file,inverse_tag)
            else
                call get_parts_onedprofile(repository,reg,filt,prof_data,tag_file)
            endif
        else
            call get_parts_onedprofile(repository,reg,filt,prof_data)
        endif
        call renormalise_bins(prof_data)
        if (logscale) prof_data%xdata = 10.**(prof_data%xdata)
    end subroutine onedprofile
end module part_profiles