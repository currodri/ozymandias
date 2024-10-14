module part_profiles
    use local
    use io_ramses
    use filtering
    use geometrical_regions
    use stats_utils
    use cosmology

    type profile_handler
        character(128) :: scaletype
        integer :: profdim,nfilter=1,zero_index
        character(128) :: xvarname
        integer :: nyvar
        character(128),dimension(:),allocatable :: yvarnames
        integer :: nbins
        integer :: nwvar
        integer :: nsubs=0
        character(128),dimension(:),allocatable :: wvarnames
        real(dbl) :: linthresh
        real(dbl),dimension(:),allocatable :: xdata
        type(pdf_handler),dimension(:),allocatable :: ydata
        type(region),dimension(:),allocatable :: subs
        type(filter),dimension(:),allocatable :: filters
    end type profile_handler

    contains
    subroutine allocate_profile_handler(prof)
        implicit none
        type(profile_handler),intent(inout) :: prof

        if (.not.allocated(prof%yvarnames)) allocate(prof%yvarnames(1:prof%nyvar))
        if (.not.allocated(prof%wvarnames)) allocate(prof%wvarnames(1:prof%nwvar))
        if (.not.allocated(prof%xdata)) allocate(prof%xdata(0:prof%nbins))
        if (.not.allocated(prof%ydata)) allocate(prof%ydata(1:prof%nbins))

        if (.not.allocated(prof%subs).and.(prof%nsubs>0)) allocate(prof%subs(1:prof%nsubs))
        if (.not.allocated(prof%filters)) allocate(prof%filters(1:prof%nfilter))
    end subroutine allocate_profile_handler

    subroutine bindata(reg,part,prof,trans_matrix,ifilt,ibin)
        use vectors
        use geometrical_regions
        implicit none
        type(region),intent(in) :: reg
        type(particle),intent(in) :: part
        type(profile_handler),intent(inout) :: prof
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix
        integer,intent(in) :: ifilt,ibin

        ! Local variables
        integer :: i,j,ipdf,index
        real(dbl) :: ytemp,wtemp,ytemp2
        character(128) :: tempvar,vartype,varname
        
        yvarloop: do i=1,prof%ydata(ibin)%nvars
            if (prof%ydata(ibin)%do_binning(i)) then
                ! Get variable
                call findbinpos_part(reg,part,ipdf,ytemp,trans_matrix,&
                                & prof%ydata(ibin)%scaletype(i),prof%ydata(ibin)%nbins,&
                                & prof%ydata(ibin)%bins(:,i),prof%ydata(ibin)%linthresh(i),&
                                & prof%ydata(ibin)%zero_index(i),prof%ydata(ibin)%varname(i))
                if (ytemp.eq.0d0) cycle
                ! Get min and max
                if (prof%ydata(ibin)%minv(i,ifilt).eq.0d0) then
                    prof%ydata(ibin)%minv(i,ifilt) = ytemp ! Just to make sure that the initial min is not zero
                else
                    prof%ydata(ibin)%minv(i,ifilt) = min(ytemp,prof%ydata(ibin)%minv(i,ifilt)) ! Min value
                end if
                prof%ydata(ibin)%maxv(i,ifilt) = max(ytemp,prof%ydata(ibin)%maxv(i,ifilt)) ! Max value
                prof%ydata(ibin)%nvalues(i,ifilt) = prof%ydata(ibin)%nvalues(i,ifilt) + 1

                if (ipdf.gt.0) then
                    wvarloop1: do j=1,prof%ydata(ibin)%nwvars
                        ! Get weights
                        tempvar = TRIM(prof%ydata(ibin)%wvarnames(j))
                        index = scan(tempvar,'/')
                        vartype = tempvar(1:index-1)
                        varname = tempvar(index+1:)
                        ytemp2 = ytemp
                        if (trim(varname)=='counts') then
                            wtemp =  1D0
                        else if (trim(varname)=='cumulative') then
                            wtemp = ytemp2
                        else
                            call getpartvalue(reg,part,prof%ydata(ibin)%wvarnames(j),&
                                            & wtemp)
                        endif

                        ! Save to PDFs
                        prof%ydata(ibin)%heights(i,ifilt,j,ipdf) = prof%ydata(ibin)%heights(i,ifilt,j,ipdf) + wtemp ! Weight to the PDF bin
                        prof%ydata(ibin)%totweights(i,ifilt,j) = prof%ydata(ibin)%totweights(i,ifilt,j) + wtemp       ! Weight

                        ! Now do it for the case of no binning (old integration method)
                        ! Get weights
                        if (trim(varname)=='counts') then
                            wtemp =  1D0
                            ytemp2 = 1D0
                        else if (trim(varname)=='cumulative') then
                            wtemp = 1D0
                        endif
                        
                        ! Save to attrs
                        prof%ydata(ibin)%total(i,ifilt,j,1) = prof%ydata(ibin)%total(i,ifilt,j,1) + ytemp2*wtemp ! Value (weighted or not)
                        prof%ydata(ibin)%total(i,ifilt,j,2) = prof%ydata(ibin)%total(i,ifilt,j,2) + wtemp       ! Weight
                    end do wvarloop1
                else
                    prof%ydata(ibin)%nout(i,ifilt) = prof%ydata(ibin)%nout(i,ifilt) + 1
                end if
            else
                ! Get variable
                call getpartvalue(reg,part,prof%ydata(ibin)%varname(i),ytemp)
                ! Get min and max
                if (prof%ydata(ibin)%minv(i,ifilt).eq.0D0) then
                    prof%ydata(ibin)%minv(i,ifilt) = ytemp ! Just to make sure that the initial min is not zero
                else
                    prof%ydata(ibin)%minv(i,ifilt) = min(ytemp,prof%ydata(ibin)%minv(i,ifilt))    ! Min value
                endif
                prof%ydata(ibin)%maxv(i,ifilt) = max(ytemp,prof%ydata(ibin)%maxv(i,ifilt))    ! Max value
                prof%ydata(ibin)%nvalues(i,ifilt) = prof%ydata(ibin)%nvalues(i,ifilt) + 1

                wvarloop2: do j=1,prof%ydata(ibin)%nwvars
                    ! Get weights
                    ytemp2 = ytemp
                    tempvar = TRIM(prof%ydata(ibin)%wvarnames(j))
                    index = scan(tempvar,'/')
                    vartype = tempvar(1:index-1)
                    varname = tempvar(index+1:)
                    if (trim(varname)=='counts') then
                        wtemp =  1D0
                        ytemp2 = 1D0
                    else if (trim(varname)=='cumulative') then
                        wtemp = 1D0
                    else
                        call getpartvalue(reg,part,prof%ydata(ibin)%wvarnames(j),wtemp)
                    endif
                    
                    ! Save to attrs
                    prof%ydata(ibin)%total(i,ifilt,j,1) = prof%ydata(ibin)%total(i,ifilt,j,1) + ytemp2*wtemp ! Value (weighted or not)
                    prof%ydata(ibin)%total(i,ifilt,j,2) = prof%ydata(ibin)%total(i,ifilt,j,2) + wtemp       ! Weight
                end do wvarloop2
            end if
        end do yvarloop
    end subroutine bindata

    subroutine renormalise_bins(prof)
        implicit none
        type(profile_handler),intent(inout) :: prof
        integer :: i,j,ibin,ifilt,index
        character(128) :: tempvar,vartype,varname

        binloop: do ibin=1,prof%nbins
            filterloop: do ifilt=1,prof%ydata(ibin)%nfilter
                varloop: do i=1,prof%ydata(ibin)%nvars
                    if (prof%ydata(ibin)%do_binning(i)) then
                        wvarloop1: do j=1,prof%ydata(ibin)%nwvars
                            tempvar = TRIM(prof%ydata(ibin)%wvarnames(j))
                            index = scan(tempvar,'/')
                            vartype = tempvar(1:index-1)
                            varname = tempvar(index+1:)
                            if (trim(varname) /= 'cumulative' .and. trim(varname) /= 'counts') then
                                prof%ydata(ibin)%heights(i,ifilt,j,:) = prof%ydata(ibin)%heights(i,ifilt,j,:) / prof%ydata(ibin)%totweights(i,ifilt,j)
                                prof%ydata(ibin)%total(i,ifilt,j,1) = prof%ydata(ibin)%total(i,ifilt,j,1) / prof%ydata(ibin)%total(i,ifilt,j,2)
                            endif
                        end do wvarloop1
                    else
                        wvarloop2: do j=1,prof%ydata(ibin)%nwvars
                            tempvar = TRIM(prof%ydata(ibin)%wvarnames(j))
                            index = scan(tempvar,'/')
                            vartype = tempvar(1:index-1)
                            varname = tempvar(index+1:)
                            if (trim(varname) /= 'cumulative' .and. trim(varname) /= 'counts') then
                                prof%ydata(ibin)%total(i,ifilt,j,1) = prof%ydata(ibin)%total(i,ifilt,j,1) / prof%ydata(ibin)%total(i,ifilt,j,2)
                            endif
                        end do wvarloop2
                    end if
                end do varloop
            end do filterloop
        end do binloop

        if (trim(prof%scaletype).eq.'log_even') prof%xdata = 10**prof%xdata
    end subroutine renormalise_bins

    subroutine get_parts_onedprofile(repository,reg,prof_data,tag_file,inverse_tag)
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
        type(profile_handler),intent(inout) :: prof_data
        character(128),intent(in),optional :: tag_file
        logical,intent(in),optional :: inverse_tag

        logical :: ok_part,ok_filter,ok_tag
        integer :: roterr
        integer :: i,j,k,itag,ifilt
        integer :: ipos,icpu,binpos
        integer :: npart,npart2,nstar,ntag
        integer :: ncpu2,ndim2
        real(dbl) :: distance,ytemp
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
            if (verbose) write(*,*)'Reading particle tags file '//TRIM(tag_file)
            read(58,'(I11)')ntag
            if (verbose) write(*,*)'Number of tagged particles in file: ',ntag
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
            if (verbose) write(*,*)'Sorting list of particle ids for binary search...'
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
            if (verbose) write(*,*)'Age simu=',sim%time_simu*sim%unit_t/(365.*24.*3600.*1d9)
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
        if (verbose) write(*,*)'Found ',npart,' particles.'
        if(nstar>0.and.verbose)then
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
                    do ifilt=1,prof_data%nfilter
                        ok_filter = filter_particle(reg,prof_data%filters(ifilt),part)
                        if (ok_filter) then
                            part%v = part%v -reg%bulk_velocity
                            call rotate_vector(part%v,trans_matrix)
                            binpos = 0
                            call findbinpos_part(reg,part,binpos,ytemp,&
                                                &trans_matrix,prof_data%scaletype,&
                                                &prof_data%nbins,prof_data%xdata,&
                                                &prof_data%linthresh,prof_data%zero_index,prof_data%xvarname)
                            if (binpos.ne.0)  call bindata(reg,part,prof_data,trans_matrix,ifilt,binpos)
                        end if
                    end do
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

    subroutine onedprofile(repository,reg,prof_data,lmax,tag_file,inverse_tag)
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        type(profile_handler),intent(inout) :: prof_data
        integer,intent(in) :: lmax
        character(128),intent(in),optional :: tag_file
        logical,intent(in),optional :: inverse_tag

        call read_hydrofile_descriptor(repository)

        call init_amr_read(repository)
        amr%lmax = lmax
        if (lmax.eq.0) amr%lmax = amr%nlevelmax
        ! Check if particle data uses family
        if (sim%dm .and. sim%hydro) call check_families(repository)

        call get_cpu_map(reg)
        if (verbose) write(*,*)'ncpu_read:',amr%ncpu_read
        if (present(tag_file)) then
            if (present(inverse_tag)) then
                call get_parts_onedprofile(repository,reg,prof_data,tag_file,inverse_tag)
            else
                call get_parts_onedprofile(repository,reg,prof_data,tag_file)
            endif
        else
            call get_parts_onedprofile(repository,reg,prof_data)
        endif
        call renormalise_bins(prof_data)
    end subroutine onedprofile
end module part_profiles