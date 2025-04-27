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
        type(part_var) :: xvar
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
        type(filter_part),dimension(:),allocatable :: filters
        type(part_var),dimension(:),allocatable :: yvars
        type(part_var),dimension(:),allocatable :: wvars
    end type profile_handler

    contains
    subroutine allocate_profile_handler(prof)
        implicit none
        type(profile_handler),intent(inout) :: prof

        if (.not.allocated(prof%yvarnames)) allocate(prof%yvarnames(1:prof%nyvar))
        if (.not.allocated(prof%wvarnames)) allocate(prof%wvarnames(1:prof%nwvar))
        if (.not.allocated(prof%yvars)) allocate(prof%yvars(1:prof%nyvar))
        if (.not.allocated(prof%wvars)) allocate(prof%wvars(1:prof%nwvar))
        if (.not.allocated(prof%xdata)) allocate(prof%xdata(0:prof%nbins))
        if (.not.allocated(prof%ydata)) allocate(prof%ydata(1:prof%nbins))

        if (.not.allocated(prof%subs).and.(prof%nsubs>0)) allocate(prof%subs(1:prof%nsubs))
        if (.not.allocated(prof%filters)) allocate(prof%filters(1:prof%nfilter))
    end subroutine allocate_profile_handler

    subroutine bindata(reg,part_data_d,part_data_i,part_data_b,prof,ibin,ifilt,trans_matrix)
        use vectors
        use geometrical_regions
        implicit none
        type(region),intent(in) :: reg
        real(dbl),dimension(1:sim%nvar_part_d),intent(in) :: part_data_d
#ifdef LONGINT
        integer(ilg),dimension(1:sim%nvar_part_i),intent(in) :: part_data_i
#else
        integer(irg),dimension(1:sim%nvar_part_i),intent(in) :: part_data_i
#endif
        integer(1),dimension(1:sim%nvar_part_b),intent(in) :: part_data_b
        type(profile_handler),intent(inout) :: prof
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix
        integer,intent(in) :: ifilt,ibin

        ! Local variables
        integer :: i,j,ipdf
        real(dbl) :: ytemp,wtemp,ytemp2
#ifdef LONGINT
        integer(ilg) :: wtemp_i,ytemp_i
#else
        integer(irg) :: wtemp_i,ytemp_i
#endif
        integer(1) :: wtemp_b,ytemp_b
        type(vector) :: dcell

        ! TODO: dcell is not used correctly, and should have some meaning in a profile
        
        yvarloop: do i=1,prof%ydata(ibin)%nvars
            if (prof%ydata(ibin)%do_binning(i)) then
                ! Get variable
                call findbinpos_part(reg,dcell,part_data_d,part_data_i,part_data_b,&
                                & ipdf,ytemp,trans_matrix,&
                                & prof%ydata(ibin)%scaletype(i),prof%ydata(ibin)%nbins,&
                                & prof%ydata(ibin)%bins(:,i),prof%ydata(ibin)%linthresh(i),&
                                & prof%ydata(ibin)%zero_index(i),prof%yvars(i))
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
                        ytemp2 = ytemp
                        if (TRIM(prof%ydata(ibin)%wvarnames(j))=='counts') then
                            wtemp =  1D0
                        else if (TRIM(prof%ydata(ibin)%wvarnames(j))=='cumulative') then
                            wtemp = ytemp2
                        else
                            if (prof%wvars(j)%vartype==1) then
                                wtemp = prof%wvars(j)%myfunction_d(amr,sim,prof%wvars(j),reg,dcell,&
                                        &part_data_d,part_data_i,part_data_b)
                            else if (prof%wvars(j)%vartype==2) then
                                wtemp_i = prof%wvars(j)%myfunction_i(amr,sim,prof%wvars(j),reg,dcell,&
                                        &part_data_d,part_data_i,part_data_b)
                                wtemp = real(wtemp_i,kind=dbl)
                            else if (prof%wvars(j)%vartype==3) then
                                wtemp_b = prof%wvars(j)%myfunction_b(amr,sim,prof%wvars(j),reg,dcell,&
                                        &part_data_d,part_data_i,part_data_b)
                                wtemp = real(wtemp_b,kind=dbl)
                            else
                                write(*,*)'Error in profile module: unknown variable type'
                                stop
                            end if
                        endif

                        ! Save to PDFs
                        prof%ydata(ibin)%heights(i,ifilt,j,ipdf) = prof%ydata(ibin)%heights(i,ifilt,j,ipdf) + wtemp ! Weight to the PDF bin
                        prof%ydata(ibin)%totweights(i,ifilt,j) = prof%ydata(ibin)%totweights(i,ifilt,j) + wtemp       ! Weight

                        ! Now do it for the case of no binning (old integration method)
                        ! Get weights
                        if (trim(prof%ydata(ibin)%wvarnames(j))=='counts') then
                            wtemp =  1D0
                            ytemp2 = 1D0
                        else if (trim(prof%ydata(ibin)%wvarnames(j))=='cumulative') then
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
                if (prof%yvars(i)%vartype==1) then
                    ytemp = prof%yvars(i)%myfunction_d(amr,sim,prof%yvars(i),reg,dcell,&
                                &part_data_d,part_data_i,part_data_b)
                else if (prof%yvars(i)%vartype==2) then
                    ytemp_i = prof%yvars(i)%myfunction_i(amr,sim,prof%yvars(i),reg,dcell,&
                                &part_data_d,part_data_i,part_data_b)
                    ytemp = real(ytemp_i,kind=dbl)
                else if (prof%yvars(i)%vartype==3) then
                    ytemp_b = prof%yvars(i)%myfunction_b(amr,sim,prof%yvars(i),reg,dcell,&
                                &part_data_d,part_data_i,part_data_b)
                    ytemp = real(ytemp_b,kind=dbl)
                else
                    write(*,*)'Error in profile module: unknown variable type'
                    stop
                end if

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
                    if (TRIM(prof%ydata(ibin)%wvarnames(j))=='counts') then
                        wtemp =  1D0
                        ytemp2 = 1D0
                    else if (TRIM(prof%ydata(ibin)%wvarnames(j))=='cumulative') then
                        wtemp = 1D0
                    else
                        if (prof%wvars(j)%vartype==1) then
                            wtemp = prof%wvars(j)%myfunction_d(amr,sim,prof%wvars(j),reg,dcell,&
                                    &part_data_d,part_data_i,part_data_b)
                        else if (prof%wvars(j)%vartype==2) then
                            wtemp_i = prof%wvars(j)%myfunction_i(amr,sim,prof%wvars(j),reg,dcell,&
                                    &part_data_d,part_data_i,part_data_b)
                            wtemp = real(wtemp_i,kind=dbl)
                        else if (prof%wvars(j)%vartype==3) then
                            wtemp_b = prof%wvars(j)%myfunction_b(amr,sim,prof%wvars(j),reg,dcell,&
                                    &part_data_d,part_data_i,part_data_b)
                            wtemp = real(wtemp_b,kind=dbl)
                        else
                            write(*,*)'Error in profile module: unknown variable type'
                            stop
                        end if
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
                            varname = TRIM(prof%ydata(ibin)%wvarnames(j))
                            if (trim(varname) /= 'cumulative' .and. trim(varname) /= 'counts') then
                                prof%ydata(ibin)%heights(i,ifilt,j,:) = prof%ydata(ibin)%heights(i,ifilt,j,:) / prof%ydata(ibin)%totweights(i,ifilt,j)
                                prof%ydata(ibin)%total(i,ifilt,j,1) = prof%ydata(ibin)%total(i,ifilt,j,1) / prof%ydata(ibin)%total(i,ifilt,j,2)
                            endif
                        end do wvarloop1
                    else
                        wvarloop2: do j=1,prof%ydata(ibin)%nwvars
                            varname = TRIM(prof%ydata(ibin)%wvarnames(j))
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

    subroutine onedprofile(repository,reg,prof_data,lmax,&
                          &part_dict,part_vtypes,tag_file,&
                          &inverse_tag)
        use geometrical_regions
        implicit none
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        type(profile_handler),intent(inout) :: prof_data
        integer,intent(in) :: lmax
        type(dictf90),intent(in),optional :: part_dict,part_vtypes
        character(128),intent(in),optional :: tag_file
        logical,intent(in),optional :: inverse_tag

        integer :: ii
        integer :: ivx,ivy,ivz

        ! Obtain details of the particle variables stored
        call read_partfile_descriptor(repository)

        ! Intialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = lmax
        if (lmax.eq.0) amr%lmax = amr%nlevelmax

        ! Compute the Hilbert curve
        call get_cpu_map(reg)

        if (prof_data%nsubs>0 .and. verbose) write(*,*)'Excluding substructures: ',prof_data%nsubs
        if (verbose) write(*,*)'ncpu_read:',amr%ncpu_read

        ! Set up the part variables quicklook tools
        if (present(part_dict).and.present(part_vtypes)) then
            ! If the user provides their own particle dictionary,
            ! use that one instead of the automatic from the
            ! particle_file_descriptor.txt (RAMSES)
            call get_partvar_tools(part_dict,part_vtypes,prof_data%nyvar,&
                                & prof_data%yvarnames,prof_data%yvars)

            ! Do it also for the weight variables
            call get_partvar_tools(part_dict,part_vtypes,prof_data%nwvar,&
                                & prof_data%wvarnames,prof_data%wvars)

            ! Do it for the xvar variable
            call set_part_var(part_dict,part_vtypes,prof_data%xvar)

            ! We also do it for the filter variables
            do ii = 1, prof_data%nfilter
                call get_filter_part_tools(part_dict,part_vtypes,prof_data%filters(ii))
            end do

            ! We always need the indexes of the velocities to perform rotation
            ! of particle velocities
            ivx = part_dict%get('velocity_x')
            ivy = part_dict%get('velocity_y')
            ivz = part_dict%get('velocity_z')

            ! If there exists a particle_file_descriptor.txt, make sure
            ! that the provided part_vtype is consistent with that one
            if (sim%isthere_part_descriptor) then
                do ii = 1, sim%nvar_part
                    if (sim%part_var_types(ii) /= part_vtypes%get(part_vtypes%keys(ii))) then
                        write(*,*)'Error: Provided particle dictionary is not consistent with the particle_file_descriptor.txt'
                        write(*,*)'sim%part_var_types: ',sim%part_var_types
                        write(*,*)sim%part_var_types(ii), part_vtypes%get(part_vtypes%keys(ii))
                        stop
                    end if
                end do
            end if

            !  Since all its alright, we just set the partIDs and partvar_types to the
            ! ones provided by the user
            partIDs = part_dict
            partvar_types = part_vtypes
        else
            ! If the user does not provide a dictionary, we just use the one
            ! provided by the particle_file_descriptor.txt (RAMSES)
            call get_partvar_tools(partIDs,partvar_types,prof_data%nyvar,&
                                & prof_data%yvarnames,prof_data%yvars)

            ! Do it also for the weight variables
            call get_partvar_tools(partIDs,partvar_types,prof_data%nwvar,&
                                & prof_data%wvarnames,prof_data%wvars)

            ! Do it for the xvar variable
            call set_part_var(partIDs,partvar_types,prof_data%xvar)

            ! We also do it for the filter variables
            do ii = 1, prof_data%nfilter
                call get_filter_part_tools(partIDs,partvar_types,prof_data%filters(ii))
            end do

            ! We always need the indexes of the velocities to perform rotation
            ! of particle velocities
            ivx = partIDs%get('velocity_x')
            ivy = partIDs%get('velocity_y')
            ivz = partIDs%get('velocity_z')
        end if

        ! Finally, renormalise the bins
        call renormalise_bins(prof_data)

        contains
            subroutine get_parts_onedprofile
                use utils, only: quick_sort_ilg, quick_sort_irg,\
                                binarysearch_ilg, binarysearch_irg
                implicit none

                logical :: ok_part,ok_filter,ok_tag,ok_sub
                integer :: roterr
                integer :: i,j,k,itag,ifilt,isub
                integer :: ipos,icpu,binpos
                integer :: npart,npart2,nstar,ntag
                integer :: ncpu2,ndim2
                integer :: npartsub
                real(dbl) :: distance,ytemp
                real(dbl),dimension(1:3,1:3) :: trans_matrix
                character(5) :: nchar,ncharcpu
                character(128) :: nomfich
                type(vector) :: xtemp,vtemp,dcell
                character(6) :: ptype
                integer,dimension(:),allocatable :: order
                integer,dimension(:),allocatable :: nparttoto
#ifndef LONGINT
                integer(irg),dimension(:),allocatable :: id,tag_id
                integer(irg),dimension(:,:),allocatable :: part_data_i
#else
                integer(ilg),dimension(:),allocatable :: id,tag_id
                integer(ilg),dimension(:,:),allocatable :: part_data_i
#endif
                real(dbl),dimension(:,:),allocatable :: part_data_d
                integer(1),dimension(:,:),allocatable :: part_data_b
                real(dbl),dimension(:,:),allocatable :: x,v

                npartsub = 0
                allocate(nparttoto(1:prof_data%nfilter))
                nparttoto = 0

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
                    
                    ! 1. Allocate particle data arrays
                    allocate(part_data_d(sim%nvar_part_d,1:npart2))
                    allocate(part_data_b(sim%nvar_part_b,1:npart2))
                    allocate(part_data_i(sim%nvar_part_i,1:npart2))

                    allocate(x(1:npart2,1:3))
                    allocate(v(1:npart2,1:3))
                    if (present(tag_file)) allocate(id(1:npart2))

                    ! 2. Loop over variables reading in the correct way as determined
                    ! by the particle_file_descriptor.txt
                    do i = 1, sim%nvar_part
                        if (sim%part_var_types(i) == 1) then
                            read(1) part_data_d(partIDs%get(partIDs%keys(i)),:)
                        elseif (sim%part_var_types(i) == 2) then
                            read(1) part_data_i(partIDs%get(partIDs%keys(i)),:)
                        elseif (sim%part_var_types(i) == 3) then
                            read(1) part_data_b(partIDs%get(partIDs%keys(i)),:)
                        endif
                    end do
                    close(1)

                    ! 3. The particle position needs to be in units of the boxlen
                    x(:,1) = part_data_d(partIDs%get('x'),:) / sim%boxlen
                    if (amr%ndim > 1) x(:,2) = part_data_d(partIDs%get('y'),:) / sim%boxlen
                    if (amr%ndim > 2) x(:,3) = part_data_d(partIDs%get('z'),:) / sim%boxlen

                    ! 4. Save also the particle velocities so they can be rotated
                    v(:,1) = part_data_d(partIDs%get('velocity_x'),:)
                    if (amr%ndim > 1) v(:,2) = part_data_d(partIDs%get('velocity_y'),:)
                    if (amr%ndim > 2) v(:,3) = part_data_d(partIDs%get('velocity_z'),:)

                    ! 5. If a tag file is used, get a hold of the particle ids
                    if (present(tag_file)) id(:) = part_data_i(partIDs%get('id'),:)

                    ! Project variables into map for particles
                    ! of interest in the region
                    partloop: do i=1,npart2
                        distance = 0D0
                        xtemp = x(i,:)
                        xtemp = xtemp - reg%centre
                        x(i,:) = xtemp
                        call checkifinside(x(i,:),reg,ok_part,distance)
                        xtemp = xtemp + reg%centre
                        x(i,:) = xtemp

                        ! If we are avoiding substructures, check if the particle is safe
                        if (prof_data%nsubs>0) then
                            ok_sub = .true.
                            do isub = 1,prof_data%nsubs
                                ok_sub = ok_sub .and. filter_sub(prof_data%subs(isub),x(i,:))
                            end do
                            if (.not.ok_sub) npartsub = npartsub + 1
                            ok_part = ok_part .and. ok_sub
                        end if

                        ! Check if tags are present for particles
                        if (present(tag_file) .and. ok_part) then
                            ok_tag = .false.      
#ifndef LONGINT
                            call binarysearch_irg(ntag,tag_id,id(i),ok_tag)
#else
                            call binarysearch_ilg(ntag,tag_id,id(i),ok_tag)
#endif
                            if (present(inverse_tag)) then
                                if (inverse_tag .and. ok_tag) ok_tag = .false.
                            end if
                            ok_part = ok_tag .and. ok_part
                        endif
                        if (ok_part) then
                            ! Rotate the particle velocity
                            vtemp = v(i,:)
                            vtemp = vtemp - reg%bulk_velocity
                            call rotate_vector(vtemp,trans_matrix)

                            ! Return the x array (now in boxlen units,
                            ! rotated for the region axis and centered)
                            xtemp = x(i,:)
                            xtemp = xtemp - reg%centre
                            call rotate_vector(xtemp,trans_matrix)
                            part_data_d(partIDs%get('x'),i) = xtemp%x
                            if (amr%ndim > 1) part_data_d(partIDs%get('y'),i) = xtemp%y
                            if (amr%ndim > 2) part_data_d(partIDs%get('z'),i) = xtemp%z

                            ! Return the velocity array (now rotated for the r
                            ! region axis and bulk velocity-corrected)
                            part_data_d(partIDs%get('velocity_x'),i) = vtemp%x
                            if (amr%ndim > 1) part_data_d(partIDs%get('velocity_y'),i) = vtemp%y
                            if (amr%ndim > 2) part_data_d(partIDs%get('velocity_z'),i) = vtemp%z

                            do ifilt=1,prof_data%nfilter
                                ok_filter = filter_particle(reg,prof_data%filters(ifilt),dcell,&
                                            &part_data_d(:,i),part_data_i(:,i),part_data_b(:,i))
                                if (ok_filter) then
                                    binpos = 0
                                    call findbinpos_part(reg,dcell,part_data_d(:,i),part_data_i(:,i),part_data_b(:,i),&
                                                        &binpos,ytemp,trans_matrix,prof_data%scaletype,&
                                                        &prof_data%nbins,prof_data%xdata,&
                                                        &prof_data%linthresh,prof_data%zero_index,prof_data%xvar)
                                    if (binpos.ne.0)  call bindata(reg,part_data_d(:,i),part_data_i(:,i),part_data_b(:,i),&
                                                        &prof_data,binpos,ifilt,trans_matrix)
                                    nparttoto(ifilt) = nparttoto(ifilt) + 1
                                end if
                            end do
                        endif
                    end do partloop
                    deallocate(x,v,part_data_d,part_data_i,part_data_b)
                    if (allocated(id)) deallocate(id)
                end do cpuloop
                if (verbose) write(*,*)'> nparttoto: ',nparttoto
                if (verbose) write(*,*)'> npartsub:  ',npartsub
                deallocate(nparttoto)
            end subroutine get_parts_onedprofile
    end subroutine onedprofile
end module part_profiles