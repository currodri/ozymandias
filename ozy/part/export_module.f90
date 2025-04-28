!--------------------------------------------------------------------------
! ozymandias:export_module.f90
!--------------------------------------------------------------------------
!
! MODULE: export_part
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> types and routines useful for extracting the raw particle data from
!> RAMSES into usable formats
!
!> @details  
!> 
! 
!
!> @date 26/07/2022   0.1.0 taking guidelines from original part2cube.f90
!--------------------------------------------------------------------------

module export_part
    use local
    use constants
    use io_ramses
    use filtering
    use cosmology

    contains

    subroutine part2skirt(repository,reg,filt,lmax,&
                        &h,smoothmethod,sedmethod,outpath,&
                        &part_dict,part_vtypes)
        use vectors
        use coordinate_systems
        use geometrical_regions

        implicit none
        ! Input/output variables
        character(128),intent(in) :: repository
        type(region),intent(in) :: reg
        type(filter_part),intent(inout) :: filt
        integer,intent(in) :: lmax
        real(dbl),intent(in) :: h
        character(100),intent(in) :: smoothmethod,sedmethod
        character(128),intent(in) :: outpath
        type(dictf90),intent(in),optional :: part_dict,part_vtypes

        ! Specific variables for this subroutine
        logical :: ok_part,ok_filter
        logical :: have_imass
        integer :: roterr
        integer :: i,j,k
        integer :: ipos,icpu,binpos
        integer :: npart,npart2,nstar,inpart=0
        integer :: ncpu2,ndim2
        integer :: nstarsaved
        real(dbl) :: distance,tempage,dx,particle_h
        character(5) :: nchar,ncharcpu
        character(6) :: ptype
        character(128) :: nomfich,varname
        character(128),dimension(1) :: skirt_varnames
        type(part_var),dimension(1) :: skirt_vars
        type(vector) :: xtemp,vtemp,dcell
        type(particle) :: part
#ifndef LONGINT
        integer(irg),dimension(:,:),allocatable :: part_data_i
#else
        integer(ilg),dimension(:,:),allocatable :: part_data_i
#endif
        real(dbl),dimension(:,:),allocatable :: part_data_d
        integer(1),dimension(:,:),allocatable :: part_data_b
        real(dbl),dimension(:,:),allocatable :: x,v

        integer :: ii
        integer :: ivx,ivy,ivz

        ! Obtain details of the hydro variables stored
        call read_hydrofile_descriptor(repository)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = lmax
        if (lmax.eq.0) amr%lmax = amr%nlevelmax

        ! Compute the Hilbert curve
        call get_cpu_map(reg)
        if (verbose) write(*,*)'ncpu_read:',amr%ncpu_read

        ! Set up the part variables quicklook tools
        ! (This is the case for SKIRT)
        skirt_varnames(1) = 'age'
        if (present(part_dict).and.present(part_vtypes)) then
            ! If the user provides their own particle dictionary,
            ! use that one instead of the automatic from the
            ! particle_file_descriptor.txt (RAMSES)
            if (part_dict%get('imass').eq.0) then
                have_imass = .false.
            else
                have_imass = .true.
            end if
            call get_partvar_tools(part_dict,part_vtypes,1,skirt_varnames,skirt_vars)

            ! Now get it for the filter
            call get_filter_part_tools(part_dict,part_vtypes,filt)

            ! We always need the indexes od the velocities to perform rotation
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
            ! from the particle_file_descriptor.txt (RAMSES)
            if (partIDs%get('imass').eq.0) then
                have_imass = .false.
            else
                have_imass = .true.
            end if
            call get_partvar_tools(partIDs,partvar_types,1,skirt_varnames,skirt_vars)

            ! Now get it for the filter
            call get_filter_part_tools(partIDs,partvar_types,filt)

            ! We always need the indexes od the velocities to perform rotation
            ivx = partIDs%get('velocity_x')
            ivy = partIDs%get('velocity_y')
            ivz = partIDs%get('velocity_z')
        end if

#ifndef IMASS
        if (sim%eta_sn .eq. -1D0) then
            write(*,*)': eta_sn=-1 and not IMASS --> should set this up!'
            stop
        end if
#endif

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

        nstarsaved = 0
        ! Open output file and add header for SKIRT format
        open(unit=7,file=TRIM(outpath),form='formatted')
        write(7,98)
        write(7,99)TRIM(sedmethod)
        write(7,97)
        write(7,101)
        write(7,97)
        97 format('#')
        98 format('# Stellar particle for simulated galaxy in RAMSES simulation')
        99 format('# SKIRT 9 import format for a particle source with the ',A,' SED family')
        101 format('# Column 1: x-coordinate (kpc)',/, &
                    '# Column 2: y-coordinate (kpc)',/,&
                    '# Column 3: z-coordinate (kpc)',/,&
                    '# Column 4: smoothing length (pc)',/,&
                    '# Column 5: x-velocity (km/s)',/, &
                    '# Column 6: y-velocity (km/s)',/,&
                    '# Column 7: z-velocity (km/s)',/,&
                    '# Column 8: initial mass (Msun)',/,&
                    '# Column 9: metallicity (1)',/,&
                    '# Column 10: age (Gyr)')

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

            ! Get variable info for particles in the 
            ! region of interest
            partloop: do i=1,npart2
                distance = 0D0
                xtemp = x(i,:)
                xtemp = xtemp - reg%centre
                x(i,:) = xtemp
                call checkifinside(x(i,:),reg,ok_part,distance)
                xtemp = xtemp + reg%centre
                x(i,:) = xtemp
                ok_filter = filter_particle(reg,filt,dcell,part_data_d(:,i),&
                                            &part_data_i(:,i),part_data_b(:,i))
                if (ok_part.and.ok_filter) then
                    nstarsaved = nstarsaved + 1
                    ! Position to kpc
                    part%x = x(i,:) * (sim%unit_l*cm2kpc)
                    ! Velocity with respect to COM of galaxy and in km/s
                    vtemp = v(i,:)
                    part%v = vtemp - reg%bulk_velocity
                    part%v = part%v * (sim%unit_v*cm2km)
                    ! Initial mass in Msun
                    if (have_imass) then
                        part%imass = part_data_d(partIDs%get('imass'),i)
                    else
                        part%imass = part_data_d(partIDs%get('mass'),i) / (1d0 - sim%eta_sn)
                    end if
                    part%imass = part%imass * (sim%unit_m*g2msun)
                    ! Particle metallicity
                    part%met = part_data_d(partIDs%get('metallicity'),i)
                    ! Particle age in Gyr
                    part%age = skirt_vars(1)%myfunction_d(amr,sim,skirt_vars(1),&
                                                &reg,dcell,part_data_d(:,i),&
                                                &part_data_i(:,i),part_data_b(:,i))
                    ! Particle level
                    part%level = part_data_i(partIDs%get('level'),i)
                    ! Get smoothing length with the choosen method
                    select case (TRIM(smoothmethod))
                        case ('constant')
                            particle_h = H
                        case ('level')
                            particle_h = (1D0/(2D0**dble(part%level)))*sim%unit_l*cm2pc
                    end select
                    write(7,100)part%x%x,part%x%y,part%x%z,&
                                particle_h,&
                                part%v%x,part%v%y,part%v%z,&
                                part%imass,&
                                part%met,&
                                part%age
                    100 format(3F10.6,F10.2,4F10.2,F10.6,F10.6)
                endif
            end do partloop
            deallocate(x,v,part_data_d,part_data_i,part_data_b)
            inpart = inpart + npart2
        end do cpuloop

        close(7)
        write(*,102)nstarsaved
        102 format('File includes ',I12,' star particles')
    end subroutine part2skirt

    subroutine part2disperse(repository,reg,filt,lmax,prob,outpath,&
                            &part_dict,part_vtypes)
#ifndef NOIFPORT
        use IFPORT
#else
#warning Compiling without IFPORT
#endif
        use vectors
        use coordinate_systems
        use geometrical_regions

        implicit none
        ! Input/output variables
        character(128),intent(in) :: repository
        type(region),intent(in) :: reg
        type(filter_part),intent(inout) :: filt
        integer,intent(in) :: lmax
        real(dbl),intent(in)    :: prob
        character(128),intent(in) :: outpath
        type(dictf90),intent(in),optional :: part_dict,part_vtypes

        ! Specific variables for this subroutine
        logical :: ok_part,ok_filter
        integer :: roterr
        integer :: i,j,k
        integer :: ipos,icpu,binpos
        integer :: npart,npart2,nstar,inpart=0
        integer :: ncpu2,ndim2
        integer :: npartsaved
        real(dbl) :: distance,tempage,dx
        character(5) :: nchar,ncharcpu
        character(6) :: ptype
        character(128) :: nomfich,varname
        type(vector) :: xtemp,vtemp,dcell
        type(particle) :: part
#ifndef LONGINT
        integer(irg),dimension(:,:),allocatable :: part_data_i
#else
        integer(ilg),dimension(:,:),allocatable :: part_data_i
#endif
        real(dbl),dimension(:,:),allocatable :: part_data_d
        integer(1),dimension(:,:),allocatable :: part_data_b
        real(dbl),dimension(:,:),allocatable :: x,v

        integer :: ii
        integer :: ivx,ivy,ivz
        real(dbl) :: partprob

        integer,parameter :: myseed = 86456

        ! Initialise random seed
        call srand(myseed)

        ! Obtain details of the hydro variables stored
        call read_hydrofile_descriptor(repository)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = lmax
        if (lmax.eq.0) amr%lmax = amr%nlevelmax

        ! Check if particle data uses family
        if (sim%dm .and. sim%hydro) call check_families(repository)

        ! Compute the Hilbert curve
        call get_cpu_map(reg)
        if (verbose) write(*,*)'ncpu_read:',amr%ncpu_read

        ! Cosmological model
        if (sim%aexp.eq.1.and.sim%h0.eq.1)sim%cosmo=.false.
        if (sim%cosmo) then
            call cosmology_model
        else
            sim%time_simu = sim%t
            write(*,*)'Age simu=',sim%time_simu*sim%unit_t/(365.*24.*3600.*1d9)
        endif

        if (present(part_dict).and.present(part_vtypes)) then
            ! If the user provides their own particle dictionary,
            ! use that one instead of the automatic from the
            ! particle_file_descriptor.txt (RAMSES)

            ! Get it for the filter
            call get_filter_part_tools(part_dict,part_vtypes,filt)

            ! We always need the indexes od the velocities to perform rotation
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
            ! Now get it for the filter
            call get_filter_part_tools(partIDs,partvar_types,filt)

            ! We always need the indexes od the velocities to perform rotation
            ivx = partIDs%get('velocity_x')
            ivy = partIDs%get('velocity_y')
            ivz = partIDs%get('velocity_z')
        end if

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

        npartsaved = 0
        ! Open output file and add header for DISPERSE
        open(unit=7,file=TRIM(outpath),form='formatted')
        write(7,97)prob
        write(7,98)
        97 format('# DM particle positions for DISPERSE using random sampling with probability ',F10.2)
        98 format('# px py pz [in kpc]')

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

            ! Get variable info for particles in the 
            ! region of interest
            partloop: do i=1,npart2
                distance = 0D0
                xtemp = x(i,:)
                xtemp = xtemp - reg%centre
                x(i,:) = xtemp
                call checkifinside(x(i,:),reg,ok_part,distance)
                xtemp = xtemp + reg%centre
                x(i,:) = xtemp
                ok_filter = filter_particle(reg,filt,dcell,part_data_d(:,i),&
                                            &part_data_i(:,i),part_data_b(:,i))
                partprob = rand()
                if (ok_part.and.ok_filter.and.(partprob.le.prob)) then
                    npartsaved = npartsaved + 1
                    ! Position to kpc
                    part%x = x(i,:) * (sim%unit_l*cm2kpc)
                    part%x = part%x * (sim%unit_l*cm2kpc)
                    write(7,100)part%x%x,part%x%y,part%x%z
                    100 format(3F16.7)
                endif
            end do partloop
            inpart = inpart + npart2
            deallocate(x,v,part_data_d,part_data_i,part_data_b)
        end do cpuloop
        
        close(7)
        write(*,102)npartsaved
        102 format('File includes ',I12,' DM particles')
    end subroutine part2disperse

    subroutine part2file(repository,reg,filt,lmax,h,smoothmethod,outpath,&
                        &part_dict,part_vtypes)
        use vectors
        use coordinate_systems
        use geometrical_regions

        implicit none
        ! Input/output variables
        character(128),intent(in) :: repository
        type(region),intent(in) :: reg
        type(filter_part),intent(inout) :: filt
        integer,intent(in) :: lmax
        real(dbl),intent(in) :: h
        character(100),intent(in) :: smoothmethod
        character(128),intent(in) :: outpath
        type(dictf90),intent(in),optional :: part_dict,part_vtypes

        ! Specific variables for this subroutine
        logical :: ok_part,ok_filter
        logical :: have_imass
        integer :: roterr
        integer :: i,j,k
        integer :: ipos,icpu,binpos
        integer :: npart,npart2,nstar,inpart=0
        integer :: ncpu2,ndim2
        integer :: nstarsaved
        real(dbl) :: distance,tempage,dx,particle_h
        character(5) :: nchar,ncharcpu
        character(6) :: ptype
        character(128) :: nomfich,varname
        character(128),dimension(1) :: out_varnames
        type(part_var),dimension(1) :: out_vars
        type(vector) :: xtemp,vtemp,dcell
        type(particle) :: part
#ifndef LONGINT
        integer(irg),dimension(:,:),allocatable :: part_data_i
#else
        integer(ilg),dimension(:,:),allocatable :: part_data_i
#endif
        real(dbl),dimension(:,:),allocatable :: part_data_d
        integer(1),dimension(:,:),allocatable :: part_data_b
        real(dbl),dimension(:,:),allocatable :: x,v

        integer :: ii
        integer :: ivx,ivy,ivz

        ! Obtain details of the hydro variables stored
        call read_hydrofile_descriptor(repository)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = lmax
        if (lmax.eq.0) amr%lmax = amr%nlevelmax

        ! Compute the Hilbert curve
        call get_cpu_map(reg)
        if (verbose) write(*,*)'ncpu_read:',amr%ncpu_read

#ifndef IMASS
        if (sim%eta_sn .eq. -1D0) then
            write(*,*)': eta_sn=-1 and not IMASS --> should set this up!'
            stop
        end if
#endif

        ! Cosmological model
        if (sim%aexp.eq.1.and.sim%h0.eq.1)sim%cosmo=.false.
        if (sim%cosmo) then
            call cosmology_model
        else
            sim%time_simu = sim%t
            write(*,*)'Age simu=',sim%time_simu*sim%unit_t/(365.*24.*3600.*1d9)
        endif

        ! Set up the part variables quicklook tools
        ! (This is the case for BASIC)
        out_varnames(1) = 'age'
        if (present(part_dict).and.present(part_vtypes)) then
            ! If the user provides their own particle dictionary,
            ! use that one instead of the automatic from the
            ! particle_file_descriptor.txt (RAMSES)
            if (part_dict%get('imass').eq.0) then
                have_imass = .false.
            else
                have_imass = .true.
            end if
            call get_partvar_tools(part_dict,part_vtypes,1,out_varnames,out_vars)

            ! Now get it for the filter
            call get_filter_part_tools(part_dict,part_vtypes,filt)

            ! We always need the indexes od the velocities to perform rotation
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
            ! from the particle_file_descriptor.txt (RAMSES)
            if (partIDs%get('imass').eq.0) then
                have_imass = .false.
            else
                have_imass = .true.
            end if
            call get_partvar_tools(partIDs,partvar_types,1,out_varnames,out_vars)

            ! Now get it for the filter
            call get_filter_part_tools(partIDs,partvar_types,filt)

            ! We always need the indexes od the velocities to perform rotation
            ivx = partIDs%get('velocity_x')
            ivy = partIDs%get('velocity_y')
            ivz = partIDs%get('velocity_z')
        end if

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

        nstarsaved = 0
        ! Open output file and add header for SKIRT format
        open(unit=7,file=TRIM(outpath),form='formatted')
        write(7,98)
        write(7,99)
        write(7,97)
        write(7,101)
        write(7,97)
        write(7, '(A11,A11,A11,A11,A11,A11,A11,A11,A11,A11)') &
                    'x', 'y', 'z', 'h', 'vx', 'vy', 'vz', 'imass', 'met', 'age'
        97 format('#')
        98 format('# Stellar particle for simulated galaxy in RAMSES simulation')
        99 format('# BASIC particle output:')
        101 format('# Column 1: x-coordinate (kpc)',/, &
                    '# Column 2: y-coordinate (kpc)',/,&
                    '# Column 3: z-coordinate (kpc)',/,&
                    '# Column 4: smoothing length (pc)',/,&
                    '# Column 5: x-velocity (km/s)',/, &
                    '# Column 6: y-velocity (km/s)',/,&
                    '# Column 7: z-velocity (km/s)',/,&
                    '# Column 8: initial mass (Msun)',/,&
                    '# Column 9: metallicity (1)',/,&
                    '# Column 10: age (Gyr)')

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

            ! Get variable info for particles in the 
            ! region of interest
            partloop: do i=1,npart2
                distance = 0D0
                xtemp = x(i,:)
                xtemp = xtemp - reg%centre
                x(i,:) = xtemp
                call checkifinside(x(i,:),reg,ok_part,distance)
                xtemp = xtemp + reg%centre
                x(i,:) = xtemp
                ok_filter = filter_particle(reg,filt,dcell,part_data_d(:,i),&
                                        &part_data_i(:,i),part_data_b(:,i))
                if (ok_part.and.ok_filter) then
                    nstarsaved = nstarsaved + 1
                    ! Position to kpc
                    part%x = x(i,:) * (sim%unit_l*cm2kpc)
                    ! Velocity with respect to COM of galaxy and in km/s
                    vtemp = v(i,:)
                    part%v = vtemp - reg%bulk_velocity
                    part%v = part%v * (sim%unit_v*cm2km)
                    ! Initial mass in Msun
                    if (have_imass) then
                        part%imass = part_data_d(partIDs%get('imass'),i)
                    else
                        part%imass = part_data_d(partIDs%get('mass'),i) / (1d0 - sim%eta_sn)
                    end if
                    part%imass = part%imass * (sim%unit_m*g2msun)
                    ! Particle metallicity
                    part%met = part_data_d(partIDs%get('metallicity'),i)
                    ! Particle age in Gyr
                    part%age = out_vars(1)%myfunction_d(amr,sim,out_vars(1),&
                                                &reg,dcell,part_data_d(:,i),&
                                                &part_data_i(:,i),part_data_b(:,i))
                    ! Particle level
                    part%level = part_data_i(partIDs%get('level'),i)
                    ! Get smoothing length with the choosen method
                    select case (TRIM(smoothmethod))
                        case ('constant')
                            particle_h = H
                        case ('level')
                            particle_h = (1D0/(2D0**dble(part%level)))*sim%unit_l*cm2pc
                    end select
                    write(7,100)part%x%x,part%x%y,part%x%z,&
                                particle_h,&
                                part%v%x,part%v%y,part%v%z,&
                                part%imass,&
                                part%met,&
                                part%age
                    100 format(3F11.6,F11.2,4F11.2,F11.6,F11.6)
                endif
            end do partloop
            deallocate(x,v,part_data_d,part_data_i,part_data_b)
            inpart = inpart + npart2
        end do cpuloop

        close(7)
        write(*,102)nstarsaved
        102 format('File includes ',I12,' star particles')
    end subroutine part2file

end module export_part