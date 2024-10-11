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
    use utils
    use io_ramses
    use filtering
    use cosmology

    contains

    subroutine read_ssp_mrelease(filename, time, Mtot)
        character(len=*), intent(in) :: filename
        real(dbl), dimension(:), allocatable, intent(out) :: time, Mtot
        integer :: num_lines
        integer :: iunit, i
        real(dbl) :: dummy_line

        ! Open the file
        write(*,*) 'Reading mass release for fake SSP particles from ',trim(filename)
        open(newunit=iunit, file=filename, status='old', action='read')

        ! Skip header lines
        do i = 1, 3
            read(iunit, *)
        end do

        ! Count the number of lines
        num_lines = 0
        do
            read(iunit, *, iostat=i)
            if (i /= 0) exit
            num_lines = num_lines + 1
        end do

        ! Rewind back to the beginning of the file
        rewind(iunit)

        ! Allocate memory
        allocate(time(num_lines), Mtot(num_lines))

        ! Skip header lines again
        do i = 1, 3
            read(iunit, *)
        end do

        ! Read time and Mtot values
        do i = 1, num_lines
            read(iunit, *) time(i), Mtot(i), dummy_line ! Ignore other columns
        end do

        ! Convert to logarithmic
        time(:) = log10(time(:))
        Mtot(:) = log10(Mtot(:))

        ! Close the file
        close(iunit)
    end subroutine read_ssp_mrelease


    subroutine part2skirt(repository,reg,filt,h,smoothmethod,sedmethod,outpath,ssp_fakestars)
#ifndef NOIFPORT
        use IFPORT
#else
#warning Compiling without IFPORT
#endif
        use utils
        use vectors
        use coordinate_systems
        use geometrical_regions

        implicit none
        ! Input/output variables
        character(128),intent(in) :: repository
        type(region),intent(in) :: reg
        type(filter),intent(in) :: filt
        real(dbl),intent(in) :: h
        character(100),intent(in) :: smoothmethod,sedmethod
        character(128),intent(in) :: outpath
        character(128),intent(in),optional :: ssp_fakestars

        ! Specific variables for this subroutine
        logical :: ok_part,ok_filter,fake_star
        integer :: roterr
        integer :: i,j,k
        integer :: ipos,icpu,binpos
        integer :: npart,npart2,nstar,inpart=0
        integer :: ncpu2,ndim2
        integer :: nstarsaved,nfakestars
        real(dbl) :: distance,tempage,dx,rcyl
        real(dbl),dimension(1:3) :: xpos
        character(5) :: nchar,ncharcpu
        character(6) :: ptype
        character(128) :: nomfich,varname
        type(vector) :: xtemp,vtemp
        type(particle) :: part
        real(dbl),dimension(:,:),allocatable :: part_data_d
        integer(1),dimension(:,:),allocatable :: part_data_b
#ifdef LONGINT
        integer(ilg),dimension(:,:),allocatable :: part_data_i
#else
        integer(irg),dimension(:,:),allocatable :: part_data_i
#endif
        real(dbl) :: particle_h,mreturn
        real(dbl),dimension(:),allocatable :: ssp_time,ssp_mreturn

        ! Obtain details of the hydro variables stored
        call read_hydrofile_descriptor(repository)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = amr%nlevelmax

        ! Check if particle data uses family
        if (sim%dm .and. sim%hydro) call check_families(repository)

        ! Read the format of the particle data stored
        call read_partfile_descriptor(repository)
        

#ifdef LONGINT
        write(*,*) 'Using LONGINT for particle IDs'
#endif
#ifndef IMASS
        if (sim%eta_sn .eq. -1D0) then
            write(*,*)': eta_sn=-1 and not IMASS --> should set this up!'
            stop
        end if
#endif

        ! Compute the Hilbert curve
        call get_cpu_map(reg)
        write(*,*)'ncpu_read:',amr%ncpu_read

        ! Cosmological model
        if (sim%aexp.eq.1.and.sim%h0.eq.1)sim%cosmo=.false.
        if (sim%cosmo) then
            call cosmology_model
        else
            sim%time_simu = sim%t
            write(*,*)'Age simu=',sim%time_simu*sim%unit_t/(365.*24.*3600.*1d9)
        endif

        ! If a fake SSP model has been given (this is for isolated sims)
        ! then we will read the SSP mass return evolution with time
        ! in order to compute "fake" initial mass for the DM particles
        ! with 0 age but mass below the minimum DM mass
        if (present(ssp_fakestars)) then
            call read_ssp_mrelease(ssp_fakestars,ssp_time,ssp_mreturn)
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
        nfakestars = 0
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

            ! Allocate particle data arrays
            allocate(part_data_d(partIDs%nd,npart2))
            allocate(part_data_i(partIDs%ni,npart2))
            allocate(part_data_b(partIDs%nb,npart2))

            ! Loop over variables reading in the correct way as determined
            ! by the pvar_info details
            do i=1,partIDs%nvar
                if (partIDs%pvar_infos(i)%variable_type=='d') then
                    read(1) part_data_d(partIDs%pvar_infos(i)%ipos,:)
                elseif (partIDs%pvar_infos(i)%variable_type=='i') then
                    read(1) part_data_i(partIDs%pvar_infos(i)%ipos,:)
                elseif (partIDs%pvar_infos(i)%variable_type=='b') then
                    read(1) part_data_b(partIDs%pvar_infos(i)%ipos,:)
                end if
            end do
            close(1)

            ! Get variable info for particles in the 
            ! region of interest
            partloop: do i=1,npart2
                distance = 0D0
                fake_star = .false.
                part%x = (/part_data_d(get_ipos(partIDs,partIDs%position_x),i),&
                            part_data_d(get_ipos(partIDs,partIDs%position_y),i),&
                            part_data_d(get_ipos(partIDs,partIDs%position_z),i)/)
                part%x = part%x / sim%boxlen
                part%x = part%x - reg%centre
                part%v = (/part_data_d(get_ipos(partIDs,partIDs%velocity_x),i),&
                            part_data_d(get_ipos(partIDs,partIDs%velocity_y),i),&
                            part_data_d(get_ipos(partIDs,partIDs%velocity_z),i)/)
                part%m = part_data_d(get_ipos(partIDs,partIDs%mass),i)
                part%id = part_data_i(get_ipos(partIDs,partIDs%identity),i)
                part%level = part_data_i(get_ipos(partIDs,partIDs%levelp),i)
                part%birth_time = part_data_d(get_ipos(partIDs,partIDs%birth_time),i)
                part%met = part_data_d(get_ipos(partIDs,partIDs%metallicity),i)
#ifdef IMASS
                part%imass = part_data_d(get_ipos(partIDs,partIDs%initial_mass),i)
                if ((part%birth_time.eq.0d0).and.(part%m<0.00102).and.(present(ssp_fakestars))) then
                    ! In the case we encounter a fake star particle (from the ICs)
                    ! we need to compute their fake age and fake initial mass
                    ! assuming a particular SSP mass return

                    ! 1. Obtain fake age given by their cylindrical radius (matching
                    ! the age of the Sun at its galactocentric distance and with a
                    ! random dispersion of 3 Gyr)
                    varname = 'dm/r_cyl'
                    call getpartvalue(reg,part,varname,tempage)
                    rcyl = tempage * (sim%unit_l*sim%boxlen*cm2kpc)
                    part%age = max(0d0,min((9d0 - 0.5d0 * rcyl) + random_gaussian(0d0,3d0),14d0)) ! in Gyr

                    ! 2. Interpolate from the age what the mass return should be
                    mreturn = interpolate_log(ssp_time,ssp_mreturn,part%age*1d9)
                    part%imass = part%m / (1d0 - mreturn)
                    ! 3. Set the metallicity fixed to the solar metallicity (Asplund 2009)
                    part%met = 0.01345d0
                    fake_star = .true.
                end if
#else
                part%imass = 0D0
                if ((part%birth_time.eq.0d0).and.(part%m<0.00102).and.(present(ssp_fakestars))) then
                    ! In the case we encounter a fake star particle (from the ICs)
                    ! we need to compute their fake age and fake initial mass
                    ! assuming a particular SSP mass return

                    ! 1. Obtain fake age given by their cylindrical radius (matching
                    ! the age of the Sun at its galactocentric distance and with a
                    ! random dispersion of 3 Gyr)
                    varname = 'dm/r_cyl'
                    call getpartvalue(reg,part,varname,tempage)
                    rcyl = tempage * (sim%unit_l*sim%boxlen*cm2kpc)
                    part%age = max(0d0,min((9d0 - 0.5d0 * rcyl) + random_gaussian(0d0,3d0),14d0)) ! in Gyr

                    ! 2. Interpolate from the age what the mass return should be
                    mreturn = interpolate_log(ssp_time,ssp_mreturn,part%age*1d9)
                    part%imass = part%m / (1d0 - mreturn)

                    ! 3. Set the metallicity fixed to the solar metallicity (Asplund 2009)
                    part%met = 0.01345d0
                    fake_star = .true.
                else
                    if (sim%family) then
                        part%family = part_data_b(get_ipos(partIDs,partIDs%family),i)
                        part%tag = part_data_b(get_ipos(partIDs,partIDs%tag),i)
                        if (part%tag==1) then
                            part%imass = part%m
                        elseif (part%tag==0.or.part%tag==-1) then
                            part%imass = part%m / (1D0 - sim%eta_sn)
                        end if
                    else
                        part%imass = part%m / (1D0 - sim%eta_sn)
                    end if
                end if
#endif
                ! Check if particle is inside the desired region
                xpos = part%x
                call checkifinside(xpos,reg,ok_part,distance)
                ! ok_filter = filter_particle(reg,filt,part)
                ! ok_part = ok_part.and.ok_filter
                ! call getparttype(part,ptype)
                ! If particle mass is less than 1e6 Msun (0.00102 in code mass units)
                ! it is consider a star particle in isolated galaxy simulations
                if (part%imass>0d0.and.ok_part) then
                    if (fake_star) then
                        nfakestars = nfakestars + 1
                    else
                        nstarsaved = nstarsaved + 1
                    end if
                    ! Position to kpc
                    part%x = part%x * (sim%unit_l*sim%boxlen*cm2kpc)
                    ! Velocity with respect to COM of galaxy and in km/s
                    part%v = part%v - reg%bulk_velocity
                    part%v = part%v * (sim%unit_v*cm2km)
                    ! Initial mass in Msun
                    part%imass = part%imass * (sim%unit_m*g2msun)
                    ! Particle age in Gyr
                    if (.not.fake_star) then
                        varname = 'star/age'
                        call getpartvalue(reg,part,varname,tempage)
                        part%age = tempage
                    end if
                    ! Get smoothing length with the choosen method
                    select case (TRIM(smoothmethod))
                        case ('constant')
                            particle_h = H
                        case ('level')
                            particle_h = (1D0/(2**part%level))*sim%unit_l*cm2pc*sim%boxlen
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
            deallocate(part_data_d,part_data_i,part_data_b)
            inpart = inpart + npart2
        end do cpuloop

        close(7)
        write(*,102)nstarsaved
        write(*,103)nfakestars
        102 format('File includes ',I12,' star particles')
        103 format('File includes ',I12,' fake star particles')
    end subroutine part2skirt

    subroutine part2disperse(repository,reg,filt,prob,outpath)
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
        type(filter),intent(in) :: filt
        real(dbl),intent(in)    :: prob
        character(128),intent(in) :: outpath

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
        type(vector) :: xtemp,vtemp
        type(particle) :: part
        integer,dimension(:),allocatable :: plevel
        real(dbl),dimension(:),allocatable :: m,age,met,imass
        real(dbl),dimension(:,:),allocatable :: x,v
        real(dbl) :: partprob

        integer,parameter :: myseed = 86456

        ! Initialise random seed
        call srand(myseed)

        ! Obtain details of the hydro variables stored
        call read_hydrofile_descriptor(repository)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = amr%nlevelmax

        ! Check if particle data uses family
        if (sim%dm .and. sim%hydro) call check_families(repository)

        ! Compute the Hilbert curve
        call get_cpu_map(reg)
        write(*,*)'ncpu_read:',amr%ncpu_read

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
            allocate(m(1:npart2))
            if(nstar>0)then
                allocate(age(1:npart2))
            endif
            allocate(x(1:npart2,1:ndim2))

            ! Read position
            do i=1,amr%ndim
                read(1)m
                x(1:npart2,i) = m/sim%boxlen
            end do

            ! Read velocity
            do i=1,amr%ndim
                read(1)m
            end do

            ! Read mass
            read(1)m
            if (nstar>0) then
                read(1) ! Skip id
                read(1) ! Skip level
                if (sim%family) then
                    read(1) ! Skip family
                    read(1) ! Skip tags
                endif
                read(1)age
                read(1) ! Skip metallicity
                read(1) ! Skip imass
            endif
            close(1)

            ! Get variable info for particles in the 
            ! region of interest
            partloop: do i=1,npart2
                distance = 0D0
                part%x = x(i,:)
                part%birth_time = age(i)
                ! Check if particle is inside the desired region
                part%x = part%x - reg%centre
                x(i,:) = part%x
                call checkifinside(x(i,:),reg,ok_part,distance)
                ok_filter = filter_particle(reg,filt,part)
                ok_part = ok_part.and.ok_filter
                call getparttype(part,ptype)
                partprob = rand()
                if (ok_part.and.(ptype.eq.'dm').and.(partprob.le.prob)) then
                    npartsaved = npartsaved + 1
                    ! Position to kpc
                    part%x = part%x * (sim%unit_l*cm2kpc)
                    write(7,100)part%x%x,part%x%y,part%x%z
                    100 format(3F16.7)
                endif
            end do partloop
            inpart = inpart + npart2
            deallocate(m,x)
            if(nstar>0)deallocate(age)
        end do cpuloop
        
        close(7)
        write(*,102)npartsaved
        102 format('File includes ',I12,' DM particles')
    end subroutine part2disperse

end module export_part