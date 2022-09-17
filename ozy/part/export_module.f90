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

    subroutine part2skirt(repository,reg,filt,h,smoothmethod,sedmethod,outpath)
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

        ! Specific variables for this subroutine
        logical :: ok_part,ok_filter
        integer :: roterr
        integer :: i,j,k
        integer :: ipos,icpu,binpos
        integer :: npart,npart2,nstar,inpart=0
        integer :: ncpu2,ndim2
        integer :: nstarsaved
        real(dbl) :: distance,tempage,dx
        character(5) :: nchar,ncharcpu
        character(6) :: ptype
        character(128) :: nomfich,varname
        type(vector) :: xtemp,vtemp
        type(particle) :: part
#ifndef IMASS
        integer(1),dimension(:), allocatable :: part_tags
#endif
        integer,dimension(:),allocatable :: plevel
        real(dbl),dimension(:),allocatable :: m,age,met,imass
        real(dbl),dimension(:,:),allocatable :: x,v

        ! Obtain details of the hydro variables stored
        call read_hydrofile_descriptor(repository)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = amr%nlevelmax

#ifndef IMASS
        if (sim%eta_sn .eq. -1D0) then
            write(*,*)': eta_sn=-1 and not IMASS --> should set this up!'
            stop
        end if
#endif
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
                    '# Column 4: smoothing length (kpc)',/,&
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
            allocate(m(1:npart2))
            if(nstar>0)then
                allocate(age(1:npart2))
                allocate(met(1:npart2))
                allocate(imass(1:npart2))
#ifndef IMASS
                allocate(part_tags(1:npart2))
#endif
            ! Settup arrays for the required smoothmethod
            if (TRIM(sedmethod).eq.'level') then
                allocate(plevel(1:npart2))
            end if
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
                read(1) ! Skip id
                if (allocated(plevel)) then
                    read(1) plevel
                else
                    read(1) ! Skip level
                end if
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
            endif
            close(1)

            ! Get variable info for particles in the 
            ! region of interest
            partloop: do i=1,npart2
                distance = 0D0
                part%x = x(i,:)
                part%v = v(i,:)
                part%m = m(i)
                if (nstar>0) then
                    part%id = 0 ! We do not care about ids here
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
                else
                    part%id = 0
                    part%age = 0D0
                    part%met = 0D0
                    part%imass = 0D0
                endif
                ! Check if particle is inside the desired region
                part%x = part%x - reg%centre
                x(i,:) = part%x
                call checkifinside(x(i,:),reg,ok_part,distance)
                ok_filter = filter_particle(reg,filt,part)
                ok_part = ok_part.and.ok_filter
                call getparttype(part,ptype)
                if (ok_part.and.(ptype.eq.'star').and.(part%m.gt.0D0)) then
                    nstarsaved = nstarsaved + 1
                    ! Position to kpc
                    part%x = part%x * (sim%unit_l*cm2kpc)
                    ! Velocity with respect to COM of galaxy and in km/s
                    part%v = part%v - reg%bulk_velocity
                    part%v = part%v * (sim%unit_v*cm2km)
                    ! Initial mass in Msun
                    part%imass = part%imass * (sim%unit_m*g2msun)
                    ! Particle age in Gyr
                    varname = 'star/age'
                    call getpartvalue(reg,part,varname,tempage)
                    part%age = tempage
                    write(7,100)part%x%x,part%x%y,part%x%z,&
                                h,&
                                part%v%x,part%v%y,part%v%z,&
                                part%imass,&
                                part%met,&
                                part%age
                    100 format(3F10.6,F10.6,4F10.2,F10.6,F10.6)
                endif
            end do partloop
            deallocate(m,x,v)
            if (nstar>0)deallocate(age,met,imass)
#ifndef IMASS
            if (nstar>0)deallocate(part_tags)
#endif
            inpart = inpart + npart2
        end do cpuloop

        close(7)
        write(*,102)nstarsaved
        102 format('File includes ',I12,' star particles')
    end subroutine part2skirt

end module export_part