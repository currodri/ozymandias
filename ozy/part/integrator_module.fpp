# 1 "integrator_module.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "integrator_module.f90"
!--------------------------------------------------------------------------
! ozymandias:integrator_module.f90
!--------------------------------------------------------------------------
!
! MODULE: part_integrator
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
!> @date 31/5/2021   0.1.1 translating from f2py code to f90wrap
!--------------------------------------------------------------------------

module part_integrator
    use local
    use io_ramses
    use filtering
    use cosmology

    type part_region_attrs
        integer :: nvars
        character(128),dimension(:),allocatable :: varnames
        integer :: nwvars
        character(128),dimension(:),allocatable :: wvarnames
        integer :: ndm,nstar
        real(dbl),dimension(:,:,:),allocatable :: data
    end type part_region_attrs

    contains

    subroutine allocate_part_regions_attrs(attrs)
        implicit none
        type(part_region_attrs),intent(inout) :: attrs

        if (.not.allocated(attrs%varnames)) allocate(attrs%varnames(attrs%nvars))
        if (.not.allocated(attrs%wvarnames)) allocate(attrs%wvarnames(attrs%nwvars))
        if (.not.allocated(attrs%data)) allocate(attrs%data(attrs%nvars,attrs%nwvars,4))
    end subroutine allocate_part_regions_attrs

    subroutine extract_data(sim,reg,part,attrs)
        use vectors
        use geometrical_regions
        implicit none

        ! Input/output variables
        type(sim_info),intent(in) :: sim
        type(region),intent(in) :: reg
        type(particle),intent(in) :: part
        type(part_region_attrs),intent(inout) :: attrs

        ! Local variables
        integer :: i,j,index
        real(dbl) :: ytemp,wtemp
        character(128) :: tempvar,vartype,varname,true_wname

        varloop: do i=1,attrs%nvars
            ! Get variable
            call getpartvalue(sim,reg,part,attrs%varnames(i),ytemp)
            wvarloop: do j=1,attrs%nwvars
                ! Get weights
                if (attrs%wvarnames(j)=='counts'.or.attrs%wvarnames(j)=='cumulative') then
                    if (ytemp .eq. 0D0) then
                        wtemp = 0D0
                    else
                        wtemp =  1D0
                    endif
                else
                    if (ytemp .eq. 0D0) then
                        wtemp = 0D0
                    else
                        index = scan(attrs%varnames(i),'/')
                        true_wname = trim(attrs%varnames(i)(1:index-1))//'/'//attrs%wvarnames(j)
                        call getpartvalue(sim,reg,part,true_wname,wtemp)
                    endif
                endif
                ! Save to attrs
                attrs%data(i,j,1) = attrs%data(i,j,1) + ytemp*wtemp ! Value (weighted or not)
                if (attrs%data(i,j,2).eq.0D0) then
                    attrs%data(i,j,2) = ytemp ! Just to make sure that the initial min is not zero
                else
                    attrs%data(i,j,2) = min(ytemp,attrs%data(i,j,2))    ! Min value
                endif
                attrs%data(i,j,3) = max(ytemp,attrs%data(i,j,3))    ! Max value
                attrs%data(i,j,4) = attrs%data(i,j,4) + wtemp       ! Weight

            end do wvarloop
        end do varloop

    end subroutine extract_data

    subroutine renormalise(sim,attrs)
        implicit none

        ! Input part_region_attrs type
        type(sim_info),intent(in) :: sim
        type(part_region_attrs),intent(inout) :: attrs

        ! Local variable
        integer :: i,j,index,index2
        character(128) :: tempvar,vartype,varname,sfrstr
        real(dbl) :: sfrind

        varloop: do i=1,attrs%nvars
            tempvar = TRIM(attrs%varnames(i))
            index = scan(tempvar,'/')
            vartype = tempvar(1:index-1)
            varname = tempvar(index+1:)
            index2 = scan(varname,'_')
            wvarloop: do j=1,attrs%nwvars
                if (attrs%wvarnames(j) .eq. 'cumulative' .and. index2.ne.0) then
                    sfrstr = varname(index2+1:)
                    read(sfrstr,'(F10.0)') sfrind
                    attrs%data(i,j,1) = attrs%data(i,j,1) * sim%unit_m/ (sfrind*1D6*2D33) ! We now have it in Msun/yr
                elseif (attrs%wvarnames(j) /= 'cumulative' .and. index2.eq.0) then
                    attrs%data(i,j,1) = attrs%data(i,j,1) / attrs%data(i,j,4)
                endif
            end do wvarloop
        end do varloop
    end subroutine renormalise

    subroutine integrate_region(repository,reg,filt,attrs)
        use vectors
        use coordinate_systems
        use geometrical_regions
        implicit none
        
        ! Input/output variables
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        type(filter),intent(in) :: filt
        type(part_region_attrs),intent(inout) :: attrs

        ! Ozymandias derived types for RAMSES
        type(hydroID) :: varIDs
        type(amr_info) :: amr
        type(sim_info) :: sim

        ! Specific variables for this subroutine
        logical :: ok_part,ok_filter
        integer :: roterr
        integer :: i,j,k
        integer :: ipos,icpu,binpos
        integer :: npart,npart2,nstar,ncpu2,ndim2
        real(dbl) :: distance
        real(dbl),dimension(1:3,1:3) :: trans_matrix
        character(5) :: nchar,ncharcpu
        character(6) :: ptype
        character(128) :: nomfich
        type(vector) :: xtemp,vtemp
        type(particle) :: part
        integer,dimension(:),allocatable :: id
        real(dbl),dimension(:),allocatable :: m,age,met,imass
        real(dbl),dimension(:,:),allocatable :: x,v

        integer :: count
        count = 0

        ! Obtain details of the hydro variables stored
        call read_hydrofile_descriptor(repository,varIDs)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository,amr,sim)
        amr%lmax = amr%nlevelmax

        ! Check if particle data uses family
        call check_families(repository,sim)

        ! Compute the Hilbert curve
        call get_cpu_map(reg,amr)
        write(*,*)'ncpu_read:',amr%ncpu_read

        ! Just make sure that initial values are zero
        attrs%data = 0D0
        attrs%nstar = 0
        attrs%ndm = 0

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
            ! write(*,*)'Processing file '//TRIM(nomfich)
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
                if (sim%family) then
                    read(1) ! Skip family
                    read(1) ! Skip tags
                endif
                read(1)age
                read(1)met
                read(1)imass
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
                part%x = part%x - reg%centre
                call rotate_vector(part%x,trans_matrix)
                x(i,:) = part%x
                call checkifinside(x(i,:),reg,ok_part,distance)
                ok_filter = filter_particle(sim,reg,filt,part)
                ok_part = ok_part.and.ok_filter
                if (ok_part) then
                    call getparttype(part,ptype)
                    if (ptype.eq.'dm') attrs%ndm = attrs%ndm + 1
                    if (ptype.eq.'star') attrs%nstar = attrs%nstar + 1
                    call rotate_vector(part%v,trans_matrix)
                    call extract_data(sim,reg,part,attrs)
                endif
            end do partloop
            deallocate(m,x,v)
            if (nstar>0)deallocate(id,age,met,imass)
        end do cpuloop

        ! Finally, just renormalise fo weighted quantities
        call renormalise(sim,attrs)

    end subroutine integrate_region

end module part_integrator
