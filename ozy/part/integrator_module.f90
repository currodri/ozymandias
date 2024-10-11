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
        real(dbl),dimension(:,:,:),allocatable :: data
        integer(irg) :: ndm,nstar,nids,ncont
#ifndef LONGINT
        integer(irg), dimension(:),allocatable :: ids
        integer(irg), dimension(:),allocatable :: cont_ids
#else
        integer(ilg), dimension(:),allocatable :: ids
        integer(ilg), dimension(:),allocatable :: cont_ids
#endif
    end type part_region_attrs

    contains

    subroutine allocate_part_regions_attrs(attrs)
        implicit none
        type(part_region_attrs),intent(inout) :: attrs

        if (.not.allocated(attrs%varnames)) allocate(attrs%varnames(attrs%nvars))
        if (.not.allocated(attrs%wvarnames)) allocate(attrs%wvarnames(attrs%nwvars))
        if (.not.allocated(attrs%data)) allocate(attrs%data(attrs%nvars,attrs%nwvars,4))
    end subroutine allocate_part_regions_attrs

    subroutine extract_data(reg,part,attrs)
        use vectors
        use geometrical_regions
        implicit none

        ! Input/output variables
        type(region),intent(in) :: reg
        type(particle),intent(in) :: part
        type(part_region_attrs),intent(inout) :: attrs

        ! Local variables
        integer :: i,j,index
        real(dbl) :: ytemp,wtemp
        character(128) :: tempvar,vartype,varname

        varloop: do i=1,attrs%nvars
            ! Get variable
            call getpartvalue(reg,part,attrs%varnames(i),ytemp)
            if (ytemp/=0D0) then
                wvarloop: do j=1,attrs%nwvars
                    tempvar = TRIM(attrs%wvarnames(j))
                    index = scan(tempvar,'/')
                    vartype = tempvar(1:index-1)
                    varname = tempvar(index+1:)
                    wtemp = 0D0
                    ! Get weights
                    if (varname=='counts'.or.varname=='cumulative') then
                        wtemp = 1D0
                    else
                        call getpartvalue(reg,part,attrs%wvarnames(j),wtemp)
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
            endif
        end do varloop

    end subroutine extract_data

    subroutine renormalise(attrs)
        implicit none

        ! Input part_region_attrs type
        type(part_region_attrs),intent(inout) :: attrs

        ! Local variable
        integer :: i,j,index,index2,indexw
        character(128) :: tempvar,vartype,varname,sfrstr
        character(128) :: tempwvar,wvartype,wvarname
        real(dbl) :: sfrind

        varloop: do i=1,attrs%nvars
            tempvar = TRIM(attrs%varnames(i))
            index = scan(tempvar,'/')
            vartype = tempvar(1:index-1)
            varname = tempvar(index+1:)
            index2 = scan(varname,'_')
            wvarloop: do j=1,attrs%nwvars
                tempwvar = TRIM(attrs%wvarnames(j))
                indexw = scan(tempwvar,'/')
                wvartype = tempwvar(1:indexw-1)
                wvarname = tempwvar(indexw+1:)
                if (wvarname .eq. 'cumulative' .and. index2.ne.0 .and. varname(1:index2-1).eq.'sfr') then
                    sfrstr = varname(index2+1:)
                    read(sfrstr,'(F10.0)') sfrind
                    attrs%data(i,j,1) = attrs%data(i,j,1) * sim%unit_m/ (sfrind*1D6*2D33) ! We now have it in Msun/yr
                elseif (wvarname /= 'cumulative' .and. index2.eq.0) then
                    attrs%data(i,j,1) = attrs%data(i,j,1) / attrs%data(i,j,4)
                endif
            end do wvarloop
        end do varloop
    end subroutine renormalise

    subroutine integrate_region(repository,reg,filt,attrs,get_ids,check_contamination)
        use utils, only:quick_sort_dp
        use vectors
        use coordinate_systems
        use geometrical_regions
        implicit none
        
        ! Input/output variables
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        type(filter),intent(in) :: filt
        type(part_region_attrs),intent(inout) :: attrs
        logical,intent(in),optional :: get_ids
        logical,intent(in),optional :: check_contamination

        ! Specific variables for this subroutine
        logical :: ok_part,ok_filter
        integer :: roterr
        integer :: i,j,k
        integer :: ipos,icpu,binpos
        integer :: npart,npart2,nstar,inpart=0
        integer :: ncpu2,ndim2
        real(dbl) :: distance
        real(dbl),dimension(1:3,1:3) :: trans_matrix
        real(dbl),dimension(1:3) :: xpos
        character(5) :: nchar,ncharcpu
        character(6) :: ptype
        character(128) :: nomfich
        type(vector) :: xtemp,vtemp
        type(particle) :: part
        real(dbl),dimension(:,:),allocatable :: part_data_d
        integer(1),dimension(:,:),allocatable :: part_data_b
#ifdef LONGINT
        integer(ilg),dimension(:,:),allocatable :: part_data_i
#else
        integer(irg),dimension(:,:),allocatable :: part_data_i
#endif
        integer :: count=0
        integer :: nmasscont=0
        integer :: npartcont=0
        real(dbl) :: ptmassmin=1d0,ptmassmax=0d0
        real(dbl) :: masshr=0d0,masslr=0d0
        real(dbl),dimension(1:100) :: massresbins=0d0
        integer,dimension(1:100) :: order=0

        
        ! Obtain details of the hydro variables stored
        call read_hydrofile_descriptor(repository)

        ! Check if particle data uses family
        if (sim%dm .and. sim%hydro) call check_families(repository)

        ! Read the format of the particle data stored
        call read_partfile_descriptor(repository)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = amr%nlevelmax

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

        ! If we're asked to get particle IDs, allocate array with maximum number of particles
        ! and initiliase to zero
        if (present(get_ids) .and. get_ids) then
            allocate(attrs%ids(1:npart))
            attrs%ids = 0
            attrs%nids = npart
        endif

        ! If we're asked to check contamination, we also need to save the IDs
        ! of the low resolution particles
        if (present(get_ids) .and. get_ids .and. &
            & present(check_contamination) .and. check_contamination) then
            allocate(attrs%cont_ids(1:npart))
            attrs%cont_ids = 0
            attrs%ncont = 0
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
                part%x = (/part_data_d(get_ipos(partIDs,partIDs%position_x),i),&
                            part_data_d(get_ipos(partIDs,partIDs%position_y),i),&
                            part_data_d(get_ipos(partIDs,partIDs%position_z),i)/)
                part%x = part%x / sim%boxlen
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
#else
                part%imass = 0D0
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
#endif
                ! Check if particle is inside the desired region
                part%x = part%x - reg%centre
                call rotate_vector(part%x,trans_matrix)
                xpos = part%x
                call checkifinside(xpos,reg,ok_part,distance)
                ok_filter = filter_particle(reg,filt,part)
                ok_part = ok_part.and.ok_filter
                if (ok_part) then
                    if (present(get_ids) .and. get_ids) attrs%ids(inpart+i) = part%id
                    call getparttype(part,ptype)
                    if (present(check_contamination) .and. check_contamination &
                        .and. ptype .eq. 'dm') then
                        ptmassmin = min(part%m,ptmassmin)
                        ptmassmax = max(part%m,ptmassmax)
                        if (.not. any(massresbins == part%m)) then
                            massresbins(nmasscont+1) = part%m
                            nmasscont = nmasscont + 1
                        endif
                        if (part%m > ptmassmin) then
                            npartcont = npartcont + 1
                            if (present(get_ids) .and. get_ids) attrs%cont_ids(inpart+i) = part%id
                            masslr = masslr + part%m
                        else
                            masshr = masshr + part%m
                        endif
                    endif
                    if (ptype.eq.'dm') attrs%ndm = attrs%ndm + 1
                    if (ptype.eq.'star') attrs%nstar = attrs%nstar + 1
                    call rotate_vector(part%v,trans_matrix)
                    call extract_data(reg,part,attrs)
                endif
            end do partloop
            deallocate(part_data_d,part_data_i,part_data_b)
            inpart = inpart + npart2
        end do cpuloop

        ! Finally, just renormalise for weighted quantities
        call renormalise(attrs)

        ! If we want to check contamination, print
        ! some global statistics
        if (present(check_contamination) .and. check_contamination) then
            write(*,*)'==== CONTAMINATION REPORT ===='
            write(*,*)'>>> Total Mass:   ',masslr+masshr
            write(*,*)'>>> % Mass in LR: ', masslr/(masslr+masshr)
            write(*,*)'>>> # of LR:      ', npartcont
            write(*,*)'>>> # of mass res:', nmasscont
            if (nmasscont .gt. 0) then
                call quick_sort_dp(massresbins(1:nmasscont),order(1:nmasscont),nmasscont)
                write(*,*)'>>> Mass resolutions >>>'
                do i=1,nmasscont
                    write(*,*) massresbins(i)
                end do
            endif
        endif
    end subroutine integrate_region

end module part_integrator