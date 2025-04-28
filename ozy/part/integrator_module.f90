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
    use geometrical_regions
    use io_ramses
    use filtering
    use cosmology
    use stats_utils

    type part_region_attrs
        integer :: nvars=1,nwvars=1
        character(128),dimension(:),allocatable :: varnames
        character(128),dimension(:),allocatable :: wvarnames
        integer :: nfilter=1,nsubs=0
        type(filter_part),dimension(:),allocatable :: filters
        type(region),dimension(:),allocatable :: subs
        type(pdf_handler) :: result
        type(part_var),dimension(:),allocatable :: vars
        type(part_var),dimension(:),allocatable :: wvars
    end type part_region_attrs

    contains

    subroutine allocate_part_regions_attrs(attrs)
        implicit none
        type(part_region_attrs),intent(inout) :: attrs

        if (.not.allocated(attrs%varnames)) allocate(attrs%varnames(attrs%nvars))
        if (.not.allocated(attrs%wvarnames)) allocate(attrs%wvarnames(attrs%nwvars))
        if (.not.allocated(attrs%filters)) allocate(attrs%filters(attrs%nfilter))
        if (.not.allocated(attrs%subs).and.attrs%nsubs>0) allocate(attrs%subs(attrs%nsubs))
        if (.not.allocated(attrs%vars)) allocate(attrs%vars(attrs%nvars))
        if (.not.allocated(attrs%wvars)) allocate(attrs%wvars(attrs%nwvars))
    end subroutine allocate_part_regions_attrs

    subroutine extract_data(reg,part_data_d,part_data_i,part_data_b,attrs,ifilt,trans_matrix)
        use vectors
        use geometrical_regions
        implicit none

        ! Input/output variables
        type(region),intent(in) :: reg
        real(dbl),dimension(1:sim%nvar_part_d),intent(in) :: part_data_d
#ifdef LONGINT
        integer(ilg),dimension(1:sim%nvar_part_i),intent(in) :: part_data_i
#else
        integer(irg),dimension(1:sim%nvar_part_i),intent(in) :: part_data_i
#endif
        integer(1),dimension(1:sim%nvar_part_b),intent(in) :: part_data_b
        integer, intent(in) :: ifilt
        type(part_region_attrs),intent(inout) :: attrs
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix

        ! Local variables
        integer :: i,j,index,ibin
        real(dbl) :: ytemp,wtemp,ytemp2
#ifdef LONGINT
        integer(ilg) :: wtemp_i,ytemp_i
#else
        integer(irg) :: wtemp_i,ytemp_i
#endif
        integer(1) :: wtemp_b,ytemp_b
        type(vector) :: dcell

        ! TODO: dcell is not used. Think about the meaning of dcell in an integration

        varloop: do i=1,attrs%result%nvars
            if (attrs%result%do_binning(i)) then
                ! Get variable
                call findbinpos_part(reg,dcell,part_data_d,part_data_i,part_data_b,&
                                    & ibin,ytemp,trans_matrix,attrs%result%scaletype(i),&
                                    & attrs%result%nbins,attrs%result%bins(:,i),&
                                    & attrs%result%linthresh(i),attrs%result%zero_index(i),&
                                    & attrs%vars(i))
                if (ytemp.eq.0d0) cycle
                ! Get min and max
                if (attrs%result%minv(i,ifilt).eq.0D0) then
                    attrs%result%minv(i,ifilt) = ytemp ! Just to make sure that the initial min is not zero
                else
                    attrs%result%minv(i,ifilt) = min(ytemp,attrs%result%minv(i,ifilt))    ! Min value
                endif
                attrs%result%maxv(i,ifilt) = max(ytemp,attrs%result%maxv(i,ifilt))    ! Max value
                attrs%result%nvalues(i,ifilt) = attrs%result%nvalues(i,ifilt) + 1
                
                if (ibin.gt.0) then
                    wvarloop1: do j=1,attrs%result%nwvars
                        ! Get weights
                        ytemp2 = ytemp
                        if (TRIM(attrs%result%wvarnames(j))=='counts') then
                            wtemp =  1D0
                        else if (TRIM(attrs%result%wvarnames(j))=='cumulative') then
                            wtemp = ytemp2
                        else
                            if (attrs%wvars(j)%vartype == 1) then
                                wtemp = attrs%wvars(j)%myfunction_d(amr,sim,attrs%wvars(j),reg,dcell,&
                                                                & part_data_d,part_data_i,part_data_b)
                            else if (attrs%wvars(j)%vartype == 2) then
                                wtemp_i = attrs%wvars(j)%myfunction_i(amr,sim,attrs%wvars(j),reg,dcell,&
                                                                & part_data_d,part_data_i,part_data_b)
                                wtemp = real(wtemp_i,kind=dbl)
                            else if (attrs%wvars(j)%vartype == 3) then
                                wtemp_b = attrs%wvars(j)%myfunction_b(amr,sim,attrs%wvars(j),reg,dcell,&
                                                                & part_data_d,part_data_i,part_data_b)
                                wtemp = real(wtemp_b,kind=dbl)
                            else
                                write(*,*)'Error: unknown variable type in extract_data'
                                stop
                            end if
                        endif
                        
                        ! Save to PDFs
                        attrs%result%heights(i,ifilt,j,ibin) = attrs%result%heights(i,ifilt,j,ibin) + wtemp ! Weight to the PDF bin
                        attrs%result%totweights(i,ifilt,j) = attrs%result%totweights(i,ifilt,j) + wtemp       ! Weight

                        ! Now do it for the case of no binning (old integration method)
                        ! Get weights
                        if (TRIM(attrs%result%wvarnames(j))=='counts') then
                            wtemp =  1D0
                            ytemp2 = 1D0
                        else if (TRIM(attrs%result%wvarnames(j))=='cumulative') then
                            wtemp = 1D0
                        endif
                        
                        ! Save to attrs
                        attrs%result%total(i,ifilt,j,1) = attrs%result%total(i,ifilt,j,1) + ytemp2*wtemp ! Value (weighted or not)
                        attrs%result%total(i,ifilt,j,2) = attrs%result%total(i,ifilt,j,2) + wtemp       ! Weight
                    end do wvarloop1
                else
                    attrs%result%nout(i,ifilt) = attrs%result%nout(i,ifilt) + 1
                end if
            else
                ! Get variable
                if (attrs%vars(i)%vartype == 1) then
                    ytemp = attrs%vars(i)%myfunction_d(amr,sim,attrs%vars(i),reg,dcell,&
                                                        & part_data_d,part_data_i,part_data_b)
                else if (attrs%vars(i)%vartype == 2) then
                    ytemp_i = attrs%vars(i)%myfunction_i(amr,sim,attrs%vars(i),reg,dcell,&
                                                        & part_data_d,part_data_i,part_data_b)
                    ytemp = real(ytemp_i,kind=dbl)
                else if (attrs%vars(i)%vartype == 3) then
                    ytemp_b = attrs%vars(i)%myfunction_b(amr,sim,attrs%vars(i),reg,dcell,&
                                                        & part_data_d,part_data_i,part_data_b)
                    ytemp = real(ytemp_b,kind=dbl)
                else
                    write(*,*)'Error: unknown variable type in extract_data'
                    stop
                end if
                ! Get min and max
                if (attrs%result%minv(i,ifilt).eq.0D0) then
                    attrs%result%minv(i,ifilt) = ytemp ! Just to make sure that the initial min is not zero
                else
                    attrs%result%minv(i,ifilt) = min(ytemp,attrs%result%minv(i,ifilt))    ! Min value
                endif
                attrs%result%maxv(i,ifilt) = max(ytemp,attrs%result%maxv(i,ifilt))    ! Max value
                attrs%result%nvalues(i,ifilt) = attrs%result%nvalues(i,ifilt) + 1

                wvarloop2: do j=1,attrs%result%nwvars
                    ! Get weights
                    ytemp2 = ytemp
                    if (TRIM(attrs%result%wvarnames(j))=='counts') then
                        wtemp =  1D0
                        ytemp2 = 1D0
                    else if (TRIM(attrs%result%wvarnames(j))=='cumulative') then
                        wtemp = 1D0
                    else
                        if (attrs%wvars(j)%vartype == 1) then
                            wtemp = attrs%wvars(j)%myfunction_d(amr,sim,attrs%wvars(j),reg,dcell,&
                                                                & part_data_d,part_data_i,part_data_b)
                        else if (attrs%wvars(j)%vartype == 2) then
                            wtemp_i = attrs%wvars(j)%myfunction_i(amr,sim,attrs%wvars(j),reg,dcell,&
                                                                & part_data_d,part_data_i,part_data_b)
                            wtemp = real(wtemp_i,kind=dbl)
                        else if (attrs%wvars(j)%vartype == 3) then
                            wtemp_b = attrs%wvars(j)%myfunction_b(amr,sim,attrs%wvars(j),reg,dcell,&
                                                                & part_data_d,part_data_i,part_data_b)
                            wtemp = real(wtemp_b,kind=dbl)
                        else
                            write(*,*)'Error: unknown variable type in extract_data'
                            stop
                        end if
                    endif
                    
                    ! Save to attrs
                    attrs%result%total(i,ifilt,j,1) = attrs%result%total(i,ifilt,j,1) + ytemp2*wtemp ! Value (weighted or not)
                    attrs%result%total(i,ifilt,j,2) = attrs%result%total(i,ifilt,j,2) + wtemp       ! Weight
                end do wvarloop2
            end if
        end do varloop

    end subroutine extract_data

    subroutine renormalise(attrs)
        implicit none

        ! Input part_region_attrs type
        type(part_region_attrs),intent(inout) :: attrs

        ! Local variable
        integer :: i,j,index2,ifilt
        character(128) :: varname,sfrstr,wvarname
        real(dbl) :: sfrind

        filterloop: do ifilt=1,attrs%nfilter
            varloop: do i=1,attrs%result%nvars
                varname = TRIM(attrs%result%varname(i))
                index2 = scan(varname,'_')
                if (attrs%result%do_binning(i)) then
                    wvarloop1: do j=1,attrs%result%nwvars
                        wvarname = TRIM(attrs%result%wvarnames(j))
                        if (trim(wvarname) .eq. 'cumulative' .and. index2.ne.0 .and. trim(varname(1:index2-1)).eq.'sfr') then
                            sfrstr = varname(index2+1:)
                            read(sfrstr,'(F10.0)') sfrind
                            attrs%result%heights(i,ifilt,j,:) = attrs%result%heights(i,ifilt,j,:) / attrs%result%totweights(i,ifilt,j) * sim%unit_m/ (sfrind*1D6*2D33) ! We now have it in Msun/yr
                            attrs%result%total(i,ifilt,j,1) = attrs%result%total(i,ifilt,j,1) / attrs%result%total(i,ifilt,j,2) * sim%unit_m/ (sfrind*1D6*2D33) ! We now have it in Msun/yr
                        elseif (trim(wvarname) /= 'cumulative' .and. index2.eq.0) then
                            attrs%result%heights(i,ifilt,j,:) = attrs%result%heights(i,ifilt,j,:) / attrs%result%totweights(i,ifilt,j)
                            attrs%result%total(i,ifilt,j,1) = attrs%result%total(i,ifilt,j,1) / attrs%result%total(i,ifilt,j,2)
                        endif
                    end do wvarloop1
                else
                    wvarloop2: do j=1,attrs%result%nwvars
                        wvarname = TRIM(attrs%result%wvarnames(j))
                        if (trim(wvarname) .eq. 'cumulative' .and. index2.ne.0 .and. trim(varname(1:index2-1)).eq.'sfr') then
                            sfrstr = varname(index2+1:)
                            read(sfrstr,'(F10.0)') sfrind
                            attrs%result%total(i,ifilt,j,1) = attrs%result%total(i,ifilt,j,1) / attrs%result%total(i,ifilt,j,2) * sim%unit_m/ (sfrind*1D6*2D33) ! We now have it in Msun/yr
                        elseif (trim(wvarname) /= 'cumulative' .and. index2.eq.0) then
                            attrs%result%total(i,ifilt,j,1) = attrs%result%total(i,ifilt,j,1) / attrs%result%total(i,ifilt,j,2)
                        endif
                    end do wvarloop2
                end if
            end do varloop
        end do filterloop
    end subroutine renormalise

    subroutine integrate_region(repository,reg,attrs,lmax,&
                                &part_dict,part_vtypes)
        use utils, only:quick_sort_dp
        use vectors
        use coordinate_systems
        use geometrical_regions
        implicit none
        
        ! Input/output variables
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        type(part_region_attrs),intent(inout) :: attrs
        integer,intent(in) :: lmax
        type(dictf90),intent(in),optional :: part_dict,part_vtypes

        ! Specific variables for this subroutine
        logical :: ok_part,ok_filter,ok_sub
        integer :: roterr
        integer :: i,j,k
        integer :: ipos,icpu,binpos,ifilt,isub
        integer :: npart,npart2,nstar,inpart=0
        integer :: ncpu2,ndim2
        real(dbl) :: distance
        real(dbl),dimension(1:3,1:3) :: trans_matrix
        character(5) :: nchar,ncharcpu
        character(6) :: ptype
        character(128) :: nomfich
        type(vector) :: xtemp,vtemp,dcell
#ifndef LONGINT
        integer(irg),dimension(:,:),allocatable :: part_data_i
#else
        integer(ilg),dimension(:,:),allocatable :: part_data_i
#endif
        real(dbl),dimension(:,:),allocatable :: part_data_d
        integer(1),dimension(:,:),allocatable :: part_data_b
        real(dbl),dimension(:,:),allocatable :: x,v

        integer :: count=0
        integer :: nmasscont=0
        integer :: npartcont=0
        real(dbl) :: ptmassmin=1d0,ptmassmax=0d0
        real(dbl) :: masshr=0d0,masslr=0d0
        real(dbl),dimension(1:100) :: massresbins=0d0
        integer,dimension(1:100) :: order=0

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

        if (attrs%nsubs>0 .and. verbose) write(*,*)'Exluding substructures: ',attrs%nsubs
        if (verbose) write(*,*)'ncpu_read:',amr%ncpu_read

        ! Set up the part variables quicklook tools
        if (present(part_dict).and.present(part_vtypes)) then
            ! If the user provides their own particle dictionary,
            ! use that one instead of the automatic from the
            ! particle_file_descriptor.txt (RAMSES)
            call get_partvar_tools(part_dict,part_vtypes,attrs%nvars,&
                                    &attrs%varnames,attrs%vars)

            ! Do the same for the weight variables
            call get_partvar_tools(part_dict,part_vtypes,attrs%nwvars,&
                                    &attrs%wvarnames,attrs%wvars)
            
            ! Do it for the filters
            do ii=1,attrs%nfilter
                call get_filter_part_tools(part_dict,part_vtypes,attrs%filters(ii))
            end do

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
            call get_partvar_tools(partIDs,partvar_types,attrs%nvars,&
                                    &attrs%varnames,attrs%vars)
            
            ! Now for the weight variables
            call get_partvar_tools(partIDs,partvar_types,attrs%nwvars,&
                                    &attrs%wvarnames,attrs%wvars)
            ! And for the filters
            do ii=1,attrs%nfilter
                call get_filter_part_tools(partIDs,partvar_types,attrs%filters(ii))
            end do
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

                ! If we are avoiding substructure, check whether we are safe
                if (attrs%nsubs>0) then
                    ok_sub = .true.
                    do isub=1,attrs%nsubs
                        ok_sub = ok_sub .and. filter_sub(attrs%subs(isub),x(i,:))
                    end do
                    ok_part = ok_part .and. ok_sub
                end if
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

                    filterloop: do ifilt=1,attrs%nfilter
                        ok_filter = filter_particle(reg,attrs%filters(ifilt),dcell,&
                                                &part_data_d(:,i),part_data_i(:,i),part_data_b(:,i))
                        if (ok_filter) call extract_data(reg,part_data_d(:,i),part_data_i(:,i),&
                                                        &part_data_b(:,i),attrs,ifilt,trans_matrix)
                    end do filterloop  
                endif
            end do partloop
            deallocate(x,v,part_data_d,part_data_i,part_data_b)
            inpart = inpart + npart2
        end do cpuloop

        ! Finally, just renormalise for weighted quantities
        call renormalise(attrs)
    end subroutine integrate_region

end module part_integrator