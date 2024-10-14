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
        integer :: nvars=1
        character(128),dimension(:),allocatable :: varnames
        integer :: nfilter=1,nsubs=0
        type(filter),dimension(:),allocatable :: filters
        type(region),dimension(:),allocatable :: subs
        type(pdf_handler) :: result
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
        if (.not.allocated(attrs%filters)) allocate(attrs%filters(attrs%nfilter))
        if (.not.allocated(attrs%subs).and.attrs%nsubs>0) allocate(attrs%subs(attrs%nsubs))
    end subroutine allocate_part_regions_attrs

    subroutine extract_data(reg,part,attrs,ifilt,trans_matrix)
        use vectors
        use geometrical_regions
        implicit none

        ! Input/output variables
        type(region),intent(in) :: reg
        type(particle),intent(in) :: part
        integer, intent(in) :: ifilt
        type(part_region_attrs),intent(inout) :: attrs
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix

        ! Local variables
        integer :: i,j,index,ibin
        real(dbl) :: ytemp,wtemp,ytemp2
        character(128) :: tempvar,vartype,varname

        varloop: do i=1,attrs%result%nvars
            if (attrs%result%do_binning(i)) then
                ! Get variable
                call findbinpos_part(reg,part,ibin,ytemp,trans_matrix,attrs%result%scaletype(i),&
                                    & attrs%result%nbins,attrs%result%bins(:,i),&
                                    & attrs%result%linthresh(i),attrs%result%zero_index(i),&
                                    & attrs%result%varname(i))
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
                        tempvar = TRIM(attrs%result%wvarnames(j))
                        index = scan(tempvar,'/')
                        vartype = tempvar(1:index-1)
                        varname = tempvar(index+1:)
                        if (trim(varname)=='counts') then
                            wtemp =  1D0
                        else if (trim(varname)=='cumulative') then
                            wtemp = ytemp2
                        else
                            call getpartvalue(reg,part,attrs%result%wvarnames(j),wtemp)
                        endif
                        
                        ! Save to PDFs
                        attrs%result%heights(i,ifilt,j,ibin) = attrs%result%heights(i,ifilt,j,ibin) + wtemp ! Weight to the PDF bin
                        attrs%result%totweights(i,ifilt,j) = attrs%result%totweights(i,ifilt,j) + wtemp       ! Weight

                        ! Now do it for the case of no binning (old integration method)
                        ! Get weights
                        if (trim(varname)=='counts') then
                            wtemp =  1D0
                            ytemp2 = 1D0
                        else if (trim(varname)=='cumulative') then
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
                call getpartvalue(reg,part,attrs%result%varname(i),ytemp)
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
                    tempvar = TRIM(attrs%result%wvarnames(j))
                    index = scan(tempvar,'/')
                    vartype = tempvar(1:index-1)
                    varname = tempvar(index+1:)
                    if (trim(varname)=='counts') then
                        wtemp =  1D0
                        ytemp2 = 1D0
                    else if (trim(varname)=='cumulative') then
                        wtemp = 1D0
                    else
                        call getpartvalue(reg,part,attrs%result%wvarnames(j),wtemp)
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
        integer :: i,j,index,index2,indexw,ifilt
        character(128) :: tempvar,vartype,varname,sfrstr
        character(128) :: tempwvar,wvartype,wvarname
        real(dbl) :: sfrind

        filterloop: do ifilt=1,attrs%nfilter
            varloop: do i=1,attrs%result%nvars
                tempvar = TRIM(attrs%result%varname(i))
                index = scan(tempvar,'/')
                vartype = tempvar(1:index-1)
                varname = tempvar(index+1:)
                index2 = scan(varname,'_')
                if (attrs%result%do_binning(i)) then
                    wvarloop1: do j=1,attrs%result%nwvars
                        tempwvar = TRIM(attrs%result%wvarnames(j))
                        indexw = scan(tempwvar,'/')
                        wvartype = tempwvar(1:indexw-1)
                        wvarname = tempwvar(indexw+1:)
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
                        tempwvar = TRIM(attrs%result%wvarnames(j))
                        indexw = scan(tempwvar,'/')
                        wvartype = tempwvar(1:indexw-1)
                        wvarname = tempwvar(indexw+1:)
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

    subroutine integrate_region(repository,reg,attrs,get_ids,check_contamination)
        use utils, only:quick_sort_dp
        use vectors
        use coordinate_systems
        use geometrical_regions
        implicit none
        
        ! Input/output variables
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        type(part_region_attrs),intent(inout) :: attrs
        logical,intent(in),optional :: get_ids
        logical,intent(in),optional :: check_contamination

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
        type(vector) :: xtemp,vtemp
        type(particle) :: part
#ifndef LONGINT
        integer(irg),dimension(:),allocatable :: id
#else
        integer(ilg),dimension(:),allocatable :: id
#endif
#ifndef IMASS
        integer(1),dimension(:), allocatable :: part_tags
#endif
        real(dbl),dimension(:),allocatable :: m,age,met,imass
        real(dbl),dimension(:,:),allocatable :: x,v

        integer :: count=0
        integer :: nmasscont=0
        integer :: npartcont=0
        real(dbl) :: ptmassmin=1d0,ptmassmax=0d0
        real(dbl) :: masshr=0d0,masslr=0d0
        real(dbl),dimension(1:100) :: massresbins=0d0
        integer,dimension(1:100) :: order=0

#ifdef LONGINT
        write(*,*) 'Using LONGINT for particle IDs'
#endif
#ifndef IMASS
        if (sim%eta_sn .eq. -1D0) then
            write(*,*)': eta_sn=-1 and not IMASS --> should set this up!'
            stop
        end if
#endif
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

        ! Just make sure that initial values are zero
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
        if (present(get_ids)) then
            if (get_ids) then
                allocate(attrs%ids(1:npart))
                attrs%ids = 0
                attrs%nids = npart
            end if
        endif

        ! If we're asked to check contamination, we also need to save the IDs
        ! of the low resolution particles
        if (present(get_ids) .and. present(check_contamination)) then
            if (get_ids .and. check_contamination) then
                allocate(attrs%cont_ids(1:npart))
                attrs%cont_ids = 0
                attrs%ncont = 0
            end if
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
            if (present(get_ids)) then
                if (get_ids .and. (.not. allocated(id))) allocate(id(1:npart2))
            end if
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
            elseif (present(get_ids)) then
                if(get_ids .and. nstar .eq. 0) read(1)id
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

                elseif (present(get_ids)) then
                    if (get_ids) then
                        part%id = id(i)
                        part%age = 0D0
                        part%met = 0D0
                        part%imass = 0D0
                    end if
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

                ! If we are avoiding substructure, check whether we are safe
                if (attrs%nsubs>0) then
                    ok_sub = .true.
                    do isub=1,attrs%nsubs
                        ok_sub = ok_sub .and. filter_sub(attrs%subs(isub),x(i,:))
                    end do
                    ok_part = ok_part .and. ok_sub
                end if
                if (ok_part) then
                    if (present(get_ids)) then
                        if (get_ids) attrs%ids(inpart+i) = part%id
                    end if
                    call getparttype(part,ptype)
                    if (present(check_contamination)) then
                        if (check_contamination .and. ptype .eq. 'dm') then
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
                        end if
                    endif
                    if (ptype.eq.'dm') attrs%ndm = attrs%ndm + 1
                    if (ptype.eq.'star') attrs%nstar = attrs%nstar + 1
                    call rotate_vector(part%v,trans_matrix)
                    filterloop: do ifilt=1,attrs%nfilter
                        ok_filter = filter_particle(reg,attrs%filters(ifilt),part)
                        if (ok_filter) call extract_data(reg,part,attrs,ifilt,trans_matrix)
                    end do filterloop  
                endif
            end do partloop
            deallocate(m,x,v)
            if (allocated(id))deallocate(id)
            if (nstar>0)deallocate(age,met,imass)
#ifndef IMASS
            if (nstar>0)deallocate(part_tags)
#endif
            inpart = inpart + npart2
        end do cpuloop

        ! Finally, just renormalise for weighted quantities
        call renormalise(attrs)
        ! If we want to check contamination, print
        ! some global statistics
        if (present(check_contamination)) then
            if (check_contamination) then
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
            end if
        endif
    end subroutine integrate_region

end module part_integrator