!--------------------------------------------------------------------------
! ozymandias:integrator_module.f90
!--------------------------------------------------------------------------
!
! MODULE: amr_integrator
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
!> @date 27/5/2021   0.1.1 translating from f2py code to f90wrap
!--------------------------------------------------------------------------

module amr_integrator
    use local
    use io_ramses
    use filtering

    type amr_region_attrs
        integer :: nvars
        character(128),dimension(:),allocatable :: varnames
        integer :: nwvars
        character(128),dimension(:),allocatable :: wvarnames
        real(dbl),dimension(:,:,:),allocatable :: data
    end type amr_region_attrs

    contains

    subroutine allocate_amr_regions_attrs(attrs)
        implicit none
        type(amr_region_attrs),intent(inout) :: attrs

        if (.not.allocated(attrs%varnames)) allocate(attrs%varnames(attrs%nvars))
        if (.not.allocated(attrs%wvarnames)) allocate(attrs%wvarnames(attrs%nwvars))
        if (.not.allocated(attrs%data)) allocate(attrs%data(attrs%nvars,attrs%nwvars,4))
    end subroutine allocate_amr_regions_attrs

    subroutine extract_data(reg,varIDs,pos,cellvars,cellsize,attrs)
        use vectors
        use geometrical_regions
        implicit none

        ! Input/output variables
        type(region),intent(in) :: reg
        type(hydroID),intent(in) :: varIDs
        real(dbl),dimension(1:3),intent(in) :: pos
        real(dbl),dimension(1:varIDs%nvar),intent(in) :: cellvars
        real(dbl),intent(in) :: cellsize
        type(amr_region_attrs),intent(inout) :: attrs

        ! Local variables
        integer :: i,j
        real(dbl) :: ytemp,wtemp
        type(vector) :: x

        x = pos
        varloop: do i=1,attrs%nvars
            ! Get variable
            call getvarvalue(varIDs,reg,cellsize,x,cellvars,attrs%varnames(i),ytemp)
            wvarloop: do j=1,attrs%nwvars
                ! Get weights
                if (attrs%wvarnames(j)=='counts'.or.attrs%wvarnames(j)=='cumulative') then
                    wtemp =  1D0
                else
                    call getvarvalue(varIDs,reg,cellsize,x,cellvars,attrs%wvarnames(j),wtemp)
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

    subroutine renormalise(attrs)
        implicit none

        ! Input/output variables
        type(amr_region_attrs),intent(inout) :: attrs

        ! Local variables
        integer :: i,j

        varloop: do i=1,attrs%nvars
            wvarloop: do j=1,attrs%nwvars
                if (attrs%wvarnames(j) /= 'cumulative') then
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
        type(amr_region_attrs),intent(inout) :: attrs

        ! Ozymandias derived types for RAMSES
        type(hydroID) :: varIDs
        type(amr_info) :: amr
        type(sim_info) :: sim

        ! Specific variables for this subroutine
        integer :: i,j,k
        integer :: ipos,icpu,ilevel,ind,idim,ivar
        integer :: ix,iy,iz,ngrida,nx_full,ny_full,nz_full
        integer :: nvarh
        integer :: roterr
        character(5) :: nchar,ncharcpu
        character(128) :: nomfich
        real(dbl) :: distance,dx
        type(vector) :: xtemp,vtemp
        logical :: ok_cell,ok_filter
        integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
        real(dbl),dimension(1:8,1:3) :: xc
        real(dbl),dimension(3,3) :: trans_matrix
        real(dbl),dimension(:,:),allocatable :: xg,x
        real(dbl),dimension(:,:,:),allocatable :: var
        integer,dimension(:,:),allocatable :: son
        logical,dimension(:),allocatable :: ref

        ! Obtain details of the hydro variables stored
        call read_hydrofile_descriptor(repository,varIDs)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository,amr,sim)
        amr%lmax = amr%nlevelmax

        ! Compute the Hilbert curve
        call get_cpu_map(reg,amr)
        write(*,*)'ncpu_read:',amr%ncpu_read

        ! Just make sure that initial values are zero
        attrs%data = 0D0

        ! Allocate grids
        allocate(ngridfile(1:amr%ncpu+amr%nboundary,1:amr%nlevelmax))
        allocate(ngridlevel(1:amr%ncpu,1:amr%nlevelmax))
        if(amr%nboundary>0)allocate(ngridbound(1:amr%nboundary,1:amr%nlevelmax))

        ! Compute linear transformation
        trans_matrix = 0D0
        call new_z_coordinates(reg%axis,trans_matrix,roterr)
        if (roterr.eq.1) then
            write(*,*) 'Incorrect CS transformation!'
            stop
        endif

        ipos=INDEX(repository,'output_')
        nchar=repository(ipos+7:ipos+13)

        ! Loop over processor files
        cpuloop: do k=1,amr%ncpu_read
            icpu = amr%cpu_list(k)
            call title(icpu,ncharcpu)

            ! Open AMR file and skip header
            nomfich = TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
            open(unit=10,file=nomfich,status='old',form='unformatted')
            do i=1,21
                read(10) ! Skip header
            end do
            ! Read grid numbers
            read(10)ngridlevel
            ngridfile(1:amr%ncpu,1:amr%nlevelmax) = ngridlevel
            read(10) ! Skip
            if(amr%nboundary>0) then
                do i=1,2
                    read(10)
                end do
                read(10)ngridbound
                ngridfile(amr%ncpu+1:amr%ncpu+amr%nboundary,1:amr%nlevelmax) = ngridbound
            endif
            read(10) ! Skip
            ! R. Teyssier: comment the single following line for old stuff
            read(10)
            if(TRIM(amr%ordering).eq.'bisection')then
                do i=1,5
                    read(10)
                end do
            else
                read(10)
            endif
            read(10)
            read(10)
            read(10)

            ! Make sure that we are not trying to access to far in the refinement map…
            call check_lmax(ngridfile,amr)
            ! Open HYDRO file and skip header
            nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
            open(unit=11,file=nomfich,status='old',form='unformatted')
            read(11)
            read(11)nvarh
            read(11)
            read(11)
            read(11)
            read(11)

            ! Loop over levels
            levelloop: do ilevel=1,amr%lmax
                ! Geometry
                dx = 0.5**ilevel
                nx_full = 2**ilevel
                ny_full = 2**ilevel
                nz_full = 2**ilevel
                do ind=1,amr%twotondim
                    iz=(ind-1)/4
                    iy=(ind-1-4*iz)/2
                    ix=(ind-1-2*iy-4*iz)
                    xc(ind,1)=(dble(ix)-0.5D0)*dx
                    xc(ind,2)=(dble(iy)-0.5D0)*dx
                    xc(ind,3)=(dble(iz)-0.5D0)*dx
                end do
                
                ! Allocate work arrays
                ngrida = ngridfile(icpu,ilevel)
                if(ngrida>0) then
                    allocate(xg (1:ngrida,1:amr%ndim))
                    allocate(son(1:ngrida,1:amr%twotondim))
                    allocate(var(1:ngrida,1:amr%twotondim,1:nvarh))
                    allocate(x  (1:ngrida,1:amr%ndim))
                    allocate(ref(1:ngrida))
                endif

                ! Loop over domains
                domloop: do j=1,amr%nboundary+amr%ncpu
                    
                    ! Read AMR data
                    if (ngridfile(j,ilevel)>0) then
                        read(10) ! Skip grid index
                        read(10) ! Skip next index
                        read(10) ! Skip prev index

                        ! Read grid center
                        do idim=1,amr%ndim
                            if(j.eq.icpu)then
                                read(10)xg(:,idim)
                            else
                                read(10)
                            endif
                        end do

                        read(10) ! Skip father index
                        do ind=1,2*amr%ndim
                            read(10) ! Skip nbor index
                        end do

                        ! Read son index
                        do ind=1,amr%twotondim
                            if(j.eq.icpu)then
                                read(10)son(:,ind)
                            else
                                read(10)
                            end if
                        end do

                        ! Skip cpu map
                        do ind=1,amr%twotondim
                            read(10)
                        end do

                        ! Skip refinement map
                        do ind=1,amr%twotondim
                            read(10)
                        end do
                    endif

                    ! Read HYDRO data
                    read(11)
                    read(11)
                    if(ngridfile(j,ilevel)>0)then
                        ! Read hydro variables
                        tndimloop: do ind=1,amr%twotondim
                            varloop: do ivar=1,nvarh
                                if (j.eq.icpu) then
                                    read(11)var(:,ind,ivar)
                                else
                                    read(11)
                                endif
                            end do varloop
                        end do tndimloop
                    endif
                end do domloop

                ! Finally, get to every cell
                if (ngrida>0) then
                    ! Loop over cells
                    cellloop: do ind=1,amr%twotondim

                        ! Compute cell center
                        do i=1,ngrida
                            x(i,1)=(xg(i,1)+xc(ind,1)-amr%xbound(1))
                            x(i,2)=(xg(i,2)+xc(ind,2)-amr%xbound(2))
                            x(i,3)=(xg(i,3)+xc(ind,3)-amr%xbound(3))
                        end do

                        ! Check if cell is refined
                        do i=1,ngrida
                            ref(i) = son(i,ind)>0.and.ilevel<amr%lmax
                        end do
                        ngridaloop: do i=1,ngrida
                            ! Check if cell is inside the desired region
                            distance = 0D0
                            xtemp = x(i,:)
                            xtemp = xtemp - reg%centre
                            call rotate_vector(xtemp,trans_matrix)
                            x(i,:) = xtemp
                            call checkifinside(x(i,:),reg,ok_cell,distance)
                            ok_filter = filter_cell(varIDs,reg,filt,xtemp,dx,var(i,ind,:))
                            ok_cell= ok_cell.and..not.ref(i).and.ok_filter
                            if (ok_cell) then
                                vtemp = var(i,ind,varIDs%vx:varIDs%vz)
                                call rotate_vector(vtemp,trans_matrix)
                                var(i,ind,varIDs%vx:varIDs%vz) = vtemp
                                call extract_data(reg,varIDs,x(i,:),var(i,ind,:),dx,attrs)
                            endif
                        end do ngridaloop
                    end do cellloop

                    deallocate(xg,son,var,ref,x)
                endif
            end do levelloop
            close(10)
            close(11)
        end do cpuloop

        ! Finally just renormalise for weighted quantities
        call renormalise(attrs)

    end subroutine integrate_region

end module amr_integrator