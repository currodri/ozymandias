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
    use geometrical_regions

    type amr_region_attrs
        integer :: nvars,nfilter=1,nsubs=0
        character(128),dimension(:),allocatable :: varnames
        integer :: nwvars
        type(filter),dimension(:),allocatable :: filters
        type(region),dimension(:),allocatable :: subs
        character(128),dimension(:),allocatable :: wvarnames
        real(dbl),dimension(:,:,:,:),allocatable :: data
    end type amr_region_attrs

    contains

    subroutine allocate_amr_regions_attrs(attrs)
        implicit none
        type(amr_region_attrs),intent(inout) :: attrs

        if (.not.allocated(attrs%varnames)) allocate(attrs%varnames(attrs%nvars))
        if (.not.allocated(attrs%wvarnames)) allocate(attrs%wvarnames(attrs%nwvars))
        if (.not.allocated(attrs%data)) allocate(attrs%data(attrs%nfilter,attrs%nvars,attrs%nwvars,4))
        if (.not.allocated(attrs%filters)) allocate(attrs%filters(attrs%nfilter))
        if (.not.allocated(attrs%subs).and.attrs%nsubs>0) allocate(attrs%subs(attrs%nsubs))
    end subroutine allocate_amr_regions_attrs

    subroutine extract_data(reg,pos,cellvars,cellsons,cellsize,attrs,ifilt,trans_matrix,grav_var)
        use vectors
        implicit none

        ! Input/output variables
        type(region),intent(in) :: reg
        real(dbl),dimension(1:3),intent(in) :: pos
        real(dbl),dimension(0:amr%twondim,1:varIDs%nvar),intent(in) :: cellvars
        integer,dimension(0:amr%twondim),intent(in) :: cellsons
        real(dbl),intent(in) :: cellsize
        type(amr_region_attrs),intent(inout) :: attrs
        integer, intent(in) :: ifilt
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix
        real(dbl),dimension(0:amr%twondim,1:4),optional,intent(in) :: grav_var

        ! Local variables
        integer :: i,j
        real(dbl) :: ytemp,wtemp
        type(vector) :: x

        x = pos
        varloop: do i=1,attrs%nvars
            ! Get variable
            if (present(grav_var)) then
                call getvarvalue(reg,cellsize,x,cellvars,cellsons,attrs%varnames(i),ytemp,trans_matrix,grav_var)
            else
                call getvarvalue(reg,cellsize,x,cellvars,cellsons,attrs%varnames(i),ytemp,trans_matrix)
            end if
            wvarloop: do j=1,attrs%nwvars
                ! Get weights
                if (attrs%wvarnames(j)=='counts'.or.attrs%wvarnames(j)=='cumulative') then
                    wtemp =  1D0
                else
                    if (present(grav_var)) then
                        call getvarvalue(reg,cellsize,x,cellvars,cellsons,attrs%wvarnames(j),wtemp,trans_matrix,grav_var)
                    else
                        call getvarvalue(reg,cellsize,x,cellvars,cellsons,attrs%wvarnames(j),wtemp,trans_matrix)
                    endif
                endif
                
                ! Save to attrs
                attrs%data(ifilt,i,j,1) = attrs%data(ifilt,i,j,1) + ytemp*wtemp ! Value (weighted or not)
                if (attrs%data(ifilt,i,j,2).eq.0D0) then
                    attrs%data(ifilt,i,j,2) = ytemp ! Just to make sure that the initial min is not zero
                else
                    attrs%data(ifilt,i,j,2) = min(ytemp,attrs%data(ifilt,i,j,2))    ! Min value
                endif
                attrs%data(ifilt,i,j,3) = max(ytemp,attrs%data(ifilt,i,j,3))    ! Max value
                attrs%data(ifilt,i,j,4) = attrs%data(ifilt,i,j,4) + wtemp       ! Weight

            end do wvarloop
        end do varloop
    end subroutine extract_data

    subroutine renormalise(attrs)
        implicit none

        ! Input/output variables
        type(amr_region_attrs),intent(inout) :: attrs

        ! Local variables
        integer :: i,j,ifilt
        filterloop: do ifilt=1,attrs%nfilter
            varloop: do i=1,attrs%nvars
                wvarloop: do j=1,attrs%nwvars
                    if (attrs%wvarnames(j) /= 'cumulative') then
                        attrs%data(ifilt,i,j,1) = attrs%data(ifilt,i,j,1) / attrs%data(ifilt,i,j,4)
                    endif
                end do wvarloop
            end do varloop
        end do filterloop
    end subroutine renormalise

    subroutine integrate_region(repository,reg,use_neigh,attrs)
        use vectors
        use coordinate_systems
        implicit none

        ! Input/output variables
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        logical, intent(in) :: use_neigh
        type(amr_region_attrs),intent(inout) :: attrs

        ! Choose type if integrator
        if (use_neigh) then
            write(*,*)'Loading neighbours...'
            call integrate_region_neigh
        else
            call integrate_region_fast
        end if
    
    contains


        subroutine integrate_region_fast
            use vectors
            use coordinate_systems
            implicit none

            ! Specific variables for this subroutine
            integer :: i,j,k
            integer :: ipos,icpu,ilevel,ind,idim,ivar,ifilt,isub
            integer :: ix,iy,iz,ngrida,nx_full,ny_full,nz_full
            integer :: tot_pos,tot_ref,total_ncell
            integer :: nvarh
            integer :: roterr
            character(5) :: nchar,ncharcpu
            character(128) :: nomfich
            real(dbl) :: distance,dx
            type(vector) :: xtemp,vtemp,gtemp
            logical :: ok_cell,ok_filter,ok_cell_each,ok_sub,read_gravity
            integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
            real(dbl),dimension(1:8,1:3) :: xc
            real(dbl),dimension(3,3) :: trans_matrix
            real(dbl),dimension(:,:),allocatable :: xg,x,xorig
            real(dbl),dimension(:,:,:),allocatable :: var,grav_var
            real(dbl),dimension(:,:),allocatable :: tempvar
            real(dbl),dimension(:,:),allocatable :: tempgrav_var
            integer,dimension(:,:),allocatable :: son
            integer,dimension(:),allocatable :: tempson
            logical,dimension(:),allocatable :: ref

            total_ncell = 0
            tot_pos = 0
            tot_ref = 0

            ! Obtain details of the hydro variables stored
            call read_hydrofile_descriptor(repository)

            ! Initialise parameters of the AMR structure and simulation attributes
            call init_amr_read(repository)
            amr%lmax = amr%nlevelmax

            ! Compute the Hilbert curve
            call get_cpu_map(reg)
            write(*,*)'ncpu_read:',amr%ncpu_read

            ! Check whether we need to read the gravity files
            read_gravity = .false.
            do ivar=1,attrs%nvars
                if (attrs%varnames(ivar)(1:4) .eq. 'grav') then
                    read_gravity = .true.
                    write(*,*)'Reading gravity files...'
                    exit
                endif
            end do

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
                call check_lmax(ngridfile)
                ! Open HYDRO file and skip header
                nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
                open(unit=11,file=nomfich,status='old',form='unformatted')
                read(11)
                read(11)nvarh
                read(11)
                read(11)
                read(11)
                read(11)

                if (read_gravity) then
                    ! Open GRAV file and skip header
                    nomfich=TRIM(repository)//'/grav_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
                    open(unit=12,file=nomfich,status='old',form='unformatted')
                    read(12) !ncpu
                    read(12) !ndim
                    read(12) !nlevelmax
                    read(12) !nboundary 
                endif

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
                        if (read_gravity)allocate(grav_var(1:ngrida,1:amr%twotondim,1:4))
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

                        if (read_gravity) then
                            ! Read GRAV data
                            read(12)
                            read(12)
                            if(ngridfile(j,ilevel)>0)then
                                do ind=1,amr%twotondim
                                    if (j.eq.icpu) then
                                        read(12)grav_var(:,ind,1)
                                    else
                                        read(12)
                                    end if
                                    do ivar=1,amr%ndim
                                        if (j.eq.icpu) then
                                            read(12)grav_var(:,ind,ivar+1)
                                        else
                                            read(12)
                                        end if
                                    end do
                                end do
                            end if
                        end if
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
                            xorig  = x
                            ngridaloop: do i=1,ngrida
                                ! Check if cell is inside the desired region
                                distance = 0D0
                                xtemp = x(i,:)
                                xtemp = xtemp - reg%centre
                                call rotate_vector(xtemp,trans_matrix)
                                x(i,:) = xtemp
                                call checkifinside(x(i,:),reg,ok_cell,distance)

                                ! If we are avoiding substructure, check whether we are safe
                                if (attrs%nsubs>0) then
                                    ok_sub = .true.
                                    do isub=1,attrs%nsubs
                                        ok_sub = ok_sub .and. filter_sub(attrs%subs(isub),xorig(i,:))
                                    end do
                                    ok_cell = ok_cell .and. ok_sub
                                end if
                                if (ok_cell) then
                                    ! Transform position to galaxy frame
                                        xtemp = xorig(i,:)
                                        xtemp = xtemp - reg%centre
                                        call rotate_vector(xtemp,trans_matrix)
                                        ! Velocity transformed
                                        vtemp = var(i,ind,varIDs%vx:varIDs%vz)
                                        vtemp = vtemp - reg%bulk_velocity
                                        call rotate_vector(vtemp,trans_matrix)

                                        ! Gravitational acc
                                        if (read_gravity) then
                                            gtemp = grav_var(i,ind,2:4)
                                            call rotate_vector(gtemp,trans_matrix)
                                        endif
                                        allocate(tempvar(0:amr%twondim,nvarh))
                                        allocate(tempson(0:amr%twondim))
                                        if (read_gravity) allocate(tempgrav_var(0:amr%twondim,1:4))
                                        ! Just add central cell as we do not want neighbours
                                        tempvar(0,:) = var(i,ind,:)
                                        tempson(0)       = son(i,ind)
                                        if (read_gravity) tempgrav_var(0,:) = grav_var(i,ind,:)
                                        tempvar(0,varIDs%vx:varIDs%vz) = vtemp
                                        if (read_gravity) tempgrav_var(0,2:4) = gtemp
                                        if (ok_cell)tot_pos = tot_pos + 1
                                        if (.not.ref(i))tot_ref = tot_ref + 1
                                        filterloop: do ifilt=1,attrs%nfilter
                                            if (read_gravity) then
                                                ok_filter = filter_cell(reg,attrs%filters(ifilt),xtemp,dx,tempvar,&
                                                                        &tempson,trans_matrix,tempgrav_var)
                                            else
                                                ok_filter = filter_cell(reg,attrs%filters(ifilt),xtemp,dx,tempvar,&
                                                                        &tempson,trans_matrix)
                                            end if
                                            ok_cell_each= ok_cell.and..not.ref(i).and.ok_filter
                                            if (ok_cell_each) then
                                                total_ncell = total_ncell + 1
                                                if (read_gravity) then
                                                    call extract_data(reg,x(i,:),tempvar,tempson,dx,attrs,ifilt,trans_matrix,tempgrav_var)
                                                else
                                                    call extract_data(reg,x(i,:),tempvar,tempson,dx,attrs,ifilt,trans_matrix)
                                                endif
                                            endif
                                        end do filterloop
                                        deallocate(tempvar,tempson)
                                        if (read_gravity) deallocate(tempgrav_var)
                                end if
                            end do ngridaloop
                        end do cellloop
                        deallocate(xg,son,var,ref,x)
                        if (read_gravity) then
                            deallocate(grav_var)
                        end if
                    endif
                end do levelloop
                close(10)
                close(11)
            end do cpuloop

            ! Finally just renormalise for weighted quantities
            call renormalise(attrs)

            write(*,*)'Total number of cells used: ', total_ncell
            write(*,*)'Total number of cells refined: ', tot_ref
            write(*,*)'Total number of cells in region: ', tot_pos

        end subroutine integrate_region_fast

        subroutine integrate_region_neigh
            use vectors
            use coordinate_systems
            implicit none

            ! Specific variables for this subroutine
            integer :: i,j,k
            integer :: ipos,icpu,ilevel,ind,idim,ivar,ifilt,iskip,inbor,ison
            integer :: ix,iy,iz,ngrida,nx_full,ny_full,nz_full
            integer :: tot_pos,tot_ref,total_ncell
            integer :: nvarh
            integer :: roterr
            character(5) :: nchar,ncharcpu
            character(128) :: nomfich
            real(dbl) :: distance,dx
            type(vector) :: xtemp,vtemp,gtemp
            logical :: ok_cell,ok_filter,ok_cell_each,read_gravity
            integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
            real(dbl),dimension(1:8,1:3) :: xc
            real(dbl),dimension(1:3,1:3) :: trans_matrix
            real(dbl),dimension(:),allocatable :: xxg,son_dens
            real(dbl),dimension(:,:),allocatable :: x
            real(dbl),dimension(:,:),allocatable :: var
            real(dbl),dimension(:,:),allocatable :: grav_var
            real(dbl),dimension(:,:),allocatable :: tempvar
            real(dbl),dimension(:,:),allocatable :: tempgrav_var
            real(dbl),dimension(:,:),allocatable :: cellpos
            integer,dimension(:,:),allocatable :: nbor
            integer,dimension(:),allocatable :: son,tempson,iig
            integer,dimension(:),allocatable :: ind_cell,ind_cell2
            integer ,dimension(0:amr%twondim) :: ind_nbor
            logical,dimension(:),allocatable :: ref
            type(level),dimension(1:100) :: grid

            total_ncell = 0
            tot_pos = 0
            tot_ref = 0

            ! Obtain details of the hydro variables stored
            call read_hydrofile_descriptor(repository)

            ! Initialise parameters of the AMR structure and simulation attributes
            call init_amr_read(repository)
            amr%lmax = amr%nlevelmax

            ! Compute the Hilbert curve
            call get_cpu_map(reg)
            write(*,*)'ncpu_read:',amr%ncpu_read

            ! Check whether we need to read the gravity files
            read_gravity = .false.
            do ivar=1,attrs%nvars
                if (attrs%varnames(ivar)(1:4) .eq. 'grav' .or.&
                & trim(attrs%varnames(ivar)) .eq. 'neighbour_accuracy') then
                    read_gravity = .true.
                    write(*,*)'Reading gravity files...'
                    exit
                endif
            end do

            ! Just make sure that initial values are zero
            attrs%data = 0D0

            ! Allocate grids
            allocate(ngridfile(1:amr%ncpu+amr%nboundary,1:amr%nlevelmax))
            allocate(ngridlevel(1:amr%ncpu,1:amr%nlevelmax))
            if(amr%nboundary>0)allocate(ngridbound(1:amr%nboundary,1:amr%nlevelmax))
            ! Compute hierarchy
            do ilevel=1,amr%lmax
                grid(ilevel)%ngrid = 0
            end do

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

                allocate(nbor(1:amr%ngridmax,1:amr%twondim))
                allocate(son(1:amr%ncoarse+amr%twotondim*amr%ngridmax))
                son = 0; nbor = 0

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
                read(10)son(1:amr%ncoarse)
                read(10)
                read(10)

                ! Open HYDRO file and skip header
                nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
                open(unit=11,file=nomfich,status='old',form='unformatted')
                read(11)
                read(11)nvarh
                read(11)
                read(11)
                read(11)
                read(11)
                allocate(var(1:amr%ncoarse+amr%twotondim*amr%ngridmax,1:nvarh))
                allocate(cellpos(1:amr%ncoarse+amr%twotondim*amr%ngridmax,1:3))
                cellpos = 0d0; var = 0d0
                
                if (read_gravity) then
                    ! Open GRAV file and skip header
                    nomfich=TRIM(repository)//'/grav_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
                    open(unit=12,file=nomfich,status='old',form='unformatted')
                    read(12) !ncpu
                    read(12) !ndim
                    read(12) !nlevelmax
                    read(12) !nboundary 
                    allocate(grav_var(1:amr%ncoarse+amr%twotondim*amr%ngridmax,1:4))
                endif

                ! Loop over levels
                levelloop1: do ilevel=1,amr%lmax
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
                    grid(ilevel)%ngrid = 0
                    ! Loop over domains
                    domloop: do j=1,amr%nboundary+amr%ncpu
                        ! Allocate work arrays
                        ngrida = ngridfile(j,ilevel)
                        if(ngrida>0)then
                            if (allocated(grid(ilevel)%ind_grid)) deallocate(grid(ilevel)%ind_grid)
                            if (allocated(grid(ilevel)%xg)) deallocate(grid(ilevel)%xg)
                            allocate(grid(ilevel)%ind_grid(1:ngrida))
                            allocate(grid(ilevel)%xg (1:ngrida,1:amr%ndim))
                            allocate(iig(1:ngrida))
                            allocate(xxg(1:ngrida))
                            
                            ! Read AMR data
                            read(10) grid(ilevel)%ind_grid
                            read(10) ! Skip next index
                            read(10) ! Skip prev index
                            if(j.eq.icpu) then
                                if (allocated(grid(ilevel)%real_ind)) deallocate(grid(ilevel)%real_ind)
                                allocate(grid(ilevel)%real_ind(1:ngrida))
                                grid(ilevel)%real_ind = grid(ilevel)%ind_grid
                                grid(ilevel)%ngrid = ngridfile(j,ilevel)
                            end if

                            ! Read grid center
                            do idim=1,amr%ndim
                                read(10)xxg
                                do i=1,ngrida
                                    grid(ilevel)%xg(i,idim) = xxg(i)
                                end do
                            end do

                            read(10) ! Skip father index
                            ! Read nbor index
                            do ind=1,amr%twondim
                                read(10)iig
                                do i=1,ngrida
                                    nbor(grid(ilevel)%ind_grid(i),ind) = iig(i)
                                end do
                            end do
                            ! Read son index
                            do ind=1,amr%twotondim
                                iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                                read(10)iig
                                do i=1,ngrida
                                    son(grid(ilevel)%ind_grid(i)+iskip) = iig(i)
                                end do
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
                        if(ngrida>0)then
                            ! Read hydro variables
                            tndimloop: do ind=1,amr%twotondim
                                iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                                varloop: do ivar=1,nvarh
                                    read(11)xxg
                                    do i=1,ngrida
                                        var(grid(ilevel)%ind_grid(i)+iskip,ivar) = xxg(i)
                                    end do
                                end do varloop
                            end do tndimloop
                        endif

                        if (read_gravity) then
                            ! Read GRAV data
                            read(12)
                            read(12)
                            if(ngrida>0)then
                                do ind=1,amr%twotondim
                                    iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                                    read(12)xxg
                                    do i=1,ngrida
                                        grav_var(grid(ilevel)%ind_grid(i)+iskip,1) = xxg(i)
                                    end do
                                    do ivar=1,amr%ndim
                                        read(12)xxg
                                        do i=1,ngrida
                                            grav_var(grid(ilevel)%ind_grid(i)+iskip,ivar+1) = xxg(i)
                                        end do
                                    end do
                                end do
                            end if
                        end if

                        !Compute positions
                        if(ngrida>0)then
                            do ind=1,amr%twotondim
                                iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                                do i=1,ngrida
                                    do ivar=1,amr%ndim
                                        cellpos(grid(ilevel)%ind_grid(i)+iskip,ivar)=(grid(ilevel)%xg(i,ivar)+xc(ind,ivar)-amr%xbound(ivar))
                                    end do
                                end do
                            end do
                        end if
                        if (ngrida>0) deallocate(iig,xxg)
                    end do domloop
                end do levelloop1
                close(10)
                close(11)
                if (read_gravity) then
                    close(12)
                end if
                ! Loop over levels again now with arrays fully filled
                levelloop2: do ilevel=1,amr%lmax
                    ! Geometry
                    dx = 0.5**ilevel
                    nx_full = 2**ilevel
                    ny_full = 2**ilevel

                    ! Allocate work arrays
                    ngrida = grid(ilevel)%ngrid
                    if(ngrida>0)then
                        allocate(ind_cell(1:ngrida))
                        allocate(x  (1:ngrida,1:amr%ndim))
                        allocate(ref(1:ngrida))
                    endif
                    if (ngrida>0) then
                        ! Loop over cells
                        cellloop: do ind=1,amr%twotondim

                            ! Get cell indexes
                            iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                            do i=1,ngrida
                                ind_cell(i) = iskip+grid(ilevel)%real_ind(i)
                            end do
                            ! Compute cell center
                            do i=1,ngrida
                                x(i,:)=cellpos(ind_cell(i),:)
                            end do

                            ! Check if cell is refined
                            do i=1,ngrida
                                ref(i) = son(ind_cell(i))>0.and.ilevel<amr%lmax
                                ! Look for cells that just got refined and de-refine them
                                if (ref(i)) then
                                    allocate(son_dens(1:amr%twotondim))
                                    do inbor=1,amr%twotondim
                                        ison = son(ind_cell(i)) + amr%ncoarse+(inbor-1)*amr%ngridmax
                                        son_dens(inbor) = var(ison,1)
                                    end do
                                    if (all(son_dens == var(ind_cell(i),1))) then
                                        do inbor=1,amr%twotondim
                                            ison = son(ind_cell(i)) + amr%ncoarse+(inbor-1)*amr%ngridmax
                                            son(ison) = son(ind_cell(i))
                                        end do
                                        son(ind_cell(i)) = 0
                                    end if
                                    deallocate(son_dens)
                                end if
                            end do

                            ngridaloop: do i=1,ngrida
                                ! Check if cell is inside the desired region
                                distance = 0D0
                                xtemp = x(i,:)
                                xtemp = xtemp - reg%centre
                                call rotate_vector(xtemp,trans_matrix)
                                x(i,:) = xtemp
                                call checkifinside(x(i,:),reg,ok_cell,distance)
                                
                                ! Velocity transformed --> ONLY FOR CENTRAL CELL
                                vtemp = var(ind_cell(i),varIDs%vx:varIDs%vz)
                                vtemp = vtemp - reg%bulk_velocity
                                call rotate_vector(vtemp,trans_matrix)

                                ! Gravitational acc --> ONLY FOR CENTRAL CELL
                                if (read_gravity) then
                                    gtemp = grav_var(ind_cell(i),2:4)
                                    call rotate_vector(gtemp,trans_matrix)
                                endif

                                ! Get neighbours
                                allocate(ind_cell2(1))
                                ind_cell2(1) = ind_cell(i)
                                call getnbor(son,nbor,ind_cell2,ind_nbor,1)
                                deallocate(ind_cell2)
                                allocate(tempvar(0:amr%twondim,nvarh))
                                allocate(tempson(0:amr%twondim))
                                if (read_gravity) allocate(tempgrav_var(0:amr%twondim,1:4))
                                ! Just correct central cell vectors for the region
                                tempvar(0,:) = var(ind_nbor(0),:)
                                tempson(0)       = son(ind_nbor(0))
                                if (read_gravity) tempgrav_var(0,:) = grav_var(ind_nbor(0),:)
                                tempvar(0,varIDs%vx:varIDs%vz) = vtemp
                                if (read_gravity) tempgrav_var(0,2:4) = gtemp
                                do inbor=0,amr%twondim
                                    tempvar(inbor,:) = var(ind_nbor(inbor),:)
                                    tempson(inbor)       = son(ind_nbor(inbor))
                                    if (read_gravity) tempgrav_var(inbor,:) = grav_var(ind_nbor(inbor),:)
                                end do
                                if (ok_cell)tot_pos = tot_pos + 1
                                if (.not.ref(i))tot_ref = tot_ref + 1
                                filterloop: do ifilt=1,attrs%nfilter
                                    if (read_gravity) then
                                        ok_filter = filter_cell(reg,attrs%filters(ifilt),xtemp,dx,tempvar,&
                                                                &tempson,trans_matrix,tempgrav_var)
                                    else
                                        ok_filter = filter_cell(reg,attrs%filters(ifilt),xtemp,dx,tempvar,&
                                                                &tempson,trans_matrix)
                                    end if
                                    ok_cell_each= ok_cell.and..not.ref(i).and.ok_filter
                                    if (ok_cell_each) then
                                        total_ncell = total_ncell + 1
                                        if (read_gravity) then
                                            call extract_data(reg,x(i,:),tempvar,tempson,dx,attrs,ifilt,trans_matrix,tempgrav_var)
                                        else
                                            call extract_data(reg,x(i,:),tempvar,tempson,dx,attrs,ifilt,trans_matrix)
                                        endif
                                    endif
                                end do filterloop
                                deallocate(tempvar,tempson)
                                if (read_gravity) deallocate(tempgrav_var)
                            end do ngridaloop
                        end do cellloop
                        deallocate(ref,x,ind_cell)
                    endif
                end do levelloop2
                deallocate(nbor,son,var,cellpos)
                close(10)
                close(11)
                if (read_gravity) then
                    close(12)
                    deallocate(grav_var)
                end if
            end do cpuloop

            ! Finally just renormalise for weighted quantities
            call renormalise(attrs)

            write(*,*)'Total number of cells used: ', total_ncell
            write(*,*)'Total number of cells refined: ', tot_ref
            write(*,*)'Total number of cells in region: ', tot_pos

        end subroutine integrate_region_neigh
    end subroutine integrate_region

end module amr_integrator