!--------------------------------------------------------------------------
! ozymandias:export_module.f90
!--------------------------------------------------------------------------
!
! MODULE: export_amr
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> types and routines useful for extracting the raw AMR data in RAMSES.
!
!> @details  
!> 
! 
!
!> @date 18/01/2022   0.1.0 taking guidelines from original amr2cube.f90
!--------------------------------------------------------------------------

module export_amr
    use local
    use io_ramses
    use filtering

    type chunk_handler
        integer :: nvars,nx=100,ny=100,nz=100
        character(128),dimension(:),allocatable :: varnames
        type(filter) :: filt
        real(dbl),dimension(:,:,:,:),allocatable :: data
    end type chunk_handler

    contains

    subroutine allocate_chunk_handler(chunk)
        implicit none
        type(chunk_handler),intent(inout) :: chunk

        if (.not.allocated(chunk%varnames)) allocate(chunk%varnames(1:chunk%nvars))
        if (.not.allocated(chunk%data)) allocate(chunk%data(1:chunk%nvars,1:chunk%nx,1:chunk%ny,1:chunk%nz))
    end subroutine allocate_chunk_handler

    subroutine get_unigrid_old(repository,reg,chunk)
        use vectors
        use coordinate_systems
        use geometrical_regions
        implicit none

        ! Input/output variables
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        type(chunk_handler),intent(inout) :: chunk

        ! Specific variables for this subroutine
        integer :: i,j,k
        integer :: ipos,icpu,ilevel,ind,idim,ivar,ifilt
        integer :: ix,iy,iz,ixp1,iyp1,izp1
        integer :: ngrida,nx_full,ny_full,nz_full,ndom
        integer :: imin,imax,jmin,jmax,kmin,kmax
        integer :: nvarh
        integer :: roterr
        character(5) :: nchar,ncharcpu
        character(128) :: nomfich
        real(dbl) :: distance,dx,value
        real(dbl) :: xmin,ymin,zmin
        real(dbl) :: ddx,ddy,ddz,dex,dey,dez,xx,yy,zz
        type(vector) :: xtemp,vtemp
        logical :: ok_cell,ok_filter,ok_cell_each
        integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
        real(dbl),dimension(1:8,1:3) :: xc
        real(dbl),dimension(3,3) :: trans_matrix
        real(dbl),dimension(:,:),allocatable :: xg,x
        real(dbl),dimension(:,:,:),allocatable :: var
        integer,dimension(:,:),allocatable :: son
        logical,dimension(:),allocatable :: ref

        type(level),dimension(1:100) :: grid

        ! Obtain details of the hydro variables stored
        call read_hydrofile_descriptor(repository)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = amr%nlevelmax

        ! Compute the Hilbert curve
        if (trim(reg%name).ne.'cube') then
            write(*,*) 'amr_chunk only supports the extraction of a cube! Check your region!'
            stop
        endif
        call get_cpu_map(reg)
        write(*,*)'ncpu_read:',amr%ncpu_read

        ! Just make sure that initial values are zero
        chunk%data = 0D0

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

        ! Compute hierarchy
        do ilevel=1,amr%lmax
            nx_full = 2**ilevel
            ny_full = 2**ilevel
            imin = int(0D0*dble(nx_full))+1
            imax = int((reg%xmax-reg%xmin)*dble(nx_full))+1
            jmin = int(0D0*dble(ny_full))+1
            jmax = int((reg%ymax-reg%ymin)*dble(ny_full))+1
            kmin = int(0D0*dble(nz_full))+1
            kmax = int((reg%zmax-reg%zmin)*dble(nz_full))+1
            allocate(grid(ilevel)%cube(1:chunk%nvars,imin:imax,jmin:jmax,kmin:kmax))
            grid(ilevel)%cube = 0D0
            grid(ilevel)%imin = imin
            grid(ilevel)%imax = imax
            grid(ilevel)%jmin = jmin
            grid(ilevel)%jmax = jmax
            grid(ilevel)%kmin = kmin
            grid(ilevel)%kmax = kmax
        end do

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
                grid(ilevel)%ngrid = ngrida
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
                            vtemp = var(i,ind,varIDs%vx:varIDs%vz)
                            call rotate_vector(vtemp,trans_matrix)
                            var(i,ind,varIDs%vx:varIDs%vz) = vtemp
                            ! Check filter
                            ok_filter = filter_cell(reg,chunk%filt,xtemp,dx,var(i,ind,:))
                            ok_cell_each= ok_cell.and..not.ref(i).and.ok_filter
                            if (ok_cell_each) then
                                ix = int((x(i,1)+0.5*(reg%xmax-reg%xmin))*dble(nx_full)) + 1
                                iy = int((x(i,2)+0.5*(reg%ymax-reg%ymin))*dble(ny_full)) + 1
                                iz = int((x(i,3)+0.5*(reg%zmax-reg%zmin))*dble(nz_full)) + 1
                                if( ix>=grid(ilevel)%imin.and.&
                                    & iy>=grid(ilevel)%jmin.and.&
                                    & iz>=grid(ilevel)%kmin.and.&
                                    & ix<=grid(ilevel)%imax.and.&
                                    & iy<=grid(ilevel)%jmax.and.&
                                    & iz<=grid(ilevel)%kmax) then
                                    amrvarloop: do ivar=1,chunk%nvars
                                        call getvarvalue(reg,dx,xtemp,var(i,ind,:),chunk%varnames(ivar),value)
                                        grid(ilevel)%cube(ivar,ix,iy,iz) = value
                                    end do amrvarloop
                                endif
                            endif
                        end do ngridaloop
                    end do cellloop

                    deallocate(xg,son,var,ref,x)
                endif
            end do levelloop
            close(10)
            close(11)
        end do cpuloop

        ! Upload to maximum level (lmax)
        nx_full = 2**amr%lmax
        ny_full = 2**amr%lmax
        nz_full = 2**amr%lmax
        imin = int(0D0*dble(nx_full))+1
        imax = int((reg%xmax-reg%xmin)*dble(nx_full))
        jmin = int(0D0*dble(ny_full))+1
        jmax = int((reg%ymax-reg%ymin)*dble(ny_full))
        kmin = int(0D0*dble(nz_full))+1
        kmax = int((reg%zmax-reg%zmin)*dble(nz_full))
        xloop: do ix = imin,imax
            xmin = ((ix-0.5)/2**amr%lmax)
            yloop: do iy=jmin,jmax
                ymin=((iy-0.5)/2**amr%lmax)
                zloop: do iz=kmin,kmax
                    zmin=((iz-0.5)/2**amr%lmax)
                    ilevelloop: do ilevel=1,amr%lmax-1
                        ndom = 2**ilevel
                        i = int(xmin*ndom)+1
                        j = int(ymin*ndom)+1
                        k = int(zmin*ndom)+1
                        projvarlooplmax: do ivar=1,chunk%nvars
                            grid(amr%lmax)%cube(ivar,ix,iy,iz) = grid(amr%lmax)%cube(ivar,ix,iy,iz) + &
                                                                & grid(ilevel)%cube(ivar,ix,iy,iz)
                        end do projvarlooplmax
                    end do ilevelloop
                end do zloop
           end do yloop
        end do xloop

        ! Adapt lmax grid to the required chunk resolution
        do i=1,chunk%nx
            xx = (dble(i)-0.5)/dble(chunk%nx)*dble(imax-imin+1)
            ix = int(xx)
            ddx = xx-ix
            dex = 1d0-ddx
            ix = ix+imin
            ix = min(ix,imax)
            ixp1 = min(ix+1,imax)
            do j=1,chunk%ny
                yy = (dble(j)-0.5)/dble(chunk%ny)*dble(jmax-jmin+1)
                iy = int(yy)
                ddy = yy-iy
                dey = 1d0-ddy
                iy = iy+jmin
                iy = min(iy,jmax)
                iyp1 = min(iy+1,jmax)
                do k=1,chunk%nz
                    zz = (dble(k)-0.5)/dble(chunk%nz)*dble(kmax-kmin+1)
                    iz = int(zz)
                    ddz = zz-iz
                    dez = 1d0-ddz
                    iz = iz+kmin
                    iz = min(iz,kmax)
                    izp1 = min(iz+1,kmax)
                    do ivar=1,chunk%nvars
                        chunk%data(ivar,i,j,k)=dex*dey*dez*grid(amr%lmax)%cube(ivar,ix  ,iy  ,iz  )+ &
                                        &      ddx*dey*dez*grid(amr%lmax)%cube(ivar,ixp1,iy  ,iz  )+ &
                                        &      dex*ddy*dez*grid(amr%lmax)%cube(ivar,ix  ,iyp1,iz  )+ &
                                        &      dex*dey*ddz*grid(amr%lmax)%cube(ivar,ix  ,iy  ,izp1)+ &
                                        &      ddx*ddy*dez*grid(amr%lmax)%cube(ivar,ixp1,iyp1,iz  )+ &
                                        &      dex*ddy*ddz*grid(amr%lmax)%cube(ivar,ix  ,iyp1,izp1)+ &
                                        &      ddx*dey*ddz*grid(amr%lmax)%cube(ivar,ixp1,iy  ,izp1)+ &
                                        &      ddx*ddy*ddz*grid(amr%lmax)%cube(ivar,ixp1,iyp1,izp1)
                    end do
                end do
            end do
        end do
    end subroutine get_unigrid_old

    subroutine get_unigrid(repository,reg,lmax,symlog,chunk)
        use vectors
        use coordinate_systems
        use geometrical_regions
        implicit none

        ! Input/output variables
        character(128),intent(in) :: repository
        type(region),intent(inout) :: reg
        integer,intent(in) :: lmax
        logical,intent(in) :: symlog
        type(chunk_handler),intent(inout) :: chunk

        ! Specific variables for this subroutine
        integer :: i,j,k
        integer :: ipos,icpu,ilevel,ind,idim,ivar,ifilt
        integer :: ix,iy,iz,ixp1,iyp1,izp1
        integer :: ngrida,nx_full,ny_full,nz_full,ndom
        integer :: imin,imax,jmin,jmax,kmin,kmax
        integer :: nvarh
        integer :: roterr
        character(5) :: nchar,ncharcpu
        character(128) :: nomfich
        real(dbl) :: distance,dx,value,signto
        real(dbl) :: xmin,ymin,zmin
        real(dbl) :: ddx,ddy,ddz,dex,dey,dez,xx,yy,zz
        real(dbl),dimension(1:3) :: newx 
        type(vector) :: xtemp,vtemp
        logical :: ok_cell,ok_filter,ok_cell_each
        integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
        real(dbl),dimension(1:8,1:3) :: xc
        real(dbl),dimension(3,3) :: trans_matrix
        real(dbl),dimension(:,:),allocatable :: xg,x
        real(dbl),dimension(:,:,:),allocatable :: var
        integer,dimension(:,:),allocatable :: son
        logical,dimension(:),allocatable :: ref

        type(level),dimension(1:100) :: grid

        ! Obtain details of the hydro variables stored
        call read_hydrofile_descriptor(repository)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = lmax
        if (lmax.eq.0) amr%lmax = amr%nlevelmax

        ! Compute the Hilbert curve
        if (trim(reg%name).ne.'cube') then
            write(*,*) 'amr_chunk only supports the extraction of a cube! Check your region!'
            stop
        endif
        call get_cpu_map(reg)
        write(*,*)'ncpu_read:',amr%ncpu_read

        ! Just make sure that initial values are zero
        chunk%data = 0D0

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

        ! Compute hierarchy
        do ilevel=1,amr%lmax
            nx_full = 2**ilevel
            ny_full = 2**ilevel
            nz_full = 2**ilevel
            imin = int(reg%xmin*dble(nx_full))+1
            imax = int(reg%xmax*dble(nx_full))+1
            jmin = int(reg%ymin*dble(ny_full))+1
            jmax = int(reg%ymax*dble(ny_full))+1
            kmin = int(reg%zmin*dble(nz_full))+1
            kmax = int(reg%zmax*dble(nz_full))+1
            allocate(grid(ilevel)%cube(1:chunk%nvars,imin:imax,jmin:jmax,kmin:kmax))
            grid(ilevel)%cube = 0D0
            grid(ilevel)%imin = imin
            grid(ilevel)%imax = imax
            grid(ilevel)%jmin = jmin
            grid(ilevel)%jmax = jmax
            grid(ilevel)%kmin = kmin
            grid(ilevel)%kmax = kmax
        end do

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
                grid(ilevel)%ngrid = ngrida
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
                            newx = xtemp
                            call rotate_vector(xtemp,trans_matrix)
                            call checkifinside(newx,reg,ok_cell,distance)
                            vtemp = var(i,ind,varIDs%vx:varIDs%vz)
                            call rotate_vector(vtemp,trans_matrix)
                            var(i,ind,varIDs%vx:varIDs%vz) = vtemp
                            ! Check filter
                            ok_filter = filter_cell(reg,chunk%filt,xtemp,dx,var(i,ind,:))
                            ok_cell_each= ok_cell.and..not.ref(i).and.ok_filter
                            if (ok_cell_each) then
                                ix=int(x(i,1)*dble(nx_full))+1
                                iy=int(x(i,2)*dble(ny_full))+1
                                iz=int(x(i,3)*dble(nz_full))+1
                                if( ix>=grid(ilevel)%imin.and.&
                                    & iy>=grid(ilevel)%jmin.and.&
                                    & iz>=grid(ilevel)%kmin.and.&
                                    & ix<=grid(ilevel)%imax.and.&
                                    & iy<=grid(ilevel)%jmax.and.&
                                    & iz<=grid(ilevel)%kmax) then
                                    amrvarloop: do ivar=1,chunk%nvars
                                        call getvarvalue(reg,dx,xtemp,var(i,ind,:),chunk%varnames(ivar),value)
                                        grid(ilevel)%cube(ivar,ix,iy,iz) = value
                                    end do amrvarloop
                                endif
                            endif
                        end do ngridaloop
                    end do cellloop

                    deallocate(xg,son,var,ref,x)
                endif
            end do levelloop
            close(10)
            close(11)
        end do cpuloop

        ! Upload to maximum level (lmax)
        nx_full = 2**amr%lmax
        ny_full = 2**amr%lmax
        nz_full = 2**amr%lmax
        imin = int(reg%xmin*dble(nx_full))+1
        imax = int(reg%xmax*dble(nx_full))
        jmin = int(reg%ymin*dble(ny_full))+1
        jmax = int(reg%ymax*dble(ny_full))
        kmin = int(reg%zmin*dble(nz_full))+1
        kmax = int(reg%zmax*dble(nz_full))
        xloop: do ix = imin,imax
            xmin = ((ix-0.5)/2**amr%lmax)
            yloop: do iy=jmin,jmax
                ymin=((iy-0.5)/2**amr%lmax)
                zloop: do iz=kmin,kmax
                    zmin=((iz-0.5)/2**amr%lmax)
                    ilevelloop: do ilevel=1,amr%lmax-1
                        ndom = 2**ilevel
                        i = int(xmin*ndom)+1
                        j = int(ymin*ndom)+1
                        k = int(zmin*ndom)+1
                        gridvarlooplmax: do ivar=1,chunk%nvars
                            grid(amr%lmax)%cube(ivar,ix,iy,iz) = grid(amr%lmax)%cube(ivar,ix,iy,iz) + &
                                                                & grid(ilevel)%cube(ivar,i,j,k)
                        end do gridvarlooplmax
                    end do ilevelloop
                end do zloop
           end do yloop
        end do xloop
        
        write(*,*)'Min value:',minval(grid(amr%lmax)%cube(1,imin:imax,jmin:jmax,kmin:kmax))
        write(*,*)'Max value:',maxval(grid(amr%lmax)%cube(1,imin:imax,jmin:jmax,kmin:kmax))
        write(*,*)'Norm:     ',sum   (dble(grid(amr%lmax)%cube(1,imin:imax,jmin:jmax,kmin:kmax))) &
       & /(imax-imin+1)/(jmax-jmin+1)/(kmax-kmin+1)
        ! Adapt lmax grid to the required chunk resolution
        do i=1,chunk%nx
            xx=(dble(i)-0.5)/dble(chunk%nx)*dble(imax-imin+1)
            ix=int(xx)
            ddx=xx-ix
            dex=1d0-ddx
            ix=ix+imin
            ix=min(ix,imax)
            ixp1=min(ix+1,imax)
            do j=1,chunk%ny
                yy=(dble(j)-0.5)/dble(chunk%ny)*dble(jmax-jmin+1)
                iy=int(yy)
                ddy=yy-iy
                dey=1d0-ddy
                iy=iy+jmin
                iy=min(iy,jmax)
                iyp1=min(iy+1,jmax)
                do k=1,chunk%nz
                    zz=(dble(k)-0.5)/dble(chunk%nz)*dble(kmax-kmin+1)
                    iz=int(zz)
                    ddz=zz-iz
                    dez=1d0-ddz
                    iz=iz+kmin
                    iz=min(iz,kmax)
                    izp1=min(iz+1,kmax)
                    do ivar=1,chunk%nvars
                        SignConserve: if (symlog) then
                            if (dex*dey*dez*(grid(amr%lmax)%cube(ivar,ix  ,iy  ,iz  ))+ &
                                 &      ddx*dey*dez*(grid(amr%lmax)%cube(ivar,ixp1,iy  ,iz  ))+ &
                                 &      dex*ddy*dez*(grid(amr%lmax)%cube(ivar,ix  ,iyp1,iz  ))+ &
                                 &      dex*dey*ddz*(grid(amr%lmax)%cube(ivar,ix  ,iy  ,izp1))+ &
                                 &      ddx*ddy*dez*(grid(amr%lmax)%cube(ivar,ixp1,iyp1,iz  ))+ &
                                 &      dex*ddy*ddz*(grid(amr%lmax)%cube(ivar,ix  ,iyp1,izp1))+ &
                                 &      ddx*dey*ddz*(grid(amr%lmax)%cube(ivar,ixp1,iy  ,izp1))+ &
                                 &      ddx*ddy*ddz*(grid(amr%lmax)%cube(ivar,ixp1,iyp1,izp1)) .ge. 0) then
                               signto = 1D0
                            else
                               signto = -1D0
                            end if
                         end if SignConserve
                        chunk%data(ivar,i,j,k)=dex*dey*dez*log10(abs(grid(amr%lmax)%cube(ivar,ix  ,iy  ,iz  )))+ &
                                        &      ddx*dey*dez*log10(abs(grid(amr%lmax)%cube(ivar,ixp1,iy  ,iz  )))+ &
                                        &      dex*ddy*dez*log10(abs(grid(amr%lmax)%cube(ivar,ix  ,iyp1,iz  )))+ &
                                        &      dex*dey*ddz*log10(abs(grid(amr%lmax)%cube(ivar,ix  ,iy  ,izp1)))+ &
                                        &      ddx*ddy*dez*log10(abs(grid(amr%lmax)%cube(ivar,ixp1,iyp1,iz  )))+ &
                                        &      dex*ddy*ddz*log10(abs(grid(amr%lmax)%cube(ivar,ix  ,iyp1,izp1)))+ &
                                        &      ddx*dey*ddz*log10(abs(grid(amr%lmax)%cube(ivar,ixp1,iy  ,izp1)))+ &
                                        &      ddx*ddy*ddz*log10(abs(grid(amr%lmax)%cube(ivar,ixp1,iyp1,izp1)))
                        chunk%data(ivar,i,j,k) = signto * 10.**(chunk%data(ivar,i,j,k))
                    end do
                end do
            end do
        end do
    end subroutine get_unigrid
end module export_amr