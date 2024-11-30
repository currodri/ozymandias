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
    use dictionary_commons
    use io_ramses
    use hydro_commons
    use filtering

    type chunk_handler
        integer :: nvars,nx=100,ny=100,nz=100
        character(128),dimension(:),allocatable :: varnames
        type(filter) :: filt
        real(dbl),dimension(:,:,:,:),allocatable :: data
        type(hydro_var),dimension(:),allocatable :: vars
    end type chunk_handler

    contains

    subroutine allocate_chunk_handler(chunk)
        implicit none
        type(chunk_handler),intent(inout) :: chunk

        if (.not.allocated(chunk%varnames)) allocate(chunk%varnames(1:chunk%nvars))
        if (.not.allocated(chunk%vars)) allocate(chunk%vars(1:chunk%nvars))
        if (.not.allocated(chunk%data)) allocate(chunk%data(1:chunk%nvars,1:chunk%nx,1:chunk%ny,1:chunk%nz))
    end subroutine allocate_chunk_handler

    ! subroutine get_unigrid_old(repository,reg,chunk)
    !     use vectors
    !     use coordinate_systems
    !     use geometrical_regions
    !     implicit none

    !     ! Input/output variables
    !     character(128),intent(in) :: repository
    !     type(region),intent(inout) :: reg
    !     type(chunk_handler),intent(inout) :: chunk

    !     ! Specific variables for this subroutine
    !     integer :: i,j,k
    !     integer :: ipos,icpu,ilevel,ind,idim,ivar,ifilt
    !     integer :: ix,iy,iz,ixp1,iyp1,izp1
    !     integer :: ngrida,nx_full,ny_full,nz_full,ndom
    !     integer :: imin,imax,jmin,jmax,kmin,kmax
    !     integer :: nvarh
    !     integer :: roterr
    !     character(5) :: nchar,ncharcpu
    !     character(128) :: nomfich
    !     real(dbl) :: distance,dx,value
    !     real(dbl) :: xmin,ymin,zmin
    !     real(dbl) :: ddx,ddy,ddz,dex,dey,dez,xx,yy,zz
    !     type(vector) :: xtemp,vtemp
    !     logical :: ok_cell,ok_filter,ok_cell_each
    !     integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
    !     real(dbl),dimension(1:8,1:3) :: xc
    !     real(dbl),dimension(3,3) :: trans_matrix
    !     real(dbl),dimension(:,:),allocatable :: xg,x
    !     real(dbl),dimension(:,:,:),allocatable :: var
    !     integer,dimension(:,:),allocatable :: son
    !     logical,dimension(:),allocatable :: ref

    !     type(level),dimension(1:100) :: grid

    !     ! Obtain details of the hydro variables stored
    !     call read_hydrofile_descriptor(repository)

    !     ! Initialise parameters of the AMR structure and simulation attributes
    !     call init_amr_read(repository)
    !     amr%lmax = amr%nlevelmax

    !     ! Compute the Hilbert curve
    !     if (trim(reg%name).ne.'cube') then
    !         write(*,*) 'amr_chunk only supports the extraction of a cube! Check your region!'
    !         stop
    !     endif
    !     call get_cpu_map(reg)
    !     write(*,*)'ncpu_read:',amr%ncpu_read

    !     ! Just make sure that initial values are zero
    !     chunk%data = 0D0

    !     ! Allocate grids
    !     allocate(ngridfile(1:amr%ncpu+amr%nboundary,1:amr%nlevelmax))
    !     allocate(ngridlevel(1:amr%ncpu,1:amr%nlevelmax))
    !     if(amr%nboundary>0)allocate(ngridbound(1:amr%nboundary,1:amr%nlevelmax))

    !     ! Compute linear transformation
    !     trans_matrix = 0D0
    !     call new_z_coordinates(reg%axis,trans_matrix,roterr)
    !     if (roterr.eq.1) then
    !         write(*,*) 'Incorrect CS transformation!'
    !         stop
    !     endif

    !     ! Compute hierarchy
    !     do ilevel=1,amr%lmax
    !         nx_full = 2**ilevel
    !         ny_full = 2**ilevel
    !         imin = int(0D0*dble(nx_full))+1
    !         imax = int((reg%xmax-reg%xmin)*dble(nx_full))+1
    !         jmin = int(0D0*dble(ny_full))+1
    !         jmax = int((reg%ymax-reg%ymin)*dble(ny_full))+1
    !         kmin = int(0D0*dble(nz_full))+1
    !         kmax = int((reg%zmax-reg%zmin)*dble(nz_full))+1
    !         allocate(grid(ilevel)%cube(1:chunk%nvars,imin:imax,jmin:jmax,kmin:kmax))
    !         grid(ilevel)%cube = 0D0
    !         grid(ilevel)%imin = imin
    !         grid(ilevel)%imax = imax
    !         grid(ilevel)%jmin = jmin
    !         grid(ilevel)%jmax = jmax
    !         grid(ilevel)%kmin = kmin
    !         grid(ilevel)%kmax = kmax
    !     end do

    !     ipos=INDEX(repository,'output_')
    !     nchar=repository(ipos+7:ipos+13)

    !     ! Loop over processor files
    !     cpuloop: do k=1,amr%ncpu_read
    !         icpu = amr%cpu_list(k)
    !         call title(icpu,ncharcpu)

    !         ! Open AMR file and skip header
    !         nomfich = TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
    !         open(unit=10,file=nomfich,status='old',form='unformatted')
    !         do i=1,21
    !             read(10) ! Skip header
    !         end do
    !         ! Read grid numbers
    !         read(10)ngridlevel
    !         ngridfile(1:amr%ncpu,1:amr%nlevelmax) = ngridlevel
    !         read(10) ! Skip
    !         if(amr%nboundary>0) then
    !             do i=1,2
    !                 read(10)
    !             end do
    !             read(10)ngridbound
    !             ngridfile(amr%ncpu+1:amr%ncpu+amr%nboundary,1:amr%nlevelmax) = ngridbound
    !         endif
    !         read(10) ! Skip
    !         ! R. Teyssier: comment the single following line for old stuff
    !         read(10)
    !         if(TRIM(amr%ordering).eq.'bisection')then
    !             do i=1,5
    !                 read(10)
    !             end do
    !         else
    !             read(10)
    !         endif
    !         read(10)
    !         read(10)
    !         read(10)

    !         ! Make sure that we are not trying to access to far in the refinement map…
    !         call check_lmax(ngridfile)
    !         ! Open HYDRO file and skip header
    !         nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
    !         open(unit=11,file=nomfich,status='old',form='unformatted')
    !         read(11)
    !         read(11)nvarh
    !         read(11)
    !         read(11)
    !         read(11)
    !         read(11)

    !         ! Loop over levels
    !         levelloop: do ilevel=1,amr%lmax
    !             ! Geometry
    !             dx = 0.5**ilevel
    !             nx_full = 2**ilevel
    !             ny_full = 2**ilevel
    !             nz_full = 2**ilevel
    !             do ind=1,amr%twotondim
    !                 iz=(ind-1)/4
    !                 iy=(ind-1-4*iz)/2
    !                 ix=(ind-1-2*iy-4*iz)
    !                 xc(ind,1)=(dble(ix)-0.5D0)*dx
    !                 xc(ind,2)=(dble(iy)-0.5D0)*dx
    !                 xc(ind,3)=(dble(iz)-0.5D0)*dx
    !             end do
                
    !             ! Allocate work arrays
    !             ngrida = ngridfile(icpu,ilevel)
    !             grid(ilevel)%ngrid = ngrida
    !             if(ngrida>0) then
    !                 allocate(xg (1:ngrida,1:amr%ndim))
    !                 allocate(son(1:ngrida,1:amr%twotondim))
    !                 allocate(var(1:ngrida,1:amr%twotondim,1:nvarh))
    !                 allocate(x  (1:ngrida,1:amr%ndim))
    !                 allocate(ref(1:ngrida))
    !             endif

    !             ! Loop over domains
    !             domloop: do j=1,amr%nboundary+amr%ncpu
                    
    !                 ! Read AMR data
    !                 if (ngridfile(j,ilevel)>0) then
    !                     read(10) ! Skip grid index
    !                     read(10) ! Skip next index
    !                     read(10) ! Skip prev index

    !                     ! Read grid center
    !                     do idim=1,amr%ndim
    !                         if(j.eq.icpu)then
    !                             read(10)xg(:,idim)
    !                         else
    !                             read(10)
    !                         endif
    !                     end do

    !                     read(10) ! Skip father index
    !                     do ind=1,2*amr%ndim
    !                         read(10) ! Skip nbor index
    !                     end do

    !                     ! Read son index
    !                     do ind=1,amr%twotondim
    !                         if(j.eq.icpu)then
    !                             read(10)son(:,ind)
    !                         else
    !                             read(10)
    !                         end if
    !                     end do

    !                     ! Skip cpu map
    !                     do ind=1,amr%twotondim
    !                         read(10)
    !                     end do

    !                     ! Skip refinement map
    !                     do ind=1,amr%twotondim
    !                         read(10)
    !                     end do
    !                 endif

    !                 ! Read HYDRO data
    !                 read(11)
    !                 read(11)
    !                 if(ngridfile(j,ilevel)>0)then
    !                     ! Read hydro variables
    !                     tndimloop: do ind=1,amr%twotondim
    !                         varloop: do ivar=1,nvarh
    !                             if (j.eq.icpu) then
    !                                 read(11)var(:,ind,ivar)
    !                             else
    !                                 read(11)
    !                             endif
    !                         end do varloop
    !                     end do tndimloop
    !                 endif
    !             end do domloop

    !             ! Finally, get to every cell
    !             if (ngrida>0) then
    !                 ! Loop over cells
    !                 cellloop: do ind=1,amr%twotondim

    !                     ! Compute cell center
    !                     do i=1,ngrida
    !                         x(i,1)=(xg(i,1)+xc(ind,1)-amr%xbound(1))
    !                         x(i,2)=(xg(i,2)+xc(ind,2)-amr%xbound(2))
    !                         x(i,3)=(xg(i,3)+xc(ind,3)-amr%xbound(3))
    !                     end do

    !                     ! Check if cell is refined
    !                     do i=1,ngrida
    !                         ref(i) = son(i,ind)>0.and.ilevel<amr%lmax
    !                     end do
    !                     ngridaloop: do i=1,ngrida
    !                         ! Check if cell is inside the desired region
    !                         distance = 0D0
    !                         xtemp = x(i,:)
    !                         xtemp = xtemp - reg%centre
    !                         call rotate_vector(xtemp,trans_matrix)
    !                         x(i,:) = xtemp
    !                         call checkifinside(x(i,:),reg,ok_cell,distance)
    !                         vtemp = var(ind_cell(i),varIDs%vx:varIDs%vz)
    !                         call rotate_vector(vtemp,trans_matrix)
    !                         var(ind_cell(i),varIDs%vx:varIDs%vz) = vtemp
    !                         ! Check filter
    !                         ok_filter = filter_cell(reg,chunk%filt,xtemp,dx,var(i,ind,:))
    !                         ok_cell_each= ok_cell.and..not.ref(i).and.ok_filter
    !                         if (ok_cell_each) then
    !                             ix = int((x(i,1)+0.5*(reg%xmax-reg%xmin))*dble(nx_full)) + 1
    !                             iy = int((x(i,2)+0.5*(reg%ymax-reg%ymin))*dble(ny_full)) + 1
    !                             iz = int((x(i,3)+0.5*(reg%zmax-reg%zmin))*dble(nz_full)) + 1
    !                             if( ix>=grid(ilevel)%imin.and.&
    !                                 & iy>=grid(ilevel)%jmin.and.&
    !                                 & iz>=grid(ilevel)%kmin.and.&
    !                                 & ix<=grid(ilevel)%imax.and.&
    !                                 & iy<=grid(ilevel)%jmax.and.&
    !                                 & iz<=grid(ilevel)%kmax) then
    !                                 amrvarloop: do ivar=1,chunk%nvars
    !                                     call getvarvalue(reg,dx,xtemp,var(i,ind,:),chunk%varnames(ivar),value)
    !                                     grid(ilevel)%cube(ivar,ix,iy,iz) = value
    !                                 end do amrvarloop
    !                             endif
    !                         endif
    !                     end do ngridaloop
    !                 end do cellloop

    !                 deallocate(xg,son,var,ref,x)
    !             endif
    !         end do levelloop
    !         close(10)
    !         close(11)
    !     end do cpuloop

    !     ! Upload to maximum level (lmax)
    !     nx_full = 2**amr%lmax
    !     ny_full = 2**amr%lmax
    !     nz_full = 2**amr%lmax
    !     imin = int(0D0*dble(nx_full))+1
    !     imax = int((reg%xmax-reg%xmin)*dble(nx_full))
    !     jmin = int(0D0*dble(ny_full))+1
    !     jmax = int((reg%ymax-reg%ymin)*dble(ny_full))
    !     kmin = int(0D0*dble(nz_full))+1
    !     kmax = int((reg%zmax-reg%zmin)*dble(nz_full))
    !     xloop: do ix = imin,imax
    !         xmin = ((ix-0.5)/2**amr%lmax)
    !         yloop: do iy=jmin,jmax
    !             ymin=((iy-0.5)/2**amr%lmax)
    !             zloop: do iz=kmin,kmax
    !                 zmin=((iz-0.5)/2**amr%lmax)
    !                 ilevelloop: do ilevel=1,amr%lmax-1
    !                     ndom = 2**ilevel
    !                     i = int(xmin*ndom)+1
    !                     j = int(ymin*ndom)+1
    !                     k = int(zmin*ndom)+1
    !                     projvarlooplmax: do ivar=1,chunk%nvars
    !                         grid(amr%lmax)%cube(ivar,ix,iy,iz) = grid(amr%lmax)%cube(ivar,ix,iy,iz) + &
    !                                                             & grid(ilevel)%cube(ivar,ix,iy,iz)
    !                     end do projvarlooplmax
    !                 end do ilevelloop
    !             end do zloop
    !        end do yloop
    !     end do xloop

    !     ! Adapt lmax grid to the required chunk resolution
    !     do i=1,chunk%nx
    !         xx = (dble(i)-0.5)/dble(chunk%nx)*dble(imax-imin+1)
    !         ix = int(xx)
    !         ddx = xx-ix
    !         dex = 1d0-ddx
    !         ix = ix+imin
    !         ix = min(ix,imax)
    !         ixp1 = min(ix+1,imax)
    !         do j=1,chunk%ny
    !             yy = (dble(j)-0.5)/dble(chunk%ny)*dble(jmax-jmin+1)
    !             iy = int(yy)
    !             ddy = yy-iy
    !             dey = 1d0-ddy
    !             iy = iy+jmin
    !             iy = min(iy,jmax)
    !             iyp1 = min(iy+1,jmax)
    !             do k=1,chunk%nz
    !                 zz = (dble(k)-0.5)/dble(chunk%nz)*dble(kmax-kmin+1)
    !                 iz = int(zz)
    !                 ddz = zz-iz
    !                 dez = 1d0-ddz
    !                 iz = iz+kmin
    !                 iz = min(iz,kmax)
    !                 izp1 = min(iz+1,kmax)
    !                 do ivar=1,chunk%nvars
    !                     chunk%data(ivar,i,j,k)=dex*dey*dez*grid(amr%lmax)%cube(ivar,ix  ,iy  ,iz  )+ &
    !                                     &      ddx*dey*dez*grid(amr%lmax)%cube(ivar,ixp1,iy  ,iz  )+ &
    !                                     &      dex*ddy*dez*grid(amr%lmax)%cube(ivar,ix  ,iyp1,iz  )+ &
    !                                     &      dex*dey*ddz*grid(amr%lmax)%cube(ivar,ix  ,iy  ,izp1)+ &
    !                                     &      ddx*ddy*dez*grid(amr%lmax)%cube(ivar,ixp1,iyp1,iz  )+ &
    !                                     &      dex*ddy*ddz*grid(amr%lmax)%cube(ivar,ix  ,iyp1,izp1)+ &
    !                                     &      ddx*dey*ddz*grid(amr%lmax)%cube(ivar,ixp1,iy  ,izp1)+ &
    !                                     &      ddx*ddy*ddz*grid(amr%lmax)%cube(ivar,ixp1,iyp1,izp1)
    !                 end do
    !             end do
    !         end do
    !     end do
    ! end subroutine get_unigrid_old

    subroutine get_unigrid(repository,reg,lmax,symlog,chunk,vardict)
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
        type(dictf90),intent(in),optional :: vardict

        ! Specific variables for this subroutine
        integer :: i,j,k
        integer :: ipos,icpu,ilevel,ind,idim,ivar,ifilt,iskip,inbor
        integer :: ix,iy,iz,ixp1,iyp1,izp1,cumngrida
        integer :: ngrida,nx_full,ny_full,nz_full,ndom
        integer :: imin,imax,jmin,jmax,kmin,kmax
        integer :: nvarh
        integer :: roterr
        integer :: ivx,ivy,ivz
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
        real(dbl),dimension(:,:),allocatable :: var
        real(dbl),dimension(:,:),allocatable :: tempvar
        integer,dimension(:,:),allocatable :: nbor
        integer,dimension(:),allocatable :: son,tempson
        integer,dimension(:),allocatable :: ind_grid,ind_cell,ind_cell2
        integer ,dimension(0:amr%twondim)::ind_nbor
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

        ! Set up hydro variables quicklook tools
        if (present(vardict)) then
            ! If the user provides their own variable dictionary,
            ! use that one instead of the automatic from the 
            ! hydro descriptor file (RAMSES)
            call get_var_tools(vardict,chunk%nvars,chunk%varnames,chunk%vars)

            ! We always need the indexes of the velocities
            ! to perform rotations of gas velocities
            ivx = vardict%get('velocity_x')
            ivy = vardict%get('velocity_y')
            ivz = vardict%get('velocity_z')
        else
            call get_var_tools(varIDs,chunk%nvars,chunk%varnames,chunk%vars)

            ! We always need the indexes of the velocities
            ! to perform rotations of gas velocities
            ivx = varIDs%get('velocity_x')
            ivy = varIDs%get('velocity_y')
            ivz = varIDs%get('velocity_z')
        end if

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

            allocate(nbor(1:amr%ngridmax,1:amr%twondim))
            allocate(son(1:amr%ncoarse+amr%twotondim*amr%ngridmax))
            cumngrida = 0

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
                    allocate(ind_grid(1:ngrida))
                    allocate(ind_cell(1:ngrida))
                    allocate(xg (1:ngrida,1:amr%ndim))
                    allocate(x  (1:ngrida,1:amr%ndim))
                    allocate(ref(1:ngrida))
                endif
                ! Loop over domains
                domloop: do j=1,amr%nboundary+amr%ncpu
                    ! Read AMR data
                    if (ngridfile(j,ilevel)>0) then
                        if(j.eq.icpu)then
                            read(10) ind_grid
                        else
                            read(10)
                        end if
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
                        do ind=1,amr%twondim
                            if(j.eq.icpu)then
                                read(10)nbor(ind_grid,ind)
                            else
                                read(10)
                            end if
                        end do

                        ! Read son index
                        do ind=1,amr%twotondim
                            iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                            if(j.eq.icpu)then
                                read(10)son(ind_grid+iskip)
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
                            iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                            varloop: do ivar=1,nvarh
                                if (j.eq.icpu) then
                                    read(11)var(ind_grid+iskip,ivar)
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
                        iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                        do i=1,ngrida
                            ref(i) = son(ind_grid(i)+iskip)>0.and.ilevel<amr%lmax
                        end do
                        
                        ! Get cell indexes
                        do i=1,ngrida
                            ind_cell(i) = iskip+ind_grid(i)
                        end do
                        ngridaloop: do i=1,ngrida
                            ! Check if cell is inside the desired region
                            distance = 0D0
                            xtemp = x(i,:)
                            xtemp = xtemp - reg%centre
                            newx = xtemp
                            call rotate_vector(xtemp,trans_matrix)
                            call checkifinside(newx,reg,ok_cell,distance)

                            ! Velocity transformed --> ONLY FOR CENTRAL CELL
                            vtemp = var(ind_cell(i),ivx:ivz)
                            vtemp = vtemp - reg%bulk_velocity
                            call rotate_vector(vtemp,trans_matrix)
                            var(ind_cell(i),ivx:ivz) = vtemp

                            ! Get neighbours
                            allocate(ind_cell2(1))
                            ind_cell2(1) = ind_cell(i)
                            call getnbor(son,nbor,ind_cell2,ind_nbor,1)
                            deallocate(ind_cell2)
                            allocate(tempvar(0:amr%twondim,nvarh))
                            allocate(tempson(0:amr%twondim))
                            do inbor=0,amr%twondim
                                tempvar(inbor,:) = var(ind_nbor(inbor),:)
                                tempson(inbor)       = son(ind_nbor(inbor))
                            end do

                            ! Check filter
                            ok_filter = filter_cell(reg,chunk%filt,xtemp,dx,tempvar,&
                                                    &tempson,trans_matrix)
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
                                        grid(ilevel)%cube(ivar,ix,iy,iz) = chunk%vars(ivar)%myfunction(amr,sim,chunk%vars(ivar),&
                                                                                            reg,dx,xtemp,tempvar,tempson,trans_matrix)
                                    end do amrvarloop
                                endif
                            endif
                            deallocate(tempvar,tempson)
                        end do ngridaloop
                    end do cellloop
                    deallocate(xg,ref,x,ind_grid,ind_cell)
                endif
                cumngrida = cumngrida + ngrida
            end do levelloop
            deallocate(nbor,son,var)
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

    subroutine amr2skirt(repository,reg,filt,varname,outpath,vardict)
        use vectors
        use coordinate_systems
        use geometrical_regions
        implicit none
        
        ! Input/output variables
        character(128),intent(in) :: repository
        type(region),intent(in) :: reg
        type(filter),intent(inout) :: filt
        character(128),intent(in) :: varname
        character(128),intent(in) :: outpath
        type(dictf90),intent(in),optional :: vardict

        ! Specific variables for this subroutine
        integer :: i,j,k
        integer :: ipos,icpu,ilevel,ind,idim,ivar,ifilt,iskip,inbor
        integer :: ix,iy,iz,ngrida,nx_full,ny_full,nz_full,cumngrida
        integer :: tot_pos,tot_ref,total_ncell,cpu_ncell
        integer :: nvarh
        integer :: roterr
        integer :: ivx,ivy,ivz
        character(5) :: nchar,ncharcpu
        character(128) :: nomfich
        real(dbl) :: distance,dx,dx2kpc
        real(dbl) :: myval
        type(vector) :: xtemp,vtemp,xtempmin,xtempmax
        logical :: ok_cell,ok_filter,ok_cell_each,read_gravity
        integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
        real(dbl),dimension(1:8,1:3) :: xc
        real(dbl),dimension(3,3) :: trans_matrix
        real(dbl),dimension(:,:),allocatable :: xg,x
        real(dbl),dimension(:,:),allocatable :: var
        real(dbl),dimension(:,:),allocatable :: grav_var
        real(dbl),dimension(:,:),allocatable :: tempvar
        real(dbl),dimension(:,:),allocatable :: tempgrav_var
        integer,dimension(:,:),allocatable :: nbor
        integer,dimension(:),allocatable :: son,tempson
        integer,dimension(:),allocatable :: ind_grid,ind_cell,ind_cell2
        integer ,dimension(:),allocatable :: ind_nbor
        logical,dimension(:),allocatable :: ref
        character(128),dimension(1) :: vname
        type(hydro_var),dimension(1) :: hvar

        total_ncell = 0
        cpu_ncell = 0

        ! Obtain details of the hydro variables stored
        call read_hydrofile_descriptor(repository)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = amr%nlevelmax

        ! Compute the Hilbert curve
        call get_cpu_map(reg)
        write(*,*)'ncpu_read:',amr%ncpu_read

        ! Set up hydro variables quicklook tools
        if (present(vardict)) then
            ! If the user provides their own variable dictionary,
            ! use that one instead of the automatic from the 
            ! hydro descriptor file (RAMSES)
            call get_var_tools(vardict,1,vname,hvar)

            ! We also do it for the filter variables
            call get_filter_var_tools(vardict,filt)

            ! We always need the indexes of the velocities
            ! to perform rotations of gas velocities
            ivx = vardict%get('velocity_x')
            ivy = vardict%get('velocity_y')
            ivz = vardict%get('velocity_z')
        else
            call get_var_tools(varIDs,1,vname,hvar)

            ! We also do it for the filter variables
            call get_filter_var_tools(varIDs,filt)

            ! We always need the indexes of the velocities
            ! to perform rotations of gas velocities
            ivx = varIDs%get('velocity_x')
            ivy = varIDs%get('velocity_y')
            ivz = varIDs%get('velocity_z')
        end if

        ! Check whether we need to read the gravity files
        read_gravity = .false.
        if (varname(1:4) .eq. 'grav') then
            read_gravity = .true.
            write(*,*)'Reading gravity files...'
        end if

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

        ! Open output file and add header for SKIRT format
        open(unit=7,file=TRIM(outpath),form='formatted')
        write(7,98)
        write(7,99)TRIM(varname)
        write(7,97)
        write(7,101)
        write(7,97)
        97 format('#')
        98 format('# Gas cell data for simulated galaxy in RAMSES simulation')
        99 format('# SKIRT 9 import format for a medium source using ',A,' method')
        101 format('# Column 1: x-min (kpc)',/, &
                    '# Column 2: y-min (kpc)',/,&
                    '# Column 3: z-min (kpc)',/,&
                    '# Column 4: x-max (kpc)',/, &
                    '# Column 5: y-max (kpc)',/,&
                    '# Column 6: z-max (kpc)',/,&
                    '# Column 7: dust mass density (Msun/pc3)')

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
            ! call check_lmax(ngridfile)

            allocate(nbor(1:amr%ngridmax,1:amr%twondim))
            allocate(son(1:amr%ncoarse+amr%twotondim*amr%ngridmax))
            son = 0; nbor = 0
            cumngrida = 0

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
            levelloop: do ilevel=1,amr%lmax
                ! Geometry
                dx = 0.5**ilevel
                dx2kpc = dx * (sim%unit_l*sim%boxlen*cm2kpc)
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
                    allocate(ind_grid(1:ngrida))
                    allocate(ind_cell(1:ngrida))
                    allocate(xg (1:ngrida,1:amr%ndim))
                    allocate(x  (1:ngrida,1:amr%ndim))
                    allocate(ref(1:ngrida))
                endif
                !write(*,*)'Level allocation fine'
                ! Loop over domains
                domloop: do j=1,amr%nboundary+amr%ncpu
                    
                    ! Read AMR data
                    if (ngridfile(j,ilevel)>0) then
                        if(j.eq.icpu)then
                            read(10) ind_grid
                        else
                            read(10)
                        end if
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
                        do ind=1,amr%twondim
                            if(j.eq.icpu)then
                                read(10)nbor(ind_grid,ind)
                            else
                                read(10)
                            end if
                        end do

                        ! Read son index
                        do ind=1,amr%twotondim
                            iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                            if(j.eq.icpu)then
                                read(10)son(ind_grid+iskip)
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
                            iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                            varloop: do ivar=1,nvarh
                                if (j.eq.icpu) then
                                    read(11)var(ind_grid+iskip,ivar)
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
                                iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                                if (j.eq.icpu) then
                                    read(12)grav_var(ind_grid+iskip,1)
                                    do ivar=1,amr%ndim
                                        read(12)grav_var(ind_grid+iskip,ivar+1)
                                    end do
                                else
                                    read(12)
                                    do ivar=1,amr%ndim
                                        read(12)
                                    end do
                                end if
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
                        iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                        do i=1,ngrida
                            ref(i) = son(ind_grid(i)+iskip)>0.and.ilevel<amr%lmax
                        end do

                        ! Get cell indexes
                        do i=1,ngrida
                            ind_cell(i) = iskip+ind_grid(i)
                        end do

                        ngridaloop: do i=1,ngrida
                            ! Check if cell is inside the desired region
                            distance = 0D0
                            xtemp = x(i,:)
                            xtemp = xtemp - reg%centre
                            x(i,:) = xtemp
                            call checkifinside(x(i,:),reg,ok_cell,distance)
                            ! Velocity transformed --> ONLY FOR CENTRAL CELL
                            vtemp = var(ind_cell(i),ivx:ivz)
                            vtemp = vtemp - reg%bulk_velocity
                            var(ind_cell(i),ivx:ivz) = vtemp

                            ! Gravitational acc --> ONLY FOR CENTRAL CELL
                            if (read_gravity) then
                                vtemp = grav_var(ind_cell(i),2:4)
                                call rotate_vector(vtemp,trans_matrix)
                                grav_var(ind_cell(i),2:4) = vtemp
                            endif

                            ! Get neighbours
                            allocate(ind_cell2(1))
                            ind_cell2(1) = ind_cell(i)
                            allocate(ind_nbor(0:amr%twondim))
                            call getnbor(son,nbor,ind_cell2,ind_nbor,1)
                            deallocate(ind_cell2)
                            allocate(tempvar(0:amr%twondim,nvarh))
                            allocate(tempson(0:amr%twondim))
                            if (read_gravity) allocate(tempgrav_var(0:amr%twondim,1:4))
                            do inbor=0,amr%twondim
                                tempvar(inbor,:) = var(ind_nbor(inbor),:)
                                tempson(inbor)       = son(ind_nbor(inbor))
                                if (read_gravity) tempgrav_var(inbor,:) = grav_var(ind_nbor(inbor),:)
                            end do
                            deallocate(ind_nbor)

                            if (read_gravity) then
                                ok_filter = filter_cell(reg,filt,xtemp,dx,tempvar,tempson,&
                                                        &trans_matrix,tempgrav_var)
                            else
                                ok_filter = filter_cell(reg,filt,xtemp,dx,tempvar,tempson,&
                                                        &trans_matrix)
                            end if
                            ok_cell_each= ok_cell.and..not.ref(i).and.ok_filter
                            cpu_ncell = cpu_ncell + 1
                            if (ok_cell_each) then
                                total_ncell = total_ncell + 1
                                if (read_gravity) then
                                    myval = hvar(1)%myfunction(amr,sim,hvar(1),reg,dx,xtemp,tempvar,&
                                                            tempson,trans_matrix,tempgrav_var)
                                else
                                    myval = hvar(1)%myfunction(amr,sim,hvar(1),reg,dx,xtemp,tempvar,&
                                                            tempson,trans_matrix)
                                endif
                                ! Position to kpc
                                xtemp = xtemp * (sim%unit_l*cm2kpc)
                                x(i,:) = xtemp
                                xtempmin = x(i,:) - dx2kpc/2D0
                                xtempmax = x(i,:) + dx2kpc/2D0
                                ! TODO: Change for the different methods
                                myval = myval * (sim%unit_d*gcm32msunpc3)
                                write(7,100)xtempmin%x,xtempmin%y,xtempmin%z,&
                                            xtempmax%x,xtempmax%y,xtempmax%z,&
                                            myval
                                100 format(6F10.6,F16.12)
                            endif
                            deallocate(tempvar,tempson)
                            if (read_gravity) deallocate(tempgrav_var)
                        end do ngridaloop
                    end do cellloop
                    deallocate(xg,ref,x,ind_grid,ind_cell)
                endif
                cumngrida = cumngrida + ngrida
            end do levelloop
            deallocate(nbor,son,var)
            close(10)
            close(11)
            if (read_gravity) then
                close(12)
                deallocate(grav_var)
            end if
        end do cpuloop

        close(7)
        write(*,*)'Number of cells in the CPU read: ', cpu_ncell
        write(*,*)'Total number of cells used:      ', total_ncell
    end subroutine amr2skirt
end module export_amr