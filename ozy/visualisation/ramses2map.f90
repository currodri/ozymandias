module obs_instruments
    use local
    use vectors

    type camera
        type(vector) :: centre,los_axis,up_vector
        real(dbl),dimension(1:2) :: region_size
        real(dbl) :: distance,far_cut_depth
        integer :: map_max_size=1024
    end type camera

    type(vector),private :: x_axis,y_axis,z_axis

    contains

    real function log2(x)
        implicit none
        real(dbl), intent(in) :: x

        log2 = log(x) / log(2D0)
    end function

    type(camera) function init_camera(centre,los_axis,up_vector,region_size,distance,far_cut_depth,map_max_size)
        implicit none
        type(vector),intent(in) :: centre,los_axis,up_vector
        real(dbl),dimension(1:2),intent(in) :: region_size
        real(dbl),intent(in) :: distance,far_cut_depth
        integer,intent(in) :: map_max_size

        init_camera%centre = centre
        
        init_camera%los_axis = los_axis
        init_camera%up_vector = up_vector

        init_camera%region_size = region_size
        init_camera%distance = distance
        init_camera%far_cut_depth = far_cut_depth
        init_camera%map_max_size = map_max_size
    end function init_camera

    integer function get_required_resolution(cam)
        implicit none
        type(camera),intent(in) :: cam
        get_required_resolution = int(ceiling(log2(cam%map_max_size/maxval(cam%region_size))))
    end function get_required_resolution

    subroutine get_map_size(cam,n_map)
        implicit none
        type(camera),intent(in) :: cam
        integer,dimension(1:2),intent(inout) :: n_map
        real(dbl) :: aspect_ratio

        aspect_ratio = dble(cam%region_size(1)) / dble(cam%region_size(2))
        if (aspect_ratio > 1D0) then
            n_map(1) = cam%map_max_size
            n_map(2) = int(nint(cam%map_max_size/aspect_ratio))
        else
            n_map(2) = cam%map_max_size
            n_map(1) = int(nint(cam%map_max_size*aspect_ratio))
        endif
    end subroutine get_map_size

    subroutine get_map_box(cam,box)
        use geometrical_regions
        implicit none
        type(camera),intent(in) :: cam
        type(region),intent(inout) :: box
        real(dbl) :: dx,dy

        dx = cam%region_size(1); dy = cam%region_size(2)
        box%centre = cam%centre
        box%axis = cam%los_axis
        box%xmin = -dx/2D0; box%ymin = -dy/2D0; box%zmin = -cam%far_cut_depth
        box%xmax = dx/2D0; box%ymax = dy/2D0; box%zmax = cam%distance
    end subroutine get_map_box

    subroutine get_pixel_edges(cam,x_edges,y_edges)
        implicit none
        type(camera),intent(in) :: cam
        real(dbl),dimension(:),allocatable,intent(inout) :: x_edges,y_edges
        real(dbl) :: dx,dy,step
        integer,dimension(1:2) :: n_map
        integer :: i

        dx = cam%region_size(1); dy = cam%region_size(2)
        call get_map_size(cam,n_map)
        allocate(x_edges(1:n_map(1)+1))
        step = dx/(n_map(1)+1)
        do i=0,n_map(1)+1
            x_edges(i) = -dx/2D0 + step*i
        end do 
        allocate(y_edges(1:n_map(2)+1))
        step = dy/(n_map(1)+1)
        do i=0,n_map(2)+1
            y_edges(i) = -dy/2D0 + step*i
        end do
    end subroutine get_pixel_edges

    subroutine get_camera_basis(cam,cam_basis)
        use basis_representations
        implicit none
        type(camera),intent(in) :: cam
        type(basis),intent(inout) :: cam_basis

        cam_basis%u(1) = cam%up_vector * cam%los_axis
        cam_basis%u(1) = cam_basis%u(1) / magnitude(cam_basis%u(1))
        cam_basis%u(2) = cam%up_vector
        cam_basis%u(3) = cam%los_axis
    end subroutine get_camera_basis

    subroutine los_transformation(cam,trans_matrix)
        use basis_representations
        implicit none
        type(camera),intent(in) :: cam
        real(dbl),dimension(1:3,1:3),intent(inout) :: trans_matrix
        type(basis) :: cam_basis
        integer :: i

        trans_matrix = 0D0
        call get_camera_basis(cam,cam_basis)
        do i=1,3
            trans_matrix(i,:) = cam_basis%u(i)
        end do
    end subroutine los_transformation

    subroutine get_bounding_box(cam,bbox)
        use geometrical_regions
        implicit none
        type(camera),intent(in) :: cam
        type(region),intent(inout) :: bbox
        type(region) :: box
        real(dbl),dimension(1:2,1:3) :: box_bounds
        real(dbl),dimension(1:3,1:2) :: coords_by_axis
        real(dbl),dimension(1:2) :: coord_array
        real(dbl),dimension(1:8,1:3) :: points, xform_corners
        real(dbl),dimension(1:3) :: xform_min,xform_max
        integer,dimension(1:3) :: shapes,nperiods,nrepeats
        integer :: i,npoints
        integer :: ipoint,iperiod,irepeat,icoord,idim

        call get_map_box(cam,box)
        box_bounds(1,:) = (/box%xmin,box%ymin,box%zmin/)
        box_bounds(2,:) = (/box%xmax,box%ymax,box%zmax/)
        coords_by_axis = transpose(box_bounds)
        shapes = (/2,2,2/)
        nperiods = (/2,4,8/)
        nperiods = (/ (int(nperiods(i)/shapes(i)),i=1,3) /)
        npoints = product(shapes)
        nrepeats = (/ (int(npoints/(nperiods(i)*shapes(i))), i=1,3) /)
        points = 0D0
        do idim=1,3
            coord_array = coords_by_axis(idim,:)
            ipoint = 1
            icoord = 1
            do iperiod=1,nperiods(idim)
                do icoord=1,shapes(idim)
                    do irepeat=1,nrepeats(idim)
                        points(ipoint,idim) = coord_array(icoord)
                        ipoint = ipoint + 1
                    end do
                end do
            end do
        end do
        call deproject_points(cam,npoints,points)
        xform_corners = points
        xform_min = (/ (minval(xform_corners(:,i)), i=1,3) /)
        xform_max = (/ (maxval(xform_corners(:,i)), i=1,3) /)
        xform_min = (/ (max(0D0,xform_min(i)), i=1,3) /)
        xform_max = (/ (min(1D0,xform_max(i)), i=1,3) /)

        bbox%xmin = xform_min(1);bbox%ymin = xform_min(2);bbox%zmin = xform_min(3)
        bbox%xmax = xform_max(1);bbox%ymax = xform_max(2);bbox%zmax = xform_max(3)
    end subroutine get_bounding_box

    subroutine deproject_points(cam,npoints,points)
        use rotations
        use basis_representations
        implicit none
        type(camera),intent(in) :: cam
        integer,intent(in) :: npoints
        real(dbl),dimension(1:npoints,1:3),intent(inout) :: points
        type(vector) :: temp_vec
        integer :: i
        real(dbl),dimension(1:3,1:3) :: trans_matrix
        type(basis) :: cam_basis

        call los_transformation(cam,trans_matrix)

        do i=1,npoints
            temp_vec = points(i,:)
            call rotate_vector(temp_vec,transpose(trans_matrix))
            temp_vec = temp_vec + cam%centre
            points(i,:) = temp_vec
        end do
    end subroutine deproject_points

    subroutine project_points(cam,npoints,points)
        use rotations
        use basis_representations
        implicit none
        type(camera),intent(in) :: cam
        integer,intent(in) :: npoints
        real(dbl),dimension(1:npoints,1:3),intent(inout) :: points
        type(vector) :: temp_vec
        integer :: i,j
        real(dbl),dimension(1:3,1:3) :: trans_matrix
        type(basis) :: cam_basis

        call los_transformation(cam,trans_matrix)

        do i=1,npoints
            temp_vec = points(i,:)
            temp_vec = temp_vec - cam%centre
            call rotate_vector(temp_vec,trans_matrix)
            points(i,:) = temp_vec
        end do
    end subroutine project_points

end module obs_instruments

module maps
    use local
    use vectors
    use rotations
    use io_ramses
    use geometrical_regions
    use obs_instruments

    type projection_handler
        character(128) :: pov
        integer :: nvars
        character(128),dimension(:),allocatable :: varnames
        character(128) :: weightvar
        real(dbl),dimension(:,:,:),allocatable :: toto
    end type projection_handler

    contains

    subroutine allocate_projection_handler(proj)
        implicit none
        type(projection_handler),intent(inout) :: proj

        if (.not.allocated(proj%varnames)) allocate(proj%varnames(1:proj%nvars))
    end subroutine allocate_projection_handler

    subroutine projection_hydro(repository,cam,bulk_velocity,proj)
        implicit none
        character(128),intent(in) :: repository
        type(camera),intent(in) :: cam
        type(vector),intent(in) :: bulk_velocity
        type(projection_handler),intent(inout) :: proj

        type(hydroID) :: varIDs
        type(amr_info) :: amr
        type(sim_info) :: sim
        type(region) :: bbox
        real(dbl),dimension(:,:),allocatable :: toto
        integer,dimension(1:2) :: n_sample
        character(128) :: nomfich
        integer :: i,j
        real(dbl) :: xx,yy

        call read_hydrofile_descriptor(repository,varIDs)

        call init_amr_read(repository,amr,sim)
        amr%lmax = min(get_required_resolution(cam),amr%nlevelmax)
        write(*,*)'Maximum resolution level: ',amr%nlevelmax
        write(*,*)'Using: ',amr%lmax
        call get_bounding_box(cam,bbox)
        bbox%name = 'cube'
        bbox%bulk_velocity = bulk_velocity
        bbox%criteria_name = 'd_euclid'
        call get_cpu_map(bbox,amr)
        call get_map_box(cam,bbox)

        ! Perform projections
        call project_cells(repository,amr,bbox,varIDs,cam,proj)
        
    end subroutine projection_hydro

    subroutine project_cells(repository,amr,bbox,varIDs,cam,proj)
        implicit none
        character(128),intent(in) :: repository
        type(amr_info),intent(inout) :: amr
        type(region),intent(in) :: bbox
        type(hydroID),intent(in) :: varIDs
        type(camera),intent(in) :: cam
        type(projection_handler),intent(inout) :: proj

        logical :: ok_cell
        integer :: i,j,k
        integer :: ipos,icpu,ilevel,ind,idim,iidim,ivar
        integer :: ix,iy,iz,ngrida,nx_full,ny_full,nz_full
        integer :: imin,imax,jmin,jmax
        integer :: nvarh
        integer :: roterr
        character(5) :: nchar,ncharcpu
        character(128) :: nomfich
        real(dbl) :: distance,dx
        type(vector) :: xtemp,vtemp
        integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
        real(dbl),dimension(1:8,1:3) :: xc
        real(dbl),dimension(1:3,1:3) :: trans_matrix
        real(dbl),dimension(:,:),allocatable :: xg,x
        real(dbl),dimension(:,:,:),allocatable :: var
        integer,dimension(:,:),allocatable :: son
        logical,dimension(:),allocatable :: ref
        real(dbl) :: rho,map,weight
        real(dbl) :: xmin,ymin
        integer :: ndom
        integer,dimension(1:2) :: n_sample
        integer :: ncells
        

        type(level),dimension(1:100) :: grid

        ! Check that between the required variables are hydro variables



        ncells = 0

        ! Compute hierarchy
        do ilevel=1,amr%lmax
            nx_full = 2**ilevel
            ny_full = 2**ilevel
            imin = int(0D0*dble(nx_full))+1
            imax = int((bbox%xmax-bbox%xmin)*dble(nx_full))+1
            jmin = int(0D0*dble(ny_full))+1
            jmax = int((bbox%ymax-bbox%ymin)*dble(ny_full))+1
            allocate(grid(ilevel)%map(1:proj%nvars,imin:imax,jmin:jmax))
            allocate(grid(ilevel)%rho(imin:imax,jmin:jmax))
            grid(ilevel)%map(:,:,:) = 0D0
            grid(ilevel)%rho(:,:) = 0D0
            grid(ilevel)%imin = imin
            grid(ilevel)%imax = imax
            grid(ilevel)%jmin = jmin
            grid(ilevel)%jmax = jmax
        end do

        call los_transformation(cam,trans_matrix)


        allocate(ngridfile(1:amr%ncpu+amr%nboundary,1:amr%nlevelmax))
        allocate(ngridlevel(1:amr%ncpu,1:amr%nlevelmax))
        if(amr%nboundary>0)allocate(ngridbound(1:amr%nboundary,1:amr%nlevelmax))

        ipos=INDEX(repository,'output_')
        nchar=repository(ipos+7:ipos+13)
        ! Loop over processor files
        cpuloop: do k=1,amr%ncpu_read
            icpu = amr%cpu_list(k)
            call title(icpu,ncharcpu)

            ! Open AMR file and skip header
            nomfich = TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
            open(unit=10,file=nomfich,status='old',form='unformatted')
            write(*,*)'Processing file '//TRIM(nomfich)
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
                if(ngrida>0)then
                    allocate(xg(1:ngrida,1:amr%ndim))
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
                        do iidim=1,amr%ndim
                            if(j.eq.icpu)then
                                read(10)xg(:,iidim)
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

                !Compute map
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

                        ! Project positions onto the camera frame
                        call project_points(cam,ngrida,x)
                        ngridaloop: do i=1,ngrida
                            ! Check if cell is inside the desired region
                            distance = 0D0
                            ! write(*,*)x(i,:)
                            ! call checkifinside(x(i,:),bbox,ok_cell,distance)
                            ok_cell = .true.
                            ok_cell= ok_cell.and..not.ref(i)
                            ! write(*,*)'x:',x(i,:)
                            if (ok_cell) then
                                ix = int((x(i,1)+0.5*(bbox%xmax-bbox%xmin))*dble(nx_full)) + 1
                                iy = int((x(i,2)+0.5*(bbox%ymax-bbox%ymin))*dble(ny_full)) + 1
                                ! write(*,*)'ix,iy:',ix,iy
                                weight = (min(x(i,3)+dx/2.,bbox%zmax)-max(x(i,3)-dx/2.,bbox%zmin))/dx
                                weight = min(1.0d0,max(weight,0.0d0))
                                if( ix>=grid(ilevel)%imin.and.&
                                    & iy>=grid(ilevel)%jmin.and.&
                                    & ix<=grid(ilevel)%imax.and.&
                                    & iy<=grid(ilevel)%jmax) then
                                    ! write(*,*)'Cell is inside'
                                    xtemp = x(i,:)
                                    vtemp = var(i,ind,varIDs%vx:varIDs%vz)
                                    call rotate_vector(vtemp,trans_matrix)
                                    var(i,ind,varIDs%vx:varIDs%vz) = vtemp
                                    
                                    call getvarvalue(varIDs,bbox,dx,xtemp,var(i,ind,:),proj%weightvar,rho)
                                    grid(ilevel)%rho(ix,iy)=grid(ilevel)%rho(ix,iy)+rho*dx*weight/(bbox%zmax-bbox%zmin)

                                    projvarloop: do ivar=1,proj%nvars
                                        call getvarvalue(varIDs,bbox,dx,xtemp,var(i,ind,:),proj%varnames(ivar),map)
                                        ! if (TRIM(proj%varnames(ivar)).eq.proj%weightvar) map = map**2
                                        grid(ilevel)%map(ivar,ix,iy)=grid(ilevel)%map(ivar,ix,iy)+map*rho*dx*weight/(bbox%zmax-bbox%zmin)
                                    end do projvarloop
                                    
                                    
                                    ncells = ncells + 1
                                endif

                            endif
                        end do ngridaloop
                    end do cellloop
                    deallocate(xg,son,var,ref,x)
                endif
            end do levelloop
        end do cpuloop
        write(*,*)'ncells:',ncells
        ! Upload to maximum level (lmax)
        nx_full = 2**amr%lmax
        ny_full = 2**amr%lmax
        imin = int(0D0*dble(nx_full))+1
        imax = int((bbox%xmax-bbox%xmin)*dble(nx_full))
        jmin = int(0D0*dble(ny_full))+1
        jmax = int((bbox%ymax-bbox%ymin)*dble(ny_full))
        xloop: do ix = imin,imax
            xmin = ((ix-0.5)/2**amr%lmax)
            yloop: do iy=jmin,jmax
                ymin=((iy-0.5)/2**amr%lmax)
                ilevelloop: do ilevel=1,amr%lmax-1
                    ndom = 2**ilevel
                    i = int(xmin*ndom)+1
                    j = int(ymin*ndom)+1
                    
                    ! Smoothing: each cell contributes to its pixel and the 4
                    ! inmediate ones to it
                    projvarlooplmax: do ivar=1,proj%nvars
                        if (.false.) then
                        ! if (ix<imax.and.iy<jmax.and.&
                        !     &ix>imin.and.iy>imin) then
                            grid(amr%lmax)%map(ivar,ix,iy)=grid(amr%lmax)%map(ivar,ix,iy) + &
                                                        & grid(ilevel)%map(ivar,i,j)
                            grid(amr%lmax)%rho(ix,iy)=grid(amr%lmax)%rho(ix,iy) + &
                                                        & grid(ilevel)%rho(i,j)
                            grid(amr%lmax)%map(ivar,ix-1,iy)=grid(amr%lmax)%map(ivar,ix-1,iy) + &
                                                        & grid(ilevel)%map(ivar,i,j)
                            grid(amr%lmax)%rho(ix-1,iy)=grid(amr%lmax)%rho(ix-1,iy) + &
                                                        & grid(ilevel)%rho(i,j)
                            grid(amr%lmax)%map(ivar,ix+1,iy)=grid(amr%lmax)%map(ivar,ix+1,iy) + &
                                                        & grid(ilevel)%map(ivar,i,j)
                            grid(amr%lmax)%rho(ix+1,iy)=grid(amr%lmax)%rho(ix+1,iy) + &
                                                        & grid(ilevel)%rho(i,j)
                            grid(amr%lmax)%map(ivar,ix,iy-1)=grid(amr%lmax)%map(ivar,ix,iy-1) + &
                                                        & grid(ilevel)%map(ivar,i,j)
                            grid(amr%lmax)%rho(ix,iy-1)=grid(amr%lmax)%rho(ix,iy-1) + &
                                                        & grid(ilevel)%rho(i,j)
                            grid(amr%lmax)%map(ivar,ix,iy+1)=grid(amr%lmax)%map(ivar,ix,iy+1) + &
                                                        & grid(ilevel)%map(ivar,i,j)
                            grid(amr%lmax)%rho(ix,iy+1)=grid(amr%lmax)%rho(ix,iy+1) + &
                                                        & grid(ilevel)%rho(i,j)
                        else
                            grid(amr%lmax)%map(ivar,ix,iy)=grid(amr%lmax)%map(ivar,ix,iy) + &
                                                        & grid(ilevel)%map(ivar,i,j)
                            grid(amr%lmax)%rho(ix,iy)=grid(amr%lmax)%rho(ix,iy) + &
                                                        & grid(ilevel)%rho(i,j)
                        endif
                    end do projvarlooplmax
                end do ilevelloop
           end do yloop
        end do xloop

        call get_map_size(cam,n_sample)
        allocate(proj%toto(1:proj%nvars,0:n_sample(1),0:n_sample(2)))
        do i=0,n_sample(1)
            ix = int(dble(i)/dble(n_sample(1))*dble(imax-imin+1))+imin
            ix = min(ix,imax)
            do j=0,n_sample(2)
                iy = int(dble(j)/dble(n_sample(2))*dble(jmax-jmin+1))+jmin
                iy = min(iy,jmax)
                projvarlooptoto: do ivar=1,proj%nvars
                    proj%toto(ivar,i,j)=grid(amr%lmax)%map(ivar,ix,iy)/grid(amr%lmax)%rho(ix,iy)
                end do projvarlooptoto
            end do
         end do
        
    end subroutine project_cells

    subroutine projection_parts(repository,cam,bulk_velocity,proj)
        implicit none
        character(128),intent(in) :: repository
        type(camera),intent(in) :: cam
        type(vector),intent(in) :: bulk_velocity
        type(projection_handler),intent(inout) :: proj

        type(amr_info) :: amr
        type(sim_info) :: sim
        type(region) :: bbox

        call init_amr_read(repository,amr,sim)
        amr%lmax = min(get_required_resolution(cam),amr%nlevelmax)
        write(*,*)'Maximum resolution level: ',amr%nlevelmax
        write(*,*)'Using: ',amr%lmax
        call get_bounding_box(cam,bbox)
        bbox%name = 'cube'
        bbox%bulk_velocity = bulk_velocity
        bbox%criteria_name = 'd_euclid'
        call get_cpu_map(bbox,amr)
        call get_map_box(cam,bbox)

        call project_particles(repository,amr,sim,bbox,cam,proj)
    end subroutine projection_parts

    subroutine project_particles(repository,amr,sim,bbox,cam,proj)
        use cosmology
        implicit none
        character(128),intent(in) :: repository
        type(amr_info),intent(in) :: amr
        type(sim_info),intent(inout) :: sim
        type(region),intent(in) :: bbox
        type(camera),intent(in) :: cam
        type(projection_handler),intent(inout) :: proj

        logical :: ok_part
        integer :: i,j,k
        integer :: ipos,icpu,ix,iy,ixp1,iyp1,ivar
        integer :: npart,npart2,nstar,ncpu2,ndim2
        real(dbl) :: weight,distance,mapvalue
        real(dbl) :: dx,dy,ddx,ddy,dex,dey
        real(dbl),dimension(1:3,1:3) :: trans_matrix
        character(5) :: nchar,ncharcpu
        character(128) :: nomfich
        type(vector) :: xtemp,vtemp
        integer,dimension(:),allocatable :: id
        integer,dimension(1:2) :: n_map
        real(dbl),dimension(:),allocatable :: m,age,met,imass
        real(dbl),dimension(:,:),allocatable :: x,v

        ! Compute transformation matrix for camera LOS
        call los_transformation(cam,trans_matrix)
        
        ! Get camera resolution
        call get_map_size(cam,n_map)
        dx = (bbox%xmax-bbox%xmin)/n_map(1)
        dy = (bbox%ymax-bbox%ymin)/n_map(2)

        ! Allocate toto

        allocate(proj%toto(1:proj%nvars,0:n_map(1)-1,0:n_map(2)-1))

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
            write(*,*)'Processing file '//TRIM(nomfich)
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

        ! Compute projected variables using CIC smoothing
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
                read(1)age
                read(1)met
                read(1)imass
            endif
            close(1)

            ! Project positions onto the camera frame
            call project_points(cam,npart2,x)

            ! Project variables into map for particles
            ! of interest in the region
            partloop: do i=1,npart2
                weight = 1D0
                distance = 0D0
                mapvalue = 0D0
                call checkifinside(x(i,:),bbox,ok_part,distance)
                if (ok_part) then
                    xtemp = x(i,:)
                    vtemp = v(i,:)
                    call rotate_vector(vtemp,trans_matrix)
                    if (nstar>0) then
                        call getpartvalue(sim,bbox,xtemp,vtemp,id(i),m(i),age(i),met(i),imass(i),proj%weightvar,weight)
                    else
                        call getpartvalue(sim,bbox,xtemp,vtemp,0,m(i),0D0,0D0,0D0,proj%weightvar,weight)
                    endif

                    projvarloop: do ivar=1,proj%nvars
                        if (nstar>0) then
                            call getpartvalue(sim,bbox,xtemp,vtemp,id(i),m(i),age(i),met(i),imass(i),proj%varnames(ivar),mapvalue)
                        else
                            call getpartvalue(sim,bbox,xtemp,vtemp,0,m(i),0D0,0D0,0D0,proj%varnames(ivar),mapvalue)
                        endif
                        if (weight.ne.0D0) then
                            ! TODO: Properly understand WOH is going on here
                            ddx = (x(i,1)-bbox%xmin)/dx
                            ddy = (x(i,2)-bbox%ymin)/dy
                            ix = int(ddx)
                            iy = int(ddy)
                            ddx = ddx - dble(ix)
                            ddy = ddy - dble(iy)
                            dex = 1D0 - ddx
                            dey = 1D0 - ddy
                            ixp1 = ix + 1
                            iyp1 = iy + 1

                            if (ix>=0.and.ix<(n_map(1)-1).and.&
                                &iy>=0.and.iy<(n_map(2)-1).and.&
                                &ddx>0.and.ddy>0) then
                                proj%toto(ivar,ix  ,iy  ) = proj%toto(ivar,ix  ,iy  ) + mapvalue*dex*dey*weight
                                proj%toto(ivar,ix  ,iyp1) = proj%toto(ivar,ix  ,iyp1) + mapvalue*dex*ddy*weight
                                proj%toto(ivar,ixp1,iy  ) = proj%toto(ivar,ixp1,iy  ) + mapvalue*ddx*dey*weight
                                proj%toto(ivar,ixp1,iyp1) = proj%toto(ivar,ixp1,iyp1) + mapvalue*ddx*ddy*weight
                            endif
                        endif
                    end do projvarloop
                endif
            end do partloop
            deallocate(m,x,v)
            if (nstar>0)deallocate(id,age,met,imass)
        end do cpuloop
    end subroutine project_particles

end module maps