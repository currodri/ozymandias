module obs_instruments
    use local
    use vectors
    use coordinate_systems
    use geometrical_regions
    use filtering

    type camera
        type(vector) :: centre,los_axis,up_vector
        type(vector) :: region_axis, region_velocity
        real(dbl),dimension(1:2) :: region_size
        real(dbl) :: distance,far_cut_depth
        integer :: map_max_size=1024
        integer :: nfilter,nsubs
        integer :: lmin,lmax
        type(filter),dimension(:),allocatable :: filters
        type(region),dimension(:),allocatable :: subs
    end type camera

    type(vector),private :: x_axis,y_axis,z_axis

    contains

    real function log2(x)
        implicit none
        real(dbl), intent(in) :: x

        log2 = log(x) / log(2D0)
    end function

    type(camera) function init_camera(centre,los_axis,up_vector,region_size,region_axis,&
                                    & region_velocity,distance,far_cut_depth,map_max_size,&
                                    & nfilter,nsubs)
        implicit none
        type(vector),intent(in) :: centre,los_axis,up_vector,region_axis,region_velocity
        real(dbl),dimension(1:2),intent(in) :: region_size
        real(dbl),intent(in) :: distance,far_cut_depth
        integer,intent(in) :: map_max_size
        integer,intent(in) :: nfilter
        integer, intent(in) :: nsubs

        init_camera%centre = centre
        
        init_camera%los_axis = los_axis
        init_camera%up_vector = up_vector
        init_camera%region_axis = region_axis

        init_camera%region_velocity = region_velocity

        init_camera%region_size = region_size
        init_camera%distance = distance
        init_camera%far_cut_depth = far_cut_depth
        init_camera%map_max_size = map_max_size

        init_camera%nfilter = nfilter
        if (.not.allocated(init_camera%filters)) allocate(init_camera%filters(nfilter))

        init_camera%nsubs = nsubs
        if (.not.allocated(init_camera%subs).and.nsubs>0) allocate(init_camera%subs(nsubs))
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
        box%axis = cam%region_axis
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
    use utils
    use vectors
    use rotations
    use io_ramses
    use geometrical_regions
    use obs_instruments

    type projection_handler
        character(128) :: pov
        integer :: nvars,nfilter
        character(128),dimension(:),allocatable :: varnames
        character(128) :: weightvar
        real(dbl),dimension(:,:,:,:),allocatable :: toto
    end type projection_handler

    contains

    subroutine allocate_projection_handler(proj)
        implicit none
        type(projection_handler),intent(inout) :: proj

        if (.not.allocated(proj%varnames)) allocate(proj%varnames(1:proj%nvars))
    end subroutine allocate_projection_handler

    subroutine projection_hydro(repository,cam,use_neigh,proj,lmax,lmin)
        implicit none
        character(128),intent(in) :: repository
        type(camera),intent(inout) :: cam
        logical,intent(in) :: use_neigh
        type(projection_handler),intent(inout) :: proj
        integer,intent(in),optional :: lmax,lmin

        type(region) :: bbox

        call read_hydrofile_descriptor(repository)

        call init_amr_read(repository)
        amr%lmax = min(get_required_resolution(cam),amr%nlevelmax)
        write(*,*)'Maximum resolution level: ',amr%nlevelmax
        write(*,*)'Using: ',amr%lmax
        cam%lmin = 1;cam%lmax = amr%lmax
        if(present(lmin)) cam%lmin = max(1,min(lmin,amr%lmax))
        if(present(lmax)) cam%lmax = min(amr%lmax,max(lmax,1))
        write(*,*)'Camera using lmin, lmax: ',cam%lmin,cam%lmax
        call get_bounding_box(cam,bbox)
        bbox%name = 'cube'
        bbox%bulk_velocity = cam%region_velocity
        bbox%criteria_name = 'd_euclid'
        call get_cpu_map(bbox)
        write(*,*)'ncpu: ',amr%ncpu_read
        call get_map_box(cam,bbox)

        if (cam%nsubs>0)write(*,*)'Excluding substructure: ',cam%nsubs

        ! Perform projections
        if (use_neigh) then
            write(*,*)'Loading neighbours...'
            call project_cells_neigh(repository,bbox,cam,proj)
        else
            call project_cells(repository,bbox,cam,proj)
        end if

        contains

        subroutine project_cells(repository,bbox,cam,proj)
            implicit none
            character(128),intent(in) :: repository
            type(region),intent(in) :: bbox
            type(camera),intent(in) :: cam
            type(projection_handler),intent(inout) :: proj

            logical :: ok_cell,ok_filter,ok_sub,read_gravity
            integer :: i,j,k
            integer :: ipos,icpu,ilevel,ind,idim,iidim,ivar,ifilt,isub
            integer :: ix,iy,iz,ngrida,nx_full,ny_full,nz_full
            integer :: imin,imax,jmin,jmax
            integer :: nvarh
            integer :: roterr
            character(5) :: nchar,ncharcpu
            character(128) :: nomfich
            real(dbl) :: distance,dx
            type(vector) :: xtemp,vtemp,gtemp
            integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
            real(dbl),dimension(1:8,1:3) :: xc
            real(dbl),dimension(1:3,1:3) :: trans_matrix
            real(dbl),dimension(:,:),allocatable :: xg,x,xorig
            real(dbl),dimension(:,:,:),allocatable :: var,grav_var
            real(dbl),dimension(:,:),allocatable :: tempvar
            real(dbl),dimension(:,:),allocatable :: tempgrav_var
            integer,dimension(:,:),allocatable :: son
            integer,dimension(:),allocatable :: tempson
            logical,dimension(:),allocatable :: ref
            real(dbl) :: rho,map,weight
            real(dbl) :: xmin,ymin
            integer :: ndom
            integer,dimension(1:2) :: n_sample
            integer :: ncells
            

            type(level),dimension(1:100) :: grid


            ncells = 0

            ! Check whether we need to read the gravity files
            read_gravity = .false.
            do ivar=1,proj%nvars
                if (proj%varnames(ivar)(1:4) .eq. 'grav') then
                    read_gravity = .true.
                    write(*,*)'Reading gravity files...'
                    exit
                endif
            end do

            ! Compute hierarchy
            do ilevel=1,amr%lmax
                nx_full = 2**ilevel
                ny_full = 2**ilevel
                imin = int(0D0*dble(nx_full))+1
                imax = int((bbox%xmax-bbox%xmin)*dble(nx_full))+1
                jmin = int(0D0*dble(ny_full))+1
                jmax = int((bbox%ymax-bbox%ymin)*dble(ny_full))+1
                allocate(grid(ilevel)%cube(1:cam%nfilter,1:proj%nvars,imin:imax,jmin:jmax))
                allocate(grid(ilevel)%map(1:cam%nfilter,imin:imax,jmin:jmax))
                grid(ilevel)%cube(:,:,:,:) = 0D0
                grid(ilevel)%map(:,:,:) = 0D0
                grid(ilevel)%imin = imin
                grid(ilevel)%imax = imax
                grid(ilevel)%jmin = jmin
                grid(ilevel)%jmax = jmax
            end do

            trans_matrix = 0D0
            call new_z_coordinates(bbox%axis,trans_matrix,roterr)
            if (roterr.eq.1) then
                write(*,*) 'Incorrect CS transformation!'
                stop
            endif


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
                !write(*,*)'Processing file '//TRIM(nomfich)
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
                        if(read_gravity) allocate(grav_var(1:ngrida,1:amr%twotondim,1:4))
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
                            xorig = x
                            call project_points(cam,ngrida,x)
                            ngridaloop: do i=1,ngrida
                                ! Check if cell is inside the desired region
                                distance = 0D0
                                xtemp = x(i,:)
                                xtemp = xtemp - bbox%centre
                                x(i,:) = xtemp
                                call checkifinside(x(i,:),bbox,ok_cell,distance)
                                xtemp = xtemp + bbox%centre
                                x(i,:) = xtemp
                                ! TODO: This is a quick fix for the cells in the limits of the box
                                ! for axis-aligned projections, but should be done in a better way
                                ! using proper weights
                                if (trim(proj%pov)=='x'.or.trim(proj%pov)=='y'.or.trim(proj%pov)=='z') then
                                    ok_cell  = (bbox%xmin <= x(i,1)+dx/2.and.x(i,1)-dx/2 <= bbox%xmax.and.&
                                                &bbox%ymin <= x(i,2)+dx/2.and.x(i,2)-dx/2 <= bbox%ymax.and.&
                                                &bbox%zmin <= x(i,3)+dx/2.and.x(i,3)-dx/2 <= bbox%zmax)
                                end if
                                ok_cell= ok_cell.and..not.ref(i)
                                ! If we are avoiding substructure, check whether we are safe
                                if (cam%nsubs>0) then
                                    ok_sub = .true.
                                    do isub=1,cam%nsubs
                                        ok_sub = ok_sub .and. filter_sub(cam%subs(isub),xorig(i,:))
                                    end do
                                    ok_cell = ok_cell .and. ok_sub
                                end if
                                
                                if (ok_cell) then
                                    ix = int((x(i,1)+0.5*(bbox%xmax-bbox%xmin))*dble(nx_full)) + 1
                                    iy = int((x(i,2)+0.5*(bbox%ymax-bbox%ymin))*dble(ny_full)) + 1
                                    weight = (min(x(i,3)+dx/2.,bbox%zmax)-max(x(i,3)-dx/2.,bbox%zmin))/dx
                                    weight = min(1.0d0,max(weight,0.0d0))
                                    if( ix>=grid(ilevel)%imin.and.&
                                        & iy>=grid(ilevel)%jmin.and.&
                                        & ix<=grid(ilevel)%imax.and.&
                                        & iy<=grid(ilevel)%jmax) then
                                        ! Transform position to galaxy frame
                                        xtemp = xorig(i,:)
                                        xtemp = xtemp - bbox%centre
                                        call rotate_vector(xtemp,trans_matrix)
                                        ! Velocity transformed
                                        vtemp = var(i,ind,varIDs%vx:varIDs%vz)
                                        vtemp = vtemp - bbox%bulk_velocity
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
                                        filterloop: do ifilt=1,cam%nfilter
                                            if (read_gravity) then
                                                ok_filter = filter_cell(bbox,cam%filters(ifilt),xtemp,dx,tempvar,&
                                                                        &tempson,trans_matrix,tempgrav_var)
                                            else
                                                ok_filter = filter_cell(bbox,cam%filters(ifilt),xtemp,dx,tempvar,&
                                                                        &tempson,trans_matrix)
                                            end if
                                            ! Finally, get hydro data
                                            if (ok_filter) then
                                                weight = 1D0
                                                if (trim(proj%weightvar) /= 'counts') then
                                                    if (read_gravity) then 
                                                        call getvarvalue(bbox,dx,xtemp,tempvar,tempson,proj%weightvar,rho,trans_matrix,tempgrav_var)
                                                    else
                                                        call getvarvalue(bbox,dx,xtemp,tempvar,tempson,proj%weightvar,rho,trans_matrix)
                                                    end if
                                                    weight = MAX(rho*dx*weight/(bbox%zmax-bbox%zmin),0D0)
                                                end if
                                                grid(ilevel)%map(ifilt,ix,iy)=grid(ilevel)%map(ifilt,ix,iy)+weight

                                                projvarloop: do ivar=1,proj%nvars
                                                    if (read_gravity) then
                                                        call getvarvalue(bbox,dx,xtemp,tempvar,tempson,proj%varnames(ivar),map,trans_matrix,tempgrav_var)
                                                    else
                                                        call getvarvalue(bbox,dx,xtemp,tempvar,tempson,proj%varnames(ivar),map,trans_matrix)
                                                    end if
                                                    grid(ilevel)%cube(ifilt,ivar,ix,iy)=grid(ilevel)%cube(ifilt,ivar,ix,iy)+map*weight
                                                end do projvarloop
                                            end if
                                        end do filterloop
                                        ncells = ncells + 1
                                        deallocate(tempvar,tempson)
                                        if (read_gravity) deallocate(tempgrav_var)
                                    endif

                                endif
                            end do ngridaloop
                        end do cellloop
                        deallocate(xg,son,var,ref,x)
                        if (read_gravity) then
                            deallocate(grav_var)
                        end if
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
            filtlooplmax: do ifilt=1,cam%nfilter
                xloop: do ix = imin,imax
                    xmin = ((ix-0.5)/2**amr%lmax)
                    yloop: do iy=jmin,jmax
                        ymin=((iy-0.5)/2**amr%lmax)
                        ilevelloop: do ilevel=1,amr%lmax-1
                            ndom = 2**ilevel
                            i = int(xmin*ndom)+1
                            j = int(ymin*ndom)+1
                                projvarlooplmax: do ivar=1,proj%nvars
                                    grid(amr%lmax)%cube(ifilt,ivar,ix,iy)=grid(amr%lmax)%cube(ifilt,ivar,ix,iy) + &
                                                                & grid(ilevel)%cube(ifilt,ivar,i,j)
                                end do projvarlooplmax
                                grid(amr%lmax)%map(ifilt,ix,iy)=grid(amr%lmax)%map(ifilt,ix,iy) + &
                                                                & grid(ilevel)%map(ifilt,i,j)
                        end do ilevelloop
                end do yloop
                end do xloop
            end do filtlooplmax

            call get_map_size(cam,n_sample)
            proj%nfilter = cam%nfilter
            allocate(proj%toto(1:proj%nfilter,1:proj%nvars,0:n_sample(1),0:n_sample(2)))
            proj%toto = 0D0
            do i=0,n_sample(1)
                ix = int(dble(i)/dble(n_sample(1))*dble(imax-imin+1))+imin
                ix = min(ix,imax)
                do j=0,n_sample(2)
                    iy = int(dble(j)/dble(n_sample(2))*dble(jmax-jmin+1))+jmin
                    iy = min(iy,jmax)
                    filtlooptoto: do ifilt=1,proj%nfilter
                        projvarlooptoto: do ivar=1,proj%nvars
                            proj%toto(ifilt,ivar,i,j)=grid(amr%lmax)%cube(ifilt,ivar,ix,iy)/grid(amr%lmax)%map(ifilt,ix,iy)
                        end do projvarlooptoto
                    end do filtlooptoto
                end do
            end do        
        end subroutine project_cells

        subroutine project_cells_neigh(repository,bbox,cam,proj)
            implicit none
            character(128),intent(in) :: repository
            type(region),intent(in) :: bbox
            type(camera),intent(in) :: cam
            type(projection_handler),intent(inout) :: proj

            logical :: ok_cell,read_gravity,ok_filter,ok_cell_each,ok_sub
            integer :: i,j,k
            integer :: ipos,icpu,ilevel,ind,idim,iidim,ivar,iskip,inbor,ison,isub
            integer :: ix,iy,iz,ngrida,cumngrida,nx_full,ny_full,nz_full
            integer :: imin,imax,jmin,jmax
            integer :: nvarh
            integer :: roterr
            character(5) :: nchar,ncharcpu
            character(128) :: nomfich
            real(dbl) :: distance,dx
            type(vector) :: xtemp,vtemp,gtemp
            integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
            real(dbl),dimension(:),allocatable :: xxg,son_dens
            real(dbl),dimension(1:8,1:3) :: xc
            real(dbl),dimension(1:3,1:3) :: trans_matrix
            real(dbl),dimension(:,:),allocatable :: x,xorig
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
            real(dbl) :: rho,map,weight
            real(dbl) :: xmin,ymin
            integer :: ndom
            integer,dimension(1:2) :: n_sample
            integer :: ncells
            integer :: ngrid_current
            integer ::test,idebug
            integer :: ifilt
            

            type(level),dimension(1:100) :: grid

            ncells = 0
            test = 0
            idebug = 0

            ! Check whether we need to read the gravity files
            read_gravity = .false.
            do ivar=1,proj%nvars
                if (proj%varnames(ivar)(1:4) .eq. 'grav' .or.&
                & trim(proj%varnames(ivar)) .eq. 'neighbour_accuracy') then
                    read_gravity = .true.
                    write(*,*)'Reading gravity files...'
                    exit
                endif
            end do

            ! Compute hierarchy
            do ilevel=1,amr%lmax
                nx_full = 2**ilevel
                ny_full = 2**ilevel
                imin = int(0D0*dble(nx_full))+1
                imax = int((bbox%xmax-bbox%xmin)*dble(nx_full))+1
                jmin = int(0D0*dble(ny_full))+1
                jmax = int((bbox%ymax-bbox%ymin)*dble(ny_full))+1
                allocate(grid(ilevel)%cube(1:cam%nfilter,1:proj%nvars,imin:imax,jmin:jmax))
                allocate(grid(ilevel)%map(1:cam%nfilter,imin:imax,jmin:jmax))
                grid(ilevel)%cube(:,:,:,:) = 0D0
                grid(ilevel)%map(:,:,:) = 0D0
                grid(ilevel)%imin = imin
                grid(ilevel)%imax = imax
                grid(ilevel)%jmin = jmin
                grid(ilevel)%jmax = jmax
                grid(ilevel)%ngrid = 0
            end do

            !call los_transformation(cam,trans_matrix)
            trans_matrix = 0D0
            call new_z_coordinates(bbox%axis,trans_matrix,roterr)
            if (roterr.eq.1) then
                write(*,*) 'Incorrect CS transformation!'
                stop
            endif


            allocate(ngridfile(1:amr%ncpu+amr%nboundary,1:amr%nlevelmax))
            allocate(ngridlevel(1:amr%ncpu,1:amr%nlevelmax))
            if(amr%nboundary>0)allocate(ngridbound(1:amr%nboundary,1:amr%nlevelmax))

            ipos=INDEX(repository,'output_')
            nchar=repository(ipos+7:ipos+13)
            ! Loop over processor files
            cpuloop: do k=1,amr%ncpu_read
                icpu = amr%cpu_list(k)
                call title(icpu,ncharcpu)

                allocate(nbor(1:amr%ngridmax,1:amr%twondim))
                allocate(son(1:amr%ncoarse+amr%twotondim*amr%ngridmax))
                nbor = 0
                son = 0

                ! Open AMR file and skip header
                nomfich = TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
                open(unit=10,file=nomfich,status='old',form='unformatted')
                !write(*,*)'Processing file '//TRIM(nomfich)
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
                cellpos = 0d0
                var = 0d0

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
                            do iidim=1,amr%ndim
                                read(10)xxg
                                do i=1,ngrida
                                    grid(ilevel)%xg(i,iidim) = xxg(i)
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
                levelloop2: do ilevel=cam%lmin,cam%lmax
                    ! Geometry
                    dx = 0.5**ilevel
                    nx_full = 2**ilevel
                    ny_full = 2**ilevel

                    ! Allocate work arrays
                    ngrida = grid(ilevel)%ngrid
                    if(ngrida>0)then
                        allocate(ind_cell(1:ngrida))
                        allocate(x  (1:ngrida,1:amr%ndim))
                        allocate(xorig(1:ngrida,1:amr%ndim))
                        allocate(ref(1:ngrida))
                    endif
                    !Compute map
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
                                        test = test + 1
                                    end if
                                    deallocate(son_dens)
                                end if
                            end do

                            ! Project positions onto the camera frame
                            xorig = x
                            call project_points(cam,ngrida,x)
                            ! print*,'son',icpu,ilevel,son(ind_cell)
                            ngridaloop: do i=1,ngrida
                                ! Check if cell is inside the desired region
                                distance = 0D0
                                xtemp = x(i,:)
                                xtemp = xtemp - bbox%centre
                                x(i,:) = xtemp
                                call checkifinside(x(i,:),bbox,ok_cell,distance)
                                xtemp = xtemp + bbox%centre
                                x(i,:) = xtemp
                                ! TODO: This is a quick fix for the cells in the limits of the box
                                ! for axis-aligned projections, but should be done in a better way
                                ! using proper weights
                                if (trim(proj%pov)=='x'.or.trim(proj%pov)=='y'.or.trim(proj%pov)=='z') then
                                    ok_cell  = (bbox%xmin <= x(i,1)+dx/2.and.x(i,1)-dx/2 <= bbox%xmax.and.&
                                                &bbox%ymin <= x(i,2)+dx/2.and.x(i,2)-dx/2 <= bbox%ymax.and.&
                                                &bbox%zmin <= x(i,3)+dx/2.and.x(i,3)-dx/2 <= bbox%zmax)
                                end if
                                ok_cell= ok_cell.and..not.ref(i)
                                ! If we are avoiding substructure, check whether we are safe
                                if (cam%nsubs>0) then
                                    ok_sub = .true.
                                    do isub=1,cam%nsubs
                                        ok_sub = ok_sub .and. filter_sub(cam%subs(isub),xorig(i,:))
                                    end do
                                    ok_cell = ok_cell .and. ok_sub
                                end if
                                if (ok_cell) then
                                    ix = int((x(i,1)+0.5*(bbox%xmax-bbox%xmin))*dble(nx_full)) + 1
                                    iy = int((x(i,2)+0.5*(bbox%ymax-bbox%ymin))*dble(ny_full)) + 1
                                    weight = (min(x(i,3)+dx/2.,bbox%zmax)-max(x(i,3)-dx/2.,bbox%zmin))/dx
                                    weight = min(1.0d0,max(weight,0.0d0))
                                    if( ix>=grid(ilevel)%imin.and.&
                                        & iy>=grid(ilevel)%jmin.and.&
                                        & ix<=grid(ilevel)%imax.and.&
                                        & iy<=grid(ilevel)%jmax) then
                                        ! print*,ilevel,ind_cell(i),son(ind_cell(i)),var(ind_cell(i),1)
                                        xtemp = xorig(i,:)
                                        xtemp = xtemp - bbox%centre
                                        call rotate_vector(xtemp,trans_matrix)

                                        ! Velocity transformed --> ONLY FOR CENTRAL CELL
                                        vtemp = var(ind_cell(i),varIDs%vx:varIDs%vz)
                                        vtemp = vtemp - bbox%bulk_velocity
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
                                        
                                        do inbor=1,amr%twondim
                                            tempvar(inbor,:) = var(ind_nbor(inbor),:)
                                            tempson(inbor)       = son(ind_nbor(inbor))
                                            if (read_gravity) tempgrav_var(inbor,:) = grav_var(ind_nbor(inbor),:)
                                        end do

                                        filterloop: do ifilt=1,cam%nfilter
                                            if (read_gravity) then
                                                ok_filter = filter_cell(bbox,cam%filters(ifilt),xtemp,dx,tempvar,&
                                                                        &tempson,trans_matrix,tempgrav_var)
                                            else
                                                ok_filter = filter_cell(bbox,cam%filters(ifilt),xtemp,dx,tempvar,&
                                                                        &tempson,trans_matrix)
                                            end if
                                            ! Finally, get hydro data
                                            if (ok_filter) then
                                                weight = 1D0
                                                if (trim(proj%weightvar) /= 'counts') then
                                                    if (read_gravity) then 
                                                        call getvarvalue(bbox,dx,xtemp,tempvar,tempson,proj%weightvar,rho,trans_matrix,tempgrav_var)
                                                    else
                                                        call getvarvalue(bbox,dx,xtemp,tempvar,tempson,proj%weightvar,rho,trans_matrix)
                                                    end if
                                                    weight = rho*dx*weight/(bbox%zmax-bbox%zmin)
                                                end if
                                                grid(ilevel)%map(ifilt,ix,iy)=grid(ilevel)%map(ifilt,ix,iy)+weight

                                                projvarloop: do ivar=1,proj%nvars
                                                    if (read_gravity) then
                                                        call getvarvalue(bbox,dx,xtemp,tempvar,tempson,proj%varnames(ivar),map,trans_matrix,tempgrav_var)
                                                    else
                                                        call getvarvalue(bbox,dx,xtemp,tempvar,tempson,proj%varnames(ivar),map,trans_matrix)
                                                    end if
                                                    grid(ilevel)%cube(ifilt,ivar,ix,iy)=grid(ilevel)%cube(ifilt,ivar,ix,iy)+map*weight
                                                end do projvarloop
                                            end if
                                        end do filterloop
                                        
                                        
                                        ncells = ncells + 1
                                        deallocate(tempvar,tempson)
                                        if (read_gravity) deallocate(tempgrav_var)
                                    endif
                                endif
                            end do ngridaloop
                        end do cellloop
                        deallocate(ref,x,xorig,ind_cell)
                    endif
                end do levelloop2
                deallocate(nbor,son,var,cellpos)
                if (read_gravity) then
                    deallocate(grav_var)
                end if
            end do cpuloop
            write(*,*)'ncells: ',ncells
            write(*,*)'de-refine: ',test
            ! Upload to maximum level (lmax)
            nx_full = 2**amr%lmax
            ny_full = 2**amr%lmax
            imin = int(0D0*dble(nx_full))+1
            imax = int((bbox%xmax-bbox%xmin)*dble(nx_full))
            jmin = int(0D0*dble(ny_full))+1
            jmax = int((bbox%ymax-bbox%ymin)*dble(ny_full))
            filtlooplmax: do ifilt=1,cam%nfilter
                xloop: do ix = imin,imax
                    xmin = ((ix-0.5)/2**amr%lmax)
                    yloop: do iy=jmin,jmax
                        ymin=((iy-0.5)/2**amr%lmax)
                        ilevelloop: do ilevel=1,amr%lmax-1
                            ndom = 2**ilevel
                            i = int(xmin*ndom)+1
                            j = int(ymin*ndom)+1
                                projvarlooplmax: do ivar=1,proj%nvars
                                    grid(amr%lmax)%cube(ifilt,ivar,ix,iy)=grid(amr%lmax)%cube(ifilt,ivar,ix,iy) + &
                                                                & grid(ilevel)%cube(ifilt,ivar,i,j)
                                end do projvarlooplmax
                                grid(amr%lmax)%map(ifilt,ix,iy)=grid(amr%lmax)%map(ifilt,ix,iy) + &
                                                                & grid(ilevel)%map(ifilt,i,j)
                        end do ilevelloop
                end do yloop
                end do xloop
            end do filtlooplmax

            call get_map_size(cam,n_sample)
            proj%nfilter = cam%nfilter
            allocate(proj%toto(1:proj%nfilter,1:proj%nvars,0:n_sample(1),0:n_sample(2)))
            proj%toto = 0D0
            do i=0,n_sample(1)
                ix = int(dble(i)/dble(n_sample(1))*dble(imax-imin+1))+imin
                ix = min(ix,imax)
                do j=0,n_sample(2)
                    iy = int(dble(j)/dble(n_sample(2))*dble(jmax-jmin+1))+jmin
                    iy = min(iy,jmax)
                    filtlooptoto: do ifilt=1,proj%nfilter
                        projvarlooptoto: do ivar=1,proj%nvars
                            proj%toto(ifilt,ivar,i,j)=grid(amr%lmax)%cube(ifilt,ivar,ix,iy)/grid(amr%lmax)%map(ifilt,ix,iy)
                        end do projvarlooptoto
                    end do filtlooptoto
                end do
            end do        
        end subroutine project_cells_neigh

        
    end subroutine projection_hydro

    subroutine projection_parts(repository,cam,proj,tag_file,inverse_tag)
        implicit none
        character(128),intent(in) :: repository
        type(camera),intent(in) :: cam
        type(projection_handler),intent(inout) :: proj
        character(128),intent(in),optional :: tag_file
        logical,intent(in),optional :: inverse_tag

        type(region) :: bbox

        call init_amr_read(repository)
        amr%lmax = amr%nlevelmax !min(get_required_resolution(cam),amr%nlevelmax)
        write(*,*)'Maximum resolution level: ',amr%nlevelmax
        write(*,*)'Using: ',amr%lmax
        if (sim%dm .and. sim%hydro) call check_families(repository)
        call get_bounding_box(cam,bbox)
        bbox%name = 'cube'
        bbox%bulk_velocity = cam%region_velocity
        bbox%criteria_name = 'd_euclid'
        call get_cpu_map(bbox)
        call get_map_box(cam,bbox)
        if (cam%nsubs>0)write(*,*)'Excluding substructure: ',cam%nsubs
        if (present(tag_file)) then
            if (present(inverse_tag)) then
                call project_particles(repository,bbox,cam,proj,tag_file,inverse_tag)
            else
                call project_particles(repository,bbox,cam,proj,tag_file)
            endif
        else
            call project_particles(repository,bbox,cam,proj)
        endif
        write(*,*)minval(proj%toto),maxval(proj%toto)
    end subroutine projection_parts

    subroutine project_particles(repository,bbox,cam,proj,tag_file,inverse_tag)
#ifndef LONGINT
        use utils, only:quick_sort_irg,binarysearch_irg
#else
        use utils, only:quick_sort_ilg,binarysearch_ilg
#endif
        use cosmology
        implicit none
        character(128),intent(in) :: repository
        type(region),intent(in) :: bbox
        type(camera),intent(in) :: cam
        type(projection_handler),intent(inout) :: proj
        character(128),intent(in),optional :: tag_file
        logical,intent(in),optional :: inverse_tag

        logical :: ok_part,ok_tag,ok_filter,ok_sub
        integer :: i,j,k,itag,ifilt,isub
        integer :: ipos,icpu,ix,iy,ixp1,iyp1,ivar
        integer(irg) :: npart,npart2,nstar,ntag
        integer :: npartsub
        integer :: ncpu2,ndim2
        real(dbl) :: weight,distance,mapvalue
        real(dbl) :: dx,dy,ddx,ddy,dex,dey
        real(dbl),dimension(1:3,1:3) :: trans_matrix
        character(5) :: nchar,ncharcpu
        character(6) :: ptype
        character(128) :: nomfich
        type(vector) :: xtemp,vtemp,dcell
        type(particle) :: part
        integer,dimension(:),allocatable :: order
        integer,dimension(:),allocatable :: nparttoto
#ifndef LONGINT
        integer(irg),dimension(:),allocatable :: id,tag_id
#else
        integer(ilg),dimension(:),allocatable :: id,tag_id
#endif
#ifndef IMASS
        integer(1),dimension(:), allocatable :: part_tags
#endif
        integer,dimension(1:2) :: n_map
        real(dbl),dimension(:),allocatable :: m,age,met
        real(dbl),dimension(:),allocatable :: imass
        real(dbl),dimension(:,:),allocatable :: x,v

        npartsub = 0
#ifndef IMASS
        if (sim%eta_sn .eq. -1D0) then
            write(*,*)': eta_sn=-1 and not IMASS --> should set this up!'
            stop
        end if
#endif

        ! If tagged particles file exists, read and allocate array
        if (present(tag_file)) then
            open(unit=58,file=TRIM(tag_file),status='old',form='formatted')
            write(*,*)'Reading particle tags file '//TRIM(tag_file)
            read(58,'(I11)')ntag
            write(*,*)'Number of tagged particles in file: ',ntag
            if (allocated(tag_id)) then
                deallocate(tag_id)
                allocate(tag_id(1:ntag))
            else
                allocate(tag_id(1:ntag))
            endif
            do itag=1,ntag
                read(58,'(I11)')tag_id(itag)
            end do
            allocate(order(1:ntag))
            write(*,*)'Sorting list of particle ids for binary search...'
#ifndef LONGINT
            call quick_sort_irg(tag_id,order,ntag)
#else
            call quick_sort_ilg(tag_id,order,ntag)
#endif
            deallocate(order)
            close(58)
        endif
        
        ! Compute transformation matrix for camera LOS
        call los_transformation(cam,trans_matrix)
        
        ! Get camera resolution
        call get_map_size(cam,n_map)
        dx = (bbox%xmax-bbox%xmin)/n_map(1)
        dy = (bbox%ymax-bbox%ymin)/n_map(2)
        dcell = (/dx,dy,0D0/)

        ! Allocate toto
        proj%nfilter = cam%nfilter
        allocate(proj%toto(1:proj%nfilter,1:proj%nvars,0:n_map(1)-1,0:n_map(2)-1))

        proj%toto = 0D0

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
        allocate(nparttoto(1:proj%nfilter))
        nparttoto = 0
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
#ifndef IMASS
                allocate(part_tags(1:npart2))
#endif
            endif
            if (present(tag_file) .and. (.not. allocated(id))) allocate(id(1:npart2))
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
            read(1)id
            read(1) ! Skip level
            if (nstar>0) then
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
            elseif (present(tag_file) .and. nstar .eq. 0) then
                read(1)id
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
                elseif (present(tag_file)) then
                    part%id = id(i)
                    part%age = 0D0
                    part%met = 0D0
                    part%imass = 0D0
                else
                    part%id = 0
                    part%age = 0D0
                    part%met = 0D0
                    part%imass = 0D0
                endif
                xtemp = x(i,:)
                xtemp = xtemp - bbox%centre
                x(i,:) = xtemp
                call checkifinside(x(i,:),bbox,ok_part,distance)
                xtemp = xtemp + bbox%centre
                x(i,:) = xtemp

                ! If we are avoiding substructure, check whether we are safe
                if (cam%nsubs>0) then
                    ok_sub = .true.
                    xtemp = xtemp + bbox%centre
                    do isub=1,cam%nsubs
                        ok_sub = ok_sub .and. filter_sub(cam%subs(isub),(/xtemp%x,xtemp%y,xtemp%z/))
                    end do
                    if (.not.ok_sub)npartsub = npartsub + 1
                    ok_part = ok_part .and. ok_sub
                end if
                ! Check if tags are present for particles
                if (present(tag_file) .and. ok_part) then
                    ok_tag = .false.      
#ifndef LONGINT
                    call binarysearch_irg(ntag,tag_id,part%id,ok_tag)
#else
                    call binarysearch_ilg(ntag,tag_id,part%id,ok_tag)
#endif
                    if (present(inverse_tag) .and. inverse_tag .and. ok_tag) ok_tag = .false.
                    ok_part = ok_tag .and. ok_part
                endif
                if (ok_part) then
                    part%v = part%v - bbox%bulk_velocity
                    call rotate_vector(part%v,trans_matrix)
                    if (nstar>0) then
                        if (TRIM(proj%weightvar).eq.'star/cumulative'.or.&
                            &TRIM(proj%weightvar).eq.'dm/cumulative') then
                            weight = 1D0
                        else
                            call getparttype(part,ptype)
                            if (ptype.eq.'star') then
                                call getpartvalue(bbox,part,proj%weightvar,weight,dcell)
                            else
                                weight = 1D0
                            endif
                        endif
                    else
                        weight = 1D0
                        ! call getpartvalue(bbox,xtemp,vtemp,0,m(i),0D0,0D0,0D0,proj%weightvar,weight)
                    endif

                    projvarloop: do ivar=1,proj%nvars
                        call getpartvalue(bbox,part,proj%varnames(ivar),mapvalue,dcell)
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
                                filtlooptoto: do ifilt=1,proj%nfilter
                                    ok_filter = filter_particle(bbox,cam%filters(ifilt),part)
                                    if (ok_filter) then
                                        proj%toto(ifilt,ivar,ix  ,iy  ) = proj%toto(ifilt,ivar,ix  ,iy  ) + mapvalue*dex*dey*weight
                                        proj%toto(ifilt,ivar,ix  ,iyp1) = proj%toto(ifilt,ivar,ix  ,iyp1) + mapvalue*dex*ddy*weight
                                        proj%toto(ifilt,ivar,ixp1,iy  ) = proj%toto(ifilt,ivar,ixp1,iy  ) + mapvalue*ddx*dey*weight
                                        proj%toto(ifilt,ivar,ixp1,iyp1) = proj%toto(ifilt,ivar,ixp1,iyp1) + mapvalue*ddx*ddy*weight
                                        nparttoto(ifilt) = nparttoto(ifilt) + 1
                                    end if
                                end do filtlooptoto
                                
                            endif
                        endif
                    end do projvarloop
                endif
            end do partloop
            deallocate(m,x,v)
            if (allocated(id))deallocate(id)
            if (nstar>0)deallocate(age,met,imass)
#ifndef IMASS
            if (nstar>0)deallocate(part_tags)
#endif
        end do cpuloop
        write(*,*)'> nparttoto: ',nparttoto
        write(*,*)'> npartsub:  ',npartsub
        deallocate(nparttoto)
    end subroutine project_particles

    subroutine healpix_hydro(repository,cam,use_neigh,proj,nside,lmax,lmin)
        implicit none
        character(128),intent(in) :: repository
        type(camera),intent(inout) :: cam
        logical,intent(in) :: use_neigh
        type(projection_handler),intent(inout) :: proj
        integer,intent(in) :: nside
        integer,intent(in),optional :: lmax,lmin

        type(region) :: bsphere

        call read_hydrofile_descriptor(repository)

        call init_amr_read(repository)
        amr%lmax = amr%nlevelmax
        write(*,*)'Maximum resolution level: ',amr%nlevelmax
        write(*,*)'Using: ',amr%lmax
        cam%lmin = 1;cam%lmax = amr%lmax
        if(present(lmin)) cam%lmin = max(1,min(lmin,amr%lmax))
        if(present(lmax)) cam%lmax = min(amr%lmax,max(lmax,1))
        write(*,*)'Camera using lmin, lmax: ',cam%lmin,cam%lmax
        bsphere%name = 'sphere'
        bsphere%bulk_velocity = cam%region_velocity
        bsphere%criteria_name = 'd_euclid'
        bsphere%centre = cam%centre
        bsphere%rmin = cam%distance
        bsphere%rmax = cam%far_cut_depth
        bsphere%axis = cam%region_axis
        call get_cpu_map(bsphere)
        write(*,*)'ncpu: ',amr%ncpu_read

        if (cam%nsubs>0)write(*,*)'Excluding substructure: ',cam%nsubs

        ! Perform projections
        call project_cells_hpix

        contains

        subroutine project_cells_hpix
            use healpix_modules
            implicit none

            logical :: ok_cell,ok_filter,ok_sub
            integer :: i,j,k
            integer :: ipos,icpu,ilevel,ind,idim,iidim,ivar,iskip
            integer :: isub,inbor,ifilt
            integer :: ix,iy,iz,ngrida,ns,cumngrida
            integer :: imin,imax
            integer :: nvarh
            integer :: roterr
            character(5) :: nchar,ncharcpu
            character(128) :: nomfich
            real(dbl) :: distance,dx
            type(vector) :: xtemp,vtemp,los,x_axis
            integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
            real(dbl),dimension(1:3) :: xvec
            real(dbl),dimension(1:8,1:3) :: xc
            real(dbl),dimension(1:3,1:3) :: trans_matrix
            real(dbl),dimension(:,:),allocatable :: xg,x,xorig
            real(dbl),dimension(:,:,:),allocatable :: var
            real(dbl),dimension(:,:),allocatable :: tempvar
            integer,dimension(:,:),allocatable :: son
            integer,dimension(:),allocatable :: tempson
            logical,dimension(:),allocatable :: ref
            real(dbl),dimension(:,:),allocatable :: proj_rho
            integer,dimension(:),allocatable :: listpix,listpix_clean
            real(dbl) :: rho,map,weight,aperture
            real(dbl) :: xmin,ymin
            integer :: ndom
            integer :: ncells
            integer :: nlist
            
            type(basis) :: hpix_basis

            x_axis = (/1D0,0D0,0D0/)

            ncells = 0

            ns = nside2npix(nside)-1
            write(*,*)'Total number of healpix map: ',ns+1
            ! TODO: We should also include in here the filters
            proj%nfilter = cam%nfilter
            allocate(proj%toto(1:proj%nfilter,1:proj%nvars,1:1,0:ns))
            allocate(proj_rho(1:proj%nfilter,0:ns))
            allocate(listpix(0:ns))
            proj%toto = 0D0
            proj_rho = 0D0
            listpix = -1

            ! Get transformation matrix. Default: Region axis is z axis, LOS is x axis
            los = x_axis - (x_axis.DOT.bsphere%axis)*bsphere%axis
            los = los / magnitude(los)
            hpix_basis%u(1) = los
            hpix_basis%u(2) = bsphere%axis * los
            hpix_basis%u(2) = hpix_basis%u(2) / magnitude(hpix_basis%u(2))
            hpix_basis%u(3) = bsphere%axis
            
            trans_matrix = 0D0
            do i=1,3
                trans_matrix(i,:) = hpix_basis%u(i)
            end do

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
                ! write(*,*)'Processing file '//TRIM(nomfich)
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

                            xorig = x
                            ngridaloop: do i=1,ngrida
                                ! Check if cell is inside the desired region
                                distance = 0D0
                                xtemp = x(i,:)

                                ! Move to center of galaxy
                                x(i,:) = xtemp - bsphere%centre

                                ! Check if cell is inside the desired region
                                call checkifinside(x(i,:),bsphere,ok_cell,distance)
                                ok_cell= ok_cell.and.(.not.ref(i))
                                listpix = -1
                                nlist = 0

                                ! If we are avoiding substructure, check whether we are safe
                                if (cam%nsubs>0) then
                                    ok_sub = .true.
                                    do isub=1,cam%nsubs
                                        ok_sub = ok_sub .and. filter_sub(cam%subs(isub),xorig(i,:))
                                    end do
                                    ok_cell = ok_cell .and. ok_sub
                                end if

                                if (ok_cell) then
                                    xtemp = x(i,:)
                                    ! Rotate position such that we have cells in the galaxy frame
                                    call rotate_vector(xtemp,trans_matrix)
                                    x(i,:) = xtemp

                                    ! Get pixels to which the cell contributes
                                    aperture = datan(0.707*dx/distance)
                                    call query_disc(nside,x(i,:),aperture,listpix,nlist)
                                    listpix_clean = pack(listpix,listpix.ge.0)
                                    if (nlist.eq.0) then
                                        listpix = -1
                                        nlist = 0
                                        call query_disc(nside,x(i,:),aperture,listpix,nlist,nest=0,inclusive=1)
                                        listpix_clean = pack(listpix,listpix.ge.0)
                                    end if

                                    ! Compute a weight for cells that are partly outside region
                                    weight = (min(distance+dx/2.,bsphere%rmax)-max(distance-dx/2.,bsphere%rmin))/dx
                                    weight = min(1.0d0,max(weight,0.0d0))

                                    ! If the cell contributes to at least one pixel, project
                                    if(nlist>0) then
                                        ! Rotate velocity with respect to galaxy frame
                                        vtemp = var(i,ind,varIDs%vx:varIDs%vz)
                                        vtemp = vtemp - bsphere%bulk_velocity
                                        call rotate_vector(vtemp,trans_matrix)
                                        var(i,ind,varIDs%vx:varIDs%vz) = vtemp

                                        ! Get neighbours
                                        allocate(tempvar(0:amr%twondim,nvarh))
                                        allocate(tempson(0:amr%twondim))
                                        ! Just add central cell as we do not want neighbours
                                        tempvar(0,:) = var(i,ind,:)
                                        tempson(0)       = son(i,ind)
                                        tempvar(0,varIDs%vx:varIDs%vz) = vtemp

                                        ! Do loop over filters
                                        filterloop: do ifilt=1,cam%nfilter
                                            ok_filter = filter_cell(bsphere,cam%filters(ifilt),xtemp,dx,tempvar,&
                                                                        &tempson,trans_matrix)
                                            if (ok_filter) then
                                                ! Get weight
                                                call getvarvalue(bsphere,dx,xtemp,tempvar,tempson,proj%weightvar,rho)
                                                do j=1,size(listpix_clean)
                                                    ix = listpix_clean(j)
                                                    proj_rho(ifilt,ix) = proj_rho(ifilt,ix)+rho*dx*weight/(bsphere%rmax-bsphere%rmin)
                                                end do

                                                ! Get variable values
                                                projvarloop: do ivar=1,proj%nvars
                                                    call getvarvalue(bsphere,dx,xtemp,tempvar,tempson,proj%varnames(ivar),map)
                                                    do j=1,size(listpix_clean)
                                                        ix = listpix_clean(j)
                                                        proj%toto(ifilt,ivar,1,ix) = proj%toto(ifilt,ivar,1,ix)+map*rho*dx*weight/(bsphere%rmax-bsphere%rmin)
                                                    end do
                                                end do projvarloop
                                            end if
                                        end do filterloop
                                        
                                        ncells = ncells + 1
                                        deallocate(tempvar,tempson)
                                    endif

                                endif
                            end do ngridaloop
                        end do cellloop
                        deallocate(xg,son,var,ref,x)
                    endif
                end do levelloop
            end do cpuloop
            write(*,*)'ncells:',ncells

            ! Renormalise cells to compute weighted values
            filtlooptoto: do ifilt=1,proj%nfilter
                do i=0,ns
                    projvarlooptoto: do ivar=1,proj%nvars
                        proj%toto(ifilt,ivar,1,i) = proj%toto(ifilt,ivar,1,i)/proj_rho(ifilt,i)
                    end do projvarlooptoto
                end do
            end do filtlooptoto

        end subroutine project_cells_hpix

    end subroutine healpix_hydro
end module maps