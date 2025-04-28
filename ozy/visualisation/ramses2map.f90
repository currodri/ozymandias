module obs_instruments
    use local
    use dictionary_commons
    use vectors
    use coordinate_systems
    use geometrical_regions
    use filtering

    type camera
        type(vector) :: centre,los_axis,up_vector
        type(vector) :: region_axis, region_velocity
        real(dbl),dimension(1:2) :: region_size
        real(dbl),dimension(1:2) :: dx_pixel
        real(dbl) :: distance,far_cut_depth
        integer :: map_max_size=1024
        integer :: nsubs
        integer :: lmin,lmax
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
                                    & nsubs)
        implicit none
        type(vector),intent(in) :: centre,los_axis,up_vector,region_axis,region_velocity
        real(dbl),dimension(1:2),intent(in) :: region_size
        real(dbl),intent(in) :: distance,far_cut_depth
        integer,intent(in) :: map_max_size
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

        init_camera%dx_pixel(:) = dble(init_camera%region_size(:)) / init_camera%map_max_size

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
        else if (aspect_ratio < 1D0) then
            n_map(2) = cam%map_max_size
            n_map(1) = int(nint(cam%map_max_size*aspect_ratio))
        else
            n_map(2) = cam%map_max_size
            n_map(1) = cam%map_max_size
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
    use dictionary_commons
    use utils
    use vectors
    use rotations
    use io_ramses
    use geometrical_regions
    use obs_instruments
    use hydro_commons
    use part_commons

    type hydro_projection_handler
        character(128) :: pov
        integer :: nvars,nfilter,nwvars
        integer, dimension(1:2) :: n_sample
        character(128),dimension(:),allocatable :: varnames
        character(128),dimension(:),allocatable :: weightvars
        real(dbl),dimension(:,:,:,:,:),allocatable :: map
        real(dbl),dimension(:,:,:,:),allocatable :: weights
        type(hydro_var),dimension(:),allocatable :: vars
        type(hydro_var),dimension(:),allocatable :: wvars
        type(filter_hydro),dimension(:),allocatable :: filters
    end type hydro_projection_handler

    type part_projection_handler
        character(128) :: pov
        integer :: nvars,nfilter,nwvars
        integer, dimension(1:2) :: n_sample
        character(128),dimension(:),allocatable :: varnames
        character(128),dimension(:),allocatable :: weightvars
        real(dbl),dimension(:,:,:,:,:),allocatable :: map
        real(dbl),dimension(:,:,:,:),allocatable :: weights
        type(part_var),dimension(:),allocatable :: vars
        type(part_var),dimension(:),allocatable :: wvars
        type(filter_part),dimension(:),allocatable :: filters
    end type part_projection_handler

    ! Gaussian kernel configuration and normalisation
    integer, parameter::rmax_npix_kernel = 50 ! Max number of pixels for radius (for computational speed-up)
    integer, parameter::kdx_expand = 25
    integer, parameter::nexpand_gauss = 30
    integer, parameter::ngauss = 10
    real(sgl), dimension(0:ngauss)::gauss_norm_values = &
        & (/1.0, 1.16825, 2.28371, 2.49311, 2.50631, 2.50663,&
            2.50663,2.50663,2.50663,2.50663,2.50663/)

    contains

    subroutine allocate_hydro_projection_handler(proj)
        implicit none
        type(hydro_projection_handler),intent(inout) :: proj

        if (.not.allocated(proj%varnames)) allocate(proj%varnames(1:proj%nvars))
        if (.not.allocated(proj%weightvars)) allocate(proj%weightvars(1:proj%nwvars))
        if (.not.allocated(proj%vars)) allocate(proj%vars(1:proj%nvars))
        if (.not.allocated(proj%wvars)) allocate(proj%wvars(1:proj%nwvars))
        if (.not.allocated(proj%filters)) allocate(proj%filters(1:proj%nfilter))
    end subroutine allocate_hydro_projection_handler

    subroutine allocate_part_projection_handler(proj)
        implicit none
        type(part_projection_handler),intent(inout) :: proj

        if (.not.allocated(proj%varnames)) allocate(proj%varnames(1:proj%nvars))
        if (.not.allocated(proj%weightvars)) allocate(proj%weightvars(1:proj%nwvars))
        if (.not.allocated(proj%vars)) allocate(proj%vars(1:proj%nvars))
        if (.not.allocated(proj%wvars)) allocate(proj%wvars(1:proj%nwvars))
        if (.not.allocated(proj%filters)) allocate(proj%filters(1:proj%nfilter))
    end subroutine allocate_part_projection_handler

    subroutine get_pos_map(n_sample,grid,ix,iy,ii_map)
        implicit none

        integer, dimension(1:2),intent(in) :: n_sample
        type(level), intent(in) :: grid
        integer, intent(in) :: ix,iy
        integer, dimension(1:2), intent(inout) :: ii_map

        real(dbl) :: xconv,yconv

        ! Grid to map transformation factors
        xconv = dble(n_sample(1))/dble(grid%imax - grid%imin + 1)
        yconv = dble(n_sample(2))/dble(grid%jmax - grid%jmin + 1)

        ! Compute map pixel
        ii_map(1) = int(dble(ix-grid%imin)*xconv)
        ii_map(2) = int(dble(iy-grid%jmin)*yconv)
    end subroutine get_pos_map

    subroutine grid_projection(proj,cam,grid,nexp_factor)
        use constants, only: halfsqrt2
        implicit none
    
        type(hydro_projection_handler), intent(inout) :: proj
        type(camera), intent(in)                :: cam
        type(level), dimension(100), intent(in)  :: grid
        real(dbl), intent(in), optional :: nexp_factor

        integer     ::nexpand
        integer     ::dx_i
        integer     ::ix,iy,ivar,ilevel,ifilt,iweight
        
        real(dbl)   ::kdx, kdx2, dx, ksize
        real(dbl)   ::nexp_f
        real(dbl)   ::kernel_norm
        integer     ::ipix, jpix
        real(dbl)   ::ipix2, jpix2
        integer     ::ixh,iyh
        real(dbl)   ::rpixel2
        integer     ::ixmin, ixmax
        integer     ::jymin, jymax
        integer     ::kxmin, kxmax
        integer     ::kymin, kymax
        real(dbl)   ::rcut2
        real(dbl)   ::xconv, yconv
        ! real(dbl)   ::kfilter
        real(dbl), dimension(:,:), allocatable :: kfilter
        integer,dimension(1:2) :: xpix

        if (verbose) write(*,*)'Using gaussian smoothing with nexp_factor = ',nexp_factor

        ! Initialise map
        call get_map_size(cam,proj%n_sample)
        allocate(proj%map(1:proj%nfilter,1:proj%nvars,1:proj%nwvars,1:proj%n_sample(1),1:proj%n_sample(2)))
        allocate(proj%weights(1:proj%nfilter,1:proj%nwvars,1:proj%n_sample(1),1:proj%n_sample(2)))
        proj%map     = 0D0
        proj%weights = 0D0

        if (present(nexp_factor)) then
            nexp_f = nexp_factor
        else
            nexp_f = 1.0d0
        end if
        ! Loop over levels
        do ilevel=cam%lmin,cam%lmax
            if (.not. grid(ilevel)%active) cycle
            dx = 0.5**ilevel
            ksize = (dx/min(cam%dx_pixel(1), cam%dx_pixel(2)))
            if (ksize .lt. 1.0d0) then
                ksize = 0.5d0 * ksize
            end if
            xconv = dble(proj%n_sample(1))/dble(grid(ilevel)%imax - grid(ilevel)%imin + 1)
            yconv = dble(proj%n_sample(2))/dble(grid(ilevel)%jmax - grid(ilevel)%jmin + 1)
            if (ksize.lt.halfsqrt2) then
                ! Loop over projected pixels in ilevel
                do iy = grid(ilevel)%jmin, grid(ilevel)%jmax
                    xpix(2) = int(dble(iy-grid(ilevel)%jmin+1)*xconv)
                    if (xpix(2) .lt. 1) cycle
                    if (xpix(2) .gt. proj%n_sample(2)) cycle
                    do ix = grid(ilevel)%imin, grid(ilevel)%imax
                        xpix(1) = int(dble(ix-grid(ilevel)%imin+1)*yconv)
                        if (xpix(1) .lt. 1) cycle
                        if (xpix(1) .gt. proj%n_sample(1)) cycle
                        ! Get the index of cell in the map grid
                        ixh = int(xpix(1)); iyh = int(xpix(2))
                        ! Project to map
                        do ifilt = 1, proj%nfilter
                            do iweight = 1, proj%nwvars
                                proj%weights(ifilt,iweight,ixh,iyh) = proj%weights(ifilt,iweight,ixh,iyh) + grid(ilevel)%map(ifilt,iweight,ix,iy)
                                do ivar = 1, proj%nvars
                                    proj%map(ifilt,ivar,iweight,ixh,iyh) = proj%map(ifilt,ivar,iweight,ixh,iyh) + grid(ilevel)%cube(ifilt,ivar,iweight,ix,iy)
                                end do
                            end do
                        end do
                    end do
                end do
            else
                ! Compute kernel translation to map properties
                kdx = ksize
                kdx = max(kdx, 0.25d0) ! Kernels should not be kdx < dx_pixel/4
                kdx2 = kdx**2.0d0

                nexpand = nexp_f*min(int(kdx_expand*kdx), nexpand_gauss)
                if (nexpand .gt. ngauss) then
                    kernel_norm = gauss_norm_values(ngauss)
                else if (nexpand .ge. 0) then
                    kernel_norm = gauss_norm_values(nexpand)
                else
                    write (*, *) "ERROR: nexpand < 0 for gauss"
                    stop
                end if
                kernel_norm = 1.0d0/(kernel_norm*twopi*kdx)
                dx_i = min(max(int(nexpand*kdx), 1), rmax_npix_kernel)
                rcut2 = dble(dx_i)**2.0d0
                call create_kfilter
                ! Loop over projected pixels in ilevel
                do iy = grid(ilevel)%jmin, grid(ilevel)%jmax
                    xpix(2) = int(dble(iy-grid(ilevel)%jmin+1)*xconv)
                    if (xpix(2) .lt. 1) cycle
                    if (xpix(2) .gt. proj%n_sample(2)) cycle
                    do ix = grid(ilevel)%imin, grid(ilevel)%imax
                        xpix(1) = int(dble(ix-grid(ilevel)%imin+1)*yconv)
                        if (xpix(1) .lt. 1) cycle
                        if (xpix(1) .gt. proj%n_sample(1)) cycle
                        ! Compute kernel limits
                        ! print*,'Checking ', ix,iy
                        ixmin = int(xpix(1) - dx_i)
                        jymin = int(xpix(2) - dx_i)
                        ixmax = int(xpix(1) + dx_i)
                        jymax = int(xpix(2) + dx_i)
                        ! print*,xpix,dx_i,ixmin,ixmax,jymin,jymax
                        kxmin = 1 - min(0,ixmin)
                        kymin = 1 - min(0,jymin)

                        ixmin = max(ixmin, 1); ixmax = min(ixmax, proj%n_sample(1))
                        jymin = max(jymin, 1); jymax = min(jymax, proj%n_sample(2))
                        kxmax = kxmin + (ixmax-ixmin)
                        kymax = kymin + (jymax-jymin)

                        ! print*,ixmin,ixmax,jymin,jymax,kxmin,kymin,kxmax,kymax
                        do ifilt = 1, proj%nfilter
                            ! if (size(proj%weights(ifilt,ixmin:ixmax,jymin:jymax),1).ne.size(kfilter(kxmin:kxmax,kymin:kymax),1)) then
                            !     print*,shape(proj%weights(ifilt,ixmin:ixmax,jymin:jymax)),shape(kfilter(kxmin:kxmax,kymin:kymax))
                            !     stop
                            ! end if
                            do iweight = 1, proj%nwvars
                                proj%weights(ifilt,iweight,ixmin:ixmax,jymin:jymax) = proj%weights(ifilt,iweight,ixmin:ixmax,jymin:jymax) + &
                                                                                &grid(ilevel)%map(ifilt,iweight,ix,iy)*kfilter(kxmin:kxmax,kymin:kymax)
                                do ivar = 1, proj%nvars
                                    proj%map(ifilt, ivar,iweight,ixmin:ixmax,jymin:jymax) = proj%map(ifilt, ivar,iweight,ixmin:ixmax,jymin:jymax) + &
                                                                        &grid(ilevel)%cube(ifilt,ivar,iweight,ix,iy)*kfilter(kxmin:kxmax,kymin:kymax)
                                end do
                            end do
                        end do
                        ! call circular_projection
                    end do
                end do
                deallocate(kfilter)
            end if
        end do

        ! Finally, normalise using the saved weights
        filtloopmap: do ifilt=1,proj%nfilter
            projvarloopmap: do ivar=1,proj%nvars
                weightvarloop: do iweight=1,proj%nwvars
                    proj%map(ifilt,ivar,iweight,:,:) = proj%map(ifilt,ivar,iweight,:,:)/proj%weights(ifilt,iweight,:,:)
                end do weightvarloop
            end do projvarloopmap
        end do filtloopmap
        
    contains

        subroutine create_kfilter
            implicit none

            allocate(kfilter(1:2*dx_i+1, 1:2*dx_i+1))
            kfilter = 0d0
            do jpix = 1, 2*dx_i+1
                jpix2 = dble(jpix - (dx_i + 1))**2
                do ipix = 1, 2*dx_i+1
                    ipix2 = dble(ipix - (dx_i + 1))**2
                    rpixel2 = ipix2 + jpix2
                    if (rpixel2 .gt. rcut2) cycle
                    ! Compute kernel weight at set distance
                    kfilter(ipix,jpix) = kernel_norm*exp(-0.5d0*rpixel2/kdx2)
                end do
            end do

        end subroutine create_kfilter

    end subroutine grid_projection

    subroutine point_projection(proj,cam,grid)
        implicit none

        type(hydro_projection_handler), intent(inout) :: proj
        type(camera), intent(in)                :: cam
        type(level), dimension(100), intent(in)  :: grid

        integer :: ix, iy, ixh, iyh
        integer :: ilevel, ivar, ifilt, iweight
        integer, dimension(1:2) :: xpix

        ! Initialise map
        allocate(proj%map(1:proj%nfilter,1:proj%nvars,1:proj%nwvars,1:proj%n_sample(1),1:proj%n_sample(2)))
        allocate(proj%weights(1:proj%nfilter,1:proj%nwvars,1:proj%n_sample(1),1:proj%n_sample(2)))
        proj%map     = 0D0
        proj%weights = 0D0

        do ilevel=cam%lmin,cam%lmax
            ! Loop over projected pixels in ilevel
            do ix = grid(ilevel)%imin, grid(ilevel)%imax
                do iy = grid(ilevel)%jmin, grid(ilevel)%jmax
                    ! Get the index of cell in the map grid
                    call get_pos_map(proj%n_sample,grid(ilevel),ix,iy,xpix)
                    ixh = int(xpix(1)); iyh = int(xpix(2))
                    ! Project to map
                    do ifilt = 1, proj%nfilter
                        do iweight = 1, proj%nwvars
                            proj%weights(ifilt,iweight,ixh,iyh) = proj%weights(ifilt,iweight,ixh,iyh) + grid(ilevel)%map(ifilt,iweight,ix,iy)
                            do ivar = 1, proj%nvars
                                proj%map(ifilt,ivar,iweight,ixh,iyh) = proj%map(ifilt,ivar,iweight,ixh,iyh) + grid(ilevel)%cube(ifilt,ivar,iweight,ix,iy)
                            end do
                        end do
                    end do
                end do
            end do
        end do

        ! Finally, normalise using the saved weights
        filtloopmap: do ifilt=1,proj%nfilter
            projvarloopmap: do ivar=1,proj%nvars
                weightvarloop: do iweight=1,proj%nwvars
                    proj%map(ifilt,ivar,iweight,:,:)=proj%map(ifilt,ivar,iweight,:,:)/proj%weights(ifilt,iweight,:,:)
                end do weightvarloop
            end do projvarloopmap
        end do filtloopmap
    end subroutine point_projection

    subroutine upload_projection(proj,cam,bbox,grid)
        implicit none

        type(hydro_projection_handler), intent(inout) :: proj
        type(camera), intent(in)                :: cam
        type(region), intent(in)                :: bbox
        type(level), dimension(100), intent(in) :: grid

        integer :: ix, iy, ixh, iyh
        integer :: ilevel, ivar, ifilt, iweight
        integer, dimension(1:2) :: xpix
        integer :: nx_full,ny_full
        integer :: imin, imax, jmin, jmax
        integer :: i,j,ndom
        real(dbl) :: xmin, ymin

        ! Initialise map
        allocate(proj%map(1:proj%nfilter,1:proj%nvars,1:proj%nwvars,1:proj%n_sample(1),1:proj%n_sample(2)))
        allocate(proj%weights(1:proj%nfilter,1:proj%nwvars,1:proj%n_sample(1),1:proj%n_sample(2)))
        proj%map     = 0D0
        proj%weights = 0D0

        ! Upload to maximum level (lmax)
        nx_full = 2**cam%lmax
        ny_full = 2**cam%lmax
        imin = 1
        imax = int((bbox%xmax-bbox%xmin)*dble(nx_full))
        jmin = 1
        jmax = int((bbox%ymax-bbox%ymin)*dble(ny_full))
        filtlooplmax: do ifilt=1,proj%nfilter
            xloop: do ix = imin,imax
                xmin = ((ix-0.5)/2**cam%lmax)
                yloop: do iy=jmin,jmax
                    ymin=((iy-0.5)/2**cam%lmax)
                    ilevelloop: do ilevel=cam%lmin,cam%lmax-1
                        ndom = 2**ilevel
                        i = int(xmin*ndom)+1
                        j = int(ymin*ndom)+1
                        weightvarlooplmax: do iweight=1,proj%nwvars
                            projvarlooplmax: do ivar=1,proj%nvars
                                grid(cam%lmax)%cube(ifilt,ivar,iweight,ix,iy)=grid(cam%lmax)%cube(ifilt,ivar,iweight,ix,iy) + &
                                                            & grid(ilevel)%cube(ifilt,ivar,iweight,i,j)
                            end do projvarlooplmax
                            grid(cam%lmax)%map(ifilt,iweight,ix,iy)=grid(cam%lmax)%map(ifilt,iweight,ix,iy) + &
                                                            & grid(ilevel)%map(ifilt,iweight,i,j)
                        end do weightvarlooplmax
                    end do ilevelloop
            end do yloop
            end do xloop
        end do filtlooplmax

        do i=1,proj%n_sample(1)
            ix = int(dble(i)/dble(proj%n_sample(1))*dble(imax-imin+1))+imin
            ix = min(ix,imax)
            do j=1,proj%n_sample(2)
                iy = int(dble(j)/dble(proj%n_sample(2))*dble(jmax-jmin+1))+jmin
                iy = min(iy,jmax)
                filtloopmap: do ifilt=1,proj%nfilter
                    weightvarloopmap: do iweight=1,proj%nwvars
                        proj%weights(ifilt,iweight,i,j)=grid(cam%lmax)%map(ifilt,iweight,ix,iy)
                        projvarloopmap: do ivar=1,proj%nvars
                            proj%map(ifilt,ivar,iweight,i,j)=grid(cam%lmax)%cube(ifilt,ivar,iweight,ix,iy)/grid(cam%lmax)%map(ifilt,iweight,ix,iy)
                        end do projvarloopmap
                    end do weightvarloopmap
                end do filtloopmap
            end do
        end do 

    end subroutine upload_projection

    subroutine projection_hydro(repository,type_projection,cam,use_neigh,proj,&
                                &lmax,lmin,nexp_factor,vardict)
        implicit none

        ! Input/output variables
        character(128),intent(in) :: repository
        character(128),intent(in) :: type_projection
        type(camera),intent(inout) :: cam
        logical,intent(in) :: use_neigh
        type(hydro_projection_handler),intent(inout) :: proj
        integer,intent(in),optional :: lmax,lmin
        real(dbl),intent(in),optional :: nexp_factor
        type(dictf90),intent(in),optional :: vardict

        type(region) :: bbox
        integer :: ivx,ivy,ivz
        integer :: ii

        ! Obtain details of the hydro variables stored
        if (.not.present(vardict)) call read_hydrofile_descriptor(repository)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = min(get_required_resolution(cam),amr%nlevelmax)
        if (verbose) write(*,*)'Maximum resolution level: ',amr%nlevelmax
        if (verbose) write(*,*)'Using: ',amr%lmax
        cam%lmin = 1;cam%lmax = amr%lmax
        if(present(lmax)) cam%lmax = max(min(lmax,amr%nlevelmax),1)
        amr%lmax = cam%lmax
        if(present(lmin)) cam%lmin = max(1,min(lmin,amr%lmax))
        if (verbose) write(*,*)'Camera using lmin, lmax: ',cam%lmin,cam%lmax
        call get_bounding_box(cam,bbox)
        bbox%name = 'cube'
        bbox%bulk_velocity = cam%region_velocity
        bbox%criteria_name = 'd_euclid'

        ! Compute the Hilbert curve for the bounding box
        call get_cpu_map(bbox)
        if (verbose) write(*,*)'ncpu: ',amr%ncpu_read
        call get_map_box(cam,bbox)

        if (cam%nsubs>0 .and. verbose)write(*,*)'Excluding substructure: ',cam%nsubs
        call get_map_size(cam,proj%n_sample)

        ! Set up hydro variables quicklook tools
        if (present(vardict)) then
            ! If the user provides their own variable dictionary,
            ! use that one instead of the automatic from the 
            ! hydro descriptor file (RAMSES)
            call get_var_tools(vardict,proj%nvars,proj%varnames,proj%vars)

            ! Do it also for the weight variable
            call get_var_tools(vardict,proj%nwvars,proj%weightvars,proj%wvars)
            
            ! We also do it for the filter variables
            do ii = 1, proj%nfilter
                call get_filter_var_tools(vardict,proj%filters(ii))
            end do

            if (verbose) then
                write(*,*) 'Using variable dicionary from user!'
                write(*,*) 'Number of variables: ',proj%nvars
                write(*,*) 'Number of weight variables: ',proj%nwvars
                write(*,*) 'Number of filters: ',proj%nfilter
                do ii = 1, proj%nvars
                    write(*,*) 'Variable ',ii,' : ',trim(proj%varnames(ii)),' at index ',vardict%get(proj%varnames(ii))
                end do
                do ii = 1, proj%nwvars
                    write(*,*) 'Weight variable ',ii,' : ',trim(proj%weightvars(ii)),' at index ',vardict%get(proj%weightvars(ii))
                end do
            end if

            ! We always need the indexes of the velocities
            ! to perform rotations of gas velocities
            ivx = vardict%get('velocity_x')
            ivy = vardict%get('velocity_y')
            ivz = vardict%get('velocity_z')
        else
            call get_var_tools(varIDs,proj%nvars,proj%varnames,proj%vars)

            ! Do it also for the weight variable
            call get_var_tools(varIDs,proj%nwvars,proj%weightvars,proj%wvars)

            ! We also do it for the filter variables
            do ii = 1, proj%nfilter
                call get_filter_var_tools(varIDs,proj%filters(ii))
            end do

            ! We always need the indexes of the velocities
            ! to perform rotations of gas velocities
            ivx = varIDs%get('velocity_x')
            ivy = varIDs%get('velocity_y')
            ivz = varIDs%get('velocity_z')
        end if

        if (cam%nsubs>0.and.verbose)write(*,*)'Excluding substructure: ',cam%nsubs

        ! Perform projections
        if (use_neigh) then
            if (verbose) write(*,*)'Loading neighbours...'
            call project_cells_neigh
        else
            call project_cells
        end if
        contains

        subroutine project_cells
            implicit none

            logical :: ok_cell,ok_filter,ok_sub,read_gravity
            integer :: i,j,k
            integer :: ipos,icpu,ilevel,ind,idim,iidim,ivar,ifilt,isub,iweight
            integer :: ix,iy,iz,ngrida,nx_full,ny_full,nz_full
            integer :: imin,imax,jmin,jmax
            integer :: nvarh
            integer :: roterr
            character(5) :: nchar,ncharcpu
            character(128) :: nomfich
            real(dbl) :: distance,dx,ksize
            type(vector) :: xtemp,vtemp,gtemp
            integer,dimension(:,:),allocatable :: ngridfile,ngridlevel,ngridbound
            real(dbl),dimension(1:8,1:3) :: xc
            real(dbl),dimension(1:3,1:3) :: trans_matrix
            real(dbl),dimension(1:proj%nvars) :: hvalues
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
            integer,dimension(1:2) :: ii_map
            integer :: ncells
            

            type(level),dimension(1:100) :: grid

            ncells = 0

            ! Check whether we need to read the gravity files
            read_gravity = .false.
            do ivar=1,proj%nvars
                if (proj%varnames(ivar)(1:4) .eq. 'grav') then
                    read_gravity = .true.
                    if (verbose) write(*,*)'Reading gravity files...'
                    exit
                endif
            end do

            ! Compute hierarchy
            do ilevel=1,amr%lmax
                nx_full = 2**ilevel
                ny_full = 2**ilevel
                ! Test this imin whether it should be +1 or not
                imin = int(0D0*dble(nx_full))+1
                imax = int((bbox%xmax-bbox%xmin)*dble(nx_full))+1
                jmin = int(0D0*dble(ny_full))+1
                jmax = int((bbox%ymax-bbox%ymin)*dble(ny_full))+1
                allocate(grid(ilevel)%cube(1:proj%nfilter,1:proj%nvars,1:proj%nwvars,imin:imax,jmin:jmax))
                allocate(grid(ilevel)%map(1:proj%nfilter,1:proj%nwvars,imin:imax,jmin:jmax))
                grid(ilevel)%cube(:,:,:,:,:) = 0D0
                grid(ilevel)%map(:,:,:,:) = 0D0
                grid(ilevel)%imin = imin
                grid(ilevel)%imax = imax
                grid(ilevel)%jmin = jmin
                grid(ilevel)%jmax = jmax
                grid(ilevel)%active = .false.
            end do

            trans_matrix = 0D0
            call new_z_coordinates(bbox%axis,trans_matrix,roterr)
            if (roterr.eq.1) then
                if (verbose) write(*,*) 'Incorrect CS transformation!'
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
                !if (verbose) write(*,*)'Processing file '//TRIM(nomfich)
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
                    if (ngrida>0.and.(ilevel.ge.cam%lmin)&
                        &.and.(ilevel.le.cam%lmax)) then
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
                                    !weight = (min(x(i,3)+dx/2.,bbox%zmax)-max(x(i,3)-dx/2.,bbox%zmin))/dx
                                    !weight = min(1.0d0,max(weight,0.0d0))
                                    if( ix>=grid(ilevel)%imin.and.&
                                        & iy>=grid(ilevel)%jmin.and.&
                                        & ix<=grid(ilevel)%imax.and.&
                                        & iy<=grid(ilevel)%jmax) then
                                        ! Transform position to galaxy frame
                                        xtemp = xorig(i,:)
                                        xtemp = xtemp - bbox%centre
                                        call rotate_vector(xtemp,trans_matrix)
                                        ! Velocity transformed
                                        vtemp = var(i,ind,ivx:ivz)
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
                                        tempvar(0,ivx:ivz) = vtemp
                                        if (read_gravity) tempgrav_var(0,2:4) = gtemp
                                        
                                        filterloop: do ifilt=1,proj%nfilter
                                            if (read_gravity) then
                                                ok_filter = filter_cell(bbox,proj%filters(ifilt),xtemp,dx,tempvar,&
                                                                        &tempson,trans_matrix,tempgrav_var)
                                            else
                                                ok_filter = filter_cell(bbox,proj%filters(ifilt),xtemp,dx,tempvar,&
                                                                        &tempson,trans_matrix)
                                            end if
                                            ! Finally, get hydro data
                                            if (ok_filter) then
                                                if (.not.grid(ilevel)%active) grid(ilevel)%active = .true.
                                                weightvarloop: do iweight=1,proj%nwvars
                                                    if (trim(proj%weightvars(iweight)) == 'counts') then
                                                        weight = 1D0
                                                    else
                                                        !MAX(rho*dx*weight/(bbox%zmax-bbox%zmin),0D0)
                                                        if (read_gravity) then
                                                            weight = proj%wvars(iweight)%myfunction(amr,sim,proj%wvars(iweight),bbox,dx*sim%boxlen,xtemp&
                                                                & ,tempvar,tempson,trans_matrix,tempgrav_var)
                                                        else
                                                            weight = proj%wvars(iweight)%myfunction(amr,sim,proj%wvars(iweight),bbox,dx*sim%boxlen,xtemp&
                                                                & ,tempvar,tempson,trans_matrix,tempgrav_var)
                                                        end if
                                                        ! weight = MAX(weight*dx/(bbox%zmax-bbox%zmin),0D0)
                                                    end if
                                                    grid(ilevel)%map(ifilt,iweight,ix,iy)=grid(ilevel)%map(ifilt,iweight,ix,iy)+weight
                                                    projvarloop: do ivar=1,proj%nvars
                                                        if (read_gravity) then
                                                            map = proj%vars(ivar)%myfunction(amr,sim,proj%vars(ivar),bbox,dx*sim%boxlen,xtemp&
                                                                & ,tempvar,tempson,trans_matrix,tempgrav_var)
                                                        else
                                                            map = proj%vars(ivar)%myfunction(amr,sim,proj%vars(ivar),bbox,dx*sim%boxlen,xtemp&
                                                                & ,tempvar,tempson,trans_matrix)
                                                        end if
                                                        grid(ilevel)%cube(ifilt,ivar,iweight,ix,iy)=grid(ilevel)%cube(ifilt,ivar,iweight,ix,iy)+map*weight
                                                    end do projvarloop
                                                end do weightvarloop
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
                    else if (ngrida>0) then
                        deallocate(xg,son,var,ref,x)
                        if (read_gravity) then
                            deallocate(grav_var)
                        end if
                    endif
                end do levelloop
            end do cpuloop
            if (verbose) write(*,*)'ncells:',ncells

            ! Select the type of projection function needed
            select case (trim(type_projection))
            case ('gauss_deposition')
                call grid_projection(proj,cam,grid,nexp_factor)
            case ('point_deposition')
                call point_projection(proj,cam,grid)
            case ('upload_deposition')
                call upload_projection(proj,cam,bbox,grid)
            case default
                if (verbose) write(*,*)'Deposition type ',trim(type_projection),' not recognised'
                if (verbose) write(*,*)'Falling back to default (upload_projection)'
                call upload_projection(proj,cam,bbox,grid)
            end select
        end subroutine project_cells

        subroutine project_cells_neigh
            implicit none

            logical :: ok_cell,read_gravity,ok_filter,ok_cell_each,ok_sub
            integer :: i,j,k
            integer :: ipos,icpu,ilevel,ind,idim,iidim,ivar,iskip,inbor,ison,isub,iweight
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
            integer ,dimension(1,0:amr%twondim) :: ind_nbor
            logical,dimension(:),allocatable :: ref
            real(dbl) :: rho,map,weight
            real(dbl) :: xmin,ymin
            integer :: ndom
            integer,dimension(1:2) :: n_sample
            integer :: ncells
            integer :: ngrid_current
            integer :: idebug
            integer :: ifilt
            

            type(level),dimension(1:100) :: grid

            ncells = 0
            idebug = 0

            ! Check whether we need to read the gravity files
            read_gravity = .false.
            do ivar=1,proj%nvars
                if (proj%varnames(ivar)(1:4) .eq. 'grav' .or.&
                & trim(proj%varnames(ivar)) .eq. 'neighbour_accuracy') then
                    read_gravity = .true.
                    if (verbose) write(*,*)'Reading gravity files...'
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
                allocate(grid(ilevel)%cube(1:proj%nfilter,1:proj%nvars,1:proj%nwvars,imin:imax,jmin:jmax))
                allocate(grid(ilevel)%map(1:proj%nfilter,1:proj%nwvars,imin:imax,jmin:jmax))
                grid(ilevel)%cube(:,:,:,:,:) = 0D0
                grid(ilevel)%map(:,:,:,:) = 0D0
                grid(ilevel)%imin = imin
                grid(ilevel)%imax = imax
                grid(ilevel)%jmin = jmin
                grid(ilevel)%jmax = jmax
                grid(ilevel)%ngrid = 0
                grid(ilevel)%active = .false.
            end do

            !call los_transformation(cam,trans_matrix)
            trans_matrix = 0D0
            call new_z_coordinates(bbox%axis,trans_matrix,roterr)
            if (roterr.eq.1) then
                if (verbose) write(*,*) 'Incorrect CS transformation!'
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
                !if (verbose) write(*,*)'Processing file '//TRIM(nomfich)
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
                                grid(ilevel)%xg(:,iidim) = xxg(:)
                            end do
                            
                            read(10) ! Skip father index
                            ! Read nbor index
                            do ind=1,amr%twondim
                                read(10)iig
                                nbor(grid(ilevel)%ind_grid(:),ind) = iig(:)
                            end do
                            ! Read son index
                            do ind=1,amr%twotondim
                                iskip = amr%ncoarse+(ind-1)*amr%ngridmax
                                read(10)iig
                                son(grid(ilevel)%ind_grid(:)+iskip) = iig(:)
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
                                    var(grid(ilevel)%ind_grid(:)+iskip,ivar) = xxg(:)
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
                                    grav_var(grid(ilevel)%ind_grid(:)+iskip,1) = xxg(:)
                                    do ivar=1,amr%ndim
                                        read(12)xxg
                                        grav_var(grid(ilevel)%ind_grid(:)+iskip,ivar+1) = xxg(:)
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
                                        vtemp = var(ind_cell(i),ivx:ivz)
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
                                        tempvar(0,:) = var(ind_nbor(1,0),:)
                                        tempson(0)       = son(ind_nbor(1,0))
                                        if (read_gravity) tempgrav_var(0,:) = grav_var(ind_nbor(1,0),:)
                                        tempvar(0,ivx:ivz) = vtemp
                                        if (read_gravity) tempgrav_var(0,2:4) = gtemp
                                        
                                        do inbor=1,amr%twondim
                                            tempvar(inbor,:) = var(ind_nbor(1,inbor),:)
                                            tempson(inbor)       = son(ind_nbor(1,inbor))
                                            if (read_gravity) tempgrav_var(inbor,:) = grav_var(ind_nbor(1,inbor),:)
                                        end do
                                        filterloop: do ifilt=1,proj%nfilter
                                            if (read_gravity) then
                                                ok_filter = filter_cell(bbox,proj%filters(ifilt),xtemp,dx*sim%boxlen,tempvar,&
                                                                        &tempson,trans_matrix,tempgrav_var)
                                            else
                                                ok_filter = filter_cell(bbox,proj%filters(ifilt),xtemp,dx*sim%boxlen,tempvar,&
                                                                        &tempson,trans_matrix)
                                            end if
                                            ! Finally, get hydro data
                                            if (ok_filter) then
                                                if (.not.grid(ilevel)%active) grid(ilevel)%active = .true.
                                                weightvarloop: do iweight=1,proj%nwvars
                                                    if (trim(proj%weightvars(iweight)) == 'counts') then
                                                        weight = 1D0
                                                    else
                                                        !MAX(rho*dx*weight/(bbox%zmax-bbox%zmin),0D0)
                                                        if (read_gravity) then
                                                            weight = proj%wvars(iweight)%myfunction(amr,sim,proj%wvars(iweight),bbox,dx*sim%boxlen,xtemp&
                                                                & ,tempvar,tempson,trans_matrix,tempgrav_var)
                                                        else
                                                            weight = proj%wvars(iweight)%myfunction(amr,sim,proj%wvars(iweight),bbox,dx*sim%boxlen,xtemp&
                                                                & ,tempvar,tempson,trans_matrix,tempgrav_var)
                                                        end if
                                                        ! weight = MAX(weight*dx/(bbox%zmax-bbox%zmin),0D0)
                                                    end if
                                                    grid(ilevel)%map(ifilt,iweight,ix,iy)=grid(ilevel)%map(ifilt,iweight,ix,iy)+weight
                                                    projvarloop: do ivar=1,proj%nvars
                                                        if (read_gravity) then
                                                            map = proj%vars(ivar)%myfunction(amr,sim,proj%vars(ivar),bbox,dx*sim%boxlen,xtemp&
                                                                & ,tempvar,tempson,trans_matrix,tempgrav_var)
                                                        else
                                                            map = proj%vars(ivar)%myfunction(amr,sim,proj%vars(ivar),bbox,dx*sim%boxlen,xtemp&
                                                                & ,tempvar,tempson,trans_matrix)
                                                        end if
                                                        grid(ilevel)%cube(ifilt,ivar,iweight,ix,iy)=grid(ilevel)%cube(ifilt,ivar,iweight,ix,iy)+map*weight
                                                    end do projvarloop
                                                end do weightvarloop
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
            if (verbose) write(*,*)'ncells: ',ncells

            ! Select the type of projection function needed
            select case (trim(type_projection))
            case ('gauss_deposition')
                call grid_projection(proj,cam,grid,nexp_factor)
            case ('point_deposition')
                call point_projection(proj,cam,grid)
            case ('upload_deposition')
                call upload_projection(proj,cam,bbox,grid)
            case default
                if (verbose) write(*,*)'Deposition type ',trim(type_projection),' not recognised'
                if (verbose) write(*,*)'Falling back to default (upload_projection)'
                call upload_projection(proj,cam,bbox,grid)
            end select        
        end subroutine project_cells_neigh

        
    end subroutine projection_hydro

    subroutine projection_parts(repository,cam,proj,part_dict,part_vtypes,tag_file,inverse_tag)
        implicit none
        character(128),intent(in) :: repository
        type(camera),intent(in) :: cam
        type(part_projection_handler),intent(inout) :: proj
        type(dictf90),intent(in),optional :: part_dict,part_vtypes
        character(128),intent(in),optional :: tag_file
        logical,intent(in),optional :: inverse_tag

        type(region) :: bbox
        integer :: ivx,ivy,ivz
        integer :: ii

        ! Obtain details of the particle variables stored
        call read_partfile_descriptor(repository)

        ! Intialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = amr%nlevelmax !min(get_required_resolution(cam),amr%nlevelmax)
        if (verbose) write(*,*)'Maximum resolution level: ',amr%nlevelmax
        if (verbose) write(*,*)'Using: ',amr%lmax
        if (sim%dm .and. sim%hydro) call check_families(repository)
        call get_bounding_box(cam,bbox)
        bbox%name = 'cube'
        bbox%bulk_velocity = cam%region_velocity
        bbox%criteria_name = 'd_euclid'

        ! Compute the Hilbert curve for the bounding box
        call get_cpu_map(bbox)
        call get_map_box(cam,bbox)

        if (cam%nsubs>0 .and. verbose)write(*,*)'Excluding substructure: ',cam%nsubs

        ! Set up the part variables quicklook tools
        if (present(part_dict).and.present(part_vtypes)) then
            ! If the user provides their own particle dictionary,
            ! use that one instead of the automatic from the
            ! particle_file_descriptor.txt (RAMSES)
            call get_partvar_tools(part_dict,part_vtypes,proj%nvars,proj%varnames,proj%vars)

            ! Do it also for the weight variables
            call get_partvar_tools(part_dict,part_vtypes,proj%nwvars,proj%weightvars,proj%wvars)

            ! We also do it for the filter variables
            do ii = 1, proj%nfilter
                call get_filter_part_tools(part_dict,part_vtypes,proj%filters(ii))
            end do

            ! We always need the indexes of the velocities to perform rotation
            ! of particle velocities
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

            ! Since all its alright, we just set the partIDs and partvar_types to the
            ! ones provided by the user
            partIDs = part_dict
            partvar_types = part_vtypes
        else
            ! If the user does not provide their own particle dictionary,
            ! use the automatic one from the particle_file_descriptor.txt (RAMSES)
            call get_partvar_tools(partIDs,partvar_types,proj%nvars,proj%varnames,proj%vars)

            ! Do it also for the weight variables
            call get_partvar_tools(partIDs,partvar_types,proj%nwvars,proj%weightvars,proj%wvars)

            ! We also do it for the filter variables
            do ii = 1, proj%nfilter
                call get_filter_part_tools(partIDs,partvar_types,proj%filters(ii))
            end do

            ! We always need the indexes of the velocities to perform rotation
            ! of particle velocities
            ivx = part_dict%get('velocity_x')
            ivy = part_dict%get('velocity_y')
            ivz = part_dict%get('velocity_z')
        end if

        call project_particles

        contains

        subroutine project_particles
#ifndef LONGINT
            use utils, only:quick_sort_irg,binarysearch_irg
#else
            use utils, only:quick_sort_ilg,binarysearch_ilg
#endif
            use cosmology
            implicit none
    
            logical :: ok_part,ok_tag,ok_filter,ok_sub
            integer :: i,j,k,itag,ifilt,isub,iweight
            integer :: ipos,icpu,ix,iy,ixp1,iyp1,ivar
            integer(irg) :: npart,npart2,nstar,ntag
            integer :: npartsub
            integer :: ncpu2,ndim2
            real(dbl) :: weight,distance,mapvalue,mapvalue_d
#ifdef LONGINT
            integer(ilg) :: mapvalue_i
#else
            integer(irg) :: mapvalue_i
#endif
            integer(1) :: mapvalue_b
            real(dbl) :: dx,dy,ddx,ddy,dex,dey
            real(dbl),dimension(1:3,1:3) :: trans_matrix
            character(5) :: nchar,ncharcpu
            character(6) :: ptype
            character(128) :: nomfich
            type(vector) :: xtemp,vtemp,dcell
            integer,dimension(:),allocatable :: order
            integer,dimension(:),allocatable :: nparttoto
#ifndef LONGINT
            integer(irg),dimension(:),allocatable :: id,tag_id
#else
            integer(ilg),dimension(:),allocatable :: id,tag_id
#endif
            integer,dimension(1:2) :: n_map
            real(dbl),dimension(:,:),allocatable :: part_data_d
            integer(1),dimension(:,:),allocatable :: part_data_b
#ifdef LONGINT
            integer(ilg),dimension(:,:),allocatable :: part_data_i
#else
            integer(irg),dimension(:,:),allocatable :: part_data_i
#endif
            real(dbl),dimension(:,:),allocatable :: x,v
    
            npartsub = 0
#ifndef IMASS
            if (sim%eta_sn .eq. -1D0) then
                if (verbose) write(*,*)': eta_sn=-1 and not IMASS --> should set this up!'
                stop
            end if
#endif
    
            ! If tagged particles file exists, read and allocate array
            if (present(tag_file)) then
                open(unit=58,file=TRIM(tag_file),status='old',form='formatted')
                if (verbose) write(*,*)'Reading particle tags file '//TRIM(tag_file)
                read(58,'(I11)')ntag
                if (verbose) write(*,*)'Number of tagged particles in file: ',ntag
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
                if (verbose) write(*,*)'Sorting list of particle ids for binary search...'
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
            allocate(proj%map(1:proj%nfilter,1:proj%nvars,1:proj%nwvars,0:n_map(1)-1,0:n_map(2)-1))
            allocate(proj%weights(1:proj%nfilter,1:proj%nwvars,0:n_map(1)-1,0:n_map(2)-1))
    
            proj%map = 0D0
    
            ! Cosmological model
            if (sim%aexp.eq.1.and.sim%h0.eq.1)sim%cosmo=.false.
            if (sim%cosmo) then
                call cosmology_model
            else
                sim%time_simu = sim%t
                if (verbose) write(*,*)'Age simu=',sim%time_simu*sim%unit_t/(365.*24.*3600.*1d9)
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
                ! if (verbose) write(*,*)'Processing file '//TRIM(nomfich)
                open(unit=1,file=nomfich,status='old',form='unformatted')
                read(1)ncpu2
                read(1)ndim2
                read(1)npart2
                read(1)
                read(1)nstar
                close(1)
                npart=npart+npart2
            end do
            if (verbose) write(*,*)'Found ',npart,' particles.'
            if(nstar>0)then
                if (verbose) write(*,*)'Found ',nstar,' star particles.'
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
                
                ! 1. Allocate particle data arrays
                allocate(part_data_d(sim%nvar_part_d,1:npart2))
                allocate(part_data_b(sim%nvar_part_b,1:npart2))
                allocate(part_data_i(sim%nvar_part_i,1:npart2))

                allocate(x(1:npart2,1:3))
                allocate(v(1:npart2,1:3))

                if (present(tag_file)) allocate(id(1:npart2))

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

                ! 5. If a tag file is used, get a hold of the particle ids
                if (present(tag_file)) id(:) = part_data_i(partIDs%get('id'),:)
    
                ! 6. Project positions onto the camera frame
                call project_points(cam,npart2,x)
    
                ! 7. Project variables into map for particles
                ! of interest in the region
                partloop: do i=1,npart2
                    distance = 0D0
                    xtemp = x(i,:)
                    xtemp = xtemp - bbox%centre
                    x(i,:) = xtemp
                    call checkifinside(x(i,:),bbox,ok_part,distance)
                    xtemp = xtemp + bbox%centre
                    x(i,:) = xtemp
    
                    ! If we are avoiding substructure, check whether we are safe
                    if (cam%nsubs>0) then
                        ok_sub = .true.
                        do isub=1,cam%nsubs
                            ok_sub = ok_sub .and. filter_sub(cam%subs(isub),x(i,:))
                        end do
                        if (.not.ok_sub)npartsub = npartsub + 1
                        ok_part = ok_part .and. ok_sub
                    end if
                    ! Check if tags are present for particles
                    if (present(tag_file) .and. ok_part) then
                        ok_tag = .false.      
#ifndef LONGINT
                        call binarysearch_irg(ntag,tag_id,id(i),ok_tag)
#else
                        call binarysearch_ilg(ntag,tag_id,id(i),ok_tag)
#endif
                        if (present(inverse_tag)) then
                            if (inverse_tag .and. ok_tag) ok_tag = .false.
                        end if
                        ok_part = ok_tag .and. ok_part
                    endif
                    if (ok_part) then
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
                            ! Rotate the particle velocity
                            vtemp = v(i,:)
                            vtemp = vtemp - bbox%bulk_velocity
                            call rotate_vector(vtemp,trans_matrix)
                            
                            ! Return the x array (now in boxlen units,
                            ! rotated for the camera frame and centered)
                            xtemp = x(i,:)
                            xtemp = xtemp - bbox%centre
                            call rotate_vector(xtemp,trans_matrix)
                            part_data_d(partIDs%get('x'),i) = xtemp%x
                            if (amr%ndim > 1) part_data_d(partIDs%get('y'),i) = xtemp%y
                            if (amr%ndim > 2) part_data_d(partIDs%get('z'),i) = xtemp%z
                            
                            ! Return the velocity array (now rotated to the
                            ! camera frame and corrected for the bulk velocity)
                            part_data_d(partIDs%get('velocity_x'),i) = vtemp%x
                            if (amr%ndim > 1) part_data_d(partIDs%get('velocity_y'),i) = vtemp%y
                            if (amr%ndim > 2) part_data_d(partIDs%get('velocity_z'),i) = vtemp%z

                            filtloopmap: do ifilt=1,proj%nfilter
                                ok_filter = filter_particle(bbox,proj%filters(ifilt),dcell,part_data_d(:,i),part_data_i(:,i),part_data_b(:,i))
                                if (ok_filter) then
                                    projvarloop: do ivar=1,proj%nvars
                                        if (proj%vars(ivar)%vartype==1) then
                                            mapvalue_d = proj%vars(ivar)%myfunction_d(amr,sim,proj%vars(ivar),bbox,dcell &
                                                & ,part_data_d(:,i),part_data_i(:,i),part_data_b(:,i))
                                            mapvalue = mapvalue_d
                                        else if (proj%vars(ivar)%vartype==2) then
                                            mapvalue_i = proj%vars(ivar)%myfunction_i(amr,sim,proj%vars(ivar),bbox,dcell &
                                                & ,part_data_d(:,i),part_data_i(:,i),part_data_b(:,i))
                                            mapvalue = real(mapvalue_i,kind=dbl)
                                        else if (proj%vars(ivar)%vartype==3) then
                                            mapvalue_b = proj%vars(ivar)%myfunction_b(amr,sim,proj%vars(ivar),bbox,dcell &
                                                & ,part_data_d(:,i),part_data_i(:,i),part_data_b(:,i))
                                            mapvalue = real(mapvalue_b,kind=dbl)
                                        end if
                                        weightloop: do iweight=1,proj%nwvars
                                            if (trim(proj%weightvars(iweight)).eq.'cumulative') then
                                                weight = 1d0
                                            else
                                                if (proj%wvars(iweight)%vartype==1) then
                                                    weight = proj%wvars(iweight)%myfunction_d(amr,sim,proj%wvars(iweight),bbox,dcell &
                                                        & ,part_data_d(:,i),part_data_i(:,i),part_data_b(:,i))
                                                else if (proj%wvars(iweight)%vartype==2) then
                                                    weight = proj%wvars(iweight)%myfunction_i(amr,sim,proj%wvars(iweight),bbox,dcell &
                                                        & ,part_data_d(:,i),part_data_i(:,i),part_data_b(:,i))
                                                else if (proj%wvars(iweight)%vartype==3) then
                                                    weight = proj%wvars(iweight)%myfunction_b(amr,sim,proj%wvars(iweight),bbox,dcell &
                                                        & ,part_data_d(:,i),part_data_i(:,i),part_data_b(:,i))
                                                end if
                                            end if
                                            ! Add the weights
                                            proj%weights(ifilt,iweight,ix,iy) = proj%weights(ifilt,iweight,ix,iy) + weight*dex*dey
                                            proj%weights(ifilt,iweight,ix,iyp1) = proj%weights(ifilt,iweight,ix,iyp1) + weight*dex*ddy
                                            proj%weights(ifilt,iweight,ixp1,iy) = proj%weights(ifilt,iweight,ixp1,iy) + weight*ddx*dey
                                            proj%weights(ifilt,iweight,ixp1,iyp1) = proj%weights(ifilt,iweight,ixp1,iyp1) + weight*ddx*ddy

                                            ! Add the map values
                                            proj%map(ifilt,ivar,iweight,ix  ,iy  ) = proj%map(ifilt,ivar,iweight,ix  ,iy  ) + mapvalue*dex*dey*weight
                                            proj%map(ifilt,ivar,iweight,ix  ,iyp1) = proj%map(ifilt,ivar,iweight,ix  ,iyp1) + mapvalue*dex*ddy*weight
                                            proj%map(ifilt,ivar,iweight,ixp1,iy  ) = proj%map(ifilt,ivar,iweight,ixp1,iy  ) + mapvalue*ddx*dey*weight
                                            proj%map(ifilt,ivar,iweight,ixp1,iyp1) = proj%map(ifilt,ivar,iweight,ixp1,iyp1) + mapvalue*ddx*ddy*weight    
                                        end do weightloop
                                    end do projvarloop
                                    nparttoto(ifilt) = nparttoto(ifilt) + 1
                                end if
                            end do filtloopmap
                                
                        endif
                    endif
                end do partloop
                deallocate(x,v,part_data_d,part_data_b,part_data_i)
                if (allocated(id))deallocate(id)
            end do cpuloop
            if (verbose) write(*,*)'> nparttoto: ',nparttoto
            if (verbose) write(*,*)'> npartsub:  ',npartsub
            deallocate(nparttoto)
        end subroutine project_particles
    end subroutine projection_parts

    subroutine healpix_hydro(repository,cam,use_neigh,proj,nside,lmax,lmin,vardict)
        implicit none
        character(128),intent(in) :: repository
        type(camera),intent(inout) :: cam
        logical,intent(in) :: use_neigh
        type(hydro_projection_handler),intent(inout) :: proj
        integer,intent(in) :: nside
        integer,intent(in),optional :: lmax,lmin
        type(dictf90),intent(in),optional :: vardict

        type(region) :: bsphere
        integer :: ivx,ivy,ivz
        integer :: ii

        ! Obtain details of the hydro variables stored
        call read_hydrofile_descriptor(repository)

        ! Initialise parameters of the AMR structure and simulation attributes
        call init_amr_read(repository)
        amr%lmax = amr%nlevelmax
        if (verbose) write(*,*)'Maximum resolution level: ',amr%nlevelmax
        if (verbose) write(*,*)'Using: ',amr%lmax
        cam%lmin = 1;cam%lmax = amr%lmax
        if(present(lmin)) cam%lmin = max(1,min(lmin,amr%lmax))
        if(present(lmax)) cam%lmax = min(amr%lmax,max(lmax,1))
        if (verbose) write(*,*)'Camera using lmin, lmax: ',cam%lmin,cam%lmax
        bsphere%name = 'sphere'
        bsphere%bulk_velocity = cam%region_velocity
        bsphere%criteria_name = 'd_euclid'
        bsphere%centre = cam%centre
        bsphere%rmin = cam%distance
        bsphere%rmax = cam%far_cut_depth
        bsphere%axis = cam%region_axis

        ! Compute the Hilbert curve for the sphere
        call get_cpu_map(bsphere)
        if (verbose) write(*,*)'ncpu: ',amr%ncpu_read

        ! Set up hydro variables quicklook tools
        if (present(vardict)) then
            ! If the user provides their own variable dictionary,
            ! use that one instead of the automatic from the 
            ! hydro descriptor file (RAMSES)
            call get_var_tools(vardict,proj%nvars,proj%varnames,proj%vars)

            ! Do it also for the weight variable
            call get_var_tools(vardict,proj%nwvars,proj%weightvars,proj%wvars)
            
            ! We also do it for the filter variables
            do ii = 1, proj%nfilter
                call get_filter_var_tools(vardict,proj%filters(ii))
            end do

            ! We always need the indexes of the velocities
            ! to perform rotations of gas velocities
            ivx = vardict%get('velocity_x')
            ivy = vardict%get('velocity_y')
            ivz = vardict%get('velocity_z')
        else
            call get_var_tools(varIDs,proj%nvars,proj%varnames,proj%vars)

            ! Do it also for the weight variable
            call get_var_tools(varIDs,proj%nwvars,proj%weightvars,proj%wvars)

            ! We also do it for the filter variables
            do ii = 1, proj%nfilter
                call get_filter_var_tools(varIDs,proj%filters(ii))
            end do

            ! We always need the indexes of the velocities
            ! to perform rotations of gas velocities
            ivx = varIDs%get('velocity_x')
            ivy = varIDs%get('velocity_y')
            ivz = varIDs%get('velocity_z')
        end if

        if (cam%nsubs>0.and.verbose)write(*,*)'Excluding substructure: ',cam%nsubs

        ! Perform projections
        call project_cells_hpix

        contains

        subroutine project_cells_hpix
            use healpix_modules
            implicit none

            logical :: ok_cell,ok_filter,ok_sub
            integer :: i,j,k
            integer :: ipos,icpu,ilevel,ind,idim,iidim
            integer :: ivar,iskip,iweight
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
            real(dbl),dimension(:,:,:),allocatable :: proj_rho
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
            if (verbose) write(*,*)'Total number of healpix map: ',ns+1
            allocate(proj%map(1:proj%nfilter,1:proj%nvars,1:proj%nwvars,1:1,0:ns))
            allocate(proj_rho(1:proj%nfilter,1:proj%nwvars,0:ns))
            allocate(listpix(0:ns))
            proj%map = 0D0
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
                ! if (verbose) write(*,*)'Processing file '//TRIM(nomfich)
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
                                        vtemp = var(i,ind,ivx:ivz)
                                        vtemp = vtemp - bsphere%bulk_velocity
                                        call rotate_vector(vtemp,trans_matrix)
                                        var(i,ind,ivx:ivz) = vtemp

                                        ! Get neighbours
                                        allocate(tempvar(0:amr%twondim,nvarh))
                                        allocate(tempson(0:amr%twondim))
                                        ! Just add central cell as we do not want neighbours
                                        tempvar(0,:) = var(i,ind,:)
                                        tempson(0)       = son(i,ind)
                                        tempvar(0,ivx:ivz) = vtemp

                                        ! Do loop over filters
                                        filterloop: do ifilt=1,proj%nfilter
                                            ok_filter = filter_cell(bsphere,proj%filters(ifilt),xtemp,dx*sim%boxlen,tempvar,&
                                                                        &tempson,trans_matrix)
                                            if (ok_filter) then

                                                ! Get variable values
                                                do j=1,size(listpix_clean)
                                                    ix = listpix_clean(j)
                                                    projvarloop: do ivar=1,proj%nvars
                                                        map = proj%vars(ivar)%myfunction(amr,sim,proj%vars(ivar),bsphere,dx*sim%boxlen,xtemp&
                                                                & ,tempvar,tempson,trans_matrix)
                                                            do iweight = 1, proj%nwvars
                                                                ! Get weight
                                                                if (trim(proj%weightvars(iweight)).eq.'cumulative') then
                                                                    rho = 1d0
                                                                else
                                                                    rho = proj%wvars(iweight)%myfunction(amr,sim,proj%wvars(iweight),bsphere,dx*sim%boxlen,xtemp&
                                                                            & ,tempvar,tempson,trans_matrix)
                                                                end if
                                                                proj_rho(ifilt,iweight,ix) = proj_rho(ifilt,iweight,ix)+rho*dx*sim%boxlen*weight/(bsphere%rmax-bsphere%rmin)
                                                                proj%map(ifilt,ivar,iweight,1,ix) = proj%map(ifilt,ivar,iweight,1,ix)+map*rho*dx*weight/(bsphere%rmax-bsphere%rmin)
                                                            end do
                                                    end do projvarloop
                                                end do
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
            if (verbose) write(*,*)'ncells:',ncells

            ! Renormalise cells to compute weighted values
            filtloopmap: do ifilt=1,proj%nfilter
                do i=0,ns
                    projvarloopmap: do ivar=1,proj%nvars
                        weightvarloopmap: do iweight=1,proj%nwvars
                            proj%map(ifilt,ivar,iweight,1,i) = proj%map(ifilt,ivar,iweight,1,i)/proj_rho(ifilt,iweight,i)
                        end do weightvarloopmap
                    end do projvarloopmap
                end do
            end do filtloopmap

        end subroutine project_cells_hpix

    end subroutine healpix_hydro
end module maps