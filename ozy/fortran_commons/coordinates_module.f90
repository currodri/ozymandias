!--------------------------------------------------------------------------
! ozymandias:coordinates_module.f90
!--------------------------------------------------------------------------
!
! MODULE: coordinate_systems
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> types, functions, operators and routines useful for coordinate
!> transformations.
!
!> @details  
!> define different coordinate transformations
! 
!
!> @date 28/3/2021   0.1 fortran routines for ozymandias
!--------------------------------------------------------------------------
module coordinate_systems
    use local
    use vectors
    use rotations
    use basis_representations

    implicit none
    public
    contains

    real(dbl) function r_sphere(p)
        type(vector), intent(in) :: p
        r_sphere = sqrt(p%x**2+p%y**2+p%z**2)
    end function r_sphere

    real(dbl) function theta_sphere(p)
        type(vector), intent(in) :: p
        real(dbl) :: r
        r = sqrt(p%x**2+p%y**2+p%z**2)
        theta_sphere = dacos(p%z/r)
    end function theta_sphere

    real(dbl) function phi_sphere(p)
        type(vector), intent(in) :: p
        phi_sphere = datan2(p%y,p%x)
    end function phi_sphere

    real(dbl) function r_cyl(p)
        type(vector), intent(in) :: p
        r_cyl = sqrt(p%x**2 + p%y**2)
    end function r_cyl

    real(dbl) function phi_cyl(p)
        type(vector), intent(in) :: p
        real(dbl) :: r
        r = sqrt(p%x**2 + p%y**2)
        if (p%x.eq.0d0.and.p%y.eq.0d0) then
            phi_cyl = 0d0
        elseif (p%x >= 0d0) then
            phi_cyl = asin(p%y/r)
        elseif (p%x > 0d0) then
            phi_cyl = atan2(p%y,p%x)
        elseif (p%x < 0d0) then
            phi_cyl = -asin(p%y/r) + 3.14159265359
        endif
    end function phi_cyl

    !---------------------------------------------------------------
    ! Subroutine: SPHERICAL BASIS FROM CARTESIAN
    !
    ! Obtain the local spherical coordinates basis given the
    ! position vector.
    !---------------------------------------------------------------
    subroutine spherical_basis_from_cartesian(p,spher_basis)
        type(vector), intent(in) :: p
        type(basis), intent(inout) :: spher_basis
        real(dbl), dimension(3,3) :: R
        real(dbl) :: theta,phi
        real(sgl),parameter :: thr = 1.0E-8
        integer :: i,j
        call initialise_basis(spher_basis)
        theta = theta_sphere(p)
        phi = phi_sphere(p)
        R(1,1) = dsin(theta)*dcos(phi)
        R(1,2) = dsin(theta)*dsin(phi)
        R(1,3) = dcos(theta)
        R(2,1) = dcos(theta)*dcos(phi)
        R(2,2) = dcos(theta)*dsin(phi)
        R(2,3) = -dsin(theta)
        R(3,1) = -dsin(phi)
        R(3,2) = dcos(phi)
        R(3,3) = 0d0
        do i=1,3
            do j=1,3
                if (abs(R(i,j)).lt.thr) R(i,j) = 0.0
            end do
        end do
        spher_basis%u(1) = transpose(R) * spher_basis%u(1)
        spher_basis%u(2) = transpose(R) * spher_basis%u(2)
        spher_basis%u(3) = transpose(R) * spher_basis%u(3)
    end subroutine spherical_basis_from_cartesian

    !---------------------------------------------------------------
    ! Subroutine: CYLINDRICAL BASIS FROM CARTESIAN
    !
    ! Obtain the local cylindrical coordinates basis given the
    ! position vector.
    !---------------------------------------------------------------
    subroutine cylindrical_basis_from_cartesian(p,cyl_basis)
        implicit none
        type(vector), intent(in) :: p
        type(basis), intent(inout) :: cyl_basis
        real(dbl) :: phi
        phi = phi_cyl(p)
        call initialise_basis(cyl_basis)
        cyl_basis%u(1) = dcos(phi) * cyl_basis%u(1) + dsin(phi) * cyl_basis%u(2)
        cyl_basis%u(2) = cyl_basis%u(3) * cyl_basis%u(1)
    end subroutine cylindrical_basis_from_cartesian

    !---------------------------------------------------------------
    ! Subroutine: ROTATION MATRIX FOR NEW Z AXIS
    !
    ! This routine obtains a rotation matrix which allows a unitary
    ! transformation from the originally aligned (x,y,z) basis to
    ! a basis with its z axis aligned with the provided axis vector.
    ! 
    ! Useful for computing vectors with respect to the angular
    ! momentum vector of a group object (halo, galaxy...).
    !---------------------------------------------------------------
    subroutine new_z_coordinates(axis,transformation_matrix,errormsg)
        implicit none
        type(vector), intent(in) :: axis
        type(vector) :: new_axis,temp_axis,ZZ
        real(dbl), dimension(3,3), intent(inout) :: transformation_matrix
        integer, intent(inout) :: errormsg
        real(dbl), dimension(3,3) :: R
        integer :: i
        real(dbl) :: theta
        R = 0d0
        ZZ%z = 1d0

        ! Is the axis vector pointing in the z direction?
        if ((axis.DOT.ZZ).eq.1d0) then
            errormsg = 0
            if (axis%z < 0d0) then
                do i=1,2
                    transformation_matrix(i,i) = 1d0
                end do
                transformation_matrix(3,3) = -1d0
            else
                do i=1,3
                    transformation_matrix(i,i) = 1d0
                end do
            endif
        else
            ! axis is not aligned with z axis, so we are going to apply rotations
            new_axis = axis
            temp_axis = axis
            new_axis%z = 0d0
            ! Find angle between axis and the x axis
            theta = dacos(new_axis%x/magnitude(new_axis))
            if (axis%y < 0d0) theta = -theta
            ! Now rotate by this much around z axis
            call euler_matrix(R,3,theta,'a')
            call rotate_vector(temp_axis,R)
            transformation_matrix = R
            ! Now find the angle between axis and the z axis
            theta = dacos(temp_axis%z/magnitude(temp_axis))
            ! And rotate around the y axis
            R = 0d0
            call euler_matrix(R,2,theta,'a')
            transformation_matrix = matmul(R,transformation_matrix)
            call rotate_vector(temp_axis,R)
            if (abs((temp_axis.DOT.ZZ)-1D0)<1D-8) then
                errormsg = 0
            else
                write(*,*)(temp_axis.DOT.ZZ)
                errormsg = 1
            endif
        endif
    end subroutine new_z_coordinates
end module coordinate_systems
!--------------------------------------------------------------------------
!
! MODULE: geometrical_regions
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> types and subroutines useful for the selection of geometrical
!> regions
!
!> @details  
!> define derived type region and the ways it can be handled
! 
!
!> @date 28/3/2021   0.1 fortran routines for ozymandias
!--------------------------------------------------------------------------
module geometrical_regions
    use local
    use vectors
    use basis_representations
    
    type region
        character(128) :: name
        type(vector) ::      centre,axis,bulk_velocity
        real(dbl)                ::       xmin,xmax,ymin,ymax,zmin,zmax
        real(dbl)                ::       rmin,rmax
        real(dbl)                ::       angle
        character(128)          ::       criteria_name
    end type region

    contains

    !---------------------------------------------------------------
    ! Subroutine: BOUNDING BOX LIMITS
    !
    ! Since the cpu mapping of RAMSES requires cartesian limits for
    ! a box inside the simulation domain, this routine computes the
    ! minimum bounding box for the arbitrary geometrical region of
    ! interest.
    !---------------------------------------------------------------
    subroutine limits(reg,lim)
        use coordinate_systems
        implicit none
        type(region), intent(in) :: reg
        real(dbl), dimension(1:3,1:2), intent(inout) :: lim
        real(dbl),dimension(1:3)     ::       u_axis,v,w,array
        real(dbl),dimension(1:8,1:3) ::       box_cyl
        real(dbl) :: cyl_rmax
        real(dbl)                    ::       cone_rmax
        type(basis) :: init_basis, out_basis
        type(vector) :: temp_vec
        integer                         ::       i,j
        real(dbl),dimension(3,3) :: trans_matrix
        integer :: roterr
        box_cyl = 0D0
        lim = 0D0
        select case (TRIM(reg%name))
        case ('cube')
            lim(1,1) = reg%xmin; lim(1,2) = reg%xmax
            lim(2,1) = reg%ymin; lim(2,2) = reg%ymax
            lim(3,1) = reg%zmin; lim(3,2) = reg%zmax
            ! if (reg%axis%z .ne. 1D0) then
            !     trans_matrix = 0D0
            !     call new_z_coordinates(reg%axis,trans_matrix,roterr)
            !     if (roterr.eq.1) then
            !         write(*,*) 'Incorrect CS transformation!'
            !         stop
            !     endif
            !     temp_vec = (/reg%xmin,reg%ymin,reg%zmin/)
            !     call rotate_vector(temp_vec,transpose(trans_matrix))
            !     lim(1,1) = min(temp_vec%x,lim(1,1)); lim(1,2) = max(temp_vec%x,lim(1,2))
            !     lim(2,1) = min(temp_vec%y,lim(2,1)); lim(3,2) = max(temp_vec%y,lim(2,2))
            !     lim(3,1) = min(temp_vec%z,lim(3,1)); lim(3,2) = max(temp_vec%z,lim(3,2))

            !     temp_vec = (/reg%xmin,reg%ymin,reg%zmax/)
            !     call rotate_vector(temp_vec,transpose(trans_matrix))
            !     lim(1,1) = min(temp_vec%x,lim(1,1)); lim(1,2) = max(temp_vec%x,lim(1,2))
            !     lim(2,1) = min(temp_vec%y,lim(2,1)); lim(3,2) = max(temp_vec%y,lim(2,2))
            !     lim(3,1) = min(temp_vec%z,lim(3,1)); lim(3,2) = max(temp_vec%z,lim(3,2))

            !     temp_vec = (/reg%xmin,reg%ymax,reg%zmin/)
            !     call rotate_vector(temp_vec,transpose(trans_matrix))
            !     lim(1,1) = min(temp_vec%x,lim(1,1)); lim(1,2) = max(temp_vec%x,lim(1,2))
            !     lim(2,1) = min(temp_vec%y,lim(2,1)); lim(3,2) = max(temp_vec%y,lim(2,2))
            !     lim(3,1) = min(temp_vec%z,lim(3,1)); lim(3,2) = max(temp_vec%z,lim(3,2))

            !     temp_vec = (/reg%xmin,reg%ymax,reg%zmax/)
            !     call rotate_vector(temp_vec,transpose(trans_matrix))
            !     lim(1,1) = min(temp_vec%x,lim(1,1)); lim(1,2) = max(temp_vec%x,lim(1,2))
            !     lim(2,1) = min(temp_vec%y,lim(2,1)); lim(3,2) = max(temp_vec%y,lim(2,2))
            !     lim(3,1) = min(temp_vec%z,lim(3,1)); lim(3,2) = max(temp_vec%z,lim(3,2))

            !     temp_vec = (/reg%xmax,reg%ymin,reg%zmin/)
            !     call rotate_vector(temp_vec,transpose(trans_matrix))
            !     lim(1,1) = min(temp_vec%x,lim(1,1)); lim(1,2) = max(temp_vec%x,lim(1,2))
            !     lim(2,1) = min(temp_vec%y,lim(2,1)); lim(3,2) = max(temp_vec%y,lim(2,2))
            !     lim(3,1) = min(temp_vec%z,lim(3,1)); lim(3,2) = max(temp_vec%z,lim(3,2))

            !     temp_vec = (/reg%xmin,reg%ymin,reg%zmax/)
            !     call rotate_vector(temp_vec,transpose(trans_matrix))
            !     lim(1,1) = min(temp_vec%x,lim(1,1)); lim(1,2) = max(temp_vec%x,lim(1,2))
            !     lim(2,1) = min(temp_vec%y,lim(2,1)); lim(3,2) = max(temp_vec%y,lim(2,2))
            !     lim(3,1) = min(temp_vec%z,lim(3,1)); lim(3,2) = max(temp_vec%z,lim(3,2))

            !     temp_vec = (/reg%xmax,reg%ymin,reg%zmax/)
            !     call rotate_vector(temp_vec,transpose(trans_matrix))
            !     lim(1,1) = min(temp_vec%x,lim(1,1)); lim(1,2) = max(temp_vec%x,lim(1,2))
            !     lim(2,1) = min(temp_vec%y,lim(2,1)); lim(3,2) = max(temp_vec%y,lim(2,2))
            !     lim(3,1) = min(temp_vec%z,lim(3,1)); lim(3,2) = max(temp_vec%z,lim(3,2))

            !     temp_vec = (/reg%xmax,reg%ymax,reg%zmax/)
            !     call rotate_vector(temp_vec,transpose(trans_matrix))
            !     lim(1,1) = min(temp_vec%x,lim(1,1)); lim(1,2) = max(temp_vec%x,lim(1,2))
            !     lim(2,1) = min(temp_vec%y,lim(2,1)); lim(3,2) = max(temp_vec%y,lim(2,2))
            !     lim(3,1) = min(temp_vec%z,lim(3,1)); lim(3,2) = max(temp_vec%z,lim(3,2))

            ! else
            !     lim(1,1) = reg%xmin; lim(1,2) = reg%xmax
            !     lim(2,1) = reg%ymin; lim(2,2) = reg%ymax
            !     lim(3,1) = reg%zmin; lim(3,2) = reg%zmax
            ! endif
        case ('sphere')
            lim(1,1) = reg%centre%x - reg%rmax; lim(1,2) = reg%centre%x + reg%rmax
            lim(2,1) = reg%centre%y - reg%rmax; lim(2,2) = reg%centre%y + reg%rmax
            lim(3,1) = reg%centre%z - reg%rmax; lim(3,2) = reg%centre%z + reg%rmax
        case ('cylinder')
            cyl_rmax = sqrt(reg%rmax**2+((reg%zmax-reg%zmin)/2D0)**2)
            lim(1,1) = reg%centre%x - cyl_rmax; lim(1,2) = reg%centre%x + cyl_rmax
            lim(2,1) = reg%centre%y - cyl_rmax; lim(2,2) = reg%centre%y + cyl_rmax
            lim(3,1) = reg%centre%z - cyl_rmax; lim(3,2) = reg%centre%z + cyl_rmax
            ! OLD METHOD - left for benchmarking
            ! call initialise_basis(init_basis)
            ! init_basis%u(1) = reg%axis
            ! array = reg%centre
            ! call mGramSchmidt(init_basis,out_basis)
            ! u_axis = out_basis%u(1)
            ! v = out_basis%u(2)
            ! w = out_basis%u(3)
            ! write(*,*)v,w
            ! do i=1,3
            !     box_cyl(1,i) = array(i) + reg%zmax*u_axis(i) + 1.4142136*reg%rmax*v(i)
            !     box_cyl(2,i) = array(i) + reg%zmax*u_axis(i) - 1.4142136*reg%rmax*v(i)
            !     box_cyl(3,i) = array(i) + reg%zmax*u_axis(i) + 1.4142136*reg%rmax*w(i)
            !     box_cyl(4,i) = array(i) + reg%zmax*u_axis(i) - 1.4142136*reg%rmax*w(i)
            !     box_cyl(5,i) = array(i) + reg%zmin*u_axis(i) + 1.4142136*reg%rmax*v(i)
            !     box_cyl(6,i) = array(i) + reg%zmin*u_axis(i) - 1.4142136*reg%rmax*v(i)
            !     box_cyl(7,i) = array(i) + reg%zmin*u_axis(i) + 1.4142136*reg%rmax*w(i)
            !     box_cyl(8,i) = array(i) + reg%zmin*u_axis(i) - 1.4142136*reg%rmax*w(i)
            ! end do
            ! do i=1,3
            !     lim(i,1) = MINVAL(box_cyl(:,i)); lim(i,2) = MAXVAL(box_cyl(:,i))
            ! end do
        case ('cone')
            call initialise_basis(init_basis)
            init_basis%u(1) = reg%axis
            array = reg%centre
            call mGramSchmidt(init_basis,out_basis)
            u_axis = out_basis%u(1)
            v = out_basis%u(2)
            w = out_basis%u(3)
            cone_rmax = reg%rmax * tan(reg%angle)
            do i=1,3
                box_cyl(1,i) = array(i) + reg%rmax*u_axis(i) + 1.4142136*cone_rmax*v(i)
                box_cyl(2,i) = array(i) + reg%rmax*u_axis(i) - 1.4142136*cone_rmax*v(i)
                box_cyl(3,i) = array(i) + reg%rmax*u_axis(i) + 1.4142136*cone_rmax*w(i)
                box_cyl(4,i) = array(i) + reg%rmax*u_axis(i) - 1.4142136*cone_rmax*w(i)
                box_cyl(5,i) = array(i) + reg%rmin*u_axis(i) + 1.4142136*cone_rmax*v(i)
                box_cyl(6,i) = array(i) + reg%rmin*u_axis(i) - 1.4142136*cone_rmax*v(i)
                box_cyl(7,i) = array(i) + reg%rmin*u_axis(i) + 1.4142136*cone_rmax*w(i)
                box_cyl(8,i) = array(i) + reg%rmin*u_axis(i) - 1.4142136*cone_rmax*w(i)
            end do
            do i=1,3
                lim(i,1) = MINVAL(box_cyl(:,i)); lim(i,2) = MAXVAL(box_cyl(:,i))
            end do
        case default
            do i=1,3
                lim(i,1) = 0; lim(i,2) = 1
            end do
        end select
    end subroutine limits

    !---------------------------------------------------------------
    ! Subroutine: CHECK IF REGION CONTAINS POINT
    !
    ! Given a particular geometrical region, it determines if it
    ! contains a given point. It also returns the distance, as
    ! defined for each region, from the origin.
    ! 
    ! WARNING: Point needs to be refered to the center of the
    ! region.
    !---------------------------------------------------------------
    subroutine checkifinside(pos,reg,ok,distance)
        implicit none
        real(dbl),dimension(1:3),intent(in) ::       pos
        type(region),intent(in)                ::       reg
        logical,intent(inout)                     ::       ok
        real(dbl),intent(inout)                ::       distance
        type(vector) :: p
        p = pos
        if (reg%name.eq.'cube') then
            call cube(p,reg,ok,distance)
        elseif (reg%name.eq.'sphere') then
            call sphere(p,reg,ok,distance)
        elseif (reg%name.eq.'cylinder') then
            call cylinder(p,reg,ok,distance)
        elseif (reg%name.eq.'cone') then
            call cone(p,reg,ok,distance)
        else
            write(*,*)'The selected region type is still not supported.'
            write(*,*)'Aborting!'
            stop
        endif
    end subroutine checkifinside

    subroutine cube(p,reg,ok,distance)
        implicit none
        type(vector),intent(in) :: p
        type(region),intent(in)                ::       reg
        logical,intent(inout)                     ::       ok
        real(dbl),intent(inout)                ::       distance
        type(vector) :: ptemp
        ptemp = p + reg%centre
        ok = (reg%xmin <= ptemp%x.and.ptemp%x <= reg%xmax.and.&
                &reg%ymin <= ptemp%y.and.ptemp%y <= reg%ymax.and.&
                &reg%zmin <= ptemp%z.and.ptemp%z <= reg%zmax)
        if (ok) then
            distance=magnitude(ptemp)
        endif
    end subroutine cube

    subroutine sphere(p,reg,ok,distance)
        implicit none
        type(vector),intent(in) :: p
        type(region),intent(in)                ::       reg
        logical,intent(inout)                     ::       ok
        real(dbl),intent(inout)                ::       distance

        distance = magnitude(p)
        ok = (reg%rmin <= distance.and.distance <= reg%rmax)
    end subroutine sphere

    subroutine cylinder(p,reg,ok,distance)
        implicit none
        type(vector),intent(in) :: p
        type(region),intent(in)                ::       reg
        logical,intent(inout)                     ::       ok
        real(dbl),intent(inout)                ::       distance

        distance = sqrt(p%x**2 + p%y**2)
        ok = (reg%rmin <= distance.and.distance <= reg%rmax.and.&
                &reg%zmin <= p%z.and.p%z <= reg%zmax)
    end subroutine cylinder

    subroutine cone(p,reg,ok,distance)
        implicit none
        type(vector),intent(in) :: p
        type(region),intent(in)                ::       reg
        logical,intent(inout)                     ::       ok
        real(dbl),intent(inout)                ::       distance
        real(dbl) :: theta

        distance = p.DOT.reg%axis
        theta = acos(distance/magnitude(p))
        ok = (reg%rmin <= distance.and.distance <= reg%rmax.and.&
                &reg%angle >= abs(theta))
    end subroutine cone

end module geometrical_regions