module math_utils
    public
    contains
        recursive function cross_product(a, b) result(cross)
            implicit none
            real(KIND=8), dimension(1:3) :: cross
            real(KIND=8), dimension(1:3), intent(in) :: a, b

            cross(1) = a(2) * b(3) - a(3) * b(2)
            cross(2) = a(3) * b(1) - a(1) * b(3)
            cross(3) = a(1) * b(2) - a(2) * b(1)
        end function cross_product
        !----------------------------------------------------------------
        ! ALL THIS FUNCTIONS ASSUME ORIGIN OF AXIS AT (0,0,0)
        recursive function r_sphere(p)
            implicit none
            real(KIND=8),dimension(1:3), intent(in) :: p
            real(KIND=8)                            :: r_sphere
            r_sphere = sqrt(sum((p**2)))
        end function r_sphere
        recursive function theta_sphere(p)
            implicit none
            real(KIND=8),dimension(1:3), intent(in) :: p
            real(KIND=8)                            :: r,theta_sphere
            r = sqrt(sum(p**2))
            theta_sphere = acos(p(3)/r)
        end function theta_sphere
        recursive function phi_sphere(p)
            implicit none
            real(KIND=8),dimension(1:3), intent(in) :: p
            real(KIND=8)                            :: phi_sphere
            phi_sphere = atan2(p(2),p(1))
        end function phi_sphere
        recursive function r_cyl(p)
            implicit none
            real(KIND=8),dimension(1:3), intent(in) :: p
            real(KIND=8)                            :: r_cyl
            r_cyl = sqrt(p(1)**2 + p(2)**2)
        end function r_cyl
        recursive function phi_cyl(p)
            implicit none
            real(KIND=8),dimension(1:3), intent(in) :: p
            real(KIND=8)                            :: r_cyl,phi_cyl
            r_cyl = sqrt(p(1)**2 + p(2)**2)
            if (p(1).eq.0d0.and.p(2).eq.0d0) then
                phi_cyl = 0d0
            elseif (p(1) >= 0d0) then
                phi_cyl = asin(p(2)/r_cyl)
            elseif (p(1) > 0d0) then
                phi_cyl = atan2(p(2),p(1))
            elseif (p(1) < 0d0) then
                phi_cyl = -asin(p(2)/r_cyl) + 3.14159265359
            endif
        end function phi_cyl
        !----------------------------------------------------------------
        recursive function matrix_vector_multiplication(A,b) result(res)
            implicit none
            real(KIND=8),dimension(1:3,1:3),intent(in) :: A
            real(KIND=8),dimension(1:3),intent(in) :: b
            real(KIND=8),dimension(1:3) :: res
            integer :: i,j
            do i=1,3
                do j=1,3
                    res(i) = 0d0
                end do
            end do
            do i=1,3
                do j=1,3
                    res(i) = res(i) + A(i,j)*b(j)
                end do
            end do
        end function matrix_vector_multiplication
        recursive function rotated_vector_3D(a,dim,angle) result(res)
            implicit none
            real(KIND=8),dimension(1:3),intent(in) :: a
            integer,intent(in)                     :: dim
            real(KIND=8),intent(in)                :: angle
            real(KIND=8),dimension(1:3)            :: res
            real(KIND=8),dimension(1:3,1:3)        :: R

            select case (dim)
            case (1)
                R(1,1) = 1; R(1,2) = 0; R(1,3) = 0
                R(2,1) = 0; R(2,2) = dcos(angle); R(2,3) = dsin(angle)
                R(3,1) = 0; R(3,2) = -dsin(angle); R(3,3) = dcos(angle)
            case(2)
                R(1,1) = dcos(angle); R(1,2) = 0; R(1,3) = -dsin(angle)
                R(2,1) = 0; R(2,2) = 1; R(2,3) = 0
                R(3,1) = dsin(angle); R(3,2) = 0; R(3,3) = dcos(angle)
            case(3)
                R(1,1) = dcos(angle); R(1,2) = dsin(angle); R(1,3) = 0
                R(2,1) = -dsin(angle); R(2,2) = dcos(angle); R(2,3) = 0
                R(3,1) = 0; R(3,2) = 0; R(3,3) = 1
            case default
                write(*,*)'dim must be 1,2 or 3.'
                write(*,*)'Aborting!'
                stop
            end select
            ! rotated_vector_3D = RESHAPE(MATMUL(R,RESHAPE(a, (/3,1/))), (/3/))
            res = matrix_vector_multiplication(R,a)
        end function rotated_vector_3D
        subroutine modify_reference_frame(na,a,centre,axis)
            implicit none
            integer                          :: na,i
            real(KIND=8),dimension(1:na,1:3) :: a
            real(KIND=8),dimension(1:3)      :: centre,axis,new_axis
            real(KIND=8)                     :: theta
            write(*,*)axis(1),axis(2),axis(3),sqrt(sum(axis**2))
            ! First translate to the new origin
            do i=1,na
                a(i,:) = a(i,:) - centre
            end do
            ! Is the axis vector pointing in the z direction?
            if (axis(1).eq.0d0.and.axis(2).eq.0d0.and.axis(3).ne.0d0) then
                if (axis(3) < 0d0) then
                    do i=1,na
                        a(i,3)=-a(i,3)
                    end do
                endif
            else
                ! axis is not aligned with z axis, so we are going to apply rotations
                new_axis = axis
                ! Find angle between axis and the x axis
                theta = dacos(axis(1))
                if (axis(2) < 0d0) theta = -theta
                ! Now rotate a by this much around z axis
                do i=1,na
                    a(i,:) = rotated_vector_3D(a(i,:),3,theta)
                end do
                new_axis = rotated_vector_3D(new_axis,3,theta)
                write(*,*)'new_axis: ',new_axis(1),new_axis(2),new_axis(3),sqrt(sum(new_axis**2))
                ! Now find the angle between new_axis and the z axis
                theta = dacos(axis(3))
                ! And rotate around the y axis
                do i=1,na
                    a(i,:) = rotated_vector_3D(a(i,:),2,theta)
                end do
                new_axis = rotated_vector_3D(new_axis,2,theta)
                write(*,*)'new_axis: ',new_axis(1),new_axis(2),new_axis(3),sqrt(sum(new_axis**2))
            endif
        end subroutine modify_reference_frame
        subroutine ortho_find(a,b,c)
            implicit none
            real(KIND=8),dimension(1:3) :: a,b,c
            real(KIND=8)                :: x1,x2,y1,y2,z1,z2,vec_norm
            integer                     :: i
            vec_norm = sqrt(sum(a**2))
            do i=1,3
                a(i) = a(i)/vec_norm
            end do
            x1 = a(1)
            y1 = a(2)
            z1 = a(3)
            if (z1.ne.0d0) then
                x2 = 1d0
                y2 = 0d0
                z2 = -(x1/z1)
            elseif (y1.ne.0d0) then
                x2 = 0d0
                z2 = 1d0
                y2 = -(z1/y1)
            else
                y2 = 1d0
                z2 = 0d0
                x2 = -(y1/x1)
            endif
            b = (/x2/vec_norm,y2/vec_norm,z2/vec_norm/)
            c = cross_product(a,b)
        end subroutine
        subroutine sphere_basis(pos,e_r,e_theta,e_phi)
            implicit none
            real(KIND=8),dimension(1:3) :: pos,e_r,e_theta,e_phi
            real(KIND=8) :: theta

            e_r = pos/sqrt(sum(pos**2))
            theta = theta_sphere(pos)
            e_phi = cross_product((/0d0,0d0,1d0/),e_r)/sin(theta)
            e_theta = cross_product(e_phi,e_r)
        end subroutine
        subroutine cyl_basis(pos,e_r,e_z,e_phi)
            implicit none
            real(KIND=8),dimension(1:3) :: pos,e_r,e_z,e_phi

            e_r = (/pos(1),pos(2),0d0/)/sqrt(pos(1)**2+pos(2)**2)
            e_z = (/0d0,0d0,1d0/)
            e_phi = cross_product(e_z,e_r)
        end subroutine
end module math_utils