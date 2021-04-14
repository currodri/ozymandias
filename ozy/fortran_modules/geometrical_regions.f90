module geometrical_regions
    use math_utils
    public
    !--------------------------------------------------------------------------
    ! Useful types
    !--------------------------------------------------------------------------

    type region
        character(128)          ::       name
        real(KIND=8),dimension(1:3) ::       centre,axis,bulk_velocity
        real(KIND=8)                ::       xmin,xmax,ymin,ymax,zmin,zmax
        real(KIND=8)                ::       rmin,rmax
        real(KIND=8)                ::       angle
        character(128)          ::       criteria_name
        ! contains
        !     procedure               ::       limits
    end type region

    contains
        subroutine limits(self,lim)
            implicit none
            type(region)                   ::       self
            real(KIND=8),dimension(1:3,1:2) ::       lim
            real(KIND=8),dimension(1:3)     ::       v,w
            real(KIND=8),dimension(1:8,1:3) ::       box_cyl
            real(KIND=8)                    ::       cone_rmax
            integer                         ::       i

            select case (TRIM(self%name))
            case ('cube')
                lim(1,1) = self%xmin; lim(1,2) = self%xmax
                lim(2,1) = self%ymin; lim(2,2) = self%ymax
                lim(3,1) = self%zmin; lim(3,2) = self%zmax
            case ('sphere')
                lim(1,1) = self%centre(1) - self%rmax; lim(1,2) = self%centre(1) + self%rmax
                lim(2,1) = self%centre(2) - self%rmax; lim(2,2) = self%centre(2) + self%rmax
                lim(3,1) = self%centre(3) - self%rmax; lim(3,2) = self%centre(3) + self%rmax
            case ('cylinder')
                call ortho_find(self%axis,v,w)
                do i=1,3
                    box_cyl(1,i) = self%centre(i) + self%zmax*self%axis(i) + 1.4142136*self%rmax*v(i)
                    box_cyl(2,i) = self%centre(i) + self%zmax*self%axis(i) - 1.4142136*self%rmax*v(i)
                    box_cyl(3,i) = self%centre(i) + self%zmax*self%axis(i) + 1.4142136*self%rmax*w(i)
                    box_cyl(4,i) = self%centre(i) + self%zmax*self%axis(i) - 1.4142136*self%rmax*w(i)
                    box_cyl(5,i) = self%centre(i) + self%zmin*self%axis(i) + 1.4142136*self%rmax*v(i)
                    box_cyl(6,i) = self%centre(i) + self%zmin*self%axis(i) - 1.4142136*self%rmax*v(i)
                    box_cyl(7,i) = self%centre(i) + self%zmin*self%axis(i) + 1.4142136*self%rmax*w(i)
                    box_cyl(8,i) = self%centre(i) + self%zmin*self%axis(i) - 1.4142136*self%rmax*w(i)
                end do
                lim(1,1) = MINVAL(box_cyl(:,1)); lim(1,2) = MAXVAL(box_cyl(:,1))
                lim(2,1) = MINVAL(box_cyl(:,2)); lim(2,2) = MAXVAL(box_cyl(:,2))
                lim(3,1) = MINVAL(box_cyl(:,3)); lim(3,2) = MAXVAL(box_cyl(:,3))
            case ('cone')
                call ortho_find(self%axis,v,w)
                cone_rmax = self%rmax * tan(self%angle)
                do i=1,3
                    box_cyl(1,i) = self%centre(i) + self%rmax*self%axis(i) + 1.4142136*cone_rmax*v(i)
                    box_cyl(2,i) = self%centre(i) + self%rmax*self%axis(i) - 1.4142136*cone_rmax*v(i)
                    box_cyl(3,i) = self%centre(i) + self%rmax*self%axis(i) + 1.4142136*cone_rmax*w(i)
                    box_cyl(4,i) = self%centre(i) + self%rmax*self%axis(i) - 1.4142136*cone_rmax*w(i)
                    box_cyl(5,i) = self%centre(i) + self%rmin*self%axis(i) + 1.4142136*cone_rmax*v(i)
                    box_cyl(6,i) = self%centre(i) + self%rmin*self%axis(i) - 1.4142136*cone_rmax*v(i)
                    box_cyl(7,i) = self%centre(i) + self%rmin*self%axis(i) + 1.4142136*cone_rmax*w(i)
                    box_cyl(8,i) = self%centre(i) + self%rmin*self%axis(i) - 1.4142136*cone_rmax*w(i)
                end do
                lim(1,1) = MINVAL(box_cyl(:,1)); lim(1,2) = MAXVAL(box_cyl(:,1))
                lim(2,1) = MINVAL(box_cyl(:,2)); lim(2,2) = MAXVAL(box_cyl(:,2))
                lim(3,1) = MINVAL(box_cyl(:,3)); lim(3,2) = MAXVAL(box_cyl(:,3))
            case default
                lim(1,1) = 0; lim(1,2) = 1
                lim(2,1) = 0; lim(2,2) = 1
                lim(3,1) = 0; lim(3,2) = 1
            end select
        end subroutine limits
        subroutine setupregion(region_type,region_attrs, reg)
            implicit none
            character(128)           ::       region_type
            real(KIND=8),dimension(1:13) ::       region_attrs
            type(region)                 ::       reg
            
            reg%name = TRIM(region_type)
            select case (reg%name)
            case ('cube')
                reg%xmin = region_attrs(1)
                reg%xmax = region_attrs(2)
                reg%ymin = region_attrs(3)
                reg%ymax = region_attrs(4)
                reg%zmin = region_attrs(5)
                reg%xmax = region_attrs(6)
                reg%axis = (/region_attrs(7),region_attrs(8),region_attrs(9)/)
                reg%bulk_velocity = (/region_attrs(10),region_attrs(11),region_attrs(12)/)
                reg%centre = (/0.5*(reg%xmax-reg%xmin),0.5*(reg%ymax-reg%ymin),0.5*(reg%zmax-reg%zmin)/)
                reg%criteria_name = 'd_euclid'
            case ('sphere')
                reg%centre = (/region_attrs(1),region_attrs(2),region_attrs(3)/)
                reg%rmin = region_attrs(4)
                reg%rmax = region_attrs(5)
                reg%axis = (/region_attrs(6),region_attrs(7),region_attrs(8)/)
                reg%bulk_velocity = (/region_attrs(9),region_attrs(10),region_attrs(11)/)
                reg%criteria_name = 'r_sphere'
            case ('cylinder')
                reg%centre = (/region_attrs(1),region_attrs(2),region_attrs(3)/)
                reg%rmin = region_attrs(4)
                reg%rmax = region_attrs(5)
                reg%zmin = region_attrs(6)
                reg%zmax = region_attrs(7)
                reg%axis = (/region_attrs(8),region_attrs(9),region_attrs(10)/)
                reg%bulk_velocity = (/region_attrs(11),region_attrs(12),region_attrs(13)/)
                reg%criteria_name = 'r_cyl'
            case ('cone')
                reg%centre = (/region_attrs(1),region_attrs(2),region_attrs(3)/)
                reg%rmin = region_attrs(4)
                reg%rmax = region_attrs(5)
                reg%angle = region_attrs(6)
                reg%axis = (/region_attrs(7),region_attrs(8),region_attrs(9)/)
                reg%bulk_velocity = (/region_attrs(10),region_attrs(11),region_attrs(12)/)
                reg%criteria_name = 'r_sphere'
            case default
                write(*,*)'The selected region type is still not supported.'
                write(*,*)'Aborting!'
                stop
            end select
            reg%axis = reg%axis/sqrt(sum(reg%axis**2))
        end subroutine setupregion
        subroutine checkifinside(pos,reg,ok,distance)
            implicit none
            
            real(kind=8),dimension(1:3) ::       pos
            type(region)                ::       reg
            logical                     ::       ok
            real(kind=8)                ::       distance
            if (reg%name.eq.'cube') then
                call cube(pos,reg,ok,distance)
            elseif (reg%name.eq.'sphere') then
                call sphere(pos,reg,ok,distance)
            elseif (reg%name.eq.'cylinder') then
                call cylinder(pos,reg,ok,distance)
            elseif (reg%name.eq.'cone') then
                call cone(pos,reg,ok,distance)
            else
                write(*,*)'The selected region type is still not supported.'
                write(*,*)'Aborting!'
                stop
            endif
        end subroutine
        subroutine cube(pos,reg,ok,distance)
            implicit none
            
            real(kind=8),dimension(1:3) ::       pos
            type(region)                ::       reg
            logical                     ::       ok
            real(kind=8)                ::       distance
            ok = (reg%xmin <= pos(1) <= reg%xmax.and.&
                    &reg%ymin <= pos(2) <= reg%ymax.and.&
                    &reg%zmin <= pos(3) <= reg%zmax)
            if (ok) then
                distance=sqrt(sum(pos**2))
            endif
        end subroutine cube

        subroutine sphere(pos,reg,ok,distance)
            implicit none
            real(kind=8),dimension(1:3) ::       pos
            type(region)                ::       reg
            logical                     ::       ok
            real(kind=8)                ::       distance

            distance = sqrt(sum(pos**2))
            ok = (reg%rmin <= distance <= reg%rmax)
        end subroutine sphere

        subroutine cylinder(pos,reg,ok,distance)
            implicit none
            real(kind=8),dimension(1:3) ::       pos
            type(region)                ::       reg
            logical                     ::       ok
            real(kind=8)                ::       distance

            distance = sqrt(pos(1)**2 + pos(2)**2)
            ok = (reg%rmin <= distance <= reg%rmax.and.reg%zmin <= pos(3) <= reg%zmax)
        end subroutine cylinder

        subroutine cone(pos,reg,ok,distance)
            implicit none
            real(kind=8),dimension(1:3) ::       pos
            type(region)                ::       reg
            logical                     ::       ok
            real(kind=8)                ::       distance,theta

            distance = dot_product(pos,reg%axis)
            theta = acos(distance/sqrt(sum(pos**2)))
            ok = (reg%rmin <= distance <= reg%rmax.and.reg%angle >= abs(theta))
        end subroutine cone

end module geometrical_regions