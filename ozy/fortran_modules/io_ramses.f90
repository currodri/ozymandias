module io_ramses
    public

    type hydroID
        integer :: nvar
        integer :: density,vx,vy,vz,thermal_pressure,metallicity
        integer :: Blx,Bly,Blz,Brx,Bry,Brz
        integer :: eCR
        integer :: xHII,xHeII,xHeIII
    end type hydroID
    

    contains
        subroutine title(n,nchar)
        !=======================================================================
            implicit none
            integer::n
            character(5)::nchar

            character::nchar1
            character(2)::nchar2
            character(3)::nchar3
            character(4)::nchar4
            character(5)::nchar5

            if(n.ge.10000)then
                write(nchar5,'(i5)') n
                nchar = nchar5
            elseif(n.ge.1000)then
                write(nchar4,'(i4)') n
                nchar = '0'//nchar4
            elseif(n.ge.100)then
                write(nchar3,'(i3)') n
                nchar = '00'//nchar3
            elseif(n.ge.10)then
                write(nchar2,'(i2)') n
                nchar = '000'//nchar2
            else
                write(nchar1,'(i1)') n
                nchar = '0000'//nchar1
            endif
        end subroutine title
! subroutine hilbert3d(x,y,z,order,bit_length,npoint)
!   implicit none

!   integer     ,INTENT(IN)                     ::bit_length,npoint
!   integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
!   real(kind=8),INTENT(OUT),dimension(1:npoint)::order

!   logical,dimension(0:3*bit_length-1)::i_bit_mask
!   logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
!   integer,dimension(0:7,0:1,0:11)::state_diagram
!   integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

!   if(bit_length>bit_size(bit_length))then
!      write(*,*)'Maximum bit length=',bit_size(bit_length)
!      write(*,*)'stop in hilbert3d'
!      stop
!   endif

!   state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
!                             &   0, 1, 3, 2, 7, 6, 4, 5,&
!                             &   2, 6, 0, 7, 8, 8, 0, 7,&
!                             &   0, 7, 1, 6, 3, 4, 2, 5,&
!                             &   0, 9,10, 9, 1, 1,11,11,&
!                             &   0, 3, 7, 4, 1, 2, 6, 5,&
!                             &   6, 0, 6,11, 9, 0, 9, 8,&
!                             &   2, 3, 1, 0, 5, 4, 6, 7,&
!                             &  11,11, 0, 7, 5, 9, 0, 7,&
!                             &   4, 3, 5, 2, 7, 0, 6, 1,&
!                             &   4, 4, 8, 8, 0, 6,10, 6,&
!                             &   6, 5, 1, 2, 7, 4, 0, 3,&
!                             &   5, 7, 5, 3, 1, 1,11,11,&
!                             &   4, 7, 3, 0, 5, 6, 2, 1,&
!                             &   6, 1, 6,10, 9, 4, 9,10,&
!                             &   6, 7, 5, 4, 1, 0, 2, 3,&
!                             &  10, 3, 1, 1,10, 3, 5, 9,&
!                             &   2, 5, 3, 4, 1, 6, 0, 7,&
!                             &   4, 4, 8, 8, 2, 7, 2, 3,&
!                             &   2, 1, 5, 6, 3, 0, 4, 7,&
!                             &   7, 2,11, 2, 7, 5, 8, 5,&
!                             &   4, 5, 7, 6, 3, 2, 0, 1,&
!                             &  10, 3, 2, 6,10, 3, 4, 4,&
!                             &   6, 1, 7, 0, 5, 2, 4, 3 /), &
!                             & (/8 ,2, 12 /) )

!   do ip=1,npoint

!      ! convert to binary
!      do i=0,bit_length-1
!         x_bit_mask(i)=btest(x(ip),i)
!         y_bit_mask(i)=btest(y(ip),i)
!         z_bit_mask(i)=btest(z(ip),i)
!      enddo

!      ! interleave bits
!      do i=0,bit_length-1
!         i_bit_mask(3*i+2)=x_bit_mask(i)
!         i_bit_mask(3*i+1)=y_bit_mask(i)
!         i_bit_mask(3*i  )=z_bit_mask(i)
!      end do

!      ! build Hilbert ordering using state diagram
!      cstate=0
!      do i=bit_length-1,0,-1
!         b2=0 ; if(i_bit_mask(3*i+2))b2=1
!         b1=0 ; if(i_bit_mask(3*i+1))b1=1
!         b0=0 ; if(i_bit_mask(3*i  ))b0=1
!         sdigit=b2*4+b1*2+b0
!         nstate=state_diagram(sdigit,0,cstate)
!         hdigit=state_diagram(sdigit,1,cstate)
!         i_bit_mask(3*i+2)=btest(hdigit,2)
!         i_bit_mask(3*i+1)=btest(hdigit,1)
!         i_bit_mask(3*i  )=btest(hdigit,0)
!         cstate=nstate
!      enddo

!      ! save Hilbert key as double precision real
!      order(ip)=0.
!      do i=0,3*bit_length-1
!         b0=0 ; if(i_bit_mask(i))b0=1
!         order(ip)=order(ip)+dble(b0)*dble(2)**i
!      end do

!   end do

! end subroutine hilbert3d
        subroutine check_lmax(ngridfile,ncpu,nboundary,nlevelmax,lmax)
            ! This simple subroutines checks for the actual maximum
            ! level of refinement active in thAlgo e simulation.
            implicit none

            integer :: ncpu,nboundary,nlevelmax,lmax
            integer,dimension(1:ncpu+nboundary,1:nlevelmax) :: ngridfile
            integer :: ngridilevel,i

            do i = nlevelmax, 0, -1
                ngridilevel=sum(ngridfile(:,i))
                if (ngridilevel .gt. 0) then
                    if (lmax .gt. i) then
                        lmax=i
                    endif
                    exit
                endif
            end do

        end subroutine check_lmax
        subroutine read_hydrofile_descriptor(repository,varIDs)
            implicit none

            character(128) ::  repository,nomfich
            logical            ::  ok
            integer            ::  nvar,i
            character(25)  ::  newVar
            character(9)   ::  igr9
            character(8)   ::  igr8
            character   ::  igr1
            character(2)   ::  igr2
            character(3)   ::  igr3
            integer            ::  newID,statn

            type(hydroID)      ::  varIDs

            nomfich=TRIM(repository)//'/hydro_file_descriptor.txt'
            inquire(file=nomfich, exist=ok) ! verify input file
            if ( ok ) then
                write(*,'(": Reading variables IDs from hydro_descriptor")')
                open(unit=10,file=nomfich,status='old',form='formatted')
                read(10,'("nvar        =",I11)')nvar
                varIDs%nvar = nvar
                do i=1,nvar
                    read(10,*) igr9,igr1,igr1,newVar
                    read(igr1,*,iostat=statn) newID
                    call select_from_descriptor_IDs(varIDs,newVar,newID)
                end do
                close(10)
            else
                write(*,'(": ",A," not found. Initializing variables to default IDs.")') trim(nomfich)
                varIDs%density = 1
                varIDs%vx = 2; varIDs%vy  = 3; varIDs%vz  = 4
                varIDs%Blx  = 5; varIDs%Bly = 6; varIDs%Blz = 7
                varIDs%Brx  = 8; varIDs%Bry = 9; varIDs%Brz = 10
                varIDs%thermal_pressure = 11; varIDs%metallicity = 12;
                varIDs%xHII = 13; varIDs%xHeII = 14; varIDs%xHeIII = 15
            end if
        end subroutine read_hydrofile_descriptor

        subroutine select_from_descriptor_IDs(varIDs,newVar,newID)
            implicit none
            integer,intent(in)           :: newID
            character(25),intent(in) :: newvar
            type(hydroID)                :: varIDs
            select case (TRIM(newVar))
            case ('density')
                varIDs%density = newID
            case ('velocity_x')
                varIDs%vx = newID
            case ('velocity_y')
                varIDs%vy = newID
            case ('velocity_z')
                varIDs%vz = newID
            case ('B_left_x')
                varIDs%Blx = newID
            case ('B_left_y')
                varIDs%Bly = newID
            case ('B_left_z')
                varIDs%Blz = newID
            case ('B_right_x')
                varIDs%Brx = newID
            case ('B_right_y')
                varIDs%Bry = newID
            case ('B_right_z')
                varIDs%Brz = newID
            case ('thermal_pressure')
                varIDs%thermal_pressure = newID
            case ('non_thermal_pressure_1')
                write(*,'(": Using non_thermal_pressure_1 as cosmic rays energy density (variable ",I2,")")') newID
                varIDs%eCR = newID
            case ('passive_scalar_1')
                write(*,'(": Using passive_scalar_1 as metallicity (variable ",I2,")")') newID
                varIDs%metallicity = newID
            ! if (pasrefine) then
            !     write(*,'(": Using passive_scalar_1 as refine field (variable ",I2,")")') newID
            !     refineID = newID
            ! else
            !     write(*,'(": Using passive_scalar_1 as metallicity (variable ",I2,")")') newID
            !     metalsID = newID
            !     has_metals=.true.
            ! end if
            ! case ('passive_scalar_2')
            ! if (pasrefine) then
            !     write(*,'(": Using passive_scalar_2 as metallicity (variable ",I2,")")') newID
            !     metalsID = newID
            !     has_metals=.true.
            ! else
            !     if (magTracers.eq.0) then
            !     write(*,'(": Using passive_scalar_2 as xHII (variable ",I2,")")') newID
            !     xHII_ID = newID
            !     has_rt=.true.
            ! end if
            ! end if
            ! case ('passive_scalar_3')
            ! if (pasrefine) then
            !     write(*,'(": Using passive_scalar_3 as xHII (variable ",I2,")")') newID
            !     xHII_ID = newID
            !     has_rt=.true.
            ! else
            !     if (magTracers.eq.0) then
            !     write(*,'(": Using passive_scalar_3 as xHeII (variable ",I2,")")') newID
            !     xHeII_ID = newID
            !     has_rt=.true.
            !     end if
            ! end if
            ! case ('passive_scalar_4')
            ! if (pasrefine) then
            !     write(*,'(": Using passive_scalar_4 as xHeII (variable ",I2,")")') newID
            !     xHeII_ID = newID
            !     has_rt=.true.
            ! else
            !     if (magTracers.eq.0) then
            !     write(*,'(": Using passive_scalar_4 as xHeIII (variable ",I2,")")') newID
            !     xHeIII_ID = newID
            !     has_rt=.true.
            !     end if
            ! end if
            ! case ('passive_scalar_5')
            ! if (pasrefine) then
            !     write(*,'(": Using passive_scalar_4 as xHeIII (variable ",I2,")")') newID
            !     xHeIII_ID = newID
            !     has_rt=.true.
            ! else
            !     write(*,'(": I have no idea what passive_scalar 4 is (variable ",I2,")")') newID
            ! end if
            case ('xHII')
                varIDs%xHII = newID
            case ('xHeII')
                varIDs%xHeII = newID
            case ('xHeIII')
                varIDs%xHeIII = newID
            end select
        end subroutine select_from_descriptor_IDs

        subroutine getvarvalue(varIDs,reg,dx,x,var,varname,value)
            use math_utils
            use geometrical_regions
            implicit none
            type(hydroID)                      :: varIDs
            type(region)                       :: reg
            real(KIND=8)                       :: dx
            real(KIND=8),dimension(1:3)        :: x,v_corrected,L,B,e_theta,e_r,e_phi,e_z
            real(KIND=8),dimension(1:varIDs%nvar) :: var
            character(128)                 :: varname
            real(KIND=8)                       :: value

            select case (TRIM(varname))
            case ('d_euclid')
                value = sqrt(sum(x**2))
            case ('r_sphere')
                value = r_sphere(x)
            case ('theta_sphere')
                value = theta_sphere(x)
            case ('phi_sphere')
                value = phi_sphere(x)
            case('r_cyl')
                write(*,*)x(1),x(2),x(3)
                value = r_cyl(x)
            case ('phi_cyl')
                value = phi_cyl(x)
            case ('v_rot')
                v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/) - reg%bulk_velocity
                e_phi = cross_product(x, (/0d0,0d0,1d0/))
                e_phi = e_theta/sqrt(sum(e_theta**2))
                value = dot_product(v_corrected,e_phi)
            case ('v_sphere_r')
                v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/) - reg%bulk_velocity
                call sphere_basis(x,e_r,e_theta,e_phi)
                value = dot_product(v_corrected,e_r)
            case ('v_sphere_phi')
                v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/) - reg%bulk_velocity
                call sphere_basis(x,e_r,e_theta,e_phi)
                value = dot_product(v_corrected,e_phi)
            case ('v_sphere_theta')
                v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/) - reg%bulk_velocity
                call sphere_basis(x,e_r,e_theta,e_phi)
                value = dot_product(v_corrected,e_theta)
            case ('v_cyl_r')
                v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/) - reg%bulk_velocity
                call cyl_basis(x,e_r,e_z,e_phi)
                value = dot_product(v_corrected,e_r)
            case ('v_cyl_z')
                v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/) - reg%bulk_velocity
                call cyl_basis(x,e_r,e_z,e_phi)
                value = dot_product(v_corrected,e_z)
            case ('v_cyl_phi')
                v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/) - reg%bulk_velocity
                call cyl_basis(x,e_r,e_z,e_phi)
                value = dot_product(v_corrected,e_phi)
            case ('density')
                value = var(varIDs%density)
            case ('mass')
                value = (var(varIDs%density) * (dx*dx)) * dx
            case ('volume')
                value = (dx * dx) * dx
            case ('metallicity')
                value = var(varIDs%metallicity)
            case ('thermal_pressure')
                value = var(varIDs%thermal_pressure)
            case ('thermal_energy')
                value = ((5d0/3d0 - 1d0) * var(varIDs%thermal_pressure) * (dx * dx)) * dx
            case ('kinetic_energy')
                value = (0.5 * (var(varIDs%density) * (dx*dx)) * dx) * sqrt(var(varIDs%vx)**2 + var(varIDs%vy)**2 + var(varIDs%vz)**2)
            case ('magnetic_energy')
                B = 0.5 *(/(var(varIDs%Blx)-var(varIDs%Brx)),(var(varIDs%Bly)-var(varIDs%Bry)),(var(varIDs%Blz)-var(varIDs%Brz))/)
                value = (0.5 * sum(B**2) * (dx*dx)) * dx
            case ('cr_energy')
                value = (var(varIDs%eCR) * (dx*dx)) * dx
            case ('momentum_x')
                value = ((var(varIDs%density) * (dx*dx)) * dx) * (var(varIDs%vx) - reg%bulk_velocity(1))
            case ('momentum_y')
                value = ((var(varIDs%density) * (dx*dx)) * dx) * (var(varIDs%vy) - reg%bulk_velocity(2))
            case ('momentum_z')
                value = ((var(varIDs%density) * (dx*dx)) * dx) * (var(varIDs%vz) - reg%bulk_velocity(3))
            case ('momentum')
                v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/) - reg%bulk_velocity
                value = ((var(varIDs%density) * (dx*dx)) * dx) * &
                        &sqrt(sum(v_corrected**2))
            case ('ang_momentum_x')
                v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/) - reg%bulk_velocity
                value = ((var(varIDs%density) * (dx*dx)) * dx) * (x(2)*v_corrected(3) &
                            &- v_corrected(2)*x(3))
            case ('ang_momentum_y')
                v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/) - reg%bulk_velocity
                value = ((var(varIDs%density) * (dx*dx)) * dx) * (x(3)*v_corrected(1) &
                            &- v_corrected(3)*x(1))
            case ('ang_momentum_z')
                v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/) - reg%bulk_velocity
                value = ((var(varIDs%density) * (dx*dx)) * dx) * (x(1)*v_corrected(2) &
                            &- v_corrected(1)*x(2))
            case ('ang_momentum')
                v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/) - reg%bulk_velocity
                L = cross_product(x,v_corrected)
                value = ((var(varIDs%density) * (dx*dx)) * dx) * sqrt(sum(L**2))
            end select

        end subroutine getvarvalue
end module io_ramses