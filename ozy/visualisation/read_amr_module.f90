!--------------------------------------------------------------------------
! ozymandias:read_amr_module.f90
!--------------------------------------------------------------------------
!
! MODULE: io_ramses
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> types and routines useful for the reading of AMR and HYDRO data from
!> RAMSES.
!
!> @details  
!> define hydroID, amr_info, sim_info and level derived types, as well
!> as functions for the initial setup of AMR and HYDRO data.
! 
!
!> @date 29/3/2021   0.1 fortran routines for ozymandias
!--------------------------------------------------------------------------
module io_ramses
    use local

    type hydroID
        integer :: nvar
        integer :: density,vx,vy,vz,thermal_pressure,metallicity
        integer :: Blx,Bly,Blz,Brx,Bry,Brz
        integer :: eCR
        integer :: xHII,xHeII,xHeIII
    end type hydroID

    type amr_info
        integer :: ncpu,ndim,nlevelmax,nboundary,twotondim,ndom
        integer :: levelmin,levelmax,lmax
        integer :: ncpu_read
        character(80) :: ordering
        integer,dimension(:),allocatable :: cpu_list
        real(dbl),dimension(:),allocatable :: bound_key
        logical,dimension(:),allocatable :: cpu_read
        real(dbl),dimension(1:3) :: xbound=(/0d0,0d0,0d0/)
    end type amr_info

    type sim_info
        real(sgl) :: t,aexp,omega_m,omega_l,omega_k,omega_b
        real(dbl) :: h0,unit_l,unit_d,unit_t
    end type sim_info

    type level
        integer::ilevel
        integer::ngrid
        real(sgl),dimension(:,:,:),pointer::cube
        real(dbl),dimension(:,:),pointer::map
        real(dbl),dimension(:,:),pointer::rho
        integer::imin
        integer::imax
        integer::jmin
        integer::jmax
        integer::kmin
        integer::kmax
    end type level
    
    type data_handler
        character(80) :: name
        logical :: x_data=.false.,y_data=.false.,z_data=.false.
        integer,dimension(2) :: nx,ny,nz
        real(dbl),dimension(:,:,:),allocatable :: x,y,z
    end type data_handler

    contains

    !---------------------------------------------------------------
    ! Subroutine: TITLE
    !
    ! # From RAMSES utils.
    !---------------------------------------------------------------
    subroutine title(n,nchar)
        implicit none
        integer, intent(in) :: n
        character(5),intent(inout) :: nchar

        character :: nchar1
        character(2) :: nchar2
        character(3) :: nchar3
        character(4) :: nchar4
        character(5) :: nchar5

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
    !---------------------------------------------------------------
    ! Subroutine: HILBERT 3D CURVE
    !
    ! # From RAMSES utils.
    !---------------------------------------------------------------
    subroutine hilbert3d(x,y,z,order,bit_length,npoint)
        implicit none
      
        integer,intent(in)                    :: bit_length,npoint
        integer,dimension(1:npoint),intent(in) :: x,y,z
        real(dbl),dimension(1:npoint),intent(out) :: order
      
        logical,dimension(0:3*bit_length-1) :: i_bit_mask
        logical,dimension(0:1*bit_length-1) :: x_bit_mask,y_bit_mask,z_bit_mask
        integer,dimension(0:7,0:1,0:11) :: state_diagram
        integer :: i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit
      
        if(bit_length>bit_size(bit_length))then
           write(*,*)'Maximum bit length=',bit_size(bit_length)
           write(*,*)'stop in hilbert3d'
           stop
        endif
      
        state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                                  &   0, 1, 3, 2, 7, 6, 4, 5,&
                                  &   2, 6, 0, 7, 8, 8, 0, 7,&
                                  &   0, 7, 1, 6, 3, 4, 2, 5,&
                                  &   0, 9,10, 9, 1, 1,11,11,&
                                  &   0, 3, 7, 4, 1, 2, 6, 5,&
                                  &   6, 0, 6,11, 9, 0, 9, 8,&
                                  &   2, 3, 1, 0, 5, 4, 6, 7,&
                                  &  11,11, 0, 7, 5, 9, 0, 7,&
                                  &   4, 3, 5, 2, 7, 0, 6, 1,&
                                  &   4, 4, 8, 8, 0, 6,10, 6,&
                                  &   6, 5, 1, 2, 7, 4, 0, 3,&
                                  &   5, 7, 5, 3, 1, 1,11,11,&
                                  &   4, 7, 3, 0, 5, 6, 2, 1,&
                                  &   6, 1, 6,10, 9, 4, 9,10,&
                                  &   6, 7, 5, 4, 1, 0, 2, 3,&
                                  &  10, 3, 1, 1,10, 3, 5, 9,&
                                  &   2, 5, 3, 4, 1, 6, 0, 7,&
                                  &   4, 4, 8, 8, 2, 7, 2, 3,&
                                  &   2, 1, 5, 6, 3, 0, 4, 7,&
                                  &   7, 2,11, 2, 7, 5, 8, 5,&
                                  &   4, 5, 7, 6, 3, 2, 0, 1,&
                                  &  10, 3, 2, 6,10, 3, 4, 4,&
                                  &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                                  & (/8 ,2, 12 /) )
      
        do ip=1,npoint
      
           ! convert to binary
           do i=0,bit_length-1
              x_bit_mask(i)=btest(x(ip),i)
              y_bit_mask(i)=btest(y(ip),i)
              z_bit_mask(i)=btest(z(ip),i)
           enddo
      
           ! interleave bits
           do i=0,bit_length-1
              i_bit_mask(3*i+2)=x_bit_mask(i)
              i_bit_mask(3*i+1)=y_bit_mask(i)
              i_bit_mask(3*i  )=z_bit_mask(i)
           end do
      
           ! build Hilbert ordering using state diagram
           cstate=0
           do i=bit_length-1,0,-1
              b2=0 ; if(i_bit_mask(3*i+2))b2=1
              b1=0 ; if(i_bit_mask(3*i+1))b1=1
              b0=0 ; if(i_bit_mask(3*i  ))b0=1
              sdigit=b2*4+b1*2+b0
              nstate=state_diagram(sdigit,0,cstate)
              hdigit=state_diagram(sdigit,1,cstate)
              i_bit_mask(3*i+2)=btest(hdigit,2)
              i_bit_mask(3*i+1)=btest(hdigit,1)
              i_bit_mask(3*i  )=btest(hdigit,0)
              cstate=nstate
           enddo
      
           ! save Hilbert key as double precision real
           order(ip)=0.
           do i=0,3*bit_length-1
              b0=0 ; if(i_bit_mask(i))b0=1
              order(ip)=order(ip)+dble(b0)*dble(2)**i
           end do
      
        end do
      
    end subroutine hilbert3d
    !---------------------------------------------------------------
    ! Subroutine: CHECK ACTIVE LEVELMAX
    !
    ! This simple subroutines checks for the actual maximum
    ! level of refinement active in the simulation.
    !---------------------------------------------------------------
    subroutine check_lmax(ngridfile,amr)

        implicit none

        type(amr_info),intent(inout) :: amr
        integer,dimension(1:amr%ncpu+amr%nboundary,1:amr%nlevelmax),intent(in) :: ngridfile
        integer :: ngridilevel,i

        do i = amr%nlevelmax, 0, -1
            ngridilevel=sum(ngridfile(:,i))
            if (ngridilevel.gt.0) then
                if (amr%lmax.gt.i) then
                    amr%lmax=i
                endif
                exit
            endif
        end do
    end subroutine check_lmax

    !---------------------------------------------------------------
    ! Subroutine: READ HYDRO IDs
    !
    ! This routine extracts the HYDRO variable IDs for a given
    ! simulation snapshot such that their order is known when
    ! hydro_*.out* files are read.
    !---------------------------------------------------------------    
    subroutine read_hydrofile_descriptor(repository,varIDs)
        implicit none

        character(128),intent(in) ::  repository
        type(hydroID),intent(inout)      ::  varIDs
        character(128) :: nomfich
        logical            ::  ok
        integer            ::  nvar,i
        character(25)  ::  newVar
        character(9)   ::  igr9
        character   ::  igr1
        integer            ::  newID,statn
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
        type(hydroID),intent(inout)                :: varIDs
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

    !---------------------------------------------------------------
    ! Subroutine: GET VARIABLE VALUE
    !
    ! For a given cell, a particular variable is computed.
    ! This is a very important subroutine, which will be modified
    ! extensively.
    !---------------------------------------------------------------
    subroutine getvarvalue(varIDs,reg,dx,x,var,varname,value)
        use vectors
        use basis_representations
        use geometrical_regions
        use coordinate_systems
        implicit none
        type(hydroID),intent(in)                      :: varIDs
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(1:varIDs%nvar),intent(in) :: var
        character(128),intent(in)                 :: varname
        real(dbl),intent(inout)                       :: value
        type(vector) :: v_corrected,L,B
        type(basis) :: temp_basis

        select case (TRIM(varname))
        case ('d_euclid')
            ! Euclidean distance
            value = magnitude(x)
        case ('r_sphere')
            ! Radius from center of sphere
            value = r_sphere(x)
        case ('theta_sphere')
            ! Value of spherical theta angle measured from the z axis
            value = theta_sphere(x)
        case ('phi_sphere')
            ! Value of spherical phi angle measure in the x-y plane 
            ! from the x axis
            value = phi_sphere(x)
        case('r_cyl')
            ! Value of cylindrical radius
            value = r_cyl(x)
        case ('phi_cyl')
            ! Value of spherical phi angle measure in the x-y plane 
            ! from the x axis
            value = phi_cyl(x)
        case ('v_sphere_r')
            ! Velocity component in the spherical radial direction
            ! 1. Correct velocity for bulk velocity of region
            ! 2. Dot product of velocity vector with spherical radial
            !    unit vector
            v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/)
            v_corrected = v_corrected - reg%bulk_velocity
            call spherical_basis_from_cartesian(x,temp_basis)
            value = v_corrected.DOT.temp_basis%u(1)
        case ('v_sphere_phi')
            ! Velocity component in the spherical azimutal (phi) direction
            ! 1. Correct velocity for bulk velocity of region
            ! 2. Dot product of velocity vector with spherical phi
            !    unit vector
            v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/)
            v_corrected = v_corrected - reg%bulk_velocity
            call spherical_basis_from_cartesian(x,temp_basis)
            value = v_corrected .DOT. temp_basis%u(3)
        case ('v_sphere_theta')
            ! Velocity component in the spherical theta direction
            ! 1. Correct velocity for bulk velocity of region
            ! 2. Dot product of velocity vector with spherical theta
            !    unit vector
            v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/)
            v_corrected = v_corrected - reg%bulk_velocity
            call spherical_basis_from_cartesian(x,temp_basis)
            value = v_corrected .DOT. temp_basis%u(2)
        case ('v_cyl_r')
            ! Velocity component in the cylindrical radial direction
            ! 1. Correct velocity for bulk velocity of region
            ! 2. Dot product of velocity vector with cylindrical
            !    radial unit vector
            v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/)
            v_corrected = v_corrected - reg%bulk_velocity
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = v_corrected.DOT.temp_basis%u(1)
        case ('v_cyl_z')
            ! Velocity component in the cylindrical z direction
            ! 1. Correct velocity for bulk velocity of region
            ! 2. Dot product of velocity vector with cylindrical
            !    z unit vector
            v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/)
            v_corrected = v_corrected - reg%bulk_velocity
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = v_corrected.DOT.temp_basis%u(3)
        case ('v_cyl_phi')
            ! Velocity component in the cylyndrical azimutal (phi) direction
            ! 1. Correct velocity for bulk velocity of region
            ! 2. Dot product of velocity vector with cylindrical
            !    phi unit vector
            v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/)
            v_corrected = v_corrected - reg%bulk_velocity
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = v_corrected.DOT.temp_basis%u(2)
        case ('density')
            ! Density
            value = var(varIDs%density)
        case ('mass')
            ! Total mass, computed as density times volume of cell (dx**3)
            value = (var(varIDs%density) * (dx*dx)) * dx
        case ('volume')
            ! Volume of cell (dx**3)
            value = (dx * dx) * dx
        case ('metallicity')
            ! Metallicity
            value = var(varIDs%metallicity)
        case ('temperature')
            ! Gas temperature
            value = var(varIDs%thermal_pressure) / var(varIDs%density) / 1.38d-16*1.66d-24
        case ('thermal_pressure')
            ! Thermal pressure
            value = var(varIDs%thermal_pressure)
        case ('thermal_energy')
            ! Thermal energy, computed as (gamma - 1)*thermal_pressure*volume
            value = ((5D0/3d0 - 1d0) * var(varIDs%thermal_pressure) * (dx * dx)) * dx
        case ('kinetic_energy')
            ! Kinetic energy, computed as 1/2*density*volume*magnitude(velocity)
            ! DISCLAIMER: Velocity not corrected for bulk velocity
            value = (0.5 * (var(varIDs%density) * (dx*dx)) * dx) * sqrt(var(varIDs%vx)**2 + var(varIDs%vy)**2 + var(varIDs%vz)**2)
        case ('magnetic_energy')
            ! Magnetic energy as magnitude(B)**2/2
            B = 0.5 *(/(var(varIDs%Blx)-var(varIDs%Brx)),(var(varIDs%Bly)-var(varIDs%Bry)),(var(varIDs%Blz)-var(varIDs%Brz))/)
            value = (0.5 * (B.DOT.B) * (dx*dx)) * dx
        case ('cr_energy')
            ! CR energy, computed as CR_energydensity*volume
            value = (var(varIDs%eCR) * (dx*dx)) * dx
        case ('momentum_x')
            ! Linear momentum in the x direction as density*volume*corrected_velocity_x
            value = ((var(varIDs%density) * (dx*dx)) * dx) * (var(varIDs%vx) - reg%bulk_velocity%x)
        case ('momentum_y')
            ! Linear momentum in the y direction density*volume*corrected_velocity_y
            value = ((var(varIDs%density) * (dx*dx)) * dx) * (var(varIDs%vy) - reg%bulk_velocity%y)
        case ('momentum_z')
            ! Linear momentum in the z direction density*volume*corrected_velocity_z
            value = ((var(varIDs%density) * (dx*dx)) * dx) * (var(varIDs%vz) - reg%bulk_velocity%z)
        case ('momentum')
            ! Magnitude of linear momentum, using corrected velocity
            v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/)
            v_corrected = v_corrected - reg%bulk_velocity
            value = ((var(varIDs%density) * (dx*dx)) * dx) * &
                    & magnitude(v_corrected)
        case ('momentum_sphere_r')
            ! Linear momentum in the spherical radial direction
            ! 1. Correct velocity for bulk velocity of region
            ! 2. Dot product of velocity vector with spherical phi
            !    unit vector
            ! 3. Multiply by mass of cell
            v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/)
            v_corrected = v_corrected - reg%bulk_velocity
            call spherical_basis_from_cartesian(x,temp_basis)
            value = (var(varIDs%density) * (dx*dx)) * dx * (v_corrected .DOT. temp_basis%u(3))
        case ('ang_momentum_x')
            ! Corrected angular momentum in the x direction
            v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/)
            v_corrected = v_corrected - reg%bulk_velocity
            value = ((var(varIDs%density) * (dx*dx)) * dx) * (x%y * v_corrected%z &
                        &- v_corrected%y * x%z)
        case ('ang_momentum_y')
            ! Corrected angular momentum in the y direction
            v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/)
            v_corrected = v_corrected - reg%bulk_velocity
            value = ((var(varIDs%density) * (dx*dx)) * dx) * (x%z*v_corrected%x &
                        &- v_corrected%z*x%x)
        case ('ang_momentum_z')
            ! Corrected angular momentum in the z direction
            v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/)
            v_corrected = v_corrected - reg%bulk_velocity
            value = ((var(varIDs%density) * (dx*dx)) * dx) * (x%x*v_corrected%y &
                        &- v_corrected%x*x%y)
        case ('ang_momentum')
            ! Corrected magnitude of angular momentum
            v_corrected = (/var(varIDs%vx),var(varIDs%vy),var(varIDs%vz)/)
            v_corrected = v_corrected - reg%bulk_velocity
            L = x * v_corrected
            value = ((var(varIDs%density) * (dx*dx)) * dx) * magnitude(L)
        case default
            write(*,*)'Variable not supported: ',TRIM(varname)
            write(*,*)'Aborting!'
            stop
        end select

    end subroutine getvarvalue

    !---------------------------------------------------------------
    ! Subroutine: INITIAL AMR SETUP
    !
    ! This routine has been adapted from the RAMSES utils, in which
    ! the initial read of the amr data is performed and saved to 
    ! amr_info and sim_info derived types.
    !---------------------------------------------------------------
    subroutine init_amr_read(repository,amr,sim)
        implicit none
        character(128), intent(in) :: repository
        type(amr_info), intent(inout) :: amr
        type(sim_info), intent(inout) :: sim
        character(5) :: nchar
        integer :: ipos,impi,i,nx,ny,nz
        character(128) :: nomfich
        logical :: ok

        ipos = index(repository,'output_')
        nchar=repository(ipos+7:ipos+13)
        nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
        ! Verify input hydro file
        inquire(file=nomfich, exist=ok)
        if ( .not. ok ) then
            print *,TRIM(nomfich)//' not found.'
            stop
        endif
        nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
        ! Verify input amr file
        inquire(file=nomfich, exist=ok)
        if ( .not. ok ) then
            print *,TRIM(nomfich)//' not found.'
            stop
        endif

        nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
        open(unit=10,file=nomfich,status='old',form='unformatted')
        read(10)amr%ncpu
        read(10)amr%ndim
        read(10)nx,ny,nz
        read(10)amr%nlevelmax
        read(10)!ngridmax
        read(10)amr%nboundary
        read(10)!ngrid_current
        read(10)!boxlen
        close(10)
        amr%xbound = (/dble(nx/2),dble(ny/2),dble(nz/2)/)
        amr%twotondim = 2**amr%ndim

        if(amr%ndim==2)then
            write(*,*)'Output file contains 2D data'
            write(*,*)'Aborting!'
            stop
        endif

        nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
        open(unit=10,file=nomfich,form='formatted',status='old')
        read(10,*)!ncpu
        read(10,*)!ndim
        read(10,'("levelmin    =",I11)')amr%levelmin
        read(10,'("levelmax    =",I11)')amr%levelmax
        read(10,*)
        read(10,*)
        read(10,*)
    
        read(10,*)
        read(10,'("time        =",E23.15)')sim%t
        read(10,'("aexp        =",E23.15)')sim%aexp
        read(10,'("H0          =",E23.15)')sim%h0
        read(10,'("omega_m     =",E23.15)')sim%omega_m
        read(10,'("omega_l     =",E23.15)')sim%omega_l
        read(10,'("omega_k     =",E23.15)')sim%omega_k
        read(10,'("omega_b     =",E23.15)')sim%omega_b
        read(10,'("unit_l      =",E23.15)')sim%unit_l
        read(10,'("unit_d      =",E23.15)')sim%unit_d
        read(10,'("unit_t      =",E23.15)')sim%unit_t
        read(10,*)
    
        read(10,'("ordering type=",A80)')amr%ordering
        read(10,*)
        allocate(amr%cpu_list(1:amr%ncpu))
        if(TRIM(amr%ordering).eq.'hilbert')then
            allocate(amr%bound_key(0:amr%ncpu))
            allocate(amr%cpu_read(1:amr%ncpu))
            amr%cpu_read=.false.
            do impi=1,amr%ncpu
                read(10,'(I8,1X,E23.15,1X,E23.15)')i,amr%bound_key(impi-1),amr%bound_key(impi)
            end do
        endif
        close(10)
    end subroutine init_amr_read

    !---------------------------------------------------------------
    ! Subroutine: COMPUTE CPU MAP
    !
    ! This has been adapted from the RAMSES utils, in which the cpu
    ! map is obtained for a desired region.
    !---------------------------------------------------------------
    subroutine get_cpu_map(reg,amr)
        use geometrical_regions
        implicit none
        type(region),intent(in) :: reg
        type(amr_info),intent(inout) :: amr
        real(dbl), dimension(1:3,1:2) :: box_limits
        real(dbl) :: xxmin,xxmax,yymin,yymax,zzmin,zzmax
        real(dbl) :: dxmax,dx,dkey
        real(dbl),dimension(1:1) :: order_min
        integer :: i,impi,j,ilevel,lmin,imin,imax,jmin,jmax,kmin,kmax
        integer :: bit_length,maxdom,ndom=1
        real(dbl),dimension(1:8) :: bounding_min,bounding_max
        integer,dimension(1:8) :: idom,jdom,kdom,cpu_min,cpu_max

        call limits(reg,box_limits)
        xxmin=box_limits(1,1) ; xxmax=box_limits(1,2)
        yymin=box_limits(2,1) ; yymax=box_limits(2,2)
        zzmin=box_limits(3,1) ; zzmax=box_limits(3,2)

        write(*,*)'limits:',xxmin,xxmax,yymin,yymax,zzmin,zzmax
        
        if(TRIM(amr%ordering).eq.'hilbert')then

            dxmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
            do ilevel=1,amr%lmax
               dx=0.5d0**ilevel
               if(dx.lt.dxmax)exit
            end do
            lmin=ilevel
            bit_length=lmin-1
            maxdom=2**bit_length
            imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
            if(bit_length>0)then
               imin=int(xxmin*dble(maxdom))
               imax=imin+1
               jmin=int(yymin*dble(maxdom))
               jmax=jmin+1
               kmin=int(zzmin*dble(maxdom))
               kmax=kmin+1
            endif
       
            dkey=(dble(2**(amr%nlevelmax+1)/dble(maxdom)))**amr%ndim
            amr%ndom=1
            if(bit_length>0)amr%ndom=8
            idom(1)=imin; idom(2)=imax
            idom(3)=imin; idom(4)=imax
            idom(5)=imin; idom(6)=imax
            idom(7)=imin; idom(8)=imax
            jdom(1)=jmin; jdom(2)=jmin
            jdom(3)=jmax; jdom(4)=jmax
            jdom(5)=jmin; jdom(6)=jmin
            jdom(7)=jmax; jdom(8)=jmax
            kdom(1)=kmin; kdom(2)=kmin
            kdom(3)=kmin; kdom(4)=kmin
            kdom(5)=kmax; kdom(6)=kmax
            kdom(7)=kmax; kdom(8)=kmax
            
            do i=1,amr%ndom
               if(bit_length>0)then
                  call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
               else
                  order_min=0.0d0
               endif
               bounding_min(i)=(order_min(1))*dkey
               bounding_max(i)=(order_min(1)+1.0D0)*dkey
            end do
            
            cpu_min=0; cpu_max=0
            do impi=1,amr%ncpu
               do i=1,amr%ndom
                  if (   amr%bound_key(impi-1).le.bounding_min(i).and.&
                       & amr%bound_key(impi  ).gt.bounding_min(i))then
                     cpu_min(i)=impi
                  endif
                  if (   amr%bound_key(impi-1).lt.bounding_max(i).and.&
                       & amr%bound_key(impi  ).ge.bounding_max(i))then
                     cpu_max(i)=impi
                  endif
               end do
            end do
            
            amr%ncpu_read=0
            do i=1,amr%ndom
               do j=cpu_min(i),cpu_max(i)
                  if(.not. amr%cpu_read(j))then
                     amr%ncpu_read=amr%ncpu_read+1
                     amr%cpu_list(amr%ncpu_read)=j
                     amr%cpu_read(j)=.true.
                  endif
               enddo
            enddo
         else
            amr%ncpu_read=amr%ncpu
            do j=1,amr%ncpu
               amr%cpu_list(j)=j
            end do
         end  if
    end subroutine get_cpu_map
end module io_ramses

module filtering
    use local
    use io_ramses

    type filter
        character(128) :: name
        integer :: ncond
        character(128), dimension(:), allocatable :: cond_vars
        character(2), dimension(:), allocatable :: cond_ops
        real(dbl), dimension(:), allocatable :: cond_vals
    end type filter

    contains

    subroutine allocate_filter(filt)
        implicit none
        type(filter), intent(inout) :: filt

        if (.not.allocated(filt%cond_vars)) allocate(filt%cond_vars(filt%ncond))
        if (.not.allocated(filt%cond_ops)) allocate(filt%cond_ops(filt%ncond))
        if (.not.allocated(filt%cond_vals)) allocate(filt%cond_vals(filt%ncond))
    end subroutine allocate_filter

    subroutine cond_string_to_filter(str, filt)
        implicit none
        character(128), intent(in) :: str
        type(filter), intent(inout) :: filt

        ! TODO: Finish this subroutine
    end subroutine cond_string_to_filter

    logical function filter_cell(varIDs,reg,filt,cell_x,cell_dx,cell_var)
        use vectors
        use geometrical_regions
        type(hydroID), intent(in) :: varIDs
        type(region), intent(in) :: reg
        type(filter), intent(in) :: filt
        real(dbl), intent(in) :: cell_dx
        type(vector), intent(in) :: cell_x
        real(dbl), dimension(1:varIDs%nvar), intent(in) :: cell_var
        integer :: i
        real(dbl) :: value

        filter_cell = .true.

        if (filt%ncond == 0) return

        do i=1,filt%ncond
            call getvarvalue(varIDs,reg,cell_dx,cell_x,cell_var,filt%cond_vars(i),value)
            select case (TRIM(filt%cond_ops(i)))
            case('/=')
                filter_cell = filter_cell .and. (value /= filt%cond_vals(i))
            case('==')
                filter_cell = filter_cell .and. (value == filt%cond_vals(i))
            case('<')
                filter_cell = filter_cell .and. (value < filt%cond_vals(i))
            case('<=')
                filter_cell = filter_cell .and. (value <= filt%cond_vals(i))
            case('>')
                filter_cell = filter_cell .and. (value > filt%cond_vals(i))
            case('>=')
                filter_cell = filter_cell .and. (value >= filt%cond_vals(i))
            case default
                write(*,*)'Relation operator not supported: ',TRIM(filt%cond_ops(i))
                write(*,*)'Aborting!'
                stop
            end select
        end do
    end function filter_cell

end module filtering