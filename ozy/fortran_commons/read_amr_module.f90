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
    use constants
    use vectors
    use cooling_module

    type hydroID
        integer :: nvar
        integer :: density=0,vx=0,vy=0,vz=0,thermal_pressure=0,metallicity=0
        integer :: Blx=0,Bly=0,Blz=0,Brx=0,Bry=0,Brz=0
        integer :: cr_pressure=0
        integer :: xHII=0,xHeII=0,xHeIII=0
        integer :: dust_density=0
        integer :: sigma2=0
    end type hydroID

    type amr_info
        integer :: ncpu,ndim,nlevelmax,nboundary,twotondim,ndom
        integer :: twondim,ngridmax,ncoarse
        integer :: levelmin,levelmax,lmax,lmin, active_lmax
        integer :: ncpu_read
        character(80) :: ordering
        integer,dimension(:),allocatable :: cpu_list
        real(dbl),dimension(:),allocatable :: bound_key
        logical,dimension(:),allocatable :: cpu_read
        real(dbl),dimension(1:3) :: xbound=(/0d0,0d0,0d0/)
    end type amr_info

    type sim_info
        logical :: cosmo=.true.,family=.false.
        logical :: dm=.false.,hydro=.false.,mhd=.false.
        logical :: cr=.false.,rt=.false.,bh=.false.
        logical :: cr_st=.false.,cr_heat=.false.,dust=.false.
        real(dbl) :: h0,t,aexp,unit_l,unit_d,unit_t,unit_m,unit_v,unit_p
        real(dbl) :: boxlen,omega_m,omega_l,omega_k,omega_b
        real(dbl) :: time_tot,time_simu,redshift,T2,nH
        integer :: n_frw
        real(dbl),dimension(:),allocatable :: aexp_frw,hexp_frw,tau_frw,t_frw
        real(dbl) :: eta_sn=-1D0
        real(dbl) :: Dcr=3D28
    end type sim_info

    type level
        logical::active
        integer::ilevel
        integer::ngrid
        integer,dimension(:),allocatable :: ind_grid
        integer,dimension(:),allocatable :: real_ind
        real(dbl),dimension(:,:),allocatable :: xg
        real(dbl),dimension(:,:,:,:),pointer::cube
        real(dbl),dimension(:,:,:),pointer::map
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

    type particle
#ifndef LONGINT
        integer(irg) :: id
#else
        integer(ilg) :: id
#endif
        type(vector) :: x,v
        real(dbl) :: m,met,imass,age,tform
    end type particle

    ! Define global variables
    type(sim_info) :: sim
    type(amr_info) :: amr
    type(hydroID)  :: varIDs
    logical :: verbose=.false.
    logical :: fix_neg_temp=.true.

    contains

    subroutine deactivate_verbose
        implicit none
        verbose = .false.
    end subroutine deactivate_verbose

    subroutine activate_verbose
        implicit none
        verbose = .true.
    end subroutine activate_verbose

    subroutine deactivate_fix_neg_temp
        implicit none
        fix_neg_temp = .false.
    end subroutine deactivate_fix_neg_temp

    subroutine activate_fix_neg_temp
        implicit none
        fix_neg_temp = .true.
    end subroutine activate_fix_neg_temp

    subroutine retrieve_vars(repository,myvars)
        implicit none

        character(128),intent(in) ::  repository
        type(hydroID), intent(inout) :: myvars

        call read_hydrofile_descriptor(repository)
        myvars = varIDs
    end subroutine retrieve_vars

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
    subroutine check_lmax(ngridfile)

        implicit none

        integer,dimension(1:amr%ncpu+amr%nboundary,1:amr%nlevelmax),intent(in) :: ngridfile
        integer :: ngridilevel,i

        do i = amr%nlevelmax, 0, -1
            ngridilevel=sum(ngridfile(:,i))
            if (ngridilevel.gt.0) then
                if (amr%lmax.gt.i) then
                    amr%active_lmax=i
                endif
                exit
            endif
        end do
    end subroutine check_lmax

    !---------------------------------------------------------------
    ! Subroutine: CHECK FAMILY
    !
    ! This routine reads the header_xxxxx.txt file from a snapshot
    ! in order to determine if particle files follow or not the
    ! family/tag format
    !---------------------------------------------------------------
    
    subroutine check_families(repository)
        implicit none

        character(len=128),intent(in) ::  repository

        character(len=128) :: nomfich,line,fline
        character(len=10) :: var1,var2,var3,var4,var5,var6
        character(5) :: nchar
        integer :: ipos,status

        ipos = index(repository,'output_')
        nchar = repository(ipos+7:ipos+13)
        nomfich = TRIM(repository)//'/header_'//TRIM(nchar)//'.txt'
        open(unit=10,file=nomfich,form='formatted',status='old')
        do
            read(10,'(A)',iostat=status)line
            if (status /= 0) exit
            fline = line
        end do
        read(fline,*)var1,var2,var3,var4,var5,var6
        if (trim(var6) .eq. 'family') then
            if (verbose) write(*,*)': This simulation uses particle families'
            sim%family = .true.
        else if (trim(var6) .eq. 'tform') then
            if (verbose) write(*,*)': This simulation uses the old particle format'
        else
            if (verbose) write(*,*)': This simulation format for particles is not recognised!'
            stop
        endif
    end subroutine check_families

    !---------------------------------------------------------------
    ! Subroutine: READ HYDRO IDs
    !
    ! This routine extracts the HYDRO variable IDs for a given
    ! simulation snapshot such that their order is known when
    ! hydro_*.out* files are read.
    !---------------------------------------------------------------    
    subroutine read_hydrofile_descriptor(repository)
        implicit none

        character(128),intent(in) ::  repository
        character(128) :: nomfich
        logical            ::  ok
        character(8)   ::  igr8
        character   ::  igr1,igrstart
        character(2)::igr2
        character(3)::igr3
        integer            ::  newID,statn

        nomfich=TRIM(repository)//'/hydro_file_descriptor.txt'
        inquire(file=nomfich, exist=ok) ! verify input file
        if ( ok ) then
            open(unit=10,file=nomfich,status='old',form='formatted')
            read(10,*) igrstart,igr8,igr1,igr2
            close(10)
            if (igrstart .eq. "n") then
                if (verbose) write(*,*) "This is an old RAMSES simulation"
                call read_hydrofile_descriptor_old(repository)
            elseif (igrstart .eq. "#") then
                if (verbose) write(*,*) "This is a new RAMSES simulation"
                call read_hydrofile_descriptor_new(repository)
            else
                if (verbose) write(*,*)" I do not recognise this sim format. Check!"
                stop
            endif
        else
            if (verbose) write(*,'(": ",A," not found. Initializing variables to default IDs.")') trim(nomfich)
            varIDs%density = 1
            varIDs%vx = 2; varIDs%vy  = 3; varIDs%vz  = 4
            varIDs%Blx  = 5; varIDs%Bly = 6; varIDs%Blz = 7
            varIDs%Brx  = 8; varIDs%Bry = 9; varIDs%Brz = 10
            varIDs%thermal_pressure = 11; varIDs%metallicity = 12;
            varIDs%xHII = 13; varIDs%xHeII = 14; varIDs%xHeIII = 15
        end if
    end subroutine read_hydrofile_descriptor

    !---------------------------------------------------------------
    ! Subroutine: READ HYDRO IDs OLD
    !
    ! This routine extracts the HYDRO variable IDs for a given
    ! simulation snapshot such that their order is known when
    ! hydro_*.out* files are read. This is for the old RAMSES
    ! format
    !---------------------------------------------------------------    
    subroutine read_hydrofile_descriptor_old(repository)
        implicit none

        character(128),intent(in) ::  repository
        character(128) :: nomfich
        logical            ::  ok
        integer            ::  nvhydro,nvloop,i !nvar,i
        character(25)  ::  newVar
        character(9)   ::  igr9
        character(8)   ::  igr8
        character   ::  igr1
        character(2)::igr2
        character(3)::igr3
        integer            ::  newID,statn
        nomfich=TRIM(repository)//'/hydro_file_descriptor.txt'
        inquire(file=nomfich, exist=ok) ! verify input file
        if ( ok ) then
            if (verbose) write(*,'(": Reading variables IDs from hydro_descriptor")')
            open(unit=10,file=nomfich,status='old',form='formatted')
            ! read(10,'("nvar        =",I11)')nvar
            read(10,*) igr9,igr1,igr2
            read(igr2,*,iostat=statn) nvhydro
            if (verbose) write(*,*)'nvar=',nvhydro
            varIDs%nvar = nvhydro
            if (nvhydro > 9) then
                nvloop = 9
            else
                nvloop = nvhydro
            end if
            do i=1,nvloop
                read(10,*) igr9,igr1,igr1,newVar
                read(igr1,*,iostat=statn) newID
                call select_from_descriptor_IDs(newVar,newID)
            end do
            if (nvhydro > 10) then
                do i=nvloop+1,nvhydro
                    read(10,*) igr8,igr3,newVar
                    igr3 = igr3(2:3);
                    read(igr3,*,iostat=statn) newID
                    call select_from_descriptor_IDs(newVar,newID)
                end do
            end if
            
            ! varIDs%nvar = nvar
            ! do i=1,nvar
            !     read(10,*) igr9,igr1,igr1,newVar
            !     read(igr1,*,iostat=statn) newID
            !     call select_from_descriptor_IDs(newVar,newID)
            ! end do
            close(10)
        else
            if (verbose) write(*,'(": ",A," not found. Initializing variables to default IDs.")') trim(nomfich)
            varIDs%density = 1
            varIDs%vx = 2; varIDs%vy  = 3; varIDs%vz  = 4
            varIDs%Blx  = 5; varIDs%Bly = 6; varIDs%Blz = 7
            varIDs%Brx  = 8; varIDs%Bry = 9; varIDs%Brz = 10
            varIDs%thermal_pressure = 11; varIDs%metallicity = 12;
            varIDs%xHII = 13; varIDs%xHeII = 14; varIDs%xHeIII = 15
        end if
    end subroutine read_hydrofile_descriptor_old


    subroutine select_from_descriptor_IDs(newVar,newID)
        implicit none
        integer,intent(in)           :: newID
        character(25),intent(in) :: newvar
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
            sim%mhd = .true.
        case ('B_left_y')
            varIDs%Bly = newID
            sim%mhd = .true.
        case ('B_left_z')
            varIDs%Blz = newID
            sim%mhd = .true.
        case ('B_right_x')
            varIDs%Brx = newID
            sim%mhd = .true.
        case ('B_right_y')
            varIDs%Bry = newID
            sim%mhd = .true.
        case ('B_right_z')
            varIDs%Brz = newID
            sim%mhd = .true.
        case ('B_x_left')
            varIDs%Blx = newID
            sim%mhd = .true.
        case ('B_y_left')
            varIDs%Bly = newID
            sim%mhd = .true.
        case ('B_z_left')
            varIDs%Blz = newID
            sim%mhd = .true.
        case ('B_x_right')
            varIDs%Brx = newID
            sim%mhd = .true.
        case ('B_y_right')
            varIDs%Bry = newID
            sim%mhd = .true.
        case ('B_z_right')
            varIDs%Brz = newID
            sim%mhd = .true.
        case ('thermal_pressure')
            varIDs%thermal_pressure = newID
        case ('pressure')
            varIDs%thermal_pressure = newID
        case ('non_thermal_pressure_1')
            if (verbose) write(*,'(": Using non_thermal_pressure_1 as cosmic ray pressure (variable ",I2,")")') newID
            varIDs%cr_pressure = newID
            sim%cr = .true.
        case ('cosmic_ray_01')
            if (verbose) write(*,'(": Using cosmic_ray_01 as cosmic ray pressure (variable ",I2,")")') newID
            varIDs%cr_pressure = newID
            sim%cr = .true.
        case ('passive_scalar_1')
            if (verbose) write(*,'(": Using passive_scalar_1 as metallicity (variable ",I2,")")') newID
            varIDs%metallicity = newID
        case ('metallicity')
            varIDs%metallicity = newID
        case ('passive_scalar_2')
            if (verbose) write(*,'(": Using passive_scalar_2 as xHII (variable ",I2,")")') newID
            varIDs%xHII = newID
        case ('scalar_01')
            if (verbose) write(*,'(": Using scalar_01 as xHII (variable ",I2,")")') newID
            varIDs%xHII = newID
        case ('passive_scalar_3')
            if (verbose) write(*,'(": Using passive_scalar_3 as xHeII (variable ",I2,")")') newID
            varIDs%xHeII = newID
        case ('scalar_02')
            if (verbose) write(*,'(": Using scalar_02 as xHeII (variable ",I2,")")') newID
            varIDs%xHeII = newID
        case ('passive_scalar_4')
            if (verbose) write(*,'(": Using passive_scalar_4 as xHeIII (variable ",I2,")")') newID
            varIDs%xHeIII = newID
        case ('scalar_03')
            if (verbose) write(*,'(": Using scalar_03 as xHeIII (variable ",I2,")")') newID
            varIDs%xHeIII = newID
        case ('xHII')
            varIDs%xHII = newID
        case ('xHeII')
            varIDs%xHeII = newID
        case ('xHeIII')
            varIDs%xHeIII = newID
        case ('H_p1_fraction')
            varIDs%xHII = newID
        case ('He_p1_fraction')
            varIDs%xHeII = newID
        case ('He_p2_fraction')
            varIDs%xHeIII = newID
        case ('dust')
            sim%dust = .true.
            varIDs%dust_density = newID
        case ('sigma2')
            varIDs%sigma2 = newID
        end select
    end subroutine select_from_descriptor_IDs

    !---------------------------------------------------------------
    ! Subroutine: READ HYDRO IDs NEW
    !
    ! This routine extracts the HYDRO variable IDs for a given
    ! simulation snapshot such that their order is known when
    ! hydro_*.out* files are read. This is for the new RAMSES
    ! format
    !---------------------------------------------------------------    
    subroutine read_hydrofile_descriptor_new(repository)
        implicit none

        character(128),intent(in) ::  repository
        character(128) :: nomfich
        logical            ::  ok
        character(25)  ::  newVar,newType
        integer            ::  newID,status,nvar=0

        nomfich=TRIM(repository)//'/hydro_file_descriptor.txt'
        inquire(file=nomfich, exist=ok) ! verify input file
        if ( ok ) then
            if (verbose) write(*,'(": Reading variables IDs from hydro_descriptor")')
            open(unit=10,file=nomfich,status='old',form='formatted')
            read(10,*)
            read(10,*)
            do
                read(10,*,iostat=status)newID,newVar,newType
                if (status /= 0) exit
                nvar = nvar + 1
                call select_from_descriptor_IDs(newVar,newID)
            end do
            close(10)
            if (verbose) write(*,*)'nvar=',nvar
            varIDs%nvar = nvar
        else
            if (verbose) write(*,'(": ",A," not found. Initializing variables to default IDs.")') trim(nomfich)
            varIDs%density = 1
            varIDs%vx = 2; varIDs%vy  = 3; varIDs%vz  = 4
            varIDs%Blx  = 5; varIDs%Bly = 6; varIDs%Blz = 7
            varIDs%Brx  = 8; varIDs%Bry = 9; varIDs%Brz = 10
            varIDs%thermal_pressure = 11; varIDs%metallicity = 12;
            varIDs%xHII = 13; varIDs%xHeII = 14; varIDs%xHeIII = 15
        end if
    end subroutine read_hydrofile_descriptor_new

    !---------------------------------------------------------------
    ! Subroutine: GET VARIABLE VALUE
    !
    ! For a given cell, a particular variable is computed.
    ! This is a very important subroutine, which will be modified
    ! extensively.
    !---------------------------------------------------------------
    subroutine getvarvalue(reg,dx,x,var,son,varname,value,trans_matrix,grav_var)
        use vectors
        use basis_representations
        use geometrical_regions
        use coordinate_systems
        implicit none
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:amr%twondim,1:varIDs%nvar),intent(in) :: var
        integer,dimension(0:amr%twondim),intent(in) :: son
        character(128),intent(in)                 :: varname
        real(dbl),intent(inout)                       :: value
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:amr%twondim,1:4),optional,intent(in) :: grav_var
        real(dbl),dimension(0:amr%twondim,1:varIDs%nvar) :: tempvar
        real(dbl),dimension(0:amr%twondim) :: totP
        integer :: i
        type(vector) :: v,L,B,vst
        type(basis) :: temp_basis
        real(dbl) :: T,rho,cV,lambda,lambda_prime,ne,ecr,nH,Tmin,Dcr,vA,cs,sigma
        real(dbl) :: dxleft,dxright
        real(dbl) :: P_r,fg_r,r
        real(dbl) :: F_a,F_s,F_d,Fcr
        real(dbl) :: bsign,gamma
        real(dbl) :: volume
        real(dbl) :: tcool,tcomp
        real(dbl) :: lambda_co, lambda_st, lambda_cr
        character(128) :: star_maker

        Tmin = 15d0 / sim%T2
        Dcr = sim%Dcr / (sim%unit_l**2 / sim%unit_t)

        select case (TRIM(varname))
        case ('d_euclid')
            ! Euclidean distance
            value = magnitude(x)
        case ('x')
            ! x - coordinate
            value = x%x
        case ('y')
            ! y - coordinate
            value = x%y
        case ('z')
            ! z - coordinate
            value = x%z
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
        case ('vx_proj')
            ! X-velocity in the projected plane
            ! v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            ! call rotate_vector(v,trans_matrix)
            value = var(0,varIDs%vx)
        case ('vy_proj')
            ! Y-velocity in the projected plane
            ! v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            ! call rotate_vector(v,trans_matrix)
            value = var(0,varIDs%vy)
        case ('vz_proj')
            ! Z-velocity in the projected plane
            ! v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            ! call rotate_vector(v,trans_matrix)
            value = var(0,varIDs%vz)
        case ('v_tangential')
            ! Tangential velocity
            ! Consider as if one substracts the radial velocity
            ! to the current Velocity vector
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call spherical_basis_from_cartesian(x,temp_basis)
            v = v - (v.DOT.temp_basis%u(1)) * temp_basis%u(1)
            value = magnitude(v)
        case ('centripetal_acc')
            ! Radial centripetal acceleration, obtained from the
            ! tangential velocity
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call spherical_basis_from_cartesian(x,temp_basis)
            v = v - (v.DOT.temp_basis%u(1)) * temp_basis%u(1)
            value = magnitude(v)**2d0 / r_sphere(x)
        case ('v_sphere_r')
            ! Velocity component in the spherical radial direction
            ! Dot product of velocity vector with spherical radial
            !    unit vector
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = v.DOT.temp_basis%u(1)
        case ('absv_sphere_r')
            ! Absolute elocity component in the spherical radial direction
            ! Dot product of velocity vector with spherical radial
            !    unit vector
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = abs(v.DOT.temp_basis%u(1))
        case ('v_sphere_phi')
            ! Velocity component in the spherical azimutal (phi) direction
            ! Dot product of velocity vector with spherical phi
            !    unit vector
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = v .DOT. temp_basis%u(3)
        case ('v_sphere_theta')
            ! Velocity component in the spherical theta direction
            ! Dot product of velocity vector with spherical theta
            !    unit vector
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = v .DOT. temp_basis%u(2)
        case ('v_cyl_r')
            ! Velocity component in the cylindrical radial direction
            ! Dot product of velocity vector with cylindrical
            !    radial unit vector
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = v.DOT.temp_basis%u(1)
        case ('v_cyl_z')
            ! Velocity component in the cylindrical z direction
            ! Dot product of velocity vector with cylindrical
            !    z unit vector
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = v.DOT.temp_basis%u(3)
        case ('v_cyl_phi')
            ! Velocity component in the cylyndrical azimutal (phi) direction
            ! Dot product of velocity vector with cylindrical
            !    phi unit vector
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = v.DOT.temp_basis%u(2)
        case ('vrcyl_overv')
            ! Velocity component in the cylindrical radial direction over total velocity magnitude
            ! Dot product of velocity vector with cylindrical
            !    radial unit vector
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = abs(v.DOT.temp_basis%u(1)) / magnitude(v)
        case ('vzcyl_overv')
            ! Velocity component in the cylindrical z direction over total velocity magnitude
            ! Dot product of velocity vector with cylindrical
            !    z unit vector
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = abs(v.DOT.temp_basis%u(3)) / magnitude(v)
        case ('vphicyl_overv')
            ! Velocity component in the cylyndrical azimutal (phi) direction over total velocity magnitude
            ! Dot product of velocity vector with cylindrical
            !    phi unit vector
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = abs(v.DOT.temp_basis%u(2)) / magnitude(v)
        case ('v_magnitude')
            ! Velocity magnitude from galaxy coordinates
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            value = magnitude(v)
        case ('v_squared')
            ! Velocity magnitude squared from galaxy coordinates
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            value = v.DOT.v
        case ('div_v')
            ! Velocity divergence
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%vx) - var(1,varIDs%vx)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%vy) - var(3,varIDs%vy)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%vz) - var(5,varIDs%vz)) / (dxright + dxleft)
            value = v%x + v%y + v%z
        case ('density')
            ! Density
            value = var(0,varIDs%density)
        case ('dust_mass')
            ! Dust mass
            value = (var(0,varIDs%dust_density) * var(0,varIDs%density) * (dx*dx)) * dx
        case ('DTM')
            ! Dust-to-metal ratio
            value = var(0,varIDs%dust_density) / (var(0,varIDs%metallicity)*var(0,varIDs%density))
        case ('mass')
            ! Total mass, computed as density times volume of cell (dx**3)
            value = (var(0,varIDs%density) * (dx*dx)) * dx
        case ('volume')
            ! Volume of cell (dx**3)
            value = (dx * dx) * dx
        case ('metallicity')
            ! Metallicity
            value = var(0,varIDs%metallicity)/0.02
        case ('metal_mass')
            ! Total metal mass
            value = (var(0,varIDs%metallicity) * var(0,varIDs%density) * (dx*dx)) * dx
        case ('dust_density')
            ! Dust density
            ! TODO: For dust simulation it should be updated
            if (.not. sim%dust) then
                value = var(0,varIDs%metallicity) * (0.4d0 *var(0,varIDs%density))
            else
                value = var(0,varIDs%dust_density) * var(0,varIDs%density)
            end if
        case ('temperature')
            ! Gas temperature
            value = var(0,varIDs%thermal_pressure) / var(0,varIDs%density)
            ! TODO: This is a quick fix for negative T in CRMHD sims
            if (value < Tmin .and. fix_neg_temp) value = Tmin
        case ('thermal_pressure')
            ! Thermal pressure
            value = max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density))
        case ('torque_therp_specific')
            ! Magnitude of the specific torque due to thermal pressure gradient
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (max(var(2,varIDs%thermal_pressure), Tmin*var(2,varIDs%density)) - &
                    & max(var(1,varIDs%thermal_pressure),Tmin*var(1,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (max(var(4,varIDs%thermal_pressure), Tmin*var(4,varIDs%density)) - &
                    & max(var(3,varIDs%thermal_pressure),Tmin*var(3,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (max(var(6,varIDs%thermal_pressure), Tmin*var(6,varIDs%density)) - &
                    & max(var(5,varIDs%thermal_pressure),Tmin*var(5,varIDs%density))) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = (r_sphere(x) / var(0,varIDs%density)) * magnitude(temp_basis%u(1) * v)
        case ('grad_thermalpressure')
            ! Magnitude of thermal pressure gradient
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (max(var(2,varIDs%thermal_pressure), Tmin*var(2,varIDs%density)) - &
                    & max(var(1,varIDs%thermal_pressure),Tmin*var(1,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (max(var(4,varIDs%thermal_pressure), Tmin*var(4,varIDs%density)) - &
                    & max(var(3,varIDs%thermal_pressure),Tmin*var(3,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (max(var(6,varIDs%thermal_pressure), Tmin*var(6,varIDs%density)) - &
                    & max(var(5,varIDs%thermal_pressure),Tmin*var(5,varIDs%density))) / (dxright + dxleft)
            value = magnitude(v)
        case ('grad_therprsphere')
            ! Thermal pressure gradient in the radial direction
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (max(var(2,varIDs%thermal_pressure), Tmin*var(2,varIDs%density)) - &
                    & max(var(1,varIDs%thermal_pressure),Tmin*var(1,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (max(var(4,varIDs%thermal_pressure), Tmin*var(4,varIDs%density)) - &
                    & max(var(3,varIDs%thermal_pressure),Tmin*var(3,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (max(var(6,varIDs%thermal_pressure), Tmin*var(6,varIDs%density)) - &
                    & max(var(5,varIDs%thermal_pressure),Tmin*var(5,varIDs%density))) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = v.DOT.temp_basis%u(1)
        case ('grad_therpz')
            ! Thermal pressure gradient in the z direction
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (max(var(2,varIDs%thermal_pressure), Tmin*var(2,varIDs%density)) - &
                    & max(var(1,varIDs%thermal_pressure),Tmin*var(1,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (max(var(4,varIDs%thermal_pressure), Tmin*var(4,varIDs%density)) - &
                    & max(var(3,varIDs%thermal_pressure),Tmin*var(3,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (max(var(6,varIDs%thermal_pressure), Tmin*var(6,varIDs%density)) - &
                    & max(var(5,varIDs%thermal_pressure),Tmin*var(5,varIDs%density))) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            value = v%z
        case ('thermal_energy')
            ! Thermal energy, computed as thermal_pressure*volume/(gamma - 1)
            value = ((max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density)) / (gamma_gas - 1d0)) * (dx * dx)) * dx
        case ('thermal_energy_specific')
            ! Specific thermal energy as E_ther/cell mass
            value = (max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density)) / (gamma_gas - 1d0)) / var(0,varIDs%density)
        case ('thermal_energy_density')
            ! Thermal energy density  as E_ther/cell volume
            value = max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density)) / (gamma_gas - 1d0)
        case ('entropy_specific')
            ! Specific entropy, following Gent 2012 equation
            ! cV = cVHydrogen * mHydrogen / kBoltzmann / mu
            ! T = (var(0,varIDs%thermal_pressure)*((sim%unit_l/sim%unit_t)**2) / var(0,varIDs%density) / kBoltzmann * mHydrogen)
            ! !TODO: This is a fix to the low temperature in CRMHD
            ! if (T<15d0 .and. fix_neg_temp) T = 15d0
            ! rho = (var(0,varIDs%density) * sim%unit_d / mHydrogen )
            ! value = cV * (log(T) - (gamma_gas - 1d0) * log(rho))
            T = max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density))/var(0,varIDs%density) * sim%T2 * mu
            rho = var(0,varIDs%density) * sim%unit_d / mHydrogen
            value = log(T) - (gamma_gas-1d0) * log(rho)
        case ('pseudo_entropy')
            ! Pseudo-entropy K, as usually defined for clusters
            if (fix_neg_temp) then
                T = max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density))
            else
                T = var(0,varIDs%thermal_pressure)
            end if
            value = T / (var(0,varIDs%density)**gamma_gas)
        case ('grad_entropy_cylz')
            ! Gradient of the pseudo entropy in the z direction
            totP(:) = 0d0
            do i=1,amr%twondim
                if (fix_neg_temp) then
                    T = max(var(i,varIDs%thermal_pressure), Tmin*var(i,varIDs%density))
                else
                    T = var(i,varIDs%thermal_pressure)
                end if
                totP(i) = T / (var(i,varIDs%density)**gamma_gas)
            end do
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (totP(2) - totP(1)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (totP(4) - totP(3)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (totP(6) - totP(5)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = v.DOT.temp_basis%u(3)
        case ('grad_entropy_rsphere')
            ! Gradient of the pseudo entropy in the r direction
            totP(:) = 0d0
            do i=1,amr%twondim
                if (fix_neg_temp) then
                    T = max(var(i,varIDs%thermal_pressure), Tmin*var(i,varIDs%density))
                else
                    T = var(i,varIDs%thermal_pressure)
                end if
                totP(i) = T / (var(i,varIDs%density)**gamma_gas)
            end do
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (totP(2) - totP(1)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (totP(4) - totP(3)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (totP(6) - totP(5)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = v.DOT.temp_basis%u(1)
        case ('sound_speed')
            ! Thermal sound speed, ideal gas
            value = sqrt(gamma_gas * (max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density)) / var(0,varIDs%density)))
        case ('effective_sound_speed')
            ! Effective sound speed of the thermal and CR fluid mixture
            value = sqrt((gamma_gas * var(0,varIDs%thermal_pressure) + gamma_crs * var(0,varIDs%cr_pressure)) / &
                        & var(0,varIDs%density))
        case ('rms_speed')
            ! RMS speed, ideal gas (Maxwellian distribution)
            value = sqrt(3d0 * (max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density)) / var(0,varIDs%density)))
        case ('kinetic_energy')
            ! Kinetic energy, computed as 1/2*density*volume*magnitude(velocity)
            value = (0.5 * (var(0,varIDs%density) * (dx*dx)) * dx) * sqrt(var(0,varIDs%vx)**2 + var(0,varIDs%vy)**2 + var(0,varIDs%vz)**2)
        case ('magnetic_energy')
            ! Magnetic energy as magnitude(B)**2/2
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            value = (0.5 * (B.DOT.B) * (dx*dx)) * dx
        case ('Bx')
            value = 0.5 *(var(0,varIDs%Blx)+var(0,varIDs%Brx))
        case ('Bx_proj')
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            v = B
            call rotate_vector(v,trans_matrix)
            value = v%x
        case ('By_proj')
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            v = B
            call rotate_vector(v,trans_matrix)
            value = v%y
        case ('Bz_proj')
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            v = B
            call rotate_vector(v,trans_matrix)
            value = v%z
        case ('By')
            value = 0.5 *(var(0,varIDs%Bly)+var(0,varIDs%Bry))
        case ('Bz')
            value = 0.5 *(var(0,varIDs%Blz)+var(0,varIDs%Brz))
        case ('magnetic_magnitude')
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            value = magnitude(B)
        case ('Bz_cyl')
            ! Magnetic field in the cylindrical z direction
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = B.DOT.temp_basis%u(3)
        case ('Bphi_cyl')
            ! Magnetic field in the cylindrical phi direction
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = B.DOT.temp_basis%u(2)
        case ('Br_cyl')
            ! Magnetic field in the cylindrical radial direction
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = B.DOT.temp_basis%u(1)
        case ('Br_sphere')
            ! Magnetic field in the spherical r direction
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = B.DOT.temp_basis%u(1)
        case ('Btheta_sphere')
            ! Magnetic field in the spherical theta direction
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = B.DOT.temp_basis%u(2)
        case ('Bphi_sphere')
            ! Magnetic field in the spherical phi direction
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = B.DOT.temp_basis%u(3)
        case ('Bzcyl_overB')
            ! Ratio of magnetic field in the cylindrical z direction to total B-field
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = abs(B.DOT.temp_basis%u(3)) / magnitude(B)
        case ('Bphicyl_overB')
            ! Ratio of magnetic field in the cylindrical phi direction to total B-field
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = abs(B.DOT.temp_basis%u(2)) / magnitude(B)
        case ('Brcyl_overB')
            ! Ratio of magnetic field in the cylindrical radial direction to total B-field
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = abs(B.DOT.temp_basis%u(1)) / magnitude(B)
        case ('magnetic_energy_specific')
            ! Specific magnetic energy as E_mag/cell mass
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            value = (0.5 * (B.DOT.B)) / var(0,varIDs%density)
        case ('magnetic_energy_density')
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            value = 0.5 * (B.DOT.B)
        case ('magnetic_pressure')
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            value = 0.5 * (B.DOT.B)
        case ('plasma_beta')
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            value = max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density)) / (0.5 * (B.DOT.B))
        case ('alfven_speed')
            ! Alfven speed defined as B / sqrt(rho)
            B = 0.5*(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            value = magnitude(B) / sqrt(var(0,varIDs%density))
        case ('cr_GH08heat')
            ! Cosmic rays hadronic and Coulomb heating from Guo&Ho(2008)
            ! (Assume fully ionised gas)
            ! TODO: Update for RT! 
            lambda = 2.63d-16 * ((sim%unit_t**3)/(sim%unit_d*(sim%unit_l**2)))
            ne = var(0,varIDs%density) * sim%unit_d / mHydrogen 
            ecr = var(0,varIDs%cr_pressure) / (4D0/3d0 - 1d0)
            ecr = ecr * (sim%unit_d * ((sim%unit_l/sim%unit_t)**2))
            value = lambda * ne * ecr
        case ('net_cooling')
            ! Net cooling rate taken from the cooling table in RAMSES output
            T = var(0,varIDs%thermal_pressure) / var(0,varIDs%density) * sim%T2
            ! TODO: This a fix only for some messed up CRMHD simulations!
            if (T<15d0) T = 15d0
            nH = var(0,varIDs%density) * sim%nH
            call solve_net_cooling(nH,T,var(0,varIDs%metallicity)/2D-2,lambda,lambda_prime)
            value = ((lambda * nH) * nH) * ((sim%unit_t**3)/(sim%unit_d*(sim%unit_l**2)))
        case ('cooling_rate')
            ! Cooling rate taken from the cooling table in RAMSES output
            T = var(0,varIDs%thermal_pressure) / var(0,varIDs%density) * sim%T2
            ! TODO: This a fix only for some messed up CRMHD simulations!
            if (T<15d0) T = 15d0
            nH = var(0,varIDs%density) * sim%nH
            call solve_cooling(nH,T,var(0,varIDs%metallicity)/2D-2,lambda,lambda_prime)
            value = ((lambda * nH) * nH) * ((sim%unit_t**3)/(sim%unit_d*(sim%unit_l**2)))
        case ('heating_rate')
            ! Heating rate taken from the cooling table in RAMSES output
            T = var(0,varIDs%thermal_pressure) / var(0,varIDs%density) * sim%T2
            ! TODO: This a fix only for some messed up CRMHD simulations!
            if (T<15d0) T = 15d0
            nH = var(0,varIDs%density) * sim%nH
            call solve_heating(nH,T,var(0,varIDs%metallicity)/2D-2,lambda,lambda_prime)
            value = ((lambda * nH) * nH) * ((sim%unit_t**3)/(sim%unit_d*(sim%unit_l**2)))
        case ('B_left_x')
            value = var(0,varIDs%Blx)
        case ('B_left_y')
            value = var(0,varIDs%Bly)
        case ('B_left_z')
            value = var(0,varIDs%Blz)
        case ('B_right_x')
            value = var(0,varIDs%Brx)
        case ('B_right_y')
            value = var(0,varIDs%Bry)
        case ('B_right_z')
            value = var(0,varIDs%Brz)
        case ('cr_energy')
            ! CR energy, computed as CR_energydensity*volume
            value = ((var(0,varIDs%cr_pressure) / (gamma_crs - 1d0)) * (dx*dx)) * dx
        case ('cr_energy_density')
            value = var(0,varIDs%cr_pressure) / (gamma_crs - 1d0)
        case ('cr_pressure')
            value = var(0,varIDs%cr_pressure)
            ! if(value<1d-22)print*,value
        case ('norm_crflux_advection')
            ! Magnitude of the advection flux of CR energy density
            ! 1. We need to reset the local gas velocity to the frame of reference of the box
            v = var(0,varIDs%vx:varIDs%vz)
            v = v + reg%bulk_velocity
            ! 2. Compute the flux as simply e_CR * magnitude(u_gas)
            value = var(0,varIDs%cr_pressure) / (gamma_crs - 1d0) / magnitude(v)
        case ('norm_crflux_streaming')
            ! Magnitude of the advection-diffusion streaming flux of CR energy density
            ! 1. Compute magnetic field vector
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            ! 2. Compute the flux magnitude as (e_CR+P_CR) * v_a
            value = var(0,varIDs%cr_pressure) * (gamma_crs/(gamma_crs - 1d0)) * magnitude(B) / sqrt(var(0,varIDs%density))
        case ('norm_crflux_diffusion')
            ! Magnitude of the diffusion flux of CR energy density
            ! 1. Compute the CR energy density gradient
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            v = v / (gamma_crs - 1d0)
            ! 2. Compute the magnetic unit vector
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            B = B / magnitude(B)
            ! 3. Compute the flux magnitude as Dcr * (b .dot. grad(e_CR))
            value = Dcr * abs(B.DOT.v)
        case ('streamflux_diffflux_ratio')
            ! 1. Compute magnetic field vector
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            ! 2. Compute the flux magnitude as (e_CR+P_CR) * v_a
            F_s = var(0,varIDs%cr_pressure) * gamma_crs * magnitude(B) / sqrt(var(0,varIDs%density))
            ! 1. Compute the CR energy density gradient
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            ! 2. Compute the magnetic unit vector
            B = B / magnitude(B)
            F_d = Dcr * abs(B.dot.v)
            value = F_s / (F_d + F_s)
        case ('total_pressure')
            ! Total pressure (thermal, turbulent, magnetic and CR)
            ! TODO: Add radiation pressure
            totP(:) = 0d0
            if (sim%hydro) then
                totP(0) = var(0,varIDs%thermal_pressure)
                ! Go back to box coordinates for the central cell, which is transformed usually
                ! before sent to read_amr
                tempvar(:,:) = var(:,:)
                v = tempvar(0,varIDs%vx:varIDs%vz)
                call rotate_vector(v,transpose(trans_matrix))
                v = v + reg%bulk_velocity
                tempvar(0,varIDs%vx:varIDs%vz) = v
                call cmp_sigma_turb(tempvar,sigma)
                totP(0) = totP(0) + var(0,varIDs%density) * sigma**2d0
            end if
            if (sim%mhd) then
                B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),&
                    & (var(0,varIDs%Bly)+var(0,varIDs%Bry)),&
                    & (var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
                totP(0) = totP(0) + 0.5 * (B.DOT.B)
            end if
            if (sim%cr) totP(0) = totP(0) + var(0,varIDs%cr_pressure)
            value = totP(0)
        case ('fP_thermal')
            ! Ratio of thermal over total pressure (thermal, turbulent, magnetic and CR)
            ! TODO: Add radiation pressure
            totP(:) = 0d0
            if (sim%hydro) then
                totP(0) =  max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density))
                ! Go back to box coordinates for the central cell, which is transformed usually
                ! before sent to read_amr
                tempvar(:,:) = var(:,:)
                v = tempvar(0,varIDs%vx:varIDs%vz)
                call rotate_vector(v,transpose(trans_matrix))
                v = v + reg%bulk_velocity
                tempvar(0,varIDs%vx:varIDs%vz) = v
                call cmp_sigma_turb(tempvar,sigma)
                totP(0) = totP(0) + var(0,varIDs%density) * sigma**2d0
            end if
            if (sim%mhd) then
                B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),&
                    & (var(0,varIDs%Bly)+var(0,varIDs%Bry)),&
                    & (var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
                totP(0) = totP(0) + 0.5 * (B.DOT.B)
            end if
            if (sim%cr) totP(0) = totP(0) + var(0,varIDs%cr_pressure)
            value =  max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density)) / totP(0)
        case ('fP_turbulent')
            ! Ratio of turbulent over total pressure (thermal, turbulent, magnetic and CR)
            ! TODO: Add radiation pressure
            totP(:) = 0d0
            sigma = 0d0
            if (sim%hydro) then
                totP(0) =  max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density))
                ! Go back to box coordinates for the central cell, which is transformed usually
                ! before sent to read_amr
                tempvar(:,:) = var(:,:)
                v = tempvar(0,varIDs%vx:varIDs%vz)
                call rotate_vector(v,transpose(trans_matrix))
                v = v + reg%bulk_velocity
                tempvar(0,varIDs%vx:varIDs%vz) = v
                call cmp_sigma_turb(tempvar,sigma)
                totP(0) = totP(0) + var(0,varIDs%density) * sigma**2d0
            end if
            if (sim%mhd) then
                B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),&
                    & (var(0,varIDs%Bly)+var(0,varIDs%Bry)),&
                    & (var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
                totP(0) = totP(0) + 0.5 * (B.DOT.B)
            end if
            if (sim%cr) totP(0) = totP(0) + var(0,varIDs%cr_pressure)
            value = var(0,varIDs%density) * sigma**2d0 / totP(0)
        case ('fP_magnetic')
            ! Ratio of magnetic over total pressure (thermal, turbulent, magnetic and CR)
            ! TODO: Add radiation pressure
            totP(:) = 0d0
            sigma = 0d0
            if (sim%hydro) then
                totP(0) =  max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density))
                ! Go back to box coordinates for the central cell, which is transformed usually
                ! before sent to read_amr
                tempvar(:,:) = var(:,:)
                v = tempvar(0,varIDs%vx:varIDs%vz)
                call rotate_vector(v,transpose(trans_matrix))
                v = v + reg%bulk_velocity
                tempvar(0,varIDs%vx:varIDs%vz) = v
                call cmp_sigma_turb(tempvar,sigma)
                totP(0) = totP(0) + var(0,varIDs%density) * sigma**2d0
            end if
            if (sim%mhd) then
                B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),&
                    & (var(0,varIDs%Bly)+var(0,varIDs%Bry)),&
                    & (var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
                totP(0) = totP(0) + 0.5 * (B.DOT.B)
            end if
            if (sim%cr) totP(0) = totP(0) + var(0,varIDs%cr_pressure)
            value = 0.5 * (B.DOT.B) / totP(0)   
        case ('fP_cr')
            ! Ratio of CR over total pressure (thermal, turbulent, magnetic and CR)
            ! TODO: Add radiation pressure
            totP(:) = 0d0
            sigma = 0d0
            if (sim%hydro) then
                totP(0) =  max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density))
                ! Go back to box coordinates for the central cell, which is transformed usually
                ! before sent to read_amr
                tempvar(:,:) = var(:,:)
                v = tempvar(0,varIDs%vx:varIDs%vz)
                call rotate_vector(v,transpose(trans_matrix))
                v = v + reg%bulk_velocity
                tempvar(0,varIDs%vx:varIDs%vz) = v
                call cmp_sigma_turb(tempvar,sigma)
                totP(0) = totP(0) + var(0,varIDs%density) * sigma**2d0
            end if
            if (sim%mhd) then
                B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),&
                    & (var(0,varIDs%Bly)+var(0,varIDs%Bry)),&
                    & (var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
                totP(0) = totP(0) + 0.5 * (B.DOT.B)
            end if
            if (sim%cr) then
                totP(0) = totP(0) + var(0,varIDs%cr_pressure)
                value = var(0,varIDs%cr_pressure) / totP(0)
            else
                value = 0d0
            end if
        case ('total_energy_density')
            ! Total energy density, considering thermal, kinetic, magnetic and CR
            totP(:) = 0d0
            if (sim%hydro) then
                ! Thermal energy density
                totP(0) = var(0,varIDs%thermal_pressure) / (5D0/3d0 - 1d0)

                ! Kinetic energy density
                ! Go back to box coordinates for the central cell, which is transformed usually
                ! before sent to read_amr
                tempvar(:,:) = var(:,:)
                v = tempvar(0,varIDs%vx:varIDs%vz)
                call rotate_vector(v,transpose(trans_matrix))
                v = v + reg%bulk_velocity
                tempvar(0,varIDs%vx:varIDs%vz) = v
                totP(0) = totP(0) + 0.5d0 * var(0,varIDs%density) *sum(tempvar(0,varIDs%vx:varIDs%vz)**2d0)
            end if

            if (sim%mhd) then
                ! Magnetic energy density
                B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),&
                    & (var(0,varIDs%Bly)+var(0,varIDs%Bry)),&
                    & (var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
                totP(0) = totP(0) +  0.5 * (B.DOT.B)
            end if
            
            if (sim%cr) then
                ! CR energy density
                totP(0) = totP(0) + var(0,varIDs%cr_pressure) / (gamma_crs - 1d0)
            end if
            value = totP(0)
        case ('fe_thermal')
            ! Thermal energy density over total energy density, considering thermal, kinetic, magnetic and CR
            totP(:) = 0d0
            if (sim%hydro) then
                ! Thermal energy density
                totP(0) = var(0,varIDs%thermal_pressure) / (gamma_gas - 1d0)

                ! Kinetic energy density
                ! Go back to box coordinates for the central cell, which is transformed usually
                ! before sent to read_amr
                tempvar(:,:) = var(:,:)
                v = tempvar(0,varIDs%vx:varIDs%vz)
                call rotate_vector(v,transpose(trans_matrix))
                v = v + reg%bulk_velocity
                tempvar(0,varIDs%vx:varIDs%vz) = v
                totP(0) = totP(0) + 0.5d0 * var(0,varIDs%density) *sum(tempvar(0,varIDs%vx:varIDs%vz)**2d0)
            end if

            if (sim%mhd) then
                ! Magnetic energy density
                B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),&
                    & (var(0,varIDs%Bly)+var(0,varIDs%Bry)),&
                    & (var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
                totP(0) = totP(0) +  0.5 * (B.DOT.B)
            end if
            
            if (sim%cr) then
                ! CR energy density
                totP(0) = totP(0) + var(0,varIDs%cr_pressure) / (gamma_crs - 1d0)
            end if
            value = (abs(var(0,varIDs%thermal_pressure)) / (gamma_gas - 1d0)) / totP(0)
        case ('chi_cr')
            value = var(0,varIDs%cr_pressure) / max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density))
        case ('cr_energy_specific')
            ! Specific CR energy, computed as CR_energydensity*volume/cell mass
            value = (var(0,varIDs%cr_pressure) / (4D0/3d0 - 1d0)) / var(0,varIDs%density)
        case ('cr_temperature_eff')
            value = var(0,varIDs%cr_pressure) / var(0,varIDs%density)
        case ('gamma_ray_luminosity')
            ! Integrate hadronic gamma-ray volumetric luminosity. This assumes the CR spectral
            ! index at 5 GeV of the MW measured by Fermi LAT (Casandjian 2015)
            ! and the CR energy density in the solar neighbourhood measure by Voyager 2
            ! (Boschini et al. 2020)
            nH  = (1d0 - var(0,varIDs%metallicity)) * var(0,varIDs%density) * sim%nH ![H/cm^3]
            ecr = (var(0,varIDs%cr_pressure) / (gamma_crs - 1d0) ) * sim%unit_p ! [erg/cm^3]
            value = LgammaH * nH * (ecr / ecr_sun) / (sim%unit_p / sim%unit_t) * (dx * (dx * dx))
        case ('grad_crp')
            ! Magnitude of CR pressure gradient
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            value = magnitude(v)
        case ('torque_crp_specific')
            ! Magnitude of torque due to CR pressure gradient
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = (r_sphere(x) / var(0,varIDs%density)) * magnitude(temp_basis%u(1) * v)
        case ('grad_crp_dotmag')
            ! Dot product of CR pressure gradient and magnetic field unit vector
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)

            B = 0.5d0 * (/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            B = B / magnitude(B)
            v = v / magnitude(v)
            value = dacos(v.DOT.B)
        case ('grad_density_dotmag')
            ! Dot product of density gradient and magnetic field unit vector
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%density) - var(1,varIDs%density)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%density) - var(3,varIDs%density)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%density) - var(5,varIDs%density)) / (dxright + dxleft)

            B = 0.5d0 * (/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            B = B / magnitude(B)
            v = v / magnitude(v)
            value = dacos(v.dot.B)
        case ('absgrad_density_dotmag')
            ! Dot product of density gradient and magnetic field unit vector
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%density) - var(1,varIDs%density)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%density) - var(3,varIDs%density)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%density) - var(5,varIDs%density)) / (dxright + dxleft)

            B = 0.5d0 * (/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            B = B / magnitude(B)
            v = v / magnitude(v)
            value = abs(v.dot.B)
        case ('grad_crprsphere')
            ! CR pressure gradient in the radial direction
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = (v.DOT.temp_basis%u(1))
        case ('gradscale_crprsphere')
            ! CR pressure gradient scale in the radial direction
            ! This is defined as Pcr/grad(Pcr)
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = (v.DOT.temp_basis%u(1))
            value = abs(var(0,varIDs%cr_pressure)/value)
        case ('gradscale_crp')
            ! CR pressure gradient scale
            ! This is defined as Pcr/grad(Pcr)
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = magnitude(v)
            value = abs(var(0,varIDs%cr_pressure)/value)
        case ('diffusion_speed')
            ! CR diffusion speed
            ! This is defined as Dcr/Lcr, with Lcr the CR pressure gradient scale
            ! CR pressure gradient scale
            ! This is defined as Pcr/grad(Pcr)
            tempvar(:,:) = var(:,:)
            if (tempvar(0,varIDs%cr_pressure).le.1d-10) then
                tempvar(0,varIDs%cr_pressure) = 1d-10
            end if
            if (any(tempvar(:,varIDs%cr_pressure).le.1d-10)) then
                do i = 1, 6
                    if (tempvar(i,varIDs%cr_pressure).le.1d-10) tempvar(i,varIDs%cr_pressure) = 1d-10
                end do
            end if
            B = 0.5d0 * (/(tempvar(0,varIDs%Blx)+tempvar(0,varIDs%Brx)),(tempvar(0,varIDs%Bly)+tempvar(0,varIDs%Bry)),(tempvar(0,varIDs%Blz)+tempvar(0,varIDs%Brz))/)

            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (tempvar(2,varIDs%cr_pressure) - tempvar(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (tempvar(4,varIDs%cr_pressure) - tempvar(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (tempvar(6,varIDs%cr_pressure) - tempvar(5,varIDs%cr_pressure)) / (dxright + dxleft)
            B = B / magnitude(B)
            value = Dcr / abs(tempvar(0,varIDs%cr_pressure)/(abs(v.DOT.B))) !abs(tempvar(0,varIDs%cr_pressure)/(magnitude(v))) ! !
            ! if (value>1d10) print*,value,tempvar(:,varIDs%cr_pressure),tempvar(:,varIDs%thermal_pressure),abs(v.DOT.B)
        case ('alfvendiff_ratio')
            ! Ratio of Alfven to diffusion speed
            tempvar(:,:) = var(:,:)
            if (tempvar(0,varIDs%cr_pressure).le.1d-10) then
                tempvar(0,varIDs%cr_pressure) = 1d-10
            end if
            if (any(tempvar(:,varIDs%cr_pressure).le.1d-10)) then
                do i = 1, 6
                    if (tempvar(i,varIDs%cr_pressure).le.1d-10) tempvar(i,varIDs%cr_pressure) = 1d-10
                end do
            end if
            B = 0.5d0 * (/(tempvar(0,varIDs%Blx)+tempvar(0,varIDs%Brx)),(tempvar(0,varIDs%Bly)+tempvar(0,varIDs%Bry)),(tempvar(0,varIDs%Blz)+tempvar(0,varIDs%Brz))/)
            vA = magnitude(B) / sqrt(tempvar(0,varIDs%density)) * (gamma_crs / (gamma_crs-1d0))

            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (tempvar(2,varIDs%cr_pressure) - tempvar(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (tempvar(4,varIDs%cr_pressure) - tempvar(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (tempvar(6,varIDs%cr_pressure) - tempvar(5,varIDs%cr_pressure)) / (dxright + dxleft)
            B = B / magnitude(B)
            value = Dcr * abs(v.DOT.B) / tempvar(0,varIDs%cr_pressure)
            value = vA / (value + vA)
        case ('grad_crpx')
            ! Gradient of CR pressure in the x direction
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            value = v%x
        case ('grad_crpy')
            ! Gradient of CR pressure in the y direction
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            value = v%y
        case ('grad_crpz')
            ! Gradient of CR pressure in the z direction
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            value = v%z
        case ('streaming_heating')
            ! CR streaming heating
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            vst = (B / sqrt(var(0,varIDs%density)))
            value = abs(vst .DOT. v)
        case ('stheatcooling_ratio')
            ! Ratio of streaming heating rate to gas cooling rate
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
            vst = (B / sqrt(var(0,varIDs%density)))
            value = abs(vst .DOT. v)

            ! Net cooling rate taken from the cooling table in RAMSES output
            T = var(0,varIDs%thermal_pressure) / var(0,varIDs%density) * sim%T2
            ! TODO: This a fix only for some messed up CRMHD simulations!
            if (T<15d0) T = 15d0
            nH = var(0,varIDs%density) * sim%nH
            call solve_net_cooling(nH,T,var(0,varIDs%metallicity)/2D-2,lambda,lambda_prime)
            lambda_co = ((lambda * nH) * nH) * ((sim%unit_t**3)/(sim%unit_d*(sim%unit_l**2)))
            value = abs(value / lambda_co)
        case ('total_coolingtime')
            !TODO: Check units!
            ! Net cooling rate taken from the cooling table in RAMSES output
            T = var(0,varIDs%thermal_pressure) / var(0,varIDs%density) * sim%T2 ! This is actually T/mu
            if (T<15) T = 15
            nH = var(0,varIDs%density) * sim%nH
            call solve_net_cooling(nH,T,var(0,varIDs%metallicity)/2D-2,lambda,lambda_prime)
            lambda_co = (lambda * nH) * nH ! [erg/s/cm^3]

            if (sim%cr .and. sim%cr_st .and. sim%cr_heat) then
                ! CR streaming heating
                dxright = dx; dxleft = dx
                if (son(1) .eq. 0) dxright = dxright * 1.5D0
                if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
                v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
                dxright = dx; dxleft = dx
                if (son(3) .eq. 0) dxright = dxright * 1.5D0
                if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
                v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
                dxright = dx; dxleft = dx
                if (son(5) .eq. 0) dxright = dxright * 1.5D0
                if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
                v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)

                B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
                vst = (B / sqrt(var(0,varIDs%density)))
                lambda_st = abs(vst .DOT. v) * (sim%unit_p / sim%unit_t) ! [erg/s/cm^3]
            else
                lambda_st = 0D0
            end if

            if (sim%cr) then
                ! Cosmic rays hadronic and Coulomb heating from Guo&Ho(2008)
                ! (Assume fully ionised gas)
                ! TODO: Update for RT! 
                lambda = 2.63d-16 ! [erg/s/cm^3]
                ne = var(0,varIDs%density) * sim%unit_d / mHydrogen 
                ecr = var(0,varIDs%cr_pressure) / (gamma_crs - 1d0)
                ecr = ecr * (sim%unit_d * ((sim%unit_l/sim%unit_t)**2))
                lambda_cr = lambda * ne * ecr ! [erg/s/cm^3]
            else
                lambda_cr = 0D0
            end if

            ! Thermal energy
            ! TODO: This a fix only for some messed up CRMHD simulations!
            if (T<15) then
                value = (Tmin * var(0,varIDs%density)) / (gamma_gas - 1d0)
            else
                value = var(0,varIDs%thermal_pressure) / (gamma_gas - 1d0)
            end if

            if ((lambda_co - lambda_st - lambda_cr)<0d0) then
                ! In the case the gas is effectively heated, the cooling time is basically large
                value = 4.34d18 / sim%unit_t
            else
                ! Convert to code time unit
                value = value / ((lambda_co - lambda_st - lambda_cr)/(sim%unit_p / sim%unit_t))
            end if
        case ('xHII')
            ! Hydrogen ionisation fraction
            value = var(0,varIDs%xHII)
        case ('xHeII')
            ! Helium first ionisation fraction
            value = var(0,varIDs%xHeII)
        case ('xHeIII')
            ! Helium second ionisation fraction
            value = var(0,varIDs%xHeIII)
        case ('momentum_x')
            ! Linear momentum in the x direction as density*volume*corrected_velocity_x
            value = ((var(0,varIDs%density) * (dx*dx)) * dx) * var(0,varIDs%vx)
        case ('momentum_y')
            ! Linear momentum in the y direction density*volume*corrected_velocity_y
            value = ((var(0,varIDs%density) * (dx*dx)) * dx) * var(0,varIDs%vy)
        case ('momentum_z')
            ! Linear momentum in the z direction density*volume*corrected_velocity_z
            value = ((var(0,varIDs%density) * (dx*dx)) * dx) * var(0,varIDs%vz)
        case ('momentum')
            ! Magnitude of linear momentum, using corrected velocity
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            value = ((var(0,varIDs%density) * (dx*dx)) * dx) * &
                    & magnitude(v)
        case ('momentum_sphere_r')
            ! Linear momentum in the spherical radial direction
            ! 1. Dot product of velocity vector with spherical r
            !    unit vector
            ! 2. Multiply by mass of cell
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = (var(0,varIDs%density) * (dx*dx)) * dx * (v .DOT. temp_basis%u(1))
        case ('absmomentum_sphere_r')
            ! Linear momentum in the spherical radial direction (absolute value)
            ! 1. Dot product of velocity vector with spherical r
            !    unit vector
            ! 2. Multiply by mass of cell
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = (var(0,varIDs%density) * (dx*dx)) * dx * (v .DOT. temp_basis%u(1))
            value = abs(value)
        case ('momentum_cyl_z')
            ! Linear momentum in the cylindrical z direction
            ! 1. Dot product of velocity vector with cylindrical z
            !    unit vector
            ! 2. Multiply by mass of cell
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call cylindrical_basis_from_cartesian(x,temp_basis)
            value = (var(0,varIDs%density) * (dx*dx)) * dx * (v .DOT. temp_basis%u(3))
        case ('ang_momentum_x')
            ! Corrected angular momentum in the x direction
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            value = ((var(0,varIDs%density) * (dx*dx)) * dx) * (x%y * v%z &
                        &- v%y * x%z)
        case ('ang_momentum_y')
            ! Corrected angular momentum in the y direction
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            value = ((var(0,varIDs%density) * (dx*dx)) * dx) * (x%z*v%x &
                        &- v%z*x%x)
        case ('ang_momentum_z')
            ! Corrected angular momentum in the z direction
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            value = ((var(0,varIDs%density) * (dx*dx)) * dx) * (x%x*v%y &
                        &- v%x*x%y)
        case ('ang_momentum_specific_x')
            ! Corrected specific angular momentum in the x direction
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            value = x%y * v%z - v%y * x%z
        case ('ang_momentum_specific_y')
            ! Corrected specific angular momentum in the y direction
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            value = x%z*v%x - v%z*x%x
        case ('ang_momentum_specific_z')
            ! Corrected specific angular momentum in the z direction
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            value = x%x*v%y - v%x*x%y
        case ('ang_momentum')
            ! Corrected magnitude of angular momentum
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            L = x * v
            value = ((var(0,varIDs%density) * (dx*dx)) * dx) * magnitude(L)
        case ('ang_momentum_specific')
            ! Corrected magnitude of specific angular momentum
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            L = x * v
            value = magnitude(L)
        case ('massflow_rate_sphere_r')
            ! Mass flow rate through the cell in the radial direction
            ! Mass per unit time
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = (var(0,varIDs%density) * (dx*dx)) * (v .DOT. temp_basis%u(1))
        case ('massflux_rate_sphere_r')
            ! Mass flux through the cell in the radial direction
            ! Mass per unit time per unit surface
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = var(0,varIDs%density) * (v .DOT. temp_basis%u(1))
        case ('grav_potential')
            ! Gravitational potential
            value = grav_var(0,1)
        case ('grav_gx')
            ! Gravitational acceleration in the x direction
            value = grav_var(0,2)
        case ('grav_gy')
            ! Gravitational acceleration in the y direction
            value = grav_var(0,3)
        case ('grav_gz')
            ! Gravitational acceleration in the z direction
            value = grav_var(0,4)
        case ('grav_grsphere')
            ! Gravitational acceleration in the radial direction
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = abs(B .DOT. temp_basis%u(1))
        case ('grav_frsphere')
            ! Total gravitational force in the radial direction
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = var(0,varIDs%density) * (B .DOT. temp_basis%u(1))
        case ('grav_fz')
            ! Gravitational force in the z direction
            B = grav_var(0,2:4)
            call rotate_vector(B,trans_matrix)
            value = var(0,varIDs%density) * B%z
        case ('escape_velocity')
            ! Local gravitational escape velocity
            value = sqrt(2d0*abs(grav_var(0,1)))
        case ('circular_velocity')
            ! Local circular velocity
            value = sqrt(abs(grav_var(0,1)))
        case ('grav_centfrsphere')
            ! Ratio of centripetal acceleration to gravitational acceleration
            ! in the spherical r direction
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call spherical_basis_from_cartesian(x,temp_basis)
            v = v - (v.DOT.temp_basis%u(1)) * temp_basis%u(1)
            value = magnitude(v)**2d0 / r_sphere(x)
            B = grav_var(0,2:4)
            value = (magnitude(v)**2d0 / r_sphere(x)) / (B .DOT. temp_basis%u(1))
        case ('torque_grav_specific')
            ! Magnitude of the specific torque due to gravity
            B = grav_var(0,2:4)
            call rotate_vector(B,trans_matrix)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = r_sphere(x) * magnitude(temp_basis%u(1) * B)
        case ('grav_crpf')
            ! Ratio of CR pressure gradient and gravitational acceleration
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            B = grav_var(0,2:4)
            value = magnitude(v) / (var(0,varIDs%density) * magnitude(B))
        case ('grav_crpfz')
            ! Ratio of CR pressure gradient and gravitational acceleration
            ! in the z direction
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            value = -v%z / (var(0,varIDs%density) * grav_var(0,4))
        case ('grav_therpfz')
            ! Ratio of thermal pressure gradient and gravitational acceleration
            ! in the z direction
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (max(var(2,varIDs%thermal_pressure), Tmin*var(2,varIDs%density)) - &
                    & max(var(1,varIDs%thermal_pressure),Tmin*var(1,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (max(var(4,varIDs%thermal_pressure), Tmin*var(4,varIDs%density)) - &
                    & max(var(3,varIDs%thermal_pressure),Tmin*var(3,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (max(var(6,varIDs%thermal_pressure), Tmin*var(6,varIDs%density)) - &
                    & max(var(5,varIDs%thermal_pressure),Tmin*var(5,varIDs%density))) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            value = -v%z / (var(0,varIDs%density) * grav_var(0,4))
        case ('grav_therpfrsphere')
            ! Ratio of thermal pressure gradient and gravitational acceleration
            ! in the spherical r direction
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (max(var(2,varIDs%thermal_pressure), Tmin*var(2,varIDs%density)) - &
                    & max(var(1,varIDs%thermal_pressure),Tmin*var(1,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (max(var(4,varIDs%thermal_pressure), Tmin*var(4,varIDs%density)) - &
                    & max(var(3,varIDs%thermal_pressure),Tmin*var(3,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (max(var(6,varIDs%thermal_pressure), Tmin*var(6,varIDs%density)) - &
                    & max(var(5,varIDs%thermal_pressure),Tmin*var(5,varIDs%density))) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = -(v.DOT.temp_basis%u(1)) / (var(0,varIDs%density) * (B .DOT. temp_basis%u(1)))
        case ('grav_therpfrspherepos')
            ! Ratio of thermal pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY POSITIVE
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (max(var(2,varIDs%thermal_pressure), Tmin*var(2,varIDs%density)) - &
                    & max(var(1,varIDs%thermal_pressure),Tmin*var(1,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (max(var(4,varIDs%thermal_pressure), Tmin*var(4,varIDs%density)) - &
                    & max(var(3,varIDs%thermal_pressure),Tmin*var(3,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (max(var(6,varIDs%thermal_pressure), Tmin*var(6,varIDs%density)) - &
                    & max(var(5,varIDs%thermal_pressure),Tmin*var(5,varIDs%density))) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = -(v.DOT.temp_basis%u(1)) / (var(0,varIDs%density) * (B .DOT. temp_basis%u(1)))
            if (value < 0) value = 0d0
        case ('grav_therpfrsphereneg')
            ! Ratio of thermal pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY NEGATIVE
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (max(var(2,varIDs%thermal_pressure), Tmin*var(2,varIDs%density)) - &
                    & max(var(1,varIDs%thermal_pressure),Tmin*var(1,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (max(var(4,varIDs%thermal_pressure), Tmin*var(4,varIDs%density)) - &
                    & max(var(3,varIDs%thermal_pressure),Tmin*var(3,varIDs%density))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (max(var(6,varIDs%thermal_pressure), Tmin*var(6,varIDs%density)) - &
                    & max(var(5,varIDs%thermal_pressure),Tmin*var(5,varIDs%density))) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = -(v.DOT.temp_basis%u(1)) / (var(0,varIDs%density) * (B .DOT. temp_basis%u(1)))
            if (value > 0) then
                value = 0d0
            else
                value = - value
            end if
        case ('grav_crpfrsphere')
            ! Ratio of CR pressure gradient and gravitational acceleration
            ! in the spherical r direction
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = -(v.DOT.temp_basis%u(1)) / (var(0,varIDs%density) * (B .DOT. temp_basis%u(1)))
        case ('grav_crpfrspherepos')
            ! Ratio of CR pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY POSITIVE
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = -(v.DOT.temp_basis%u(1)) / (var(0,varIDs%density) * (B .DOT. temp_basis%u(1)))
            if (value < 0) value = 0d0
        case ('grav_crpfrsphereneg')
            ! Ratio of CR pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY NEGATIVE
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = -(v.DOT.temp_basis%u(1)) / (var(0,varIDs%density) * (B .DOT. temp_basis%u(1)))
            if (value > 0) then
                value = 0d0
            else
                value = - value
            end if
        case ('torque_magp_specific')
            ! Magnitude of the specific torque due to magnetic pressure gradient
            totP(:) = 0d0
            do i=1,amr%twondim
                B = 0.5 *(/(var(i,varIDs%Blx)+var(i,varIDs%Brx)),&
                    & (var(i,varIDs%Bly)+var(i,varIDs%Bry)),&
                    & (var(i,varIDs%Blz)+var(i,varIDs%Brz))/)
                value = 0.5 * (B.DOT.B)
                totP(i) = value
            end do
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (totP(2) - totP(1)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (totP(4) - totP(3)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (totP(6) - totP(5)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = (r_sphere(x) / var(0,varIDs%density)) * magnitude(temp_basis%u(1) * v)
        case ('grad_magprsphere')
            ! Magnetic pressure gradient
            ! in the spherical r direction
            totP(:) = 0d0
            do i=1,amr%twondim
                B = 0.5 *(/(var(i,varIDs%Blx)+var(i,varIDs%Brx)),&
                    & (var(i,varIDs%Bly)+var(i,varIDs%Bry)),&
                    & (var(i,varIDs%Blz)+var(i,varIDs%Brz))/)
                value = 0.5 * (B.DOT.B)
                totP(i) = value
            end do
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (totP(2) - totP(1)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (totP(4) - totP(3)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (totP(6) - totP(5)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = v.DOT.temp_basis%u(1)
        case ('grav_magpfrsphere')
            ! Ratio of magnetic pressure gradient and gravitational acceleration
            ! in the spherical r direction
            totP(:) = 0d0
            do i=1,amr%twondim
                B = 0.5 *(/(var(i,varIDs%Blx)+var(i,varIDs%Brx)),&
                    & (var(i,varIDs%Bly)+var(i,varIDs%Bry)),&
                    & (var(i,varIDs%Blz)+var(i,varIDs%Brz))/)
                value = 0.5 * (B.DOT.B)
                totP(i) = value
            end do
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (totP(2) - totP(1)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (totP(4) - totP(3)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (totP(6) - totP(5)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = -(v.DOT.temp_basis%u(1)) / (var(0,varIDs%density) * (B .DOT. temp_basis%u(1)))
        case ('grav_magpfrspherepos')
            ! Ratio of magnetic pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY POSITIVE
            totP(:) = 0d0
            do i=1,amr%twondim
                B = 0.5 *(/(var(i,varIDs%Blx)+var(i,varIDs%Brx)),&
                    & (var(i,varIDs%Bly)+var(i,varIDs%Bry)),&
                    & (var(i,varIDs%Blz)+var(i,varIDs%Brz))/)
                value = 0.5 * (B.DOT.B)
                totP(i) = value
            end do
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (totP(2) - totP(1)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (totP(4) - totP(3)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (totP(6) - totP(5)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = -(v.DOT.temp_basis%u(1)) / (var(0,varIDs%density) * (B .DOT. temp_basis%u(1)))
            if (value < 0) value = 0d0
        case ('grav_magpfrsphereneg')
            ! Ratio of magnetic pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY NEGATIVE
            totP(:) = 0d0
            do i=1,amr%twondim
                B = 0.5 *(/(var(i,varIDs%Blx)+var(i,varIDs%Brx)),&
                    & (var(i,varIDs%Bly)+var(i,varIDs%Bry)),&
                    & (var(i,varIDs%Blz)+var(i,varIDs%Brz))/)
                value = 0.5 * (B.DOT.B)
                totP(i) = value
            end do
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (totP(2) - totP(1)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (totP(4) - totP(3)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (totP(6) - totP(5)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = -(v.DOT.temp_basis%u(1)) / (var(0,varIDs%density) * (B .DOT. temp_basis%u(1)))
            if (value > 0) then
                value = 0d0
            else
                value = - value
            end if
        case ('torque_totp_specific')
            ! Magnitude of the specific torque due to total pressure gradient
            totP(:) = 0d0
            do i=1,amr%twondim
                if (sim%hydro) totP(i) = max(var(i,varIDs%thermal_pressure), Tmin*var(i,varIDs%density))
                if (sim%mhd) then
                    B = 0.5 *(/(var(i,varIDs%Blx)+var(i,varIDs%Brx)),&
                        & (var(i,varIDs%Bly)+var(i,varIDs%Bry)),&
                        & (var(i,varIDs%Blz)+var(i,varIDs%Brz))/)
                    totP(i) = totP(i) + 0.5 * (B.DOT.B)
                end if
                if (sim%cr) totP(i) = totP(i) + var(i,varIDs%cr_pressure)
            end do
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (totP(2) - totP(1)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (totP(4) - totP(3)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (totP(6) - totP(5)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = (r_sphere(x) / var(0,varIDs%density)) * magnitude(temp_basis%u(1) * v)
        case ('grav_totpfrsphere')
            ! Ratio of total pressure gradient and gravitational acceleration
            ! in the spherical r direction
            totP(:) = 0d0
            do i=1,amr%twondim
                if (sim%hydro) totP(i) = var(i,varIDs%thermal_pressure)
                if (sim%mhd) then
                    B = 0.5 *(/(var(i,varIDs%Blx)+var(i,varIDs%Brx)),&
                        & (var(i,varIDs%Bly)+var(i,varIDs%Bry)),&
                        & (var(i,varIDs%Blz)+var(i,varIDs%Brz))/)
                    totP(i) = totP(i) + 0.5 * (B.DOT.B)
                end if
                if (sim%cr) totP(i) = totP(i) + var(i,varIDs%cr_pressure)
            end do
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (totP(2) - totP(1)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (totP(4) - totP(3)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (totP(6) - totP(5)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = -(v.DOT.temp_basis%u(1)) / (var(0,varIDs%density) * (B .DOT. temp_basis%u(1)))
        case ('grav_totpfrspherepos')
            ! Ratio of total pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY POSITIVE
            totP(:) = 0d0
            do i=1,amr%twondim
                if (sim%hydro) totP(i) = var(i,varIDs%thermal_pressure)
                if (sim%mhd) then
                    B = 0.5 *(/(var(i,varIDs%Blx)+var(i,varIDs%Brx)),&
                        & (var(i,varIDs%Bly)+var(i,varIDs%Bry)),&
                        & (var(i,varIDs%Blz)+var(i,varIDs%Brz))/)
                    totP(i) = totP(i) + 0.5 * (B.DOT.B)
                end if
                if (sim%cr) totP(i) = totP(i) + var(i,varIDs%cr_pressure)
            end do
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (totP(2) - totP(1)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (totP(4) - totP(3)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (totP(6) - totP(5)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = -(v.DOT.temp_basis%u(1)) / (var(0,varIDs%density) * (B .DOT. temp_basis%u(1)))
            if (value < 0) value = 0d0
        case ('grav_totpfrsphereneg')
            ! Ratio of total pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY NEGATIVE
            totP(:) = 0d0
            do i=1,amr%twondim
                if (sim%hydro) totP(i) = var(i,varIDs%thermal_pressure)
                if (sim%mhd) then
                    B = 0.5 *(/(var(i,varIDs%Blx)+var(i,varIDs%Brx)),&
                        & (var(i,varIDs%Bly)+var(i,varIDs%Bry)),&
                        & (var(i,varIDs%Blz)+var(i,varIDs%Brz))/)
                    totP(i) = totP(i) + 0.5 * (B.DOT.B)
                end if
                if (sim%cr) totP(i) = totP(i) + var(i,varIDs%cr_pressure)
            end do
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (totP(2) - totP(1)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (totP(4) - totP(3)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (totP(6) - totP(5)) / (dxright + dxleft)
            call rotate_vector(v,trans_matrix)
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = -(v.DOT.temp_basis%u(1)) / (var(0,varIDs%density) * (B .DOT. temp_basis%u(1)))
            if (value > 0) then
                value = 0d0
            else
                value = - value
            end if
        case ('freefall_time')
            ! Free-fall time
            B = grav_var(0,2:4)
            call spherical_basis_from_cartesian(x,temp_basis)
            fg_r = B .DOT. temp_basis%u(1)
            r = r_sphere(x)
            ! if (fg_r > 0d0) then
            !     ! In the case the acceleration is outward, we just set the free-fall
            !     ! timesclae to a very large value
            !     value = 4.34d18 / sim%unit_t
            ! else
            !     value = sqrt((2d0*r)/abs(fg_r))
            ! end if
            value = sqrt((2d0*r)/abs(fg_r))
        case ('inflow_time')
            ! Radial inflow timescale, defined as t_inflow = r / v_r(inflow)
            v = (/var(0,varIDs%vx),var(0,varIDs%vy),var(0,varIDs%vz)/)
            call spherical_basis_from_cartesian(x,temp_basis)
            value = v.DOT.temp_basis%u(1)
            if (value > 0d0) then
                ! In the case the gas is outflowing, the timescale is just larger
                ! than the age of the Universe in code time units
                value = 4.34d18 / sim%unit_t
            else
                r = r_sphere(x)
                value = r / abs(value)
            end if
        case ('compression_time')
            ! Compression timescale, as defined in Birnboim and Dekel (2003)
            ! 1. Velocity divergence
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%vx) - var(1,varIDs%vx)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%vy) - var(3,varIDs%vy)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%vz) - var(5,varIDs%vz)) / (dxright + dxleft)
            value = v%x + v%y + v%z
            if (value > 0d0) then
                ! In case the gas is diverging instead of converging, the timescale is
                ! just very large
                value = 4.34d18 / sim%unit_t
            else
                value = (21d0/5d0) / abs(value)
            end if
        case ('compcooltime_ratio')
            ! Ratio of compression time to cooling time
            ! 1. Velocity divergence
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(2,varIDs%vx) - var(1,varIDs%vx)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(4,varIDs%vy) - var(3,varIDs%vy)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(6,varIDs%vz) - var(5,varIDs%vz)) / (dxright + dxleft)
            tcomp = v%x + v%y + v%z
            if (tcomp > 0d0) then
                ! In case the gas is diverging instead of converging, the timescale is
                ! just very large
                tcomp = 4.34d18 / sim%unit_t
            else
                tcomp = (21d0/5d0) / abs(tcomp)
            end if
            ! 2. Cooling time
            T = var(0,varIDs%thermal_pressure) / var(0,varIDs%density) * sim%T2 ! This is actually T/mu
            if (T<15) T = 15
            nH = var(0,varIDs%density) * sim%nH
            call solve_net_cooling(nH,T,var(0,varIDs%metallicity)/2D-2,lambda,lambda_prime)
            lambda_co = (lambda * nH) * nH ! [erg/s/cm^3]

            if (sim%cr .and. sim%cr_st .and. sim%cr_heat) then
                ! CR streaming heating
                dxright = dx; dxleft = dx
                if (son(1) .eq. 0) dxright = dxright * 1.5D0
                if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
                v%x = (var(2,varIDs%cr_pressure) - var(1,varIDs%cr_pressure)) / (dxright + dxleft)
                dxright = dx; dxleft = dx
                if (son(3) .eq. 0) dxright = dxright * 1.5D0
                if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
                v%y = (var(4,varIDs%cr_pressure) - var(3,varIDs%cr_pressure)) / (dxright + dxleft)
                dxright = dx; dxleft = dx
                if (son(5) .eq. 0) dxright = dxright * 1.5D0
                if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
                v%z = (var(6,varIDs%cr_pressure) - var(5,varIDs%cr_pressure)) / (dxright + dxleft)

                B = 0.5 *(/(var(0,varIDs%Blx)+var(0,varIDs%Brx)),(var(0,varIDs%Bly)+var(0,varIDs%Bry)),(var(0,varIDs%Blz)+var(0,varIDs%Brz))/)
                vst = (B / sqrt(var(0,varIDs%density)))
                lambda_st = abs(vst .DOT. v) * (sim%unit_p / sim%unit_t) ! [erg/s/cm^3]
            else
                lambda_st = 0D0
            end if

            if (sim%cr) then
                ! Cosmic rays hadronic and Coulomb heating from Guo&Ho(2008)
                ! (Assume fully ionised gas)
                ! TODO: Update for RT! 
                lambda = 2.63d-16 ! [erg/s/cm^3]
                ne = var(0,varIDs%density) * sim%unit_d / mHydrogen 
                ecr = var(0,varIDs%cr_pressure) / (gamma_crs - 1d0)
                ecr = ecr * (sim%unit_d * ((sim%unit_l/sim%unit_t)**2))
                lambda_cr = lambda * ne * ecr ! [erg/s/cm^3]
            else
                lambda_cr = 0D0
            end if

            ! Thermal energy
            ! TODO: This a fix only for some messed up CRMHD simulations!
            if (T<15) then
                tcool = (Tmin * var(0,varIDs%density)) / (gamma_gas - 1d0)
            else
                tcool = var(0,varIDs%thermal_pressure) / (gamma_gas - 1d0)
            end if

            if ((lambda_co - lambda_st - lambda_cr)<0d0) then
                ! In the case the gas is effectively heated, the cooling time is basically large
                tcool = 4.34d18 / sim%unit_t
            else
                ! Convert to code time unit
                tcool = tcool / ((lambda_co - lambda_st - lambda_cr)/(sim%unit_p / sim%unit_t))
            end if
            value =  tcomp / tcool
        case ('eff_FKmag')
            ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
            star_maker = 'FKmag'
            ! Go back to box coordinates for the central cell, which is transformed usually
            ! before sent to read_amr
            ! Converging flow check
            tempvar(:,:) = var(:,:)
            v = tempvar(0,varIDs%vx:varIDs%vz)
            call rotate_vector(v,transpose(trans_matrix))
            v = v + reg%bulk_velocity
            tempvar(0,varIDs%vx:varIDs%vz) = v
            value  = sf_eff(reg,dx,x,tempvar,star_maker,.true.)
        case ('eff_FKmagnocr')
            ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
            ! In this case the contribution of CR pressure is ignore in the computation
            ! of the effective sound speed
            star_maker = 'FKmag'
            ! Go back to box coordinates for the central cell, which is transformed usually
            ! before sent to read_amr
            ! Converging flow check
            tempvar(:,:) = var(:,:)
            v = tempvar(0,varIDs%vx:varIDs%vz)
            call rotate_vector(v,transpose(trans_matrix))
            v = v + reg%bulk_velocity
            tempvar(0,varIDs%vx:varIDs%vz) = v
            value  = sf_eff(reg,dx,x,tempvar,star_maker,.false.)
        case ('eff_FK2')
            ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
            star_maker = 'FK2'
            ! Go back to box coordinates for the central cell, which is transformed usually
            ! before sent to read_amr
            ! Converging flow check
            tempvar(:,:) = var(:,:)
            v = tempvar(0,varIDs%vx:varIDs%vz)
            call rotate_vector(v,transpose(trans_matrix))
            v = v + reg%bulk_velocity
            tempvar(0,varIDs%vx:varIDs%vz) = v
            value  = sf_eff(reg,dx,x,tempvar,star_maker,.false.)
        case ('jeans_mtt')
            ! MTT model Jeans length
            star_maker = 'jeans_mtt'
            ! Go back to box coordinates for the central cell, which is transformed usually
            ! before sent to read_amr
            ! Converging flow check
            tempvar(:,:) = var(:,:)
            v = tempvar(0,varIDs%vx:varIDs%vz)
            call rotate_vector(v,transpose(trans_matrix))
            v = v + reg%bulk_velocity
            tempvar(0,varIDs%vx:varIDs%vz) = v
            value  = sf_eff(reg,dx,x,tempvar,star_maker,.true.)
        case ('jeansmtt_dx')
            ! MTT model Jeans length in cell units
            star_maker = 'jeansmtt_dx'
            ! Go back to box coordinates for the central cell, which is transformed usually
            ! before sent to read_amr
            ! Converging flow check
            tempvar(:,:) = var(:,:)
            v = tempvar(0,varIDs%vx:varIDs%vz)
            call rotate_vector(v,transpose(trans_matrix))
            v = v + reg%bulk_velocity
            tempvar(0,varIDs%vx:varIDs%vz) = v
            value  = sf_eff(reg,dx,x,tempvar,star_maker,.true.)
        case ('jeans_mttnocr')
            ! MTT model Jeans length without CR pressure support
            star_maker = 'jeans_mtt'
            ! Go back to box coordinates for the central cell, which is transformed usually
            ! before sent to read_amr
            ! Converging flow check
            tempvar(:,:) = var(:,:)
            v = tempvar(0,varIDs%vx:varIDs%vz)
            call rotate_vector(v,transpose(trans_matrix))
            v = v + reg%bulk_velocity
            tempvar(0,varIDs%vx:varIDs%vz) = v
            value  = sf_eff(reg,dx,x,tempvar,star_maker,.false.)
        case ('jeansmtt_dxnocr')
            ! MTT model Jeans length in cell units without CR pressure support
            star_maker = 'jeansmtt_dx'
            ! Go back to box coordinates for the central cell, which is transformed usually
            ! before sent to read_amr
            ! Converging flow check
            tempvar(:,:) = var(:,:)
            v = tempvar(0,varIDs%vx:varIDs%vz)
            call rotate_vector(v,transpose(trans_matrix))
            v = v + reg%bulk_velocity
            tempvar(0,varIDs%vx:varIDs%vz) = v
            value  = sf_eff(reg,dx,x,tempvar,star_maker,.false.)
        case ('neighbour_accuracy')
            ! This variable is used as a debugging method to check
            ! whether the nearest neighbour implementation is able to 
            ! recover correct gravitational acceleration
            
            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (grav_var(2,1) - grav_var(1,1)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (grav_var(4,1) - grav_var(3,1)) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (grav_var(6,1) - grav_var(5,1)) / (dxright + dxleft)
            B = grav_var(0,2:4)
            call rotate_vector(B,transpose(trans_matrix))
            value = -(B.DOT.v)/(magnitude(v)*magnitude(B))
        case ('sigma')
            ! Go back to box coordinates for the central cell, which is transformed usually
            ! before sent to read_amr
            ! Converging flow check
            tempvar(:,:) = var(:,:)
            v = tempvar(0,varIDs%vx:varIDs%vz)
            call rotate_vector(v,transpose(trans_matrix))
            v = v + reg%bulk_velocity
            tempvar(0,varIDs%vx:varIDs%vz) = v
            call cmp_sigma_turb(tempvar,value)
        case ('turbulent_pressure')
            ! Go back to box coordinates for the central cell, which is transformed usually
            ! before sent to read_amr
            tempvar(:,:) = var(:,:)
            v = tempvar(0,varIDs%vx:varIDs%vz)
            call rotate_vector(v,transpose(trans_matrix))
            v = v + reg%bulk_velocity
            tempvar(0,varIDs%vx:varIDs%vz) = v
            call cmp_sigma_turb(tempvar,sigma)
            value = var(0,varIDs%density) * sigma**2d0
        case ('mach_number')
            ! Mach number (velocity/sound speed) in the reference frame of the galaxy
            cs = sqrt(5D0/3d0 * (max(var(0,varIDs%thermal_pressure), Tmin*var(0,varIDs%density)) / var(0,varIDs%density)))
            v = tempvar(0,varIDs%vx:varIDs%vz)
            value = magnitude(v)/cs
        case default
            write(*,*)'Variable not supported: ',TRIM(varname)
            write(*,*)'Aborting!'
            stop
        end select

    end subroutine getvarvalue

    !---------------------------------------------------------------
    ! Subroutine: GET ETA SN
    !
    ! When the initial particle mass is not saved, the initial
    ! star particle mass needs to be computed by using the tags
    ! and the value of the parameter eta_sn in the namelist of the
    ! simulation
    !---------------------------------------------------------------
    subroutine get_eta_sn(repository)
        implicit none
        character(128), intent(in) :: repository
        integer :: i,status,ipos,jpos
        character(128) :: nomfich
        character(6)  ::  param
        real(dbl) ::  pvalue
        character(200) :: line
        logical :: ok,ok_found

        nomfich=TRIM(repository)//'/namelist.txt'
        
        inquire(file=nomfich, exist=ok) ! verify input file
        if (ok) then
            if (verbose) write(*,'(": Reading eta_sn from namelist.txt")')
            ok_found = .false.
            open(unit=15,file=nomfich,status='old',form='formatted')
            do
                read(15,'(A)',iostat=status)line
                if (status /= 0) exit
                param = line(1:6)
                if (param=='eta_sn') then
                    ipos = index(line,'=')
                    jpos = index(line,'!')
                    read(line(ipos+1:jpos-1),'(F10.0)')pvalue
                    sim%eta_sn = pvalue
                    if (verbose) write(*,*)': Found eta_sn=',sim%eta_sn
                    ok_found = .true.
                    exit
                end if
            end do
            if (.not.ok_found) then
                if (verbose) write(*,'(": Namelist does not have eta_sn: eta_sn set to 0.317")')
                sim%eta_sn = 0.317D0
            end if
        else
            if (verbose) write(*,'(": Namelist not found: eta_sn set to 0.2")')
            sim%eta_sn = 0.317D0 !2D-1
        end if
    end subroutine get_eta_sn
    !---------------------------------------------------------------
    ! Subroutine: GET Dcr
    !
    ! When we have a CRMHD simulation, look for the diffusion
    ! constant value as defined in the namelist
    ! TODO: This only works for NUT sims!
    !---------------------------------------------------------------
    subroutine get_Dcr(repository)
        implicit none
        character(128), intent(in) :: repository
        integer :: i,status,ipos,jpos
        character(128) :: nomfich
        character(6)  ::  param
        real(dbl) ::  pvalue
        character(200) :: line
        logical :: ok
        nomfich=TRIM(repository)//'/../CRiMHD+SfFb.nml'
        
        inquire(file=nomfich, exist=ok) ! verify input file
        if (ok) then
            if (verbose) write(*,'(": Reading Dcr from namelist.txt")')
            open(unit=15,file=nomfich,status='old',form='formatted')
            do
                read(15,'(A)',iostat=status)line
                if (status /= 0) exit
                param = line(1:3)
                if (param=='Dcr') then
                    ipos = index(line,'=')
                    jpos = index(line,'!')
                    read(line(ipos+1:jpos-1),'(F10.0)')pvalue
                    sim%Dcr = pvalue
                    if (verbose) write(*,*)': Found Dcr=',sim%Dcr
                    exit
                end if
            end do
        else
            if (verbose) write(*,'(": Namelist not found: Dcr set to 3.0d28")')
            sim%Dcr = 3.0d28
        end if
    end subroutine get_Dcr
    !---------------------------------------------------------------
    ! Function: CMP SIGMA TURB
    !
    ! This function obtains the local dispersion velocity of the gas
    ! in the same way it is used for the MTT models of star formation
    !----------------------------------------------------------------
    subroutine cmp_sigma_turb(var,sigma)
        use vectors
        use basis_representations
        use geometrical_regions
        use coordinate_systems

        implicit none
        real(dbl),dimension(0:amr%twondim,1:varIDs%nvar),intent(in) :: var
        real(dbl),intent(out)           :: sigma

        logical :: isConvergent
        real(dbl) :: d
        real(dbl),dimension(0:amr%twondim) ::darr,uarr,varr,warr
        real(dbl) :: divv,uavg,vavg,wavg,dtot
        real(dbl) :: px_div,py_div,pz_div,Jx,Jy,Jz,rho_local
        real(dbl) :: trgv,ul,ur

        ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
        ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor. 
        ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation 
        ! from neighbouring cell values and differentiate. 
        ! Get neighbor cells if they exist, otherwise use straight injection from local cell

        d = var(0,varIDs%density)
        darr = var(:,varIDs%density)

        uarr = var(:,varIDs%vx)
        varr = var(:,varIDs%vy)
        warr = var(:,varIDs%vz)
        divv  = (uarr(2)*darr(2)-uarr(1)*darr(1)) &
            & + (varr(4)*darr(4)-varr(3)*darr(3)) & 
            & + (warr(6)*darr(6)-warr(5)*darr(5))
        
        ! Average velocity
        dtot  = sum(darr)
        uavg  = sum(darr*uarr)/dtot
        vavg  = sum(darr*varr)/dtot
        wavg  = sum(darr*warr)/dtot
        ! Subtract the mean velocity field
        uarr(:) = uarr(:) - uavg
        varr(:) = varr(:) - vavg
        warr(:) = warr(:) - wavg
        ! Subtract the symmetric divergence field                    
        ! ex)  (---->,<--): only subtract (-->,<--): result (-->,0) 
        ! ex)  (<----,-->): only subtract (<--,-->): result (<--,0)
        px_div = min( abs(darr(1)*uarr(1)),abs(darr(2)*uarr(2)))
        py_div = min( abs(darr(3)*varr(3)),abs(darr(4)*varr(4)))
        pz_div = min( abs(darr(5)*warr(5)),abs(darr(6)*warr(6)))

        isConvergent = darr(2)*uarr(2) - darr(1)*uarr(1) < 0 
        if (isConvergent) then
            uarr(1) = uarr(1) - px_div/darr(1)
            uarr(2) = uarr(2) + px_div/darr(2)
        else ! comment out if you do not want to subtract outflows
            uarr(1) = uarr(1) + px_div/darr(1)
            uarr(2) = uarr(2) - px_div/darr(2)
        end if 

        isConvergent = darr(4)*varr(4) - darr(3)*varr(3) < 0
        if (isConvergent) then
            varr(3) = varr(3) - py_div/darr(3)
            varr(4) = varr(4) + py_div/darr(4)
        else ! comment out if you do not want to subtract outflows
            varr(3) = varr(3) + py_div/darr(3)
            varr(4) = varr(4) - py_div/darr(4)
        end if

        isConvergent = darr(6)*warr(6) - darr(5)*warr(5) < 0
        if (isConvergent) then 
            warr(5) = warr(5) - pz_div/darr(5)
            warr(6) = warr(6) + pz_div/darr(6)
        else ! comment out if you do not want to subtract outflows
            warr(5) = warr(5) + pz_div/darr(5)
            warr(6) = warr(6) - pz_div/darr(6)
        end if

        ! subtract the rotational velocity field (x-y) plane
        ! ^y       <-        |4|        |-u|
        ! |       |  |     |1| |2|   |-v|  |+v|
        ! --->x    ->        |3|        |+u|
        Jz  = - varr(1)*darr(1) + varr(2)*darr(2) &
            &   + uarr(3)*darr(3) - uarr(4)*darr(4)
        Jz  = Jz / 4.0

        varr(1) = varr(1) + Jz/darr(1) 
        varr(2) = varr(2) - Jz/darr(2) 
        uarr(3) = uarr(3) - Jz/darr(3)
        uarr(4) = uarr(4) + Jz/darr(4)

        ! subtract the rotational velocity field (y-z) plane
        ! ^z       <-        |6|        |-v|  
        ! |       |  |     |3| |4|   |-w|  |+w|
        ! --->y    ->        |5|        |+v|
        Jx  = - warr(3)*darr(3) + warr(4)*darr(4) &
            &   + varr(5)*darr(5) - varr(6)*darr(6)
        Jx  = Jx / 4.0

        warr(3) = warr(3) + Jx/darr(3) 
        warr(4) = warr(4) - Jx/darr(4) 
        varr(5) = varr(5) - Jx/darr(5)
        varr(6) = varr(6) + Jx/darr(6)

        ! subtract the rotational velocity field (x-z) plane
        ! ^z       ->        |6|        |+u|  
        ! |       |  |     |1| |2|   |+w|  |-w|
        ! --->x    <-        |5|        |-u|
        Jy  = + warr(1)*darr(1) - warr(2)*darr(2) &
            &   - uarr(5)*darr(5) + uarr(6)*darr(6)
        Jy  = Jy / 4.0

        warr(1) = warr(1) - Jy/darr(1) 
        warr(2) = warr(2) + Jy/darr(2) 
        uarr(5) = uarr(5) + Jy/darr(5)
        uarr(6) = uarr(6) - Jy/darr(6)

        ! From this point, uarr,varr,warr is just the turbulent velocity
        trgv  = 0.0

        !x-direc
        ul    = (darr(2)*uarr(2) + d*uarr(0))/(darr(2)+d)
        ur    = (darr(1)*uarr(1) + d*uarr(0))/(darr(1)+d)
        trgv  = trgv + (ur-ul)**2
        !y-direc
        ul    = (darr(4)*varr(4) + d*varr(0))/(darr(4)+d)
        ur    = (darr(3)*varr(3) + d*varr(0))/(darr(3)+d)
        trgv  = trgv + (ur-ul)**2
        !z-direc
        ul    = (darr(6)*warr(6) + d*warr(0))/(darr(6)+d)
        ur    = (darr(5)*warr(5) + d*warr(0))/(darr(5)+d)
        trgv  = trgv + (ur-ul)**2
        !z-direc; tangential component - y
        ul    = (darr(6)*varr(6) + d*varr(0))/(darr(6)+d)
        ur    = (darr(5)*varr(5) + d*varr(0))/(darr(5)+d)
        trgv  = trgv + (ur-ul)**2
        !y-direc; tangential component - z
        ul    = (darr(4)*warr(4) + d*warr(0))/(darr(4)+d)
        ur    = (darr(3)*warr(3) + d*warr(0))/(darr(3)+d)
        trgv  = trgv + (ur-ul)**2
        !z-direc; tangential component - x
        ul    = (darr(6)*uarr(6) + d*uarr(0))/(darr(6)+d)
        ur    = (darr(5)*uarr(5) + d*uarr(0))/(darr(5)+d)
        trgv  = trgv + (ur-ul)**2
        !x-direc; tangential component - z
        ul    = (darr(2)*warr(2) + d*warr(0))/(darr(2)+d)
        ur    = (darr(1)*warr(1) + d*warr(0))/(darr(1)+d)
        trgv  = trgv + (ur-ul)**2
        !y-direc; tangential component - x
        ul    = (darr(4)*uarr(4) + d*uarr(0))/(darr(4)+d)
        ur    = (darr(3)*uarr(3) + d*uarr(0))/(darr(3)+d)
        trgv  = trgv + (ur-ul)**2
        !x-direc; tangential component - y
        ul    = (darr(2)*varr(2) + d*varr(0))/(darr(2)+d)
        ur    = (darr(1)*varr(1) + d*varr(0))/(darr(1)+d)
        trgv  = trgv + (ur-ul)**2
        
        ! we want the square root of this
        sigma = sqrt(trgv)
    end subroutine
    !---------------------------------------------------------------
    ! Function: SF EFF
    !
    ! Magneto-thermo-turbulent models of star formation commonly used
    ! in RAMSES simulations have non-constant star formation
    ! efficiency, so this routines obtains the actual value in the
    ! same way it is obtained in the RAMSES version from
    ! S. Martin-Alvarez
    !----------------------------------------------------------------
    function sf_eff(reg,dx,x,var,star_maker,use_crs)
        use vectors
        use basis_representations
        use geometrical_regions
        use coordinate_systems
        implicit none
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:amr%twondim,1:varIDs%nvar),intent(in) :: var
        character(128),intent(in)                 :: star_maker
        logical,intent(in)  :: use_crs
        real(dbl) :: sf_eff
        real(dbl) :: Tmin
        ! SF variables
        logical :: isConvergent
        real(dbl) :: n_gmc,n_dc,nCOM,d_gmc,d_dc
        real(dbl) :: nISM,n_star
        real(dbl) :: t0,t_star,del_star
        real(dbl) :: d,d0
        real(dbl) :: factG
        real(dbl),dimension(0:amr%twondim) ::darr,uarr,varr,warr
        real(dbl) :: divv,uavg,vavg,wavg,dtot
        real(dbl) :: px_div,py_div,pz_div,Jx,Jy,Jz,rho_local
        real(dbl) :: trgv,ul,ur
        real(dbl) :: temp,c_s2,temperature
        real(dbl) :: Bx,By,Bz,Pmag,invbeta
        real(dbl) :: lamjt,sf_lam
        real(dbl) :: e_cts,phi_t,theta,alpha0
        real(dbl) :: betafunc,scrit,sigs

        ! TODO: These parameters should NOT be hardcoded!
        n_gmc = 10D0
        n_dc = 1D10
        t_star = 0.632456 ! Star formation time scale (in Gyr) at the density threshold. Def: 0.632456
        del_star=200.     ! EOS density scale
        sf_lam = 1D0
        Tmin = 15d0/sim%T2

        ! 1. VARIABLE INITIALISATION
        sf_eff = 0D0
        ! Star formation time sclae from Gyr to code units
        t0    = t_star*(1D0/s2Gyr)/sim%unit_t
        ! ISM density threshold from H/cc to code units
        nISM = n_star
        nCOM  = del_star*sim%omega_b*rhoc*(sim%h0/100.)**2/sim%aexp**3*XH/mHydrogen
        nISM  = max(nCOM,nISM)
        d0    = nISM/sim%nH
        factG = 1d0
        if (sim%cosmo) factG = 3d0/4d0/twopi*sim%omega_m*sim%aexp
        d_gmc = max(nCOM,n_gmc)/sim%nH
        d_dc  = max(nCOM,n_dc)/sim%nH    ! dark cloud

        ! 2. COMPUTE SF MODEL
        if (TRIM(star_maker)=='FKmag') then
            ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
            d = var(0,varIDs%density)
            if (d<d_gmc) return
            ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
            ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor. 
            ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation 
            ! from neighbouring cell values and differentiate. 
            ! Get neighbor cells if they exist, otherwise use straight injection from local cell

            darr = var(:,varIDs%density)

            uarr = var(:,varIDs%vx)
            varr = var(:,varIDs%vy)
            warr = var(:,varIDs%vz)
            divv  = (uarr(2)*darr(2)-uarr(1)*darr(1)) &
                & + (varr(4)*darr(4)-varr(3)*darr(3)) & 
                & + (warr(6)*darr(6)-warr(5)*darr(5))
            if (divv>0) return
            
            ! Average velocity
            dtot  = sum(darr)
            uavg  = sum(darr*uarr)/dtot
            vavg  = sum(darr*varr)/dtot
            wavg  = sum(darr*warr)/dtot
            ! Subtract the mean velocity field
            uarr(:) = uarr(:) - uavg
            varr(:) = varr(:) - vavg
            warr(:) = warr(:) - wavg
            ! Subtract the symmetric divergence field                    
            ! ex)  (---->,<--): only subtract (-->,<--): result (-->,0) 
            ! ex)  (<----,-->): only subtract (<--,-->): result (<--,0)
            px_div = min( abs(darr(1)*uarr(1)),abs(darr(2)*uarr(2)))
            py_div = min( abs(darr(3)*varr(3)),abs(darr(4)*varr(4)))
            pz_div = min( abs(darr(5)*warr(5)),abs(darr(6)*warr(6)))

            isConvergent = darr(2)*uarr(2) - darr(1)*uarr(1) < 0 
            if (isConvergent) then
                uarr(1) = uarr(1) - px_div/darr(1)
                uarr(2) = uarr(2) + px_div/darr(2)
            else ! comment out if you do not want to subtract outflows
                uarr(1) = uarr(1) + px_div/darr(1)
                uarr(2) = uarr(2) - px_div/darr(2)
            end if 

            isConvergent = darr(4)*varr(4) - darr(3)*varr(3) < 0
            if (isConvergent) then
                varr(3) = varr(3) - py_div/darr(3)
                varr(4) = varr(4) + py_div/darr(4)
            else ! comment out if you do not want to subtract outflows
                varr(3) = varr(3) + py_div/darr(3)
                varr(4) = varr(4) - py_div/darr(4)
            end if

            isConvergent = darr(6)*warr(6) - darr(5)*warr(5) < 0
            if (isConvergent) then 
                warr(5) = warr(5) - pz_div/darr(5)
                warr(6) = warr(6) + pz_div/darr(6)
            else ! comment out if you do not want to subtract outflows
                warr(5) = warr(5) + pz_div/darr(5)
                warr(6) = warr(6) - pz_div/darr(6)
            end if

            ! subtract the rotational velocity field (x-y) plane
            ! ^y       <-        |4|        |-u|
            ! |       |  |     |1| |2|   |-v|  |+v|
            ! --->x    ->        |3|        |+u|
            Jz  = - varr(1)*darr(1) + varr(2)*darr(2) &
                &   + uarr(3)*darr(3) - uarr(4)*darr(4)
            Jz  = Jz / 4.0

            varr(1) = varr(1) + Jz/darr(1) 
            varr(2) = varr(2) - Jz/darr(2) 
            uarr(3) = uarr(3) - Jz/darr(3)
            uarr(4) = uarr(4) + Jz/darr(4)

            ! subtract the rotational velocity field (y-z) plane
            ! ^z       <-        |6|        |-v|  
            ! |       |  |     |3| |4|   |-w|  |+w|
            ! --->y    ->        |5|        |+v|
            Jx  = - warr(3)*darr(3) + warr(4)*darr(4) &
                &   + varr(5)*darr(5) - varr(6)*darr(6)
            Jx  = Jx / 4.0

            warr(3) = warr(3) + Jx/darr(3) 
            warr(4) = warr(4) - Jx/darr(4) 
            varr(5) = varr(5) - Jx/darr(5)
            varr(6) = varr(6) + Jx/darr(6)

            ! subtract the rotational velocity field (x-z) plane
            ! ^z       ->        |6|        |+u|  
            ! |       |  |     |1| |2|   |+w|  |-w|
            ! --->x    <-        |5|        |-u|
            Jy  = + warr(1)*darr(1) - warr(2)*darr(2) &
                &   - uarr(5)*darr(5) + uarr(6)*darr(6)
            Jy  = Jy / 4.0

            warr(1) = warr(1) - Jy/darr(1) 
            warr(2) = warr(2) + Jy/darr(2) 
            uarr(5) = uarr(5) + Jy/darr(5)
            uarr(6) = uarr(6) - Jy/darr(6)

            ! From this point, uarr,varr,warr is just the turbulent velocity
            trgv  = 0.0

            !x-direc
            ul    = (darr(2)*uarr(2) + d*uarr(0))/(darr(2)+d)
            ur    = (darr(1)*uarr(1) + d*uarr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc
            ul    = (darr(4)*varr(4) + d*varr(0))/(darr(4)+d)
            ur    = (darr(3)*varr(3) + d*varr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc
            ul    = (darr(6)*warr(6) + d*warr(0))/(darr(6)+d)
            ur    = (darr(5)*warr(5) + d*warr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc; tangential component - y
            ul    = (darr(6)*varr(6) + d*varr(0))/(darr(6)+d)
            ur    = (darr(5)*varr(5) + d*varr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc; tangential component - z
            ul    = (darr(4)*warr(4) + d*warr(0))/(darr(4)+d)
            ur    = (darr(3)*warr(3) + d*warr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc; tangential component - x
            ul    = (darr(6)*uarr(6) + d*uarr(0))/(darr(6)+d)
            ur    = (darr(5)*uarr(5) + d*uarr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !x-direc; tangential component - z
            ul    = (darr(2)*warr(2) + d*warr(0))/(darr(2)+d)
            ur    = (darr(1)*warr(1) + d*warr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc; tangential component - x
            ul    = (darr(4)*uarr(4) + d*uarr(0))/(darr(4)+d)
            ur    = (darr(3)*uarr(3) + d*uarr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !x-direc; tangential component - y
            ul    = (darr(2)*varr(2) + d*varr(0))/(darr(2)+d)
            ur    = (darr(1)*varr(1) + d*varr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2

            ! Now compute sound speed squared
            temp = var(0,varIDs%thermal_pressure) / var(0,varIDs%density)
            ! TODO: This is a quick fix for negative T in CRMHD sims
            if (temp < Tmin .and. fix_neg_temp) temp = Tmin
            temp = temp * gamma_gas
            ! TODO: This should be change to add also radiation pressure
            if (sim%cr.and.use_crs) then
                temp = temp + var(0,varIDs%cr_pressure) / var(0,varIDs%density) * gamma_crs
            end if
            temp = max(temp,smallc**2)
            c_s2 = temp
            if (sim%mhd) then
                ! Added magnetic pressure to the support as contribution to c_s (assuming isothermal gas within the cell)
                Bx = var(0,varIDs%Blx)+var(0,varIDs%Brx)
                By = var(0,varIDs%Bly)+var(0,varIDs%Bry)
                Bz = var(0,varIDs%Blz)+var(0,varIDs%Brz)
                ! 1/2 factor in Pmag comes from (Br+Bl)/2
                ! Pmag = 0.03978873577*0.25*(Bx**2 + By**2 + Bz**2) ! Pmag = B^2/(8*Pi)
                Pmag = 0.125d0*(Bx**2 + By**2 + Bz**2) ! Pmag = B^2/2
                invbeta = Pmag/(c_s2*d) ! beta = c_s^2 * rho / Pmag
                c_s2 = c_s2*(1.+invbeta) ! (c_s)_eff = c_s*(1.+beta^-1)**0.5
            end if
            ! Calculate "turbulent" Jeans length in cell units, lamjt 
            ! (see e.g. Chandrasekhar 51, Bonazzola et al 87, Federrath & Klessen 2012 eq 36)
            rho_local = d
            lamjt = (pi*trgv + sqrt(pi*pi*trgv*trgv + 36.0*pi*c_s2*factG*rho_local*dx**2))/(6.0*factG*rho_local*dx**2)
            if (lamjt > sf_lam) return ! Jeans length resolved: gas is stable
            ! print*,'Jeans length notk resolved'
            ! Jeans length not resolved --> form stars to lower density and stabilise gas
            ! corresponding virial parameter for homogeneous sphere <= 1.5 in turbulent dominated 
            ! limit and <= 0.5 in pressure dominated limit (in good agreement with observations,
            ! at least for massive (>= 10^5 M_sun) clouds see Kauffmann, Pillai, Goldsmith 2013, Fig.1)
            alpha0 = 5d0/(pi*factG*rho_local)*(trgv + c_s2)/dx**2
            ! Compute star formation efficiency per free-fall time (Federrath & Klessen 2012 eq 41)
            if (sim%mhd) then
                ! e_cts is the unresolved (for us) proto-stellar feedback: i.e. the dense gas core-to-star efficiency
                e_cts = 0.5  ! would be 1.0 without feedback (Federrath & Klessen 2012)
                phi_t = 0.57 ; theta = 0.33 ! best fit values of the Padoan & Nordlund multi-scale sf model to GMC simulation data
                betafunc = ((1.+0.925*invbeta**1.5)**(2./3.)) / ((1.+invbeta)**2) ! Magnetic field reduction of the critical log density
                scrit = log(0.067*betafunc/(theta**2)*alpha0*trgv/c_s2) ! best fit from Padoan & Nordlund MS model again
            else
                ! e_cts is the unresolved (for us) proto-stellar feedback: i.e. the dense gas core-to-star efficiency
                e_cts = 0.5  ! would be 1.0 without feedback (Federrath & Klessen 2012)
                phi_t = 0.57 ; theta = 0.33 ! best fit values of the Padoan & Nordlund multi-scale sf model to GMC simulation data 
                scrit = log(0.067/theta**2*alpha0*trgv/c_s2) ! best fit from Padoan & Nordlund MS model again
            end if
            sigs  = log(1.0+0.16*trgv/c_s2) ! factor 0.16 is b^2 where b=0.4 for a mixture of turbulence forcing modes
            
            sf_eff = e_cts/2.0*phi_t*exp(3.0/8.0*sigs)*(2.0-erfc((sigs-scrit)/sqrt(2.0*sigs)))
        else if (TRIM(star_maker)=='FK2') then
            ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
            d = var(0,varIDs%density)
            if (d<d_gmc) return
            ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
            ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor. 
            ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation 
            ! from neighbouring cell values and differentiate. 
            ! Get neighbor cells if they exist, otherwise use straight injection from local cell

            darr = var(:,varIDs%density)

            uarr = var(:,varIDs%vx)
            varr = var(:,varIDs%vy)
            warr = var(:,varIDs%vz)
            divv  = (uarr(2)*darr(2)-uarr(1)*darr(1)) &
                & + (varr(4)*darr(4)-varr(3)*darr(3)) & 
                & + (warr(6)*darr(6)-warr(5)*darr(5))
            if (divv>0) return
            
            ! Average velocity
            dtot  = sum(darr)
            uavg  = sum(darr*uarr)/dtot
            vavg  = sum(darr*varr)/dtot
            wavg  = sum(darr*warr)/dtot
            ! Subtract the mean velocity field
            uarr(:) = uarr(:) - uavg
            varr(:) = varr(:) - vavg
            warr(:) = warr(:) - wavg
            ! Subtract the symmetric divergence field                    
            ! ex)  (---->,<--): only subtract (-->,<--): result (-->,0) 
            ! ex)  (<----,-->): only subtract (<--,-->): result (<--,0)
            px_div = min( abs(darr(1)*uarr(1)),abs(darr(2)*uarr(2)))
            py_div = min( abs(darr(3)*varr(3)),abs(darr(4)*varr(4)))
            pz_div = min( abs(darr(5)*warr(5)),abs(darr(6)*warr(6)))

            isConvergent = darr(2)*uarr(2) - darr(1)*uarr(1) < 0 
            if (isConvergent) then
                uarr(1) = uarr(1) - px_div/darr(1)
                uarr(2) = uarr(2) + px_div/darr(2)
            else ! comment out if you do not want to subtract outflows
                uarr(1) = uarr(1) + px_div/darr(1)
                uarr(2) = uarr(2) - px_div/darr(2)
            end if 

            isConvergent = darr(4)*varr(4) - darr(3)*varr(3) < 0
            if (isConvergent) then
                varr(3) = varr(3) - py_div/darr(3)
                varr(4) = varr(4) + py_div/darr(4)
            else ! comment out if you do not want to subtract outflows
                varr(3) = varr(3) + py_div/darr(3)
                varr(4) = varr(4) - py_div/darr(4)
            end if

            isConvergent = darr(6)*warr(6) - darr(5)*warr(5) < 0
            if (isConvergent) then 
                warr(5) = warr(5) - pz_div/darr(5)
                warr(6) = warr(6) + pz_div/darr(6)
            else ! comment out if you do not want to subtract outflows
                warr(5) = warr(5) + pz_div/darr(5)
                warr(6) = warr(6) - pz_div/darr(6)
            end if

            ! subtract the rotational velocity field (x-y) plane
            ! ^y       <-        |4|        |-u|
            ! |       |  |     |1| |2|   |-v|  |+v|
            ! --->x    ->        |3|        |+u|
            Jz  = - varr(1)*darr(1) + varr(2)*darr(2) &
                &   + uarr(3)*darr(3) - uarr(4)*darr(4)
            Jz  = Jz / 4.0

            varr(1) = varr(1) + Jz/darr(1) 
            varr(2) = varr(2) - Jz/darr(2) 
            uarr(3) = uarr(3) - Jz/darr(3)
            uarr(4) = uarr(4) + Jz/darr(4)

            ! subtract the rotational velocity field (y-z) plane
            ! ^z       <-        |6|        |-v|  
            ! |       |  |     |3| |4|   |-w|  |+w|
            ! --->y    ->        |5|        |+v|
            Jx  = - warr(3)*darr(3) + warr(4)*darr(4) &
                &   + varr(5)*darr(5) - varr(6)*darr(6)
            Jx  = Jx / 4.0

            warr(3) = warr(3) + Jx/darr(3) 
            warr(4) = warr(4) - Jx/darr(4) 
            varr(5) = varr(5) - Jx/darr(5)
            varr(6) = varr(6) + Jx/darr(6)

            ! subtract the rotational velocity field (x-z) plane
            ! ^z       ->        |6|        |+u|  
            ! |       |  |     |1| |2|   |+w|  |-w|
            ! --->x    <-        |5|        |-u|
            Jy  = + warr(1)*darr(1) - warr(2)*darr(2) &
                &   - uarr(5)*darr(5) + uarr(6)*darr(6)
            Jy  = Jy / 4.0

            warr(1) = warr(1) - Jy/darr(1) 
            warr(2) = warr(2) + Jy/darr(2) 
            uarr(5) = uarr(5) + Jy/darr(5)
            uarr(6) = uarr(6) - Jy/darr(6)

            ! From this point, uarr,varr,warr is just the turbulent velocity
            trgv  = 0.0

            !x-direc
            ul    = (darr(2)*uarr(2) + d*uarr(0))/(darr(2)+d)
            ur    = (darr(1)*uarr(1) + d*uarr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc
            ul    = (darr(4)*varr(4) + d*varr(0))/(darr(4)+d)
            ur    = (darr(3)*varr(3) + d*varr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc
            ul    = (darr(6)*warr(6) + d*warr(0))/(darr(6)+d)
            ur    = (darr(5)*warr(5) + d*warr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc; tangential component - y
            ul    = (darr(6)*varr(6) + d*varr(0))/(darr(6)+d)
            ur    = (darr(5)*varr(5) + d*varr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc; tangential component - z
            ul    = (darr(4)*warr(4) + d*warr(0))/(darr(4)+d)
            ur    = (darr(3)*warr(3) + d*warr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc; tangential component - x
            ul    = (darr(6)*uarr(6) + d*uarr(0))/(darr(6)+d)
            ur    = (darr(5)*uarr(5) + d*uarr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !x-direc; tangential component - z
            ul    = (darr(2)*warr(2) + d*warr(0))/(darr(2)+d)
            ur    = (darr(1)*warr(1) + d*warr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc; tangential component - x
            ul    = (darr(4)*uarr(4) + d*uarr(0))/(darr(4)+d)
            ur    = (darr(3)*uarr(3) + d*uarr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !x-direc; tangential component - y
            ul    = (darr(2)*varr(2) + d*varr(0))/(darr(2)+d)
            ur    = (darr(1)*varr(1) + d*varr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2

            ! Now compute sound speed squared
            temp = var(0,varIDs%thermal_pressure) / var(0,varIDs%density)
            ! TODO: This is a quick fix for negative T in CRMHD sims
            if (temp < Tmin .and. fix_neg_temp) temp = Tmin
            temp = temp * gamma_gas
            temp = max(temp,smallc**2)
            c_s2 = temp
            ! Calculate "turbulent" Jeans length in cell units, lamjt 
            ! (see e.g. Chandrasekhar 51, Bonazzola et al 87, Federrath & Klessen 2012 eq 36)
            rho_local = d
            lamjt = (pi*trgv + sqrt(pi*pi*trgv*trgv + 36.0*pi*c_s2*factG*rho_local*dx**2))/(6.0*factG*rho_local*dx**2)
            if (lamjt > sf_lam) return ! Jeans length resolved: gas is stable
            ! print*,'Jeans length notk resolved'
            ! Jeans length not resolved --> form stars to lower density and stabilise gas
            ! corresponding virial parameter for homogeneous sphere <= 1.5 in turbulent dominated 
            ! limit and <= 0.5 in pressure dominated limit (in good agreement with observations,
            ! at least for massive (>= 10^5 M_sun) clouds see Kauffmann, Pillai, Goldsmith 2013, Fig.1)
            alpha0 = 5d0/(pi*factG*rho_local)*(trgv + c_s2)/dx**2
            ! Compute star formation efficiency per free-fall time (Federrath & Klessen 2012 eq 41)
            ! e_cts is the unresolved (for us) proto-stellar feedback: i.e. the dense gas core-to-star efficiency
            e_cts = 0.5  ! would be 1.0 without feedback (Federrath & Klessen 2012)
            phi_t = 0.57 ; theta = 0.33 ! best fit values of the Padoan & Nordlund multi-scale sf model to GMC simulation data 
            sigs  = log(1.0+0.16*trgv/c_s2) ! factor 0.16 is b^2 where b=0.4 for a mixture of turbulence forcing modes
            scrit = log(0.067/theta**2*alpha0*trgv/c_s2) ! best fit from Padoan & Nordlund MS model again 
            sf_eff = e_cts/2.0*phi_t*exp(3.0/8.0*sigs)*(2.0-erfc((sigs-scrit)/sqrt(2.0*sigs)))
        end if

        ! 3. COMPUTE JEANS LENGTH
        if (TRIM(star_maker)=='jeans_mtt') then
            ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
            ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor. 
            ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation 
            ! from neighbouring cell values and differentiate. 
            ! Get neighbor cells if they exist, otherwise use straight injection from local cell
            d = var(0,varIDs%density)
            darr = var(:,varIDs%density)

            uarr = var(:,varIDs%vx)
            varr = var(:,varIDs%vy)
            warr = var(:,varIDs%vz)
            divv  = (uarr(2)*darr(2)-uarr(1)*darr(1)) &
                & + (varr(4)*darr(4)-varr(3)*darr(3)) & 
                & + (warr(6)*darr(6)-warr(5)*darr(5))
            
            ! Average velocity
            dtot  = sum(darr)
            uavg  = sum(darr*uarr)/dtot
            vavg  = sum(darr*varr)/dtot
            wavg  = sum(darr*warr)/dtot
            ! Subtract the mean velocity field
            uarr(:) = uarr(:) - uavg
            varr(:) = varr(:) - vavg
            warr(:) = warr(:) - wavg
            ! Subtract the symmetric divergence field                    
            ! ex)  (---->,<--): only subtract (-->,<--): result (-->,0) 
            ! ex)  (<----,-->): only subtract (<--,-->): result (<--,0)
            px_div = min( abs(darr(1)*uarr(1)),abs(darr(2)*uarr(2)))
            py_div = min( abs(darr(3)*varr(3)),abs(darr(4)*varr(4)))
            pz_div = min( abs(darr(5)*warr(5)),abs(darr(6)*warr(6)))

            isConvergent = darr(2)*uarr(2) - darr(1)*uarr(1) < 0 
            if (isConvergent) then
                uarr(1) = uarr(1) - px_div/darr(1)
                uarr(2) = uarr(2) + px_div/darr(2)
            else ! comment out if you do not want to subtract outflows
                uarr(1) = uarr(1) + px_div/darr(1)
                uarr(2) = uarr(2) - px_div/darr(2)
            end if 

            isConvergent = darr(4)*varr(4) - darr(3)*varr(3) < 0
            if (isConvergent) then
                varr(3) = varr(3) - py_div/darr(3)
                varr(4) = varr(4) + py_div/darr(4)
            else ! comment out if you do not want to subtract outflows
                varr(3) = varr(3) + py_div/darr(3)
                varr(4) = varr(4) - py_div/darr(4)
            end if

            isConvergent = darr(6)*warr(6) - darr(5)*warr(5) < 0
            if (isConvergent) then 
                warr(5) = warr(5) - pz_div/darr(5)
                warr(6) = warr(6) + pz_div/darr(6)
            else ! comment out if you do not want to subtract outflows
                warr(5) = warr(5) + pz_div/darr(5)
                warr(6) = warr(6) - pz_div/darr(6)
            end if

            ! subtract the rotational velocity field (x-y) plane
            ! ^y       <-        |4|        |-u|
            ! |       |  |     |1| |2|   |-v|  |+v|
            ! --->x    ->        |3|        |+u|
            Jz  = - varr(1)*darr(1) + varr(2)*darr(2) &
                &   + uarr(3)*darr(3) - uarr(4)*darr(4)
            Jz  = Jz / 4.0

            varr(1) = varr(1) + Jz/darr(1) 
            varr(2) = varr(2) - Jz/darr(2) 
            uarr(3) = uarr(3) - Jz/darr(3)
            uarr(4) = uarr(4) + Jz/darr(4)

            ! subtract the rotational velocity field (y-z) plane
            ! ^z       <-        |6|        |-v|  
            ! |       |  |     |3| |4|   |-w|  |+w|
            ! --->y    ->        |5|        |+v|
            Jx  = - warr(3)*darr(3) + warr(4)*darr(4) &
                &   + varr(5)*darr(5) - varr(6)*darr(6)
            Jx  = Jx / 4.0

            warr(3) = warr(3) + Jx/darr(3) 
            warr(4) = warr(4) - Jx/darr(4) 
            varr(5) = varr(5) - Jx/darr(5)
            varr(6) = varr(6) + Jx/darr(6)

            ! subtract the rotational velocity field (x-z) plane
            ! ^z       ->        |6|        |+u|  
            ! |       |  |     |1| |2|   |+w|  |-w|
            ! --->x    <-        |5|        |-u|
            Jy  = + warr(1)*darr(1) - warr(2)*darr(2) &
                &   - uarr(5)*darr(5) + uarr(6)*darr(6)
            Jy  = Jy / 4.0

            warr(1) = warr(1) - Jy/darr(1) 
            warr(2) = warr(2) + Jy/darr(2) 
            uarr(5) = uarr(5) + Jy/darr(5)
            uarr(6) = uarr(6) - Jy/darr(6)

            ! From this point, uarr,varr,warr is just the turbulent velocity
            trgv  = 0.0

            !x-direc
            ul    = (darr(2)*uarr(2) + d*uarr(0))/(darr(2)+d)
            ur    = (darr(1)*uarr(1) + d*uarr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc
            ul    = (darr(4)*varr(4) + d*varr(0))/(darr(4)+d)
            ur    = (darr(3)*varr(3) + d*varr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc
            ul    = (darr(6)*warr(6) + d*warr(0))/(darr(6)+d)
            ur    = (darr(5)*warr(5) + d*warr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc; tangential component - y
            ul    = (darr(6)*varr(6) + d*varr(0))/(darr(6)+d)
            ur    = (darr(5)*varr(5) + d*varr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc; tangential component - z
            ul    = (darr(4)*warr(4) + d*warr(0))/(darr(4)+d)
            ur    = (darr(3)*warr(3) + d*warr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc; tangential component - x
            ul    = (darr(6)*uarr(6) + d*uarr(0))/(darr(6)+d)
            ur    = (darr(5)*uarr(5) + d*uarr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !x-direc; tangential component - z
            ul    = (darr(2)*warr(2) + d*warr(0))/(darr(2)+d)
            ur    = (darr(1)*warr(1) + d*warr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc; tangential component - x
            ul    = (darr(4)*uarr(4) + d*uarr(0))/(darr(4)+d)
            ur    = (darr(3)*uarr(3) + d*uarr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !x-direc; tangential component - y
            ul    = (darr(2)*varr(2) + d*varr(0))/(darr(2)+d)
            ur    = (darr(1)*varr(1) + d*varr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2

            ! Now compute sound speed squared
            temp = var(0,varIDs%thermal_pressure) / var(0,varIDs%density)
            ! TODO: This is a quick fix for negative T in CRMHD sims
            if (temp < Tmin .and. fix_neg_temp) temp = Tmin
            temp = temp * gamma_gas
            ! TODO: This should be change to add also radiation pressure
            if (sim%cr.and.use_crs) then
                temp = temp + var(0,varIDs%cr_pressure) / var(0,varIDs%density) * gamma_crs
            end if
            temp = max(temp,smallc**2)
            c_s2 = temp
            if (sim%mhd) then
                ! Added magnetic pressure to the support as contribution to c_s (assuming isothermal gas within the cell)
                Bx = var(0,varIDs%Blx)+var(0,varIDs%Brx)
                By = var(0,varIDs%Bly)+var(0,varIDs%Bry)
                Bz = var(0,varIDs%Blz)+var(0,varIDs%Brz)
                ! 1/2 factor in Pmag comes from (Br+Bl)/2
                ! Pmag = 0.03978873577*0.25*(Bx**2 + By**2 + Bz**2) ! Pmag = B^2/(8*Pi)
                Pmag = 0.125d0*(Bx**2 + By**2 + Bz**2) ! Pmag = B^2/2
                invbeta = Pmag/(c_s2*d) ! beta = c_s^2 * rho / Pmag
                c_s2 = c_s2*(1.+invbeta) ! (c_s)_eff = c_s*(1.+beta^-1)**0.5
            end if
            ! Calculate "turbulent" Jeans length in code length
            ! (see e.g. Chandrasekhar 51, Bonazzola et al 87, Federrath & Klessen 2012 eq 36)
            rho_local = d
            lamjt = (pi*trgv + sqrt(pi*pi*trgv*trgv + 36.0*pi*c_s2*factG*rho_local*dx**2))/(6.0*factG*rho_local*dx)
            sf_eff = lamjt ! TODO: Don't like this change of variable names...
        elseif (TRIM(star_maker)=='jeansmtt_dx') then
            ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
            ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor. 
            ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation 
            ! from neighbouring cell values and differentiate. 
            ! Get neighbor cells if they exist, otherwise use straight injection from local cell

            d = var(0,varIDs%density)
            darr = var(:,varIDs%density)

            uarr = var(:,varIDs%vx)
            varr = var(:,varIDs%vy)
            warr = var(:,varIDs%vz)
            divv  = (uarr(2)*darr(2)-uarr(1)*darr(1)) &
                & + (varr(4)*darr(4)-varr(3)*darr(3)) & 
                & + (warr(6)*darr(6)-warr(5)*darr(5))
            
            ! Average velocity
            dtot  = sum(darr)
            uavg  = sum(darr*uarr)/dtot
            vavg  = sum(darr*varr)/dtot
            wavg  = sum(darr*warr)/dtot
            ! Subtract the mean velocity field
            uarr(:) = uarr(:) - uavg
            varr(:) = varr(:) - vavg
            warr(:) = warr(:) - wavg
            ! Subtract the symmetric divergence field                    
            ! ex)  (---->,<--): only subtract (-->,<--): result (-->,0) 
            ! ex)  (<----,-->): only subtract (<--,-->): result (<--,0)
            px_div = min( abs(darr(1)*uarr(1)),abs(darr(2)*uarr(2)))
            py_div = min( abs(darr(3)*varr(3)),abs(darr(4)*varr(4)))
            pz_div = min( abs(darr(5)*warr(5)),abs(darr(6)*warr(6)))

            isConvergent = darr(2)*uarr(2) - darr(1)*uarr(1) < 0 
            if (isConvergent) then
                uarr(1) = uarr(1) - px_div/darr(1)
                uarr(2) = uarr(2) + px_div/darr(2)
            else ! comment out if you do not want to subtract outflows
                uarr(1) = uarr(1) + px_div/darr(1)
                uarr(2) = uarr(2) - px_div/darr(2)
            end if 

            isConvergent = darr(4)*varr(4) - darr(3)*varr(3) < 0
            if (isConvergent) then
                varr(3) = varr(3) - py_div/darr(3)
                varr(4) = varr(4) + py_div/darr(4)
            else ! comment out if you do not want to subtract outflows
                varr(3) = varr(3) + py_div/darr(3)
                varr(4) = varr(4) - py_div/darr(4)
            end if

            isConvergent = darr(6)*warr(6) - darr(5)*warr(5) < 0
            if (isConvergent) then 
                warr(5) = warr(5) - pz_div/darr(5)
                warr(6) = warr(6) + pz_div/darr(6)
            else ! comment out if you do not want to subtract outflows
                warr(5) = warr(5) + pz_div/darr(5)
                warr(6) = warr(6) - pz_div/darr(6)
            end if

            ! subtract the rotational velocity field (x-y) plane
            ! ^y       <-        |4|        |-u|
            ! |       |  |     |1| |2|   |-v|  |+v|
            ! --->x    ->        |3|        |+u|
            Jz  = - varr(1)*darr(1) + varr(2)*darr(2) &
                &   + uarr(3)*darr(3) - uarr(4)*darr(4)
            Jz  = Jz / 4.0

            varr(1) = varr(1) + Jz/darr(1) 
            varr(2) = varr(2) - Jz/darr(2) 
            uarr(3) = uarr(3) - Jz/darr(3)
            uarr(4) = uarr(4) + Jz/darr(4)

            ! subtract the rotational velocity field (y-z) plane
            ! ^z       <-        |6|        |-v|  
            ! |       |  |     |3| |4|   |-w|  |+w|
            ! --->y    ->        |5|        |+v|
            Jx  = - warr(3)*darr(3) + warr(4)*darr(4) &
                &   + varr(5)*darr(5) - varr(6)*darr(6)
            Jx  = Jx / 4.0

            warr(3) = warr(3) + Jx/darr(3) 
            warr(4) = warr(4) - Jx/darr(4) 
            varr(5) = varr(5) - Jx/darr(5)
            varr(6) = varr(6) + Jx/darr(6)

            ! subtract the rotational velocity field (x-z) plane
            ! ^z       ->        |6|        |+u|  
            ! |       |  |     |1| |2|   |+w|  |-w|
            ! --->x    <-        |5|        |-u|
            Jy  = + warr(1)*darr(1) - warr(2)*darr(2) &
                &   - uarr(5)*darr(5) + uarr(6)*darr(6)
            Jy  = Jy / 4.0

            warr(1) = warr(1) - Jy/darr(1) 
            warr(2) = warr(2) + Jy/darr(2) 
            uarr(5) = uarr(5) + Jy/darr(5)
            uarr(6) = uarr(6) - Jy/darr(6)

            ! From this point, uarr,varr,warr is just the turbulent velocity
            trgv  = 0.0

            !x-direc
            ul    = (darr(2)*uarr(2) + d*uarr(0))/(darr(2)+d)
            ur    = (darr(1)*uarr(1) + d*uarr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc
            ul    = (darr(4)*varr(4) + d*varr(0))/(darr(4)+d)
            ur    = (darr(3)*varr(3) + d*varr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc
            ul    = (darr(6)*warr(6) + d*warr(0))/(darr(6)+d)
            ur    = (darr(5)*warr(5) + d*warr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc; tangential component - y
            ul    = (darr(6)*varr(6) + d*varr(0))/(darr(6)+d)
            ur    = (darr(5)*varr(5) + d*varr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc; tangential component - z
            ul    = (darr(4)*warr(4) + d*warr(0))/(darr(4)+d)
            ur    = (darr(3)*warr(3) + d*warr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc; tangential component - x
            ul    = (darr(6)*uarr(6) + d*uarr(0))/(darr(6)+d)
            ur    = (darr(5)*uarr(5) + d*uarr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !x-direc; tangential component - z
            ul    = (darr(2)*warr(2) + d*warr(0))/(darr(2)+d)
            ur    = (darr(1)*warr(1) + d*warr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc; tangential component - x
            ul    = (darr(4)*uarr(4) + d*uarr(0))/(darr(4)+d)
            ur    = (darr(3)*uarr(3) + d*uarr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !x-direc; tangential component - y
            ul    = (darr(2)*varr(2) + d*varr(0))/(darr(2)+d)
            ur    = (darr(1)*varr(1) + d*varr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2

            ! Now compute sound speed squared
            temp = var(0,varIDs%thermal_pressure) / var(0,varIDs%density)
            ! TODO: This is a quick fix for negative T in CRMHD sims
            if (temp < Tmin .and. fix_neg_temp) temp = Tmin
            temp = temp * gamma_gas
            ! TODO: This should be change to add also radiation pressure
            if (sim%cr.and.use_crs) then
                temp = temp + var(0,varIDs%cr_pressure) / var(0,varIDs%density) * gamma_crs
            end if
            temp = max(temp,smallc**2)
            c_s2 = temp
            if (sim%mhd) then
                ! Added magnetic pressure to the support as contribution to c_s (assuming isothermal gas within the cell)
                Bx = var(0,varIDs%Blx)+var(0,varIDs%Brx)
                By = var(0,varIDs%Bly)+var(0,varIDs%Bry)
                Bz = var(0,varIDs%Blz)+var(0,varIDs%Brz)
                ! 1/2 factor in Pmag comes from (Br+Bl)/2
                ! Pmag = 0.03978873577*0.25*(Bx**2 + By**2 + Bz**2) ! Pmag = B^2/(8*Pi)
                Pmag = 0.125d0*(Bx**2 + By**2 + Bz**2) ! Pmag = B^2/2
                invbeta = Pmag/(c_s2*d) ! beta = c_s^2 * rho / Pmag
                c_s2 = c_s2*(1.+invbeta) ! (c_s)_eff = c_s*(1.+beta^-1)**0.5
            end if
            ! Calculate "turbulent" Jeans length in cell units
            ! (see e.g. Chandrasekhar 51, Bonazzola et al 87, Federrath & Klessen 2012 eq 36)
            rho_local = d
            lamjt = (pi*trgv + sqrt(pi*pi*trgv*trgv + 36.0*pi*c_s2*factG*rho_local*dx**2))/(6.0*factG*rho_local*dx**2d0)
            sf_eff = lamjt ! TODO: Don't like this change of variable names...
        end if
        
    end function sf_eff

    !---------------------------------------------------------------
    ! Subroutine: INITIAL AMR SETUP
    !
    ! This routine has been adapted from the RAMSES utils, in which
    ! the initial read of the amr data is performed and saved to 
    ! amr_info and sim_info derived types.
    !---------------------------------------------------------------
    subroutine init_amr_read(repository)
        implicit none
        character(128), intent(in) :: repository
        character(5) :: nchar
        character(13) :: namestr13
        character(14) :: namestr14
        integer :: ipos,impi,i,nx,ny,nz
        character(128) :: nomfich
        character(128)  :: cooling_file
        logical :: ok

        ipos = index(repository,'output_')
        nchar=repository(ipos+7:ipos+13)
        nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
        ! Verify input hydro file
        inquire(file=nomfich, exist=ok)
        if ( .not. ok ) then
            print *,TRIM(nomfich)//' not found.'
            if (verbose) write(*,*)': No hydro data in this simulation.'
        else
            sim%hydro = .true.
        endif
        nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
        ! Verify input part file
        inquire(file=nomfich, exist=ok)
        if ( .not. ok ) then
            print *,TRIM(nomfich)//' not found.'
            if (verbose) write(*,*)': No particles in this simulation.'
        else
            sim%dm = .true.
        endif
        if (sim%dm .and. .not. sim%hydro .and. verbose) write(*,*)': This is a DM only simulation.'
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
        read(10)amr%ngridmax
        read(10)amr%nboundary
        read(10)!ngrid_current
        read(10)!boxlen
        close(10)
        amr%xbound = (/dble(nx/2),dble(ny/2),dble(nz/2)/)
        amr%twotondim = 2**amr%ndim
        amr%twondim = 2*amr%ndim
        amr%ncoarse = nx*ny*nz

        if(amr%ndim==2)then
            if (verbose) write(*,*)'Output file contains 2D data'
            if (verbose) write(*,*)'Aborting!'
            stop
        endif

        nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
        open(unit=10,file=nomfich,form='formatted',status='old')
        read(10,*)!ncpu
        read(10,*)!ndim
        read(10,'(A13,I11)')namestr13,amr%levelmin
        read(10,'(A13,I11)')namestr13,amr%levelmax
        read(10,*)
        read(10,*)
        read(10,*)
    
        read(10,'(A13,E23.15)')namestr13,sim%boxlen
        read(10,'(A13,E23.15)')namestr13,sim%t
        read(10,'(A13,E23.15)')namestr13,sim%aexp
        read(10,'(A13,E23.15)')namestr13,sim%h0
        read(10,'(A13,E23.15)')namestr13,sim%omega_m
        read(10,'(A13,E23.15)')namestr13,sim%omega_l
        read(10,'(A13,E23.15)')namestr13,sim%omega_k
        read(10,'(A13,E23.15)')namestr13,sim%omega_b
        read(10,'(A13,E23.15)')namestr13,sim%unit_l
        read(10,'(A13,E23.15)')namestr13,sim%unit_d
        read(10,'(A13,E23.15)')namestr13,sim%unit_t
        read(10,*)
        sim%unit_m = ((sim%unit_d*sim%unit_l)*sim%unit_l)*sim%unit_l
        read(10,'(A14,A80)')namestr14,amr%ordering
        read(10,*)
        if(allocated(amr%cpu_list))deallocate(amr%cpu_list)
        allocate(amr%cpu_list(1:amr%ncpu))
        if(TRIM(amr%ordering).eq.'hilbert')then
            if(allocated(amr%bound_key))deallocate(amr%bound_key,amr%cpu_read)
            allocate(amr%bound_key(0:amr%ncpu))
            allocate(amr%cpu_read(1:amr%ncpu))
            amr%cpu_read=.false.
            do impi=1,amr%ncpu
                read(10,'(I8,1X,E23.15,1X,E23.15)')i,amr%bound_key(impi-1),amr%bound_key(impi)
            end do
        endif
        close(10)

#ifndef IMASS
#warning: I am compiling without IMASS
        call get_eta_sn(repository)
#endif
        sim%redshift = 1.0d0/sim%aexp - 1.0d0                          ! Current redshift
        sim%T2 = mHydrogen / kBoltzmann * ((sim%unit_l/sim%unit_t)**2) ! Temperature conversion factor
        sim%nH = XH / mHydrogen * sim%unit_d                           ! nH conversion factor
        sim%unit_v = (sim%unit_l/sim%unit_t)                           ! Velocity unit
        sim%unit_p = sim%unit_d*((sim%unit_l/sim%unit_t)**2)           ! Pressure unit

        if (sim%hydro) then
            ! Also initialise the cooling table
            cooling_file=TRIM(repository)//'/cooling_'//TRIM(nchar)//'.out'
            inquire(file=cooling_file, exist=ok)
            if(ok) call read_cool(cooling_file)
            call get_Dcr(repository)
        endif
    end subroutine init_amr_read

    !---------------------------------------------------------------
    ! Subroutine: COMPUTE CPU MAP
    !
    ! This has been adapted from the RAMSES utils, in which the cpu
    ! map is obtained for a desired region.
    !---------------------------------------------------------------
    subroutine get_cpu_map(reg)
        use geometrical_regions
        implicit none
        type(region),intent(in) :: reg
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

    subroutine getnborgrids(son,nbor,igrid,igridn,ngrid)
        implicit none
        integer::ngrid
        integer,dimension(1:amr%ncoarse+amr%twotondim*amr%ngridmax)::son
        integer,dimension(1:amr%ngridmax,1:amr%twondim)::nbor
        integer,dimension(1:ngrid)::igrid
        integer,dimension(1:ngrid,0:amr%twondim)::igridn
        !---------------------------------------------------------
        ! This routine computes the index of the 6 neighboring 
        ! grids for grid igrid(:). The index for the central 
        ! grid is stored in igridn(:,0). If for some reasons
        ! the neighboring grids don't exist, then igridn(:,j) = 0.
        !---------------------------------------------------------
        integer::i,j
    
        ! Store central grid
        do i=1,ngrid
            igridn(i,0)=igrid(i)
        end do
        ! Store neighboring grids
        do j=1,amr%twondim
            do i=1,ngrid
                if (son(nbor(igrid(i),j))>0) then
                    igridn(i,j)=son(nbor(igrid(i),j))
                else
                    igridn(i,j)=-nbor(igrid(i),j)
                end if
            end do
        end do
    end subroutine getnborgrids

    subroutine getnborcells(igridn,ind,icelln,ng)
        implicit none
        integer::ng,ind
        integer,dimension(1:ng,0:amr%twondim)::igridn
        integer,dimension(1:ng,1:amr%twondim)::icelln
        !--------------------------------------------------------------
        ! This routine computes the index of 6-neighboring cells
        ! The user must provide igridn = index of the 6 neighboring
        ! grids and the cell's grid (see routine getnborgrids). 
        ! ind is the cell index in the grid.
        !--------------------------------------------------------------
        integer::i,in,ig,ih,iskip
        integer,dimension(1:8,1:6)::ggg,hhh
        
        ggg(1:8,1)=(/1,0,1,0,1,0,1,0/); hhh(1:8,1)=(/2,1,4,3,6,5,8,7/)
        ggg(1:8,2)=(/0,2,0,2,0,2,0,2/); hhh(1:8,2)=(/2,1,4,3,6,5,8,7/)
        ggg(1:8,3)=(/3,3,0,0,3,3,0,0/); hhh(1:8,3)=(/3,4,1,2,7,8,5,6/)
        ggg(1:8,4)=(/0,0,4,4,0,0,4,4/); hhh(1:8,4)=(/3,4,1,2,7,8,5,6/)
        ggg(1:8,5)=(/5,5,5,5,0,0,0,0/); hhh(1:8,5)=(/5,6,7,8,1,2,3,4/)
        ggg(1:8,6)=(/0,0,0,0,6,6,6,6/); hhh(1:8,6)=(/5,6,7,8,1,2,3,4/)
        
        ! Reset indices
        icelln(1:ng,1:amr%twondim)=0
        ! Compute cell numbers
        do in=1,amr%twondim
            ig=ggg(ind,in)
            ih=hhh(ind,in)
            iskip=amr%ncoarse+(ih-1)*amr%ngridmax
            do i=1,ng
                if(igridn(i,ig)>0)then
                    icelln(i,in)=iskip+igridn(i,ig)
                else
                    icelln(i,in)=abs(igridn(i,ig))
                end if
            end do
        end do        
      
    end subroutine getnborcells

    subroutine getnbor(son,nbor,ind_cell,ind_father,ncell)
        implicit none
        integer::ncell
        integer,dimension(1:amr%ncoarse+amr%twotondim*amr%ngridmax)::son
        integer,dimension(1:amr%ngridmax,1:amr%twondim)::nbor
        integer,dimension(1:ncell)::ind_cell
        integer,dimension(1:ncell,0:amr%twondim)::ind_father
        !-----------------------------------------------------------------
        ! This subroutine determines the 2*ndim neighboring cells
        ! cells of the input cell (ind_cell). 
        ! If for some reasons they don't exist, the routine returns 
        ! the input cell.
        !-----------------------------------------------------------------
        integer::nxny,i,idim,j,iok,ind
        integer,dimension(1:3)::ibound,iskip1,iskip2
        integer,dimension(1:ncell,1:3)::ix
        integer,dimension(1:ncell)::ind_grid_father,pos
        integer,dimension(1:ncell,0:amr%twondim)::igridn,igridn_ok
        integer,dimension(1:ncell,1:amr%twondim)::icelln_ok

        ! if (verbose) write(*,*)'ncoarse,ngridmax,twondim: ',amr%ncoarse,amr%ngridmax,amr%twondim
        ! Get father cell
        do i=1,ncell
           ind_father(i,0)=ind_cell(i)
        end do
        
        ! Get father cell position in the grid
        do i=1,ncell
           pos(i)=(ind_father(i,0)-amr%ncoarse-1)/amr%ngridmax+1
        end do
        
        ! Get father grid
        do i=1,ncell
           ind_grid_father(i)=ind_father(i,0)-amr%ncoarse-(pos(i)-1)*amr%ngridmax
        end do
        
        ! Get neighboring father grids
        call getnborgrids(son,nbor,ind_grid_father,igridn,ncell)
        
        ! Loop over position
        do ind=1,amr%twotondim
           
           ! Select father cells that sit at position ind
           do j=0,amr%twondim
              iok=0
              do i=1,ncell
                 if(pos(i)==ind)then
                    iok=iok+1
                    igridn_ok(iok,j)=igridn(i,j)
                 end if
              end do
           end do
           
           ! Get neighboring cells for selected cells
           if(iok>0)call getnborcells(igridn_ok,ind,icelln_ok,iok)
           
           ! Update neighboring father cells for selected cells
           do j=1,amr%twondim
              iok=0
              do i=1,ncell
                 if(pos(i)==ind)then
                    iok=iok+1
                    if(icelln_ok(iok,j)>0)then
                       ind_father(i,j)=icelln_ok(iok,j)
                    else
                       ind_father(i,j)=ind_cell(i)
                    end if
                 end if
              end do
           end do
           
        end do
           
          
    end subroutine getnbor
    
    subroutine getparttype(part,ptype)
        implicit none
        type(particle),intent(in) :: part
        character(6),intent(inout) :: ptype

        if (part%age.ne.0D0) then
            ptype = 'star'
        else
            ptype = 'dm'
        endif
    end subroutine getparttype

    subroutine getpartvalue(reg,part,var,value,dx)
        use vectors
        use basis_representations
        use geometrical_regions
        use coordinate_systems
        implicit none
        type(region),intent(in) :: reg
        type(particle),intent(in) ::part
        character(128),intent(in) :: var
        real(dbl),intent(inout) :: value
        type(vector),optional,intent(in) :: dx
        
        type(vector) :: v,L
        type(basis) :: temp_basis
        character(6) :: ptype
        character(128) :: tempvar,vartype,varname,sfrstr,sfrtype
        integer :: index,iii,index2,index3
        real(dbl) :: time,age,birth_date,sfrind,current_age_univ

        tempvar = TRIM(var)
        index = scan(tempvar,'/')
        vartype = tempvar(1:index-1)
        varname = tempvar(index+1:)

        if (sim%dm .and. .not. sim%hydro) then
            ptype = 'dm'
        else
            call getparttype(part,ptype)
        endif

        if (vartype.eq.'dm'.and.ptype.eq.'dm') then
            select case (TRIM(varname))
            case ('mass')
                ! Mass
                value = part%m
            case ('density')
                ! Density
                if (present(dx)) then
                    value = part%m / (dx%x*dx%y+dx%z)
                else
                    write(*,*)'Can not compute a particle density without cell size!'
                    stop
                endif
            case ('sdensity')
                ! Surface density
                if (present(dx)) then
                    value = part%m / (dx%x*dx%y)
                else
                    write(*,*)'Can not compute a particle surface density without cell size!'
                    stop
                endif
            case ('x')
                ! x - coordinate
                value = part%x%x
            case ('y')
                ! y - coordinate
                value = part%x%y
            case ('z')
                ! z - coordinate
                value = part%x%z
            case ('vx')
                ! velocity x - component
                value = part%v%x
            case ('vy')
                ! velocity y - component
                value = part%v%y
            case ('vz')
                ! velocity z - component
                value = part%v%z
            case ('d_euclid')
                ! Euclidean distance
                value = magnitude(part%x)
            case ('r_sphere')
                ! Radius from center of sphere
                value = r_sphere(part%x)
            case ('theta_sphere')
                ! Value of spherical theta angle measured from the z axis
                value = theta_sphere(part%x)
            case ('phi_sphere')
                ! Value of spherical phi angle measure in the x-y plane 
                ! from the x axis
                value = phi_sphere(part%x)
            case('r_cyl')
                ! Value of cylindrical radius
                value = r_cyl(part%x)
            case ('phi_cyl')
                ! Value of spherical phi angle measure in the x-y plane 
                ! from the x axis
                value = phi_cyl(part%x)
            case ('v_tangential')
                ! Tangential velocity magnitude
                ! Consider as if one substracts the radial velocity,
                ! the remaining component is the tangential velocity
                v = part%v
                call spherical_basis_from_cartesian(part%x,temp_basis)
                v = v - (v.DOT.temp_basis%u(1))*temp_basis%u(1)
                value = magnitude(v)
            case ('centripetal_acc')
                ! Radial centripetal acceleration, obtained from the
                ! tangential velocity
                v = part%v
                call spherical_basis_from_cartesian(part%x,temp_basis)
                v = v - (v.DOT.temp_basis%u(1))*temp_basis%u(1)
                value = magnitude(v)**2d0 / r_sphere(part%x)
            case ('v_sphere_r')
                ! Velocity component in the spherical radial direction
                ! Dot product of velocity vector with spherical radial
                !    unit vector
                v = part%v
                call spherical_basis_from_cartesian(part%x,temp_basis)
                value = v.DOT.temp_basis%u(1)
            case ('v_sphere_phi')
                ! Velocity component in the spherical azimutal (phi) direction
                ! Dot product of velocity vector with spherical phi
                !    unit vector
                v = part%v
                call spherical_basis_from_cartesian(part%x,temp_basis)
                value = v .DOT. temp_basis%u(3)
            case ('v_sphere_theta')
                ! Velocity component in the spherical theta direction
                ! Dot product of velocity vector with spherical theta
                !    unit vector
                v = part%v
                call spherical_basis_from_cartesian(part%x,temp_basis)
                value = v .DOT. temp_basis%u(2)
            case ('v_cyl_r')
                ! Velocity component in the cylindrical radial direction
                ! Dot product of velocity vector with cylindrical
                !    radial unit vector
                v = part%v
                call cylindrical_basis_from_cartesian(part%x,temp_basis)
                value = v.DOT.temp_basis%u(1)
            case ('v_cyl_z')
                ! Velocity component in the cylindrical z direction
                ! Dot product of velocity vector with cylindrical
                !    z unit vector
                v = part%v
                call cylindrical_basis_from_cartesian(part%x,temp_basis)
                value = v.DOT.temp_basis%u(3)
            case ('v_cyl_phi')
                ! Velocity component in the cylyndrical azimutal (phi) direction
                ! Dot product of velocity vector with cylindrical
                !    phi unit vector
                v = part%v
                call cylindrical_basis_from_cartesian(part%x,temp_basis)
                value = v.DOT.temp_basis%u(2)
            case ('momentum_x')
                ! Linear momentum in the x direction as mass*corrected_velocity_x
                value = part%m * (part%v%x)
            case ('momentum_y')
                ! Linear momentum in the y direction mass*corrected_velocity_y
                value = part%m * (part%v%y)
            case ('momentum_z')
                ! Linear momentum in the z direction mass*corrected_velocity_z
                value = part%m * (part%v%z)
            case ('momentum')
                ! Magnitude of linear momentum, using corrected velocity
                v = part%v
                value = part%m * magnitude(v)
            case ('momentum_sphere_r')
                ! Linear momentum in the spherical radial direction
                ! Dot product of velocity vector with spherical phi
                !    unit vector
                ! 3. Multiply by mass of particle
                v = part%v
                call spherical_basis_from_cartesian(part%x,temp_basis)
                value = part%m * (v .DOT. temp_basis%u(1))
            case ('ang_momentum_x')
                ! Corrected angular momentum in the x direction
                v = part%v
                value = part%m * (part%x%y * v%z &
                            &- v%y * part%x%z)
            case ('ang_momentum_y')
                ! Corrected angular momentum in the y direction
                v = part%v
                value = part%m * (part%x%z * v%x &
                            &- v%z * part%x%x)
            case ('ang_momentum_z')
                ! Corrected angular momentum in the z direction
                v = part%v
                value = part%m * (part%x%x * v%y &
                            &- v%x * part%x%y)
            case ('ang_momentum')
                ! Corrected magnitude of angular momentum
                v = part%v
                L = part%x * v
                value = part%m * magnitude(L)
            case ('ang_momentum_specific_x')
                ! Corrected specific angular momentum in the x direction
                v = part%v
                value = part%x%y * v%z - v%y * part%x%z
            case ('ang_momentum_specific_y')
                ! Corrected specific angular momentum in the y direction
                v = part%v
                value = part%x%z * v%x - v%z * part%x%x
            case ('ang_momentum_specific_z')
                ! Corrected specific angular momentum in the z direction
                v = part%v
                value = part%x%x * v%y - v%x * part%x%y
            case ('ang_momentum_specific')
                ! Corrected magnitude of specific angular momentum
                v = part%v
                L = part%x * v
                value = magnitude(L)
            end select
        elseif (vartype.eq.'star'.and.ptype.eq.'star') then
            index2 = scan(varname,'_')
            sfrstr = varname(1:index2-1)
            if (trim(sfrstr).ne.'sfr') then
                select case (TRIM(varname))
                case ('mass')
                    ! Mass
                    value = part%m
                case ('density')
                    ! Density
                    if (present(dx)) then
                        value = part%m / (dx%x*dx%y+dx%z)
                    else
                        write(*,*)'Can not compute a particle density without cell size!'
                        stop
                    endif
                case ('sdensity')
                    ! Surface density
                    if (present(dx)) then
                        value = part%m / (dx%x*dx%y)
                    else
                        write(*,*)'Can not compute a particle surface density without cell size!'
                        stop
                    endif
                case ('x')
                    ! x - coordinate
                    value = part%x%x
                case ('y')
                    ! y - coordinate
                    value = part%x%y
                case ('z')
                    ! z - coordinate
                    value = part%x%z
                case ('vx')
                    ! velocity x - component
                    value = part%v%x
                case ('vy')
                    ! velocity y - component
                    value = part%v%y
                case ('vz')
                    ! velocity z - component
                    value = part%v%z
                case ('d_euclid')
                    ! Euclidean distance
                    value = magnitude(part%x)
                case ('r_sphere')
                    ! Radius from center of sphere
                    value = r_sphere(part%x)
                case ('theta_sphere')
                    ! Value of spherical theta angle measured from the z axis
                    value = theta_sphere(part%x)
                case ('phi_sphere')
                    ! Value of spherical phi angle measure in the x-y plane 
                    ! from the x axis
                    value = phi_sphere(part%x)
                case('r_cyl')
                    ! Value of cylindrical radius
                    value = r_cyl(part%x)
                case ('phi_cyl')
                    ! Value of spherical phi angle measure in the x-y plane 
                    ! from the x axis
                    value = phi_cyl(part%x)
                case ('v_sphere_r')
                    ! Velocity component in the spherical radial direction
    
                    ! Dot product of velocity vector with spherical radial
                    !    unit vector
                    v = part%v
                    call spherical_basis_from_cartesian(part%x,temp_basis)
                    value = v.DOT.temp_basis%u(1)
                case ('v_sphere_phi')
                    ! Velocity component in the spherical azimutal (phi) direction
                    ! Dot product of velocity vector with spherical phi
                    !    unit vector
                    v = part%v
                    call spherical_basis_from_cartesian(part%x,temp_basis)
                    value = v .DOT. temp_basis%u(3)
                case ('v_sphere_theta')
                    ! Velocity component in the spherical theta direction
                    ! Dot product of velocity vector with spherical theta
                    !    unit vector
                    v = part%v
                    call spherical_basis_from_cartesian(part%x,temp_basis)
                    value = v .DOT. temp_basis%u(2)
                case ('v_cyl_r')
                    ! Velocity component in the cylindrical radial direction
                    ! Dot product of velocity vector with cylindrical
                    !    radial unit vector
                    v = part%v
                    call cylindrical_basis_from_cartesian(part%x,temp_basis)
                    value = v.DOT.temp_basis%u(1)
                case ('v_cyl_z')
                    ! Velocity component in the cylindrical z direction
                    ! Dot product of velocity vector with cylindrical
                    !    z unit vector
                    v = part%v
                    call cylindrical_basis_from_cartesian(part%x,temp_basis)
                    value = v.DOT.temp_basis%u(3)
                case ('v_cyl_phi')
                    ! Velocity component in the cylyndrical azimutal (phi) direction
                    ! Dot product of velocity vector with cylindrical
                    !    phi unit vector
                    v = part%v
                    call cylindrical_basis_from_cartesian(part%x,temp_basis)
                    value = v.DOT.temp_basis%u(2)
                case ('momentum_x')
                    ! Linear momentum in the x direction as mass*corrected_velocity_x
                    value = part%m * part%v%x
                case ('momentum_y')
                    ! Linear momentum in the y direction mass*corrected_velocity_y
                    value = part%m * part%v%y
                case ('momentum_z')
                    ! Linear momentum in the z direction mass*corrected_velocity_z
                    value = part%m * part%v%z
                case ('momentum')
                    ! Magnitude of linear momentum, using corrected velocity
                    v = part%v
                    value = part%m * magnitude(v)
                case ('momentum_sphere_r')
                    ! Linear momentum in the spherical radial direction
                    ! Dot product of velocity vector with spherical phi
                    !    unit vector
                    ! 3. Multiply by mass of particle
                    v = part%v
                    call spherical_basis_from_cartesian(part%x,temp_basis)
                    value = part%m * (v .DOT. temp_basis%u(1))
                case ('ang_momentum_x')
                    ! Corrected angular momentum in the x direction
                    v = part%v
                    value = part%m * (part%x%y * v%z &
                                &- v%y * part%x%z)
                case ('ang_momentum_y')
                    ! Corrected angular momentum in the y direction
                    v = part%v
                    value = part%m * (part%x%z * v%x &
                                &- v%z * part%x%x)
                case ('ang_momentum_z')
                    ! Corrected angular momentum in the z direction
                    v = part%v
                    value = part%m * (part%x%x * v%y &
                                &- v%x * part%x%y)
                case ('ang_momentum')
                    ! Corrected magnitude of angular momentum
                    v = part%v
                    L = part%x * v
                    value = part%m * magnitude(L)
                case ('metallicity')
                    ! Metallicity
                    value = part%met
                case ('age')
                    ! Age
                    if (sim%cosmo) then
                        iii = 1
                        do while(sim%tau_frw(iii)>part%age.and.iii<sim%n_frw)
                            iii = iii + 1
                        end do
                        ! Interpolate time
#ifdef AGEPROPER
                        time = part%age
#else
                        time = sim%t_frw(iii)*(part%age-sim%tau_frw(iii-1))/(sim%tau_frw(iii)-sim%tau_frw(iii-1))+ &
                                & sim%t_frw(iii-1)*(part%age-sim%tau_frw(iii))/(sim%tau_frw(iii-1)-sim%tau_frw(iii))
#endif
                        value = (sim%time_simu-time)/(sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                        if (value<0d0) value = 0d0
                    endif
                case ('birth_date')
                    ! Birth date
                    if (sim%cosmo) then
                        iii = 1
                        do while(sim%tau_frw(iii)>part%age.and.iii<sim%n_frw)
                            iii = iii + 1
                        end do
                        ! Interpolate time
                        time = sim%t_frw(iii)*(part%age-sim%tau_frw(iii-1))/(sim%tau_frw(iii)-sim%tau_frw(iii-1))+ &
                                & sim%t_frw(iii-1)*(part%age-sim%tau_frw(iii))/(sim%tau_frw(iii-1)-sim%tau_frw(iii))
                        age = (sim%time_simu-value)/(sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                        value = (sim%time_tot+age)/(sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                    endif
                end select
            else
                ! This is for the case in which the SFR using a particular time indicator is required
                ! It should be used as 'star/sfr_xxx', where xxx is some multiple of Myr
                sfrstr = varname(index2+1:)
                index3 = scan(sfrstr,'_')
                if (index3.eq.0) then
                    read(sfrstr,'(F10.0)') sfrind
                    ! We want it in units of Gyr, that's why we divide by 1e+3
                    sfrind = sfrind/1D3
                    ! Birth date
                    if (sim%cosmo) then
                        iii = 1
                        do while(sim%tau_frw(iii)>part%age.and.iii<sim%n_frw)
                            iii = iii + 1
                        end do
                        ! Interpolate time
                        current_age_univ = (sim%time_tot+sim%time_simu)/(sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                        time = sim%t_frw(iii)*(part%age-sim%tau_frw(iii-1))/(sim%tau_frw(iii)-sim%tau_frw(iii-1))+ &
                                & sim%t_frw(iii-1)*(part%age-sim%tau_frw(iii))/(sim%tau_frw(iii-1)-sim%tau_frw(iii))
                        birth_date = (sim%time_tot+time)/(sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                    endif
                    ! Compute SFR by binning star particles by age
                    if (birth_date >= (current_age_univ - sfrind)) then
                        value = part%imass
                    else
                        value = 0D0
                    endif
                else
                    ! In the case that projections or profiles are obtained, we need volumetric
                    ! information (e.g. SFR density, SFR surface density)
                    sfrtype = sfrstr(1:index3-1)
                    sfrstr = sfrstr(index3+1:)
                    read(sfrstr,'(F10.0)')sfrind
                    ! We want it in units of Gyr, that's why we divide by 1e+3
                    sfrind = sfrind/1D3
                    ! Birth date
                    if (sim%cosmo) then
                        iii = 1
                        do while(sim%tau_frw(iii)>part%age.and.iii<sim%n_frw)
                            iii = iii + 1
                        end do
                        ! Interpolate time
                        current_age_univ = (sim%time_tot+sim%time_simu)/(sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
#ifdef AGEPROPER
                        time = part%age
#else
                        time = sim%t_frw(iii)*(part%age-sim%tau_frw(iii-1))/(sim%tau_frw(iii)-sim%tau_frw(iii-1))+ &
                                & sim%t_frw(iii-1)*(part%age-sim%tau_frw(iii))/(sim%tau_frw(iii-1)-sim%tau_frw(iii))
#endif
                        birth_date = (sim%time_tot+time)/(sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                    endif
                    ! write(*,*)birth_date, current_age_univ, part%age
                    ! Compute SFR by binning star particles by age
                    if (birth_date >= (current_age_univ - sfrind)) then
                        select case (TRIM(sfrtype))
                            case ('density')
                                if (present(dx)) then
                                    value = part%imass / (dx%x*dx%y*dx%z) / sfrind
                                else
                                    write(*,*)'Can not compute a particle density without cell size!'
                                endif
                            case ('surface')
                                if (present(dx)) then
                                    value = part%imass / (dx%x*dx%y) / sfrind
                                else
                                    write(*,*)'Can not compute a particle surface density without cell size!'
                                endif
                            case default
                                write(*,*)'This type of SFR is not recognised: ',TRIM(sfrtype)
                                write(*,*)'Aborting!'
                                stop
                        end select
                    else
                        value = 0D0
                    endif
                endif
            endif
        else
            value = 0D0
        endif
        
    end subroutine getpartvalue
end module io_ramses

module filtering
    use local
    use io_ramses

    type filter
        character(128) :: name
        integer :: ncond=0
        character(128), dimension(:), allocatable :: cond_vars
        character(128), dimension(:), allocatable :: cond_vars_comp
        character(2), dimension(:), allocatable :: cond_ops
        real(dbl), dimension(:), allocatable :: cond_vals
        logical, dimension(:), allocatable :: use_var
    end type filter

    contains

    subroutine allocate_filter(filt)
        implicit none
        type(filter), intent(inout) :: filt

        if (.not.allocated(filt%cond_vars)) allocate(filt%cond_vars(filt%ncond))
        if (.not.allocated(filt%cond_vars_comp)) allocate(filt%cond_vars_comp(filt%ncond))
        if (.not.allocated(filt%cond_ops)) allocate(filt%cond_ops(filt%ncond))
        if (.not.allocated(filt%cond_vals)) allocate(filt%cond_vals(filt%ncond))
        if (.not.allocated(filt%use_var)) allocate(filt%use_var(filt%ncond))
        filt%use_var = .false.
    end subroutine allocate_filter

    subroutine cond_string_to_filter(str, filt)
        implicit none
        character(128), intent(in) :: str
        type(filter), intent(inout) :: filt

        ! TODO: Finish this subroutine
    end subroutine cond_string_to_filter

    logical function filter_cell(reg,filt,cell_x,cell_dx,cell_var,cell_son,&
                                &trans_matrix,grav_var)
        use vectors
        use geometrical_regions
        type(region), intent(in) :: reg
        type(filter), intent(in) :: filt
        real(dbl), intent(in) :: cell_dx
        type(vector), intent(in) :: cell_x
        real(dbl), dimension(0:amr%twondim,1:varIDs%nvar), intent(in) :: cell_var
        integer,dimension(0:amr%twondim),intent(in) :: cell_son
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix
        real(dbl),dimension(0:amr%twondim,1:4),intent(in),optional :: grav_var
        integer :: i
        real(dbl) :: value,filt_value

        filter_cell = .true.

        if (filt%ncond == 0) return

        do i=1,filt%ncond
            if (present(grav_var)) then
                call getvarvalue(reg,cell_dx,cell_x,cell_var,cell_son,&
                                &filt%cond_vars(i),value,trans_matrix,grav_var)
            else
                call getvarvalue(reg,cell_dx,cell_x,cell_var,cell_son,&
                                    &filt%cond_vars(i),value,trans_matrix)
            end if
            if (filt%use_var(i)) then
                if (present(grav_var)) then
                    call getvarvalue(reg,cell_dx,cell_x,cell_var,cell_son,&
                                    &filt%cond_vars_comp(i),filt_value,trans_matrix,grav_var)
                else
                    call getvarvalue(reg,cell_dx,cell_x,cell_var,cell_son,&
                                        &filt%cond_vars_comp(i),filt_value,trans_matrix)
                end if
                filt_value = filt%cond_vals(i) * filt_value
            else
                filt_value = filt%cond_vals(i)
            end if
            select case (TRIM(filt%cond_ops(i)))
            case('/=')
                filter_cell = filter_cell .and. (value /= filt_value)
            case('==')
                filter_cell = filter_cell .and. (value == filt_value)
            case('<')
                filter_cell = filter_cell .and. (value < filt_value)
            case('<=')
                filter_cell = filter_cell .and. (value <= filt_value)
            case('>')
                filter_cell = filter_cell .and. (value > filt_value)
            case('>=')
                filter_cell = filter_cell .and. (value >= filt_value)
            case default
                write(*,*)'Relation operator not supported: ',TRIM(filt%cond_ops(i))
                write(*,*)'Aborting!'
                stop
            end select
        end do
    end function filter_cell

    logical function filter_particle(reg,filt,part,dx)
        use geometrical_regions
        type(region),intent(in) :: reg
        type(filter),intent(in) :: filt
        type(particle),intent(in) :: part
        type(vector),optional,intent(in) :: dx
        integer :: i
        real(dbl) :: value

        filter_particle = .true.
        if (filt%ncond == 0) return

        do i=1,filt%ncond
            if (present(dx)) then
                call getpartvalue(reg,part,filt%cond_vars(i),value,dx)
            else
                call getpartvalue(reg,part,filt%cond_vars(i),value)
            endif
            select case (TRIM(filt%cond_ops(i)))
            case('/=')
                filter_particle = filter_particle .and. (value /= filt%cond_vals(i))
            case('==')
                filter_particle = filter_particle .and. (value == filt%cond_vals(i))
            case('<')
                filter_particle = filter_particle .and. (value < filt%cond_vals(i))
            case('<=')
                filter_particle = filter_particle .and. (value <= filt%cond_vals(i))
            case('>')
                filter_particle = filter_particle .and. (value > filt%cond_vals(i))
            case('>=')
                filter_particle = filter_particle .and. (value >= filt%cond_vals(i))
            case default
                write(*,*)'Relation operator not supported: ',TRIM(filt%cond_ops(i))
                write(*,*)'Aborting!'
                stop
            end select
        end do
    end function filter_particle

    logical function filter_sub(sub,cell_x)
        use vectors
        use geometrical_regions
        type(region),intent(in) :: sub
        real(dbl),dimension(1:3), intent(in) :: cell_x
        integer :: i
        real(dbl) :: distance
        real(dbl),dimension(1:3) :: pos

        filter_sub = .True.

        pos = cell_x - (/sub%centre%x,sub%centre%y,sub%centre%z/)

        call checkifinside(pos,sub,filter_sub,distance)
        filter_sub = .not.filter_sub
    end function filter_sub

end module filtering