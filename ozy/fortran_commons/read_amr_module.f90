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
    use dictionary_commons
    use constants
    use vectors
    use cooling_module

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
        logical :: isthere_part_descriptor=.false.
        real(dbl) :: h0,t,aexp,unit_l,unit_d,unit_t,unit_m,unit_v,unit_p
        real(dbl) :: boxlen,omega_m,omega_l,omega_k,omega_b
        real(dbl) :: time_tot,time_simu,redshift,T2,nH
        integer :: n_frw, nvar, nvar_part, nvar_part_d, nvar_part_i, nvar_part_b
        integer,dimension(:),allocatable :: part_var_types
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
        real(dbl),dimension(:,:,:,:,:),pointer::cube
        real(dbl),dimension(:,:,:,:),pointer::map
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

    type rt_info
        integer :: nRTvar=0
        integer :: nIons=0
        integer :: nGroups=0
        integer :: iIons=0
        integer :: rtdp=0
        real(dbl) :: X=0d0, Y=0d0
        real(dbl) :: scale_np=0d0, scale_pf=0d0, rt_c_fraction=0d0
        real(dbl) :: n_star=0d0, T2_star=0d0, g_star=0d0
        integer,dimension(:),allocatable :: spec2group
        real(dbl),dimension(:),allocatable :: groupL0,groupL1,group_egy
        real(dbl),dimension(:,:),allocatable :: group_csn,group_cse
    end type rt_info

    type particle
#ifndef LONGINT
        integer(irg) :: id,level
#else
        integer(ilg) :: id,level
#endif
        type(vector) :: x,v
        real(dbl) :: m,met,imass,age,tform
    end type particle

    ! Define global variables
    type(sim_info) :: sim
    type(amr_info) :: amr
    type(rt_info)  :: rtinfo
    type(dictf90)  :: varIDs
    type(dictf90)  :: partIDs,partvar_types
    type(dictf90)  :: rtIDs
    logical :: verbose=.false.
    logical :: fix_neg_temp=.true.
    real(dbl)      :: Tmin,cV
    real(dbl)      :: lambda_crGH08

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
        type(dictf90), intent(inout) :: myvars

        call read_hydrofile_descriptor(repository)
        myvars = varIDs
    end subroutine retrieve_vars

    subroutine retrieve_rtinfo(repository,nchar,myrtinfo)
        implicit none
        character(128),intent(in) :: repository
        character(5), intent(in) :: nchar
        type(rt_info), intent(inout) :: myrtinfo

        ! Read the info file into module rtinfo, then copy to the passed object
        call read_rt_info(repository,nchar)
        myrtinfo = rtinfo
    end subroutine retrieve_rtinfo

    subroutine retrieve_rtIDs(repository,nchar,myrtids)
        implicit none

        character(128),intent(in) :: repository
        character(5), intent(in) :: nchar
        type(dictf90), intent(inout) :: myrtids

        ! Ensure rt info and rtIDs are populated
        call read_rt_info(repository,nchar)
        myrtids = rtIDs
    end subroutine retrieve_rtIDs

    subroutine print_rtinfo(repository,nchar)
        implicit none
        character(128),intent(in) :: repository
        character(5), intent(in) :: nchar
        integer :: i,j

        ! Check for the presence of RT files
        call check_rt_files(repository,nchar)

        if (.not. sim%rt) then
            write(*,*) 'No RT files detected for repository:', trim(repository)
            return
        end if

        ! Read rt info into module rtinfo
        call read_rt_info(repository,nchar)

        write(*,*) '---- RT info ----'
        write(*,'("nRTvar      =",I6)') rtinfo%nRTvar
        write(*,'("nIons       =",I6)') rtinfo%nIons
        write(*,'("nGroups     =",I6)') rtinfo%nGroups
        write(*,'("iIons       =",I6)') rtinfo%iIons
        write(*,'("rtprecision =",I6)') rtinfo%rtdp
        write(*,'("X_fraction  =",E23.15)') rtinfo%X
        write(*,'("Y_fraction  =",E23.15)') rtinfo%Y
        write(*,*)
        write(*,'("unit_np     =",E23.15)') rtinfo%scale_np
        write(*,'("unit_pf     =",E23.15)') rtinfo%scale_pf
        write(*,'("rt_c_frac   =",E23.15)') rtinfo%rt_c_fraction
        write(*,*)
        write(*,'("n_star      =",E23.15)') rtinfo%n_star
        write(*,'("T2_star     =",E23.15)') rtinfo%T2_star
        write(*,'("g_star      =",E23.15)') rtinfo%g_star
        write(*,*)
        if (allocated(rtinfo%groupL0)) then
            write(*,*) 'groupL0:'
            write(*,'(20E12.3)') (rtinfo%groupL0(i), i=1,min(size(rtinfo%groupL0),20))
        end if
        if (allocated(rtinfo%groupL1)) then
            write(*,*) 'groupL1:'
            write(*,'(20E12.3)') (rtinfo%groupL1(i), i=1,min(size(rtinfo%groupL1),20))
        end if
        if (allocated(rtinfo%spec2group)) then
            write(*,*) 'spec2group:'
            write(*,'(20I6)') (rtinfo%spec2group(i), i=1,size(rtinfo%spec2group))
        end if

        if (allocated(rtinfo%group_egy)) then
            write(*,*) 'Per-group properties:'
            do i=1,rtinfo%nGroups
                write(*,'("--- Group ",I2)') i
                write(*,'("  egy =",1PE12.3)') rtinfo%group_egy(i)
                if (allocated(rtinfo%group_csn)) then
                    write(*,'("  csn =",20(1PE12.3))') (rtinfo%group_csn(i,j), j=1,min(20,size(rtinfo%group_csn,2)))
                end if
                if (allocated(rtinfo%group_cse)) then
                    write(*,'("  cse =",20(1PE12.3))') (rtinfo%group_cse(i,j), j=1,min(20,size(rtinfo%group_cse,2)))
                end if
            end do
        end if

        write(*,*) '---- end RT info ----'

    end subroutine print_rtinfo


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

        if (varIDs%count.ne.0) return 

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
            write(*,'(": ",A," not found. Stopping!")') trim(nomfich)
            stop
        end if
    end subroutine read_hydrofile_descriptor

    !---------------------------------------------------------------
    ! Subroutine: READ PART IDs
    ! This routine extracts the PART variable IDs for a given
    ! simulation snapshot such that their order is known when
    ! part_*.out* files are read. This only applies if the RAMSES
    ! version includes a standarised part_file_descriptor.txt
    !---------------------------------------------------------------
    subroutine read_partfile_descriptor(repository)
        implicit none

        character(128),intent(in) :: repository
        character(256)            :: nomfich
        logical                   :: ok
        character(25)  ::  newVar,newType
        integer            ::  newID,status,i,TypeID

        if (partIDs%count.ne.0) return

        nomfich=TRIM(repository)//'/part_file_descriptor.txt'
        inquire(file=nomfich, exist=ok) ! verify input file
        if (ok) then
            sim%nvar_part_d = 0; sim%nvar_part_i = 0; sim%nvar_part_b = 0
            sim%nvar_part = 0
            if (verbose) write(*,'(": Reading part IDs from part_file_descriptor.txt")')
            open(unit=111,file=nomfich,status='old',form='formatted')
            read(111,*) ! Skip header
            read(111,*) ! Skip header
            do
                read(111,*,iostat=status)newID,newVar,newType
                if (status /= 0) exit
                sim%nvar_part = sim%nvar_part + 1
                if (trim(newType) .eq. 'd') then
                    TypeID = 1
                    sim%nvar_part_d = sim%nvar_part_d + 1
                    call partIDs%add(newVar,sim%nvar_part_d)
                elseif (trim(newType) .eq. 'i') then
                    TypeID = 2
                    sim%nvar_part_i = sim%nvar_part_i + 1
                    call partIDs%add(newVar,sim%nvar_part_i)
                elseif (trim(newType) .eq. 'b') then
                    TypeID = 3
                    sim%nvar_part_b = sim%nvar_part_b + 1
                    call partIDs%add(newVar,sim%nvar_part_b)
                else
                    write(*,'(": ",A," particle type not found. Stopping!")')
                    stop
                end if
                sim%part_var_types(sim%nvar_part) = TypeID
                call partvar_types%add(newVar,TypeID)
            end do
            close(111)
            if (verbose) then
                write(*,*)'nvar_part=',sim%nvar_part
                write(*,*)'nvar_part_d=',sim%nvar_part_d
                write(*,*)'nvar_part_i=',sim%nvar_part_i
                write(*,*)'nvar_part_b=',sim%nvar_part_b
            end if
            sim%isthere_part_descriptor = .true.
        else
            write(*,'(": ",A," not found. Using default and information from header.txt!")') trim(nomfich)
            sim%nvar_part_d = 0; sim%nvar_part_i = 0; sim%nvar_part_b = 0
            sim%nvar_part = 0
            call read_headerfile(repository)
            sim%isthere_part_descriptor = .false.
        end if
    end subroutine read_partfile_descriptor

    !---------------------------------------------------------------
    ! Subroutine: READ HEADER
    !
    ! This routine reads the header_xxxxx.txt file from a snapshot
    ! in order to determine the ordering and number of particle
    ! variables used in the RAMSES simulation.
    !---------------------------------------------------------------
    subroutine read_headerfile(repository)
        use utils, only: get_word
        implicit none
        character(len=128), intent(in) :: repository
        integer :: iunit, ios, i, ipos, num_vars
        character(len=256) :: line
        character(256)            :: nomfich
        character(len=50), dimension(:), allocatable :: variables
        character(len=25) :: varname
        character(5) :: nchar
        logical :: found_header
        integer :: rt_nkeys, idcount
        character(len=32) :: keyname

        if (partIDs%count.ne.0) return

        ipos = index(repository,'output_')
        nchar=repository(ipos+7:ipos+13)
        nomfich=TRIM(repository)//'/header_'//TRIM(nchar)//'.txt'
    
        ! Open the file
        iunit = 24
        open(unit=iunit, file=nomfich, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print*, "Error opening file: ", nomfich
            stop
        end if
    
        found_header = .false.
    
        ! Read file line by line
        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
    
            ! Check if the line starts with "pos"
            if (trim(line(1:3)) == 'pos') then
                found_header = .true.
                exit
            end if
        end do
    
        ! If header not found, exit
        if (.not. found_header) then
            print*, "Error: No 'pos' header found in file."
            close(iunit)
            stop
        end if
    
        ! Split the header line into variable names
        num_vars = 0
        allocate(variables(50))  ! Allocate enough space for variables
    
        i = 1
        do
            call get_word(line, i, variables(num_vars+1))
            if (len_trim(variables(num_vars+1)) == 0) exit
            num_vars = num_vars + 1
        end do
    
        ! Close the file
        close(iunit)
        ! Build rtIDs dictionary: entries to access photon density and fluxes
        if (rtinfo%nGroups > 0) then
            rt_nkeys = rtinfo%nGroups * (1 + amr%ndim)
            if (rt_nkeys < 1) rt_nkeys = 1
            ! Initialize rtIDs dictionary (resize)
            call rtIDs%init(rt_nkeys)

            idcount = 0
            do i = 1, rtinfo%nGroups
                ! base key for this group, e.g. 'rt_g1'
                write(keyname, '(A,I0)') 'rt_g', i
                idcount = idcount + 1
                call rtIDs%add(trim(keyname)//'_n', idcount)  ! photon density
                if (amr%ndim >= 1) then
                    idcount = idcount + 1
                    call rtIDs%add(trim(keyname)//'_fx', idcount)
                end if
                if (amr%ndim >= 2) then
                    idcount = idcount + 1
                    call rtIDs%add(trim(keyname)//'_fy', idcount)
                end if
                if (amr%ndim >= 3) then
                    idcount = idcount + 1
                    call rtIDs%add(trim(keyname)//'_fz', idcount)
                end if
            end do
        end if
        
        ! Now loop over the variables assigning them to the particle
        ! variables dictionary
        if (.not.allocated(sim%part_var_types)) then
            if (amr%ndim>2) then
                allocate(sim%part_var_types(num_vars+4))
                call partIDs%init(num_vars+4)
                call partvar_types%init(num_vars+4)
            elseif (amr%ndim>1) then
                allocate(sim%part_var_types(num_vars+2))
                call partIDs%init(num_vars+2)
                call partvar_types%init(num_vars+2)
            else
                allocate(sim%part_var_types(num_vars))
                call partIDs%init(num_vars)
                call partvar_types%init(num_vars)
            end if
        end if
        
        do i = 1, num_vars
            if (trim(variables(i)) == 'pos') then
                sim%nvar_part = sim%nvar_part + 1
                sim%nvar_part_d = sim%nvar_part_d + 1
                call partIDs%add('x',sim%nvar_part_d)
                call partvar_types%add('x',1)
                sim%part_var_types(sim%nvar_part) = 1
                if (amr%ndim > 1) then
                    sim%nvar_part = sim%nvar_part + 1
                    sim%nvar_part_d = sim%nvar_part_d + 1
                    call partIDs%add('y',sim%nvar_part_d)
                    call partvar_types%add('y',1)
                    sim%part_var_types(sim%nvar_part) = 1
                end if
                if (amr%ndim >2) then
                    sim%nvar_part = sim%nvar_part + 1
                    sim%nvar_part_d = sim%nvar_part_d + 1
                    call partIDs%add('z',sim%nvar_part_d)
                    call partvar_types%add('z',1)
                    sim%part_var_types(sim%nvar_part) = 1
                end if
            else if (trim(variables(i)) == 'vel') then
                sim%nvar_part = sim%nvar_part + 1
                sim%nvar_part_d = sim%nvar_part_d + 1
                call partIDs%add('velocity_x',sim%nvar_part_d)
                call partvar_types%add('velocity_x',1)
                sim%part_var_types(sim%nvar_part) = 1
                if (amr%ndim > 1) then
                    sim%nvar_part = sim%nvar_part + 1
                    sim%nvar_part_d = sim%nvar_part_d + 1
                    call partIDs%add('velocity_y',sim%nvar_part_d)
                    call partvar_types%add('velocity_y',1)
                    sim%part_var_types(sim%nvar_part) = 1
                end if
                if (amr%ndim >2) then
                    sim%nvar_part = sim%nvar_part + 1
                    sim%nvar_part_d = sim%nvar_part_d + 1
                    call partIDs%add('velocity_z',sim%nvar_part_d)
                    call partvar_types%add('velocity_z',1)
                    sim%part_var_types(sim%nvar_part) = 1
                end if
            else if (trim(variables(i)) == 'mass') then
                sim%nvar_part = sim%nvar_part + 1
                sim%nvar_part_d = sim%nvar_part_d + 1
                call partIDs%add('mass',sim%nvar_part_d)
                call partvar_types%add('mass',1)
                sim%part_var_types(sim%nvar_part) = 1
            else if (trim(variables(i)) == 'iord') then
                sim%nvar_part = sim%nvar_part + 1
                sim%nvar_part_i = sim%nvar_part_i + 1
                call partIDs%add('iord',sim%nvar_part_i)
                call partvar_types%add('iord',2)
                sim%part_var_types(sim%nvar_part) = 2
            else if (trim(variables(i)) == 'level') then
                sim%nvar_part = sim%nvar_part + 1
                sim%nvar_part_i = sim%nvar_part_i + 1
                call partIDs%add('level',sim%nvar_part_i)
                call partvar_types%add('level',2)
                sim%part_var_types(sim%nvar_part) = 2
            else if (trim(variables(i)) == 'family') then
                sim%nvar_part = sim%nvar_part + 1
                sim%nvar_part_b = sim%nvar_part_b + 1
                call partIDs%add('family',sim%nvar_part_b)
                call partvar_types%add('family',3)
                sim%part_var_types(sim%nvar_part) = 3
            else if (trim(variables(i)) == 'tag') then
                sim%nvar_part = sim%nvar_part + 1
                sim%nvar_part_b = sim%nvar_part_b + 1
                call partIDs%add('tag',sim%nvar_part_b)
                call partvar_types%add('tag',3)
                sim%part_var_types(sim%nvar_part) = 3
            else if (trim(variables(i)) == 'tform') then
                sim%nvar_part = sim%nvar_part + 1
                sim%nvar_part_d = sim%nvar_part_d + 1
                call partIDs%add('birth_time',sim%nvar_part_d)
                call partvar_types%add('birth_time',1)
                sim%part_var_types(sim%nvar_part) = 1
            else if (trim(variables(i)) == 'metal') then
                sim%nvar_part = sim%nvar_part + 1
                sim%nvar_part_d = sim%nvar_part_d + 1
                call partIDs%add('metallicity',sim%nvar_part_d)
                call partvar_types%add('metallicity',1)
                sim%part_var_types(sim%nvar_part) = 1
            else if (trim(variables(i)) == 'imass') then
                sim%nvar_part = sim%nvar_part + 1
                sim%nvar_part_d = sim%nvar_part_d + 1
                call partIDs%add('initial_mass',sim%nvar_part_d)
                call partvar_types%add('initial_mass',1)
                sim%part_var_types(sim%nvar_part) = 1
            else
                sim%nvar_part = sim%nvar_part + 1
                sim%nvar_part_d = sim%nvar_part_d + 1
                call partIDs%add(trim(variables(i)),sim%nvar_part_d)
                call partvar_types%add(trim(variables(i)),1)
                sim%part_var_types(sim%nvar_part) = 1
            end if

        end do

        
        if (verbose) then
            write(*,*)'nvar_part=',sim%nvar_part
            write(*,*)'nvar_part_d=',sim%nvar_part_d
            write(*,*)'nvar_part_i=',sim%nvar_part_i
            write(*,*)'nvar_part_b=',sim%nvar_part_b
        end if

    end subroutine read_headerfile


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
        write(*,'(": Reading variables IDs from hydro_descriptor")')
        open(unit=10,file=nomfich,status='old',form='formatted')
        read(10,*) igr9,igr1,igr2
        read(igr2,*,iostat=statn) nvhydro
        write(*,*)'nvar=',nvhydro
        sim%nvar = nvhydro
        if (nvhydro > 9) then
            nvloop = 9
        else
            nvloop = nvhydro
        end if
        do i=1,nvloop
            read(10,*) igr9,igr1,igr1,newVar
            read(igr1,*,iostat=statn) newID
            call varIDs%add(newVar,newID)
        end do
        if (nvhydro > 10) then
            do i=nvloop+1,nvhydro
                read(10,*) igr8,igr3,newVar
                igr3 = igr3(2:3);
                read(igr3,*,iostat=statn) newID
                call varIDs%add(newVar,newID)
            end do
        end if
        close(10)
    end subroutine read_hydrofile_descriptor_old

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
        if (verbose) write(*,'(": Reading variables IDs from hydro_descriptor")')
        open(unit=10,file=nomfich,status='old',form='formatted')
        read(10,*)
        read(10,*)
        do
            read(10,*,iostat=status)newID,newVar,newType
            if (status /= 0) exit
            nvar = nvar + 1
            call varIDs%add(newVar,newID)
        end do
        close(10)
        if (verbose) write(*,*)'nvar=',nvar
        sim%nvar = nvar
    end subroutine read_hydrofile_descriptor_new

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
        logical :: ok

        nomfich=TRIM(repository)//'/namelist.txt'
        
        inquire(file=nomfich, exist=ok) ! verify input file
        if (ok) then
            if (verbose) write(*,'(": Reading eta_sn from namelist.txt")')
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
                    exit
                end if
            end do
        else
            if (verbose) write(*,'(": Namelist not found: eta_sn set to 0.2")')
            sim%eta_sn = 2D-1
        end if

    end subroutine get_eta_sn

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
            if (sim%cr) call get_Dcr(repository)
            Tmin = 15d0 / sim%T2
            cV = XH * cVHydrogen * mHydrogen / kBoltzmann
            lambda_crGH08 = 2.63d-16 * ((sim%unit_t**3)/(sim%unit_d*(sim%unit_l**2)))
        endif

        ! Check for the presence of RT files
        call check_rt_files(repository,nchar)
        if (sim%rt) call read_rt_info(repository,nchar)
    end subroutine init_amr_read

    !---------------------------------------------------------------
    ! Subroutine: READ RT INFO
    !
    ! Read the info_rt_XXXXX.txt file produced by `output_rtInfo`
    ! and store contents in the global `rtinfo` object.
    !---------------------------------------------------------------
    subroutine read_rt_info(repository,nchar)
        use utils, only: get_word
        implicit none
        character(128), intent(in) :: repository
        character(5), intent(in) :: nchar
        character(128) :: nomfich
        integer :: iunit, ios, i, j, k, nTok
        character(len=512) :: line
        character(len=80) :: token
        integer :: pos
        character(len=512) :: chunk
        integer :: len_chunk, jend
        integer, dimension(20) :: tmpi
        real(dbl), dimension(20) :: tmpr

        nomfich = TRIM(repository)//'/info_rt_'//TRIM(nchar)//'.txt'

        ! Try to open the file
        iunit = 30
        open(unit=iunit, file=nomfich, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            if (verbose) write(*,*)'info_rt file not found:', TRIM(nomfich)
            return
        end if

        ! Initialize/clear any existing rtinfo arrays
        if (allocated(rtinfo%spec2group)) deallocate(rtinfo%spec2group)
        if (allocated(rtinfo%groupL0)) deallocate(rtinfo%groupL0)
        if (allocated(rtinfo%groupL1)) deallocate(rtinfo%groupL1)
        if (allocated(rtinfo%group_egy)) deallocate(rtinfo%group_egy)
        if (allocated(rtinfo%group_csn)) deallocate(rtinfo%group_csn)
        if (allocated(rtinfo%group_cse)) deallocate(rtinfo%group_cse)

        do
            read(iunit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            ! Find and parse known keys
            if (index(line,'nRTvar') /= 0) then
                pos = index(line,'=')
                read(line(pos+1:),*) rtinfo%nRTvar
            else if (index(line,'nIons') /= 0) then
                pos = index(line,'=')
                read(line(pos+1:),*) rtinfo%nIons
                ! if groups already known, (re)allocate cs arrays to have correct second dim
                if (rtinfo%nIons > 0 .and. rtinfo%nGroups > 0) then
                    if (allocated(rtinfo%group_csn)) then
                        deallocate(rtinfo%group_csn)
                    end if
                    if (allocated(rtinfo%group_cse)) then
                        deallocate(rtinfo%group_cse)
                    end if
                    allocate(rtinfo%group_csn(rtinfo%nGroups, rtinfo%nIons))
                    allocate(rtinfo%group_cse(rtinfo%nGroups, rtinfo%nIons))
                end if
                else if (index(line,'nGroups') /= 0) then
                pos = index(line,'=')
                read(line(pos+1:),*) rtinfo%nGroups
                ! allocate group arrays once nGroups known
                if (rtinfo%nGroups > 0) then
                    if (.not. allocated(rtinfo%groupL0)) allocate(rtinfo%groupL0(rtinfo%nGroups))
                    if (.not. allocated(rtinfo%groupL1)) allocate(rtinfo%groupL1(rtinfo%nGroups))
                    if (.not. allocated(rtinfo%group_egy)) allocate(rtinfo%group_egy(rtinfo%nGroups))
                    ! allocate cs arrays when nIons is known; for now allocate with 1 col if nIons not yet known
                    if (rtinfo%nIons > 0) then
                        if (.not. allocated(rtinfo%group_csn)) allocate(rtinfo%group_csn(rtinfo%nGroups,rtinfo%nIons))
                        if (.not. allocated(rtinfo%group_cse)) allocate(rtinfo%group_cse(rtinfo%nGroups,rtinfo%nIons))
                    else
                        if (.not. allocated(rtinfo%group_csn)) allocate(rtinfo%group_csn(rtinfo%nGroups,1))
                        if (.not. allocated(rtinfo%group_cse)) allocate(rtinfo%group_cse(rtinfo%nGroups,1))
                    end if
                end if
            else if (index(line,'iIons') /= 0) then
                pos = index(line,'=')
                read(line(pos+1:),*) rtinfo%iIons
            else if (index(line,'rtprecision') /= 0) then
                pos = index(line,'=')
                read(line(pos+1:),*) rtinfo%rtdp
            else if (index(line,'X_fraction') /= 0) then
                pos = index(line,'=')
                read(line(pos+1:),*) rtinfo%X
            else if (index(line,'Y_fraction') /= 0) then
                pos = index(line,'=')
                read(line(pos+1:),*) rtinfo%Y
            else if (index(line,'unit_np') /= 0) then
                pos = index(line,'=')
                read(line(pos+1:),*) rtinfo%scale_np
            else if (index(line,'unit_pf') /= 0) then
                pos = index(line,'=')
                read(line(pos+1:),*) rtinfo%scale_pf
            else if (index(line,'rt_c_frac') /= 0) then
                pos = index(line,'=')
                read(line(pos+1:),*) rtinfo%rt_c_fraction
            else if (index(line,'n_star') /= 0) then
                pos = index(line,'=')
                read(line(pos+1:),*) rtinfo%n_star
            else if (index(line,'T2_star') /= 0) then
                pos = index(line,'=')
                read(line(pos+1:),*) rtinfo%T2_star
            else if (index(line,'g_star') /= 0) then
                pos = index(line,'=')
                read(line(pos+1:),*) rtinfo%g_star
            else if (index(line,'groupL0') /= 0) then
                pos = index(line,'=')
                if (rtinfo%nGroups > 0) then
                    ! writer uses 20F12.3 after the label; parse fixed-width F12.3 fields
                    k = 0
                    chunk = line(pos+1:)
                    len_chunk = len_trim(chunk)
                    do i = 1, min(20, rtinfo%nGroups)
                        k = k + 1
                        j = (k-1)*12 + 1
                        if (j > len_chunk) exit
                        jend = min(j+11, len_chunk)
                        token = adjustl(chunk(j:jend))
                        read(token, '(F12.3)', iostat=ios) tmpr(k)
                        if (ios /= 0) then
                            exit
                        end if
                        rtinfo%groupL0(k) = tmpr(k)
                    end do
                end if
            else if (index(line,'groupL1') /= 0) then
                pos = index(line,'=')
                if (rtinfo%nGroups > 0) then
                    k = 0
                    chunk = line(pos+1:)
                    len_chunk = len_trim(chunk)
                    do i = 1, min(20, rtinfo%nGroups)
                        k = k + 1
                        j = (k-1)*12 + 1
                        if (j > len_chunk) exit
                        jend = min(j+11, len_chunk)
                        token = adjustl(chunk(j:jend))
                        read(token, '(F12.3)', iostat=ios) tmpr(k)
                        if (ios /= 0) then
                            exit
                        end if
                        rtinfo%groupL1(k) = tmpr(k)
                    end do
                end if
            else if (index(line,'spec2group') /= 0) then
                pos = index(line,'=')
                ! spec2group is written as 20I12; only first nIons are relevant
                if (rtinfo%nIons <= 0) then
                    cycle
                end if
                if (.not. allocated(rtinfo%spec2group)) allocate(rtinfo%spec2group(rtinfo%nIons))
                chunk = line(pos+1:)
                len_chunk = len_trim(chunk)
                do j = 1, rtinfo%nIons
                    k = (j-1)*12 + 1
                    if (k > len_chunk) then
                        rtinfo%spec2group(j) = 0
                    else
                        jend = min(k+11, len_chunk)
                        token = adjustl(chunk(k:jend))
                        read(token, '(I12)', iostat=ios) tmpi(j)
                        if (ios /= 0) then
                            rtinfo%spec2group(j) = 0
                        else
                            rtinfo%spec2group(j) = tmpi(j)
                        end if
                    end if
                end do
            else if (index(line,'egy') /= 0 .and. index(line,'group_egy')==0) then
                ! This line label in write_group_props is 'egy      [eV] ='
                pos = index(line,'=')
                ! Skip: per-group energies are read in the group blocks below
                cycle
            else if (index(line,'csn') /= 0) then
                ! group_csn line: multiple lines per group; handle in loop below
                ! The writer outputs group_csn(ip,:) and group_cse(ip,:) in the loop.
                ! To capture these we re-read the file from current position storing per-group
                ! Move file back to current beginning isn't trivial; instead, continue reading
                ! and when encountering the per-group markers, use a simple state machine.
                cycle
            end if
        end do

        rewind(iunit)

        ! Re-scan file to capture per-group properties (group_egy, group_csn, group_cse)
        i = 0
        do
            read(iunit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            if (index(line,'---Group') /= 0) then
                ! Next lines: group_egy, group_csn, group_cse belong to group ip
                i = i + 1
                ! read next non-blank line expecting 'egy' for this group
                read(iunit,'(A)',iostat=ios) line
                if (ios /= 0) exit
                pos = index(line,'=')
                if (pos > 0) then
                    read(line(pos+1:),*) rtinfo%group_egy(i)
                end if
                ! read csn line (specifically nIons entries)
                read(iunit,'(A)',iostat=ios) line
                if (ios /= 0) exit
                pos = index(line,'=')
                if (pos > 0) then
                    do j = 1, max(1, rtinfo%nIons)
                        k = (j-1)*12 + 1
                        if (k > len_trim(line(pos+1:))) then
                            rtinfo%group_csn(i,j) = 0.0d0
                        else
                            token = adjustl(line(pos+k:pos+k+11))
                            read(token, '(1PE12.3)', iostat=ios) rtinfo%group_csn(i,j)
                            if (ios /= 0) rtinfo%group_csn(i,j) = 0.0d0
                        end if
                    end do
                end if
                ! read cse line
                read(iunit,'(A)',iostat=ios) line
                if (ios /= 0) exit
                pos = index(line,'=')
                if (pos > 0) then
                    do j = 1, max(1, rtinfo%nIons)
                        k = (j-1)*12 + 1
                        if (k > len_trim(line(pos+1:))) then
                            rtinfo%group_cse(i,j) = 0.0d0
                        else
                            token = adjustl(line(pos+k:pos+k+11))
                            read(token, '(1PE12.3)', iostat=ios) rtinfo%group_cse(i,j)
                            if (ios /= 0) rtinfo%group_cse(i,j) = 0.0d0
                        end if
                    end do
                end if
                if (i >= rtinfo%nGroups) exit
            end if
        end do

        close(iunit)
        if (verbose) write(*,*)'Read RT info from: ',TRIM(nomfich)
    end subroutine read_rt_info

    !---------------------------------------------------------------
    ! Subroutine: CHECK RT FILES
    !
    ! This routine checks the output directory for the presence of
    ! files starting with "rt_". If found, it sets the variable
    ! sim%rt to .true.
    !---------------------------------------------------------------
    subroutine check_rt_files(repository,nchar)
        implicit none
        character(128), intent(in) :: repository
        character(5), intent(in) :: nchar
        character(128) :: filename
        logical :: file_found     

        file_found = .false.

        write(filename, '(A,"/rt_",A,".out",I5.5)') trim(repository), trim(nchar), 1
        inquire(file=filename, exist=file_found)
        if (file_found) then
            sim%rt = .true.
        end if      
    end subroutine check_rt_files

    !---------------------------------------------------------------
    ! Subroutine: GET DCR
    !
    ! This routine reads the namelist.nml file from a snapshot
    ! in order to determine the value of the parameter Dcr
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

        nomfich=TRIM(repository)//'/namelist.txt'
        
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
            if (verbose) write(*,'(": Namelist not found: Dcr set to 3D28")')
            sim%Dcr = 3D28
        end if

    end subroutine get_Dcr

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
    
end module io_ramses
