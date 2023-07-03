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
        real(dbl) :: h0,t,aexp,unit_l,unit_d,unit_t,unit_m,unit_v,unit_p
        real(dbl) :: boxlen,omega_m,omega_l,omega_k,omega_b
        real(dbl) :: time_tot,time_simu,redshift,T2,nH
        integer :: n_frw, nvar
        real(dbl),dimension(:),allocatable :: aexp_frw,hexp_frw,tau_frw,t_frw
        real(dbl) :: eta_sn=-1D0
        real(dbl) :: Dcr=3D28
    end type sim_info

    type level
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
    type(dictf90)  :: varIDs

    real(dbl)      :: Tmin,cV
    real(dbl)      :: lambda_crGH08

    contains

    subroutine retrieve_vars(repository,myvars)
        implicit none

        character(128),intent(in) ::  repository
        type(dictf90), intent(inout) :: myvars

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
            write(*,*)': This simulation uses particle families'
            sim%family = .true.
        else if (trim(var6) .eq. 'tform') then
            write(*,*)': This simulation uses the old particle format'
        else
            write(*,*)': This simulation format for particles is not recognised!'
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
                write(*,*) "This is an old RAMSES simulation"
                call read_hydrofile_descriptor_old(repository)
            elseif (igrstart .eq. "#") then
                write(*,*) "This is a new RAMSES simulation"
                call read_hydrofile_descriptor_new(repository)
            else
                write(*,*)" I do not recognise this sim format. Check!"
                stop
            endif
        else
            write(*,'(": ",A," not found. Stopping!")') trim(nomfich)
            stop
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
        write(*,'(": Reading variables IDs from hydro_descriptor")')
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
        write(*,*)'nvar=',nvar
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
            write(*,'(": Reading eta_sn from namelist.txt")')
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
                    write(*,*)': Found eta_sn=',sim%eta_sn
                    exit
                end if
            end do
        else
            write(*,'(": Namelist not found: eta_sn set to 0.2")')
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
            write(*,*)': No hydro data in this simulation.'
        else
            sim%hydro = .true.
        endif
        nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
        ! Verify input part file
        inquire(file=nomfich, exist=ok)
        if ( .not. ok ) then
            print *,TRIM(nomfich)//' not found.'
            write(*,*)': No particles in this simulation.'
        else
            sim%dm = .true.
        endif
        if (sim%dm .and. .not. sim%hydro) write(*,*)': This is a DM only simulation.'
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
    
        read(10,'("boxlen      =",E23.15)')sim%boxlen
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
        sim%unit_m = ((sim%unit_d*sim%unit_l)*sim%unit_l)*sim%unit_l
        read(10,'("ordering type=",A80)')amr%ordering
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
            Tmin = 15d0 / sim%T2
            cV = XH * cVHydrogen * mHydrogen / kBoltzmann
            lambda_crGH08 = 2.63d-16 * ((sim%unit_t**3)/(sim%unit_d*(sim%unit_l**2)))
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
        write(*,*)'limits:',xxmin,xxmax,yymin,yymax,zzmin,zzmax
        write(*,*)'ordering: ',TRIM(amr%ordering)
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
                igridn(i,j)=son(nbor(igrid(i),j))
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

        ! write(*,*)'ncoarse,ngridmax,twondim: ',amr%ncoarse,amr%ngridmax,amr%twondim
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