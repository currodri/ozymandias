program HalosExtractor
  implicit none
  integer(kind=4),allocatable      :: first_part(:),nb_of_parts(:)
  integer(kind=4),allocatable      :: linked_list(:), liste_parts(:)
  integer(kind=4)                  :: nb_of_halos, nb_of_subhalos
  integer(kind=4)                  :: nbodies
  real(kind=4)                     :: massp
  real(kind=4)                     :: age_univ,Lbox_pt,Lbox_pt2,Hub_pt,omega_0,hubble,omega_lambda_0
  real(kind=4)                     :: omega_t,omega_lambda_t,omega_f,omega_lambda_f,omega_c_f
  real(kind=4)                     :: rho_crit,aexp,Lboxp,mboxp,af,ai,Lf,H_f,H_i
  type vector
     real(kind=8)      :: x,y,z
  end type vector
  type shape
     real(kind=4)      :: a,b,c
  end type shape
  type baryon
     real(kind=4)      :: rvir,mvir,tvir,cvel
  end type baryon
  type hprofile
     real(kind=4)      :: rho_0,r_c
  end type hprofile
  type halo
     type (baryon)     :: datas
     type (shape)      :: sh
     type (vector)     :: p
     type (vector)     :: v
     type (vector)     :: L
     type (hprofile)   :: halo_profile
     integer(kind=4)   :: my_number
     integer(kind=4)   :: my_timestep
     integer(kind=4)   :: nbsub
     integer(kind=4)   :: hosthalo
     integer(kind=4)   :: hostsub
     integer(kind=4)   :: level
     integer(kind=4)   :: nextsub
     real(kind=4)      :: m
     real(kind=4)      :: r
     real(kind=4)      :: spin
     real(kind=4)      :: ek,ep,et
  end type halo
  type (halo),allocatable          :: liste_halos(:)
  !======================================================================
  ! Definitions specific to input/output
  !======================================================================
  character(len=80)         :: data_dir
  character(len=3)          :: file_num
  integer(kind=4)           :: numstep
  integer(kind=4),parameter :: errunit = 0
  logical                   :: fsub         ! flag to notify whether subhaloes are included
  logical                   :: no_rewrite


  ! Filename variables
  character(len=22)         :: filenameIn
  character(len=28)         :: filenameOut
  character(len=128)        :: execFolder=""
  integer                   :: ifile
  integer,parameter         :: maxfiles = 999
  logical                   :: enquiref, enquirerf
  logical                   :: do_file
  integer                   :: ngs, ndm, nqq
  logical :: working=.false.

  ! Received folder
  integer       :: n
  integer       :: iargc, tempnum
  character(len=4)   :: opt
  character(len=128) :: arg
  
  write(*,'("###############################################")')
  write(*,'("                 Halo Extractor                ")')
  write(*,'("###############################################")')

  ! Receive folder for execution
  n = iargc()
  if (n.gt.0) then
     call getarg(1,arg)
     execFolder = trim(arg)
     if (n.gt.1) then
        call getarg(1,arg)
        read (arg,*) tempnum
        if (tempnum.eq.1) no_rewrite=.false.
        if (tempnum.ne.1) no_rewrite=.true.
     else
        no_rewrite=.true.
     end if
  end if
     
  ! Initialise number of files done
  ngs = 0; ndm = 0; nqq = 0
  
  ! Do Halos Extractor for gas+star halo files
  do ifile=1,maxfiles
     ! Generate file names for current iteration
     call title3c(ifile,file_num)
     write(filenameIn,'(a,a3)') 'tree_brick_starsub_',file_num
     write(filenameOut,'(a,a3)') 'clean_tree_brick_starsub_',file_num
     enquirerf=.false.
     do_file=.false.
     inquire(file=trim(execFolder)//trim(filenameIn), exist=enquiref)
     if (no_rewrite) inquire(file=trim(execFolder)//trim(filenameOut), exist=enquirerf)
     if ((enquiref.eq..true.).and.(enquirerf.eq..false.)) do_file=.true.
     ! Read and write halos
     if (do_file) then
        call check_working
        ngs = ngs + 1
        call read_tree_brick
        call write_tree_brick_clean
        deallocate(liste_halos,nb_of_parts)
     end if
  end do

  ! Do Halos Extractor for DM halo files
  do ifile=1,maxfiles
     ! Generate file names for current iteration
     call title3c(ifile,file_num)
     write(filenameIn,'(a,a3)') 'tree_bricks',file_num
     write(filenameOut,'(a,a3)') 'clean_tree_bricks',file_num
     enquirerf=.false.
     do_file=.false.
     inquire(file=trim(execFolder)//trim(filenameIn), exist=enquiref)
     if (no_rewrite) inquire(file=trim(execFolder)//trim(filenameOut), exist=enquirerf)
     if ((enquiref.eq..true.).and.(enquirerf.eq..false.)) do_file=.true.
     ! Read and write halos
     if (do_file) then
        call check_working
        ndm = ndm + 1
        call read_tree_brick
        call write_tree_brick_clean
        deallocate(liste_halos,nb_of_parts)
     end if
  end do

  ! Report results
  if (ngs.gt.0) write(*,'(": Number of gas+stars halos cleaned: ",I3)') ngs
  if (ndm.gt.0) write(*,'(": Number of dm halos cleaned       : ",I3)') ndm
  if (nqq.gt.0) write(*,'(": Number of other halos cleaned    : ",I3)') nqq
  if (nqq.eq.0 .and. ndm.eq.0 .and. ngs.eq.0) write(*,'(": No halos found in this folder.")')

  write(*,'(": Execution concluded, closing...")')
  
contains
!############################################################################
!############################################################################
!############################################################################
!############################################################################
  subroutine check_working
    if (working.eq..false.) then
       working=.true.
       write(*,'(": Be pacient: I have found some files and I am working...")')
    end if
  end subroutine check_working
!############################################################################
!############################################################################
!############################################################################
!############################################################################
  subroutine title3c(n,nchar)
    implicit none
    integer(kind=4) :: n
    character*3     :: nchar
    character*1     :: nchar1
    character*2     :: nchar2
    character*3     :: nchar3
    if(n.ge.100)then
       write(nchar3,'(i3)') n
       nchar = nchar3
    elseif(n.ge.10)then
       write(nchar2,'(i2)') n
       nchar = '0'//nchar2
    else
       write(nchar1,'(i1)') n
       nchar = '00'//nchar1
    endif
  end subroutine title3c
!############################################################################
!############################################################################
!############################################################################
!############################################################################
  subroutine read_tree_brick
    implicit none
    integer(kind=4)                                         :: i,unitfile,start,j  
    integer(kind=4) ,allocatable                            :: members(:)
    logical                                                 :: done
    done = .false.
    unitfile = 44
    open(unit=unitfile,file=trim(execFolder)//trim(filenameIn),form='unformatted',status='unknown')
    read(unitfile) nbodies
    read(unitfile) massp
    read(unitfile) aexp
    read(unitfile) omega_t
    read(unitfile) age_univ
    read(unitfile) nb_of_halos, nb_of_subhalos    
    allocate(nb_of_parts(1:nb_of_halos + nb_of_subhalos))
    allocate(liste_halos(1:nb_of_halos + nb_of_subhalos))
    do i=1,nb_of_halos + nb_of_subhalos
       read(unitfile) nb_of_parts(i)
       ! read list of particles in each halo
       allocate(members(1:nb_of_parts(i)))
       read(unitfile) members
       deallocate(members)
       ! read each halo properties
       call read_halo(liste_halos(i),unitfile)         
    end do
    close(unitfile)
    return
  end subroutine read_tree_brick
!############################################################################
!############################################################################
!############################################################################
!############################################################################
  subroutine write_tree_brick_clean
    implicit none
    integer(kind=4)                                         :: i,unitfile,start,j  
    integer(kind=4) ,allocatable                            :: members(:)
    logical                                                 :: done
    done = .false. 
    unitfile = 44
    open(unit=unitfile,file=trim(execFolder)//trim(filenameOut),form='unformatted',status='unknown')
    write(unitfile) nbodies
    write(unitfile) massp
    write(unitfile) aexp
    write(unitfile) omega_t
    write(unitfile) age_univ
    write(unitfile) nb_of_halos, nb_of_subhalos
    do i=1,nb_of_halos + nb_of_subhalos
       ! write list of particles in each halo
       write(unitfile) nb_of_parts(i)
       ! write each halo properties
       call write_halo(liste_halos(i),unitfile)         
    enddo
    close(unitfile)
    return
  end subroutine write_tree_brick_clean
!############################################################################
!############################################################################
!############################################################################
!############################################################################
  subroutine read_halo(h,unitfile)
    implicit none
    integer(kind=4) :: unitfile
    type (halo)     :: h
    read(unitfile) h%my_number
    read(unitfile) h%my_timestep 
    read(unitfile) h%level,h%hosthalo,h%hostsub,h%nbsub,h%nextsub
    read(unitfile) h%m
    read(unitfile) h%p%x,h%p%y,h%p%z
    read(unitfile) h%v%x,h%v%y,h%v%z
    read(unitfile) h%L%x,h%L%y,h%L%z 
    read(unitfile) h%r, h%sh%a, h%sh%b, h%sh%c
    read(unitfile) h%ek,h%ep,h%et
    read(unitfile) h%spin
    read(unitfile) h%datas%rvir,h%datas%mvir,h%datas%tvir,h%datas%cvel
    read(unitfile) h%halo_profile%rho_0,h%halo_profile%r_c
    return
  end subroutine read_halo
!############################################################################
!############################################################################
!############################################################################
!############################################################################
  subroutine write_halo(h,unitfile)
    implicit none
    integer(kind=4) :: unitfile
    type (halo)     :: h
    write(unitfile) h%my_number
    write(unitfile) h%my_timestep 
    write(unitfile) h%level,h%hosthalo,h%hostsub,h%nbsub,h%nextsub
    write(unitfile) h%m
    write(unitfile) h%p%x,h%p%y,h%p%z
    write(unitfile) h%v%x,h%v%y,h%v%z
    write(unitfile) h%L%x,h%L%y,h%L%z 
    write(unitfile) h%r, h%sh%a, h%sh%b, h%sh%c
    write(unitfile) h%ek,h%ep,h%et
    write(unitfile) h%spin
    write(unitfile) h%datas%rvir,h%datas%mvir,h%datas%tvir,h%datas%cvel
    write(unitfile) h%halo_profile%rho_0,h%halo_profile%r_c
    return
  end subroutine write_halo
!############################################################################
!############################################################################
!############################################################################
!############################################################################
end program HalosExtractor
