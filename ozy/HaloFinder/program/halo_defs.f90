module halo_defs

  public

  ! precision for particle positions
#ifdef STARS
  integer(kind=4),parameter :: dp=kind(1d0)
#else
  integer(kind=4),parameter :: dp=kind(1e0)
#endif

  ! I/O units
  integer(kind=4),parameter :: haloPropsUnit = 45


  !======================================================================
  ! useful types :
  !======================================================================
  type vector
     real(dp)      :: x,y,z
  end type vector

  type shape 
     real(dp)      :: a,b,c
  end type shape

  type baryon
     real(dp)      :: rvir,mvir,tvir,cvel
  end type baryon

  type hprofile
     real(dp)      :: rho_0,r_c
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
     real(dp)      :: m
     real(dp)      :: r
     real(dp)      :: spin
     real(dp)      :: ek,ep,et
     ! jeje 
     integer(kind=4)   :: minPartID   ! smallest ID of particles in the (sub-)halo
     ! end jeje 
#ifdef CONTAM
     integer(kind=4)   :: contaminated ! 0 if not, 1 if yes... 
#endif
  end type halo
  !======================================================================


  !======================================================================
  ! parameters relative to the simulation analysis
  !======================================================================
  character(len=8),parameter :: gravsoft  = 'cubsplin'         ! type of gravitational softening 
  integer(kind=4)            :: nbPes                          ! obsolete vars for reading treecode format
  integer(kind=4)            :: nsteps                         ! number of timesteps
  real(dp)               :: alpha,tnow                     ! 
  integer(kind=4)            :: nbodies                        ! number of particles in the simulation 
  integer(kind=4)            :: nMembers                       ! minimal number of particles in a fof halo
  real(dp)               :: b_init                         ! linking length parameter of the fof at z=0
  character(len=4),parameter :: profile   = 'TSIS'             ! type of halo profile (only isothermal_sphere yet)
  integer(kind=4),parameter  :: ninterp   = 30000              ! nb of bins for interp. of smoothed grav. field
  integer(kind=4)            :: FlagPeriod                     ! flag for periodicity of boundary conditions 
  real(dp)               :: fPeriod(3)
  !----------- For gadget format: ----------------------------------------
  integer (kind=4)           :: nhr                            ! to read only the selected HR particles
  real (kind=4)              :: minlrmrat                      ! to recognize contaminated halos
  !======================================================================


  !======================================================================
  ! Definitions specific to input/output
  !======================================================================
  character(len=80)         :: data_dir
  character(len=3)          :: file_num
  integer(kind=4)           :: numstep
  integer(kind=4),parameter :: errunit = 6
  logical(kind=4)           :: write_resim_masses              ! for writing resim_masses.dat file
  !======================================================================


  !======================================================================
  ! Constants
  !======================================================================
  real(dp),parameter    :: gravconst = 430.1d0              ! G in units of (km/s)^2 Mpc/(10^11 Msol)
  real(dp),parameter    :: pi        = 3.141592654d0
  !======================================================================


  !======================================================================
  ! Global variables 
  !======================================================================
  real(dp),allocatable             :: pos(:,:),vel(:,:)
  real(dp)                     :: massp
  real(dp),allocatable         :: epsvect(:),mass(:)
  ! array of particle fields
  ! (e.g. (/'pos','vel','mass','iord','level'/) for a DM-only run):
  character(30) :: ParticleFields(20)  
#ifdef STARS
  real(dp),allocatable         :: stellarAge(:),stellarZ(:),massInit(:)
#endif
  real(dp)                     :: omega_t,omega_lambda_t,omega_f,omega_b_f,omega_lambda_f,omega_c_f
  real(dp)                     :: rho_crit,aexp,Lboxp,mboxp,af,ai,Lf,H_f,H_i
  real(dp)                     :: age_univ,Lbox_pt,Lbox_pt2,Hub_pt,omega_0,hubble,omega_lambda_0
  real(dp)                     :: vir_overdens,rho_mean
  integer(kind=4),allocatable      :: linked_list(:), liste_parts(:)
  integer(kind=4),allocatable      :: first_part(:),nb_of_parts(:)
  type (halo),allocatable          :: liste_halos(:)
  real(dp)                     :: phsmooth(0:1+ninterp)
  integer(kind=4)                  :: nb_of_halos, nb_of_subhalos
  integer(kind=4)                  :: numero_step  
  character(len=3)                 :: type
  ! leo
  real(dp)                     :: omega_z,delta_crit_z,rho_crit_z
  ! end leo
  !======================================================================


  !======================================================================
  ! defs for Adaptahop
  !======================================================================
   integer(kind=4), parameter :: nparmax=512*512*512
   integer(kind=4),parameter  :: nparbuffer=128**3
   integer(kind=8),parameter  :: ncellbuffermin=128**3
   integer(kind=4),parameter  :: nlevelmax=30,nbodiespercell=20
   integer(kind=4)            :: lin=10,lsin=11,lin2=12
   integer(kind=4)            :: loudis=12,lounei=13,lounode=14,loupartnode=15,lounodedyn=16
   real(dp)            :: bignum=1.e30
   ! Physical constants (units : m s kg) ->
   real(kind=8), parameter :: gravitational_constant=6.6726d-11
   real(kind=8), parameter :: critical_density= 1.8788d-26
   real(kind=8), parameter :: mega_parsec=3.0857d22
   real(kind=8), parameter :: solar_mass=1.989d30  ! in kg
   real(kind=8), parameter :: convert_in_mps2=1.d6
   real(dp), allocatable    :: vxout(:),vyout(:),vzout(:),vdisout(:)
   integer(kind=4), allocatable :: mass_cell(:)
   real(kind=8), allocatable    :: tmass_cell(:)
   real(kind=8), allocatable    :: vcell(:,:)
   real(dp), allocatable    :: size_cell(:)
   real(dp), allocatable    :: pos_cell(:,:)
   integer(kind=4), allocatable :: sister(:)
   integer(kind=4), allocatable :: firstchild(:)
   integer(kind=4), allocatable :: idpart(:),idpart_tmp(:)
   integer(kind=4), allocatable :: iparneigh(:,:)
   real(kind=8), allocatable    :: distance(:)
   real(dp), allocatable    :: density(:)
   integer(kind=4), allocatable :: firstpart(:)
   integer(kind=4), allocatable :: igrouppart(:)
   integer(kind=4), allocatable :: idgroup(:),idgroup_tmp(:)
   integer(kind=4), allocatable :: igroupid(:)
   integer(kind=4), allocatable :: color(:)
   integer(kind=4), allocatable :: partnode(:)
   real(dp), allocatable    :: densityg(:)
   real(dp)    :: sizeroot
   real(kind=8)    :: xlong, ylong, zlong, boxsize, boxsize2
   real(kind=8)    :: xlongs2, ylongs2, zlongs2
   real(dp)    :: omega0,omegaL,GMphys
   real(dp)    :: aexp_max
   integer(kind=4) :: nvoisins,nhop,ntype
   integer(kind=8) :: ncellmx,ncellbuffer
   integer(kind=4) :: ngroups,nmembthresh,nnodes,nnodesmax
   integer(kind=4) :: ncpu,nmpi,niterations
   real(dp)    :: rho_threshold
   logical         :: verbose,megaverbose,periodic
   real(kind=8)    :: fgas
   real(kind=8)    :: fudge,alphap,epsilon,fudgepsilon
   real(kind=8)    :: pos_shift(3),pos_renorm,velrenorm

   type grp
      sequence
      integer(kind=4) :: nhnei
      integer(kind=4) :: njunk ! To avoid missalignement in memory
      integer(kind=4), dimension(:),pointer :: isad_gr(:)
      real(dp), dimension(:),pointer    :: rho_saddle_gr(:)
   end type grp

   type supernode
      sequence
      integer(kind=4) :: level
      integer(kind=4) :: mother
      integer(kind=4) :: firstchild
      integer(kind=4) :: nsisters
      integer(kind=4) :: sister
      real(dp)    :: rho_saddle
      real(dp)    :: density
      real(dp)    :: densmax
      real(dp)    :: radius
      integer(kind=4) :: mass
      real(dp)    :: truemass
      real(dp)    :: position(3)
   end type supernode

   type (grp), allocatable       :: group(:)
   type (supernode), allocatable :: node(:)

!======================================================================
! Flags for halo finder selection
!======================================================================
   character(len=3)    :: method       ! flag to notify which and how the halofinder is to be used
   logical             :: fsub         ! flag to notify whether subhaloes are included               
   logical             :: cdm          ! flag to select particle closest to the cdm instead of the one with the highest density
   logical             :: ssm          ! flag to activate shrinking spheres method instead of cdm or highest density 
!======================================================================
! array to build the structure tree
!======================================================================
integer(kind=4), allocatable :: first_daughter(:), mother(:), first_sister(:), level(:)
integer(kind=4) :: nstruct

! used for the merger history method
integer(kind=4), allocatable :: npfather(:),ex_liste_parts(:),removesub(:)
integer(kind=4), allocatable :: ex_level(:),ex_nb_of_parts(:)
integer(kind=4) ::  ex_nb_of_structs


contains

!///////////////////////////////////////////////////////////////////////
!***********************************************************************
  subroutine clear_halo(h)

    implicit none

    type (halo)     :: h

    h%my_number            = 0
    h%my_timestep          = 0
    h%nbsub                = 0
    h%hosthalo             = 0
    h%hostsub              = 0
    h%level                = 1
    h%nextsub              = -1
    h%m                    = 0.0
    h%p%x                  = 0.0
    h%p%y                  = 0.0    
    h%p%z                  = 0.0  
    h%v%x                  = 0.0
    h%v%y                  = 0.0
    h%v%z                  = 0.0
    h%L%x                  = 0.0
    h%L%y                  = 0.0
    h%L%z                  = 0.0
    h%spin                 = 0.0
    h%r                    = 0.0
    h%sh%a                 = 0.0
    h%sh%b                 = 0.0
    h%sh%c                 = 0.0
    h%ek                   = 0.0
    h%ep                   = 0.0
    h%et                   = 0.0
    h%datas%rvir           = 0.0
    h%datas%mvir           = 0.0  
    h%datas%tvir           = 0.0
    h%datas%cvel           = 0.0
    h%halo_profile%rho_0   = 0.0
    h%halo_profile%r_c     = 0.0
    ! jeje 
    h%minPartID            = -1
#ifdef CONTAM 
    h%contaminated         = -1
#endif

    return

  end subroutine clear_halo

!***********************************************************************
! jeje 
  
  subroutine copy_halo(h1,h2)
    
    ! copy h1 into h2
    
    implicit none

    type(halo) :: h1,h2

    h2%my_number            = h1%my_number         
    h2%my_timestep          = h1%my_timestep       
    h2%nbsub                = h1%nbsub             
    h2%hosthalo             = h1%hosthalo          
    h2%hostsub              = h1%hostsub           
    h2%level                = h1%level             
    h2%nextsub              = h1%nextsub           
    h2%m                    = h1%m                 
    h2%p%x                  = h1%p%x               
    h2%p%y                  = h1%p%y                   
    h2%p%z                  = h1%p%z                 
    h2%v%x                  = h1%v%x               
    h2%v%y                  = h1%v%y               
    h2%v%z                  = h1%v%z               
    h2%L%x                  = h1%L%x               
    h2%L%y                  = h1%L%y               
    h2%L%z                  = h1%L%z               
    h2%spin                 = h1%spin              
    h2%r                    = h1%r                 
    h2%sh%a                 = h1%sh%a              
    h2%sh%b                 = h1%sh%b              
    h2%sh%c                 = h1%sh%c              
    h2%ek                   = h1%ek                
    h2%ep                   = h1%ep                
    h2%et                   = h1%et                
    h2%datas%rvir           = h1%datas%rvir        
    h2%datas%mvir           = h1%datas%mvir          
    h2%datas%tvir           = h1%datas%tvir        
    h2%datas%cvel           = h1%datas%cvel        
    h2%halo_profile%rho_0   = h1%halo_profile%rho_0
    h2%halo_profile%r_c     = h1%halo_profile%r_c  
    ! jeje 
    h2%minPartID = h1%minPartID
#ifdef CONTAM
    h2%contaminated         = h1%contaminated
#endif

    return
     
  end subroutine copy_halo
! end jeje 
!***********************************************************************
!///////////////////////////////////////////////////////////////////////

end module halo_defs


  
