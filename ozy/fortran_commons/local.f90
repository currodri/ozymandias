!--------------------------------------------------------------------------
! ozymandias:local.f90
!--------------------------------------------------------------------------
!
! MODULE: local
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> local definitions of single and double precision, general constants and variables
!
!> @details  
!> defines the kind-parameters for short and long integers, and single/double 
!> precision reals.
! 
!
!> @date 28/3/2021   0.1 fortran routines for ozymandias
!--------------------------------------------------------------------------
module local

!> @note This module must be "use"d by every program, subroutine, and function!

! The entire ozymandias package should be processor independent.  This can
! be accomplished by the use of the "kind" parameters.

implicit none
private
public :: sgl,dbl,ish,irg,ilg,rt_real_kind,hydro_real_kind

! Define the "kind" parameters for single and double precision reals, 
!> single precision real kind parameter
  integer,parameter                     :: sgl = SELECTED_REAL_KIND(p=6,r=37)   
!> double precision real kind parameter
  integer,parameter                     :: dbl = SELECTED_REAL_KIND(p=13,r=200) 

! Define the "kind" parameters for short and regular integers,
!> short integer kind parameter 
  integer,parameter                     :: ish = SELECTED_INT_KIND(3) 
!> long integer kind parameter  
  integer,parameter                     :: irg = SELECTED_INT_KIND(9)
!> long integer kind parameter  
  integer,parameter                     :: ilg = SELECTED_INT_KIND(12)

! Define the real kind used for radiation transfer arrays (rt_var).
! Use the preprocessor macro RTPRE (set in the build) to pick single or
! double precision for RT variables. Defining this alias here keeps
! conditional preprocessing out of individual routine argument lists,
! which helps tools like f90wrap/f2py parse consistent signatures.
#if RTPRE==4
  integer,parameter :: rt_real_kind = sgl
#warning "RT variables compiled with SINGLE precision (RTPRE=4)"
#elif RTPRE==8
  integer,parameter :: rt_real_kind = dbl
#warning "RT variables compiled with DOUBLE precision (RTPRE=8)"
#else
  integer,parameter :: rt_real_kind = dbl
#warning "RT variables compiled with DOUBLE precision (default)"
#endif

! Define the real kind used for hydro arrays (hydro_var).
! Use the preprocessor macro UPRE (set in the build) to pick single or
! double precision for hydro variables. Defining this alias here keeps
! conditional preprocessing out of individual routine argument lists,
! which helps tools like f90wrap/f2py parse consistent signatures.
#if UPRE==4
  integer,parameter :: hydro_real_kind = sgl
#warning "Hydro variables compiled with SINGLE precision (UPRE=4)"
#elif UPRE==8
  integer,parameter :: hydro_real_kind = dbl
#warning "Hydro variables compiled with DOUBLE precision (UPRE=8)"
#else
  integer,parameter :: hydro_real_kind = dbl
#warning "Hydro variables compiled with DOUBLE precision (default)"
#endif

end module local

module constants
  use local
  implicit none
  ! Other constants
  real(dbl) :: pi=3.14159265359
  real(dbl) :: twopi = 6.283185307179586
  real(dbl) :: halfsqrt2 = 0.7071067811865476
  real(dbl) :: smallc = 1D-10
  real(dbl) :: smallr = 1D-10
  ! Mass
  real(dbl),parameter::g2msun=5.02739933e-034
  ! Density
  real(dbl),parameter::gcm32msunpc3=1.47755759e+22
  ! SFR
  real(dbl),parameter::gs2msunyr=1.585e-26
  real(dbl),parameter::gscm22msunyrkpc2=1.51112d17
  ! Distance
  real(dbl),parameter::cm2km=1d-5
  real(dbl),parameter::cm2pc=3.24078e-19
  real(dbl),parameter::cm2kpc=3.2408d-22
  real(dbl),parameter::cm2Mpc=3.2408d-25
  real(dbl),parameter::pc2cm=1/cm2pc
  real(dbl),parameter::Armstrong2cm=1d-8
  ! Time
  real(dbl),parameter::s2Gyr=3.1709e-17
  real(dbl),parameter::s2Myr=3.1709e-14
  real(dbl),parameter::s2kyr=3.1709e-11
  real(dbl),parameter::s2yr=3.1709e-8
  real(dbl),parameter::Myr2s=3.1536e13
  ! Energy
  real(dbl),parameter::eV2erg=1.6022e-12
  ! Fundamental constants ##################################################
  ! Radiation constants
  real(dbl),parameter::alphaB=2.6d-13 ! Case B recombination rate (cm^3 s^−1; Ferland 1992)
  ! (Physical) constants
  real(dbl),parameter::Ggrav=6.674d-8     ! Gravitational constant (cm^3 g^-1 s^-2)
  real(dbl),parameter::clight=29979245800.0d0      ! Speed of light (cm/s)
  real(dbl),parameter::m_electron=9.1093837015d-28 ! Electron mass (g)
  real(dbl),parameter :: mHydrogen=1.673532784796145D-24 ! Mass of hydrogen atom in g
  real(dbl),parameter :: kBoltzmann=1.380649D-16 ! Boltzmann constant in erg/K
  real(dbl),parameter :: gas_R=8.31446261815324D7 ! Ideal gas constant R in erg/(k*mol)
  real(dbl),parameter :: atoweightH=1.00784 ! Atomic weight of Hydrogen in atomic mass units
  real(dbl),parameter :: cVHydrogen=1.4D8 !123746.76463754028 ! Specific heat capacity at constat volume, in erg/(K*g)
  real(dbl),parameter :: XH=0.76 ! Hydrogen fraction
  real(dbl),parameter :: YHe=0.24 ! Helium fraction
  real(dbl),parameter :: mu=0.5882352941176471 ! Primordial, fully ionised mean molecular weight
  real(dbl),parameter :: gamma_gas=1.6666667   ! Always assuming monatomic adiabatic gas
  real(dbl),parameter :: gamma_cr=1.3333333   ! Cosmic rays assumed always relativistic
  ! Elemental constants
  real(dbl),parameter :: amu2g = 1.66054d-24       ! Atomic mass units in grams
  real(dbl),parameter :: mH_amu = 1.007825d0       ! Hydrogen molecular weight [amu]
  real(dbl),parameter :: mHe_amu = 4.002602d0      ! Helium molecular weight [amu]
  real(dbl),parameter :: mO_amu = 15.9994d0        ! Oxygen molecular weight [amu]
  real(dbl),parameter :: mC_amu = 12.0107d0        ! Carbon molecular weight [amu]
  real(dbl),parameter :: mMg_amu = 24.305d0        ! Magnesium molecular weight [amu]
  real(dbl),parameter :: mSi_amu = 28.0855d0       ! Silicon molecular weight [amu]
  real(dbl),parameter :: mFe_amu = 55.854d0        ! Iron molecular weight [amu]
  real(dbl),parameter :: mOovermH = 15.875d0       ! Oxygen to Hydrogen mass ratio [dimensionless]
  ! Molecular constants
  real(dbl),parameter :: CO_to_Cmass = 0.4288      ! CO to C mass [dimensionless]
  ! Cosmology
  real(dbl),parameter :: rhoc=1.8800000d-29   ! Critical density of the Universe (g cm-3)
  ! Milky Way constants
  real(dbl),parameter :: ecr_sun=1.4d-12 ! CR energy density in the Solar neighbourhood (Boschini et al. 2020, erg/cm-3)
  real(dbl),parameter :: LgammaH=1.1d-28 ! Gamma-ray luminosity per H atom in the 0.5-5 GeV range (Casandjian 2015)
  ! Constants for Yohan's version of dust simulations
  real(dbl) :: SiOverSil=0.163d0
end module constants