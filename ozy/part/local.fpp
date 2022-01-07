# 1 "local.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "local.f90"
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
public :: sgl,dbl,ish,irg,ilg

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


end module local

module constants
  use local
  implicit none
  public :: mHydrogen,kBoltzmann,gas_R,atoweightH,cVHydrogen
  ! Physical constants
  real(dbl) :: mHydrogen=1.673532784796145D-24 ! Mass of hydrogen atom in g
  real(dbl) :: kBoltzmann=1.380649D-16 ! Boltzmann constant in erg/K
  real(dbl) :: gas_R=8.31446261815324D7 ! Ideal gas constant R in erg/(k*mol)
  real(dbl) :: atoweightH=1.00784 ! Atomic weight of Hydrogen in atomic mass units
  real(dbl) :: cVHydrogen=1.4D8 !123746.76463754028 ! Specific heat capacity at constat volume, in erg/(K*g)
end module constants
