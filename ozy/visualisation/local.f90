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
public :: sgl,dbl,ish,irg

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


end module local
