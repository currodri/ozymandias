!--------------------------------------------------------------------------
! ozymandias:variables.f90
!--------------------------------------------------------------------------
!
! MODULE: variables
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> types and functions that allow the handling of hydro raw variables and
!> derived variables
!
!> @details  
!> 
! 
!
!> @date 08/5/2023   0.3 making ozymandias fully compatible with non-cosmo
!--------------------------------------------------------------------------
module hydro_commons
    use local
    use dictionary_commons

    type hydro_var
        character(128) :: name, rname
        integer, dimension(:), allocatable :: ids
    end type hydro_var

    contains

    function cleaned_vars(vardict,varlist)
        implicit none

        
    end function cleaned_vars

end module hydro_commons