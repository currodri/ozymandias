!--------------------------------------------------------------------------
! ozymandias:dictionary_module.f90
!--------------------------------------------------------------------------
!
! MODULE: dictionary_commons
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> types, functions, operators and routines to construct Python-like
!> dictionaries
!
!> @details  
!> defines the derived type dictf90
! 
!
!> @date 08/5/2023   0.3 making ozymandias fully compatible with non-cosmo
!--------------------------------------------------------------------------

module dictionary_commons
    use local

    implicit none

    private
    public :: dictf90

    type dictf90
        character(128), dimension(:), allocatable :: keys
        integer, dimension(:), allocatable            :: values
        integer            :: count = 0

        contains
            procedure :: init => allocate_dict_class
            procedure :: add => add_to_dict_class
            procedure :: get => find_key_class
    end type dictf90

    contains

    subroutine allocate_dict(self,nkeys)
        implicit none

        type(dictf90), intent(inout) :: self
        integer, intent(in)          :: nkeys

        if (.not.allocated(self%keys)) allocate(self%keys(nkeys))
        if (.not.allocated(self%values)) allocate(self%values(nkeys))
        self%count = 0
    end subroutine allocate_dict

    subroutine allocate_dict_class(self,nkeys)
        implicit none

        class(dictf90), intent(inout) :: self
        integer, intent(in)          :: nkeys

        call allocate_dict(self,nkeys)
    end subroutine allocate_dict_class

    subroutine add_to_dict(self, key, value)
        implicit none

        type(dictf90), intent(inout) :: self
        character(len=*), intent(in) :: key
        integer, intent(in)          :: value

        integer                      :: i
        
        ! Check if the key already exists in the dictionary
        do i = 1, self%count
            if (trim(self%keys(i)) == trim(key)) then
                self%values(i) = value
                return
            end if
        end do

        ! If the key does not exist, add it to the dictionary
        self%count              = self%count + 1
        self%keys(self%count)   = key
        self%values(self%count) = value

    end subroutine add_to_dict

    subroutine add_to_dict_class(self, key, value)
        implicit none

        class(dictf90), intent(inout) :: self
        character(len=*), intent(in) :: key
        integer, intent(in)          :: value

        call add_to_dict(self,key,value)
    end subroutine add_to_dict_class

    integer function find_key(self,key) 
        implicit none

        type(dictf90), intent(in)  :: self
        character(len=*), intent(in) :: key

        integer :: i

        do i = 1, self%count
            if (trim(self%keys(i)) == trim(key)) then
                find_key = self%values(i)
                return
            end if
        end do

        ! If key not found, return 0 index
        find_key = 0

    end function find_key

    integer function find_key_class(self,key) 
        implicit none

        class(dictf90), intent(in)  :: self
        character(len=*), intent(in) :: key

        find_key_class = find_key(self,key)
    end function find_key_class

end module dictionary_commons