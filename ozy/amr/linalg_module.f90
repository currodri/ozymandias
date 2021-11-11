!--------------------------------------------------------------------------
! ozymandias:linalg_module.f90
!--------------------------------------------------------------------------
!
! MODULE: vectors
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> types, functions, operators and routines useful for linear algebra
!
!> @details  
!> defines the derived type vector, with its attributes and bound
!> procedures
! 
!
!> @date 28/3/2021   0.1 fortran routines for ozymandias
!--------------------------------------------------------------------------
module vectors
    use local, only: sgl, dbl
    
    implicit none

    type vector
        real(dbl) :: x = 0d0
        real(dbl) :: y = 0d0
        real(dbl) :: z = 0d0
    end type vector

    type array_vectors
        integer :: n
        type(vector),dimension(:),allocatable :: list
    end type array_vectors

    interface assignment (=)
        module procedure array_to_vector
        module procedure vector_to_array
    end interface

    interface operator (+)
        module procedure vector_add
    end interface

    interface operator (-)
        module procedure vector_subtract
    end interface

    interface operator (*)
        module procedure vector_times_real
        module procedure real_times_vector
        module procedure vector_times_int
        module procedure int_times_vector
        module procedure cross_product
        module procedure matrix_times_vector
    end interface

    interface operator (/)
        module procedure vector_div_real
        module procedure vector_div_int
    end interface

    interface operator (.DOT.)
        module procedure my_dot_product
    end interface

    interface operator (.eq.)
        module procedure compare_vectors
    end interface

    contains
    subroutine array_to_vector(vec_result,array)
        type(vector), intent(inout) :: vec_result
        real(dbl), dimension(3), intent(in) :: array
        vec_result%x = array(1)
        vec_result%y = array(2)
        vec_result%z = array(3)
    end subroutine array_to_vector

    subroutine vector_to_array(array_result,vec_1)
        type(vector), intent(in) :: vec_1
        real(dbl), dimension(3), intent(inout) :: array_result
        array_result(1) = vec_1%x
        array_result(2) = vec_1%y
        array_result(3) = vec_1%z
    end subroutine vector_to_array

    type(vector) function vector_add(vec_1,vec_2)
        type(vector), intent(in) :: vec_1,vec_2
        vector_add%x = vec_1%x + vec_2%x
        vector_add%y = vec_1%y + vec_2%y
        vector_add%z = vec_1%z + vec_2%Z
    end function vector_add

    type(vector) function vector_subtract(vec_1,vec_2)
        type(vector), intent(in) :: vec_1,vec_2
        vector_subtract%x = vec_1%x - vec_2%x
        vector_subtract%y = vec_1%y - vec_2%y
        vector_subtract%z = vec_1%z - vec_2%Z
    end function vector_subtract

    type(vector) function vector_times_real(vec_1,real_2)
        type(vector), intent(in) :: vec_1
        real(dbl), intent(in) :: real_2
        vector_times_real%x = vec_1%x * real_2
        vector_times_real%y = vec_1%y * real_2
        vector_times_real%z = vec_1%z * real_2
    end function vector_times_real

    type(vector) function real_times_vector(real_1,vec_2)
        type(vector), intent(in) :: vec_2
        real(dbl), intent(in) :: real_1
        real_times_vector%x = vec_2%x * real_1
        real_times_vector%y = vec_2%y * real_1
        real_times_vector%z = vec_2%z * real_1
    end function real_times_vector

    type(vector) function vector_times_int(vec_1,int_2)
        type(vector), intent(in) :: vec_1
        integer, intent(in) :: int_2
        vector_times_int%x = vec_1%x * real(int_2)
        vector_times_int%y = vec_1%y * real(int_2)
        vector_times_int%z = vec_1%z * real(int_2)
    end function vector_times_int

    type(vector) function int_times_vector(int_1,vec_2)
        type(vector), intent(in) :: vec_2
        integer, intent(in) :: int_1
        int_times_vector%x = vec_2%x * real(int_1)
        int_times_vector%y = vec_2%y * real(int_1)
        int_times_vector%z = vec_2%z * real(int_1)
    end function int_times_vector

    type(vector) function vector_div_real(vec_1,real_2)
        type(vector), intent(in) :: vec_1
        real(dbl), intent(in) :: real_2
        vector_div_real%x = vec_1%x / real_2
        vector_div_real%y = vec_1%y / real_2
        vector_div_real%z = vec_1%z / real_2
    end function vector_div_real

    type(vector) function vector_div_int(vec_1,int_2)
        type(vector), intent(in) :: vec_1
        integer, intent(in) :: int_2
        vector_div_int%x = vec_1%x / real(int_2)
        vector_div_int%y = vec_1%y / real(int_2)
        vector_div_int%z = vec_1%z / real(int_2)
    end function vector_div_int 

    real(dbl) function my_dot_product(vec_1,vec_2)
        type(vector), intent(in) :: vec_1,vec_2
        my_dot_product = vec_1%x*vec_2%x + vec_1%y*vec_2%y &
                        + vec_1%z*vec_2%z
    end function my_dot_product

    type(vector) function cross_product(vec_1,vec_2)
        type(vector), intent(in) :: vec_1,vec_2
        cross_product%x = vec_1%y*vec_2%z - vec_1%z*vec_2%y
        cross_product%y = vec_1%z*vec_2%x - vec_1%x*vec_2%z
        cross_product%z = vec_1%x*vec_2%y - vec_1%y*vec_2%x
    end function cross_product

    type(vector) function matrix_times_vector(matrix_1,vec_2)
        real(dbl), dimension(1:3,1:3),intent(in) :: matrix_1
        type(vector), intent(in) :: vec_2
        real(dbl),dimension(1:3) :: array_2
        array_2 = vec_2
        matrix_times_vector = matmul(matrix_1,array_2)
    end function matrix_times_vector

    real(dbl) function magnitude(vec_1)
        type(vector), intent(in) :: vec_1
        magnitude = sqrt(vec_1.DOT.vec_1)
    end function magnitude

    logical function compare_vectors(vec_1,vec_2)
        type(vector), intent(in) :: vec_1,vec_2
        compare_vectors = (vec_1%x.eq.vec_2%x) .and.&
                            &(vec_1%y.eq.vec_2%y).and.&
                            &(vec_1%z.eq.vec_2%z)
    end function compare_vectors
end module vectors
!--------------------------------------------------------------------------
!
! MODULE: rotations
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> subroutines useful for rotations
!
!> @details  
!> define basic rotation matrices and how they are applied to vectors
! 
!
!> @date 28/3/2021   0.1 fortran routines for ozymandias
!--------------------------------------------------------------------------
module rotations
    use local
    use vectors
    implicit none
    interface rotate_vector
        module procedure rotate_vector_single
        module procedure rotate_vector_array
    end interface rotate_vector
    contains
    
    !---------------------------------------------------------------
    ! Subroutine: EULER MATRIX
    !
    ! Computes the rotation matrix for a particular cartesian axis
    ! and a given angle of rotation.
    ! Options for passive or active rotation are especified in the
    ! argument `ap'.
    !---------------------------------------------------------------
    subroutine euler_matrix(R,dim,angle,ap)
        implicit none
        integer, intent(in) :: dim
        real(dbl), intent(in) :: angle
        character(1), intent(in) :: ap
        real(dbl), dimension(1:3,1:3),intent(inout) :: R
        real(sgl),parameter :: thr = 1.0E-8
        integer(irg) :: i,j
        select case (dim)
        case (1)
            R(1,1) = 1; R(1,2) = 0; R(1,3) = 0
            R(2,1) = 0; R(2,2) = dcos(angle); R(2,3) = -dsin(angle)
            R(3,1) = 0; R(3,2) = dsin(angle); R(3,3) = dcos(angle)
        case(2)
            R(1,1) = dcos(angle); R(1,2) = 0; R(1,3) = dsin(angle)
            R(2,1) = 0; R(2,2) = 1; R(2,3) = 0
            R(3,1) = -dsin(angle); R(3,2) = 0; R(3,3) = dcos(angle)
        case(3)
            R(1,1) = dcos(angle); R(1,2) = -dsin(angle); R(1,3) = 0
            R(2,1) = dsin(angle); R(2,2) = dcos(angle); R(2,3) = 0
            R(3,1) = 0; R(3,2) = 0; R(3,3) = 1
        case default
            write(*,*)'dim must be 1,2 or 3.'
            write(*,*)'Aborting!'
            stop
        end select
        do i=1,3
            do j=1,3
                if (abs(R(i,j)).lt.thr) R(i,j) = 0.0
            end do
        end do
        if (ap.ne.'p') then
            R= transpose(R)
        end if
    end subroutine euler_matrix

    !---------------------------------------------------------------
    ! Subroutine: ROTATE VECTOR
    !
    ! Rotates vector derived types by matrix multiplication with
    ! a given rotation matrix.
    ! This is built as an interface, such that single and array of
    ! vectors can be multiplied with the single generic procedure
    ! `rotate_vector'.
    !---------------------------------------------------------------
    subroutine rotate_vector_single(vec,rotation_matrix)
        type(vector),intent(inout) :: vec
        real(dbl),dimension(1:3,1:3),intent(in) :: rotation_matrix
        vec = rotation_matrix * vec
    end subroutine rotate_vector_single

    subroutine rotate_vector_array(vec,rotation_matrix)
        type(array_vectors),intent(inout) :: vec
        real(dbl),dimension(1:3,1:3),intent(in) :: rotation_matrix
        integer :: i
        do i=1,vec%n
            vec%list(i) = rotation_matrix * vec%list(i)
        end do
    end subroutine rotate_vector_array
end module rotations
!--------------------------------------------------------------------------
!
! MODULE: basis_representations
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> types, functions and routines useful for basis representations
!
!> @details  
!> defines the derived type basis and how a orthonormal one can be obtained
! 
!
!> @date 28/3/2021   0.1 fortran routines for ozymandias
!--------------------------------------------------------------------------
module basis_representations
    use local
    use vectors

    type basis
        type(vector) :: u(3)
    end type

    contains
    subroutine initialise_basis(this)
        type(basis),intent(inout) :: this
        this%u(1)%x = 1D0
        this%u(2)%y = 1D0
        this%u(3)%z = 1D0
    end subroutine initialise_basis

    !---------------------------------------------------------------
    ! Subroutine: MODIFIED GRAM-SCHMIDT METHOD
    !
    ! This procedure computes an orthonormal basis for a provided 
    ! basis derived type using the numerically modified Gram-Schmidt
    ! method.
    !---------------------------------------------------------------
    subroutine mGramSchmidt(vecs,e)
        implicit none
        type(basis), intent(in) :: vecs
        type(basis), intent(inout) :: e
        integer :: i,j
        e%u(1) = vecs%u(1) / magnitude(vecs%u(1))
        do i=2,3
            e%u(i) = vecs%u(i)
            do j=1,i-1
                e%u(i) = e%u(i) - (e%u(i).DOT.e%u(j))/magnitude(e%u(j)) * e%u(j)
            end do
            e%u(i) = e%u(i) / magnitude(e%u(i))
        end do
    end subroutine mGramSchmidt
end module basis_representations

