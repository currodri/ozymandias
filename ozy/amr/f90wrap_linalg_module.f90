! Module vectors defined in file linalg_module.fpp

subroutine f90wrap_vector__get__x(this, f90wrap_x)
    use vectors, only: vector
    implicit none
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(vector_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_x
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_x = this_ptr%p%x
end subroutine f90wrap_vector__get__x

subroutine f90wrap_vector__set__x(this, f90wrap_x)
    use vectors, only: vector
    implicit none
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(vector_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_x
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%x = f90wrap_x
end subroutine f90wrap_vector__set__x

subroutine f90wrap_vector__get__y(this, f90wrap_y)
    use vectors, only: vector
    implicit none
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(vector_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_y
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_y = this_ptr%p%y
end subroutine f90wrap_vector__get__y

subroutine f90wrap_vector__set__y(this, f90wrap_y)
    use vectors, only: vector
    implicit none
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(vector_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_y
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%y = f90wrap_y
end subroutine f90wrap_vector__set__y

subroutine f90wrap_vector__get__z(this, f90wrap_z)
    use vectors, only: vector
    implicit none
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(vector_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_z
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_z = this_ptr%p%z
end subroutine f90wrap_vector__get__z

subroutine f90wrap_vector__set__z(this, f90wrap_z)
    use vectors, only: vector
    implicit none
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in)   :: this(2)
    type(vector_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_z
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%z = f90wrap_z
end subroutine f90wrap_vector__set__z

subroutine f90wrap_vector_initialise(this)
    use vectors, only: vector
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type(vector_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_vector_initialise

subroutine f90wrap_vector_finalise(this)
    use vectors, only: vector
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type(vector_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_vector_finalise

subroutine f90wrap_array_vectors__get__n(this, f90wrap_n)
    use vectors, only: array_vectors
    implicit none
    type array_vectors_ptr_type
        type(array_vectors), pointer :: p => NULL()
    end type array_vectors_ptr_type
    integer, intent(in)   :: this(2)
    type(array_vectors_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_n
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_n = this_ptr%p%n
end subroutine f90wrap_array_vectors__get__n

subroutine f90wrap_array_vectors__set__n(this, f90wrap_n)
    use vectors, only: array_vectors
    implicit none
    type array_vectors_ptr_type
        type(array_vectors), pointer :: p => NULL()
    end type array_vectors_ptr_type
    integer, intent(in)   :: this(2)
    type(array_vectors_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_n
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%n = f90wrap_n
end subroutine f90wrap_array_vectors__set__n

subroutine f90wrap_array_vectors__array_getitem__list(f90wrap_this, f90wrap_i, listitem)
    
    use vectors, only: array_vectors, vector
    implicit none
    
    type array_vectors_ptr_type
        type(array_vectors), pointer :: p => NULL()
    end type array_vectors_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(array_vectors_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: listitem(2)
    type(vector_ptr_type) :: list_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%list)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%list)) then
            call f90wrap_abort("array index out of range")
        else
            list_ptr%p => this_ptr%p%list(f90wrap_i)
            listitem = transfer(list_ptr,listitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_array_vectors__array_getitem__list

subroutine f90wrap_array_vectors__array_setitem__list(f90wrap_this, f90wrap_i, listitem)
    
    use vectors, only: array_vectors, vector
    implicit none
    
    type array_vectors_ptr_type
        type(array_vectors), pointer :: p => NULL()
    end type array_vectors_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(array_vectors_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: listitem(2)
    type(vector_ptr_type) :: list_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%list)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%list)) then
            call f90wrap_abort("array index out of range")
        else
            list_ptr = transfer(listitem,list_ptr)
            this_ptr%p%list(f90wrap_i) = list_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_array_vectors__array_setitem__list

subroutine f90wrap_array_vectors__array_len__list(f90wrap_this, f90wrap_n)
    
    use vectors, only: array_vectors, vector
    implicit none
    
    type array_vectors_ptr_type
        type(array_vectors), pointer :: p => NULL()
    end type array_vectors_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(array_vectors_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%list)) then
        f90wrap_n = size(this_ptr%p%list)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_array_vectors__array_len__list

subroutine f90wrap_array_vectors_initialise(this)
    use vectors, only: array_vectors
    implicit none
    
    type array_vectors_ptr_type
        type(array_vectors), pointer :: p => NULL()
    end type array_vectors_ptr_type
    type(array_vectors_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_array_vectors_initialise

subroutine f90wrap_array_vectors_finalise(this)
    use vectors, only: array_vectors
    implicit none
    
    type array_vectors_ptr_type
        type(array_vectors), pointer :: p => NULL()
    end type array_vectors_ptr_type
    type(array_vectors_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_array_vectors_finalise

subroutine f90wrap_magnitude(ret_magnitude, vec_1)
    use vectors, only: vector, magnitude
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    real(8), intent(out) :: ret_magnitude
    type(vector_ptr_type) :: vec_1_ptr
    integer, intent(in), dimension(2) :: vec_1
    vec_1_ptr = transfer(vec_1, vec_1_ptr)
    ret_magnitude = magnitude(vec_1=vec_1_ptr%p)
end subroutine f90wrap_magnitude

subroutine f90wrap_array_to_vector(vec_result, array)
    use vectors, only: assignment(=), vector
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type(vector_ptr_type) :: vec_result_ptr
    integer, intent(in), dimension(2) :: vec_result
    real(8), dimension(3), intent(in) :: array
    vec_result_ptr = transfer(vec_result, vec_result_ptr)
    vec_result_ptr%p = array
end subroutine f90wrap_array_to_vector

subroutine f90wrap_vector_to_array(array_result, vec_1)
    use vectors, only: assignment(=), vector
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    real(8), dimension(3), intent(inout) :: array_result
    type(vector_ptr_type) :: vec_1_ptr
    integer, intent(in), dimension(2) :: vec_1
    vec_1_ptr = transfer(vec_1, vec_1_ptr)
    array_result = vec_1_ptr%p
end subroutine f90wrap_vector_to_array

! End of module vectors defined in file linalg_module.fpp

! Module rotations defined in file linalg_module.fpp

subroutine f90wrap_euler_matrix(r, dim, angle, ap, n0, n1)
    use rotations, only: euler_matrix
    implicit none
    
    real(8), intent(inout), dimension(n0,n1) :: r
    integer, intent(in) :: dim
    real(8), intent(in) :: angle
    character(1), intent(in) :: ap
    integer :: n0
    !f2py intent(hide), depend(r) :: n0 = shape(r,0)
    integer :: n1
    !f2py intent(hide), depend(r) :: n1 = shape(r,1)
    call euler_matrix(R=r, dim=dim, angle=angle, ap=ap)
end subroutine f90wrap_euler_matrix

subroutine f90wrap_rotate_vector_single(vec, rotation_matrix, n0, n1)
    use rotations, only: rotate_vector
    use vectors, only: vector
    implicit none
    
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    type(vector_ptr_type) :: vec_ptr
    integer, intent(in), dimension(2) :: vec
    real(8), intent(in), dimension(n0,n1) :: rotation_matrix
    integer :: n0
    !f2py intent(hide), depend(rotation_matrix) :: n0 = shape(rotation_matrix,0)
    integer :: n1
    !f2py intent(hide), depend(rotation_matrix) :: n1 = shape(rotation_matrix,1)
    vec_ptr = transfer(vec, vec_ptr)
    call rotate_vector(vec=vec_ptr%p, rotation_matrix=rotation_matrix)
end subroutine f90wrap_rotate_vector_single

subroutine f90wrap_rotate_vector_array(vec, rotation_matrix, n0, n1)
    use rotations, only: rotate_vector
    use vectors, only: array_vectors
    implicit none
    
    type array_vectors_ptr_type
        type(array_vectors), pointer :: p => NULL()
    end type array_vectors_ptr_type
    type(array_vectors_ptr_type) :: vec_ptr
    integer, intent(in), dimension(2) :: vec
    real(8), intent(in), dimension(n0,n1) :: rotation_matrix
    integer :: n0
    !f2py intent(hide), depend(rotation_matrix) :: n0 = shape(rotation_matrix,0)
    integer :: n1
    !f2py intent(hide), depend(rotation_matrix) :: n1 = shape(rotation_matrix,1)
    vec_ptr = transfer(vec, vec_ptr)
    call rotate_vector(vec=vec_ptr%p, rotation_matrix=rotation_matrix)
end subroutine f90wrap_rotate_vector_array

! End of module rotations defined in file linalg_module.fpp

! Module basis_representations defined in file linalg_module.fpp

subroutine f90wrap_basis__array_getitem__u(f90wrap_this, f90wrap_i, uitem)
    
    use basis_representations, only: basis
    use vectors, only: vector
    implicit none
    
    type basis_ptr_type
        type(basis), pointer :: p => NULL()
    end type basis_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(basis_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: uitem(2)
    type(vector_ptr_type) :: u_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%u)) then
        call f90wrap_abort("array index out of range")
    else
        u_ptr%p => this_ptr%p%u(f90wrap_i)
        uitem = transfer(u_ptr,uitem)
    endif
end subroutine f90wrap_basis__array_getitem__u

subroutine f90wrap_basis__array_setitem__u(f90wrap_this, f90wrap_i, uitem)
    
    use basis_representations, only: basis
    use vectors, only: vector
    implicit none
    
    type basis_ptr_type
        type(basis), pointer :: p => NULL()
    end type basis_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(basis_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: uitem(2)
    type(vector_ptr_type) :: u_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%u)) then
        call f90wrap_abort("array index out of range")
    else
        u_ptr = transfer(uitem,u_ptr)
        this_ptr%p%u(f90wrap_i) = u_ptr%p
    endif
end subroutine f90wrap_basis__array_setitem__u

subroutine f90wrap_basis__array_len__u(f90wrap_this, f90wrap_n)
    
    use basis_representations, only: basis
    use vectors, only: vector
    implicit none
    
    type basis_ptr_type
        type(basis), pointer :: p => NULL()
    end type basis_ptr_type
    type vector_ptr_type
        type(vector), pointer :: p => NULL()
    end type vector_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(basis_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    f90wrap_n = size(this_ptr%p%u)
end subroutine f90wrap_basis__array_len__u

subroutine f90wrap_basis_initialise(this)
    use basis_representations, only: basis
    implicit none
    
    type basis_ptr_type
        type(basis), pointer :: p => NULL()
    end type basis_ptr_type
    type(basis_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_basis_initialise

subroutine f90wrap_basis_finalise(this)
    use basis_representations, only: basis
    implicit none
    
    type basis_ptr_type
        type(basis), pointer :: p => NULL()
    end type basis_ptr_type
    type(basis_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_basis_finalise

subroutine f90wrap_initialise_basis(this)
    use basis_representations, only: basis, initialise_basis
    implicit none
    
    type basis_ptr_type
        type(basis), pointer :: p => NULL()
    end type basis_ptr_type
    type(basis_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call initialise_basis(this=this_ptr%p)
end subroutine f90wrap_initialise_basis

subroutine f90wrap_mgramschmidt(vecs, e)
    use basis_representations, only: mgramschmidt, basis
    implicit none
    
    type basis_ptr_type
        type(basis), pointer :: p => NULL()
    end type basis_ptr_type
    type(basis_ptr_type) :: vecs_ptr
    integer, intent(in), dimension(2) :: vecs
    type(basis_ptr_type) :: e_ptr
    integer, intent(in), dimension(2) :: e
    vecs_ptr = transfer(vecs, vecs_ptr)
    e_ptr = transfer(e, e_ptr)
    call mgramschmidt(vecs=vecs_ptr%p, e=e_ptr%p)
end subroutine f90wrap_mgramschmidt

! End of module basis_representations defined in file linalg_module.fpp

