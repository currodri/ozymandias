module utils
    use local

    contains

    ! Subroutine to extract words from a line
    subroutine get_word(line, start_pos, word)
        implicit none
        character(len=*), intent(in) :: line
        integer, intent(inout) :: start_pos
        character(len=50), intent(out) :: word
        integer :: i, len_line

        len_line = len_trim(line)
        word = ""

        ! Skip leading spaces
        do while (start_pos <= len_line .and. line(start_pos:start_pos) == ' ')
            start_pos = start_pos + 1
        end do

        ! Extract the word
        i = 1
        do while (start_pos <= len_line .and. line(start_pos:start_pos) /= ' ')
            word(i:i) = line(start_pos:start_pos)
            start_pos = start_pos + 1
            i = i + 1
        end do

        ! Skip trailing spaces
        do while (start_pos <= len_line .and. line(start_pos:start_pos) == ' ')
            start_pos = start_pos + 1
        end do

    end subroutine get_word

    function check_number_suffix(str) result(is_number)
        implicit none
        character(len=*), intent(in) :: str
        logical :: is_number
        integer :: i, len_str

        is_number = .false.
        len_str = len_trim(str)

        ! Find last underscore ('_')
        do i = len_str, 1, -1
            if (str(i:i) == '_') then
                ! Check if remaining substring contains only digits
                if (i < len_str) then
                    is_number = all_digits(str(i+1:len_str))
                end if
                exit
            end if
        end do
    end function check_number_suffix

    function all_digits(substr) result(is_digit)
        implicit none
        character(len=*), intent(in) :: substr
        logical :: is_digit
        integer :: j, len_sub

        is_digit = .true.
        len_sub = len_trim(substr)

        ! Check each character is a digit
        do j = 1, len_sub
            if (index('0123456789', substr(j:j)) == 0) then
                is_digit = .false.
                exit
            end if
        end do
    end function all_digits

    function get_cleaned_string(str) result(cleaned_str)
        implicit none
        character(len=*), intent(in) :: str
        character(len=len(str)) :: cleaned_str
        integer :: i, len_str

        cleaned_str = str  ! Default: same as input
        len_str = len_trim(str)

        ! Find last underscore ('_')
        do i = len_str, 1, -1
            if (str(i:i) == '_') then
                ! Check if remaining substring contains only digits
                if (i < len_str) then
                    if (all_digits(str(i+1:len_str))) then
                        cleaned_str = str(1:i-1)  ! Remove '_number' part
                    end if
                end if
                exit
            end if
        end do
    end function get_cleaned_string

    SUBROUTINE quick_sort_irg(list, order, n)
        IMPLICIT NONE
        ! Quick sort routine from:
        ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
        ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
        ! Modified by Alan Miller to include an associated integer array which gives
        ! the positions of the elements in the original order.
        
        
        INTEGER :: n
        INTEGER(irg), DIMENSION (1:n), INTENT(INOUT)  :: list
        INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order
        
        ! Local variable
        INTEGER :: i
        
        DO i = 1, n
            order(i) = i
        END DO
        
        CALL quick_sort_1_int(1, n)
        
        CONTAINS
        
        RECURSIVE SUBROUTINE quick_sort_1_int(left_end, right_end)
        
            INTEGER, INTENT(IN) :: left_end, right_end
        
            !     Local variables
            INTEGER             :: i, j, itemp
            INTEGER(irg)        :: reference, temp
            INTEGER, PARAMETER  :: max_simple_sort_size = 6
        
            IF (right_end < left_end + max_simple_sort_size) THEN
                ! Use interchange sort for small lists
                CALL interchange_sort_int(left_end, right_end)
        
            ELSE
                ! Use partition ("quick") sort
                reference = list((left_end + right_end)/2)
                i = left_end - 1; j = right_end + 1
        
            DO
                ! Scan list from left end until element >= reference is found
                DO
                    i = i + 1
                    IF (list(i) >= reference) EXIT
                END DO
                ! Scan list from right end until element <= reference is found
                DO
                    j = j - 1
                    IF (list(j) <= reference) EXIT
                END DO
        
        
                IF (i < j) THEN
                    ! Swap two out-of-order elements
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                ELSE IF (i == j) THEN
                    i = i + 1
                    EXIT
                ELSE
                    EXIT
                END IF
                END DO
                IF (left_end < j) CALL quick_sort_1_int(left_end, j)
                IF (i < right_end) CALL quick_sort_1_int(i, right_end)
            END IF
        
        END SUBROUTINE quick_sort_1_int
        
        
        SUBROUTINE interchange_sort_int(left_end, right_end)
        
            INTEGER, INTENT(IN) :: left_end, right_end
        
            !     Local variables                                                                                                                                                                           
            INTEGER             :: i, j, itemp
            INTEGER(irg)        :: temp
        
            DO i = left_end, right_end - 1
                DO j = i+1, right_end
                IF (list(i) > list(j)) THEN
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                END IF
                END DO
            END DO
        
        END SUBROUTINE interchange_sort_int
        
    END SUBROUTINE quick_sort_irg

    SUBROUTINE quick_sort_ilg(list, order, n)
        IMPLICIT NONE
        ! Quick sort routine from:
        ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
        ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
        ! Modified by Alan Miller to include an associated integer array which gives
        ! the positions of the elements in the original order.
        
        
        INTEGER :: n
        INTEGER(ilg), DIMENSION (1:n), INTENT(INOUT)  :: list
        INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order
        
        ! Local variable
        INTEGER :: i
        
        DO i = 1, n
            order(i) = i
        END DO
        
        CALL quick_sort_1_int(1, n)
        
        CONTAINS
        
        RECURSIVE SUBROUTINE quick_sort_1_int(left_end, right_end)
        
            INTEGER, INTENT(IN) :: left_end, right_end
        
            !     Local variables
            INTEGER             :: i, j, itemp
            INTEGER(ilg)        :: reference, temp
            INTEGER, PARAMETER  :: max_simple_sort_size = 6
        
            IF (right_end < left_end + max_simple_sort_size) THEN
                ! Use interchange sort for small lists
                CALL interchange_sort_int(left_end, right_end)
        
            ELSE
                ! Use partition ("quick") sort
                reference = list((left_end + right_end)/2)
                i = left_end - 1; j = right_end + 1
        
            DO
                ! Scan list from left end until element >= reference is found
                DO
                    i = i + 1
                    IF (list(i) >= reference) EXIT
                END DO
                ! Scan list from right end until element <= reference is found
                DO
                    j = j - 1
                    IF (list(j) <= reference) EXIT
                END DO
        
        
                IF (i < j) THEN
                    ! Swap two out-of-order elements
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                ELSE IF (i == j) THEN
                    i = i + 1
                    EXIT
                ELSE
                    EXIT
                END IF
                END DO
                IF (left_end < j) CALL quick_sort_1_int(left_end, j)
                IF (i < right_end) CALL quick_sort_1_int(i, right_end)
            END IF
        
        END SUBROUTINE quick_sort_1_int
        
        
        SUBROUTINE interchange_sort_int(left_end, right_end)
        
            INTEGER, INTENT(IN) :: left_end, right_end
        
            !     Local variables                                                                                                                                                                           
            INTEGER             :: i, j, itemp
            INTEGER(ilg)        :: temp
        
            DO i = left_end, right_end - 1
                DO j = i+1, right_end
                IF (list(i) > list(j)) THEN
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                END IF
                END DO
            END DO
        
        END SUBROUTINE interchange_sort_int
        
    END SUBROUTINE quick_sort_ilg

    SUBROUTINE quick_sort_dp(list, order, n)
        IMPLICIT NONE
        ! Quick sort routine from:
        ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
        ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
        ! Modified by Alan Miller to include an associated integer array which gives
        ! the positions of the elements in the original order.
        
        
        INTEGER :: n
        REAL(dbl), DIMENSION (1:n), INTENT(INOUT)  :: list
        INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order
        
        ! Local variable
        INTEGER :: i
        
        DO i = 1, n
            order(i) = i
        END DO
        
        CALL quick_sort_1_dp(1, n)
        
        CONTAINS
        
        RECURSIVE SUBROUTINE quick_sort_1_dp(left_end, right_end)
        
            INTEGER, INTENT(IN) :: left_end, right_end
        
            !     Local variables
            INTEGER             :: i, j, itemp
            REAL(kind=8)              :: reference, temp
            INTEGER, PARAMETER  :: max_simple_sort_size = 6
        
            IF (right_end < left_end + max_simple_sort_size) THEN
                ! Use interchange sort for small lists
                CALL interchange_sort_dp(left_end, right_end)
        
            ELSE
                ! Use partition ("quick") sort
                reference = list((left_end + right_end)/2)
                i = left_end - 1; j = right_end + 1
        
            DO
                ! Scan list from left end until element >= reference is found
                DO
                    i = i + 1
                    IF (list(i) >= reference) EXIT
                END DO
                ! Scan list from right end until element <= reference is found
                DO
                    j = j - 1
                    IF (list(j) <= reference) EXIT
                END DO
        
        
                IF (i < j) THEN
                    ! Swap two out-of-order elements
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                ELSE IF (i == j) THEN
                    i = i + 1
                    EXIT
                ELSE
                    EXIT
                END IF
                END DO
                IF (left_end < j) CALL quick_sort_1_dp(left_end, j)
                IF (i < right_end) CALL quick_sort_1_dp(i, right_end)
            END IF
        
        END SUBROUTINE quick_sort_1_dp
        
        
        SUBROUTINE interchange_sort_dp(left_end, right_end)
        
            INTEGER, INTENT(IN) :: left_end, right_end
        
            !     Local variables                                                                                                                                                                           
            INTEGER             :: i, j, itemp
            REAL(kind=8)           :: temp
        
            DO i = left_end, right_end - 1
                DO j = i+1, right_end
                IF (list(i) > list(j)) THEN
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                END IF
                END DO
            END DO
        
        END SUBROUTINE interchange_sort_dp
        
    END SUBROUTINE quick_sort_dp

    subroutine binarysearch_irg(n,list,goal,found)
        implicit none
        
        integer,intent(in) :: n
        integer(irg),dimension(1:n),intent(in) :: list
        integer(irg),intent(in) :: goal
        logical,intent(inout) :: found

        integer :: low
        integer :: high
        integer :: middle

        low = 1
        high = n
        found = .false.
        
        do while(low <= high .and. .not. found)
            middle = (low + high)/2
            if (goal == list(middle)) then
                found = .true.
            elseif (goal < list(middle)) then
                high = middle - 1
            else
                low = middle +1
            endif
        end do

    end subroutine binarysearch_irg

    subroutine binarysearch_ilg(n,list,goal,found)
        implicit none
        
        integer,intent(in) :: n
        integer(ilg),dimension(1:n),intent(in) :: list
        integer(ilg),intent(in) :: goal
        logical,intent(inout) :: found

        integer :: low
        integer :: high
        integer :: middle

        low = 1
        high = n
        found = .false.
        
        do while(low <= high .and. .not. found)
            middle = (low + high)/2
            if (goal == list(middle)) then
                found = .true.
            elseif (goal < list(middle)) then
                high = middle - 1
            else
                low = middle +1
            endif
        end do

    end subroutine binarysearch_ilg

    subroutine binarysearch_dp(n,list,goal,ifound)
        implicit none
        
        integer,intent(in) :: n
        real(dbl),dimension(1:n),intent(in) :: list
        real(dbl),intent(in) :: goal
        logical :: found
        integer,intent(inout) :: ifound

        integer :: low
        integer :: high
        integer :: middle

        low = 1
        high = n
        found = .false.
        ifound = 0
        
        do while(low <= high .and. .not. found)
            middle = (low + high)/2
            if (goal == list(middle)) then
                found = .true.
                ifound = middle
            elseif (goal < list(middle)) then
                high = middle - 1
            else
                low = middle +1
            endif
        end do

    end subroutine binarysearch_dp

    !***********************************************************************
    subroutine indexx(n,arr,indx)
    
    implicit none
    
    integer(kind=4)           :: n,indx(n)
    real(dbl)              :: arr(n)
    integer(kind=4),parameter :: m=7,nstack=50
    integer(kind=4)           :: i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
    real(dbl)              :: a
    
    do j = 1,n
        indx(j) = j
    enddo
    
    jstack = 0
    l      = 1
    ir     = n
    1 if (ir-l .lt. m) then
        do j = l+1,ir
            indxt = indx(j)
            a     = arr(indxt)
            do i = j-1,1,-1
            if (arr(indx(i)) .le. a) goto 2
            indx(i+1) = indx(i)
            enddo
            i         = 0
    2       indx(i+1) = indxt
        enddo
        if (jstack .eq. 0) return
        ir     = istack(jstack)
        l      = istack(jstack-1)
        jstack = jstack-2
    else
        k         = (l+ir)/2
        itemp     = indx(k)
        indx(k)   = indx(l+1)
        indx(l+1) = itemp
        if (arr(indx(l+1)) .gt. arr(indx(ir))) then
            itemp     = indx(l+1)
            indx(l+1) = indx(ir)
            indx(ir)  = itemp
        endif
        if (arr(indx(l)) .gt. arr(indx(ir))) then
            itemp    = indx(l)
            indx(l)  = indx(ir)
            indx(ir) = itemp
        endif
        if (arr(indx(l+1)) .gt. arr(indx(l))) then
            itemp     = indx(l+1)
            indx(l+1) = indx(l)
            indx(l)   = itemp
        endif
        i     = l+1
        j     = ir
        indxt = indx(l)
        a     = arr(indxt)
    3    continue
        i     = i+1
        if (arr(indx(i)) .lt. a) goto 3
    4    continue
        j     = j-1
        if (arr(indx(j)) .gt. a) goto 4
        if (j .lt. i) goto 5
        itemp   = indx(i)
        indx(i) = indx(j)
        indx(j) = itemp
        goto 3
    5    continue
        indx(l) = indx(j)
        indx(j) = indxt
        jstack  = jstack+2
        if (jstack .gt. nstack) stop 'nstack too small in indexx'
        if (ir-i+1 .ge. j-l) then
            istack(jstack)   = ir
            istack(jstack-1) = i
            ir               = j-1
        else
            istack(jstack)   = j-1
            istack(jstack-1) = l
            l                = i
        endif
    endif
    goto 1
    
    end subroutine indexx   

    subroutine locate(xx,n,x,j)
 
    implicit none
    
    integer(kind=4) ::  n
    real(dbl)    ::  xx(n),x
    integer(kind=4) ::  j,jl,ju,jm
    
    jl = 0
    ju = n+1
    do while (ju-jl .gt. 1) 
        jm = (ju+jl)/2
        if ((xx(n) .gt. xx(1)) .eqv. (x .gt. xx(jm))) then
            jl = jm
        else
            ju = jm
        endif
    enddo
    j = jl
    
    return
    
    end subroutine locate
    ! subroutine bspline_regridding(n,x,y,z,fcn,o,new_n,new_x,new_y,new_z,new_fcn,extrap)
    !     use bspline_module

    !     implicit none
    !     integer,dimension(1:3),intent(in) :: n
    !     real(dbl),dimension(1:n(1)),intent(in) :: x
    !     real(dbl),dimension(1:n(2)),intent(in) :: y
    !     real(dbl),dimension(1:n(3)),intent(in) :: z
    !     real(dbl),dimension(1:n(1),1:n(2),1:n(3)), intent(in) :: fcn
    !     integer,dimension(1:3),intent(in) :: o
    !     integer,dimension(1:3),intent(in) :: new_n
    !     real(dbl),dimension(1:new_n(1)),intent(in) :: new_x
    !     real(dbl),dimension(1:new_n(2)),intent(in) :: new_y
    !     real(dbl),dimension(1:new_n(3)),intent(in) :: new_z
    !     real(dbl),dimension(1:new_n(1),1:new_n(2),1:new_n(3)), intent(inout) :: new_fcn
    !     logical,optional :: extrap

    !     integer :: i,j,k
    !     integer :: iflag  !! status flag
    !     integer :: inbvx,inbvy,inbvz,iloy,iloz
    !     integer :: idx=0,idy=0,idz=0
    !     real(dbl) :: val
    !     real(dbl),dimension(1:n(1),1:n(2),1:n(3)) :: bcoeff
    !     integer,parameter :: iknot  = 0    !! automatically select the knots
    !     real(dbl),dimension(n(1)+o(1))    :: tx       !! x knots
    !     real(dbl),dimension(n(2)+o(2))    :: ty       !! y knots
    !     real(dbl),dimension(n(3)+o(3))    :: tz       !! z knots
    !     real(dbl),dimension(o(2),o(3))    :: w2
    !     real(dbl),dimension(o(3))         :: w1
    !     real(dbl),dimension(3*maxval(o))     :: w0

    !     inbvx = 1
    !     inbvy = 1
    !     inbvz = 1
    !     iloy = 1
    !     iloz = 1

    !     call db3ink(x,n,y,n,z,n,fcn,o(1),o(2),o(3),iknot,tx,ty,tz,bcoeff,iflag)

    !     if (iflag/=0) error stop 'error calling db2ink'
    !     do i=1,new_n(1)
    !         do j=1,new_n(2)
    !             do k=1,new_n(3)
    !                 call db3val(new_x(i),new_y(j),new_z(k),idx,idy,idz,tx,ty,tz,n(1),n(2),n(3),o(1),o(2),o(3),bcoeff,val,iflag,&
    !                             inbvx,inbvy,inbvz,iloy,iloz,w1,w2,extrap)
    !                 if (iflag/=0) error stop 'error calling db2val'
    !                 new_fcn(i,j,k) = val
    !             end do
    !         end do
    !     end do

    ! end subroutine bspline_regridding
end module utils