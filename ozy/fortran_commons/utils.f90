module utils
    use local

    contains

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