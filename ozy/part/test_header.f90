program read_header
    implicit none

    character(len=128) :: repository,line, final_line
    character(len=10) :: a,b,c,d,e,f,g,h,i
    integer :: status

    call read_params

    open(unit=10,file=repository,form='formatted',status='old')

    do
        read(10,'(A)',iostat=status)line
        if (status /= 0) exit
        final_line = line
    end do
    read(final_line,*)a,b,c,d,e,f
    write(*,*)f
    if (trim(f) .eq. 'family') then
        write(*,*)'This simulation uses particle families'
    else if (trim(f) .eq. 'tform') then
        write(*,*)'This simulation uses the old particle format'
    endif


    contains
    subroutine read_params
        implicit none

        integer :: i,n
        integer :: iargc
        character(len=10) :: opt
        character(len=128) :: arg

        n = iargc()

        do i = 1,n,2
            call getarg(i,opt)
            if (i == n) then
                print '("option ",a4," has no argument")', opt
                stop 2
            end if
            call getarg(i+1,arg)
            select case (trim(opt))
            case ('-inp')
                repository = trim(arg)
            case default
                print '("unknown option ",a2," ignored")', opt
            end select
        end do

        return
    end subroutine read_params
end program read_header