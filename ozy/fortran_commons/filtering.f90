module filtering
    use local
    use io_ramses
    use hydro_commons

    type filter
        character(128) :: name
        integer :: ncond
        type(hydro_var), dimension(:), allocatable :: cond_vars
        type(hydro_var), dimension(:), allocatable :: cond_vars_comp
        character(2), dimension(:), allocatable :: cond_ops
        real(dbl), dimension(:), allocatable :: cond_vals
        logical, dimension(:), allocatable :: use_var
    end type filter

    contains

    subroutine allocate_filter(filt)
        implicit none
        type(filter), intent(inout) :: filt

        if (.not.allocated(filt%cond_vars)) allocate(filt%cond_vars(filt%ncond))
        if (.not.allocated(filt%cond_vars_comp)) allocate(filt%cond_vars_comp(filt%ncond))
        if (.not.allocated(filt%cond_ops)) allocate(filt%cond_ops(filt%ncond))
        if (.not.allocated(filt%cond_vals)) allocate(filt%cond_vals(filt%ncond))
        if (.not.allocated(filt%use_var)) allocate(filt%use_var(filt%ncond))
        filt%use_var = .false.
    end subroutine allocate_filter

    subroutine get_filter_var_tools(vardict,filt)
        implicit none

        type(dictf90),intent(in) :: vardict
        type(filter),intent(inout):: filt

        integer :: i

        if (filt%ncond == 0) return

        do i=1,filt%ncond
            call set_hydro_var(vardict,filt%cond_vars(i))
            if (filt%use_var(i)) call set_hydro_var(vardict,filt%cond_vars_comp(i))
        end do
    end subroutine get_filter_var_tools

    logical function filter_cell(reg,filt,cell_x,cell_dx,cell_var,cell_son,&
                                &trans_matrix,grav_var)
        use vectors
        use geometrical_regions
        type(region), intent(in) :: reg
        type(filter), intent(in) :: filt
        real(dbl), intent(in) :: cell_dx
        type(vector), intent(in) :: cell_x
        real(dbl), dimension(0:amr%twondim,1:sim%nvar), intent(in) :: cell_var
        integer,dimension(0:amr%twondim),intent(in) :: cell_son
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix
        real(dbl),dimension(0:amr%twondim,1:4),intent(in),optional :: grav_var
        integer :: i
        real(dbl) :: value,filt_value

        filter_cell = .true.

        if (filt%ncond == 0) return

        do i=1,filt%ncond
            if (present(grav_var)) then
                value = filt%cond_vars(i)%myfunction(amr,sim,filt%cond_vars(i),reg,cell_dx,&
                                                    cell_x,cell_var,cell_son,trans_matrix,&
                                                    grav_var)
            else
                value = filt%cond_vars(i)%myfunction(amr,sim,filt%cond_vars(i),reg,cell_dx,&
                                                    cell_x,cell_var,cell_son,trans_matrix)
            end if
            if (filt%use_var(i)) then
                if (present(grav_var)) then
                    value = filt%cond_vars_comp(i)%myfunction(amr,sim,filt%cond_vars_comp(i),reg,cell_dx,&
                                                            cell_x,cell_var,cell_son,trans_matrix,&
                                                            grav_var)
                else
                    value = filt%cond_vars_comp(i)%myfunction(amr,sim,filt%cond_vars_comp(i),reg,cell_dx,&
                                                            cell_x,cell_var,cell_son,trans_matrix)
                end if
                filt_value = filt%cond_vals(i) * filt_value
            else
                filt_value = filt%cond_vals(i)
            end if
            select case (TRIM(filt%cond_ops(i)))
            case('/=')
                filter_cell = filter_cell .and. (value /= filt_value)
            case('==')
                filter_cell = filter_cell .and. (value == filt_value)
            case('<')
                filter_cell = filter_cell .and. (value < filt_value)
            case('<=')
                filter_cell = filter_cell .and. (value <= filt_value)
            case('>')
                filter_cell = filter_cell .and. (value > filt_value)
            case('>=')
                filter_cell = filter_cell .and. (value >= filt_value)
            case default
                write(*,*)'Relation operator not supported: ',TRIM(filt%cond_ops(i))
                write(*,*)'Aborting!'
                stop
            end select
        end do
    end function filter_cell

    logical function filter_particle(reg,filt,part,dx)
        use geometrical_regions
        type(region),intent(in) :: reg
        type(filter),intent(in) :: filt
        type(particle),intent(in) :: part
        type(vector),optional,intent(in) :: dx
        integer :: i
        real(dbl) :: value

        filter_particle = .true.
        if (filt%ncond == 0) return

        do i=1,filt%ncond
            if (present(dx)) then
                call getpartvalue(reg,part,filt%cond_vars(i)%name,value,dx)
            else
                call getpartvalue(reg,part,filt%cond_vars(i)%name,value)
            endif
            select case (TRIM(filt%cond_ops(i)))
            case('/=')
                filter_particle = filter_particle .and. (value /= filt%cond_vals(i))
            case('==')
                filter_particle = filter_particle .and. (value == filt%cond_vals(i))
            case('<')
                filter_particle = filter_particle .and. (value < filt%cond_vals(i))
            case('<=')
                filter_particle = filter_particle .and. (value <= filt%cond_vals(i))
            case('>')
                filter_particle = filter_particle .and. (value > filt%cond_vals(i))
            case('>=')
                filter_particle = filter_particle .and. (value >= filt%cond_vals(i))
            case default
                write(*,*)'Relation operator not supported: ',TRIM(filt%cond_ops(i))
                write(*,*)'Aborting!'
                stop
            end select
        end do
    end function filter_particle

    logical function filter_sub(sub,cell_x)
        use vectors
        use geometrical_regions
        type(region),intent(in) :: sub
        real(dbl),dimension(1:3), intent(in) :: cell_x
        integer :: i
        real(dbl) :: distance
        real(dbl),dimension(1:3) :: pos

        filter_sub = .True.

        pos = cell_x - (/sub%centre%x,sub%centre%y,sub%centre%z/)

        call checkifinside(pos,sub,filter_sub,distance)
        filter_sub = .not.filter_sub
    end function filter_sub

end module filtering