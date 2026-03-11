module filtering_hydro
    use local
    use io_ramses
    use hydro_commons

    type filter_hydro
        character(128) :: name
        integer :: ncond
        type(hydro_var), dimension(:), allocatable :: cond_vars
        type(hydro_var), dimension(:), allocatable :: cond_vars_comp
        character(128), dimension(:), allocatable :: cond_vars_name
        character(128), dimension(:), allocatable :: cond_vars_comp_name
        character(2), dimension(:), allocatable :: cond_ops
        real(dbl), dimension(:), allocatable :: cond_vals
        logical, dimension(:), allocatable :: use_var
    end type filter_hydro

    contains

    subroutine allocate_filter_hydro(filt)
        implicit none
        type(filter_hydro), intent(inout) :: filt

        if (.not.allocated(filt%cond_vars)) allocate(filt%cond_vars(filt%ncond))
        if (.not.allocated(filt%cond_vars_name)) allocate(filt%cond_vars_name(filt%ncond))
        if (.not.allocated(filt%cond_vars_comp)) allocate(filt%cond_vars_comp(filt%ncond))
        if (.not.allocated(filt%cond_vars_comp_name)) allocate(filt%cond_vars_comp_name(filt%ncond))
        if (.not.allocated(filt%cond_ops)) allocate(filt%cond_ops(filt%ncond))
        if (.not.allocated(filt%cond_vals)) allocate(filt%cond_vals(filt%ncond))
        if (.not.allocated(filt%use_var)) allocate(filt%use_var(filt%ncond))
        filt%use_var = .false.
    end subroutine allocate_filter_hydro

    subroutine get_filter_var_tools(my_sim,vardict,filt)
        implicit none

        type(sim_info), intent(in) :: my_sim
        type(dictf90),intent(in) :: vardict
        type(filter_hydro),intent(inout) :: filt

        logical :: ok_check
        integer :: i, ivar

        ! Loop over the conditions
        if (filt%ncond == 0) return
        do i=1,filt%ncond
            ! 1. If the condition is 'none' just ignore this filter and set to 0 conds
            if (TRIM(filt%cond_vars_name(i)) == 'none') then
                filt%ncond = 0
                cycle
            end if

            ! 2. Set the variable
            filt%cond_vars(i)%name = filt%cond_vars_name(i)
            call set_hydro_var(my_sim,vardict,filt%cond_vars(i))
            if (filt%use_var(i)) then
                filt%cond_vars_comp(i)%name = filt%cond_vars_comp_name(i)
                call set_hydro_var(my_sim,vardict,filt%cond_vars_comp(i))
            end if            
        
        end do ! i

    end subroutine get_filter_var_tools

    logical function filter_cell(reg,filt,cell_x,cell_dx,cell_var,cell_son,&
                                &trans_matrix,grav_var,rt_var)
        use vectors
        use geometrical_regions
        type(region), intent(in) :: reg
        type(filter_hydro), intent(in) :: filt
        real(dbl), intent(in) :: cell_dx
        type(vector), intent(in) :: cell_x
        real(dbl), dimension(0:amr%twondim,1:sim%nvar), intent(in) :: cell_var
        integer,dimension(0:amr%twondim),intent(in) :: cell_son
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix
        real(dbl),dimension(0:amr%twondim,1:4),intent(in),optional :: grav_var
        real(dbl),dimension(0:amr%twondim,1:rtinfo%nRTvar),intent(in),optional :: rt_var

        integer :: i
        real(dbl) :: value,filt_value

        filter_cell = .true.

        if (filt%ncond == 0) return

        do i=1,filt%ncond
            if (present(grav_var).and.present(rt_var)) then
                value = filt%cond_vars(i)%myfunction(amr,sim,rtinfo,filt%cond_vars(i),reg,cell_dx,&
                                                    cell_x,cell_var,cell_son,trans_matrix,&
                                                    grav_var,rt_var)
            elseif (present(grav_var)) then
                value = filt%cond_vars(i)%myfunction(amr,sim,rtinfo,filt%cond_vars(i),reg,cell_dx,&
                                                    cell_x,cell_var,cell_son,trans_matrix,&
                                                    grav_var)
            elseif (present(rt_var)) then
                value = filt%cond_vars(i)%myfunction(amr,sim,rtinfo,filt%cond_vars(i),reg,cell_dx,&
                                                    cell_x,cell_var,cell_son,trans_matrix, rt_var=rt_var)
            else
                value = filt%cond_vars(i)%myfunction(amr,sim,rtinfo,filt%cond_vars(i),reg,cell_dx,&
                                                    cell_x,cell_var,cell_son,trans_matrix)
            end if
            if (filt%use_var(i)) then
                if (present(grav_var).and.present(rt_var)) then
                    filt_value = filt%cond_vars_comp(i)%myfunction(amr,sim,rtinfo,filt%cond_vars_comp(i),reg,cell_dx,&
                                                            cell_x,cell_var,cell_son,trans_matrix,&
                                                            grav_var,rt_var)
                elseif (present(grav_var)) then
                    filt_value = filt%cond_vars_comp(i)%myfunction(amr,sim,rtinfo,filt%cond_vars_comp(i),reg,cell_dx,&
                                                            cell_x,cell_var,cell_son,trans_matrix,&
                                                            grav_var)
                elseif (present(rt_var)) then
                    filt_value = filt%cond_vars_comp(i)%myfunction(amr,sim,rtinfo,filt%cond_vars_comp(i),reg,cell_dx,&
                                                            cell_x,cell_var,cell_son,trans_matrix, rt_var=rt_var)
                else
                    filt_value = filt%cond_vars_comp(i)%myfunction(amr,sim,rtinfo,filt%cond_vars_comp(i),reg,cell_dx,&
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

end module filtering_hydro