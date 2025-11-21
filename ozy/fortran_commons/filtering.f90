module filtering
    use local
    use io_ramses
    use hydro_commons
    use part_commons

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

    type filter_part
        character(128) :: name
        integer :: ncond_d,ncond_i,ncond_b,ncond
        integer, dimension(:), allocatable :: cond_vtype
        type(part_var), dimension(:), allocatable :: cond_vars
        type(part_var), dimension(:), allocatable :: cond_vars_comp
        character(128), dimension(:), allocatable :: cond_vars_name
        character(128), dimension(:), allocatable :: cond_vars_comp_name
        character(2), dimension(:), allocatable :: cond_ops
        real(dbl), dimension(:), allocatable :: cond_vals_d
#ifdef LONGINT
        integer(ilg), dimension(:), allocatable :: cond_vals_i
#else
        integer(irg), dimension(:), allocatable :: cond_vals_i
#endif
        integer(1), dimension(:), allocatable :: cond_vals_b
        logical, dimension(:), allocatable :: use_var
    end type filter_part

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

    subroutine allocate_filter_part(filt)
        implicit none
        type(filter_part), intent(inout) :: filt

        if (.not.allocated(filt%cond_vars)) allocate(filt%cond_vars(filt%ncond))
        if (.not.allocated(filt%cond_vars_name)) allocate(filt%cond_vars_name(filt%ncond))
        if (.not.allocated(filt%cond_vars_comp)) allocate(filt%cond_vars_comp(filt%ncond))
        if (.not.allocated(filt%cond_vars_comp_name)) allocate(filt%cond_vars_comp_name(filt%ncond))
        if (.not.allocated(filt%cond_ops)) allocate(filt%cond_ops(filt%ncond))
        if (.not.allocated(filt%cond_vtype)) allocate(filt%cond_vtype(filt%ncond))
        if (.not.allocated(filt%cond_vals_d)) allocate(filt%cond_vals_d(filt%ncond))
        if (.not.allocated(filt%cond_vals_i)) allocate(filt%cond_vals_i(filt%ncond))
        if (.not.allocated(filt%cond_vals_b)) allocate(filt%cond_vals_b(filt%ncond))
        if (.not.allocated(filt%use_var)) allocate(filt%use_var(filt%ncond))
        filt%use_var = .false.
    end subroutine allocate_filter_part

    subroutine get_filter_var_tools(vardict,filt)
        implicit none

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
            call set_hydro_var(vardict,filt%cond_vars(i))
            if (filt%use_var(i)) then
                filt%cond_vars_comp(i)%name = filt%cond_vars_comp_name(i)
                call set_hydro_var(vardict,filt%cond_vars_comp(i))
            end if            
        
        end do ! i

    end subroutine get_filter_var_tools

    subroutine get_filter_part_tools(vardict,vtypedict,filt)
        implicit none

        type(dictf90),intent(in) :: vardict, vtypedict
        type(filter_part),intent(inout) :: filt

        logical :: ok_check
        integer :: i, ivar


        ! Loop over the conditions
        if (filt%ncond == 0) return
        do i = 1, filt%ncond
            ! 1. If the condition is 'none' just ignore this filter and set to 0 conds
            if (TRIM(filt%cond_vars_name(i)) == 'none') then
                filt%ncond = 0
                cycle
            end if

            ! 2. Set the variable
            filt%cond_vars(i)%name = filt%cond_vars_name(i)
            call set_part_var(vardict,vtypedict,filt%cond_vars(i))
            if (filt%cond_vars(i)%vartype==1) then
                filt%ncond_d = filt%ncond_d + 1
                filt%cond_vtype(i) = 1
            elseif (filt%cond_vars(i)%vartype==2) then
                filt%ncond_i = filt%ncond_i + 1
                filt%cond_vtype(i) = 2
            elseif (filt%cond_vars(i)%vartype==3) then
                filt%ncond_b = filt%ncond_b + 1
                filt%cond_vtype(i) = 3
            end if
            if (filt%use_var(i)) then
                filt%cond_vars_comp(i)%name = filt%cond_vars_comp_name(i)
                call set_part_var(vardict,vtypedict,filt%cond_vars_comp(i))
            end if
        end do ! i
    end subroutine get_filter_part_tools

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
#if RTPRE==4
        real(sgl),dimension(0:amr%twondim,1:rtinfo%nRTvar),intent(in),optional :: rt_var
#elif RTPRE==8
        real(dbl),dimension(0:amr%twondim,1:rtinfo%nRTvar),intent(in),optional :: rt_var
#endif
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

    logical function filter_particle(reg,filt,dx,part_var_d,part_var_i,part_var_b)
        
        use vectors
        use geometrical_regions
        type(region), intent(in) :: reg
        type(filter_part), intent(in) :: filt
        type(vector), intent(in) :: dx
        real(dbl), dimension(1:sim%nvar_part_d), intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg), dimension(1:sim%nvar_part_i), intent(in) :: part_var_i
#else
        integer(irg), dimension(1:sim%nvar_part_i), intent(in) :: part_var_i
#endif
        integer(1), dimension(1:sim%nvar_part_b), intent(in) :: part_var_b

        integer :: i, counter_d, counter_i, counter_b
        real(dbl) :: value_d,filt_value_d
#ifdef LONGINT
        integer(ilg) :: value_i,filt_value_i
#else
        integer(irg) :: value_i,filt_value_i
#endif
        integer(1) :: value_b,filt_value_b

        filter_particle = .true.

        if (filt%ncond == 0) return

        do i = 1, filt%ncond
            if (filt%cond_vtype(i) == 1) then
                value_d = filt%cond_vars(i)%myfunction_d(amr,sim,filt%cond_vars(i),reg,dx,&
                                                    part_var_d,part_var_i,part_var_b)
                if (filt%use_var(i)) then
                    filt_value_d = filt%cond_vars_comp(i)%myfunction_d(amr,sim,filt%cond_vars_comp(i),reg,dx,&
                                                            part_var_d,part_var_i,part_var_b)
                    filt_value_d = filt%cond_vals_d(i) * filt_value_d
                else
                    filt_value_d = filt%cond_vals_d(i)
                end if
                select case (TRIM(filt%cond_ops(i)))
                case('/=')
                    filter_particle = filter_particle .and. (value_d /= filt_value_d)
                case('==')
                    filter_particle = filter_particle .and. (value_d == filt_value_d)
                case('<')
                    filter_particle = filter_particle .and. (value_d < filt_value_d)
                case('<=')
                    filter_particle = filter_particle .and. (value_d <= filt_value_d)
                case('>')
                    filter_particle = filter_particle .and. (value_d > filt_value_d)
                case('>=')
                    filter_particle = filter_particle .and. (value_d >= filt_value_d)
                case default
                    write(*,*)'Relation operator not supported: ',TRIM(filt%cond_ops(i))
                    write(*,*)'Aborting!'
                    stop
                end select
            else if (filt%cond_vtype(i) == 2) then
                value_i = filt%cond_vars(i)%myfunction_i(amr,sim,filt%cond_vars(i),reg,dx,&
                                                    part_var_d,part_var_i,part_var_b)
                if (filt%use_var(i)) then
                    filt_value_i = filt%cond_vars_comp(i)%myfunction_i(amr,sim,filt%cond_vars_comp(i),reg,dx,&
                                                            part_var_d,part_var_i,part_var_b)
                    filt_value_i = filt%cond_vals_i(i) * filt_value_i
                else
                    filt_value_i = filt%cond_vals_i(i)
                end if
                select case (TRIM(filt%cond_ops(i)))
                case('/=')
                    filter_particle = filter_particle .and. (value_i /= filt_value_i)
                case('==')
                    filter_particle = filter_particle .and. (value_i == filt_value_i)
                case('<')
                    filter_particle = filter_particle .and. (value_i < filt_value_i)
                case('<=')
                    filter_particle = filter_particle .and. (value_i <= filt_value_i)
                case('>')
                    filter_particle = filter_particle .and. (value_i > filt_value_i)
                case('>=')
                    filter_particle = filter_particle .and. (value_i >= filt_value_i)
                case default
                    write(*,*)'Relation operator not supported: ',TRIM(filt%cond_ops(i))
                    write(*,*)'Aborting!'
                    stop
                end select
            else if (filt%cond_vtype(i) == 3) then
                value_b = filt%cond_vars(i)%myfunction_b(amr,sim,filt%cond_vars(i),reg,dx,&
                                                    part_var_d,part_var_i,part_var_b)
                if (filt%use_var(i)) then
                    filt_value_b = filt%cond_vars_comp(i)%myfunction_b(amr,sim,filt%cond_vars_comp(i),reg,dx,&
                                                            part_var_d,part_var_i,part_var_b)
                    filt_value_b = filt%cond_vals_b(i) * filt_value_b
                else
                    filt_value_b = filt%cond_vals_b(i)
                end if
                select case (TRIM(filt%cond_ops(i)))
                case('/=')
                    filter_particle = filter_particle .and. (value_b /= filt_value_b)
                case('==')
                    filter_particle = filter_particle .and. (value_b == filt_value_b)
                case('<')
                    filter_particle = filter_particle .and. (value_b < filt_value_b)
                case('<=')
                    filter_particle = filter_particle .and. (value_b <= filt_value_b)
                case('>')
                    filter_particle = filter_particle .and. (value_b > filt_value_b)
                case('>=')
                    filter_particle = filter_particle .and. (value_b >= filt_value_b)
                case default
                    write(*,*)'Relation operator not supported: ',TRIM(filt%cond_ops(i))
                    write(*,*)'Aborting!'
                    stop
                end select
            end if
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