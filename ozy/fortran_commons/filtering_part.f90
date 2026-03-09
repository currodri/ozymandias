module filtering_part
    use local
    use io_ramses
    use part_commons

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
            else
                write(*,*)'Variable type not supported in filtering: ',filt%cond_vars(i)%name,filt%cond_vars(i)%vartype
                write(*,*)'Aborting!'
                stop
            end if
            if (filt%use_var(i)) then
                filt%cond_vars_comp(i)%name = filt%cond_vars_comp_name(i)
                call set_part_var(vardict,vtypedict,filt%cond_vars_comp(i))
            end if
        end do ! i
    end subroutine get_filter_part_tools

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
            else
                write(*,*)'Variable type not supported in filtering: ',filt%cond_vtype(i)
                write(*,*)'Aborting!'
                stop
            end if
        end do

        
    end function filter_particle

end module filtering_part