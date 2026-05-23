!--------------------------------------------------------------------------
! ozymandias:variables.f90
!--------------------------------------------------------------------------
!
! MODULE: variables
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> types and functions that allow the handling of hydro raw variables and
!> derived variables
!
!> @details  
!> 
! 
!
!> @date 08/5/2023   0.3 making ozymandias fully compatible with non-cosmo
!--------------------------------------------------------------------------
module part_commons
    use local
    use constants
    use dictionary_commons
    use vectors
    use basis_representations
    use coordinate_systems
    use geometrical_regions
    use io_ramses, only: amr_info, sim_info, partIDs, partvar_types

    type part_var
        character(128) :: name,type
        integer :: vartype
        integer, dimension(:), allocatable :: ids,vtypes
        real(dbl) :: num_suffix
        procedure(myinterface_d),pointer,nopass :: myfunction_d
        procedure(myinterface_i),pointer,nopass :: myfunction_i
        procedure(myinterface_b),pointer,nopass :: myfunction_b
    end type part_var

    contains
    ! BINNING FUNCTIONS
    subroutine findbinpos_part(my_amr,my_sim,reg,dcell,part_data_d,&
                            & part_data_i,part_data_b,&
                            & ibin,value,trans_matrix,&
                            & scaletype,nbins,bins,linthresh,&
                            & zero_index,xvar)
        use vectors
        use geometrical_regions
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dcell
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_data_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_data_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_data_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_data_b
        integer,intent(inout) :: ibin
        real(dbl),intent(inout) :: value
#ifdef LONGINT
        integer(ilg) :: value_i
#else
        integer(irg) :: value_i
#endif
        integer(1) :: value_b
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix
        character(128),intent(in) :: scaletype
        integer,intent(in) :: nbins
        real(dbl),intent(in) :: linthresh
        integer,intent(in) :: zero_index
        real(dbl) :: origvalue
        real(dbl),dimension(0:nbins) :: bins
        type(part_var),intent(in) :: xvar

        ! Get variable value
        if (xvar%vartype==1) then
            value = xvar%myfunction_d(my_amr,my_sim,xvar,reg,dcell,part_data_d,part_data_i,part_data_b)
        else if (xvar%vartype==2) then
            value_i = xvar%myfunction_i(my_amr,my_sim,xvar,reg,dcell,part_data_d,part_data_i,part_data_b)
            value = real(value_i,kind=dbl)
        else if (xvar%vartype==3) then
            value_b = xvar%myfunction_b(my_amr,my_sim,xvar,reg,dcell,part_data_d,part_data_i,part_data_b)
            value = real(value_b,kind=dbl)
        else
            write(*,*)'ERROR: Unknown variable type in findbinpos_part'
            stop
        end if
        origvalue = value

        ! Make sure we are not out of the outer boundaries
        if (trim(scaletype).eq.'log_even') then
            value = log10(value)
        end if
        if (value .eq. bins(nbins)) then
            ibin = nbins
            value = origvalue
            return
        else if (value<bins(0).or.value>bins(nbins)) then
            ibin = 0
            value = origvalue
            return
        end if

        value = origvalue

        ! Transform value depending how the bins are provided
        if (trim(scaletype).eq.'log_even') then
            value = log10(value)
            ibin = int(dble(nbins)*(value-bins(0))/(bins(nbins)-bins(0))) + 1
        else if (trim(scaletype).eq.'linear_even') then
            ibin = int(dble(nbins)*(value-bins(0))/(bins(nbins)-bins(0))) + 1
        else if (trim(scaletype).eq.'symlog') then
            if (value <= -linthresh) then
                ! Negative logarithmic region
                value = log10(-value)
                ibin = -int(dble(zero_index-1) * (value - log10(-bins(0))) / (log10(-bins(0)) - log10(-bins(zero_index-1)))) + 1
            elseif (value > linthresh) then
                ! Positive logarithmic region
                value = log10(value)
                ibin = int(dble(nbins - zero_index + 1) * (value - log10(bins(zero_index))) / (log10(bins(nbins)) - log10(bins(zero_index)))) + zero_index
            else
                ! Linear region around zero
                ibin = zero_index
            endif
        else
            ibin = 1
            do while (ibin .lt. nbins)
                if (value .le. bins(ibin)) exit
                ibin = ibin + 1
            end do
        end if
        
        ! Last check
        if (ibin < 0) ibin = 0

        value = origvalue
    end subroutine findbinpos_part

    ! RAW VARIABLES
    function raw_part_d(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: raw_part_d

        raw_part_d = part_var_d(pvar%ids(1))
    end function raw_part_d

    function raw_part_i(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

#ifdef LONGINT
        integer(ilg) :: raw_part_i
#else
        integer(irg) :: raw_part_i
#endif
        raw_part_i = part_var_i(pvar%ids(1))
    end function raw_part_i

    function raw_part_b(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        integer(1) :: raw_part_b

        raw_part_b = part_var_b(pvar%ids(1))
    end function raw_part_b

    ! GEOMETRICAL VARIABLES
    function myinterface_d(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: myinterface_d

        myinterface_d = part_var_d(pvar%ids(1))
    end function myinterface_d

    function myinterface_i(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

#ifdef LONGINT
        integer(ilg) :: myinterface_i
#else
        integer(irg) :: myinterface_i
#endif
        myinterface_i = part_var_i(pvar%ids(1))
    end function myinterface_i

    function myinterface_b(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        integer(1) :: myinterface_b

        myinterface_b = part_var_b(pvar%ids(1))
    end function myinterface_b

    function d_euclid_wrap(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: d_euclid_wrap
        type(vector) :: pos

        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))
        d_euclid_wrap = magnitude(pos)
    end function d_euclid_wrap

    function r_sphere_wrap(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: r_sphere_wrap
        type(vector) :: pos

        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        ! Radius from center of the sphere
        r_sphere_wrap = r_sphere(pos)
    end function r_sphere_wrap

    function theta_sphere_wrap(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: theta_sphere_wrap
        type(vector) :: pos

        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        ! Value of spherical theta angle measured from the z axis
        theta_sphere_wrap = theta_sphere(pos)
    end function theta_sphere_wrap

    function phi_sphere_wrap(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: phi_sphere_wrap
        type(vector) :: pos

        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        ! Value of spherical phi angle measure in the x-y plane 
        ! from the x axis
        phi_sphere_wrap = phi_sphere(pos)
    end function phi_sphere_wrap

    function r_cyl_wrap(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else   
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: r_cyl_wrap
        type(vector) :: pos

        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        ! Radius from the z axis
        r_cyl_wrap = r_cyl(pos)
    end function r_cyl_wrap

    function phi_cyl_wrap(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT 
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: phi_cyl_wrap
        type(vector) :: pos

        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        ! Value of cylindrical phi angle measure in the x-y plane 
        ! from the x axis
        phi_cyl_wrap = phi_cyl(pos)
    end function phi_cyl_wrap

    subroutine check_geovar(vardict,varname,pvar,ok)
        implicit none

        type(dictf90), intent(in)      :: vardict
        character(128),intent(in) :: varname
        type(part_var),intent(inout) :: pvar
        logical,intent(inout) :: ok

        ok = .true.

        select case (trim(varname))
        case ('d_euclid')
            ! Euclidean distance
            pvar%type = 'geometric'
            pvar%name = 'd_euclid'
            pvar%vartype = 1
            allocate(pvar%ids(3))
            allocate(pvar%vtypes(3))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%myfunction_d => d_euclid_wrap
        case ('r_sphere')
            ! Radius from center of sphere
            pvar%type = 'geometric'
            pvar%name = 'r_sphere'
            pvar%vartype = 1
            allocate(pvar%ids(3))
            allocate(pvar%vtypes(3))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%myfunction_d => r_sphere_wrap
        case ('theta_sphere')
            ! Spherical theta angle
            pvar%type = 'geometric'
            pvar%name = 'theta_sphere'
            pvar%vartype = 1
            allocate(pvar%ids(3))
            allocate(pvar%vtypes(3))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%myfunction_d => theta_sphere_wrap
        case ('phi_sphere')
            ! Spherical phi angle
            pvar%type = 'geometric'
            pvar%name = 'phi_sphere'
            pvar%vartype = 1
            allocate(pvar%ids(3))
            allocate(pvar%vtypes(3))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%myfunction_d => phi_sphere_wrap
        case ('r_cyl')
            ! Radius from z axis
            pvar%type = 'geometric'
            pvar%name = 'r_cyl'
            pvar%vartype = 1
            allocate(pvar%ids(3))
            allocate(pvar%vtypes(3))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%myfunction_d => r_cyl_wrap
        case ('phi_cyl')
            ! Cylindrical phi angle
            pvar%type = 'geometric'
            pvar%name = 'phi_cyl'
            pvar%vartype = 1
            allocate(pvar%ids(3))
            allocate(pvar%vtypes(3))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%myfunction_d => phi_cyl_wrap
        case default
            ok = .false.
        end select
    end subroutine

    ! DERIVED VARIABLES
    function v_sphere_r(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: v_sphere_r
        type(vector) :: pos,v
        type(basis) :: temp_basis

        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        ! Velocity component in the spherical radial direction
        ! Dot product of velocity vector with spherical radial
        !    unit vector
        call spherical_basis_from_cartesian(pos,temp_basis)
        v_sphere_r = v.DOT.temp_basis%u(1)
    end function v_sphere_r

    function v_sphere_phi(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: v_sphere_phi
        type(vector) :: pos,v
        type(basis) :: temp_basis

        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        ! Velocity component in the spherical phi direction
        ! Dot product of velocity vector with spherical phi
        !    unit vector
        call spherical_basis_from_cartesian(pos,temp_basis)
        v_sphere_phi = v.DOT.temp_basis%u(3)
    end function v_sphere_phi

    function v_sphere_theta(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: v_sphere_theta
        type(vector) :: pos,v
        type(basis) :: temp_basis

        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        ! Velocity component in the spherical theta direction
        ! Dot product of velocity vector with spherical theta
        !    unit vector
        call spherical_basis_from_cartesian(pos,temp_basis)
        v_sphere_theta = v.DOT.temp_basis%u(2)
    end function v_sphere_theta

    function v_cyl_z(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: v_cyl_z
        type(vector) :: pos,v
        type(basis) :: temp_basis

        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        ! Velocity component in the cylindrical z direction
        ! Dot product of velocity vector with cylindrical z
        !    unit vector
        call cylindrical_basis_from_cartesian(pos,temp_basis)
        v_cyl_z = v.DOT.temp_basis%u(3)
    end function v_cyl_z

    function v_cyl_phi(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: v_cyl_phi
        type(vector) :: pos,v
        type(basis) :: temp_basis

        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        ! Velocity component in the cylindrical phi direction
        ! Dot product of velocity vector with cylindrical phi
        !    unit vector
        call cylindrical_basis_from_cartesian(pos,temp_basis)
        v_cyl_phi = v.DOT.temp_basis%u(2)
    end function v_cyl_phi

    function v_cyl_r(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: v_cyl_r
        type(vector) :: pos,v
        type(basis) :: temp_basis

        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        ! Velocity component in the cylindrical r direction
        ! Dot product of velocity vector with cylindrical r
        !    unit vector
        call cylindrical_basis_from_cartesian(pos,temp_basis)
        v_cyl_r = v.DOT.temp_basis%u(1)
    end function v_cyl_r

    function v_magnitude(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: v_magnitude
        type(vector) :: v

        v%x = part_var_d(pvar%ids(1))
        v%y = part_var_d(pvar%ids(2))
        v%z = part_var_d(pvar%ids(3))

        ! Magnitude of the velocity vector
        v_magnitude = magnitude(v)
    end function v_magnitude

    function v_squared(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif 
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: v_squared
        type(vector) :: v

        v%x = part_var_d(pvar%ids(1))
        v%y = part_var_d(pvar%ids(2))
        v%z = part_var_d(pvar%ids(3))

        ! Square of the magnitude of the velocity vector
        v_squared = v.DOT.v
    end function v_squared

    function v_tangential(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif 
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: v_tangential
        type(vector) :: pos,v
        type(basis) :: temp_basis

        ! Tangential velocity magnitude
        ! Consider as if one substracts the radial velocity,
        ! the remaining component is the tangential velocity
        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        call spherical_basis_from_cartesian(pos,temp_basis)
        v_tangential = sqrt((v .DOT. v) - (v .DOT. temp_basis%u(1))**2d0)
    end function v_tangential

    function centripetal_acc(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif 
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: centripetal_acc
        type(vector) :: pos,v
        type(basis) :: temp_basis

        ! Centripetal acceleration
        ! The centripetal acceleration is the acceleration
        !    required to keep an object moving in a circle
        !    at a constant speed
        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        call spherical_basis_from_cartesian(pos,temp_basis)
        centripetal_acc = ((v.DOT.v) - (v.DOT.temp_basis%u(1))**2)/r_sphere(pos)
    end function centripetal_acc

    function momentum_x(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif 
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: momentum_x
        type(vector) :: v
        real(dbl) :: m
        
        v%x = part_var_d(pvar%ids(1))
        v%y = part_var_d(pvar%ids(2))
        v%z = part_var_d(pvar%ids(3))
        m = part_var_d(pvar%ids(4))

        ! x-component of momentum
        momentum_x = m*v%x
    end function momentum_x

    function momentum_y(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else 
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: momentum_y
        type(vector) :: v
        real(dbl) :: m
        
        v%x = part_var_d(pvar%ids(1))
        v%y = part_var_d(pvar%ids(2))
        v%z = part_var_d(pvar%ids(3))
        m = part_var_d(pvar%ids(4))

        ! y-component of momentum
        momentum_y = m*v%y
    end function momentum_y

    function momentum_z(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: momentum_z
        type(vector) :: v
        real(dbl) :: m
        
        v%x = part_var_d(pvar%ids(1))
        v%y = part_var_d(pvar%ids(2))
        v%z = part_var_d(pvar%ids(3))
        m = part_var_d(pvar%ids(4))

        ! z-component of momentum
        momentum_z = m*v%z
    end function momentum_z

    function momentum(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif 
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: momentum
        type(vector) :: v
        real(dbl) :: m

        v%x = part_var_d(pvar%ids(1))
        v%y = part_var_d(pvar%ids(2))
        v%z = part_var_d(pvar%ids(3))
        m = part_var_d(pvar%ids(4))

        ! Magnitude of momentum
        momentum = m*magnitude(v)
    end function momentum

    function momentum_sphere_r(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif 
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b
    
        real(dbl) :: momentum_sphere_r
        type(vector) :: pos,v
        type(basis) :: temp_basis

        ! Linear momentum in the spherical radial direction
        ! 1. Dot product of velocity vector with spherical r
        !    unit vector
        ! 2. Multiply by mass of particle
        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        call spherical_basis_from_cartesian(pos,temp_basis)
        momentum_sphere_r = part_var_d(pvar%ids(7))*(v.DOT.temp_basis%u(1))
    end function momentum_sphere_r

    function momentum_cyl_z(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif 
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: momentum_cyl_z
        type(vector) :: pos,v
        type(basis) :: temp_basis

        ! Linear momentum in the cylindrical z direction
        ! 1. Dot product of velocity vector with cylindrical z
        !    unit vector
        ! 2. Multiply by mass of particle
        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        call cylindrical_basis_from_cartesian(pos,temp_basis)
        momentum_cyl_z = part_var_d(pvar%ids(7))*(v.DOT.temp_basis%u(3))
    end function momentum_cyl_z

    function ang_momentum_x(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif 
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: ang_momentum_x
        type(vector) :: pos,v
        real(dbl) :: m

        ! Angular momentum x-component
        ! 1. Cross product of position vector with linear momentum
        ! 2. Multiply by mass of particle
        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        m = part_var_d(pvar%ids(7))

        ang_momentum_x = m*(pos%y*v%z - pos%z*v%y)
    end function ang_momentum_x

    function ang_momentum_y(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else 
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: ang_momentum_y
        type(vector) :: pos,v
        real(dbl) :: m

        ! Angular momentum y-component
        ! 1. Cross product of position vector with linear momentum
        ! 2. Multiply by mass of particle
        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        m = part_var_d(pvar%ids(7))

        ang_momentum_y = m*(pos%z*v%x - pos%x*v%z)
    end function ang_momentum_y

    function ang_momentum_z(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: ang_momentum_z
        type(vector) :: pos,v
        real(dbl) :: m

        ! Angular momentum z-component
        ! 1. Cross product of position vector with linear momentum
        ! 2. Multiply by mass of particle
        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        m = part_var_d(pvar%ids(7))

        ang_momentum_z = m*(pos%x*v%y - pos%y*v%x)
    end function ang_momentum_z

    function ang_momentum(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: ang_momentum
        type(vector) :: pos,v
        real(dbl) :: m

        ! Magnitude of angular momentum
        ! 1. Cross product of position vector with linear momentum
        ! 2. Multiply by mass of particle
        pos%x = part_var_d(pvar%ids(1))
        pos%y = part_var_d(pvar%ids(2))
        pos%z = part_var_d(pvar%ids(3))

        v%x = part_var_d(pvar%ids(4))
        v%y = part_var_d(pvar%ids(5))
        v%z = part_var_d(pvar%ids(6))

        m = part_var_d(pvar%ids(7))

        ang_momentum = m*magnitude(pos*v)
    end function ang_momentum

    function density(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: density
        real(dbl) :: m

        ! Density of the particle
        ! Mass of the particle divided by the volume of the cell
        m = part_var_d(pvar%ids(1))
        density = m / (dx%x*dx%y*dx%z)
    end function density

    function sdensity(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif 
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: sdensity
        real(dbl) :: m

        ! Surface density of the particle
        ! Mass of the particle divided by the surface of the cell
        m = part_var_d(pvar%ids(1))
        sdensity = m / (dx%x*dx%y)
        ! if (m*my_sim%unit_m*g2msun<1980d0) print*,'mass,surface: ',m*my_sim%unit_m*g2msun,(dx%x*dx%y)*(my_sim%unit_l*cm2kpc*my_sim%boxlen)**2
    end function sdensity

    subroutine check_dervar(vardict,varname,pvar,ok)
        implicit none

        type(dictf90),intent(in) :: vardict
        character(128),intent(in) :: varname
        type(part_var),intent(inout) :: pvar
        logical,intent(out) :: ok

        ok = .true.
        select case(varname)
        case('v_sphere_r')
            ! Velocity component in the spherical radial direction
            pvar%type = 'derived'
            pvar%name = 'v_sphere_r'
            pvar%vartype = 1
            allocate(pvar%ids(6))
            allocate(pvar%vtypes(6))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%myfunction_d => v_sphere_r
        case('v_sphere_phi')
            ! Velocity component in the spherical phi direction
            pvar%type = 'derived'
            pvar%name = 'v_sphere_phi'
            pvar%vartype = 1
            allocate(pvar%ids(6))
            allocate(pvar%vtypes(6))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%myfunction_d => v_sphere_phi
        case('v_sphere_theta')
            ! Velocity component in the spherical theta direction
            pvar%type = 'derived'
            pvar%name = 'v_sphere_theta'
            pvar%vartype = 1
            allocate(pvar%ids(6))
            allocate(pvar%vtypes(6))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%myfunction_d => v_sphere_theta
        case('v_cyl_z')
            ! Velocity component in the cylindrical z direction
            pvar%type = 'derived'
            pvar%name = 'v_cyl_z'
            pvar%vartype = 1
            allocate(pvar%ids(6))
            allocate(pvar%vtypes(6))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%myfunction_d => v_cyl_z
        case('v_cyl_phi')
            ! Velocity component in the cylindrical phi direction
            pvar%type = 'derived'
            pvar%name = 'v_cyl_phi'
            pvar%vartype = 1
            allocate(pvar%ids(6))
            allocate(pvar%vtypes(6))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%myfunction_d => v_cyl_phi
        case('v_cyl_r')
            ! Velocity component in the cylindrical r direction
            pvar%type = 'derived'
            pvar%name = 'v_cyl_r'
            pvar%vartype = 1
            allocate(pvar%ids(6))
            allocate(pvar%vtypes(6))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%myfunction_d => v_cyl_r
        case('v_magnitude')
            ! Magnitude of the velocity vector
            pvar%type = 'derived'
            pvar%name = 'v_magnitude'
            pvar%vartype = 1
            allocate(pvar%ids(3))
            allocate(pvar%vtypes(3))
            pvar%ids(1) = vardict%get('velocity_x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('velocity_y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('velocity_z')
            pvar%vtypes(3) = 1
            pvar%myfunction_d => v_magnitude
        case('v_squared')
            ! Square of the magnitude of the velocity vector
            pvar%type = 'derived'
            pvar%name = 'v_squared'
            pvar%vartype = 1
            allocate(pvar%ids(3))
            allocate(pvar%vtypes(3))
            pvar%ids(1) = vardict%get('velocity_x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('velocity_y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('velocity_z')
            pvar%vtypes(3) = 1
            pvar%myfunction_d => v_squared
        case('v_tangential')
            ! Tangential velocity magnitude
            pvar%type = 'derived'
            pvar%name = 'v_tangential'
            pvar%vartype = 1
            allocate(pvar%ids(6))
            allocate(pvar%vtypes(6))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%myfunction_d => v_tangential
        case('centripetal_acc')
            ! Centripetal acceleration
            pvar%type = 'derived'
            pvar%name = 'centripetal_acc'
            pvar%vartype = 1
            allocate(pvar%ids(7))
            allocate(pvar%vtypes(7))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%ids(7) = vardict%get('r_sphere')
            pvar%vtypes(7) = 1
            pvar%myfunction_d => centripetal_acc
        case('momentum_x')
            ! x-component of momentum
            pvar%type = 'derived'
            pvar%name = 'momentum_x'
            pvar%vartype = 1
            allocate(pvar%ids(4))
            allocate(pvar%vtypes(4))
            pvar%ids(1) = vardict%get('velocity_x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('velocity_y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('velocity_z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('mass')
            pvar%vtypes(4) = 1
            pvar%myfunction_d => momentum_x
        case('momentum_y')
            ! y-component of momentum
            pvar%type = 'derived'
            pvar%name = 'momentum_y'
            pvar%vartype = 1
            allocate(pvar%ids(4))
            allocate(pvar%vtypes(4))
            pvar%ids(1) = vardict%get('velocity_x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('velocity_y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('velocity_z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('mass')
            pvar%vtypes(4) = 1
            pvar%myfunction_d => momentum_y
        case('momentum_z')
            ! z-component of momentum
            pvar%type = 'derived'
            pvar%name = 'momentum_z'
            pvar%vartype = 1
            allocate(pvar%ids(4))
            allocate(pvar%vtypes(4))
            pvar%ids(1) = vardict%get('velocity_x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('velocity_y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('velocity_z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('mass')
            pvar%vtypes(4) = 1
            pvar%myfunction_d => momentum_z
        case('momentum')
            ! Magnitude of momentum
            pvar%type = 'derived'
            pvar%name = 'momentum'
            pvar%vartype = 1
            allocate(pvar%ids(4))
            allocate(pvar%vtypes(4))
            pvar%ids(1) = vardict%get('velocity_x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('velocity_y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('velocity_z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('mass')
            pvar%vtypes(4) = 1
            pvar%myfunction_d => momentum
        case('momentum_sphere_r')
            ! Linear momentum in the spherical radial direction
            pvar%type = 'derived'
            pvar%name = 'momentum_sphere_r'
            pvar%vartype = 1
            allocate(pvar%ids(7))
            allocate(pvar%vtypes(7))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%ids(7) = vardict%get('mass')
            pvar%vtypes(7) = 1
            pvar%myfunction_d => momentum_sphere_r
        case('momentum_cyl_z')
            ! Linear momentum in the cylindrical z direction
            pvar%type = 'derived'
            pvar%name = 'momentum_cyl_z'
            pvar%vartype = 1
            allocate(pvar%ids(7))
            allocate(pvar%vtypes(7))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%ids(7) = vardict%get('mass')
            pvar%vtypes(7) = 1
            pvar%myfunction_d => momentum_cyl_z
        case('ang_momentum_x')
            ! Angular momentum x-component
            pvar%type = 'derived'
            pvar%name = 'ang_momentum_x'
            pvar%vartype = 1
            allocate(pvar%ids(7))
            allocate(pvar%vtypes(7))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%ids(7) = vardict%get('mass')
            pvar%vtypes(7) = 1
            pvar%myfunction_d => ang_momentum_x
        case('ang_momentum_y')
            ! Angular momentum y-component
            pvar%type = 'derived'
            pvar%name = 'ang_momentum_y'
            pvar%vartype = 1
            allocate(pvar%ids(7))
            allocate(pvar%vtypes(7))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%ids(7) = vardict%get('mass')
            pvar%vtypes(7) = 1
            pvar%myfunction_d => ang_momentum_y
        case('ang_momentum_z')
            ! Angular momentum z-component
            pvar%type = 'derived'
            pvar%name = 'ang_momentum_z'
            pvar%vartype = 1
            allocate(pvar%ids(7))
            allocate(pvar%vtypes(7))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%ids(7) = vardict%get('mass')
            pvar%vtypes(7) = 1
            pvar%myfunction_d => ang_momentum_z
        case('ang_momentum')
            ! Magnitude of angular momentum
            pvar%type = 'derived'
            pvar%name = 'ang_momentum'
            pvar%vartype = 1
            allocate(pvar%ids(7))
            allocate(pvar%vtypes(7))
            pvar%ids(1) = vardict%get('x')
            pvar%vtypes(1) = 1
            pvar%ids(2) = vardict%get('y')
            pvar%vtypes(2) = 1
            pvar%ids(3) = vardict%get('z')
            pvar%vtypes(3) = 1
            pvar%ids(4) = vardict%get('velocity_x')
            pvar%vtypes(4) = 1
            pvar%ids(5) = vardict%get('velocity_y')
            pvar%vtypes(5) = 1
            pvar%ids(6) = vardict%get('velocity_z')
            pvar%vtypes(6) = 1
            pvar%ids(7) = vardict%get('mass')
            pvar%vtypes(7) = 1
            pvar%myfunction_d => ang_momentum
        case('density')
            ! Density of the particle
            pvar%type = 'derived'
            pvar%name = 'density'
            pvar%vartype = 1
            allocate(pvar%ids(1))
            allocate(pvar%vtypes(1))
            pvar%ids(1) = vardict%get('mass')
            pvar%vtypes(1) = 1
            pvar%myfunction_d => density
        case('sdensity')
            ! Surface density of the particle
            pvar%type = 'derived'
            pvar%name = 'sdensity'
            pvar%vartype = 1
            allocate(pvar%ids(1))
            allocate(pvar%vtypes(1))
            pvar%ids(1) = vardict%get('mass')
            pvar%vtypes(1) = 1
            pvar%myfunction_d => sdensity
        case default
            ok = .false.
        end select
    end subroutine check_dervar

    ! DERIVED STAR VARIABLES
    function age(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: age
        real(dbl) :: t
        integer :: iii

        if (part_var_d(pvar%ids(1)).eq.0d0) then
            age = 0d0
            return
        end if
        
        ! Age of the star
        ! Time of the simulation minus the time of the star formation
        if (my_sim%cosmo) then
            iii = 1
            do while(my_sim%tau_frw(iii)>part_var_d(pvar%ids(1)).and.iii<my_sim%n_frw)
                iii = iii + 1
            end do
            ! Interpolate time
#ifdef AGEPROPER
            t = part_var_d(pvar%ids(1))
#else
            if (my_sim%rt) then
                ! RT simulations always force use_proper_time=.true.
                t = part_var_d(pvar%ids(1))
            else
                t = my_sim%t_frw(iii)*(part_var_d(pvar%ids(1))-my_sim%tau_frw(iii-1))/(my_sim%tau_frw(iii)-my_sim%tau_frw(iii-1))+ &
                & my_sim%t_frw(iii-1)*(part_var_d(pvar%ids(1))-my_sim%tau_frw(iii))/(my_sim%tau_frw(iii-1)-my_sim%tau_frw(iii))
            end if
#endif
            age = (my_sim%time_simu-t)/(my_sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
        else
            t = part_var_d(pvar%ids(1))
            age = (my_sim%time_simu-t)*my_sim%unit_t/(365.*24.*3600.*1d9)
        end if

    end function age

    function birth_date(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: birth_date
        real(dbl) :: t,age
        integer :: iii

        ! Birth date of the star
        ! Time of the star formation
        if (my_sim%cosmo) then
            iii = 1
            do while(my_sim%tau_frw(iii)>part_var_d(pvar%ids(1)).and.iii<my_sim%n_frw)
                iii = iii + 1
            end do
            ! Interpolate time
#ifdef AGEPROPER
            t = part_var_d(pvar%ids(1))
#else
            if (my_sim%rt) then
                ! RT simulations always force use_proper_time=.true.
                t = part_var_d(pvar%ids(1))
            else
                t = my_sim%t_frw(iii)*(part_var_d(pvar%ids(1))-my_sim%tau_frw(iii-1))/(my_sim%tau_frw(iii)-my_sim%tau_frw(iii-1))+ &
                & my_sim%t_frw(iii-1)*(part_var_d(pvar%ids(1))-my_sim%tau_frw(iii))/(my_sim%tau_frw(iii-1)-my_sim%tau_frw(iii))
            end if
#endif
            age = (my_sim%time_simu-t)
            birth_date = (my_sim%time_tot+age)/(my_sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
        else
            birth_date = part_var_d(pvar%ids(1))*my_sim%unit_t/(365.*24.*3600.*1d9)
        end if

    end function birth_date

    function sfr(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: sfr
        real(dbl) :: t,birth_date
        integer :: iii
        real(dbl) :: sfrind,current_age_univ,m

        sfr = 0d0

        if (part_var_d(pvar%ids(1)).eq.0d0) return

        ! Star formation rate
        ! This quantity only makes sense on a cumulative sense,
        ! as the SFR per individual stellar particle is not well defined
        ! 1. Get the indicator desired (timescale over which SFR is averaged)
        sfrind = pvar%num_suffix

        ! Get the initial mass of the stellar particle
#ifdef IMASS
        m = part_var_d(pvar%ids(2))
#else
        m = part_var_d(pvar%ids(2)) / (1d0 - my_sim%eta_sn)
#endif

        ! 2. Get the particle birth date
        if (my_sim%cosmo) then
            iii = 1
            do while(my_sim%tau_frw(iii)>part_var_d(pvar%ids(1)).and.iii<my_sim%n_frw)
                iii = iii + 1
            end do
            ! Interpolate time
#ifdef AGEPROPER
            t = part_var_d(pvar%ids(1))
#else
            if (my_sim%rt) then
                ! RT simulations always force use_proper_time=.true.
                t = part_var_d(pvar%ids(1))
            else
                t = my_sim%t_frw(iii)*(part_var_d(pvar%ids(1))-my_sim%tau_frw(iii-1))/(my_sim%tau_frw(iii)-my_sim%tau_frw(iii-1))+ &
                & my_sim%t_frw(iii-1)*(part_var_d(pvar%ids(1))-my_sim%tau_frw(iii))/(my_sim%tau_frw(iii-1)-my_sim%tau_frw(iii))
            end if
#endif
            current_age_univ = (my_sim%time_tot+my_sim%time_simu)/(my_sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
            birth_date = (my_sim%time_tot+t)/(my_sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
            ! 3. If the particle is older than the indicator, the SFR is zero
            if (birth_date >= (current_age_univ - sfrind)) then
                sfr = part_var_d(pvar%ids(2))
            else
                sfr = 0d0
            end if
        else
            birth_date = part_var_d(pvar%ids(1))*my_sim%unit_t/(365.*24.*3600.*1d9)
            current_age_univ = my_sim%time_simu*my_sim%unit_t/(365.*24.*3600.*1d9)
            ! 3. If the particle is older than the indicator, the SFR is zero
            if (birth_date >= (current_age_univ - sfrind)) then
                sfr = part_var_d(pvar%ids(2))
            else
                sfr = 0d0
            end if
        end if
    end function sfr

    function sfr_surface(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: sfr_surface
        real(dbl) :: t,birth_date
        integer :: iii
        real(dbl) :: sfrind,current_age_univ,m

        sfr_surface = 0d0

        if (part_var_d(pvar%ids(1)).eq.0d0) return

        ! Star formation rate
        ! This quantity only makes sense on a cumlative sense,
        ! as the SFR per individual stellar particle is not well defined
        ! 1. Get the indicator desired (timescale over which SFR is averaged)
        sfrind = pvar%num_suffix
    
        ! Get the initial mass of the stellar particle
#ifdef IMASS
        m = part_var_d(pvar%ids(2))
#else
        m = part_var_d(pvar%ids(2)) / (1d0 - my_sim%eta_sn)
#endif

        ! 2. Get the particle birth date
        if (my_sim%cosmo) then
            iii = 1
            do while(my_sim%tau_frw(iii)>part_var_d(pvar%ids(1)).and.iii<my_sim%n_frw)
                iii = iii + 1
            end do
            ! Interpolate time
#ifdef AGEPROPER
            t = part_var_d(pvar%ids(1))
#else
            if (my_sim%rt) then
                ! RT simulations always force use_proper_time=.true.
                t = part_var_d(pvar%ids(1))
            else
                t = my_sim%t_frw(iii)*(part_var_d(pvar%ids(1))-my_sim%tau_frw(iii-1))/(my_sim%tau_frw(iii)-my_sim%tau_frw(iii-1))+ &
                & my_sim%t_frw(iii-1)*(part_var_d(pvar%ids(1))-my_sim%tau_frw(iii))/(my_sim%tau_frw(iii-1)-my_sim%tau_frw(iii))
            end if
#endif
            current_age_univ = (my_sim%time_tot+my_sim%time_simu)/(my_sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
            birth_date = (my_sim%time_tot+t)/(my_sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
            ! 3. If the particle is older than the indicator, the SFR is zero
            if (birth_date >= (current_age_univ - sfrind)) then
                sfr_surface = part_var_d(pvar%ids(2)) / (dx%x*dx%y) / sfrind
            else
                sfr_surface = 0d0
            end if
        else
            birth_date = part_var_d(pvar%ids(1))*my_sim%unit_t/(365.*24.*3600.*1d9)
            current_age_univ = my_sim%time_simu*my_sim%unit_t/(365.*24.*3600.*1d9)
            ! 3. If the particle is older than the indicator, the SFR is zero
            if (birth_date >= (current_age_univ - sfrind)) then
                sfrind = sfrind / my_sim%unit_t * (365.*24.*3600.*1d9)
                sfr_surface = part_var_d(pvar%ids(2)) / (dx%x*dx%y) / sfrind
                ! if (part_var_d(pvar%ids(2))*my_sim%unit_m*g2msun<1980) then
                !     print*,'mass,surface,sfrind:',part_var_d(pvar%ids(2))*my_sim%unit_m*g2msun,dx%x*dx%y*(my_sim%unit_l*my_sim%boxlen)**2*cm2kpc**2,sfrind*my_sim%unit_t/(365.*24.*3600.*1d9)
                !     print*,'sfr_surface:',sfr_surface*my_sim%unit_m/((my_sim%unit_l*my_sim%boxlen)**2)/my_sim%unit_t*gscm22msunyrkpc2
                ! endif
            else
                sfr_surface = 0d0
            end if
        end if
    end function sfr_surface

    function sfr_density(my_amr,my_sim,pvar,reg,dx,part_var_d,part_var_i,part_var_b)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(part_var),intent(in) :: pvar
        type(region),intent(in) :: reg
        type(vector),intent(in) :: dx
        real(dbl),dimension(1:my_sim%nvar_part_d),intent(in) :: part_var_d
#ifdef LONGINT
        integer(ilg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#else
        integer(irg),dimension(1:my_sim%nvar_part_i),intent(in) :: part_var_i
#endif
        integer(1),dimension(1:my_sim%nvar_part_b),intent(in) :: part_var_b

        real(dbl) :: sfr_density
        real(dbl) :: t,birth_date
        integer :: iii
        real(dbl) :: sfrind,current_age_univ,m

        sfr_density = 0d0

        if (part_var_d(pvar%ids(1)).eq.0d0) return

        ! Star formation rate
        ! This quantity only makes sense on a cumlative sense,
        ! as the SFR per individual stellar particle is not well defined
        ! 1. Get the indicator desired (timescale over which SFR is averaged)
        sfrind = pvar%num_suffix

        ! Get the initial mass of the stellar particle
#ifdef IMASS
        m = part_var_d(pvar%ids(2))
#else
        m = part_var_d(pvar%ids(2)) / (1d0 - my_sim%eta_sn)
#endif

        ! 2. Get the particle birth date
        if (my_sim%cosmo) then
            iii = 1
            do while(my_sim%tau_frw(iii)>part_var_d(pvar%ids(1)).and.iii<my_sim%n_frw)
                iii = iii + 1
            end do
            ! Interpolate time
#ifdef AGEPROPER
            t = part_var_d(pvar%ids(1))
#else
            if (my_sim%rt) then
                ! RT simulations always force use_proper_time=.true.
                t = part_var_d(pvar%ids(1))
            else
                t = my_sim%t_frw(iii)*(part_var_d(pvar%ids(1))-my_sim%tau_frw(iii-1))/(my_sim%tau_frw(iii)-my_sim%tau_frw(iii-1))+ &
                & my_sim%t_frw(iii-1)*(part_var_d(pvar%ids(1))-my_sim%tau_frw(iii))/(my_sim%tau_frw(iii-1)-my_sim%tau_frw(iii))
            end if
#endif
            current_age_univ = (my_sim%time_tot+my_sim%time_simu)/(my_sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
            birth_date = (my_sim%time_tot+t)/(my_sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
            ! 3. If the particle is older than the indicator, the SFR is zero
            if (birth_date >= (current_age_univ - sfrind)) then
                sfr_density = part_var_d(pvar%ids(2)) / (dx%x*dx%y*dx%z) / sfrind
            else
                sfr_density = 0d0
            end if
        else
            birth_date = part_var_d(pvar%ids(1))
            current_age_univ = my_sim%time_simu*my_sim%unit_t/(365.*24.*3600.*1d9)
            ! 3. If the particle is older than the indicator, the SFR is zero
            if (birth_date >= (current_age_univ - sfrind)) then
                sfr_density = part_var_d(pvar%ids(2)) / (dx%x*dx%y*dx%z) / sfrind
            else
                sfr_density = 0d0
            end if
        end if
    end function sfr_density

    subroutine check_starvar(vardict,varname,pvar,ok)
        use utils, only: get_cleaned_string,get_numeric_suffix
        implicit none
        type(part_var),intent(inout) :: pvar
        type(dictf90),intent(in) :: vardict
        character(128),intent(in) :: varname
        logical,intent(out) :: ok

        character(128) :: clean_name
        integer :: index2,index3

        ok = .true.

        clean_name = get_cleaned_string(varname)

        select case(clean_name)
        case('age')
            ! Age of the star
            pvar%type = 'derived'
            pvar%name = 'age'
            pvar%vartype = 1
            allocate(pvar%ids(1))
            allocate(pvar%vtypes(1))
            pvar%ids(1) = vardict%get('birth_time')
            pvar%vtypes(1) = 1
            pvar%myfunction_d => age
        case('birth_date')
            ! Birth date of the star
            pvar%type = 'derived'
            pvar%name = 'birth_date'
            pvar%vartype = 1
            allocate(pvar%ids(1))
            allocate(pvar%vtypes(1))
            pvar%ids(1) = vardict%get('birth_time')
            pvar%vtypes(1) = 1
            pvar%myfunction_d => birth_date
        case('sfr')
            ! Star formation rate
            pvar%type = 'derived'
            pvar%name = varname
            pvar%vartype = 1
            allocate(pvar%ids(2))
            allocate(pvar%vtypes(2))
            pvar%ids(1) = vardict%get('birth_time')
            pvar%vtypes(1) = 1
#ifdef IMASS
            pvar%ids(2) = vardict%get('initial_mass')
#else
            pvar%ids(2) = vardict%get('mass')
#endif
            pvar%vtypes(2) = 1
            pvar%myfunction_d => sfr
            pvar%num_suffix = dble(get_numeric_suffix(varname))/1D3
        case('sfr_surface')
            ! Star formation rate per unit area
            pvar%type = 'derived'
            pvar%name = varname
            pvar%vartype = 1
            allocate(pvar%ids(2))
            allocate(pvar%vtypes(2))
            pvar%ids(1) = vardict%get('birth_time')
            pvar%vtypes(1) = 1
#ifdef IMASS
            pvar%ids(2) = vardict%get('initial_mass')
#else
            pvar%ids(2) = vardict%get('mass')
#endif
            pvar%vtypes(2) = 1
            pvar%myfunction_d => sfr_surface
            pvar%num_suffix = dble(get_numeric_suffix(varname))/1D3
        case('sfr_density')
            ! Star formation rate
            pvar%type = 'derived'
            pvar%name = varname
            pvar%vartype = 1
            allocate(pvar%ids(2))
            allocate(pvar%vtypes(2))
            pvar%ids(1) = vardict%get('birth_time')
            pvar%vtypes(1) = 1
#ifdef IMASS
            pvar%ids(2) = vardict%get('initial_mass')
#else
            pvar%ids(2) = vardict%get('mass')
#endif
            pvar%vtypes(2) = 1
            pvar%myfunction_d => sfr_density
            pvar%num_suffix = dble(get_numeric_suffix(varname))/1D3
        case default
            ok = .false.
        end select
    end subroutine check_starvar

    subroutine get_partvar_tools(vardict,vtypedict,nreq,reqvars,cleaned_vars)
        implicit none

        type(dictf90), intent(in) :: vardict,vtypedict
        integer, intent(in) :: nreq
        character(128),dimension(:),intent(in) :: reqvars
        type(part_var),dimension(:),intent(out) :: cleaned_vars

        logical :: ok_check
        integer :: i, ivar

        ! Loop over the requested variables
        do i = 1, nreq
            ivar = vardict%get(reqvars(i))
            if (ivar.ne.0) then
                ! 1. If variable is a raw variable defined in the
                ! particle dictionary, we just simply store it
                cleaned_vars(i)%type = 'raw'
                cleaned_vars(i)%name = reqvars(i)
                allocate(cleaned_vars(i)%ids(1))
                allocate(cleaned_vars(i)%vtypes(1))
                cleaned_vars(i)%ids(1) = ivar
                cleaned_vars(i)%vtypes(1) = vtypedict%get(reqvars(i))
                cleaned_vars(i)%vartype = vtypedict%get(reqvars(i))
                if (cleaned_vars(i)%vtypes(1).eq.1) then
                    ! Double float variable
                    cleaned_vars(i)%myfunction_d => raw_part_d
                else if (cleaned_vars(i)%vtypes(1).eq.2) then
                    ! Integer variable
                    cleaned_vars(i)%myfunction_i => raw_part_i
                elseif (cleaned_vars(i)%vtypes(1).eq.3) then
                    ! Single integer variable
                    cleaned_vars(i)%myfunction_b => raw_part_b
                end if
            else if (trim(reqvars(i)) == 'cumulative') then
                ! 2. Cumulative (just adding up counts)
                cleaned_vars(i)%type = 'cumulative'
                cleaned_vars(i)%name = 'cumulative'
                allocate(cleaned_vars(i)%ids(1))
                allocate(cleaned_vars(i)%vtypes(1))
                cleaned_vars(i)%ids(1) = 0
                cleaned_vars(i)%vtypes(1) = 1
                cleaned_vars(i)%vartype = 1
            else
                ! 3. Variable is not a raw variable, so it is either a
                ! geometrical, derived or star derived variable
                call check_geovar(vardict,reqvars(i),cleaned_vars(i),ok_check)
                if (ok_check) cycle

                ! Now check for a derived particle variable
                call check_dervar(vardict,reqvars(i),cleaned_vars(i),ok_check)
                if (ok_check) cycle

                ! And finally check for derived star variables
                call check_starvar(vardict,reqvars(i),cleaned_vars(i),ok_check)
                if (.not.ok_check) then
                    write(*,*) 'ERROR: Variable ',trim(reqvars(i)),' not found in particle dictionary'
                    write(*,*) 'vardict: ',vardict%keys
                    stop
                end if
            end if
        end do

    end subroutine get_partvar_tools

    subroutine set_part_var(vardict,vtypedict,pvar)
        implicit none

        type(dictf90), intent(in) :: vardict,vtypedict
        type(part_var),intent(inout) :: pvar

        logical :: ok_check
        integer :: i, ivar

        ivar = vardict%get(pvar%name)

        if (ivar.ne.0) then
            ! 1. If variable is a raw variable defined in the
            ! particle dictionary, we just simply store it
            pvar%type = 'raw'
            pvar%vartype = vtypedict%get(pvar%name)
            allocate(pvar%ids(1))
            allocate(pvar%vtypes(1))
            pvar%ids(1) = ivar
            pvar%vtypes(1) = vtypedict%get(pvar%name)
            if (pvar%vtypes(1).eq.1) then
                ! Double float variable
                pvar%myfunction_d => raw_part_d
            else if (pvar%vtypes(1).eq.2) then
                ! Integer variable
                pvar%myfunction_i => raw_part_i
            elseif (pvar%vtypes(1).eq.3) then
                ! Single integer variable
                pvar%myfunction_b => raw_part_b
            end if
        else if (trim(pvar%name) == 'cumulative') then
            ! 2. Cumulative (just adding up counts)
            pvar%type = 'cumulative'
            pvar%vartype = 1
            allocate(pvar%ids(1))
            allocate(pvar%vtypes(1))
            pvar%ids(1) = 0
            pvar%vtypes(1) = 1
        else
            ! 3. Variable is not a raw variable, so it is either a
            ! geometrical, derived or star derived variable
            call check_geovar(vardict,pvar%name,pvar,ok_check)
            if (ok_check) return

            ! Now check for a derived particle variable
            call check_dervar(vardict,pvar%name,pvar,ok_check)
            if (ok_check) return

            ! And finally check for derived star variables
            call check_starvar(vardict,pvar%name,pvar,ok_check)
            
            if (.not.ok_check) then
                write(*,*) 'ERROR: Variable ',trim(pvar%name),' not found in particle dictionary'
                stop
            end if
        end if            

    end subroutine set_part_var
    
end module part_commons