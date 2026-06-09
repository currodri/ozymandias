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
module hydro_commons
    use local
    use constants
    use dictionary_commons
    use vectors
    use basis_representations
    use coordinate_systems
    use geometrical_regions
    use io_ramses, only: amr_info, sim_info, rt_info, Tmin, cV, lambda_crGH08
    use cooling_module

    type hydro_var
        character(128) :: name,type
        integer, dimension(:), allocatable :: ids
        procedure(myinterface),pointer,nopass :: myfunction
    end type hydro_var

    contains

    ! BINNING FUNCTIONS
    subroutine findbinpos(my_amr,my_sim,my_rt,reg,x,cvars,csons,csize,&
                            & ibin,value,trans_matrix,&
                            & scaletype,nbins,bins,linthresh,&
                            & zero_index,xvar,&
                            & gvars,rtvars)
        implicit none

        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(region),intent(in) :: reg
        type(vector),intent(in) :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: cvars
        integer,dimension(0:my_amr%twondim),intent(in) :: csons
        real(dbl),intent(in) :: csize
        integer,intent(inout) :: ibin
        real(dbl),intent(inout) :: value
        real(dbl),dimension(1:3,1:3),intent(in) :: trans_matrix
        character(128),intent(in) :: scaletype
        integer,intent(in) :: nbins
        real(dbl),intent(in) :: linthresh
        integer,intent(in) :: zero_index
        real(dbl) :: origvalue
        real(dbl),dimension(0:nbins),intent(in) :: bins
        type(hydro_var),intent(in) :: xvar
        real(dbl),dimension(0:,:),optional,intent(in) :: gvars
        real(dbl),dimension(0:,:),optional,intent(in) :: rtvars

        ! Get variable value
        if (present(gvars) .and. present(rtvars)) then
            value = xvar%myfunction(my_amr,my_sim,my_rt,xvar,reg,csize,x,cvars,csons,trans_matrix,gvars,rtvars)
        else if (present(gvars)) then
            value = xvar%myfunction(my_amr,my_sim,my_rt,xvar,reg,csize,x,cvars,csons,trans_matrix,gvars)
        else if (present(rtvars)) then
            value = xvar%myfunction(my_amr,my_sim,my_rt,xvar,reg,csize,x,cvars,csons,trans_matrix,rt_var=rtvars)
        else
            value = xvar%myfunction(my_amr,my_sim,my_rt,xvar,reg,csize,x,cvars,csons,trans_matrix)
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
    end subroutine findbinpos

    ! RAW VARIABLES
    function raw_hydro(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var
        real(dbl) :: raw_hydro
        raw_hydro = var(0,hvar%ids(1))
    end function raw_hydro

    ! GEOMETRICAL VARIABLES

    function myinterface(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: myinterface
        myinterface = magnitude(x)
    end function myinterface

    function d_euclid_wrap(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: d_euclid_wrap

        ! Euclidean distance
        d_euclid_wrap = magnitude(x)
    end function d_euclid_wrap

    function x_coord_wrap(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: x_coord_wrap

        ! x - coordinate
        x_coord_wrap = x%x
    end function x_coord_wrap

    function y_coord_wrap(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: y_coord_wrap

        ! y - coordinate
        y_coord_wrap = x%y
    end function y_coord_wrap

    function z_coord_wrap(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: z_coord_wrap

        ! z - coordinate
        z_coord_wrap = x%z
    end function z_coord_wrap

    function r_sphere_wrap(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: r_sphere_wrap

        ! Radius from center of sphere
        r_sphere_wrap = r_sphere(x)
    end function r_sphere_wrap

    function theta_sphere_wrap(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: theta_sphere_wrap

        ! Value of spherical theta angle measured from the z axis
        theta_sphere_wrap = theta_sphere(x)
    end function theta_sphere_wrap

    function phi_sphere_wrap(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: phi_sphere_wrap

        ! Value of spherical phi angle measure in the x-y plane 
        ! from the x axis
        phi_sphere_wrap = phi_sphere(x)
    end function phi_sphere_wrap

    function r_cyl_wrap(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: r_cyl_wrap

        ! Value of cylindrical radius
        r_cyl_wrap = r_cyl(x)
    end function r_cyl_wrap

    function phi_cyl_wrap(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: phi_cyl_wrap

        ! Value of cylindrical phi angle measure in the x-y plane 
        ! from the x axis
        phi_cyl_wrap = phi_cyl(x)
    end function phi_cyl_wrap

    subroutine check_geovar(varname,hvar,ok)
        implicit none

        character(128), intent(in) :: varname
        type(hydro_var), intent(inout) :: hvar
        logical, intent(inout)        :: ok
        
        ok = .true.

        select case (trim(varname))
        case ('d_euclid')
            ! Euclidean distance
            hvar%type = 'geometric'
            hvar%name = 'd_euclid'
            hvar%myfunction => d_euclid_wrap
        case ('x')
            ! x - coordinate
            hvar%type = 'geometric'
            hvar%name = 'x'
            hvar%myfunction => x_coord_wrap
        case ('y')
            ! y - coordinate
            hvar%type = 'geometric'
            hvar%name = 'y'
            hvar%myfunction => y_coord_wrap
        case ('z')
            ! z - coordinate
            hvar%type = 'geometric'
            hvar%name = 'z'
            hvar%myfunction => z_coord_wrap
        case ('r_sphere')
            ! Radius from center of sphere
            hvar%type = 'geometric'
            hvar%name = 'r_sphere'
            hvar%myfunction => r_sphere_wrap
        case ('theta_sphere')
            ! Value of spherical theta angle measured from the z axis
            hvar%type = 'geometric'
            hvar%name = 'theta_sphere'
            hvar%myfunction => theta_sphere_wrap
        case ('phi_sphere')
            ! Value of spherical phi angle measure in the x-y plane 
            ! from the x axis
            hvar%type = 'geometric'
            hvar%name = 'phi_sphere'
            hvar%myfunction => phi_sphere_wrap
        case('r_cyl')
            ! Value of cylindrical radius
            hvar%type = 'geometric'
            hvar%name = 'r_cyl'
            hvar%myfunction => r_cyl_wrap
        case ('phi_cyl')
            ! Value of spherical phi angle measure in the x-y plane 
            ! from the x axis
            hvar%type = 'geometric'
            hvar%name = 'phi_cyl'
            hvar%myfunction => phi_cyl_wrap
        case default
            ok = .false.
        end select
    end subroutine check_geovar

    ! GRAVITATIONAL VARIABLES

    function grav_potential(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_potential

        ! Gravitational potential
        grav_potential = grav_var(0,1)
    end function grav_potential

    function grav_gx(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_gx

        ! Gravitational acceleration in the x direction
        grav_gx = grav_var(0,2)
    end function grav_gx

    function grav_gy(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_gy

        ! Gravitational acceleration in the y direction
        grav_gy = grav_var(0,3)
    end function grav_gy

    function grav_gz(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_gz

        ! Gravitational acceleration in the z direction
        grav_gz = grav_var(0,4)
    end function grav_gz

    subroutine check_gravvar(varname,hvar,ok)
        implicit none

        character(128), intent(in) :: varname
        type(hydro_var), intent(inout) :: hvar
        logical, intent(inout)        :: ok
        
        ok = .true.
        select case (trim(varname))
        case ('grav_potential')
            ! Gravitational potential
            hvar%type = 'gravity'
            hvar%name = 'grav_potential'
            hvar%myfunction => grav_potential
        case ('grav_gx')
            ! Gravitational acceleration in the x direction
            hvar%type = 'gravity'
            hvar%name = 'grav_gx'
            hvar%myfunction => grav_gx
        case ('grav_gy')
            ! Gravitational acceleration in the y direction
            hvar%type = 'gravity'
            hvar%name = 'grav_gy'
            hvar%myfunction => grav_gy
        case ('grav_gz')
            ! Gravitational acceleration in the z direction
            hvar%type = 'gravity'
            hvar%name = 'grav_gz'
            hvar%myfunction => grav_gz
        case default
            ok = .false.
        end select

    end subroutine check_gravvar

    ! DERIVED VARIABLES

    function cell_volume(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: cell_volume

        cell_volume = (dx * dx) * dx
    end function cell_volume

    function cell_mass(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: cell_mass

        cell_mass = (var(0,hvar%ids(1)) * (dx * dx)) * dx
    end function cell_mass

    function v_sphere_r(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: v_sphere_r
        type(vector) :: v
        type(basis) :: temp_basis

        ! Velocity component in the spherical radial direction
        ! Dot product of velocity vector with spherical radial
        !    unit vector
        v = (/var(0,hvar%ids(1)),var(0,hvar%ids(2)),var(0,hvar%ids(3))/)
        call spherical_basis_from_cartesian(x,temp_basis)
        v_sphere_r = v.DOT.temp_basis%u(1)
    end function v_sphere_r

    function v_sphere_phi(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: v_sphere_phi
        type(vector) :: v
        type(basis) :: temp_basis

        ! Velocity component in the spherical azimutal (phi) direction
        ! Dot product of velocity vector with spherical phi
        !    unit vector
        v = (/var(0,hvar%ids(1)),var(0,hvar%ids(2)),var(0,hvar%ids(3))/)
        call spherical_basis_from_cartesian(x,temp_basis)
        v_sphere_phi = v .DOT. temp_basis%u(3)
    end function v_sphere_phi

    function v_sphere_theta(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: v_sphere_theta
        type(vector) :: v
        type(basis) :: temp_basis

        ! Velocity component in the spherical theta direction
        ! Dot product of velocity vector with spherical theta
        !    unit vector
        v = (/var(0,hvar%ids(1)),var(0,hvar%ids(2)),var(0,hvar%ids(3))/)
        call spherical_basis_from_cartesian(x,temp_basis)
        v_sphere_theta = v .DOT. temp_basis%u(2)
    end function v_sphere_theta

    function v_cyl_r(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: v_cyl_r
        type(vector) :: v
        type(basis) :: temp_basis

        ! Velocity component in the cylindrical radial direction
        ! Dot product of velocity vector with cylindrical
        !    radial unit vector
        v = (/var(0,hvar%ids(1)),var(0,hvar%ids(2)),var(0,hvar%ids(3))/)
        call cylindrical_basis_from_cartesian(x,temp_basis)
        v_cyl_r = v .DOT. temp_basis%u(1)
    end function v_cyl_r

    function v_cyl_z(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: v_cyl_z
        type(vector) :: v
        type(basis) :: temp_basis

        ! Velocity component in the cylindrical radial direction
        ! Dot product of velocity vector with cylindrical
        !    z unit vector
        v = (/var(0,hvar%ids(1)),var(0,hvar%ids(2)),var(0,hvar%ids(3))/)
        call cylindrical_basis_from_cartesian(x,temp_basis)
        v_cyl_z = v .DOT. temp_basis%u(3)
    end function v_cyl_z

    function v_cyl_phi(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: v_cyl_phi
        type(vector) :: v
        type(basis) :: temp_basis

        ! Velocity component in the cylindrical radial direction
        ! Dot product of velocity vector with cylindrical
        !    phi unit vector
        v = (/var(0,hvar%ids(1)),var(0,hvar%ids(2)),var(0,hvar%ids(3))/)
        call cylindrical_basis_from_cartesian(x,temp_basis)
        v_cyl_phi = v .DOT. temp_basis%u(2)
    end function v_cyl_phi

    function v_magnitude(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: v_magnitude
        type(vector) :: v
        type(basis) :: temp_basis

        ! Velocity magnitude from galaxy coordinates
        v = (/var(0,hvar%ids(1)),var(0,hvar%ids(2)),var(0,hvar%ids(3))/)
        v_magnitude = magnitude(v)
    end function v_magnitude

    function v_squared(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: v_squared
        type(vector) :: v
        type(basis) :: temp_basis

        ! Velocity magnitude squared from galaxy coordinates
        v = (/var(0,hvar%ids(1)),var(0,hvar%ids(2)),var(0,hvar%ids(3))/)
        v_squared = v.DOT.v
    end function v_squared

    function momentum_x(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: momentum_x
        type(vector) :: v
        type(basis) :: temp_basis

        ! Linear momentum in the x direction as density*volume*corrected_velocity_x
        momentum_x = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * var(0,hvar%ids(2))
    end function momentum_x

    function momentum_y(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: momentum_y
        type(vector) :: v
        type(basis) :: temp_basis

        ! Linear momentum in the y direction as density*volume*corrected_velocity_y
        momentum_y = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * var(0,hvar%ids(2))
    end function momentum_y

    function momentum_z(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: momentum_z
        type(vector) :: v
        type(basis) :: temp_basis

        ! Linear momentum in the z direction as density*volume*corrected_velocity_z
        momentum_z = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * var(0,hvar%ids(2))
    end function momentum_z

    function momentum(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: momentum
        type(vector) :: v
        type(basis) :: temp_basis

        ! Magnitude of linear momentum, using corrected velocity
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        momentum = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * &
                    & magnitude(v)
    end function momentum

    function momentum_sphere_r(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: momentum_sphere_r
        type(vector) :: v
        type(basis) :: temp_basis

        ! Linear momentum in the spherical radial direction
        ! 1. Dot product of velocity vector with spherical r
        !    unit vector
        ! 2. Multiply by mass of cell
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        call spherical_basis_from_cartesian(x,temp_basis)
        momentum_sphere_r = (var(0,hvar%ids(1)) * (dx*dx)) * dx * (v .DOT. temp_basis%u(1))
    end function momentum_sphere_r

    function momentum_cyl_z(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: momentum_cyl_z
        type(vector) :: v
        type(basis) :: temp_basis

        ! Linear momentum in the cylindrical z direction
        ! 1. Dot product of velocity vector with cylindrical z
        !    unit vector
        ! 2. Multiply by mass of cell
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        call cylindrical_basis_from_cartesian(x,temp_basis)
        momentum_cyl_z = (var(0,hvar%ids(1)) * (dx*dx)) * dx * (v .DOT. temp_basis%u(3))
    end function momentum_cyl_z

    function ang_momentum_x(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: ang_momentum_x
        type(vector) :: v
        type(basis) :: temp_basis

        ! Corrected angular momentum in the x direction
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        ang_momentum_x = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * (x%y * v%z &
                            &- v%y * x%z)
    end function ang_momentum_x

    function ang_momentum_y(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: ang_momentum_y
        type(vector) :: v
        type(basis) :: temp_basis

        ! Corrected angular momentum in the y direction
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        ang_momentum_y = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * (x%z*v%x &
                            &- v%z*x%x)
    end function ang_momentum_y

    function ang_momentum_z(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: ang_momentum_z
        type(vector) :: v
        type(basis) :: temp_basis

        ! Corrected angular momentum in the z direction
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        ang_momentum_z = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * (x%x*v%y &
                        &- v%x*x%y)
    end function ang_momentum_z

    function ang_momentum(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: ang_momentum
        type(vector) :: v,L
        type(basis) :: temp_basis

        ! Corrected angular momentum in the z direction
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        L = x * v
        ang_momentum = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * magnitude(L)
    end function ang_momentum

    function massflow_rate_sphere_r(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: massflow_rate_sphere_r
        type(vector) :: v
        type(basis) :: temp_basis

        ! Mass flow rate through the cell in the radial direction
            ! Mass per unit time
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        call spherical_basis_from_cartesian(x,temp_basis)
        massflow_rate_sphere_r = (var(0,hvar%ids(1)) * (dx*dx)) * (v .DOT. temp_basis%u(1))    
    end function massflow_rate_sphere_r

    function massflux_rate_sphere_r(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: massflux_rate_sphere_r
        type(vector) :: v
        type(basis) :: temp_basis

        ! Mass flux through the cell in the radial direction
        ! Mass per unit time per unit surface
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        call spherical_basis_from_cartesian(x,temp_basis)
        massflux_rate_sphere_r = var(0,hvar%ids(1)) * (v .DOT. temp_basis%u(1))  
    end function massflux_rate_sphere_r

    function kinetic_energy(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: kinetic_energy

        ! Kinetic energy, computed as 1/2*density*volume*magnitude(velocity)
        ! NOTE: Not corrected!
        kinetic_energy = (0.5d0 * (var(0,hvar%ids(1)) * (dx*dx)) * dx) &
                                 & * sqrt(var(0,hvar%ids(2))**2 + var(0,hvar%ids(3))**2 &
                                 & + var(0,hvar%ids(4))**2)
    end function kinetic_energy

    function sigma(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: sigma
        type(vector) :: v
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar) :: tempvar

        ! Local velocity dispersion
        ! Go back to box coordinates for the central cell, which is transformed usually
        ! before sent to read_amr
        ! Converging flow check
        tempvar = var
        v = (/tempvar(0,hvar%ids(2)),tempvar(0,hvar%ids(3)),tempvar(0,hvar%ids(4))/)
        call rotate_vector(v,transpose(trans_matrix))
        v = v + reg%bulk_velocity
        tempvar(0,hvar%ids(2)) = v%x
        tempvar(0,hvar%ids(3)) = v%y
        tempvar(0,hvar%ids(4)) = v%z
        call cmp_sigma_turb(my_amr,my_sim,hvar,tempvar,sigma)
    end function sigma

    function turb_Mach_number(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var
        real(dbl) :: turb_Mach_number
        real(dbl) :: sound_speed
        real(dbl) :: sigma
        type(vector) :: v
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar) :: tempvar

        ! Local velocity dispersion
        ! Go back to box coordinates for the central cell, which is transformed usually
        ! before sent to read_amr
        ! Converging flow check
        tempvar = var
        v = (/tempvar(0,hvar%ids(2)),tempvar(0,hvar%ids(3)),tempvar(0,hvar%ids(4))/)
        call rotate_vector(v,transpose(trans_matrix))
        v = v + reg%bulk_velocity
        tempvar(0,hvar%ids(2)) = v%x
        tempvar(0,hvar%ids(3)) = v%y
        tempvar(0,hvar%ids(4)) = v%z
        call cmp_sigma_turb(my_amr,my_sim,hvar,tempvar,sigma)

        ! Thermal sound speed, ideal gas
        sound_speed = sqrt(gamma_gas * (max(var(0,hvar%ids(5)), Tmin*var(0,hvar%ids(1))) &
                            & / var(0,hvar%ids(1))))
        turb_Mach_number = sigma / sound_speed
    end function turb_Mach_number

    function temperature(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: temperature

        ! Gas temperature
        temperature = var(0,hvar%ids(2)) / var(0,hvar%ids(1))
        if (temperature < 0d0) then
            temperature = Tmin
        endif
    end function temperature

    function thermal_energy(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: thermal_energy

        ! Thermal energy, computed as thermal_pressure*volume/(gamma - 1)
        thermal_energy = ((max(var(0,hvar%ids(2)), Tmin*var(0,hvar%ids(1))) &
                & / (gamma_gas - 1d0)) * (dx * dx)) * dx
    end function thermal_energy

    function thermal_energy_specific(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: thermal_energy_specific

        ! Specific thermal energy as E_ther/cell mass
        thermal_energy_specific = (max(var(0,hvar%ids(2)), Tmin*var(0,hvar%ids(1))) &
                                    & / (gamma_gas - 1d0)) / var(0,hvar%ids(1))
    end function thermal_energy_specific

    function thermal_energy_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: thermal_energy_density

        ! Thermal energy density  as E_ther/cell volume
        thermal_energy_density = (max(var(0,hvar%ids(2)), Tmin*var(0,hvar%ids(1))) &
                                    & / (gamma_gas - 1d0))
    end function thermal_energy_density

    function entropy_specific(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: entropy_specific

        real(dbl) :: T,rho

        ! Specific entropy, following Gent 2012 equation
        T = var(0,hvar%ids(2))/var(0,hvar%ids(1)) * my_sim%T2 * mu
        rho = var(0,hvar%ids(1)) * my_sim%unit_d / mHydrogen
        entropy_specific = log(T) - (gamma_gas-1d0) * log(rho)
    end function entropy_specific

    function sound_speed(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: sound_speed

        ! Thermal sound speed, ideal gas
        sound_speed = sqrt(gamma_gas * (max(var(0,hvar%ids(2)), Tmin*var(0,hvar%ids(1))) &
                            & / var(0,hvar%ids(1))))
    end function sound_speed

    function grad_thermalpressure(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grad_thermalpressure

        real(dbl) :: dxright,dxleft
        type(vector) :: v

        ! Magnitude of thermal pressure gradient
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (max(var(1,hvar%ids(2)), Tmin*var(1,hvar%ids(1))) - &
                & max(var(2,hvar%ids(2)),Tmin*var(2,hvar%ids(1)))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (max(var(3,hvar%ids(2)), Tmin*var(3,hvar%ids(1))) - &
                & max(var(4,hvar%ids(2)),Tmin*var(4,hvar%ids(1)))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (max(var(5,hvar%ids(2)), Tmin*var(5,hvar%ids(1))) - &
                & max(var(6,hvar%ids(2)),Tmin*var(6,hvar%ids(1)))) / (dxright + dxleft)
        grad_thermalpressure = magnitude(v)
    end function grad_thermalpressure

    function grad_therprsphere(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grad_therprsphere

        type(vector) :: v
        type(basis) :: temp_basis
        real(dbl) :: dxright,dxleft

        ! Thermal pressure gradient in the radial direction
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (max(var(1,hvar%ids(2)), Tmin*var(1,hvar%ids(1))) - &
                & max(var(2,hvar%ids(2)),Tmin*var(2,hvar%ids(1)))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (max(var(3,hvar%ids(2)), Tmin*var(3,hvar%ids(1))) - &
                & max(var(4,hvar%ids(2)),Tmin*var(4,hvar%ids(1)))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (max(var(5,hvar%ids(2)), Tmin*var(5,hvar%ids(1))) - &
                & max(var(6,hvar%ids(2)),Tmin*var(6,hvar%ids(1)))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        call spherical_basis_from_cartesian(x,temp_basis)
        grad_therprsphere = v.DOT.temp_basis%u(1)
    end function grad_therprsphere

    function grad_therpz(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grad_therpz

        type(vector) :: v
        real(dbl) :: dxright,dxleft

        ! Thermal pressure gradient in the z direction
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (max(var(1,hvar%ids(2)), Tmin*var(1,hvar%ids(1))) - &
                & max(var(2,hvar%ids(2)),Tmin*var(2,hvar%ids(1)))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (max(var(3,hvar%ids(2)), Tmin*var(3,hvar%ids(1))) - &
                & max(var(4,hvar%ids(2)),Tmin*var(4,hvar%ids(1)))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (max(var(5,hvar%ids(2)), Tmin*var(5,hvar%ids(1))) - &
                & max(var(6,hvar%ids(2)),Tmin*var(6,hvar%ids(1)))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        grad_therpz = v%z
    end function grad_therpz

    function magnetic_energy(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: magnetic_energy

        type(vector) :: B

        ! Magnetic energy as magnitude(B)**2/2
        B = 0.5 *(/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                    (var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                    (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
        magnetic_energy = (0.5 * (B.DOT.B) * (dx*dx)) * dx
    end function magnetic_energy

    function magnetic_magnitude(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: magnetic_magnitude

        type(vector) :: B

        ! Magnetic field magnitude
        B = 0.5 *(/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                    (var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                    (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
        magnetic_magnitude = magnitude(B)
    end function magnetic_magnitude

    function magnetic_energy_specific(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: magnetic_energy_specific

        type(vector) :: B

        ! Magnetic field magnitude
        B = 0.5 *(/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                    (var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                    (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
        magnetic_energy_specific = (0.5 * (B.DOT.B)) / var(0,hvar%ids(7))
    end function magnetic_energy_specific

    function magnetic_energy_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: magnetic_energy_density

        type(vector) :: B

        ! Magnetic field magnitude
        B = 0.5 *(/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                    (var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                    (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
        magnetic_energy_density = 0.5 * (B.DOT.B)
    end function magnetic_energy_density

    function alfven_speed(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: alfven_speed

        type(vector) :: B

        ! Magnetic field magnitude
        B = 0.5 *(/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                    (var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                    (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
        alfven_speed = magnitude(B) / sqrt(var(0,hvar%ids(7)))
    end function alfven_speed

    function grav_magpfrsphere(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_magpfrsphere

        integer :: i
        type(vector) :: B,v,g
        type(basis) :: temp_basis
        real(dbl) :: dxright,dxleft
        real(dbl),dimension(1:my_amr%twondim) :: totP

        ! Magnetic field magnitude
        totP(:) = 0d0
        do i=1,my_amr%twondim
            B = 0.5 *(/(var(i,hvar%ids(1))+var(i,hvar%ids(4))),&
                        (var(i,hvar%ids(2))+var(i,hvar%ids(5))),&
                        (var(i,hvar%ids(3))+var(i,hvar%ids(6)))/)
            totP(i) = 0.5d0 * (B.DOT.B)
        end do
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (totP(1) - totP(2)) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (totP(3) - totP(4)) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (totP(5) - totP(6)) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)

        g = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_magpfrsphere = -(v.DOT.temp_basis%u(1)) / (var(0,hvar%ids(7)) * (B .DOT. temp_basis%u(1)))
    end function grav_magpfrsphere

    function grav_magpfrspherepos(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_magpfrspherepos

        integer :: i
        type(vector) :: B,v,g
        type(basis) :: temp_basis
        real(dbl) :: dxright,dxleft
        real(dbl),dimension(1:my_amr%twondim) :: totP

        ! Magnetic field magnitude
        totP(:) = 0d0
        do i=1,my_amr%twondim
            B = 0.5 *(/(var(i,hvar%ids(1))+var(i,hvar%ids(4))),&
                        (var(i,hvar%ids(2))+var(i,hvar%ids(5))),&
                        (var(i,hvar%ids(3))+var(i,hvar%ids(6)))/)
            totP(i) = 0.5d0 * (B.DOT.B)
        end do
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (totP(1) - totP(2)) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (totP(3) - totP(4)) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (totP(5) - totP(6)) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)

        g = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_magpfrspherepos = -(v.DOT.temp_basis%u(1)) / (var(0,hvar%ids(7)) * (B .DOT. temp_basis%u(1)))
        if (grav_magpfrspherepos < 0) grav_magpfrspherepos = 0d0
    end function grav_magpfrspherepos

    function grav_magpfrsphereneg(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_magpfrsphereneg

        integer :: i
        type(vector) :: B,v,g
        type(basis) :: temp_basis
        real(dbl) :: dxright,dxleft
        real(dbl),dimension(1:my_amr%twondim) :: totP

        ! Magnetic field magnitude
        totP(:) = 0d0
        do i=1,my_amr%twondim
            B = 0.5 *(/(var(i,hvar%ids(1))+var(i,hvar%ids(4))),&
                        (var(i,hvar%ids(2))+var(i,hvar%ids(5))),&
                        (var(i,hvar%ids(3))+var(i,hvar%ids(6)))/)
            totP(i) = 0.5d0 * (B.DOT.B)
        end do
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (totP(1) - totP(2)) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (totP(3) - totP(4)) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (totP(5) - totP(6)) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)

        g = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_magpfrsphereneg = -(v.DOT.temp_basis%u(1)) / (var(0,hvar%ids(7)) * (B .DOT. temp_basis%u(1)))
        if (grav_magpfrsphereneg > 0) grav_magpfrsphereneg = 0d0
    end function grav_magpfrsphereneg

    function cr_energy(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: cr_energy

        ! CR energy, computed as CR_energydensity*volume
        cr_energy = ((var(0,hvar%ids(1)) / (gamma_cr - 1d0)) * (dx*dx)) * dx
    end function cr_energy

    function cr_energy_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: cr_energy_density

        ! CR energy density
        cr_energy_density = var(0,hvar%ids(1)) / (gamma_cr - 1d0)
    end function cr_energy_density

    function cr_energy_specific(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: cr_energy_specific

        ! Specific CR energy, computed as CR_energydensity*volume/cell mass
        cr_energy_specific = (var(0,hvar%ids(1)) / (gamma_cr - 1d0)) / var(0,hvar%ids(2))
    end function cr_energy_specific

    function cr_temperature_eff(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: cr_temperature_eff

        ! Effective CR temperature
        cr_temperature_eff = var(0,hvar%ids(1)) / var(0,hvar%ids(2))
    end function cr_temperature_eff

    function cr_GH08heat(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: cr_GH08heat

        real(dbl) :: ne,ecr

        ! Cosmic rays hadronic and Coulomb heating from Guo&Ho(2008)
        ! (Assume fully ionised gas)
        ! TODO: Update for RT! 
        ne = var(0,hvar%ids(2)) * my_sim%unit_d / mHydrogen 
        ecr = var(0,hvar%ids(1)) / (gamma_cr - 1d0)
        ecr = ecr * (my_sim%unit_d * ((my_sim%unit_l/my_sim%unit_t)**2))
        cr_GH08heat = lambda_crGH08 * ne * ecr
    end function cr_GH08heat

    function grad_crp(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grad_crp
        type(vector) :: v
        real(dbl) :: dxright,dxleft

        ! Magnitude of CR pressure gradient
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        grad_crp = magnitude(v)
    end function grad_crp

    function grad_crpx(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grad_crpx
        type(vector) :: v
        real(dbl) :: dxright,dxleft

        ! Magnitude of CR pressure gradient
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        grad_crpx = v%x
    end function grad_crpx

    function grad_crpy(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grad_crpy
        type(vector) :: v
        real(dbl) :: dxright,dxleft

        ! Magnitude of CR pressure gradient
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        grad_crpy = v%y
    end function grad_crpy

    function grad_crpz(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grad_crpz
        type(vector) :: v
        real(dbl) :: dxright,dxleft

        ! Magnitude of CR pressure gradient
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        grad_crpz = v%z
    end function grad_crpz

    function grad_crprsphere(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grad_crprsphere
        type(vector) :: v
        real(dbl) :: dxright,dxleft
        type(basis) :: temp_basis

        ! CR pressure gradient in the radial direction
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        call spherical_basis_from_cartesian(x,temp_basis)
        grad_crprsphere = (v.DOT.temp_basis%u(1))
    end function grad_crprsphere

    function gradscale_crprsphere(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: gradscale_crprsphere
        type(vector) :: v
        real(dbl) :: dxright,dxleft
        type(basis) :: temp_basis

        ! CR pressure gradient scale in the radial direction
        ! This is defined as Pcr/grad(Pcr)
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        call spherical_basis_from_cartesian(x,temp_basis)
        gradscale_crprsphere = abs(var(0,hvar%ids(1))/(v.DOT.temp_basis%u(1)))
    end function gradscale_crprsphere

    function gradscale_crp(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: gradscale_crp
        type(vector) :: v
        real(dbl) :: dxright,dxleft
        type(basis) :: temp_basis

        ! CR pressure gradient scale
        ! This is defined as Pcr/grad(Pcr)
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        call spherical_basis_from_cartesian(x,temp_basis)
        gradscale_crp = abs(var(0,hvar%ids(1))/(magnitude(v)))
    end function gradscale_crp

    function diffusion_speed(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: diffusion_speed
        type(vector) :: v
        real(dbl) :: dxright,dxleft
        type(basis) :: temp_basis

        ! CR diffusion speed
        ! This is defined as Dcr/Lcr, with Lcr the CR pressure gradient scale
        ! CR pressure gradient scale
        ! This is defined as Pcr/grad(Pcr)
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        call spherical_basis_from_cartesian(x,temp_basis)
        diffusion_speed = my_sim%Dcr / abs(var(0,hvar%ids(1))/magnitude(v))
    end function diffusion_speed

    function alfvendiff_ratio(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: alfvendiff_ratio
        type(vector) :: v,B
        real(dbl) :: dxright,dxleft
        real(dbl) :: vA,vdiff

        ! Ratio of Alfven to diffusion speed
        B = 0.5d0 * (/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                    &(var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                    (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
        vA = magnitude(B) / sqrt(var(0,hvar%ids(8)))

        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(7)) - var(2,hvar%ids(7))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(7)) - var(4,hvar%ids(7))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(7)) - var(6,hvar%ids(7))) / (dxright + dxleft)
        B = B / magnitude(B)
        vdiff = my_sim%Dcr / abs(var(0,hvar%ids(7))/(abs(v.DOT.B)))
        alfvendiff_ratio = (5d0/3d0) * vA / (vdiff + (5d0/3d0) * vA)
    end function alfvendiff_ratio

    function streaming_heating(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: streaming_heating
        type(vector) :: v,B,vA
        real(dbl) :: dxright,dxleft

        ! Streaming heating
        B = 0.5d0 * (/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                    &(var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                    (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
        vA = B / sqrt(var(0,hvar%ids(8)))

        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(7)) - var(2,hvar%ids(7))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(7)) - var(4,hvar%ids(7))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(7)) - var(6,hvar%ids(7))) / (dxright + dxleft)
        streaming_heating = abs(vA .DOT. v)
    end function streaming_heating

    function net_cooling(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var
    
        real(dbl) :: net_cooling
        real(dbl) :: T,nH,Z,lambda,lambda_prime

        T = var(0,hvar%ids(1)) / var(0,hvar%ids(2)) * my_sim%T2
        ! TODO: This a fix only for some messed up CRMHD simulations!
        if (T<15d0) T = 15d0
        nH = var(0,hvar%ids(2)) * my_sim%nH
        Z  = var(0,hvar%ids(3)) / 2D-2
        call solve_cooling(nH,T,Z,lambda,lambda_prime)
        net_cooling = ((lambda * nH) * nH) * ((my_sim%unit_t**3)/(my_sim%unit_d*(my_sim%unit_l**2)))
    end function net_cooling

    function stheatcooling_ratio(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: stheatcooling_ratio
        type(vector) :: v,B,vA
        real(dbl) :: dxright,dxleft
        real(dbl) :: T,nH,Z,lambda,lambda_prime,lambda_st

        ! Ratio of Alfven to diffusion speed
        B = 0.5d0 * (/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                    &(var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                    (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
        vA = B / sqrt(var(0,hvar%ids(8)))

        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(7)) - var(2,hvar%ids(7))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(7)) - var(4,hvar%ids(7))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(7)) - var(6,hvar%ids(7))) / (dxright + dxleft)
        lambda_st = abs(vA .DOT. v)

        T = var(0,hvar%ids(9)) / var(0,hvar%ids(8)) * my_sim%T2
        ! TODO: This a fix only for some messed up CRMHD simulations!
        if (T<15d0) T = 15d0
        nH = var(0,hvar%ids(8)) * my_sim%nH
        Z  = var(0,hvar%ids(10)) / 2D-2
        call solve_cooling(nH,T,Z,lambda,lambda_prime)
        lambda = ((lambda * nH) * nH) * ((my_sim%unit_t**3)/(my_sim%unit_d*(my_sim%unit_l**2)))
        stheatcooling_ratio = abs(lambda_st / lambda)
    end function stheatcooling_ratio

    function total_coolingtime(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: total_coolingtime
        type(vector) :: v,B,vA
        real(dbl) :: dxright,dxleft
        real(dbl) :: ne,ecr
        real(dbl) :: T,nH,Z,lambda,lambda_prime,lambda_st,lambda_cr

        if (my_sim%cr .and. my_sim%cr_st .and. my_sim%cr_heat) then
            B = 0.5d0 * (/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                        &(var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                        (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
            vA = B / sqrt(var(0,hvar%ids(8)))

            dxright = dx; dxleft = dx
            if (son(1) .eq. 0) dxright = dxright * 1.5D0
            if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
            v%x = (var(1,hvar%ids(7)) - var(2,hvar%ids(7))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(3) .eq. 0) dxright = dxright * 1.5D0
            if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
            v%y = (var(3,hvar%ids(7)) - var(4,hvar%ids(7))) / (dxright + dxleft)
            dxright = dx; dxleft = dx
            if (son(5) .eq. 0) dxright = dxright * 1.5D0
            if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
            v%z = (var(5,hvar%ids(7)) - var(6,hvar%ids(7))) / (dxright + dxleft)
            lambda_st = abs(vA .DOT. v)
        else
            lambda_st = 0D0
        end if

        if (my_sim%cr) then
            ! Cosmic rays hadronic and Coulomb heating from Guo&Ho(2008)
            ! (Assume fully ionised gas)
            ! TODO: Update for RT!
            lambda = 2.63d-16 * ((my_sim%unit_t**3)/(my_sim%unit_d*(my_sim%unit_l**2)))
            ne = var(0,hvar%ids(8)) * my_sim%nH
            ecr = var(0,hvar%ids(7)) / (gamma_cr - 1d0)
            ecr = ecr * (my_sim%unit_d * ((my_sim%unit_l/my_sim%unit_t)**2))
            lambda_cr = lambda * ne * ecr
        else
            lambda_cr = 0D0
        end if

        T = var(0,hvar%ids(9)) / var(0,hvar%ids(8)) * my_sim%T2
        ! TODO: This a fix only for some messed up CRMHD simulations!
        if (T<15d0) T = 15d0
        nH = var(0,hvar%ids(8)) * my_sim%nH
        Z  = var(0,hvar%ids(10)) / 2D-2
        call solve_cooling(nH,T,Z,lambda,lambda_prime)
        lambda = ((lambda * nH) * nH) * ((my_sim%unit_t**3)/(my_sim%unit_d*(my_sim%unit_l**2)))

        ! Thermal energy, computed as thermal_pressure*volume/(gamma - 1)
        ! TODO: This a fix only for some messed up CRMHD simulations!
        if (T<15) then
            total_coolingtime = (Tmin * var(0,hvar%ids(8))) / (gamma_gas - 1d0)
        else
            total_coolingtime = var(0,hvar%ids(9)) / (gamma_gas - 1d0)
        end if

        total_coolingtime = total_coolingtime / (lambda - lambda_st - lambda_cr)
    end function total_coolingtime

    function grav_frsphere(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_frsphere
        type(vector) :: B
        type(basis) :: temp_basis

        B = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_frsphere = var(0,hvar%ids(1)) * (B .DOT. temp_basis%u(1))
    end function grav_frsphere

    function escape_velocity(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: escape_velocity

        escape_velocity = sqrt(2d0*grav_var(0,1))
    end function escape_velocity

    function grav_crpf(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_crpf
        type(vector) :: v,g
        real(dbl) :: dxright,dxleft

        ! Ratio of CR pressure gradient and gravitational acceleration
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        
        g = grav_var(0,2:4)
        grav_crpf = magnitude(v) / (var(0,hvar%ids(2)) * magnitude(g))
    end function grav_crpf

    function grav_crpfz(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_crpfz
        type(vector) :: v,g
        real(dbl) :: dxright,dxleft

        ! Ratio of CR pressure gradient and gravitational acceleration
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        g = grav_var(0,2:4)
        grav_crpfz = -v%z / (var(0,hvar%ids(2)) * magnitude(g))
    end function grav_crpfz

    function grav_therpfz(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_therpfz
        type(vector) :: v,g
        real(dbl) :: dxright,dxleft

        ! Ratio of thermal pressure gradient and gravitational acceleration
        ! in the z direction
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        g = grav_var(0,2:4)
        grav_therpfz = -v%z / (var(0,hvar%ids(2)) * magnitude(g))
    end function grav_therpfz

    function grav_therpfrsphere(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_therpfrsphere
        type(vector) :: v,g
        type(basis) :: temp_basis
        real(dbl) :: dxright,dxleft

        ! Ratio of thermal pressure gradient and gravitational acceleration
        ! in the z direction
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        g = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_therpfrsphere = -(v .DOT. temp_basis%u(1)) / (var(0,hvar%ids(2)) * (g.DOT.temp_basis%u(1)))
    end function grav_therpfrsphere

    function grav_therpfrspherepos(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_therpfrspherepos
        type(vector) :: v,g
        type(basis) :: temp_basis
        real(dbl) :: dxright,dxleft

        ! Ratio of thermal pressure gradient and gravitational acceleration
        ! in the z direction
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        g = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_therpfrspherepos = -(v .DOT. temp_basis%u(1)) / (var(0,hvar%ids(2)) * (g.DOT.temp_basis%u(1)))
        if (grav_therpfrspherepos < 0d0) grav_therpfrspherepos = 0d0
    end function grav_therpfrspherepos

    function grav_therpfrsphereneg(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_therpfrsphereneg
        type(vector) :: v,g
        type(basis) :: temp_basis
        real(dbl) :: dxright,dxleft

        ! Ratio of thermal pressure gradient and gravitational acceleration
        ! in the spherical r direction
        ! ONLY NEGATIVE
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        g = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_therpfrsphereneg = -(v .DOT. temp_basis%u(1)) / (var(0,hvar%ids(2)) * (g.DOT.temp_basis%u(1)))
        if (grav_therpfrsphereneg > 0d0) then
            grav_therpfrsphereneg = 0d0
        else
            grav_therpfrsphereneg = -grav_therpfrsphereneg
        end if
    end function grav_therpfrsphereneg

    function grav_crpfrsphere(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_crpfrsphere
        type(vector) :: v,g
        type(basis) :: temp_basis
        real(dbl) :: dxright,dxleft

        ! Ratio of CR pressure gradient and gravitational acceleration
        ! in the z direction
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        g = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_crpfrsphere = -(v .DOT. temp_basis%u(1)) / (var(0,hvar%ids(2)) * (g.DOT.temp_basis%u(1)))
    end function grav_crpfrsphere

    function grav_crpfrspherepos(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_crpfrspherepos
        type(vector) :: v,g
        type(basis) :: temp_basis
        real(dbl) :: dxright,dxleft

        ! Ratio of CR pressure gradient and gravitational acceleration
        ! in the z direction
        ! ONLY POSITIVE
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        g = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_crpfrspherepos = -(v .DOT. temp_basis%u(1)) / (var(0,hvar%ids(2)) * (g.DOT.temp_basis%u(1)))
        if (grav_crpfrspherepos < 0d0) grav_crpfrspherepos = 0d0
    end function grav_crpfrspherepos

    function grav_crpfrsphereneg(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: grav_crpfrsphereneg
        type(vector) :: v,g
        type(basis) :: temp_basis
        real(dbl) :: dxright,dxleft

        ! Ratio of thermal pressure gradient and gravitational acceleration
        ! in the spherical r direction
        ! ONLY NEGATIVE
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (var(1,hvar%ids(1)) - var(2,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (var(3,hvar%ids(1)) - var(4,hvar%ids(1))) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (var(5,hvar%ids(1)) - var(6,hvar%ids(1))) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        g = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_crpfrsphereneg = -(v .DOT. temp_basis%u(1)) / (var(0,hvar%ids(2)) * (g.DOT.temp_basis%u(1)))
        if (grav_crpfrsphereneg > 0d0) then
            grav_crpfrsphereneg = 0d0
        else
            grav_crpfrsphereneg = -grav_crpfrsphereneg
        end if
    end function grav_crpfrsphereneg

    function grav_totpfrsphere(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        integer :: i
        real(dbl) :: grav_totpfrsphere
        type(vector) :: v,g,B
        type(basis) :: temp_basis
        real(dbl) :: dxright,dxleft
        real(dbl),dimension(1:my_amr%twondim) :: totP

        totP(:) = 0d0
        do i = 1, my_amr%twondim
            if (my_sim%hydro) totP(i) = var(i,hvar%ids(9))
            if (my_sim%mhd) then
                B = 0.5d0 * (/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                        &(var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                        (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
                totP(i) = totP(i) + 0.5d0 * (B.DOT.B)
            end if
            if (my_sim%cr) totP(i) = totP(i) + var(i,hvar%ids(7))
        end do

        ! Ratio of CR pressure gradient and gravitational acceleration
        ! in the z direction
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (totP(1) - totP(2)) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (totP(3) - totP(4)) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (totP(5) - totP(6)) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        g = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_totpfrsphere = -(v .DOT. temp_basis%u(1)) / (var(0,hvar%ids(8)) * (g.DOT.temp_basis%u(1)))
    end function grav_totpfrsphere

    function grav_totpfrspherepos(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        integer :: i
        real(dbl) :: grav_totpfrspherepos
        type(vector) :: v,g,B
        type(basis) :: temp_basis
        real(dbl) :: dxright,dxleft
        real(dbl),dimension(1:my_amr%twondim) :: totP

        totP(:) = 0d0
        do i = 1, my_amr%twondim
            if (my_sim%hydro) totP(i) = var(i,hvar%ids(9))
            if (my_sim%mhd) then
                B = 0.5d0 * (/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                        &(var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                        (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
                totP(i) = totP(i) + 0.5d0 * (B.DOT.B)
            end if
            if (my_sim%cr) totP(i) = totP(i) + var(i,hvar%ids(7))
        end do

        ! Ratio of CR pressure gradient and gravitational acceleration
        ! in the z direction
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (totP(1) - totP(2)) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (totP(3) - totP(4)) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (totP(5) - totP(6)) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        g = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_totpfrspherepos = -(v .DOT. temp_basis%u(1)) / (var(0,hvar%ids(8)) * (g.DOT.temp_basis%u(1)))
        if (grav_totpfrspherepos < 0d0) grav_totpfrspherepos = 0d0
    end function grav_totpfrspherepos

    function grav_totpfrsphereneg(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        integer :: i
        real(dbl) :: grav_totpfrsphereneg
        type(vector) :: v,g,B
        type(basis) :: temp_basis
        real(dbl) :: dxright,dxleft
        real(dbl),dimension(1:my_amr%twondim) :: totP

        totP(:) = 0d0
        do i = 1, my_amr%twondim
            if (my_sim%hydro) totP(i) = var(i,hvar%ids(9))
            if (my_sim%mhd) then
                B = 0.5d0 * (/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                        &(var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                        (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
                totP(i) = totP(i) + 0.5d0 * (B.DOT.B)
            end if
            if (my_sim%cr) totP(i) = totP(i) + var(i,hvar%ids(7))
        end do

        ! Ratio of CR pressure gradient and gravitational acceleration
        ! in the z direction
        dxright = dx; dxleft = dx
        if (son(1) .eq. 0) dxright = dxright * 1.5D0
        if (son(2) .eq. 0) dxleft = dxleft * 1.5D0
        v%x = (totP(1) - totP(2)) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(3) .eq. 0) dxright = dxright * 1.5D0
        if (son(4) .eq. 0) dxleft = dxleft * 1.5D0
        v%y = (totP(3) - totP(4)) / (dxright + dxleft)
        dxright = dx; dxleft = dx
        if (son(5) .eq. 0) dxright = dxright * 1.5D0
        if (son(6) .eq. 0) dxleft = dxleft * 1.5D0
        v%z = (totP(5) - totP(6)) / (dxright + dxleft)
        call rotate_vector(v,trans_matrix)
        g = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_totpfrsphereneg = -(v .DOT. temp_basis%u(1)) / (var(0,hvar%ids(8)) * (g.DOT.temp_basis%u(1)))
        if (grav_totpfrsphereneg > 0d0) then
            grav_totpfrsphereneg = 0d0
        else
            grav_totpfrsphereneg = -grav_totpfrsphereneg
        end if
    end function grav_totpfrsphereneg

    function eff_FKmag(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: eff_FKmag
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar):: tempvar
        type(vector) :: v

        ! Go back to box coordinates for the central cell, which is transformed usually
        ! before sent to read_amr
        tempvar(:,:) = var(:,:)
        v = tempvar(0,hvar%ids(7):hvar%ids(9))
        call rotate_vector(v,transpose(trans_matrix))
        v = v + reg%bulk_velocity
        tempvar(0,hvar%ids(7):hvar%ids(9)) = v
        eff_FKmag = sf_eff(my_amr,my_sim,hvar,reg,dx,x,tempvar,'FKmag',.true.)

    end function eff_FKmag

    function eff_FKmagnocr(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: eff_FKmagnocr
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar):: tempvar
        type(vector) :: v

        ! Go back to box coordinates for the central cell, which is transformed usually
        ! before sent to read_amr
        tempvar(:,:) = var(:,:)
        v = tempvar(0,hvar%ids(7):hvar%ids(9))
        call rotate_vector(v,transpose(trans_matrix))
        v = v + reg%bulk_velocity
        tempvar(0,hvar%ids(7):hvar%ids(9)) = v
        eff_FKmagnocr = sf_eff(my_amr,my_sim,hvar,reg,dx,x,tempvar,'FKmag',.false.)

    end function eff_FKmagnocr

    function eff_FK2(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: eff_FK2
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar):: tempvar
        type(vector) :: v

        ! Go back to box coordinates for the central cell, which is transformed usually
        ! before sent to read_amr
        tempvar(:,:) = var(:,:)
        v = tempvar(0,hvar%ids(7):hvar%ids(9))
        call rotate_vector(v,transpose(trans_matrix))
        v = v + reg%bulk_velocity
        tempvar(0,hvar%ids(7):hvar%ids(9)) = v
        eff_FK2 = sf_eff(my_amr,my_sim,hvar,reg,dx,x,tempvar,'FK2',.false.)

    end function eff_FK2

    ! ELEMENT VARIABLES
    function hydrogen_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: hydrogen_density
        real(dbl) :: metal_mass,dust_mass
        integer :: i

        ! 1. Add up the metal mass
        metal_mass = 0d0
        do i = 2, 10
            metal_mass = metal_mass + var(0,hvar%ids(i))
        end do

        ! 2. Add up the dust mass
        dust_mass = 0d0
        do i = 11, 16
            dust_mass = dust_mass + var(0,hvar%ids(i))
        end do

        hydrogen_density = var(0,hvar%ids(1)) * (1d0 - metal_mass - dust_mass)
    end function hydrogen_density

    function hydrogen_mass(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: hydrogen_mass
        real(dbl) :: metal_mass,dust_mass
        integer :: i

        ! 1. Add up the metal mass
        metal_mass = 0d0
        do i = 2, 10
            metal_mass = metal_mass + var(0,hvar%ids(i))
        end do

        ! 2. Add up the dust mass
        dust_mass = 0d0
        do i = 11, 16
            dust_mass = dust_mass + var(0,hvar%ids(i))
        end do

        hydrogen_mass = (var(0,hvar%ids(1)) * (1d0 - metal_mass - dust_mass)) * dx**3d0
    end function hydrogen_mass

    function HI_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: HI_density
        real(dbl) :: metal_mass,dust_mass,xHI
        integer :: i

        ! 1. Add up the metal mass
        metal_mass = 0d0
        do i = 2, 10
            metal_mass = metal_mass + var(0,hvar%ids(i))
        end do

        ! 2. Add up the dust mass
        dust_mass = 0d0
        do i = 11, 16
            dust_mass = dust_mass + var(0,hvar%ids(i))
        end do

        ! 3. Compute the HI density
        xHI = var(0,hvar%ids(17))
        HI_density = var(0,hvar%ids(1)) * (1d0 - metal_mass - dust_mass) * xHI
    end function HI_density

function HII_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: HII_density
        real(dbl) :: metal_mass,dust_mass,xHII
        integer :: i

        ! 1. Add up the metal mass
        metal_mass = 0d0
        do i = 2, 10
            metal_mass = metal_mass + var(0,hvar%ids(i))
        end do

        ! 2. Add up the dust mass
        dust_mass = 0d0
        do i = 11, 16
            dust_mass = dust_mass + var(0,hvar%ids(i))
        end do

        ! 3. Compute the HI density
        xHII = var(0,hvar%ids(17))
        HII_density = var(0,hvar%ids(1)) * (1d0 - metal_mass - dust_mass) * xHII
    end function HII_density

    function H2_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: H2_density
        real(dbl) :: metal_mass,dust_mass,xH2
        integer :: i

        ! 1. Add up the metal mass
        metal_mass = 0d0
        do i = 2, 10
            metal_mass = metal_mass + var(0,hvar%ids(i))
        end do

        ! 2. Add up the dust mass
        dust_mass = 0d0
        do i = 11, 16
            dust_mass = dust_mass + var(0,hvar%ids(i))
        end do

        ! 3. Compute the xH2
        xH2 = (1.d0 - var(0,hvar%ids(17)) - var(0,hvar%ids(18))) / 2d0

        H2_density = var(0,hvar%ids(1)) * (1d0 - metal_mass - dust_mass) * xH2
    end function H2_density

function H2_mass(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: H2_mass
        real(dbl) :: metal_mass,dust_mass,xH2
        integer :: i

        ! 1. Add up the metal mass
        metal_mass = 0d0
        do i = 2, 10
            metal_mass = metal_mass + var(0,hvar%ids(i))
        end do

        ! 2. Add up the dust mass
        dust_mass = 0d0
        do i = 11, 16
            dust_mass = dust_mass + var(0,hvar%ids(i))
        end do

        ! 3. Compute the xH2
        xH2 = (1.d0 - var(0,hvar%ids(17)) - var(0,hvar%ids(18))) / 2d0

        H2_mass = (var(0,hvar%ids(1)) * (1d0 - metal_mass - dust_mass) * xH2) * (dx**3d0)
    end function H2_mass

    function xH2(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: xH2

        ! 1. Compute the xH2
        xH2 = (1.d0 - var(0,hvar%ids(1)) - var(0,hvar%ids(2))) / 2d0
    end function xH2

    function xCO(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: xCO
        real(dbl) :: mCO,mPAH,mCgrains

        ! 1. Compute the C mass in CO_fraction
        mCO = var(0,hvar%ids(2)) * CO_to_Cmass

        ! 2. Compute the C mass in PAHs
        mPAH = var(0,hvar%ids(3)) + var(0,hvar%ids(4))

        ! 3. Compute the C mass in Cgrains
        mCgrains = var(0,hvar%ids(5)) + var(0,hvar%ids(6))

        ! 4. Compute xCO
        xCO = mCO / (var(0,hvar%ids(1)) + mCO + mPAH + mCgrains)
    end function xCO

    function xCI(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: xCI
        real(dbl) :: mCO,mPAH,mCgrains

        ! 1. Compute the C mass in CO_fraction
        mCO = var(0,hvar%ids(2)) * CO_to_Cmass

        ! 2. Compute the C mass in PAHs
        mPAH = var(0,hvar%ids(3)) + var(0,hvar%ids(4))

        ! 3. Compute the C mass in Cgrains
        mCgrains = var(0,hvar%ids(5)) + var(0,hvar%ids(6))

        ! 4. Compute xCI
        xCI = (var(0,hvar%ids(1))*var(0,hvar%ids(7))) / (var(0,hvar%ids(1)) + mCO + mPAH + mCgrains)
    end function xCI

    function xCII(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: xCII
        real(dbl) :: mCO,mPAH,mCgrains

        ! 1. Compute the C mass in CO_fraction
        mCO = var(0,hvar%ids(2)) * CO_to_Cmass

        ! 2. Compute the C mass in PAHs
        mPAH = var(0,hvar%ids(3)) + var(0,hvar%ids(4))

        ! 3. Compute the C mass in Cgrains
        mCgrains = var(0,hvar%ids(5)) + var(0,hvar%ids(6))

        ! 4. Compute xCII
        xCII = (var(0,hvar%ids(1))*var(0,hvar%ids(7))) / (var(0,hvar%ids(1)) + mCO + mPAH + mCgrains)
    end function xCII

    function oxygen_abundance(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: oxygen_abundance
        real(dbl) :: metal_mass,dust_mass,Hmassfraction,Omassfraction
        integer :: i

        ! 1. Add up the metal mass
        metal_mass = 0d0
        do i = 1, 9
            metal_mass = metal_mass + var(0,hvar%ids(i))
        end do

        ! 2. Add up the dust mass
        dust_mass = 0d0
        do i = 10, 15
            dust_mass = dust_mass + var(0,hvar%ids(i))
        end do

        ! 3. Compute the mass fraction of hydrogen
        Hmassfraction = (1d0 - metal_mass - dust_mass) * XH

        ! 4. Compute the oxygen abundance
        Omassfraction = var(0,hvar%ids(2))
        oxygen_abundance = 12d0 + log10( Omassfraction / Hmassfraction / mOovermH)
    end function oxygen_abundance

    function oxygen_mass(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var
        
        real(dbl) :: oxygen_mass

        oxygen_mass = var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * (dx**3d0)
    end function oxygen_mass


    ! DUST VARIABLES
    function PAHSmall_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: PAHSmall_density

        PAHSmall_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function PAHSmall_density

    function xC_PAHSmall(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: xC_PAHSmall
        real(dbl) :: mCO,mPAH,mCgrains

        ! 1. Compute the C mass in CO_fraction
        mCO = var(0,hvar%ids(2)) * CO_to_Cmass

        ! 2. Compute the C mass in PAHs
        mPAH = var(0,hvar%ids(3)) + var(0,hvar%ids(4))

        ! 3. Compute the C mass in Cgrains
        mCgrains = var(0,hvar%ids(5)) + var(0,hvar%ids(6))

        ! 4. Compute xC_PAHSmall
        xC_PAHSmall = var(0,hvar%ids(3)) / (var(0,hvar%ids(1)) + mCO + mPAH + mCgrains)
        
    end function xC_PAHSmall

    function PAHLarge_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: PAHLarge_density

        PAHLarge_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function PAHLarge_density

    function xC_PAHLarge(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: xC_PAHLarge
        real(dbl) :: mCO,mPAH,mCgrains

        ! 1. Compute the C mass in CO_fraction
        mCO = var(0,hvar%ids(2)) * CO_to_Cmass

        ! 2. Compute the C mass in PAHs
        mPAH = var(0,hvar%ids(3)) + var(0,hvar%ids(4))

        ! 3. Compute the C mass in Cgrains
        mCgrains = var(0,hvar%ids(5)) + var(0,hvar%ids(6))

        ! 4. Compute xC_PAHLarge
        xC_PAHLarge = var(0,hvar%ids(4)) / (var(0,hvar%ids(1)) + mCO + mPAH + mCgrains)
        
    end function xC_PAHLarge

    function CSmall_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: CSmall_density

        CSmall_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function CSmall_density

    function xC_CSmall(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: xC_CSmall
        real(dbl) :: mCO,mPAH,mCgrains

        ! 1. Compute the C mass in CO_fraction
        mCO = var(0,hvar%ids(2)) * CO_to_Cmass

        ! 2. Compute the C mass in PAHs
        mPAH = var(0,hvar%ids(3)) + var(0,hvar%ids(4))

        ! 3. Compute the C mass in Cgrains
        mCgrains = var(0,hvar%ids(5)) + var(0,hvar%ids(6))

        ! 4. Compute xC_CSmall
        xC_CSmall = var(0,hvar%ids(5)) / (var(0,hvar%ids(1)) + mCO + mPAH + mCgrains)
        
    end function xC_CSmall

    function CLarge_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: CLarge_density

        CLarge_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function CLarge_density

    function xC_CLarge(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: xC_CLarge
        real(dbl) :: mCO,mPAH,mCgrains

        ! 1. Compute the C mass in CO_fraction
        mCO = var(0,hvar%ids(2)) * CO_to_Cmass

        ! 2. Compute the C mass in PAHs
        mPAH = var(0,hvar%ids(3)) + var(0,hvar%ids(4))

        ! 3. Compute the C mass in Cgrains
        mCgrains = var(0,hvar%ids(5)) + var(0,hvar%ids(6))

        ! 4. Compute xC_CLarge
        xC_CLarge = var(0,hvar%ids(6)) / (var(0,hvar%ids(1)) + mCO + mPAH + mCgrains)
        
    end function xC_CLarge

    function SilSmall_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: SilSmall_density

        SilSmall_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function SilSmall_density

    function SilLarge_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: SilLarge_density

        SilLarge_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function SilLarge_density

    function CI_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: CI_density

        CI_density = (var(0,hvar%ids(1)) * var(0,hvar%ids(2))) * var(0,hvar%ids(3))
    end function CI_density

    function CII_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: CII_density

        CII_density = (var(0,hvar%ids(1)) * var(0,hvar%ids(2))) * var(0,hvar%ids(3))
    end function CII_density

    function CO_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: CO_density

        CO_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function CO_density

    function PAHSmall_mass(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: PAHSmall_mass

        PAHSmall_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function PAHSmall_mass

    function PAHLarge_mass(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: PAHLarge_mass

        PAHLarge_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function PAHLarge_mass

    function CSmall_mass(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: CSmall_mass

        CSmall_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function CSmall_mass

    function CLarge_mass(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: CLarge_mass

        CLarge_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function CLarge_mass

    function SilSmall_mass(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: SilSmall_mass

        SilSmall_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function SilSmall_mass

    function SilLarge_mass(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: SilLarge_mass

        SilLarge_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function SilLarge_mass

    function CO_mass(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: CO_mass

        CO_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function CO_mass

    function dust_density(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: dust_density
        real(dbl) :: total_dust_fraction

        total_dust_fraction = var(0,hvar%ids(2)) + var(0,hvar%ids(3)) &
                            & + var(0,hvar%ids(4)) + var(0,hvar%ids(5)) &
                            & + var(0,hvar%ids(6)) + var(0,hvar%ids(7))
        dust_density = var(0,hvar%ids(1)) * total_dust_fraction
    end function dust_density

    function dust_mass(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: dust_mass
        real(dbl) :: total_dust_fraction

        total_dust_fraction = var(0,hvar%ids(2)) + var(0,hvar%ids(3)) &
                            & + var(0,hvar%ids(4)) + var(0,hvar%ids(5)) &
                            & + var(0,hvar%ids(6)) + var(0,hvar%ids(7))
        if (total_dust_fraction .ge. 1d0) then
            print *, "Error: Total dust fraction cannot be greater than or equal to 1."
            stop
        end if
        dust_mass = var(0,hvar%ids(1)) * total_dust_fraction * (dx*dx) * dx
    end function dust_mass

    function DTM(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: DTM
        real(dbl) :: total_dust_mass, total_metal_mass
        integer :: i

        total_dust_mass = var(0,hvar%ids(1)) + var(0,hvar%ids(2)) &
                            & + var(0,hvar%ids(3)) + var(0,hvar%ids(4)) &
                            & + var(0,hvar%ids(5)) + var(0,hvar%ids(6))
        if (hvar%ids(7).eq. 0) then
            total_metal_mass = 0d0
            do i = 8, size(hvar%ids(:))
                total_metal_mass = total_metal_mass + var(0,hvar%ids(i))
            end do
        else
            total_metal_mass = var(0,hvar%ids(7))
        end if
        DTM = total_dust_mass / (total_dust_mass + total_metal_mass)
    end function DTM

    function DTG(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: DTG

        DTG = var(0,hvar%ids(1)) + var(0,hvar%ids(2)) + var(0,hvar%ids(3)) + var(0,hvar%ids(4)) &
                & + var(0,hvar%ids(5)) + var(0,hvar%ids(6))
    end function DTG

    function STL(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: STL
        real(dbl) :: total_small_grains, total_large_grains

        total_small_grains = var(0,hvar%ids(1)) + var(0,hvar%ids(3))
        total_large_grains = var(0,hvar%ids(2)) + var(0,hvar%ids(4))
        STL = total_small_grains / total_large_grains
    end function STL

    function CSR(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: CSR
        real(dbl) :: total_carbon, total_silicon

        total_carbon = var(0,hvar%ids(1)) + var(0,hvar%ids(2))
        total_silicon = var(0,hvar%ids(3)) + var(0,hvar%ids(4))
        CSR = total_carbon / total_silicon
    end function CSR

    function qPAH(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: qPAH
        real(dbl) :: total_dust, total_pah

        total_dust = var(0,hvar%ids(1)) + var(0,hvar%ids(2)) + var(0,hvar%ids(3)) + var(0,hvar%ids(4)) + &
                     & var(0,hvar%ids(5)) + var(0,hvar%ids(6))
        total_pah = var(0,hvar%ids(5)) + var(0,hvar%ids(6))
        qPAH = total_pah / total_dust
    end function qPAH

    function fionPAH(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: fionPAH
        real(dbl) :: total_ion, total_pah

        total_pah = var(0,hvar%ids(1)) + var(0,hvar%ids(2))
        total_ion = var(0,hvar%ids(1)) * var(0,hvar%ids(3)) + var(0,hvar%ids(2)) * var(0,hvar%ids(4))
        
        fionPAH = total_ion / total_pah
    end function fionPAH

    function STL_PAH(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: STL_PAH

        STL_PAH = var(0,hvar%ids(1)) / var(0,hvar%ids(2))
    end function STL_PAH

    function CSmall_Stokes(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: CSmall_Stokes
        real(dbl) :: c_s,t_s,t_turb,sdust,adust
        real(dbl) :: vel_dis,value
        type(vector) :: v
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar) :: tempvar

        ! 1. Compute the local sound speed
        c_s = sqrt(gamma_gas * (max(var(0,hvar%ids(2)), Tmin*var(0,hvar%ids(1))) &
                            & / var(0,hvar%ids(1))))

        ! 2. Compute the stopping time
        sdust = 2.2d0 / my_sim%unit_d
        adust = 0.005d-4 / my_sim%unit_l
        t_s = sdust * adust / (var(0,hvar%ids(1)) * sqrt(8d0/pi) * c_s)

        ! 3. Compute the dynamical (turbulent) time
        ! Go back to box coordinates for the central cell, which is transformed usually
        ! before sent to read_amr
        tempvar(:,:) = var(:,:)
        v = (/tempvar(0,hvar%ids(3)),tempvar(0,hvar%ids(4)),tempvar(0,hvar%ids(5))/)
        call rotate_vector(v,transpose(trans_matrix))
        v = v + reg%bulk_velocity
        tempvar(0,hvar%ids(3)) = v%x
        tempvar(0,hvar%ids(4)) = v%y
        tempvar(0,hvar%ids(5)) = v%z
        call cmp_sigma_turb(my_amr,my_sim,hvar,tempvar,vel_dis)
        t_turb = dx / vel_dis

        ! 4. Compute the Stokes number
        CSmall_Stokes = t_s / t_turb
    end function CSmall_Stokes

    function CLarge_Stokes(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: CLarge_Stokes
        real(dbl) :: c_s,t_s,t_turb,sdust,adust
        real(dbl) :: vel_dis,value
        type(vector) :: v
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar) :: tempvar

        ! 1. Compute the local sound speed
        c_s = sqrt(gamma_gas * (max(var(0,hvar%ids(2)), Tmin*var(0,hvar%ids(1))) &
                            & / var(0,hvar%ids(1))))

        ! 2. Compute the stopping time
        sdust = 2.2d0 / my_sim%unit_d
        adust = 0.1d-4 / my_sim%unit_l
        t_s = sdust * adust / (var(0,hvar%ids(1)) * sqrt(8d0/pi) * c_s)

        ! 3. Compute the dynamical (turbulent) time
        ! Go back to box coordinates for the central cell, which is transformed usually
        ! before sent to read_amr
        tempvar(:,:) = var(:,:)
        v = (/tempvar(0,hvar%ids(3)),tempvar(0,hvar%ids(4)),tempvar(0,hvar%ids(5))/)
        call rotate_vector(v,transpose(trans_matrix))
        v = v + reg%bulk_velocity
        tempvar(0,hvar%ids(3)) = v%x
        tempvar(0,hvar%ids(4)) = v%y
        tempvar(0,hvar%ids(5)) = v%z
        call cmp_sigma_turb(my_amr,my_sim,hvar,tempvar,vel_dis)
        t_turb = dx / vel_dis
        
        ! 4. Compute the Stokes number
        CLarge_Stokes = t_s / t_turb
    end function CLarge_Stokes

    function SilSmall_Stokes(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: SilSmall_Stokes
        real(dbl) :: c_s,t_s,t_turb,sdust,adust
        real(dbl) :: vel_dis,value
        type(vector) :: v
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar) :: tempvar

        ! 1. Compute the local sound speed
        c_s = sqrt(gamma_gas * (max(var(0,hvar%ids(2)), Tmin*var(0,hvar%ids(1))) &
                            & / var(0,hvar%ids(1))))

        ! 2. Compute the stopping time
        sdust = 3.3d0 / my_sim%unit_d
        adust = 0.005d-4 / my_sim%unit_l
        t_s = sdust * adust / (var(0,hvar%ids(1)) * sqrt(8d0/pi) * c_s)

        ! 3. Compute the dynamical (turbulent) time
        ! Go back to box coordinates for the central cell, which is transformed usually
        ! before sent to read_amr
        tempvar(:,:) = var(:,:)
        v = (/tempvar(0,hvar%ids(3)),tempvar(0,hvar%ids(4)),tempvar(0,hvar%ids(5))/)
        call rotate_vector(v,transpose(trans_matrix))
        v = v + reg%bulk_velocity
        tempvar(0,hvar%ids(3)) = v%x
        tempvar(0,hvar%ids(4)) = v%y
        tempvar(0,hvar%ids(5)) = v%z
        call cmp_sigma_turb(my_amr,my_sim,hvar,tempvar,vel_dis)
        t_turb = dx / vel_dis
        
        ! 4. Compute the Stokes number
        SilSmall_Stokes = t_s / t_turb
    end function SilSmall_Stokes

    function SilLarge_Stokes(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        real(dbl) :: SilLarge_Stokes
        real(dbl) :: c_s,t_s,t_turb,sdust,adust
        real(dbl) :: vel_dis,value
        type(vector) :: v
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar) :: tempvar

        ! 1. Compute the local sound speed
        c_s = sqrt(gamma_gas * (max(var(0,hvar%ids(2)), Tmin*var(0,hvar%ids(1))) &
                            & / var(0,hvar%ids(1))))

        ! 2. Compute the stopping time
        sdust = 3.3d0 / my_sim%unit_d
        adust = 0.1d-4 / my_sim%unit_l
        t_s = sdust * adust / (var(0,hvar%ids(1)) * sqrt(8d0/pi) * c_s)

        ! 3. Compute the dynamical (turbulent) time
        ! Go back to box coordinates for the central cell, which is transformed usually
        ! before sent to read_amr
        tempvar(:,:) = var(:,:)
        v = (/tempvar(0,hvar%ids(3)),tempvar(0,hvar%ids(4)),tempvar(0,hvar%ids(5))/)
        call rotate_vector(v,transpose(trans_matrix))
        v = v + reg%bulk_velocity
        tempvar(0,hvar%ids(3)) = v%x
        tempvar(0,hvar%ids(4)) = v%y
        tempvar(0,hvar%ids(5)) = v%z
        call cmp_sigma_turb(my_amr,my_sim,hvar,tempvar,vel_dis)
        t_turb = dx / vel_dis
        
        ! 4. Compute the Stokes number
        SilLarge_Stokes = t_s / t_turb
    end function SilLarge_Stokes

    ! RT VARIABLES
    function G0(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var

        integer :: igroup,igrp
        real(dbl) :: G0
        real(dbl) :: energy_density

        ! Compute the energy density in the Habing band (6eV-13.6eV)
        energy_density = 0d0
        do igroup = 1, my_rt%nGroups
            igrp = 1 + (my_amr%ndim + 1) * (igroup - 1)
            if (my_rt%group_egy(igroup) > 5.4 .and. my_rt%group_egy(igroup) < 13.6) then
                energy_density = energy_density + rt_var(0,igrp) * my_rt%group_egy(igroup)
            end if
        end do
        
        ! Convert from code units to cgs units [erg/cm2/s]
        energy_density = energy_density * my_sim%unit_v * my_rt%scale_np * eV2erg

        ! Convert to G0 units (Habing field)
        ! 1 Habing = 1.6d-3 erg/cm2/s
        G0 = energy_density / 1.6d-3
    end function G0

    function charging_gamma(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var
        real(dbl) :: charging_gamma
        real(dbl) :: metal_mass,dust_mass,hydrogen_density
        real(dbl) :: nH,nHe,ne,T,local_mu
        integer :: igroup,igrp,i
        real(dbl) :: G0
        real(dbl) :: energy_density

        ! 1. Add up the metal mass
        metal_mass = 0d0
        do i = 2, 10
            metal_mass = metal_mass + var(0,hvar%ids(i))
        end do

        ! 2. Add up the dust mass
        dust_mass = 0d0
        do i = 11, 16
            dust_mass = dust_mass + var(0,hvar%ids(i))
        end do

        ! 3. Compute nH and nHe
        hydrogen_density = var(0,hvar%ids(1)) * (1d0 - metal_mass - dust_mass) * XH
        nH = hydrogen_density / mHydrogen * my_sim%unit_d
        nHe = 0.25d0 * nH * YHe / XH

        ! 4. Compute ne (TODO: This should also use the metals)
        ne = nH * var(0,hvar%ids(17)) + nHe * (var(0,hvar%ids(18)) + 2d0 * var(0,hvar%ids(19)))

        ! 5. Compute T in K
        local_mu = 1d0 / (XH * (1d0+var(0,hvar%ids(17))) + &
                    & 0.25d0 * YHe * (1d0 + var(0,hvar%ids(18)) + 2d0 * var(0,hvar%ids(19))))
        T = var(0,hvar%ids(20))/var(0,hvar%ids(1)) * my_sim%T2 * local_mu

        ! 6. Compute the energy density in the Habing band (6eV-13.6eV)
        energy_density = 0d0
        do igroup = 1, my_rt%nGroups
            igrp = 1 + (my_amr%ndim + 1) * (igroup - 1)
            if (my_rt%group_egy(igroup) > 5.4 .and. my_rt%group_egy(igroup) < 13.6) then
                energy_density = energy_density + rt_var(0,igrp) * my_rt%group_egy(igroup)
            end if
        end do
        
        ! 7. Convert from code units to cgs units [erg/cm2/s]
        energy_density = energy_density * my_sim%unit_v * my_rt%scale_np * eV2erg

        ! 8. Convert to G0 units (Habing field)
        ! 1 Habing = 1.6d-3 erg/cm2/s
        G0 = energy_density / 1.6d-3

        ! 9. Compute the charging gamma
        charging_gamma = G0 * sqrt(T) / ne

    end function charging_gamma

	function h2form_prism_ratio(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var
        real(dbl) :: h2form_prism_ratio

        real(dbl) :: dtg_mw,h2formation_prism,TK

        ! 1. Add up the dust mass
        dtg_mw = (var(0,hvar%ids(1)) + var(0,hvar%ids(2)) + var(0,hvar%ids(3)) + var(0,hvar%ids(4)) &
        & + var(0,hvar%ids(5)) + var(0,hvar%ids(6))) / (1.d0/162.d0)

        ! 2. Compute the H2 formation rate on dust grains using the PRISM formula
        TK = var(0,hvar%ids(7)) * ((my_sim%unit_v**2d0) * 1.66d-27 / 1.3806200d-23)
        h2formation_prism = 3.5d-17 * dtg_mw * 1d1 * sqrt(TK/100.d0)

        ! 3. Compute the ratio of the live dust H2 formation rate and the one from PRISM
        h2form_prism_ratio = var(0,hvar%ids(8)) / h2formation_prism

    end function h2form_prism_ratio

    function heatingfrac_peh(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var
        real(dbl) :: heatingfrac_peh

        heatingfrac_peh = var(0,hvar%ids(2)) / var(0,hvar%ids(1))

    end function heatingfrac_peh

    function heatingfrac_h2(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var
        real(dbl) :: heatingfrac_h2

        heatingfrac_h2 = var(0,hvar%ids(2)) / var(0,hvar%ids(1))

    end function heatingfrac_h2

    function coolingfrac_dustrecomb(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var
        real(dbl) :: coolingfrac_dustrecomb

        coolingfrac_dustrecomb = var(0,hvar%ids(2)) / var(0,hvar%ids(1))
        
    end function coolingfrac_dustrecomb

    function coolingfrac_dustcoll(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var
        real(dbl) :: coolingfrac_dustcoll

        coolingfrac_dustcoll = var(0,hvar%ids(2)) / var(0,hvar%ids(1))
        
    end function coolingfrac_dustcoll

    function coolingfrac_CO(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var
        real(dbl) :: coolingfrac_CO
        coolingfrac_CO = var(0,hvar%ids(2)) / var(0,hvar%ids(1))
        
    end function coolingfrac_CO

    function coolingfrac_CII(my_amr,my_sim,my_rt,hvar,reg,dx,x,var,son,trans_matrix,grav_var,rt_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(rt_info),intent(in) :: my_rt
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:,:),optional,intent(in) :: grav_var
        real(dbl),dimension(0:,:),optional,intent(in) :: rt_var
        real(dbl) :: coolingfrac_CII
        coolingfrac_CII = var(0,hvar%ids(2)) / var(0,hvar%ids(1))
        
    end function coolingfrac_CII

    subroutine check_dervar(vardict,varname,hvar,ok)
        implicit none

        type(dictf90), intent(in)      :: vardict
        character(128), intent(in) :: varname
        type(hydro_var), intent(inout) :: hvar
        logical, intent(inout)        :: ok
        
        ok = .true.
        select case (trim(varname))
        case ('volume')
            ! Cell volume
            hvar%type = 'derived'
            hvar%name = 'volume'
            hvar%myfunction => cell_volume
        case ('mass')
            ! Cell mass
            hvar%type = 'derived'
            hvar%name = 'mass'
            allocate(hvar%ids(1))
            hvar%ids(1) = vardict%get('density')
            hvar%myfunction => cell_mass
        case ('v_sphere_r')
            ! Velocity component in the spherical radial direction
            ! Dot product of velocity vector with spherical radial
            !    unit vector
            hvar%type = 'derived'
            hvar%name = 'v_sphere_r'
            allocate(hvar%ids(3))
            hvar%ids(1) = vardict%get('velocity_x')
            hvar%ids(2) = vardict%get('velocity_y')
            hvar%ids(3) = vardict%get('velocity_z')
            hvar%myfunction => v_sphere_r
        case ('v_sphere_phi')
            ! Velocity component in the spherical azimutal (phi) direction
            ! Dot product of velocity vector with spherical phi
            !    unit vector
            hvar%type = 'derived'
            hvar%name = 'v_sphere_phi'
            allocate(hvar%ids(3))
            hvar%ids(1) = vardict%get('velocity_x')
            hvar%ids(2) = vardict%get('velocity_y')
            hvar%ids(3) = vardict%get('velocity_z')
            hvar%myfunction => v_sphere_phi
        case ('v_sphere_theta')
            ! Velocity component in the spherical theta direction
            ! Dot product of velocity vector with spherical theta
            !    unit vector
            hvar%type = 'derived'
            hvar%name = 'v_sphere_theta'
            allocate(hvar%ids(3))
            hvar%ids(1) = vardict%get('velocity_x')
            hvar%ids(2) = vardict%get('velocity_y')
            hvar%ids(3) = vardict%get('velocity_z')
            hvar%myfunction => v_sphere_theta
        case ('v_cyl_r')
            ! Velocity component in the cylindrical radial direction
            ! Dot product of velocity vector with cylindrical
            !    radial unit vector
            hvar%type = 'derived'
            hvar%name = 'v_cyl_r'
            allocate(hvar%ids(3))
            hvar%ids(1) = vardict%get('velocity_x')
            hvar%ids(2) = vardict%get('velocity_y')
            hvar%ids(3) = vardict%get('velocity_z')
            hvar%myfunction => v_cyl_r
        case ('v_cyl_z')
            ! Velocity component in the cylindrical z direction
            ! Dot product of velocity vector with cylindrical
            !    z unit vector
            hvar%type = 'derived'
            hvar%name = 'v_cyl_z'
            allocate(hvar%ids(3))
            hvar%ids(1) = vardict%get('velocity_x')
            hvar%ids(2) = vardict%get('velocity_y')
            hvar%ids(3) = vardict%get('velocity_z')
            hvar%myfunction => v_cyl_z
        case ('v_cyl_phi')
            ! Velocity component in the cylyndrical azimutal (phi) direction
            ! Dot product of velocity vector with cylindrical
            !    phi unit vector
            hvar%type = 'derived'
            hvar%name = 'v_cyl_phi'
            allocate(hvar%ids(3))
            hvar%ids(1) = vardict%get('velocity_x')
            hvar%ids(2) = vardict%get('velocity_y')
            hvar%ids(3) = vardict%get('velocity_z')
            hvar%myfunction => v_cyl_phi
        case ('v_magnitude')
            ! Velocity magnitude from galaxy coordinates
            hvar%type = 'derived'
            hvar%name = 'v_magnitude'
            allocate(hvar%ids(3))
            hvar%ids(1) = vardict%get('velocity_x')
            hvar%ids(2) = vardict%get('velocity_y')
            hvar%ids(3) = vardict%get('velocity_z')
            hvar%myfunction => v_magnitude
        case ('v_squared')
            ! Velocity magnitude squared from galaxy coordinates
            hvar%type = 'derived'
            hvar%name = 'v_squared'
            allocate(hvar%ids(3))
            hvar%ids(1) = vardict%get('velocity_x')
            hvar%ids(2) = vardict%get('velocity_y')
            hvar%ids(3) = vardict%get('velocity_z')
            hvar%myfunction => v_squared
        case ('momentum_x')
            ! Linear momentum in the x direction as density*volume*corrected_velocity_x
            hvar%type = 'derived'
            hvar%name = 'momentum_x'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_x')
            hvar%myfunction => momentum_x
        case ('momentum_y')
            ! Linear momentum in the y direction as density*volume*corrected_velocity_y
            hvar%type = 'derived'
            hvar%name = 'momentum_y'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_y')
            hvar%myfunction => momentum_y
        case ('momentum_z')
            ! Linear momentum in the z direction as density*volume*corrected_velocity_z
            hvar%type = 'derived'
            hvar%name = 'momentum_z'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_z')
            hvar%myfunction => momentum_z
        case ('momentum')
            ! Magnitude of linear momentum, using corrected velocity
            hvar%type = 'derived'
            hvar%name = 'momentum'
            allocate(hvar%ids(4))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_x')
            hvar%ids(3) = vardict%get('velocity_y')
            hvar%ids(4) = vardict%get('velocity_z')
            hvar%myfunction => momentum
        case ('momentum_sphere_r')
            ! Linear momentum in the spherical radial direction
            ! 1. Dot product of velocity vector with spherical r
            !    unit vector
            ! 2. Multiply by mass of cell
            hvar%type = 'derived'
            hvar%name = 'momentum_sphere_r'
            allocate(hvar%ids(4))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_x')
            hvar%ids(3) = vardict%get('velocity_y')
            hvar%ids(4) = vardict%get('velocity_z')
            hvar%myfunction => momentum_sphere_r
        case ('momentum_cyl_z')
            ! Linear momentum in the cylindrical z direction
            ! 1. Dot product of velocity vector with cylindrical z
            !    unit vector
            ! 2. Multiply by mass of cell
            hvar%type = 'derived'
            hvar%name = 'momentum_cyl_z'
            allocate(hvar%ids(4))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_x')
            hvar%ids(3) = vardict%get('velocity_y')
            hvar%ids(4) = vardict%get('velocity_z')
            hvar%myfunction => momentum_cyl_z
        case ('ang_momentum_x')
            ! Corrected angular momentum in the x direction
            hvar%type = 'derived'
            hvar%name = 'ang_momentum_x'
            allocate(hvar%ids(4))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_x')
            hvar%ids(3) = vardict%get('velocity_y')
            hvar%ids(4) = vardict%get('velocity_z')
            hvar%myfunction => ang_momentum_x
        case ('ang_momentum_y')
            ! Corrected angular momentum in the y direction
            hvar%type = 'derived'
            hvar%name = 'ang_momentum_y'
            allocate(hvar%ids(4))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_x')
            hvar%ids(3) = vardict%get('velocity_y')
            hvar%ids(4) = vardict%get('velocity_z')
            hvar%myfunction => ang_momentum_y
        case ('ang_momentum_z')
            ! Corrected angular momentum in the z direction
            hvar%type = 'derived'
            hvar%name = 'ang_momentum_z'
            allocate(hvar%ids(4))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_x')
            hvar%ids(3) = vardict%get('velocity_y')
            hvar%ids(4) = vardict%get('velocity_z')
            hvar%myfunction => ang_momentum_z
        case ('ang_momentum')
            ! Corrected magnitude of angular momentum
            hvar%type = 'derived'
            hvar%name = 'ang_momentum'
            allocate(hvar%ids(4))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_x')
            hvar%ids(3) = vardict%get('velocity_y')
            hvar%ids(4) = vardict%get('velocity_z')
            hvar%myfunction => ang_momentum
        case ('massflow_rate_sphere_r')
            ! Mass flow rate through the cell in the radial direction
            ! Mass per unit time
            hvar%type = 'derived'
            hvar%name = 'massflow_rate_sphere_r'
            allocate(hvar%ids(4))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_x')
            hvar%ids(3) = vardict%get('velocity_y')
            hvar%ids(4) = vardict%get('velocity_z')
            hvar%myfunction => massflow_rate_sphere_r
        case ('massflux_rate_sphere_r')
            ! Mass flux through the cell in the radial direction
            ! Mass per unit time per unit surface
            hvar%type = 'derived'
            hvar%name = 'massflux_rate_sphere_r'
            allocate(hvar%ids(4))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_x')
            hvar%ids(3) = vardict%get('velocity_y')
            hvar%ids(4) = vardict%get('velocity_z')
            hvar%myfunction => massflux_rate_sphere_r
        case ('kinetic_energy')
            ! Kinetic energy, computed as 1/2*density*volume*magnitude(velocity)
            ! NOTE: Not corrected!
            hvar%type = 'derived'
            hvar%name = 'kinetic_energy'
            allocate(hvar%ids(4))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_x')
            hvar%ids(3) = vardict%get('velocity_y')
            hvar%ids(4) = vardict%get('velocity_z')
            hvar%myfunction => kinetic_energy
        case ('sigma')
            ! Local velocity dispersion
            hvar%type = 'derived'
            hvar%name = 'sigma'
            allocate(hvar%ids(4))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_x')
            hvar%ids(3) = vardict%get('velocity_y')
            hvar%ids(4) = vardict%get('velocity_z')
            hvar%myfunction => sigma
        case ('temperature')
            ! Gas temperature
            hvar%type = 'derived'
            hvar%name = 'temperature'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('thermal_pressure')
            hvar%myfunction => temperature
        case ('thermal_energy')
            ! Thermal energy, computed as thermal_pressure*volume/(gamma - 1)
            hvar%type = 'derived'
            hvar%name = 'thermal_energy'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('thermal_pressure')
            hvar%myfunction => thermal_energy
        case ('thermal_energy_specific')
            ! Specific thermal energy as E_ther/cell mass
            hvar%type = 'derived'
            hvar%name = 'thermal_energy_specific'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('thermal_pressure')
            hvar%myfunction => thermal_energy_specific
        case ('thermal_energy_density')
            ! Thermal energy density  as E_ther/cell volume
            hvar%type = 'derived'
            hvar%name = 'thermal_energy_density'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('thermal_pressure')
            hvar%myfunction => thermal_energy_density
        case ('entropy_specific')
            ! Specific entropy, following Gent 2012 equation
            hvar%type = 'derived'
            hvar%name = 'entropy_specific'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('thermal_pressure')
            hvar%myfunction => entropy_specific
        case ('sound_speed')
            ! Thermal sound speed, ideal gas
            hvar%type = 'derived'
            hvar%name = 'sound_speed'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('thermal_pressure')
            hvar%myfunction => sound_speed
        case ('turb_Mach_number')
            ! Mach number, using ideal gas sound speed
            hvar%type = 'derived'
            hvar%name = 'turb_Mach_number'
            allocate(hvar%ids(5))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('velocity_x')
            hvar%ids(3) = vardict%get('velocity_y')
            hvar%ids(4) = vardict%get('velocity_z')
            hvar%ids(5) = vardict%get('thermal_pressure')
            hvar%myfunction => turb_Mach_number
        case ('grad_thermalpressure')
            ! Magnitude of thermal pressure gradient
            hvar%type = 'derived'
            hvar%name = 'grad_thermalpressure'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('thermal_pressure')
            hvar%myfunction => grad_thermalpressure
        case ('grad_therprsphere')
            ! Thermal pressure gradient in the radial direction
            hvar%type = 'derived'
            hvar%name = 'grad_therprsphere'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('thermal_pressure')
            hvar%myfunction => grad_therprsphere
        case ('grad_therpz')
            ! Thermal pressure gradient in the z direction
            hvar%type = 'derived'
            hvar%name = 'grad_therpz'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('thermal_pressure')
            hvar%myfunction => grad_therpz
        case ('magnetic_energy')
            ! Magnetic energy as magnitude(B)**2/2
            hvar%type = 'derived'
            hvar%name = 'magnetic_energy'
            allocate(hvar%ids(6))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%myfunction => magnetic_energy
        case ('magnetic_magnitude')
            ! Magnetic field magnitude
            hvar%type = 'derived'
            hvar%name = 'magnetic_magnitude'
            allocate(hvar%ids(6))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%myfunction => magnetic_magnitude
        case ('magnetic_energy_specific')
            ! Specific magnetic energy as E_mag/cell mass
            hvar%type = 'derived'
            hvar%name = 'magnetic_energy_specific'
            allocate(hvar%ids(7))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('density')
            hvar%myfunction => magnetic_energy_specific
        case ('magnetic_energy_density')
            ! Magnetic energy density
            hvar%type = 'derived'
            hvar%name = 'magnetic_energy_density'
            allocate(hvar%ids(6))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%myfunction => magnetic_energy_density
        case ('magnetic_pressure')
            ! Magnetic pressure == magnetic energy density
            hvar%type = 'derived'
            hvar%name = 'magnetic_pressure'
            allocate(hvar%ids(6))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%myfunction => magnetic_energy_density
        case ('alfven_speed')
            ! Alfven speed defined as B / sqrt(rho)
            hvar%type = 'derived'
            hvar%name = 'alfven_speed'
            allocate(hvar%ids(7))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('density')
            hvar%myfunction => alfven_speed
        case ('grav_magpfrsphere')
            ! Ratio of magnetic pressure gradient and gravitational acceleration
            ! in the spherical r direction
            hvar%type = 'derived'
            hvar%name = 'grav_magpfrsphere'
            allocate(hvar%ids(7))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('density')
            hvar%myfunction => grav_magpfrsphere
        case ('grav_magpfrspherepos')
            ! Ratio of magnetic pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY POSITIVE
            hvar%type = 'derived'
            hvar%name = 'grav_magpfrspherepos'
            allocate(hvar%ids(7))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('density')
            hvar%myfunction => grav_magpfrspherepos
        case ('grav_magpfrsphereneg')
            ! Ratio of magnetic pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY NEGATIVE
            hvar%type = 'derived'
            hvar%name = 'grav_magpfrsphereneg'
            allocate(hvar%ids(7))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('density')
            hvar%myfunction => grav_magpfrsphereneg
        case ('cr_energy')
            ! CR energy, computed as CR_energydensity*volume
            hvar%type = 'derived'
            hvar%name = 'cr_energy'
            allocate(hvar%ids(1))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%myfunction => cr_energy
        case ('cr_energy_density')
            ! CR energy density
            hvar%type = 'derived'
            hvar%name = 'cr_energy_density'
            allocate(hvar%ids(1))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%myfunction => cr_energy_density
        case ('cr_energy_specific')
            ! Specific CR energy, computed as CR_energydensity*volume/cell mass
            hvar%type = 'derived'
            hvar%name = 'cr_energy_specific'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%ids(2) = vardict%get('density')
            hvar%myfunction => cr_energy_specific
        case ('cr_temperature_eff')
            ! Effective CR temperature
            hvar%type = 'derived'
            hvar%name = 'cr_temperature_eff'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%ids(2) = vardict%get('density')
            hvar%myfunction => cr_temperature_eff
        case ('cr_GH08heat')
            ! Cosmic rays hadronic and Coulomb heating from Guo&Ho(2008)
            ! (Assume fully ionised gas)
            ! TODO: Update for RT! 
            hvar%type = 'derived'
            hvar%name = 'cr_GH08heat'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%ids(2) = vardict%get('density')
            hvar%myfunction => cr_GH08heat
        case ('grad_crp')
            ! Magnitude of CR pressure gradient
            hvar%type = 'derived'
            hvar%name = 'grad_crp'
            allocate(hvar%ids(1))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%myfunction => grad_crp
        case ('grad_crpx')
            ! Gradient of CR pressure in the x direction
            hvar%type = 'derived'
            hvar%name = 'grad_crpx'
            allocate(hvar%ids(1))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%myfunction => grad_crpx
        case ('grad_crpy')
            ! Gradient of CR pressure in the y direction
            hvar%type = 'derived'
            hvar%name = 'grad_crpy'
            allocate(hvar%ids(1))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%myfunction => grad_crpy
        case ('grad_crpz')
            ! Gradient of CR pressure in the z direction
            hvar%type = 'derived'
            hvar%name = 'grad_crpz'
            allocate(hvar%ids(1))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%myfunction => grad_crpz
        case ('grad_crprsphere')
            ! CR pressure gradient in the radial direction
            hvar%type = 'derived'
            hvar%name = 'grad_crprsphere'
            allocate(hvar%ids(1))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%myfunction => grad_crprsphere
        case ('gradscale_crprsphere')
            ! CR pressure gradient scale in the radial direction
            ! This is defined as Pcr/grad(Pcr)
            hvar%type = 'derived'
            hvar%name = 'gradscale_crprsphere'
            allocate(hvar%ids(1))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%myfunction => gradscale_crprsphere
        case ('gradscale_crp')
            ! CR pressure gradient scale
            ! This is defined as Pcr/grad(Pcr)
            hvar%type = 'derived'
            hvar%name = 'gradscale_crp'
            allocate(hvar%ids(1))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%myfunction => gradscale_crp
        case ('diffusion_speed')
            ! CR diffusion speed
            ! This is defined as Dcr/Lcr, with Lcr the CR pressure gradient scale
            ! CR pressure gradient scale
            ! This is defined as Pcr/grad(Pcr)
            hvar%type = 'derived'
            hvar%name = 'diffusion_speed'
            allocate(hvar%ids(1))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%myfunction => diffusion_speed
        case ('alfvendiff_ratio')
            ! Ratio of Alfven to diffusion speed
            hvar%type = 'derived'
            hvar%name = 'alfvendiff_ratio'
            allocate(hvar%ids(7))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('cr_pressure')
            hvar%ids(8) = vardict%get('density')
            hvar%myfunction => alfvendiff_ratio
        case ('streaming_heating')
            ! CR streaming heating
            hvar%type = 'derived'
            hvar%name = 'streaming_heating'
            allocate(hvar%ids(8))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('cr_pressure')
            hvar%ids(8) = vardict%get('density')
            hvar%myfunction => streaming_heating
        case ('net_cooling')
            ! Net cooling rate taken from the cooling table in RAMSES output
            hvar%type = 'derived'
            hvar%name = 'net_cooling'
            allocate(hvar%ids(3))
            hvar%ids(1) = vardict%get('thermal_pressure')
            hvar%ids(2) = vardict%get('density')
            hvar%ids(3) = vardict%get('metallicity')
            hvar%myfunction => net_cooling
        case ('stheatcooling_ratio')
            ! Ratio of streaming heating rate to gas cooling rate
            hvar%type = 'derived'
            hvar%name = 'stheatcooling_ratio'
            allocate(hvar%ids(10))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('cr_pressure')
            hvar%ids(8) = vardict%get('density')
            hvar%ids(9) = vardict%get('thermal_pressure')
            hvar%ids(10) = vardict%get('metallicity')
            hvar%myfunction => stheatcooling_ratio
        case ('total_coolingtime')
            ! Ratio of streaming heating rate to gas cooling rate
            hvar%type = 'derived'
            hvar%name = 'total_coolingtime'
            allocate(hvar%ids(10))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('cr_pressure')
            hvar%ids(8) = vardict%get('density')
            hvar%ids(9) = vardict%get('thermal_pressure')
            hvar%ids(10) = vardict%get('metallicity')
            hvar%myfunction => total_coolingtime
        case ('grav_frsphere')
            ! Total gravitational force in the radial direction
            hvar%type = 'derived'
            hvar%name = 'grav_frsphere'
            allocate(hvar%ids(1))
            hvar%ids(1) = vardict%get('density')
            hvar%myfunction => grav_frsphere
        case ('escape_velocity')
            ! Local gravitational escape velocity
            hvar%type = 'derived'
            hvar%name = 'escape_velocity'
            hvar%myfunction => escape_velocity
        case ('grav_crpf')
            ! Ratio of CR pressure gradient and gravitational acceleration
            hvar%type = 'derived'
            hvar%name = 'grav_crpf'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%ids(2) = vardict%get('density')
            hvar%myfunction => grav_crpf
        case ('grav_crpfz')
            ! Ratio of CR pressure gradient and gravitational acceleration
            ! in the z direction
            hvar%type = 'derived'
            hvar%name = 'grav_crpfz'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%ids(2) = vardict%get('density')
            hvar%myfunction => grav_crpfz
        case ('grav_therpfz')
            ! Ratio of thermal pressure gradient and gravitational acceleration
            ! in the z direction
            hvar%type = 'derived'
            hvar%name = 'grav_therpfz'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('thermal_pressure')
            hvar%ids(2) = vardict%get('density')
            hvar%myfunction => grav_therpfz
        case ('grav_therpfrsphere')
            ! Ratio of thermal pressure gradient and gravitational acceleration
            ! in the z direction
            hvar%type = 'derived'
            hvar%name = 'grav_therpfrsphere'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('thermal_pressure')
            hvar%ids(2) = vardict%get('density')
            hvar%myfunction => grav_therpfrsphere
        case ('grav_therpfrspherepos')
            ! Ratio of thermal pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY POSITIVE
            hvar%type = 'derived'
            hvar%name = 'grav_therpfrspherepos'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('thermal_pressure')
            hvar%ids(2) = vardict%get('density')
            hvar%myfunction => grav_therpfrspherepos
        case ('grav_therpfrsphereneg')
            ! Ratio of thermal pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY NEGATIVE
            hvar%type = 'derived'
            hvar%name = 'grav_therpfrsphereneg'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('thermal_pressure')
            hvar%ids(2) = vardict%get('density')
            hvar%myfunction => grav_therpfrsphereneg
        case ('grav_crpfrsphere')
            ! Ratio of CR pressure gradient and gravitational acceleration
            ! in the z direction
            hvar%type = 'derived'
            hvar%name = 'grav_crpfrsphere'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%ids(2) = vardict%get('density')
            hvar%myfunction => grav_crpfrsphere
        case ('grav_crpfrspherepos')
            ! Ratio of CR pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY POSITIVE
            hvar%type = 'derived'
            hvar%name = 'grav_crpfrspherepos'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%ids(2) = vardict%get('density')
            hvar%myfunction => grav_crpfrspherepos
        case ('grav_crpfrsphereneg')
            ! Ratio of CR pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY NEGATIVE
            hvar%type = 'derived'
            hvar%name = 'grav_crpfrsphereneg'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('cr_pressure')
            hvar%ids(2) = vardict%get('density')
            hvar%myfunction => grav_crpfrsphereneg
        case ('grav_totpfrsphere')
            ! Ratio of total pressure gradient and gravitational acceleration
            ! in the z direction
            hvar%type = 'derived'
            hvar%name = 'grav_totpfrsphere'
            allocate(hvar%ids(9))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('cr_pressure')
            hvar%ids(8) = vardict%get('density')
            hvar%ids(9) = vardict%get('thermal_pressure')
            hvar%myfunction => grav_totpfrsphere
        case ('grav_totpfrspherepos')
            ! Ratio of total pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY POSITIVE
            hvar%type = 'derived'
            hvar%name = 'grav_totpfrspherepos'
            allocate(hvar%ids(9))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('cr_pressure')
            hvar%ids(8) = vardict%get('density')
            hvar%ids(9) = vardict%get('thermal_pressure')
            hvar%myfunction => grav_totpfrspherepos
        case ('grav_totpfrsphereneg')
            ! Ratio of total pressure gradient and gravitational acceleration
            ! in the spherical r direction
            ! ONLY NEGATIVE
            hvar%type = 'derived'
            hvar%name = 'grav_totpfrsphereneg'
            allocate(hvar%ids(9))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('cr_pressure')
            hvar%ids(8) = vardict%get('density')
            hvar%ids(9) = vardict%get('thermal_pressure')
            hvar%myfunction => grav_totpfrsphereneg
        case ('eff_FKmag')
            ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
            hvar%type = 'derived'
            hvar%name = 'eff_FKmag'
            allocate(hvar%ids(12))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('velocity_x')
            hvar%ids(8) = vardict%get('velocity_y')
            hvar%ids(9) = vardict%get('velocity_z')
            hvar%ids(10) = vardict%get('density')
            hvar%ids(11) = vardict%get('thermal_pressure')
            hvar%ids(12) = vardict%get('cr_pressure')
            hvar%myfunction => eff_FKmag
        case ('eff_FKmagnocr')
            ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
            ! In this case the contribution of CR pressure is ignore in the computation
            ! of the effective sound speed
            hvar%type = 'derived'
            hvar%name = 'eff_FKmagnocr'
            allocate(hvar%ids(12))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('velocity_x')
            hvar%ids(8) = vardict%get('velocity_y')
            hvar%ids(9) = vardict%get('velocity_z')
            hvar%ids(10) = vardict%get('density')
            hvar%ids(11) = vardict%get('thermal_pressure')
            hvar%ids(12) = vardict%get('cr_pressure')
            hvar%myfunction => eff_FKmagnocr
        case ('eff_FK2')
            ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
            hvar%type = 'derived'
            hvar%name = 'eff_FK2'
            allocate(hvar%ids(12))
            hvar%ids(1) = vardict%get('B_left_x')
            hvar%ids(2) = vardict%get('B_left_y')
            hvar%ids(3) = vardict%get('B_left_z')
            hvar%ids(4) = vardict%get('B_right_x')
            hvar%ids(5) = vardict%get('B_right_y')
            hvar%ids(6) = vardict%get('B_right_z')
            hvar%ids(7) = vardict%get('velocity_x')
            hvar%ids(8) = vardict%get('velocity_y')
            hvar%ids(9) = vardict%get('velocity_z')
            hvar%ids(10) = vardict%get('density')
            hvar%ids(11) = vardict%get('thermal_pressure')
            hvar%ids(12) = vardict%get('cr_pressure')
            hvar%myfunction => eff_FK2
        case ('hydrogen_density')
            ! Density of the hydrogen gas
            hvar%type = 'derived'
            hvar%name = 'hydrogen_density'
            allocate(hvar%ids(16))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('iron_fraction')
            hvar%ids(3) = vardict%get('oxygen_fraction')
            hvar%ids(4) = vardict%get('carbon_fraction')
            hvar%ids(5) = vardict%get('silicon_fraction')
            hvar%ids(6) = vardict%get('magnesium_fraction')
            hvar%ids(7) = vardict%get('sulfur_fraction')
            hvar%ids(8) = vardict%get('nitrogen_fraction')
            hvar%ids(9) = vardict%get('calcium_fraction')
            hvar%ids(10) = vardict%get('neon_fraction')
            hvar%ids(11) = vardict%get('PAHSmall_fraction')
            hvar%ids(12) = vardict%get('PAHLarge_fraction')
            hvar%ids(13) = vardict%get('CSmall_fraction')
            hvar%ids(14) = vardict%get('CLarge_fraction')
            hvar%ids(15) = vardict%get('SilSmall_fraction')
            hvar%ids(16) = vardict%get('SilLarge_fraction')
            hvar%myfunction => hydrogen_density
        case ('hydrogen_mass')
            ! Mass of the hydrogen gas
            hvar%type = 'derived'
            hvar%name = 'hydrogen_mass'
            allocate(hvar%ids(16))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('iron_fraction')
            hvar%ids(3) = vardict%get('oxygen_fraction')
            hvar%ids(4) = vardict%get('carbon_fraction')
            hvar%ids(5) = vardict%get('silicon_fraction')
            hvar%ids(6) = vardict%get('magnesium_fraction')
            hvar%ids(7) = vardict%get('sulfur_fraction')
            hvar%ids(8) = vardict%get('nitrogen_fraction')
            hvar%ids(9) = vardict%get('calcium_fraction')
            hvar%ids(10) = vardict%get('neon_fraction')
            hvar%ids(11) = vardict%get('PAHSmall_fraction')
            hvar%ids(12) = vardict%get('PAHLarge_fraction')
            hvar%ids(13) = vardict%get('CSmall_fraction')
            hvar%ids(14) = vardict%get('CLarge_fraction')
            hvar%ids(15) = vardict%get('SilSmall_fraction')
            hvar%ids(16) = vardict%get('SilLarge_fraction')
            hvar%myfunction => hydrogen_mass
        case ('HI_density')
            ! Density of the neutral hydrogen gas
            hvar%type = 'derived'
            hvar%name = 'HI_density'
            allocate(hvar%ids(17))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('iron_fraction')
            hvar%ids(3) = vardict%get('oxygen_fraction')
            hvar%ids(4) = vardict%get('carbon_fraction')
            hvar%ids(5) = vardict%get('silicon_fraction')
            hvar%ids(6) = vardict%get('magnesium_fraction')
            hvar%ids(7) = vardict%get('sulfur_fraction')
            hvar%ids(8) = vardict%get('nitrogen_fraction')
            hvar%ids(9) = vardict%get('calcium_fraction')
            hvar%ids(10) = vardict%get('neon_fraction')
            hvar%ids(11) = vardict%get('PAHSmall_fraction')
            hvar%ids(12) = vardict%get('PAHLarge_fraction')
            hvar%ids(13) = vardict%get('CSmall_fraction')
            hvar%ids(14) = vardict%get('CLarge_fraction')
            hvar%ids(15) = vardict%get('SilSmall_fraction')
            hvar%ids(16) = vardict%get('SilLarge_fraction')
            hvar%ids(17) = vardict%get('xHI')
            hvar%myfunction => HI_density
        case ('HII_density')
            ! Density of the ionised hydrogen gas
            hvar%type = 'derived'
            hvar%name = 'HII_density'
            allocate(hvar%ids(17))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('iron_fraction')
            hvar%ids(3) = vardict%get('oxygen_fraction')
            hvar%ids(4) = vardict%get('carbon_fraction')
            hvar%ids(5) = vardict%get('silicon_fraction')
            hvar%ids(6) = vardict%get('magnesium_fraction')
            hvar%ids(7) = vardict%get('sulfur_fraction')
            hvar%ids(8) = vardict%get('nitrogen_fraction')
            hvar%ids(9) = vardict%get('calcium_fraction')
            hvar%ids(10) = vardict%get('neon_fraction')
            hvar%ids(11) = vardict%get('PAHSmall_fraction')
            hvar%ids(12) = vardict%get('PAHLarge_fraction')
            hvar%ids(13) = vardict%get('CSmall_fraction')
            hvar%ids(14) = vardict%get('CLarge_fraction')
            hvar%ids(15) = vardict%get('SilSmall_fraction')
            hvar%ids(16) = vardict%get('SilLarge_fraction')
            hvar%ids(17) = vardict%get('xHII')
            hvar%myfunction => HII_density
        case ('H2_density')
            ! Density of the molecular hydrogen gas
            hvar%type = 'derived'
            hvar%name = 'H2_density'
            allocate(hvar%ids(18))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('iron_fraction')
            hvar%ids(3) = vardict%get('oxygen_fraction')
            hvar%ids(4) = vardict%get('carbon_fraction')
            hvar%ids(5) = vardict%get('silicon_fraction')
            hvar%ids(6) = vardict%get('magnesium_fraction')
            hvar%ids(7) = vardict%get('sulfur_fraction')
            hvar%ids(8) = vardict%get('nitrogen_fraction')
            hvar%ids(9) = vardict%get('calcium_fraction')
            hvar%ids(10) = vardict%get('neon_fraction')
            hvar%ids(11) = vardict%get('PAHSmall_fraction')
            hvar%ids(12) = vardict%get('PAHLarge_fraction')
            hvar%ids(13) = vardict%get('CSmall_fraction')
            hvar%ids(14) = vardict%get('CLarge_fraction')
            hvar%ids(15) = vardict%get('SilSmall_fraction')
            hvar%ids(16) = vardict%get('SilLarge_fraction')
            hvar%ids(17) = vardict%get('xHI')
            hvar%ids(18) = vardict%get('xHII')
            hvar%myfunction => H2_density
        case ('H2_mass')
            ! Mass of the molecular hydrogen gas
            hvar%type = 'derived'
            hvar%name = 'H2_mass'
            allocate(hvar%ids(18))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('iron_fraction')
            hvar%ids(3) = vardict%get('oxygen_fraction')
            hvar%ids(4) = vardict%get('carbon_fraction')
            hvar%ids(5) = vardict%get('silicon_fraction')
            hvar%ids(6) = vardict%get('magnesium_fraction')
            hvar%ids(7) = vardict%get('sulfur_fraction')
            hvar%ids(8) = vardict%get('nitrogen_fraction')
            hvar%ids(9) = vardict%get('calcium_fraction')
            hvar%ids(10) = vardict%get('neon_fraction')
            hvar%ids(11) = vardict%get('PAHSmall_fraction')
            hvar%ids(12) = vardict%get('PAHLarge_fraction')
            hvar%ids(13) = vardict%get('CSmall_fraction')
            hvar%ids(14) = vardict%get('CLarge_fraction')
            hvar%ids(15) = vardict%get('SilSmall_fraction')
            hvar%ids(16) = vardict%get('SilLarge_fraction')
            hvar%ids(17) = vardict%get('xHI')
            hvar%ids(18) = vardict%get('xHII')
            hvar%myfunction => H2_mass
        case ('xH2')
            ! Fraction of hydrogen in molecular form
            hvar%type = 'derived'
            hvar%name = 'xH2'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('xHI')
            hvar%ids(2) = vardict%get('xHII')
            hvar%myfunction => xH2
        case ('xCO')
            ! Fraction of carbon in CO
            hvar%type = 'derived'
            hvar%name = 'xCO'
            allocate(hvar%ids(6))
            hvar%ids(1) = vardict%get('carbon_fraction')
            hvar%ids(2) = vardict%get('CO_fraction')
            hvar%ids(3) = vardict%get('PAHSmall_fraction')
            hvar%ids(4) = vardict%get('PAHLarge_fraction')
            hvar%ids(5) = vardict%get('CSmall_fraction')
            hvar%ids(6) = vardict%get('CLarge_fraction')
            hvar%myfunction => xCO
        case ('xCI')
            ! Fraction of carbon in CI
            hvar%type = 'derived'
            hvar%name = 'xCI'
            allocate(hvar%ids(7))
            hvar%ids(1) = vardict%get('carbon_fraction')
            hvar%ids(2) = vardict%get('CO_fraction')
            hvar%ids(3) = vardict%get('PAHSmall_fraction')
            hvar%ids(4) = vardict%get('PAHLarge_fraction')
            hvar%ids(5) = vardict%get('CSmall_fraction')
            hvar%ids(6) = vardict%get('CLarge_fraction')
            hvar%ids(7) = vardict%get('carbon_01')
            hvar%myfunction => xCI
        case ('xCII')
            ! Fraction of carbon in CII
            hvar%type = 'derived'
            hvar%name = 'xCII'
            allocate(hvar%ids(7))
            hvar%ids(1) = vardict%get('carbon_fraction')
            hvar%ids(2) = vardict%get('CO_fraction')
            hvar%ids(3) = vardict%get('PAHSmall_fraction')
            hvar%ids(4) = vardict%get('PAHLarge_fraction')
            hvar%ids(5) = vardict%get('CSmall_fraction')
            hvar%ids(6) = vardict%get('CLarge_fraction')
            hvar%ids(7) = vardict%get('carbon_02')
            hvar%myfunction => xCII
        case ('CI_density')
            ! Density of the neutral carbon
            hvar%type = 'derived'
            hvar%name = 'CI_density'
            allocate(hvar%ids(3))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('carbon_fraction')
            hvar%ids(3) = vardict%get('carbon_01')
            hvar%myfunction => CI_density
        case ('CII_density')
            ! Density of the ionised carbon
            hvar%type = 'derived'
            hvar%name = 'CII_density'
            allocate(hvar%ids(3))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('carbon_fraction')
            hvar%ids(3) = vardict%get('carbon_02')
            hvar%myfunction => CII_density
        case ('PAHSmall_density')
            ! Density of the small PAHs
            hvar%type = 'derived'
            hvar%name = 'PAHSmall_density'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('PAHSmall_fraction')
            hvar%myfunction => PAHSmall_density
        case ('xC_PAHSmall')
            ! Fraction of carbon in PAHSmall
            hvar%type = 'derived'
            hvar%name = 'xC_PAHSmall'
            allocate(hvar%ids(6))
            hvar%ids(1) = vardict%get('carbon_fraction')
            hvar%ids(2) = vardict%get('CO_fraction')
            hvar%ids(3) = vardict%get('PAHSmall_fraction')
            hvar%ids(4) = vardict%get('PAHLarge_fraction')
            hvar%ids(5) = vardict%get('CSmall_fraction')
            hvar%ids(6) = vardict%get('CLarge_fraction')
            hvar%myfunction => xC_PAHSmall
        case ('PAHLarge_density')
            ! Density of the large PAHs
            hvar%type = 'derived'
            hvar%name = 'PAHLarge_density'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('PAHLarge_fraction')
            hvar%myfunction => PAHLarge_density
        case ('xC_PAHLarge')
            ! Fraction of carbon in PAHLarge
            hvar%type = 'derived'
            hvar%name = 'xC_PAHLarge'
            allocate(hvar%ids(6))
            hvar%ids(1) = vardict%get('carbon_fraction')
            hvar%ids(2) = vardict%get('CO_fraction')
            hvar%ids(3) = vardict%get('PAHSmall_fraction')
            hvar%ids(4) = vardict%get('PAHLarge_fraction')
            hvar%ids(5) = vardict%get('CSmall_fraction')
            hvar%ids(6) = vardict%get('CLarge_fraction')
            hvar%myfunction => xC_PAHLarge
        case ('CSmall_density')
            ! Density of the small C grains
            hvar%type = 'derived'
            hvar%name = 'CSmall_density'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('CSmall_fraction')
            hvar%myfunction => CSmall_density
        case ('xC_CSmall')
            ! Fraction of carbon in CSmall
            hvar%type = 'derived'
            hvar%name = 'xC_CSmall'
            allocate(hvar%ids(6))
            hvar%ids(1) = vardict%get('carbon_fraction')
            hvar%ids(2) = vardict%get('CO_fraction')
            hvar%ids(3) = vardict%get('PAHSmall_fraction')
            hvar%ids(4) = vardict%get('PAHLarge_fraction')
            hvar%ids(5) = vardict%get('CSmall_fraction')
            hvar%ids(6) = vardict%get('CLarge_fraction')
            hvar%myfunction => xC_CSmall
        case ('CLarge_density')
            ! Density of the large PAHs
            hvar%type = 'derived'
            hvar%name = 'CLarge_density'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('CLarge_fraction')
            hvar%myfunction => CLarge_density
        case ('xC_CLarge')
            ! Fraction of carbon in CLarge
            hvar%type = 'derived'
            hvar%name = 'xC_CLarge'
            allocate(hvar%ids(6))
            hvar%ids(1) = vardict%get('carbon_fraction')
            hvar%ids(2) = vardict%get('CO_fraction')
            hvar%ids(3) = vardict%get('PAHSmall_fraction')
            hvar%ids(4) = vardict%get('PAHLarge_fraction')
            hvar%ids(5) = vardict%get('CSmall_fraction')
            hvar%ids(6) = vardict%get('CLarge_fraction')
            hvar%myfunction => xC_CLarge
        case ('SilSmall_density')
            ! Density of the small PAHs
            hvar%type = 'derived'
            hvar%name = 'SilSmall_density'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('SilSmall_fraction')
            hvar%myfunction => SilSmall_density
        case ('SilLarge_density')
            ! Density of the large PAHs
            hvar%type = 'derived'
            hvar%name = 'SilLarge_density'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('SilLarge_fraction')
            hvar%myfunction => SilLarge_density
        case ('CO_density')
            ! Density of the CO molecules
            hvar%type = 'derived'
            hvar%name = 'CO_density'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('CO_fraction')
            hvar%myfunction => CO_density
        case ('PAHSmall_mass')
            ! mass of the small PAHs
            hvar%type = 'derived'
            hvar%name = 'PAHSmall_mass'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('PAHSmall_fraction')
            hvar%myfunction => PAHSmall_mass
        case ('PAHLarge_mass')
            ! mass of the large PAHs
            hvar%type = 'derived'
            hvar%name = 'PAHLarge_mass'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('PAHLarge_fraction')
            hvar%myfunction => PAHLarge_mass
        case ('CSmall_mass')
            ! mass of the small C grains
            hvar%type = 'derived'
            hvar%name = 'CSmall_mass'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('CSmall_fraction')
            hvar%myfunction => CSmall_mass
        case ('CLarge_mass')
            ! mass of the large PAHs
            hvar%type = 'derived'
            hvar%name = 'CLarge_mass'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('CLarge_fraction')
            hvar%myfunction => CLarge_mass
        case ('SilSmall_mass')
            ! mass of the small PAHs
            hvar%type = 'derived'
            hvar%name = 'SilSmall_mass'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('SilSmall_fraction')
            hvar%myfunction => SilSmall_mass
        case ('SilLarge_mass')
            ! mass of the large PAHs
            hvar%type = 'derived'
            hvar%name = 'SilLarge_mass'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('SilLarge_fraction')
            hvar%myfunction => SilLarge_mass
        case ('CO_mass')
            ! mass of the CO molecules
            hvar%type = 'derived'
            hvar%name = 'CO_mass'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('CO_fraction')
            hvar%myfunction => CO_mass
        case ('oxygen_mass')
            ! mass of the oxygen
            hvar%type = 'derived'
            hvar%name = 'oxygen_mass'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('oxygen_fraction')
            hvar%myfunction => oxygen_mass
        case ('oxygen_abundance')
            ! Oxygen abundance 12+log(O/H)
            hvar%type = 'derived'
            hvar%name = 'oxygen_abundance'
            allocate(hvar%ids(15))
            hvar%ids(1) = vardict%get('iron_fraction')
            hvar%ids(2) = vardict%get('oxygen_fraction')
            hvar%ids(3) = vardict%get('carbon_fraction')
            hvar%ids(4) = vardict%get('silicon_fraction')
            hvar%ids(5) = vardict%get('magnesium_fraction')
            hvar%ids(6) = vardict%get('sulfur_fraction')
            hvar%ids(7) = vardict%get('nitrogen_fraction')
            hvar%ids(8) = vardict%get('calcium_fraction')
            hvar%ids(9) = vardict%get('neon_fraction')
            hvar%ids(10) = vardict%get('PAHSmall_fraction')
            hvar%ids(11) = vardict%get('PAHLarge_fraction')
            hvar%ids(12) = vardict%get('CSmall_fraction')
            hvar%ids(13) = vardict%get('CLarge_fraction')
            hvar%ids(14) = vardict%get('SilSmall_fraction')
            hvar%ids(15) = vardict%get('SilLarge_fraction')
            hvar%myfunction => oxygen_abundance
        case ('dust_density')
            ! Density of the dust
            hvar%type = 'derived'
            hvar%name = 'dust_density'
            allocate(hvar%ids(7))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('CSmall_fraction')
            hvar%ids(3) = vardict%get('CLarge_fraction')
            hvar%ids(4) = vardict%get('SilSmall_fraction')
            hvar%ids(5) = vardict%get('SilLarge_fraction')
            hvar%ids(6) = vardict%get('PAHSmall_fraction')
            hvar%ids(7) = vardict%get('PAHLarge_fraction')
            hvar%myfunction => dust_density
        case ('dust_mass')
            ! Mass of the dust
            hvar%type = 'derived'
            hvar%name = 'dust_mass'
            allocate(hvar%ids(7))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('CSmall_fraction')
            hvar%ids(3) = vardict%get('CLarge_fraction')
            hvar%ids(4) = vardict%get('SilSmall_fraction')
            hvar%ids(5) = vardict%get('SilLarge_fraction')
            hvar%ids(6) = vardict%get('PAHSmall_fraction')
            hvar%ids(7) = vardict%get('PAHLarge_fraction')
            hvar%myfunction => dust_mass
        case ('DTM')
            ! dust-to-metal mass ratio
            hvar%type = 'derived'
            hvar%name = 'DTM'
            allocate(hvar%ids(17))
            hvar%ids(1) = vardict%get('CSmall_fraction')
            hvar%ids(2) = vardict%get('CLarge_fraction')
            hvar%ids(3) = vardict%get('SilSmall_fraction')
            hvar%ids(4) = vardict%get('SilLarge_fraction')
            hvar%ids(5) = vardict%get('PAHSmall_fraction')
            hvar%ids(6) = vardict%get('PAHLarge_fraction')
            hvar%ids(7) = vardict%get('metallicity')
            hvar%ids(8) = vardict%get('iron_fraction')
            hvar%ids(9) = vardict%get('oxygen_fraction')
            hvar%ids(10) = vardict%get('nitrogen_fraction')
            hvar%ids(11) = vardict%get('magnesium_fraction')
            hvar%ids(12) = vardict%get('neon_fraction')
            hvar%ids(13) = vardict%get('silicon_fraction')
            hvar%ids(14) = vardict%get('calcium_fraction')
            hvar%ids(15) = vardict%get('carbon_fraction')
            hvar%ids(16) = vardict%get('sulfur_fraction')
            hvar%ids(17) = vardict%get('CO_fraction')
            hvar%myfunction => DTM
        case ('DTG')
            ! dust-to-gas mass ratio
            hvar%type = 'derived'
            hvar%name = 'DTG'
            allocate(hvar%ids(6))
            hvar%ids(1) = vardict%get('CSmall_fraction')
            hvar%ids(2) = vardict%get('CLarge_fraction')
            hvar%ids(3) = vardict%get('SilSmall_fraction')
            hvar%ids(4) = vardict%get('SilLarge_fraction')
            hvar%ids(5) = vardict%get('PAHSmall_fraction') 
            hvar%ids(6) = vardict%get('PAHLarge_fraction')
            hvar%myfunction => DTG
        case ('STL')
            ! small-to-large grain mass ratio
            hvar%type = 'derived'
            hvar%name = 'STL'
            allocate(hvar%ids(5))
            hvar%ids(1) = vardict%get('CSmall_fraction')
            hvar%ids(2) = vardict%get('CLarge_fraction')
            hvar%ids(3) = vardict%get('SilSmall_fraction')
            hvar%ids(4) = vardict%get('SilLarge_fraction')
            hvar%ids(5) = vardict%get('metallicity')
            hvar%myfunction => STL
        case ('CSR')
            ! carbon-to-silicate grain mass ratio
            hvar%type = 'derived'
            hvar%name = 'CSR'
            allocate(hvar%ids(5))
            hvar%ids(1) = vardict%get('CSmall_fraction')
            hvar%ids(2) = vardict%get('CLarge_fraction')
            hvar%ids(3) = vardict%get('SilSmall_fraction')
            hvar%ids(4) = vardict%get('SilLarge_fraction')
            hvar%ids(5) = vardict%get('metallicity')
            hvar%myfunction => CSR
        case ('qPAH')
            ! PAH-to-dust mass ratio
            hvar%type = 'derived'
            hvar%name = 'qPAH'
            allocate(hvar%ids(6))
            hvar%ids(1) = vardict%get('CSmall_fraction')
            hvar%ids(2) = vardict%get('CLarge_fraction')
            hvar%ids(3) = vardict%get('SilSmall_fraction')
            hvar%ids(4) = vardict%get('SilLarge_fraction')
            hvar%ids(5) = vardict%get('PAHSmall_fraction')
            hvar%ids(6) = vardict%get('PAHLarge_fraction')
            hvar%myfunction => qPAH
        case ('fionPAH')
            ! PAH-to-dust mass ratio
            hvar%type = 'derived'
            hvar%name = 'fionPAH'
            allocate(hvar%ids(4))
            hvar%ids(1) = vardict%get('PAHSmall_fraction')
            hvar%ids(2) = vardict%get('PAHLarge_fraction')
            hvar%ids(3) = vardict%get('PAHSmall_fion')
            hvar%ids(4) = vardict%get('PAHLarge_fion')
            hvar%myfunction => fionPAH
        case ('STL_PAH')
            ! PAH-to-dust mass ratio
            hvar%type = 'derived'
            hvar%name = 'STL_PAH'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('PAHSmall_fraction')
            hvar%ids(2) = vardict%get('PAHLarge_fraction')
            hvar%myfunction => STL_PAH
        case ('CSmall_Stokes')
            ! Stokes number of the small C grains
            hvar%type = 'derived'
            hvar%name = 'CSmall_Stokes'
            allocate(hvar%ids(5))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('thermal_pressure')
            hvar%ids(3) = vardict%get('velocity_x')
            hvar%ids(4) = vardict%get('velocity_y')
            hvar%ids(5) = vardict%get('velocity_z')
            hvar%myfunction => CSmall_Stokes
        case ('CLarge_Stokes')
            ! Stokes number of the large C grains
            hvar%type = 'derived'
            hvar%name = 'CLarge_Stokes'
            allocate(hvar%ids(5))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('thermal_pressure')
            hvar%ids(3) = vardict%get('velocity_x')
            hvar%ids(4) = vardict%get('velocity_y')
            hvar%ids(5) = vardict%get('velocity_z')
            hvar%myfunction => CLarge_Stokes
        case ('SilSmall_Stokes')
            ! Stokes number of the small silicate grains
            hvar%type = 'derived'
            hvar%name = 'SilSmall_Stokes'
            allocate(hvar%ids(5))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('thermal_pressure')
            hvar%ids(3) = vardict%get('velocity_x')
            hvar%ids(4) = vardict%get('velocity_y')
            hvar%ids(5) = vardict%get('velocity_z')
            hvar%myfunction => SilSmall_Stokes
        case ('SilLarge_Stokes')
            ! Stokes number of the large silicate grains
            hvar%type = 'derived'
            hvar%name = 'SilLarge_Stokes'
            allocate(hvar%ids(5))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('thermal_pressure')
            hvar%ids(3) = vardict%get('velocity_x')
            hvar%ids(4) = vardict%get('velocity_y')
            hvar%ids(5) = vardict%get('velocity_z')
            hvar%myfunction => SilLarge_Stokes
        case ('G0')
            ! Interstellar radiation field strength in Habing units
            hvar%type = 'derived'
            hvar%name = 'G0'
            hvar%myfunction => G0
        case ('charging_gamma')
            ! Charging gamma parameter G0 sqrt(T)/n_e
            hvar%type = 'derived'
            hvar%name = 'charging_gamma'
            allocate(hvar%ids(20))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('iron_fraction')
            hvar%ids(3) = vardict%get('oxygen_fraction')
            hvar%ids(4) = vardict%get('carbon_fraction')
            hvar%ids(5) = vardict%get('silicon_fraction')
            hvar%ids(6) = vardict%get('magnesium_fraction')
            hvar%ids(7) = vardict%get('sulfur_fraction')
            hvar%ids(8) = vardict%get('nitrogen_fraction')
            hvar%ids(9) = vardict%get('calcium_fraction')
            hvar%ids(10) = vardict%get('neon_fraction')
            hvar%ids(11) = vardict%get('PAHSmall_fraction')
            hvar%ids(12) = vardict%get('PAHLarge_fraction')
            hvar%ids(13) = vardict%get('CSmall_fraction')
            hvar%ids(14) = vardict%get('CLarge_fraction')
            hvar%ids(15) = vardict%get('SilSmall_fraction')
            hvar%ids(16) = vardict%get('SilLarge_fraction')
            hvar%ids(17) = vardict%get('xHII')
            hvar%ids(18) = vardict%get('xHeII')
            hvar%ids(19) = vardict%get('xHeIII')
            hvar%ids(20) = vardict%get('thermal_pressure')
            hvar%myfunction => charging_gamma
        case ('h2form_prism_ratio')
            ! Ratio of H2 formation from live dust to the scaling
            ! in the original PRISM model
            hvar%type = 'derived'
            hvar%name = 'h2form_prism_ratio'
            allocate(hvar%ids(8))
            hvar%ids(1) = vardict%get('CSmall_fraction')
            hvar%ids(2) = vardict%get('CLarge_fraction')
            hvar%ids(3) = vardict%get('SilSmall_fraction')
            hvar%ids(4) = vardict%get('SilLarge_fraction')
            hvar%ids(5) = vardict%get('PAHSmall_fraction')
            hvar%ids(6) = vardict%get('PAHLarge_fraction')
            hvar%ids(7) = vardict%get('cooling_temperature')
            hvar%ids(8) = vardict%get('h2formation_rate')
            hvar%myfunction => h2form_prism_ratio
        case ('heatingfrac_peh')
            ! Fraction of total heating contributed by PEH
            hvar%type = 'derived'
            hvar%name = 'heatingfrac_peh'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('heating_rate')
            hvar%ids(2) = vardict%get('heating_pe')
            hvar%myfunction => heatingfrac_peh
        case ('heatingfrac_h2')
            ! Fraction of total heating contributed by H2
            hvar%type = 'derived'
            hvar%name = 'heatingfrac_h2'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('heating_rate')
            hvar%ids(2) = vardict%get('heating_h2')
            hvar%myfunction => heatingfrac_h2
        case ('coolingfrac_dustrecomb')
            ! Fraction of total cooling contributed by dust recombination
            hvar%type = 'derived'
            hvar%name = 'coolingfrac_dustrecomb'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('cooling_rate')
            hvar%ids(2) = vardict%get('cooling_dust_rec')
            hvar%myfunction => coolingfrac_dustrecomb
        case ('coolingfrac_dustcoll')
            ! Fraction of total cooling contributed by dust collisional cooling
            hvar%type = 'derived'
            hvar%name = 'coolingfrac_dustcoll'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('cooling_rate')
            hvar%ids(2) = vardict%get('cooling_dust')
            hvar%myfunction => coolingfrac_dustcoll
        case ('coolingfrac_CO')
            ! Fraction of total cooling contributed by CO cooling
            hvar%type = 'derived'
            hvar%name = 'coolingfrac_CO'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('cooling_rate')
            hvar%ids(2) = vardict%get('cooling_CO')
            hvar%myfunction => coolingfrac_CO
        case ('coolingfrac_CII')
            ! Fraction of total cooling contributed by CII cooling
            hvar%type = 'derived'
            hvar%name = 'coolingfrac_CII'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('cooling_rate')
            hvar%ids(2) = vardict%get('cooling_CII')
            hvar%myfunction => coolingfrac_CII
        case default
            ok = .false.
        end select

    end subroutine check_dervar

    subroutine get_var_tools(vardict,nreq,reqvars,cleaned_vars)
        implicit none

        type(dictf90), intent(in)      :: vardict
        integer, intent(in)            :: nreq
        character(128), dimension(:), intent(in):: reqvars
        type(hydro_var), dimension(:),intent(inout) :: cleaned_vars

        logical                        :: ok_check
        integer                        :: i, ivar
        
        ! Loop over requested variable names
        do i = 1, nreq
            ivar = vardict%get(trim(reqvars(i)))
            if (ivar.ne.0) then
                ! 1. If variable is a raw variable defined in the
                ! hydro descriptor dictionary, just simply store it
                cleaned_vars(i)%type = 'raw'
                cleaned_vars(i)%name = reqvars(i)
                allocate(cleaned_vars(i)%ids(1))
                cleaned_vars(i)%ids(1) = ivar
                cleaned_vars(i)%myfunction => raw_hydro
            else if (trim(reqvars(i)) == 'cumulative') then
                ! Cumulative variable
                cleaned_vars(i)%type = 'cumulative'
                cleaned_vars(i)%name = reqvars(i)
                allocate(cleaned_vars(i)%ids(1))
                cleaned_vars(i)%ids(1) = 0
            else if (trim(reqvars(i)) == 'column') then
                ! Column density weighting (special weight type)
                cleaned_vars(i)%type = 'column'
                cleaned_vars(i)%name = reqvars(i)
                allocate(cleaned_vars(i)%ids(1))
                cleaned_vars(i)%ids(1) = 0
            else if (trim(reqvars(i)) == 'counts') then
                ! Counts weighting (special weight type)
                cleaned_vars(i)%type = 'counts'
                cleaned_vars(i)%name = reqvars(i)
                allocate(cleaned_vars(i)%ids(1))
                cleaned_vars(i)%ids(1) = 0
            else
                ! 2. Variable is not a raw variable, so it is either a
                ! geometrical variable a derived variable or a gravity
                ! variable
                call check_geovar(reqvars(i),cleaned_vars(i),ok_check)
                if (ok_check) cycle

                ! 3. Then the variable may be a derived one, so we go
                ! through the defined derived variables to pre-define
                ! the indexes (for speed-up)
                call check_dervar(vardict,reqvars(i),cleaned_vars(i),ok_check)
                if (ok_check) cycle

                ! 4. The option left: we are asked for a gravity variable
                call check_gravvar(reqvars(i),cleaned_vars(i),ok_check)
                
                if (.not.ok_check) then
                    write(*,*)'Variable not supported: ',TRIM(reqvars(i))
                    write(*,*)'Aborting!',nreq,i
                    stop
                end if
            end if
        end do
    end subroutine get_var_tools

    subroutine set_hydro_var(vardict,hvar)
        implicit none

        type(dictf90), intent(in)      :: vardict
        type(hydro_var), intent(inout) :: hvar

        logical                        :: ok_check
        integer                        :: i, ivar
        
        ivar = vardict%get(hvar%name)
        
        if (ivar.ne.0) then
            ! 1. If variable is a raw variable defined in the
            ! hydro descriptor dictionary, just simply store it
            hvar%type = 'raw'
            allocate(hvar%ids(1))
            hvar%ids(1) = ivar
            hvar%myfunction => raw_hydro
        else if (trim(hvar%name) == 'cumulative') then
            ! Cumulative variable
            hvar%type = 'cumulative'
            allocate(hvar%ids(1))
            hvar%ids(1) = 0
        else if (trim(hvar%name) == 'column') then
            ! Column density weighting (special weight type)
            hvar%type = 'column'
            allocate(hvar%ids(1))
            hvar%ids(1) = 0
        else if (trim(hvar%name) == 'counts') then
            ! Counts weighting (special weight type)
            hvar%type = 'counts'
            allocate(hvar%ids(1))
            hvar%ids(1) = 0
        else
            ! 2. Variable is not a raw variable, so it is either a
            ! geometrical variable a derived variable or a gravity
            ! variable
            call check_geovar(hvar%name,hvar,ok_check)
            if (ok_check) return

            ! 3. Then the variable maybe a derived one, so we go
            ! through the defined derived variables to pre-define
            ! the indexes (for speed-up)
            call check_dervar(vardict,hvar%name,hvar,ok_check)
            if (ok_check) return

            ! 4. The option left: we are asked for a gravity variable
            call check_gravvar(hvar%name,hvar,ok_check)
            
            if (.not.ok_check) then
                write(*,*)'Variable not supported: ',TRIM(hvar%name)
                write(*,*)'Aborting!'
                stop
            end if
        end if
    end subroutine set_hydro_var

    !---------------------------------------------------------------
    ! Function: SF EFF
    !
    ! Magneto-thermo-turbulent models of star formation commonly used
    ! in RAMSES simulations have non-constant star formation
    ! efficiency, so this routines obtains the actual value in the
    ! same way it is obtained in the RAMSES version from
    ! S. Martin-Alvarez
    !----------------------------------------------------------------
    function sf_eff(my_amr,my_sim,hvar,reg,dx,x,var,star_maker,use_crs)
        use vectors
        use basis_representations
        use geometrical_regions
        use coordinate_systems
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        character(len=*),intent(in)                 :: star_maker
        logical,intent(in)  :: use_crs
        real(dbl) :: sf_eff
        ! SF variables
        logical :: isConvergent
        real(dbl) :: n_gmc,n_dc,nCOM,d_gmc,d_dc
        real(dbl) :: nISM,n_star
        real(dbl) :: t0,t_star,del_star
        real(dbl) :: d,d0
        real(dbl) :: factG
        real(dbl),dimension(0:my_amr%twondim) ::darr,uarr,varr,warr
        real(dbl) :: divv,uavg,vavg,wavg,dtot
        real(dbl) :: px_div,py_div,pz_div,Jx,Jy,Jz,rho_local
        real(dbl) :: trgv,ul,ur
        real(dbl) :: temp,c_s2
        real(dbl) :: Bx,By,Bz,Pmag,invbeta
        real(dbl) :: lamjt,sf_lam
        real(dbl) :: e_cts,phi_t,theta,alpha0
        real(dbl) :: betafunc,scrit,sigs

        ! TODO: These parameters should NOT be hardcoded!
        n_gmc = 10D0
        n_dc = 1D10
        t_star = 0.632456 ! Star formation time scale (in Gyr) at the density threshold. Def: 0.632456
        del_star=200.     ! EOS density scale
        sf_lam = 1D0

        ! 1. VARIABLE INITIALISATION
        sf_eff = 0D0
        ! Star formation time scale from Gyr to code units
        t0    = t_star*(1D0/s2Gyr)/my_sim%unit_t
        ! ISM density threshold from H/cc to code units
        nISM = n_star
        nCOM  = del_star*my_sim%omega_b*rhoc*(my_sim%h0/100.)**2/my_sim%aexp**3*XH/mHydrogen
        nISM  = max(nCOM,nISM)
        d0    = nISM/my_sim%nH
        factG = 1d0
        if (my_sim%cosmo) factG = 3d0/4d0/twopi*my_sim%omega_m*my_sim%aexp
        d_gmc = max(nCOM,n_gmc)/my_sim%nH
        d_dc  = max(nCOM,n_dc)/my_sim%nH    ! dark cloud

        ! 2. COMPUTE SF MODEL
        if (TRIM(star_maker)=='FKmag') then
            ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
            d = var(0,hvar%ids(10))
            if (d<d_gmc) return
            ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
            ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor. 
            ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation 
            ! from neighbouring cell values and differentiate. 
            ! Get neighbor cells if they exist, otherwise use straight injection from local cell

            darr = var(:,hvar%ids(10))

            uarr = var(:,hvar%ids(7))
            varr = var(:,hvar%ids(8))
            warr = var(:,hvar%ids(9))
            divv  = (uarr(2)*darr(2)-uarr(1)*darr(1)) &
                & + (varr(4)*darr(4)-varr(3)*darr(3)) & 
                & + (warr(6)*darr(6)-warr(5)*darr(5))
            if (divv>0) return
            
            ! Average velocity
            dtot  = sum(darr)
            uavg  = sum(darr*uarr)/dtot
            vavg  = sum(darr*varr)/dtot
            wavg  = sum(darr*warr)/dtot
            ! Subtract the mean velocity field
            uarr(:) = uarr(:) - uavg
            varr(:) = varr(:) - vavg
            warr(:) = warr(:) - wavg
            ! Subtract the symmetric divergence field                    
            ! ex)  (---->,<--): only subtract (-->,<--): result (-->,0) 
            ! ex)  (<----,-->): only subtract (<--,-->): result (<--,0)
            px_div = min( abs(darr(1)*uarr(1)),abs(darr(2)*uarr(2)))
            py_div = min( abs(darr(3)*varr(3)),abs(darr(4)*varr(4)))
            pz_div = min( abs(darr(5)*warr(5)),abs(darr(6)*warr(6)))

            isConvergent = darr(2)*uarr(2) - darr(1)*uarr(1) < 0 
            if (isConvergent) then
                uarr(1) = uarr(1) - px_div/darr(1)
                uarr(2) = uarr(2) + px_div/darr(2)
            else ! comment out if you do not want to subtract outflows
                uarr(1) = uarr(1) + px_div/darr(1)
                uarr(2) = uarr(2) - px_div/darr(2)
            end if 

            isConvergent = darr(4)*varr(4) - darr(3)*varr(3) < 0
            if (isConvergent) then
                varr(3) = varr(3) - py_div/darr(3)
                varr(4) = varr(4) + py_div/darr(4)
            else ! comment out if you do not want to subtract outflows
                varr(3) = varr(3) + py_div/darr(3)
                varr(4) = varr(4) - py_div/darr(4)
            end if

            isConvergent = darr(6)*warr(6) - darr(5)*warr(5) < 0
            if (isConvergent) then 
                warr(5) = warr(5) - pz_div/darr(5)
                warr(6) = warr(6) + pz_div/darr(6)
            else ! comment out if you do not want to subtract outflows
                warr(5) = warr(5) + pz_div/darr(5)
                warr(6) = warr(6) - pz_div/darr(6)
            end if

            ! subtract the rotational velocity field (x-y) plane
            ! ^y       <-        |4|        |-u|
            ! |       |  |     |1| |2|   |-v|  |+v|
            ! --->x    ->        |3|        |+u|
            Jz  = - varr(1)*darr(1) + varr(2)*darr(2) &
                &   + uarr(3)*darr(3) - uarr(4)*darr(4)
            Jz  = Jz / 4.0

            varr(1) = varr(1) + Jz/darr(1) 
            varr(2) = varr(2) - Jz/darr(2) 
            uarr(3) = uarr(3) - Jz/darr(3)
            uarr(4) = uarr(4) + Jz/darr(4)

            ! subtract the rotational velocity field (y-z) plane
            ! ^z       <-        |6|        |-v|  
            ! |       |  |     |3| |4|   |-w|  |+w|
            ! --->y    ->        |5|        |+v|
            Jx  = - warr(3)*darr(3) + warr(4)*darr(4) &
                &   + varr(5)*darr(5) - varr(6)*darr(6)
            Jx  = Jx / 4.0

            warr(3) = warr(3) + Jx/darr(3) 
            warr(4) = warr(4) - Jx/darr(4) 
            varr(5) = varr(5) - Jx/darr(5)
            varr(6) = varr(6) + Jx/darr(6)

            ! subtract the rotational velocity field (x-z) plane
            ! ^z       ->        |6|        |+u|  
            ! |       |  |     |1| |2|   |+w|  |-w|
            ! --->x    <-        |5|        |-u|
            Jy  = + warr(1)*darr(1) - warr(2)*darr(2) &
                &   - uarr(5)*darr(5) + uarr(6)*darr(6)
            Jy  = Jy / 4.0

            warr(1) = warr(1) - Jy/darr(1) 
            warr(2) = warr(2) + Jy/darr(2) 
            uarr(5) = uarr(5) + Jy/darr(5)
            uarr(6) = uarr(6) - Jy/darr(6)

            ! From this point, uarr,varr,warr is just the turbulent velocity
            trgv  = 0.0

            !x-direc
            ul    = (darr(2)*uarr(2) + d*uarr(0))/(darr(2)+d)
            ur    = (darr(1)*uarr(1) + d*uarr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc
            ul    = (darr(4)*varr(4) + d*varr(0))/(darr(4)+d)
            ur    = (darr(3)*varr(3) + d*varr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc
            ul    = (darr(6)*warr(6) + d*warr(0))/(darr(6)+d)
            ur    = (darr(5)*warr(5) + d*warr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc; tangential component - y
            ul    = (darr(6)*varr(6) + d*varr(0))/(darr(6)+d)
            ur    = (darr(5)*varr(5) + d*varr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc; tangential component - z
            ul    = (darr(4)*warr(4) + d*warr(0))/(darr(4)+d)
            ur    = (darr(3)*warr(3) + d*warr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc; tangential component - x
            ul    = (darr(6)*uarr(6) + d*uarr(0))/(darr(6)+d)
            ur    = (darr(5)*uarr(5) + d*uarr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !x-direc; tangential component - z
            ul    = (darr(2)*warr(2) + d*warr(0))/(darr(2)+d)
            ur    = (darr(1)*warr(1) + d*warr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc; tangential component - x
            ul    = (darr(4)*uarr(4) + d*uarr(0))/(darr(4)+d)
            ur    = (darr(3)*uarr(3) + d*uarr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !x-direc; tangential component - y
            ul    = (darr(2)*varr(2) + d*varr(0))/(darr(2)+d)
            ur    = (darr(1)*varr(1) + d*varr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2

            ! Now compute sound speed squared
            temp = var(0,hvar%ids(11)) / var(0,hvar%ids(10)) * gamma_gas
            ! TODO: This is a quick fix for negative T in CRMHD sims
            if (temp<0d0) temp = Tmin * gamma_gas
            ! TODO: This should be change to add also radiation pressure
            if (my_sim%cr.and.use_crs) then
                temp = temp + var(0,hvar%ids(12)) / var(0,hvar%ids(10)) * gamma_cr
            end if
            temp = max(temp,smallc**2)
            c_s2 = temp
            if (my_sim%mhd) then
                ! Added magnetic pressure to the support as contribution to c_s (assuming isothermal gas within the cell)
                Bx = var(0,hvar%ids(1))+var(0,hvar%ids(4))
                By = var(0,hvar%ids(2))+var(0,hvar%ids(5))
                Bz = var(0,hvar%ids(3))+var(0,hvar%ids(6))
                ! 1/2 factor in Pmag comes from (Br+Bl)/2
                ! Pmag = 0.03978873577*0.25*(Bx**2 + By**2 + Bz**2) ! Pmag = B^2/(8*Pi)
                Pmag = 0.125d0*(Bx**2 + By**2 + Bz**2) ! Pmag = B^2/2
                invbeta = Pmag/(c_s2*d) ! beta = c_s^2 * rho / Pmag
                c_s2 = c_s2*(1.+invbeta) ! (c_s)_eff = c_s*(1.+beta^-1)**0.5
            end if
            ! Calculate "turbulent" Jeans length in cell units, lamjt 
            ! (see e.g. Chandrasekhar 51, Bonazzola et al 87, Federrath & Klessen 2012 eq 36)
            rho_local = d
            lamjt = (pi*trgv + sqrt(pi*pi*trgv*trgv + 36.0*pi*c_s2*factG*rho_local*dx**2))/(6.0*factG*rho_local*dx**2)
            if (lamjt > sf_lam) return ! Jeans length resolved: gas is stable
            ! print*,'Jeans length notk resolved'
            ! Jeans length not resolved --> form stars to lower density and stabilise gas
            ! corresponding virial parameter for homogeneous sphere <= 1.5 in turbulent dominated 
            ! limit and <= 0.5 in pressure dominated limit (in good agreement with observations,
            ! at least for massive (>= 10^5 M_sun) clouds see Kauffmann, Pillai, Goldsmith 2013, Fig.1)
            alpha0 = 5d0/(pi*factG*rho_local)*(trgv + c_s2)/dx**2
            ! Compute star formation efficiency per free-fall time (Federrath & Klessen 2012 eq 41)
            if (my_sim%mhd) then
                ! e_cts is the unresolved (for us) proto-stellar feedback: i.e. the dense gas core-to-star efficiency
                e_cts = 0.5  ! would be 1.0 without feedback (Federrath & Klessen 2012)
                phi_t = 0.57 ; theta = 0.33 ! best fit values of the Padoan & Nordlund multi-scale sf model to GMC simulation data
                betafunc = ((1.+0.925*invbeta**1.5)**(2./3.)) / ((1.+invbeta)**2) ! Magnetic field reduction of the critical log density
                scrit = log(0.067*betafunc/(theta**2)*alpha0*trgv/c_s2) ! best fit from Padoan & Nordlund MS model again
            else
                ! e_cts is the unresolved (for us) proto-stellar feedback: i.e. the dense gas core-to-star efficiency
                e_cts = 0.5  ! would be 1.0 without feedback (Federrath & Klessen 2012)
                phi_t = 0.57 ; theta = 0.33 ! best fit values of the Padoan & Nordlund multi-scale sf model to GMC simulation data 
                scrit = log(0.067/theta**2*alpha0*trgv/c_s2) ! best fit from Padoan & Nordlund MS model again
            end if
            sigs  = log(1.0+0.16*trgv/c_s2) ! factor 0.16 is b^2 where b=0.4 for a mixture of turbulence forcing modes
            
            sf_eff = e_cts/2.0*phi_t*exp(3.0/8.0*sigs)*(2.0-erfc((sigs-scrit)/sqrt(2.0*sigs)))
        else if (TRIM(star_maker)=='FK2') then
            ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
            d = var(0,hvar%ids(10))
            if (d<d_gmc) return
            ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
            ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor. 
            ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation 
            ! from neighbouring cell values and differentiate. 
            ! Get neighbor cells if they exist, otherwise use straight injection from local cell

            darr = var(:,hvar%ids(10))

            uarr = var(:,hvar%ids(7))
            varr = var(:,hvar%ids(8))
            warr = var(:,hvar%ids(9))
            divv  = (uarr(2)*darr(2)-uarr(1)*darr(1)) &
                & + (varr(4)*darr(4)-varr(3)*darr(3)) & 
                & + (warr(6)*darr(6)-warr(5)*darr(5))
            if (divv>0) return
            
            ! Average velocity
            dtot  = sum(darr)
            uavg  = sum(darr*uarr)/dtot
            vavg  = sum(darr*varr)/dtot
            wavg  = sum(darr*warr)/dtot
            ! Subtract the mean velocity field
            uarr(:) = uarr(:) - uavg
            varr(:) = varr(:) - vavg
            warr(:) = warr(:) - wavg
            ! Subtract the symmetric divergence field                    
            ! ex)  (---->,<--): only subtract (-->,<--): result (-->,0) 
            ! ex)  (<----,-->): only subtract (<--,-->): result (<--,0)
            px_div = min( abs(darr(1)*uarr(1)),abs(darr(2)*uarr(2)))
            py_div = min( abs(darr(3)*varr(3)),abs(darr(4)*varr(4)))
            pz_div = min( abs(darr(5)*warr(5)),abs(darr(6)*warr(6)))

            isConvergent = darr(2)*uarr(2) - darr(1)*uarr(1) < 0 
            if (isConvergent) then
                uarr(1) = uarr(1) - px_div/darr(1)
                uarr(2) = uarr(2) + px_div/darr(2)
            else ! comment out if you do not want to subtract outflows
                uarr(1) = uarr(1) + px_div/darr(1)
                uarr(2) = uarr(2) - px_div/darr(2)
            end if 

            isConvergent = darr(4)*varr(4) - darr(3)*varr(3) < 0
            if (isConvergent) then
                varr(3) = varr(3) - py_div/darr(3)
                varr(4) = varr(4) + py_div/darr(4)
            else ! comment out if you do not want to subtract outflows
                varr(3) = varr(3) + py_div/darr(3)
                varr(4) = varr(4) - py_div/darr(4)
            end if

            isConvergent = darr(6)*warr(6) - darr(5)*warr(5) < 0
            if (isConvergent) then 
                warr(5) = warr(5) - pz_div/darr(5)
                warr(6) = warr(6) + pz_div/darr(6)
            else ! comment out if you do not want to subtract outflows
                warr(5) = warr(5) + pz_div/darr(5)
                warr(6) = warr(6) - pz_div/darr(6)
            end if

            ! subtract the rotational velocity field (x-y) plane
            ! ^y       <-        |4|        |-u|
            ! |       |  |     |1| |2|   |-v|  |+v|
            ! --->x    ->        |3|        |+u|
            Jz  = - varr(1)*darr(1) + varr(2)*darr(2) &
                &   + uarr(3)*darr(3) - uarr(4)*darr(4)
            Jz  = Jz / 4.0

            varr(1) = varr(1) + Jz/darr(1) 
            varr(2) = varr(2) - Jz/darr(2) 
            uarr(3) = uarr(3) - Jz/darr(3)
            uarr(4) = uarr(4) + Jz/darr(4)

            ! subtract the rotational velocity field (y-z) plane
            ! ^z       <-        |6|        |-v|  
            ! |       |  |     |3| |4|   |-w|  |+w|
            ! --->y    ->        |5|        |+v|
            Jx  = - warr(3)*darr(3) + warr(4)*darr(4) &
                &   + varr(5)*darr(5) - varr(6)*darr(6)
            Jx  = Jx / 4.0

            warr(3) = warr(3) + Jx/darr(3) 
            warr(4) = warr(4) - Jx/darr(4) 
            varr(5) = varr(5) - Jx/darr(5)
            varr(6) = varr(6) + Jx/darr(6)

            ! subtract the rotational velocity field (x-z) plane
            ! ^z       ->        |6|        |+u|  
            ! |       |  |     |1| |2|   |+w|  |-w|
            ! --->x    <-        |5|        |-u|
            Jy  = + warr(1)*darr(1) - warr(2)*darr(2) &
                &   - uarr(5)*darr(5) + uarr(6)*darr(6)
            Jy  = Jy / 4.0

            warr(1) = warr(1) - Jy/darr(1) 
            warr(2) = warr(2) + Jy/darr(2) 
            uarr(5) = uarr(5) + Jy/darr(5)
            uarr(6) = uarr(6) - Jy/darr(6)

            ! From this point, uarr,varr,warr is just the turbulent velocity
            trgv  = 0.0

            !x-direc
            ul    = (darr(2)*uarr(2) + d*uarr(0))/(darr(2)+d)
            ur    = (darr(1)*uarr(1) + d*uarr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc
            ul    = (darr(4)*varr(4) + d*varr(0))/(darr(4)+d)
            ur    = (darr(3)*varr(3) + d*varr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc
            ul    = (darr(6)*warr(6) + d*warr(0))/(darr(6)+d)
            ur    = (darr(5)*warr(5) + d*warr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc; tangential component - y
            ul    = (darr(6)*varr(6) + d*varr(0))/(darr(6)+d)
            ur    = (darr(5)*varr(5) + d*varr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc; tangential component - z
            ul    = (darr(4)*warr(4) + d*warr(0))/(darr(4)+d)
            ur    = (darr(3)*warr(3) + d*warr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !z-direc; tangential component - x
            ul    = (darr(6)*uarr(6) + d*uarr(0))/(darr(6)+d)
            ur    = (darr(5)*uarr(5) + d*uarr(0))/(darr(5)+d)
            trgv  = trgv + (ur-ul)**2
            !x-direc; tangential component - z
            ul    = (darr(2)*warr(2) + d*warr(0))/(darr(2)+d)
            ur    = (darr(1)*warr(1) + d*warr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2
            !y-direc; tangential component - x
            ul    = (darr(4)*uarr(4) + d*uarr(0))/(darr(4)+d)
            ur    = (darr(3)*uarr(3) + d*uarr(0))/(darr(3)+d)
            trgv  = trgv + (ur-ul)**2
            !x-direc; tangential component - y
            ul    = (darr(2)*varr(2) + d*varr(0))/(darr(2)+d)
            ur    = (darr(1)*varr(1) + d*varr(0))/(darr(1)+d)
            trgv  = trgv + (ur-ul)**2

            ! Now compute sound speed squared
            temp = var(0,hvar%ids(11)) / var(0,hvar%ids(10)) * gamma_gas
            temp = max(temp,smallc**2)
            c_s2 = temp
            ! Calculate "turbulent" Jeans length in cell units, lamjt 
            ! (see e.g. Chandrasekhar 51, Bonazzola et al 87, Federrath & Klessen 2012 eq 36)
            rho_local = d
            lamjt = (pi*trgv + sqrt(pi*pi*trgv*trgv + 36.0*pi*c_s2*factG*rho_local*dx**2))/(6.0*factG*rho_local*dx**2)
            if (lamjt > sf_lam) return ! Jeans length resolved: gas is stable
            ! print*,'Jeans length notk resolved'
            ! Jeans length not resolved --> form stars to lower density and stabilise gas
            ! corresponding virial parameter for homogeneous sphere <= 1.5 in turbulent dominated 
            ! limit and <= 0.5 in pressure dominated limit (in good agreement with observations,
            ! at least for massive (>= 10^5 M_sun) clouds see Kauffmann, Pillai, Goldsmith 2013, Fig.1)
            alpha0 = 5d0/(pi*factG*rho_local)*(trgv + c_s2)/dx**2
            ! Compute star formation efficiency per free-fall time (Federrath & Klessen 2012 eq 41)
            ! e_cts is the unresolved (for us) proto-stellar feedback: i.e. the dense gas core-to-star efficiency
            e_cts = 0.5  ! would be 1.0 without feedback (Federrath & Klessen 2012)
            phi_t = 0.57 ; theta = 0.33 ! best fit values of the Padoan & Nordlund multi-scale sf model to GMC simulation data 
            sigs  = log(1.0+0.16*trgv/c_s2) ! factor 0.16 is b^2 where b=0.4 for a mixture of turbulence forcing modes
            scrit = log(0.067/theta**2*alpha0*trgv/c_s2) ! best fit from Padoan & Nordlund MS model again 
            sf_eff = e_cts/2.0*phi_t*exp(3.0/8.0*sigs)*(2.0-erfc((sigs-scrit)/sqrt(2.0*sigs)))
        end if

        
    end function sf_eff

    !---------------------------------------------------------------
    ! Function: CMP SIGMA TURB
    !
    ! This function obtains the local dispersion velocity of the gas
    ! in the same way it is used for the MTT models of star formation
    !----------------------------------------------------------------
    subroutine cmp_sigma_turb(my_amr,my_sim,hvar,var,sg)
        use vectors
        use basis_representations
        use geometrical_regions
        use coordinate_systems

        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        real(dbl),intent(out)           :: sg

        logical :: isConvergent
        real(dbl) :: d
        real(dbl),dimension(0:my_amr%twondim) ::darr,uarr,varr,warr
        real(dbl) :: divv,uavg,vavg,wavg,dtot
        real(dbl) :: px_div,py_div,pz_div,Jx,Jy,Jz,rho_local
        real(dbl) :: trgv,ul,ur

        ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
        ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor. 
        ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation 
        ! from neighbouring cell values and differentiate. 
        ! Get neighbor cells if they exist, otherwise use straight injection from local cell

        d = var(0,hvar%ids(1))
        darr = var(:,hvar%ids(1))

        uarr = var(:,hvar%ids(2))
        varr = var(:,hvar%ids(3))
        warr = var(:,hvar%ids(4))
        divv  = (uarr(2)*darr(2)-uarr(1)*darr(1)) &
            & + (varr(4)*darr(4)-varr(3)*darr(3)) & 
            & + (warr(6)*darr(6)-warr(5)*darr(5))
        
        ! Average velocity
        dtot  = sum(darr)
        uavg  = sum(darr*uarr)/dtot
        vavg  = sum(darr*varr)/dtot
        wavg  = sum(darr*warr)/dtot
        ! Subtract the mean velocity field
        uarr(:) = uarr(:) - uavg
        varr(:) = varr(:) - vavg
        warr(:) = warr(:) - wavg
        ! Subtract the symmetric divergence field                    
        ! ex)  (---->,<--): only subtract (-->,<--): result (-->,0) 
        ! ex)  (<----,-->): only subtract (<--,-->): result (<--,0)
        px_div = min( abs(darr(1)*uarr(1)),abs(darr(2)*uarr(2)))
        py_div = min( abs(darr(3)*varr(3)),abs(darr(4)*varr(4)))
        pz_div = min( abs(darr(5)*warr(5)),abs(darr(6)*warr(6)))

        isConvergent = darr(2)*uarr(2) - darr(1)*uarr(1) < 0 
        if (isConvergent) then
            uarr(1) = uarr(1) - px_div/darr(1)
            uarr(2) = uarr(2) + px_div/darr(2)
        else ! comment out if you do not want to subtract outflows
            uarr(1) = uarr(1) + px_div/darr(1)
            uarr(2) = uarr(2) - px_div/darr(2)
        end if 

        isConvergent = darr(4)*varr(4) - darr(3)*varr(3) < 0
        if (isConvergent) then
            varr(3) = varr(3) - py_div/darr(3)
            varr(4) = varr(4) + py_div/darr(4)
        else ! comment out if you do not want to subtract outflows
            varr(3) = varr(3) + py_div/darr(3)
            varr(4) = varr(4) - py_div/darr(4)
        end if

        isConvergent = darr(6)*warr(6) - darr(5)*warr(5) < 0
        if (isConvergent) then 
            warr(5) = warr(5) - pz_div/darr(5)
            warr(6) = warr(6) + pz_div/darr(6)
        else ! comment out if you do not want to subtract outflows
            warr(5) = warr(5) + pz_div/darr(5)
            warr(6) = warr(6) - pz_div/darr(6)
        end if

        ! subtract the rotational velocity field (x-y) plane
        ! ^y       <-        |4|        |-u|
        ! |       |  |     |1| |2|   |-v|  |+v|
        ! --->x    ->        |3|        |+u|
        Jz  = - varr(1)*darr(1) + varr(2)*darr(2) &
            &   + uarr(3)*darr(3) - uarr(4)*darr(4)
        Jz  = Jz / 4.0

        varr(1) = varr(1) + Jz/darr(1) 
        varr(2) = varr(2) - Jz/darr(2) 
        uarr(3) = uarr(3) - Jz/darr(3)
        uarr(4) = uarr(4) + Jz/darr(4)

        ! subtract the rotational velocity field (y-z) plane
        ! ^z       <-        |6|        |-v|  
        ! |       |  |     |3| |4|   |-w|  |+w|
        ! --->y    ->        |5|        |+v|
        Jx  = - warr(3)*darr(3) + warr(4)*darr(4) &
            &   + varr(5)*darr(5) - varr(6)*darr(6)
        Jx  = Jx / 4.0

        warr(3) = warr(3) + Jx/darr(3) 
        warr(4) = warr(4) - Jx/darr(4) 
        varr(5) = varr(5) - Jx/darr(5)
        varr(6) = varr(6) + Jx/darr(6)

        ! subtract the rotational velocity field (x-z) plane
        ! ^z       ->        |6|        |+u|  
        ! |       |  |     |1| |2|   |+w|  |-w|
        ! --->x    <-        |5|        |-u|
        Jy  = + warr(1)*darr(1) - warr(2)*darr(2) &
            &   - uarr(5)*darr(5) + uarr(6)*darr(6)
        Jy  = Jy / 4.0

        warr(1) = warr(1) - Jy/darr(1) 
        warr(2) = warr(2) + Jy/darr(2) 
        uarr(5) = uarr(5) + Jy/darr(5)
        uarr(6) = uarr(6) - Jy/darr(6)

        ! From this point, uarr,varr,warr is just the turbulent velocity
        trgv  = 0.0

        !x-direc
        ul    = (darr(2)*uarr(2) + d*uarr(0))/(darr(2)+d)
        ur    = (darr(1)*uarr(1) + d*uarr(0))/(darr(1)+d)
        trgv  = trgv + (ur-ul)**2
        !y-direc
        ul    = (darr(4)*varr(4) + d*varr(0))/(darr(4)+d)
        ur    = (darr(3)*varr(3) + d*varr(0))/(darr(3)+d)
        trgv  = trgv + (ur-ul)**2
        !z-direc
        ul    = (darr(6)*warr(6) + d*warr(0))/(darr(6)+d)
        ur    = (darr(5)*warr(5) + d*warr(0))/(darr(5)+d)
        trgv  = trgv + (ur-ul)**2
        !z-direc; tangential component - y
        ul    = (darr(6)*varr(6) + d*varr(0))/(darr(6)+d)
        ur    = (darr(5)*varr(5) + d*varr(0))/(darr(5)+d)
        trgv  = trgv + (ur-ul)**2
        !y-direc; tangential component - z
        ul    = (darr(4)*warr(4) + d*warr(0))/(darr(4)+d)
        ur    = (darr(3)*warr(3) + d*warr(0))/(darr(3)+d)
        trgv  = trgv + (ur-ul)**2
        !z-direc; tangential component - x
        ul    = (darr(6)*uarr(6) + d*uarr(0))/(darr(6)+d)
        ur    = (darr(5)*uarr(5) + d*uarr(0))/(darr(5)+d)
        trgv  = trgv + (ur-ul)**2
        !x-direc; tangential component - z
        ul    = (darr(2)*warr(2) + d*warr(0))/(darr(2)+d)
        ur    = (darr(1)*warr(1) + d*warr(0))/(darr(1)+d)
        trgv  = trgv + (ur-ul)**2
        !y-direc; tangential component - x
        ul    = (darr(4)*uarr(4) + d*uarr(0))/(darr(4)+d)
        ur    = (darr(3)*uarr(3) + d*uarr(0))/(darr(3)+d)
        trgv  = trgv + (ur-ul)**2
        !x-direc; tangential component - y
        ul    = (darr(2)*varr(2) + d*varr(0))/(darr(2)+d)
        ur    = (darr(1)*varr(1) + d*varr(0))/(darr(1)+d)
        trgv  = trgv + (ur-ul)**2
        
        ! we want the square root of this
        sg = sqrt(trgv)
    end subroutine cmp_sigma_turb

end module hydro_commons