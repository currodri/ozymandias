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
    use io_ramses, only: amr_info, sim_info, Tmin, cV, lambda_crGH08
    use cooling_module

    type hydro_var
        character(128) :: name,type
        integer, dimension(:), allocatable :: ids
        procedure(myinterface),pointer,nopass :: myfunction
    end type hydro_var

    contains

    ! RAW VARIABLES
    function raw_hydro(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: raw_hydro
        raw_hydro = var(0,hvar%ids(1))
    end function raw_hydro

    ! GEOMETRICAL VARIABLES

    function myinterface(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: myinterface
        myinterface = magnitude(x)
    end function myinterface

    function d_euclid_wrap(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: d_euclid_wrap

        ! Euclidean distance
        d_euclid_wrap = magnitude(x)
    end function d_euclid_wrap

    function x_coord_wrap(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: x_coord_wrap

        ! x - coordinate
        x_coord_wrap = x%x
    end function x_coord_wrap

    function y_coord_wrap(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: y_coord_wrap

        ! y - coordinate
        y_coord_wrap = x%y
    end function y_coord_wrap

    function z_coord_wrap(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: z_coord_wrap

        ! z - coordinate
        z_coord_wrap = x%z
    end function z_coord_wrap

    function r_sphere_wrap(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: r_sphere_wrap

        ! Radius from center of sphere
        r_sphere_wrap = r_sphere(x)
    end function r_sphere_wrap

    function theta_sphere_wrap(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: theta_sphere_wrap

        ! Value of spherical theta angle measured from the z axis
        theta_sphere_wrap = theta_sphere(x)
    end function theta_sphere_wrap

    function phi_sphere_wrap(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: phi_sphere_wrap

        ! Value of spherical phi angle measure in the x-y plane 
        ! from the x axis
        phi_sphere_wrap = phi_sphere(x)
    end function phi_sphere_wrap

    function r_cyl_wrap(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: r_cyl_wrap

        ! Value of cylindrical radius
        r_cyl_wrap = r_cyl(x)
    end function r_cyl_wrap

    function phi_cyl_wrap(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_potential(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: grav_potential

        ! Gravitational potential
        grav_potential = grav_var(0,1)
    end function grav_potential

    function grav_gx(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: grav_gx

        ! Gravitational acceleration in the x direction
        grav_gx = grav_var(0,2)
    end function grav_gx

    function grav_gy(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: grav_gy

        ! Gravitational acceleration in the y direction
        grav_gy = grav_var(0,3)
    end function grav_gy

    function grav_gz(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function cell_volume(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: cell_volume

        cell_volume = (dx * dx) * dx
    end function cell_volume

    function cell_mass(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: cell_mass

        cell_mass = (var(0,hvar%ids(1)) * (dx * dx)) * dx
    end function cell_mass

    function v_sphere_r(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function v_sphere_phi(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function v_sphere_theta(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function v_cyl_r(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function v_cyl_z(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function v_cyl_phi(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function v_magnitude(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: v_magnitude
        type(vector) :: v
        type(basis) :: temp_basis

        ! Velocity magnitude from galaxy coordinates
        v = (/var(0,hvar%ids(1)),var(0,hvar%ids(2)),var(0,hvar%ids(3))/)
        v_magnitude = magnitude(v)
    end function v_magnitude

    function v_squared(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: v_squared
        type(vector) :: v
        type(basis) :: temp_basis

        ! Velocity magnitude squared from galaxy coordinates
        v = (/var(0,hvar%ids(1)),var(0,hvar%ids(2)),var(0,hvar%ids(3))/)
        v_squared = v.DOT.v
    end function v_squared

    function momentum_x(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: momentum_x
        type(vector) :: v
        type(basis) :: temp_basis

        ! Linear momentum in the x direction as density*volume*corrected_velocity_x
        momentum_x = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * var(0,hvar%ids(2))
    end function momentum_x

    function momentum_y(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: momentum_y
        type(vector) :: v
        type(basis) :: temp_basis

        ! Linear momentum in the y direction as density*volume*corrected_velocity_y
        momentum_y = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * var(0,hvar%ids(2))
    end function momentum_y

    function momentum_z(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: momentum_z
        type(vector) :: v
        type(basis) :: temp_basis

        ! Linear momentum in the z direction as density*volume*corrected_velocity_z
        momentum_z = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * var(0,hvar%ids(2))
    end function momentum_z

    function momentum(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: momentum
        type(vector) :: v
        type(basis) :: temp_basis

        ! Magnitude of linear momentum, using corrected velocity
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        momentum = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * &
                    & magnitude(v)
    end function momentum

    function momentum_sphere_r(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function momentum_cyl_z(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function ang_momentum_x(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: ang_momentum_x
        type(vector) :: v
        type(basis) :: temp_basis

        ! Corrected angular momentum in the x direction
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        ang_momentum_x = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * (x%y * v%z &
                            &- v%y * x%z)
    end function ang_momentum_x

    function ang_momentum_y(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: ang_momentum_y
        type(vector) :: v
        type(basis) :: temp_basis

        ! Corrected angular momentum in the y direction
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        ang_momentum_y = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * (x%z*v%x &
                            &- v%z*x%x)
    end function ang_momentum_y

    function ang_momentum_z(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: ang_momentum_z
        type(vector) :: v
        type(basis) :: temp_basis

        ! Corrected angular momentum in the z direction
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        ang_momentum_z = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * (x%x*v%y &
                        &- v%x*x%y)
    end function ang_momentum_z

    function ang_momentum(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: ang_momentum
        type(vector) :: v,L
        type(basis) :: temp_basis

        ! Corrected angular momentum in the z direction
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        L = x * v
        ang_momentum = ((var(0,hvar%ids(1)) * (dx*dx)) * dx) * magnitude(L)
    end function ang_momentum

    function massflow_rate_sphere_r(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: massflow_rate_sphere_r
        type(vector) :: v
        type(basis) :: temp_basis

        ! Mass flow rate through the cell in the radial direction
            ! Mass per unit time
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        call spherical_basis_from_cartesian(x,temp_basis)
        massflow_rate_sphere_r = (var(0,hvar%ids(1)) * (dx*dx)) * (v .DOT. temp_basis%u(1))    
    end function massflow_rate_sphere_r

    function massflux_rate_sphere_r(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: massflux_rate_sphere_r
        type(vector) :: v
        type(basis) :: temp_basis

        ! Mass flux through the cell in the radial direction
        ! Mass per unit time per unit surface
        v = (/var(0,hvar%ids(2)),var(0,hvar%ids(3)),var(0,hvar%ids(4))/)
        call spherical_basis_from_cartesian(x,temp_basis)
        massflux_rate_sphere_r = var(0,hvar%ids(1)) * (v .DOT. temp_basis%u(1))  
    end function massflux_rate_sphere_r

    function kinetic_energy(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: kinetic_energy

        ! Kinetic energy, computed as 1/2*density*volume*magnitude(velocity)
        ! NOTE: Not corrected!
        kinetic_energy = (0.5d0 * (var(0,hvar%ids(1)) * (dx*dx)) * dx) &
                                 & * sqrt(var(0,hvar%ids(2))**2 + var(0,hvar%ids(3))**2 &
                                 & + var(0,hvar%ids(4))**2)
    end function kinetic_energy

    function sigma(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: sigma,value
        type(vector) :: v
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar) :: tempvar

        ! Local velocity dispersion
        ! Go back to box coordinates for the central cell, which is transformed usually
        ! before sent to read_amr
        ! Converging flow check
        tempvar(:,:) = var(:,:)
        v = (/tempvar(0,hvar%ids(2)),tempvar(0,hvar%ids(3)),tempvar(0,hvar%ids(4))/)
        call rotate_vector(v,transpose(trans_matrix))
        v = v + reg%bulk_velocity
        tempvar(0,hvar%ids(2)) = v%x
        tempvar(0,hvar%ids(3)) = v%y
        tempvar(0,hvar%ids(4)) = v%z
        call cmp_sigma_turb(my_amr,my_sim,hvar,tempvar,value)
        sigma = value
    end function sigma

    function temperature(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: temperature

        ! Gas temperature
        temperature = var(0,hvar%ids(2)) / var(0,hvar%ids(1))
        if (temperature < 0d0) then
            temperature = Tmin
        endif
    end function temperature

    function thermal_energy(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: thermal_energy

        ! Thermal energy, computed as thermal_pressure*volume/(gamma - 1)
        thermal_energy = ((max(var(0,hvar%ids(2)), Tmin*var(0,hvar%ids(1))) &
                & / (gamma_gas - 1d0)) * (dx * dx)) * dx
    end function thermal_energy

    function thermal_energy_specific(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: thermal_energy_specific

        ! Specific thermal energy as E_ther/cell mass
        thermal_energy_specific = (max(var(0,hvar%ids(2)), Tmin*var(0,hvar%ids(1))) &
                                    & / (gamma_gas - 1d0)) / var(0,hvar%ids(1))
    end function thermal_energy_specific

    function thermal_energy_density(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: thermal_energy_density

        ! Thermal energy density  as E_ther/cell volume
        thermal_energy_density = (max(var(0,hvar%ids(2)), Tmin*var(0,hvar%ids(1))) &
                                    & / (gamma_gas - 1d0))
    end function thermal_energy_density

    function entropy_specific(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: entropy_specific

        real(dbl) :: T,rho

        ! Specific entropy, following Gent 2012 equation
        T = var(0,hvar%ids(2))/var(0,hvar%ids(1)) * my_sim%T2 * mu
        rho = var(0,hvar%ids(1)) * my_sim%unit_d / mHydrogen
        entropy_specific = log(T) - (gamma_gas-1d0) * log(rho)
    end function entropy_specific

    function sound_speed(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: sound_speed

        ! Thermal sound speed, ideal gas
        sound_speed = sqrt(gamma_gas * (max(var(0,hvar%ids(2)), Tmin*var(0,hvar%ids(1))) &
                            & / var(0,hvar%ids(1))))
    end function sound_speed

    function grad_thermalpressure(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grad_therprsphere(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grad_therpz(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function magnetic_energy(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: magnetic_energy

        type(vector) :: B

        ! Magnetic energy as magnitude(B)**2/2
        B = 0.5 *(/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                    (var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                    (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
        magnetic_energy = (0.5 * (B.DOT.B) * (dx*dx)) * dx
    end function magnetic_energy

    function magnetic_magnitude(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: magnetic_magnitude

        type(vector) :: B

        ! Magnetic field magnitude
        B = 0.5 *(/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                    (var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                    (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
        magnetic_magnitude = magnitude(B)
    end function magnetic_magnitude

    function magnetic_energy_specific(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: magnetic_energy_specific

        type(vector) :: B

        ! Magnetic field magnitude
        B = 0.5 *(/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                    (var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                    (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
        magnetic_energy_specific = (0.5 * (B.DOT.B)) / var(0,hvar%ids(7))
    end function magnetic_energy_specific

    function magnetic_energy_density(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: magnetic_energy_density

        type(vector) :: B

        ! Magnetic field magnitude
        B = 0.5 *(/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                    (var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                    (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
        magnetic_energy_density = 0.5 * (B.DOT.B)
    end function magnetic_energy_density

    function alfven_speed(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: alfven_speed

        type(vector) :: B

        ! Magnetic field magnitude
        B = 0.5 *(/(var(0,hvar%ids(1))+var(0,hvar%ids(4))),&
                    (var(0,hvar%ids(2))+var(0,hvar%ids(5))),&
                    (var(0,hvar%ids(3))+var(0,hvar%ids(6)))/)
        alfven_speed = magnitude(B) / sqrt(var(0,hvar%ids(7)))
    end function alfven_speed

    function grav_magpfrsphere(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_magpfrspherepos(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_magpfrsphereneg(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function cr_energy(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: cr_energy

        ! CR energy, computed as CR_energydensity*volume
        cr_energy = ((var(0,hvar%ids(1)) / (gamma_cr - 1d0)) * (dx*dx)) * dx
    end function cr_energy

    function cr_energy_density(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: cr_energy_density

        ! CR energy density
        cr_energy_density = var(0,hvar%ids(1)) / (gamma_cr - 1d0)
    end function cr_energy_density

    function cr_energy_specific(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: cr_energy_specific

        ! Specific CR energy, computed as CR_energydensity*volume/cell mass
        cr_energy_specific = (var(0,hvar%ids(1)) / (gamma_cr - 1d0)) / var(0,hvar%ids(2))
    end function cr_energy_specific

    function cr_temperature_eff(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: cr_temperature_eff

        ! Effective CR temperature
        cr_temperature_eff = var(0,hvar%ids(1)) / var(0,hvar%ids(2))
    end function cr_temperature_eff

    function cr_GH08heat(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grad_crp(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grad_crpx(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grad_crpy(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grad_crpz(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grad_crprsphere(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function gradscale_crprsphere(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function gradscale_crp(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function diffusion_speed(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function alfvendiff_ratio(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function streaming_heating(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function net_cooling(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var
    
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

    function stheatcooling_ratio(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function total_coolingtime(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_frsphere(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: grav_frsphere
        type(vector) :: B
        type(basis) :: temp_basis

        B = grav_var(0,2:4)
        call spherical_basis_from_cartesian(x,temp_basis)
        grav_frsphere = var(0,hvar%ids(1)) * (B .DOT. temp_basis%u(1))
    end function grav_frsphere

    function escape_velocity(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: escape_velocity

        escape_velocity = sqrt(2d0*grav_var(0,1))
    end function escape_velocity

    function grav_crpf(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_crpfz(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_therpfz(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_therpfrsphere(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_therpfrspherepos(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_therpfrsphereneg(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_crpfrsphere(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_crpfrspherepos(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_crpfrsphereneg(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_totpfrsphere(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_totpfrspherepos(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function grav_totpfrsphereneg(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function eff_FKmag(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function eff_FKmagnocr(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function eff_FK2(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    ! DUST VARIABLES
    function PAHSmall_density(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: PAHSmall_density

        PAHSmall_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function PAHSmall_density

    function PAHLarge_density(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: PAHLarge_density

        PAHLarge_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function PAHLarge_density

    function CSmall_density(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: CSmall_density

        CSmall_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function CSmall_density

    function CLarge_density(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: CLarge_density

        CLarge_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function CLarge_density

    function SilSmall_density(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: SilSmall_density

        SilSmall_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function SilSmall_density

    function SilLarge_density(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: SilLarge_density

        SilLarge_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function SilLarge_density

    function CO_density(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: CO_density

        CO_density = var(0,hvar%ids(1)) * var(0,hvar%ids(2))
    end function CO_density

    function PAHSmall_mass(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: PAHSmall_mass

        PAHSmall_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function PAHSmall_mass

    function PAHLarge_mass(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: PAHLarge_mass

        PAHLarge_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function PAHLarge_mass

    function CSmall_mass(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: CSmall_mass

        CSmall_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function CSmall_mass

    function CLarge_mass(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: CLarge_mass

        CLarge_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function CLarge_mass

    function SilSmall_mass(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: SilSmall_mass

        SilSmall_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function SilSmall_mass

    function SilLarge_mass(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: SilLarge_mass

        SilLarge_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function SilLarge_mass

    function CO_mass(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: CO_mass

        CO_mass = ((var(0,hvar%ids(1)) * var(0,hvar%ids(2)) * dx) * dx) * dx
    end function CO_mass

    function DTM(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: DTM
        real(dbl) :: total_dust_mass, total_metal_mass

        total_dust_mass = var(0,hvar%ids(1)) + var(0,hvar%ids(2)) + var(0,hvar%ids(3)) + var(0,hvar%ids(4))
        total_metal_mass = var(0,hvar%ids(5))
        DTM = total_dust_mass / (total_dust_mass + total_metal_mass)
    end function DTM

    function STL(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: STL
        real(dbl) :: total_small_grains, total_large_grains

        total_small_grains = var(0,hvar%ids(1)) + var(0,hvar%ids(3))
        total_large_grains = var(0,hvar%ids(2)) + var(0,hvar%ids(4))
        STL = total_small_grains / total_large_grains
    end function STL

    function CSR(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

        real(dbl) :: CSR
        real(dbl) :: total_carbon, total_silicon

        total_carbon = var(0,hvar%ids(1)) + var(0,hvar%ids(2))
        total_silicon = var(0,hvar%ids(3)) + var(0,hvar%ids(4))
        CSR = total_carbon / total_silicon
    end function CSR

    function CSmall_Stokes(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function CLarge_Stokes(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function SilSmall_Stokes(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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

    function SilLarge_Stokes(my_amr,my_sim,hvar,reg,dx,x,var,son,trans_matrix,grav_var)
        implicit none
        type(amr_info),intent(in) :: my_amr
        type(sim_info),intent(in) :: my_sim
        type(hydro_var), intent(in) :: hvar
        type(region),intent(in)                       :: reg
        real(dbl),intent(in)                       :: dx
        type(vector),intent(in)        :: x
        real(dbl),dimension(0:my_amr%twondim,1:my_sim%nvar),intent(in) :: var
        integer,dimension(0:my_amr%twondim),intent(in) :: son
        real(dbl),dimension(1:3,1:3),optional,intent(in) :: trans_matrix
        real(dbl),dimension(0:my_amr%twondim,1:4),optional,intent(in) :: grav_var

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
        case ('PAHSmall_density')
            ! Density of the small PAHs
            hvar%type = 'derived'
            hvar%name = 'PAHSmall_density'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('PAHSmall_fraction')
            hvar%myfunction => PAHSmall_density
        case ('PAHLarge_density')
            ! Density of the large PAHs
            hvar%type = 'derived'
            hvar%name = 'PAHLarge_density'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('PAHLarge_fraction')
            hvar%myfunction => PAHLarge_density
        case ('CSmall_density')
            ! Density of the small C grains
            hvar%type = 'derived'
            hvar%name = 'CSmall_density'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('CSmall_fraction')
            hvar%myfunction => CSmall_density
        case ('CLarge_density')
            ! Density of the large PAHs
            hvar%type = 'derived'
            hvar%name = 'CLarge_density'
            allocate(hvar%ids(2))
            hvar%ids(1) = vardict%get('density')
            hvar%ids(2) = vardict%get('CLarge_fraction')
            hvar%myfunction => CLarge_density
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
        case ('DTM')
            ! dust-to-metal mass ratio
            hvar%type = 'derived'
            hvar%name = 'DTM'
            allocate(hvar%ids(5))
            hvar%ids(1) = vardict%get('CSmall_fraction')
            hvar%ids(2) = vardict%get('CLarge_fraction')
            hvar%ids(3) = vardict%get('SilSmall_fraction')
            hvar%ids(4) = vardict%get('SilLarge_fraction')
            hvar%ids(5) = vardict%get('metallicity')
            hvar%myfunction => DTM
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
            birth_date = part_var_d(pvar%ids(1))
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
            birth_date = part_var_d(pvar%ids(1))
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
            birth_date = part_var_d(pvar%ids(1))
            current_age_univ = my_sim%time_simu*my_sim%unit_t/(365.*24.*3600.*1d9)
            ! 3. If the particle is older than the indicator, the SFR is zero
            if (birth_date >= (current_age_univ - sfrind)) then
                sfr_surface = part_var_d(pvar%ids(2)) / (dx%x*dx%y) / sfrind
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
            print*,get_numeric_suffix(varname),pvar%num_suffix
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