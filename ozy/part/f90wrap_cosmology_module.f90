! Module cosmology defined in file cosmology_module.fpp

subroutine f90wrap_cosmology_model(sim)
    use cosmology, only: cosmology_model
    use io_ramses, only: sim_info
    implicit none
    
    type sim_info_ptr_type
        type(sim_info), pointer :: p => NULL()
    end type sim_info_ptr_type
    type(sim_info_ptr_type) :: sim_ptr
    integer, intent(in), dimension(2) :: sim
    sim_ptr = transfer(sim, sim_ptr)
    call cosmology_model(sim=sim_ptr%p)
end subroutine f90wrap_cosmology_model

subroutine f90wrap_friedmann(o_mat_0, o_vac_0, o_k_0, alpha, axp_min, axp_out, hexp_out, tau_out, t_out, ntable, &
    age_tot, n0, n1, n2, n3)
    use cosmology, only: friedmann
    implicit none
    
    real(8), intent(in) :: o_mat_0
    real(8), intent(in) :: o_vac_0
    real(8), intent(in) :: o_k_0
    real(8), intent(in) :: alpha
    real(8), intent(in) :: axp_min
    real(8), intent(inout), dimension(n0) :: axp_out
    real(8), intent(inout), dimension(n1) :: hexp_out
    real(8), intent(inout), dimension(n2) :: tau_out
    real(8), intent(inout), dimension(n3) :: t_out
    integer, intent(in) :: ntable
    real(8), intent(inout) :: age_tot
    integer :: n0
    !f2py intent(hide), depend(axp_out) :: n0 = shape(axp_out,0)
    integer :: n1
    !f2py intent(hide), depend(hexp_out) :: n1 = shape(hexp_out,0)
    integer :: n2
    !f2py intent(hide), depend(tau_out) :: n2 = shape(tau_out,0)
    integer :: n3
    !f2py intent(hide), depend(t_out) :: n3 = shape(t_out,0)
    call friedmann(O_mat_0=o_mat_0, O_vac_0=o_vac_0, O_k_0=o_k_0, alpha=alpha, axp_min=axp_min, axp_out=axp_out, &
        hexp_out=hexp_out, tau_out=tau_out, t_out=t_out, ntable=ntable, age_tot=age_tot)
end subroutine f90wrap_friedmann

subroutine f90wrap_dadtau(axp_tau, o_mat_0, o_vac_0, ret_dadtau, o_k_0)
    use cosmology, only: dadtau
    implicit none
    
    real(8) :: axp_tau
    real(8) :: o_mat_0
    real(8) :: o_vac_0
    real(8), intent(out) :: ret_dadtau
    real(8) :: o_k_0
    ret_dadtau = dadtau(axp_tau=axp_tau, O_mat_0=o_mat_0, O_vac_0=o_vac_0, O_k_0=o_k_0)
end subroutine f90wrap_dadtau

subroutine f90wrap_dadt(axp_t, o_mat_0, o_vac_0, ret_dadt, o_k_0)
    use cosmology, only: dadt
    implicit none
    
    real(8) :: axp_t
    real(8) :: o_mat_0
    real(8) :: o_vac_0
    real(8), intent(out) :: ret_dadt
    real(8) :: o_k_0
    ret_dadt = dadt(axp_t=axp_t, O_mat_0=o_mat_0, O_vac_0=o_vac_0, O_k_0=o_k_0)
end subroutine f90wrap_dadt

! End of module cosmology defined in file cosmology_module.fpp

