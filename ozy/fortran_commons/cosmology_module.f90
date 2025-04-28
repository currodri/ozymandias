module cosmology
    use local
    use io_ramses
    contains

    !---------------------------------------------------------------
    ! Subroutine: COSMOLOGY MODEL
    !
    ! This subroutines initialises the FLRW cosmological model
    ! using the parameters read from a simulation snapshot.
    ! It creates from there the look-up tables and the age of the
    ! simulation in that snapshot.
    !---------------------------------------------------------------
    subroutine cosmology_model
        implicit none

        ! These values are hardcoded, such that consistency with other
        ! RAMSES utils can be used
        integer :: i
        real(dbl) :: alpha,axp_min

        sim%n_frw = 1000
        alpha = 1.d-6
        axp_min = 1.d-3

        ! Allocate look-up tables
        if(allocated(sim%aexp_frw))deallocate(sim%aexp_frw,sim%hexp_frw,sim%tau_frw,sim%t_frw)
        allocate(sim%aexp_frw(0:sim%n_frw),sim%hexp_frw(0:sim%n_frw))
        allocate(sim%tau_frw(0:sim%n_frw),sim%t_frw(0:sim%n_frw))

        ! Compute Friedmann model look-up tables
        if (verbose) write(*,*)'Computing Friedmann model'
        call friedmann(dble(sim%omega_m),dble(sim%omega_l),dble(sim%omega_k), &
            & alpha,axp_min,sim%aexp_frw,sim%hexp_frw,sim%tau_frw,sim%t_frw,sim%n_frw,sim%time_tot)

        ! Find neighboring expansion factors
        i=1
        do while(sim%aexp_frw(i)>sim%aexp.and.i<sim%n_frw)
            i=i+1
        end do

        ! Interpolate time
        sim%time_simu=sim%t_frw(i)*(sim%aexp-sim%aexp_frw(i-1))/(sim%aexp_frw(i)-sim%aexp_frw(i-1))+ &
            & sim%t_frw(i-1)*(sim%aexp-sim%aexp_frw(i))/(sim%aexp_frw(i-1)-sim%aexp_frw(i))
        if (verbose) write(*,*)'Age simu=',(sim%time_tot+sim%time_simu)/(sim%h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
    end subroutine cosmology_model


    !---------------------------------------------------------------
    ! Subroutine: FRIEDMANN
    !
    ! This subroutine assumes that axp = 1 at z = 0 (today) 
    ! and that t and tau = 0 at z = 0 (today).              
    ! axp is the expansion factor, hexp the Hubble constant 
    ! defined as hexp=1/axp*daxp/dtau, tau the conformal    
    ! time, and t the look-back time, both in unit of 1/H0. 
    ! alpha is the required accuracy and axp_min is the     
    ! starting expansion factor of the look-up table.       
    ! ntable is the required size of the look-up table. 
    !---------------------------------------------------------------
    subroutine friedmann(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
        & axp_out,hexp_out,tau_out,t_out,ntable,age_tot)
        implicit none
        integer,intent(in) :: ntable
        real(dbl),intent(in) :: O_mat_0, O_vac_0, O_k_0
        real(dbl),intent(in) :: alpha,axp_min
        real(dbl),intent(inout) :: age_tot
        real(dbl),dimension(0:ntable),intent(inout) :: axp_out,hexp_out,tau_out,t_out
        
        real(dbl) :: axp_tau, axp_t
        real(dbl) :: axp_tau_pre, axp_t_pre
        real(dbl) :: dtau,dt
        real(dbl) :: tau,t
        integer :: nstep,nout,nskip

        axp_tau = 1.0D0
        axp_t = 1.0D0
        tau = 0.0D0
        t = 0.0D0
        nstep = 0
        
        do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
            nstep = nstep + 1
            dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
            axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
            axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
            tau = tau - dtau
            
            dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
            axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
            axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
            t = t - dt
        end do

        age_tot=-t
        if (verbose) write(*,666)-t
        666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

        nskip=nstep/ntable
        
        axp_t = 1.d0
        t = 0.d0
        axp_tau = 1.d0
        tau = 0.d0
        nstep = 0
        nout=0
        t_out(nout)=t
        tau_out(nout)=tau
        axp_out(nout)=axp_tau
        hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

        do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
            nstep = nstep + 1
            dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
            axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
            axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
            tau = tau - dtau

            dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
            axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
            axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
            t = t - dt
            
            if(mod(nstep,nskip)==0)then
                nout=nout+1
                t_out(nout)=t
                tau_out(nout)=tau
                axp_out(nout)=axp_tau
                hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
            end if
        end do
        t_out(ntable)=t
        tau_out(ntable)=tau
        axp_out(ntable)=axp_tau
        hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

    end subroutine friedmann

    function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0) 
        real(dbl)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
        dadtau = axp_tau*axp_tau*axp_tau *  &
                &   ( O_mat_0 + &
                &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
                &     O_k_0   * axp_tau )
        dadtau = sqrt(dadtau)
        return
    end function dadtau

    function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
        real(dbl)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
        dadt   = (1.0D0/axp_t)* &
                &   ( O_mat_0 + &
                &     O_vac_0 * axp_t*axp_t*axp_t + &
                &     O_k_0   * axp_t )
        dadt = sqrt(dadt)
        return
    end function dadt
end module cosmology