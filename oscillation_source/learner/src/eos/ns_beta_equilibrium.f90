module ns_beta_equilibrium
    use ns_precision,only:wp
    use ns_status
    use ns_nuclear_toy,only:nuclear_model_t,nuclear_state_t
    use ns_electron_gas,only:electron_state_t
    implicit none
    private
    type,public :: beta_state_t
        logical::valid=.false.,charge_neutrality_closed=.false.,causal_local=.false.
        real(wp)::n_b_m3=0._wp,proton_fraction=0._wp,alpha=0._wp,n_p_m3=0._wp,n_e_m3=0._wp
        real(wp)::beta_residual_mev=0._wp,energy_j_m3=0._wp,pressure_pa=0._wp,chi_equilibrium=0._wp
        type(nuclear_state_t)::nuclear
        type(electron_state_t)::electron
        integer::iterations=0
    end type
    public :: solve_beta_equilibrium
contains
    subroutine solve_beta_equilibrium(model,n,state,status)
        type(nuclear_model_t),intent(in)::model;real(wp),intent(in)::n
        type(beta_state_t),intent(out)::state;type(status_t),intent(out)::status
        state=beta_state_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','student beta-equilibrium slot not implemented')
    end subroutine
end module ns_beta_equilibrium
