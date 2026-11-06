module ns_nuclear_toy
    use ns_precision,only:wp
    use ns_status
    implicit none
    private
    type,public :: nuclear_model_t
        logical :: ready=.false.,regulated=.false.
        real(wp)::n0_fm3=0._wp,binding_mev=0._wp,k0_mev=0._wp,s0_mev=0._wp,regulator_c=0._wp
        real(wp)::m_nucleon_kg=0._wp,t0_mev=0._wp,a_mev=0._wp,b_mev=0._wp,sigma=0._wp
    end type
    type,public :: nuclear_state_t
        logical :: valid=.false.
        real(wp)::n_b_m3=0._wp,n_b_fm3=0._wp,u=0._wp,alpha=0._wp,proton_fraction=0._wp
        real(wp)::symmetric_energy_mev=0._wp,symmetry_mev=0._wp,energy_excess_mev=0._wp
        real(wp)::exact_kinetic_mev=0._wp,exact_kinetic_excess_mev=0._wp
        real(wp)::pressure_nuclear_pa=0._wp,energy_nuclear_j_m3=0._wp,chi_nuclear=0._wp
    end type
    public :: init_nuclear_model,evaluate_nuclear,symmetry_energy_mev,exact_kinetic_per_baryon
    public :: scan_nuclear_causality
contains
    subroutine init_nuclear_model(model,regulated,regulator_c,n0_fm3,binding_mev,k0_mev,s0_mev,status)
        type(nuclear_model_t),intent(out)::model
        logical,intent(in)::regulated
        real(wp),intent(in)::regulator_c,n0_fm3,binding_mev,k0_mev,s0_mev
        type(status_t),intent(out)::status
        model=nuclear_model_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','student nuclear-model slot not implemented')
    end subroutine
    subroutine evaluate_nuclear(model,n_b_m3,alpha,state,status)
        type(nuclear_model_t),intent(in)::model;real(wp),intent(in)::n_b_m3,alpha
        type(nuclear_state_t),intent(out)::state;type(status_t),intent(out)::status
        state=nuclear_state_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','student nuclear-state slot not implemented')
    end subroutine
    subroutine symmetry_energy_mev(model,u,s,status)
        type(nuclear_model_t),intent(in)::model;real(wp),intent(in)::u;real(wp),intent(out)::s
        type(status_t),intent(out)::status
        s=0._wp;call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','student symmetry-energy slot not implemented')
    end subroutine
    subroutine exact_kinetic_per_baryon(model,n_b_m3,alpha,t,status)
        type(nuclear_model_t),intent(in)::model;real(wp),intent(in)::n_b_m3,alpha;real(wp),intent(out)::t
        type(status_t),intent(out)::status
        t=0._wp;call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','student exact-kinetic slot not implemented')
    end subroutine
    subroutine scan_nuclear_causality(model,u_max,max_chi,causal,status)
        type(nuclear_model_t),intent(in)::model;real(wp),intent(in)::u_max;real(wp),intent(out)::max_chi
        logical,intent(out)::causal;type(status_t),intent(out)::status
        max_chi=0._wp;causal=.false.
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','student nuclear-causality slot not implemented')
    end subroutine
end module ns_nuclear_toy
