! T04 learning slots. See labs/ch03/task.md. No solution import or automatic fallback.
module ns_polytrope
    use ns_precision, only: wp
    use ns_status
    use ns_eos_types
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
    implicit none
    private
    type, public :: polytrope_t
        private
        logical :: ready=.false.
    end type
    public :: init_polytrope,init_polytrope_k,eos_at_n,eos_at_pressure,eos_at_mu
    public :: polytrope_domain,polytrope_ready
contains
    pure logical function polytrope_ready(model)
        type(polytrope_t),intent(in)::model
        polytrope_ready=model%ready
    end function
    pure subroutine init_polytrope(model,n_anchor,p_anchor,gamma,c_j,n_min,n_max,status)
        type(polytrope_t),intent(out)::model
        real(wp),intent(in)::n_anchor,p_anchor,gamma,c_j,n_min,n_max
        type(status_t),intent(out)::status
        model=polytrope_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','T04: implement polytrope initialization')
    end subroutine
    pure subroutine init_polytrope_k(model,k_si,gamma,c_j,n_anchor,n_min,n_max,status)
        type(polytrope_t),intent(out)::model
        real(wp),intent(in)::k_si,gamma,c_j,n_anchor,n_min,n_max
        type(status_t),intent(out)::status
        model=polytrope_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','T04: implement K constructor')
    end subroutine
    pure subroutine polytrope_domain(model,lower,upper,status)
        type(polytrope_t),intent(in)::model
        type(eos_state_t),intent(out)::lower,upper
        type(status_t),intent(out)::status
        call clear_eos_state(lower);call clear_eos_state(upper)
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','T04: implement declared domain')
    end subroutine
    pure subroutine eos_at_n(model,n_b_m3,state,status)
        type(polytrope_t),intent(in)::model
        real(wp),intent(in)::n_b_m3
        type(eos_state_t),intent(out)::state
        type(status_t),intent(out)::status
        call clear_eos_state(state)
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','T04: implement density evaluation')
    end subroutine
    pure subroutine eos_at_pressure(model,p_pa,side,state,status)
        type(polytrope_t),intent(in)::model
        real(wp),intent(in)::p_pa
        integer,intent(in)::side
        type(eos_state_t),intent(out)::state
        type(status_t),intent(out)::status
        call clear_eos_state(state)
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','T04: implement pressure inverse')
    end subroutine
    pure subroutine eos_at_mu(model,mu_j,state,dn_dmu,status)
        type(polytrope_t),intent(in)::model
        real(wp),intent(in)::mu_j
        type(eos_state_t),intent(out)::state
        real(wp),intent(out)::dn_dmu
        type(status_t),intent(out)::status
        call clear_eos_state(state);dn_dmu=ieee_value(0._wp,ieee_quiet_nan)
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','T04: implement chemical-potential inverse')
    end subroutine
end module
