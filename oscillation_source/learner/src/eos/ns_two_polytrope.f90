module ns_two_polytrope
    use ns_precision,only:wp
    use ns_status
    use ns_polytrope
    use ns_eos_types,only:eos_state_t,clear_eos_state
    implicit none
    private
    type, public :: two_polytrope_t
        logical::ready=.false.
    end type
    public :: init_continuous_two_polytrope,two_polytrope_transition,two_polytrope_ready
contains
    pure logical function two_polytrope_ready(model)
        type(two_polytrope_t),intent(in)::model
        two_polytrope_ready=model%ready
    end function
    subroutine init_continuous_two_polytrope(model,core,gamma_crust,n_crust_min,n_transition,status)
        type(two_polytrope_t),intent(out)::model
        type(polytrope_t),intent(in)::core
        real(wp),intent(in)::gamma_crust,n_crust_min,n_transition
        type(status_t),intent(out)::status
        model=two_polytrope_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation', &
            'student two-polytrope slot not implemented')
    end subroutine
    subroutine two_polytrope_transition(model,crust_state,core_state,status)
        type(two_polytrope_t),intent(in)::model
        type(eos_state_t),intent(out)::crust_state,core_state
        type(status_t),intent(out)::status
        call clear_eos_state(crust_state);call clear_eos_state(core_state)
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation', &
            'student two-polytrope transition slot not implemented')
    end subroutine
end module ns_two_polytrope
