! Continuous thermodynamic two-polytrope teaching construction with fixed core.
module ns_two_polytrope
    use ns_precision, only: wp
    use ns_status
    use ns_eos_types
    use ns_polytrope
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    type, public :: two_polytrope_t
        private
        type(polytrope_t) :: crust, core
        real(wp) :: n_transition_m3=0._wp, p_transition_pa=0._wp
        logical :: ready=.false.
    end type
    public :: init_continuous_two_polytrope, two_polytrope_transition, two_polytrope_ready
contains
    pure logical function two_polytrope_ready(model)
        type(two_polytrope_t),intent(in)::model
        two_polytrope_ready=model%ready
    end function

    pure subroutine init_continuous_two_polytrope(model,core,gamma_crust,n_crust_min,n_transition,status)
        type(two_polytrope_t),intent(out)::model
        type(polytrope_t),intent(in)::core
        real(wp),intent(in)::gamma_crust,n_crust_min,n_transition
        type(status_t),intent(out)::status
        type(eos_state_t)::s_core,s_crust
        real(wp)::c_crust
        model=two_polytrope_t()
        if(.not.polytrope_ready(core)) then
            call fail_status(status,NS_UNINITIALIZED,'implementation','core branch is not initialized');return
        end if
        if(.not.all(ieee_is_finite([gamma_crust,n_crust_min,n_transition]))) then
            call fail_status(status,NS_NONFINITE_INPUT,'physics','nonfinite two-polytrope input');return
        end if
        if(gamma_crust<=1._wp .or. n_crust_min<=0._wp .or. n_transition<=n_crust_min) then
            call fail_status(status,NS_INVALID_INPUT,'physics','require Gamma_crust>1 and 0<n_crust_min<n_transition');return
        end if
        call eos_at_n(core,n_transition,s_core,status)
        if(.not.is_ok(status)) return
        c_crust=(s_core%energy_j_m3-s_core%pressure_pa/(gamma_crust-1._wp))/n_transition
        if(.not.ieee_is_finite(c_crust) .or. c_crust<=0._wp) then
            call fail_status(status,NS_MATCHING_INCONSISTENT,'physics', &
                'continuous energy completion requires positive crust C')
            return
        end if
        call init_polytrope(model%crust,n_transition,s_core%pressure_pa,gamma_crust,c_crust, &
                            n_crust_min,n_transition,status)
        if(.not.is_ok(status)) return
        call eos_at_n(model%crust,n_transition,s_crust,status)
        if(.not.is_ok(status)) return
        if(s_crust%pressure_pa/=s_core%pressure_pa .or. s_crust%energy_j_m3/=s_core%energy_j_m3) then
            call fail_status(status,NS_MATCHING_INCONSISTENT,'physics', &
                'continuous knot failed exact P/epsilon closure')
            return
        end if
        model%core=core;model%n_transition_m3=n_transition;model%p_transition_pa=s_core%pressure_pa
        model%ready=.true.;call clear_status(status)
    end subroutine

    pure subroutine two_polytrope_transition(model,crust_state,core_state,status)
        type(two_polytrope_t),intent(in)::model
        type(eos_state_t),intent(out)::crust_state,core_state
        type(status_t),intent(out)::status
        call clear_eos_state(crust_state);call clear_eos_state(core_state)
        if(.not.model%ready) then
            call fail_status(status,NS_UNINITIALIZED,'implementation','two-polytrope is not initialized');return
        end if
        call eos_at_n(model%crust,model%n_transition_m3,crust_state,status);if(.not.is_ok(status))return
        call eos_at_n(model%core,model%n_transition_m3,core_state,status)
    end subroutine
end module ns_two_polytrope
