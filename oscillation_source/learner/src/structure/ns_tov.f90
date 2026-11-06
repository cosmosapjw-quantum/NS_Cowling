! SI equilibrium structure. All state variables carry SI units; q is dimensionless.
module ns_tov
    use ns_precision, only: wp
    use ns_status
    use ns_physical_constants, only: g_si,c_si,pi
    use ns_numerics_types, only: ode_options_t,ode_result_t
    use ns_eos_types
    use ns_background_records
    use ns_surface_events
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
    implicit none
    private
    abstract interface
        subroutine eos_pressure_callback(p_pa,side,state,status)
            import wp,eos_state_t,status_t
            real(wp),intent(in)::p_pa
            integer,intent(in)::side
            type(eos_state_t),intent(out)::state
            type(status_t),intent(out)::status
        end subroutine
    end interface
    public :: eos_pressure_callback,tov_rhs,center_start,integrate_tov
contains
    subroutine tov_rhs(eos,r,y,dy,status)
        procedure(eos_pressure_callback)::eos
        real(wp),intent(in)::r,y(:)
        real(wp),intent(out)::dy(:)
        type(status_t),intent(out)::status
        dy=ieee_value(0._wp,ieee_quiet_nan)
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','learner SI TOV RHS')
    end subroutine
    subroutine center_start(eos,pc,r0,nu0,y,status)
        procedure(eos_pressure_callback)::eos
        real(wp),intent(in)::pc,r0,nu0
        real(wp),intent(out)::y(3)
        type(status_t),intent(out)::status
        y=ieee_value(0._wp,ieee_quiet_nan)
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','learner center series')
    end subroutine
    subroutine integrate_tov(eos,pc,options,profile,eos_events)
        procedure(eos_pressure_callback)::eos
        real(wp),intent(in)::pc
        type(star_options_t),intent(in)::options
        type(star_profile_t),intent(out)::profile
        real(wp),optional,intent(in)::eos_events(:)
        profile=star_profile_t()
        call fail_status(profile%status,NS_NOT_IMPLEMENTED,'implementation','learner TOV integration')
    end subroutine
end module
