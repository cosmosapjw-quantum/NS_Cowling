module ns_maxwell
    use ns_precision, only: wp
    use ns_status
    use ns_eos_types
    use ns_polytrope
    use ns_matching_types, only: transition_state_t
    use ns_numerics_types, only: root_options_t,root_result_t
    use ns_bisection, only: bisect
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: solve_maxwell_match
contains
    subroutine solve_maxwell_match(low,high,mu_left_j,mu_right_j,options,transition,status)
        type(polytrope_t),intent(in)::low,high
        real(wp),intent(in)::mu_left_j,mu_right_j
        type(root_options_t),intent(in)::options
        type(transition_state_t),intent(out)::transition
        type(status_t),intent(out)::status
        type(root_result_t)::root
        type(eos_state_t)::lo_state,hi_state,lleft,hleft,lright,hright
        real(wp)::dummy,fl,fr,p_scale,mu_scale,thermo_scale,thermo_resid
        transition=transition_state_t()
        if(.not.all(ieee_is_finite([mu_left_j,mu_right_j]))) then
            call fail_status(status,NS_NONFINITE_INPUT,'physics','nonfinite Maxwell bracket');return
        end if
        if(mu_left_j<=0._wp .or. mu_right_j<=mu_left_j) then
            call fail_status(status,NS_INVALID_INPUT,'physics','Maxwell bracket must satisfy 0<mu_left<mu_right');return
        end if
        call eos_at_mu(low,mu_left_j,lleft,dummy,status);if(.not.is_ok(status))return
        call eos_at_mu(high,mu_left_j,hleft,dummy,status);if(.not.is_ok(status))return
        call eos_at_mu(low,mu_right_j,lright,dummy,status);if(.not.is_ok(status))return
        call eos_at_mu(high,mu_right_j,hright,dummy,status);if(.not.is_ok(status))return
        fl=hleft%pressure_pa-lleft%pressure_pa;fr=hright%pressure_pa-lright%pressure_pa
        call bisect(residual,mu_left_j,mu_right_j,options,root)
        if(.not.is_ok(root%status)) then
            status=root%status;return
        end if
        call eos_at_mu(low,root%root,lo_state,dummy,status);if(.not.is_ok(status))return
        call eos_at_mu(high,root%root,hi_state,dummy,status);if(.not.is_ok(status))return
        transition%mu_t_j=root%root;transition%p_low_pa=lo_state%pressure_pa;transition%p_high_pa=hi_state%pressure_pa
        transition%p_t_pa=0.5_wp*(transition%p_low_pa+transition%p_high_pa)
        transition%pressure_residual_pa=transition%p_high_pa-transition%p_low_pa
        transition%mu_low_j=lo_state%mu_j;transition%mu_high_j=hi_state%mu_j
        transition%mu_mismatch_j=hi_state%mu_j-lo_state%mu_j
        transition%n_low_m3=lo_state%n_b_m3;transition%n_high_m3=hi_state%n_b_m3
        transition%delta_n_m3=hi_state%n_b_m3-lo_state%n_b_m3
        transition%energy_low_j_m3=lo_state%energy_j_m3;transition%energy_high_j_m3=hi_state%energy_j_m3
        transition%delta_e_j_m3=hi_state%energy_j_m3-lo_state%energy_j_m3
        transition%root_iterations=root%iterations;transition%root_evaluations=root%evaluations
        transition%branch_selection_valid=fl<0._wp .and. fr>0._wp .and. transition%delta_n_m3>0._wp
        p_scale=max(abs(transition%p_t_pa),tiny(1._wp));mu_scale=max(abs(transition%mu_t_j),tiny(1._wp))
        thermo_scale=max(abs(transition%delta_e_j_m3),abs(transition%mu_t_j*transition%delta_n_m3),tiny(1._wp))
        thermo_resid=transition%delta_e_j_m3-transition%mu_t_j*transition%delta_n_m3
        transition%thermodynamically_consistent= &
            abs(transition%pressure_residual_pa)<=max(options%ftol,5e-12_wp*p_scale) .and. &
            abs(transition%mu_mismatch_j)<=max(options%xtol,5e-12_wp*mu_scale) .and. &
            abs(thermo_resid)<=5e-12_wp*thermo_scale .and. transition%branch_selection_valid
        if(.not.transition%thermodynamically_consistent) then
            call fail_status(status,NS_MATCHING_INCONSISTENT,'physics', &
                'Maxwell crossing violates phase ordering or thermodynamic closure')
            return
        end if
        transition%valid=.true.;call clear_status(status)
    contains
        subroutine residual(mu,value,st)
            real(wp),intent(in)::mu
            real(wp),intent(out)::value
            type(status_t),intent(out)::st
            type(eos_state_t)::a,b
            real(wp)::d
            call eos_at_mu(low,mu,a,d,st);if(.not.is_ok(st))return
            call eos_at_mu(high,mu,b,d,st);if(.not.is_ok(st))return
            value=b%pressure_pa-a%pressure_pa
            if(.not.ieee_is_finite(value)) then
                call fail_status(st,NS_NONFINITE_RESULT,'numerics','Maxwell pressure residual nonfinite')
            end if
        end subroutine
    end subroutine
end module ns_maxwell
