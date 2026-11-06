module ns_mu_bridge
    use ns_precision, only: wp
    use ns_status
    use ns_eos_types
    use ns_polytrope
    use ns_matching_types, only: mu_bridge_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: construct_branch_bridge,construct_bridge_endpoints,evaluate_mu_bridge
contains
    subroutine construct_branch_bridge(low,high,mu_low_j,mu_high_j,bridge,status)
        type(polytrope_t),intent(in)::low,high
        real(wp),intent(in)::mu_low_j,mu_high_j
        type(mu_bridge_t),intent(out)::bridge
        type(status_t),intent(out)::status
        type(eos_state_t)::a,b
        real(wp)::dummy
        bridge=mu_bridge_t()
        if(.not.all(ieee_is_finite([mu_low_j,mu_high_j]))) then
            call fail_status(status,NS_NONFINITE_INPUT,'physics','nonfinite bridge interval');return
        end if
        if(mu_low_j<=0._wp .or. mu_high_j<=mu_low_j) then
            call fail_status(status,NS_INVALID_INPUT,'physics','bridge requires 0<mu_low<mu_high');return
        end if
        call eos_at_mu(low,mu_low_j,a,dummy,status);if(.not.is_ok(status))return
        call eos_at_mu(high,mu_high_j,b,dummy,status);if(.not.is_ok(status))return
        call construct_bridge_endpoints(mu_low_j,a%pressure_pa,a%n_b_m3,mu_high_j,b%pressure_pa,b%n_b_m3,bridge,status)
    end subroutine

    subroutine construct_bridge_endpoints(mu0,p0,n0,mu1,p1,n1,bridge,status)
        real(wp),intent(in)::mu0,p0,n0,mu1,p1,n1
        type(mu_bridge_t),intent(out)::bridge
        type(status_t),intent(out)::status
        real(wp)::h,d20,d21,g0,g1
        bridge=mu_bridge_t()
        if(.not.all(ieee_is_finite([mu0,p0,n0,mu1,p1,n1]))) then
            call fail_status(status,NS_NONFINITE_INPUT,'physics','nonfinite bridge endpoint');return
        end if
        if(mu0<=0._wp .or. mu1<=mu0 .or. min(p0,p1,n0,n1)<=0._wp) then
            call fail_status(status,NS_INVALID_INPUT,'physics','bridge endpoints require positive ordered mu,P,n');return
        end if
        h=mu1-mu0
        bridge%mu_low_j=mu0;bridge%mu_high_j=mu1;bridge%p_low_pa=p0;bridge%p_high_pa=p1
        bridge%n_low_m3=n0;bridge%n_high_m3=n1
        bridge%a0=p0;bridge%a1=h*n0
        bridge%a2=3._wp*(p1-p0)-h*(2._wp*n0+n1)
        bridge%a3=2._wp*(p0-p1)+h*(n0+n1)
        if(.not.all(ieee_is_finite([bridge%a0,bridge%a1,bridge%a2,bridge%a3]))) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','bridge coefficients not representable');return
        end if
        d20=2._wp*bridge%a2/(h*h);d21=(2._wp*bridge%a2+6._wp*bridge%a3)/(h*h)
        g0=mu0*d20-n0;g1=mu1*d21-n1
        bridge%min_d2p=min(d20,d21);bridge%min_causal_margin=min(g0,g1)
        bridge%interval_extrema_checked=.true.
        ! P'' is linear. For mu>0, d(mu P''-P')/dmu = mu P''' has fixed sign.
        ! Endpoint extrema therefore exactly determine stability and chi<=1 on the cubic interval.
        bridge%stable_interval=bridge%min_d2p>0._wp
        bridge%causal_interval=bridge%stable_interval .and. bridge%min_causal_margin>=0._wp
        if(.not.bridge%stable_interval .or. .not.bridge%causal_interval) then
            call fail_status(status,NS_MATCHING_INCONSISTENT,'physics', &
                'cubic P(mu) bridge violates stability or causality on interval')
            return
        end if
        bridge%ready=.true.;call clear_status(status)
    end subroutine

    subroutine evaluate_mu_bridge(bridge,mu,state,dn_dmu,status)
        type(mu_bridge_t),intent(in)::bridge
        real(wp),intent(in)::mu
        type(eos_state_t),intent(out)::state
        real(wp),intent(out)::dn_dmu
        type(status_t),intent(out)::status
        real(wp)::h,t,p,n,eps,chi
        call clear_eos_state(state);dn_dmu=0._wp
        if(.not.bridge%ready) then
            call fail_status(status,NS_UNINITIALIZED,'implementation','bridge is not admitted');return
        end if
        if(.not.ieee_is_finite(mu)) then
            call fail_status(status,NS_NONFINITE_INPUT,'physics','nonfinite bridge query');return
        end if
        if(mu<bridge%mu_low_j .or. mu>bridge%mu_high_j) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','bridge query outside closed mu interval');return
        end if
        h=bridge%mu_high_j-bridge%mu_low_j;t=(mu-bridge%mu_low_j)/h
        p=bridge%a0+bridge%a1*t+bridge%a2*t*t+bridge%a3*t*t*t
        n=(bridge%a1+2._wp*bridge%a2*t+3._wp*bridge%a3*t*t)/h
        dn_dmu=(2._wp*bridge%a2+6._wp*bridge%a3*t)/(h*h)
        if(.not.all(ieee_is_finite([p,n,dn_dmu])) .or. min(p,n,dn_dmu)<=0._wp) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','bridge evaluation is nonpositive or nonfinite');return
        end if
        eps=mu*n-p;chi=n/(mu*dn_dmu)
        if(.not.all(ieee_is_finite([eps,chi])) .or. eps<=0._wp .or. chi<=0._wp .or. chi>1._wp) then
            call fail_status(status,NS_MATCHING_INCONSISTENT,'physics', &
                'bridge evaluation violates positive causal thermodynamics')
            return
        end if
        state=eos_state_t(p,eps,n,mu,chi,SIDE_REGULAR,.true.);call clear_status(status)
    end subroutine
end module ns_mu_bridge
