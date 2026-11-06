! Cold single polytrope, SI: epsilon includes the declared C*n rest/linear energy.
! Analytic inverses; no bisection, interpolation, beta-equilibrium or global EOS state.
module ns_polytrope
    use ns_precision, only: wp
    use ns_status
    use ns_safe_math, only: checked_scale
    use ns_eos_types
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    private
    type, public :: polytrope_t
        private
        real(wp) :: n_anchor=0._wp,p_anchor=0._wp,gamma=0._wp,c_j=0._wp
        real(wp) :: n_min=0._wp,n_max=0._wp,thermal_anchor=0._wp
        type(eos_state_t) :: lower,upper
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
        type(polytrope_t)::candidate
        type(eos_state_t)::anchor
        model=polytrope_t()
        if (.not.all(ieee_is_finite([n_anchor,p_anchor,gamma,c_j,n_min,n_max]))) then
            call fail_status(status,NS_NONFINITE_INPUT,'physics','nonfinite polytrope parameter');return
        end if
        if (gamma<=1._wp .or. c_j<=0._wp .or. p_anchor<=0._wp .or. n_min<=0._wp &
            .or. n_max<=n_min .or. n_anchor<n_min .or. n_anchor>n_max) then
            call fail_status(status,NS_INVALID_INPUT,'physics','require Gamma>1, C>0, P_anchor>0 and 0<n_min<=n_anchor<=n_max')
            return
        end if
        candidate%n_anchor=n_anchor;candidate%p_anchor=p_anchor;candidate%gamma=gamma
        candidate%c_j=c_j;candidate%n_min=n_min;candidate%n_max=n_max
        call evaluate_density(candidate,n_min,candidate%lower,status)
        if (.not.is_ok(status)) return
        call evaluate_density(candidate,n_max,candidate%upper,status)
        if (.not.is_ok(status)) return
        call evaluate_density(candidate,n_anchor,anchor,status)
        if (.not.is_ok(status)) return
        ! With C>0, Gamma>1, chi(n) is increasing; checking the upper end bounds the interval.
        if (candidate%upper%chi>1._wp) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','specified interval violates causal chi<=1');return
        end if
        if (candidate%lower%mu_j<=c_j .or. candidate%upper%mu_j<=candidate%lower%mu_j &
            .or. candidate%upper%pressure_pa<=candidate%lower%pressure_pa) then
            call fail_status(status,NS_STAGNATION,'numerics','mu or pressure interval cannot be resolved in real64');return
        end if
        call positive_product(anchor%pressure_pa/anchor%n_b_m3,gamma/(gamma-1._wp), &
                              candidate%thermal_anchor,status)
        if (.not.is_ok(status)) return
        candidate%ready=.true.;model=candidate
        call clear_status(status)
    end subroutine

    pure subroutine init_polytrope_k(model,k_si,gamma,c_j,n_anchor,n_min,n_max,status)
        type(polytrope_t),intent(out)::model
        real(wp),intent(in)::k_si,gamma,c_j,n_anchor,n_min,n_max
        type(status_t),intent(out)::status
        real(wp)::p,lg
        model=polytrope_t()
        if (.not.all(ieee_is_finite([k_si,gamma,c_j,n_anchor,n_min,n_max]))) then
            call fail_status(status,NS_NONFINITE_INPUT,'physics','nonfinite K constructor input');return
        end if
        if(k_si<=0._wp .or. gamma<=1._wp .or. n_anchor<=0._wp) then
            call fail_status(status,NS_INVALID_INPUT,'physics','K, anchor and Gamma must be positive, Gamma>1');return
        end if
        lg=log(n_anchor)
        if(abs(lg)>1._wp) then
            if(gamma>huge(lg)/abs(lg)) then
                call fail_status(status,NS_RANGE_ERROR,'numerics','K constructor exponent overflow');return
            end if
        end if
        call positive_exp(log(k_si)+gamma*lg,p,status)
        if(.not.is_ok(status)) return
        call init_polytrope(model,n_anchor,p,gamma,c_j,n_min,n_max,status)
    end subroutine

    pure subroutine polytrope_domain(model,lower,upper,status)
        type(polytrope_t),intent(in)::model
        type(eos_state_t),intent(out)::lower,upper
        type(status_t),intent(out)::status
        call clear_eos_state(lower);call clear_eos_state(upper)
        if(.not.model%ready) then
            call fail_status(status,NS_UNINITIALIZED,'implementation','polytrope was not successfully initialized');return
        end if
        lower=model%lower;upper=model%upper;call clear_status(status)
    end subroutine

    pure subroutine eos_at_n(model,n_b_m3,state,status)
        type(polytrope_t),intent(in)::model
        real(wp),intent(in)::n_b_m3
        type(eos_state_t),intent(out)::state
        type(status_t),intent(out)::status
        call clear_eos_state(state)
        call check_query(model,n_b_m3,model%n_min,model%n_max,status)
        if(.not.is_ok(status)) return
        call evaluate_density(model,n_b_m3,state,status)
    end subroutine

    pure subroutine eos_at_pressure(model,p_pa,side,state,status)
        type(polytrope_t),intent(in)::model
        real(wp),intent(in)::p_pa
        integer,intent(in)::side
        type(eos_state_t),intent(out)::state
        type(status_t),intent(out)::status
        real(wp)::n
        call clear_eos_state(state)
        if(side/=SIDE_REGULAR .and. side/=SIDE_INNER .and. side/=SIDE_OUTER) then
            call fail_status(status,NS_INVALID_INPUT,'physics','unknown EOS side label');return
        end if
        call check_query(model,p_pa,model%lower%pressure_pa,model%upper%pressure_pa,status)
        if(.not.is_ok(status)) return
        if(p_pa==model%lower%pressure_pa) then
            state=model%lower
        else if(p_pa==model%upper%pressure_pa) then
            state=model%upper
        else
            call scaled_power(model%n_anchor,p_pa,model%p_anchor,1._wp/model%gamma,n,status)
            if(.not.is_ok(status)) return
            call inverse_roundoff_boundary(model,n,status)
            if(.not.is_ok(status)) return
            call evaluate_density(model,n,state,status)
            if(.not.is_ok(status)) return
        end if
        state%side=side
    end subroutine

    pure subroutine eos_at_mu(model,mu_j,state,dn_dmu,status)
        type(polytrope_t),intent(in)::model
        real(wp),intent(in)::mu_j
        type(eos_state_t),intent(out)::state
        real(wp),intent(out)::dn_dmu
        type(status_t),intent(out)::status
        real(wp)::n
        call clear_eos_state(state);dn_dmu=ieee_value(0._wp,ieee_quiet_nan)
        call check_query(model,mu_j,model%lower%mu_j,model%upper%mu_j,status)
        if(.not.is_ok(status)) return
        if(mu_j==model%lower%mu_j) then
            state=model%lower
        else if(mu_j==model%upper%mu_j) then
            state=model%upper
        else
            call scaled_power(model%n_anchor,mu_j-model%c_j,model%thermal_anchor, &
                              1._wp/(model%gamma-1._wp),n,status)
            if(.not.is_ok(status)) return
            call inverse_roundoff_boundary(model,n,status)
            if(.not.is_ok(status)) return
            call evaluate_density(model,n,state,status)
            if(.not.is_ok(status)) return
        end if
        ! dn/dmu=n/(mu*chi); log form avoids n*n or a tiny intermediate product.
        call positive_exp(log(state%n_b_m3)-log(state%mu_j)-log(state%chi),dn_dmu,status)
        if(.not.is_ok(status)) call clear_eos_state(state)
    end subroutine

    pure subroutine check_query(model,value,lo,hi,status)
        type(polytrope_t),intent(in)::model
        real(wp),intent(in)::value,lo,hi
        type(status_t),intent(out)::status
        if(.not.model%ready) then
            call fail_status(status,NS_UNINITIALIZED,'implementation','polytrope was not successfully initialized')
        else if(.not.ieee_is_finite(value)) then
            call fail_status(status,NS_NONFINITE_INPUT,'physics','EOS query is not finite')
        else if(value<lo .or. value>hi) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','query outside the declared closed positive-density interval')
        else
            call clear_status(status)
        end if
    end subroutine

    pure subroutine evaluate_density(model,n,state,status)
        type(polytrope_t),intent(in)::model
        real(wp),intent(in)::n
        type(eos_state_t),intent(out)::state
        type(status_t),intent(out)::status
        type(eos_state_t)::candidate
        real(wp)::p,rest,internal,thermal,mu,energy,q,chi
        call clear_eos_state(state)
        call scaled_power(model%p_anchor,n,model%n_anchor,model%gamma,p,status)
        if(.not.is_ok(status)) return
        call positive_product(n,model%c_j,rest,status)
        if(.not.is_ok(status)) return
        call positive_product(p,1._wp/(model%gamma-1._wp),internal,status)
        if(.not.is_ok(status)) return
        call positive_sum(rest,internal,energy,status)
        if(.not.is_ok(status)) return
        q=p/n
        if(.not.ieee_is_finite(q) .or. q<=0._wp) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','P/n not representable');return
        end if
        call positive_product(q,model%gamma/(model%gamma-1._wp),thermal,status)
        if(.not.is_ok(status)) return
        call positive_sum(model%c_j,thermal,mu,status)
        if(.not.is_ok(status)) return
        chi=model%gamma*(q/mu)
        if(.not.ieee_is_finite(chi) .or. chi<=0._wp) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','positive chi not representable');return
        end if
        candidate=eos_state_t(p,energy,n,mu,chi,SIDE_REGULAR,.true.)
        state=candidate;call clear_status(status)
    end subroutine

    pure subroutine inverse_roundoff_boundary(model,n,status)
        type(polytrope_t),intent(in)::model
        real(wp),intent(inout)::n
        type(status_t),intent(out)::status
        real(wp),parameter::tol=32._wp*epsilon(1._wp)
        ! Input P/mu was already INSIDE its strict interval. Only inverse arithmetic rounding is snapped.
        if(n<model%n_min) then
            if((model%n_min-n)/model%n_min<=tol) then
                n=model%n_min
            else
                call fail_status(status,NS_STAGNATION,'numerics','inverse density below bound beyond rounding');return
            end if
        else if(n>model%n_max) then
            if((n-model%n_max)/model%n_max<=tol) then
                n=model%n_max
            else
                call fail_status(status,NS_STAGNATION,'numerics','inverse density above bound beyond rounding');return
            end if
        end if
        call clear_status(status)
    end subroutine

    pure subroutine positive_sum(a,b,value,status)
        real(wp),intent(in)::a,b
        real(wp),intent(out)::value
        type(status_t),intent(out)::status
        value=ieee_value(0._wp,ieee_quiet_nan)
        if(a>huge(a)-b) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','EOS sum overflows');return
        end if
        value=a+b;call clear_status(status)
    end subroutine

    pure subroutine positive_product(a,b,value,status)
        real(wp),intent(in)::a,b
        real(wp),intent(out)::value
        type(status_t),intent(out)::status
        call checked_scale(a,b,value,status)
        if(is_ok(status)) then
            if(value<=0._wp) call fail_status(status,NS_RANGE_ERROR,'numerics','positive EOS product vanished')
        end if
    end subroutine

    pure subroutine positive_exp(exponent,value,status)
        real(wp),intent(in)::exponent
        real(wp),intent(out)::value
        type(status_t),intent(out)::status
        value=ieee_value(0._wp,ieee_quiet_nan)
        if(.not.ieee_is_finite(exponent)) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','EOS logarithmic magnitude not finite');return
        end if
        if(exponent>log(huge(value)) .or. exponent<log(nearest(0._wp,1._wp))) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','EOS value outside real64 magnitude range');return
        end if
        value=exp(exponent)
        if(.not.ieee_is_finite(value) .or. value<=0._wp) then
            value=ieee_value(0._wp,ieee_quiet_nan)
            call fail_status(status,NS_RANGE_ERROR,'numerics','EOS exponential not representable');return
        end if
        call clear_status(status)
    end subroutine

    pure subroutine scaled_power(base,x,anchor,power,value,status)
        real(wp),intent(in)::base,x,anchor,power
        real(wp),intent(out)::value
        type(status_t),intent(out)::status
        real(wp)::lr,t,factor,logvalue
        value=ieee_value(0._wp,ieee_quiet_nan)
        if(x==anchor) then
            value=base;call clear_status(status);return
        end if
        lr=log(x)-log(anchor)
        if(abs(lr)<500._wp) lr=log(x/anchor)
        if(abs(lr)>1._wp) then
            if(abs(power)>huge(power)/abs(lr)) then
                call fail_status(status,NS_RANGE_ERROR,'numerics','power exponent overflow');return
            end if
        end if
        t=power*lr;logvalue=log(base)+t
        if(abs(t)<log(huge(t))-2._wp) then
            factor=exp(t)
            call positive_product(base,factor,value,status)
        else
            call positive_exp(logvalue,value,status)
        end if
    end subroutine
end module ns_polytrope
