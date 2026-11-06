! SI equilibrium structure. All state variables carry SI units; q is dimensionless.
module ns_tov
    use ns_precision, only: wp
    use ns_status
    use ns_physical_constants, only: g_si,c_si,pi
    use ns_numerics_types, only: ode_options_t,ode_result_t
    use ns_adaptive_rk, only: adaptive_integrate
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
        real(wp),intent(in)::r,y(:) ! m[kg], P[Pa], nu[1]
        real(wp),intent(out)::dy(:)
        type(status_t),intent(out)::status
        type(eos_state_t)::s
        real(wp)::compact,phi
        dy=ieee_value(0._wp,ieee_quiet_nan)
        call fail_status(status,NS_INVALID_INPUT,'physics','TOV requires r>0,m>=0,P>0 and a three-component finite state')
        if(size(y)/=3.or.size(dy)/=3) return
        if(.not.ieee_is_finite(r).or..not.all(ieee_is_finite(y))) return
        if(r<=0._wp.or.y(1)<0._wp.or.y(2)<=0._wp) return
        compact=2._wp*(g_si/c_si**2)*(y(1)/r)
        if(.not.ieee_is_finite(compact).or.compact>=1._wp) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','TOV trapped-surface denominator is not positive');return
        end if
        call eos(y(2),SIDE_REGULAR,s,status)
        if(.not.is_ok(status)) return
        if(.not.s%valid.or..not.all(ieee_is_finite([s%energy_j_m3,s%pressure_pa]))) then
            call fail_status(status,NS_INVALID_INPUT,'physics','invalid EOS callback state');return
        end if
        if(s%energy_j_m3<=0._wp) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','positive energy density required');return
        end if
        if(abs(s%pressure_pa/y(2)-1._wp)>256._wp*epsilon(r)) then
            call fail_status(status,NS_INVALID_INPUT,'physics','EOS callback returned the wrong requested pressure');return
        end if
        phi=(g_si/c_si**2)*(y(1)/r**2+4._wp*pi*r*(y(2)/c_si**2))/(1._wp-compact)
        dy(1)=4._wp*pi*r**2*(s%energy_j_m3/c_si**2)
        dy(2)=-(s%energy_j_m3+y(2))*phi
        dy(3)=phi
        if(.not.all(ieee_is_finite(dy)).or.phi<0._wp) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','nonfinite TOV derivative');return
        end if
        call clear_status(status)
    end subroutine

    subroutine center_start(eos,pc,r0,nu0,y,status)
        procedure(eos_pressure_callback)::eos
        real(wp),intent(in)::pc,r0,nu0
        real(wp),intent(out)::y(3)
        type(status_t),intent(out)::status
        type(eos_state_t)::s,starter
        real(wp)::alpha,drop
        y=ieee_value(0._wp,ieee_quiet_nan)
        call fail_status(status,NS_INVALID_INPUT,'physics','finite Pc>0,r0>0,nu0 required')
        if(.not.all(ieee_is_finite([pc,r0,nu0]))) return
        if(pc<=0._wp.or.r0<=0._wp) return
        call eos(pc,SIDE_REGULAR,s,status)
        if(.not.is_ok(status)) return
        if(.not.s%valid.or.s%energy_j_m3<=0._wp.or..not.ieee_is_finite(s%energy_j_m3)) then
            call fail_status(status,NS_INVALID_INPUT,'physics','invalid central EOS state');return
        end if
        alpha=(2._wp*pi*g_si/3._wp)*((s%energy_j_m3+pc)/c_si**2)* &
              ((s%energy_j_m3+3._wp*pc)/c_si**2)
        drop=alpha*r0**2
        if(.not.ieee_is_finite(drop).or.drop<=0._wp.or.drop>=0.01_wp*pc) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics', &
            'center expansion requires resolved pressure drop <1 percent');return
        end if
        y=[(4._wp*pi/3._wp)*(s%energy_j_m3/c_si**2)*r0**3,pc-drop,nu0]
        if(y(2)==pc.or..not.all(ieee_is_finite(y))) then
            call fail_status(status,NS_STAGNATION,'numerics','center correction not representable');return
        end if
        call eos(y(2),SIDE_REGULAR,starter,status)
    end subroutine

    subroutine integrate_tov(eos,pc,options,profile,eos_events)
        procedure(eos_pressure_callback)::eos
        real(wp),intent(in)::pc
        type(star_options_t),intent(in)::options
        type(star_profile_t),intent(out)::profile
        real(wp),optional,intent(in)::eos_events(:) ! smooth pressure stops; no sharp phase jumps yet
        type(eos_state_t)::s
        type(status_t)::st
        type(ode_options_t)::opt
        type(ode_result_t)::traj
        real(wp)::y0(3),p0,qe,shift,dy(3),p,compact,frac
        real(wp),allocatable::stops(:)
        integer::i,n,astat
        profile=star_profile_t()
        profile%pc_pa=pc;profile%r0_m=options%r0_m;profile%p_cut_pa=options%p_cut_pa
        profile%eos_id=options%eos_id;profile%surface_policy=options%surface_policy
        call fail_status(profile%status,NS_INVALID_INPUT,'physics','invalid surface, radius limit or missing EOS identity')
        if(options%surface_policy/='finite_cutoff_schwarzschild_time') return
        if(len_trim(options%eos_id)==0.or.options%eos_id=='unspecified') return
        if(.not.all(ieee_is_finite([options%r_limit_m,options%p_cut_pa]))) return
        if(options%r_limit_m<=options%r0_m.or.options%p_cut_pa<=0._wp) return
        call center_start(eos,pc,options%r0_m,options%nu_start,y0,profile%status)
        if(.not.is_ok(profile%status)) return
        p0=y0(2)
        call pressure_span(p0,options%p_cut_pa,qe,profile%status)
        if(.not.is_ok(profile%status)) return
        call eos(options%p_cut_pa,SIDE_REGULAR,s,profile%status)
        if(.not.is_ok(profile%status)) return
        ! Validate EOS domain BEFORE integration; no extrapolated endpoint.
        if(present(eos_events)) then
            call pressure_stop_coordinates(eos_events,p0,options%p_cut_pa,stops,profile%status)
            if(.not.is_ok(profile%status)) return
        else
            allocate(stops(0))
        end if
        profile%p_start_pa=p0
        opt=options%ode
        if(.not.ieee_is_finite(opt%h_initial).or.opt%h_initial<=0._wp) then
            call fail_status(profile%status,NS_INVALID_INPUT,'numerics','initial dq must be positive finite');return
        end if
        ! dr/dq is steep at the center. Resolve the declared series before growing dq.
        frac=(pc-p0)/pc
        opt%h_initial=min(opt%h_initial,0.1_wp*frac)
        call adaptive_integrate(q_rhs,0._wp,qe,[options%r0_m,y0(1),y0(3)],opt,traj,stops)
        profile%accepted=traj%accepted;profile%rejected=traj%rejected
        profile%evaluations=traj%evaluations;profile%stops_hit=traj%stops_hit
        profile%status=traj%status
        if(.not.traj%valid.or..not.is_ok(traj%status)) return
        n=size(traj%x)
        allocate(profile%q(n),profile%r_m(n),profile%m_kg(n),profile%p_pa(n), &
                 profile%energy_j_m3(n),profile%n_b_m3(n),profile%mu_j(n),profile%chi(n), &
                 profile%nu(n),profile%lambda(n),profile%dnu_dr(n),profile%segment(n),profile%side(n),stat=astat)
        if(astat/=0) then
            call fail_status(profile%status,NS_RANGE_ERROR,'environment','background allocation failure');return
        end if
        profile%q=traj%x;profile%r_m=traj%y(1,:);profile%m_kg=traj%y(2,:)
        compact=2._wp*(g_si/c_si**2)*(profile%m_kg(n)/profile%r_m(n))
        if(compact<=0._wp.or.compact>=1._wp) then
            call fail_status(profile%status,NS_OUT_OF_DOMAIN,'physics','invalid final compactness');return
        end if
        shift=0.5_wp*log(1._wp-compact)-traj%y(3,n)
        profile%nu=traj%y(3,:)+shift;profile%lapse_shift=shift
        do i=1,n
            call pressure_at_q(traj%x(i),p0,options%p_cut_pa,qe,p,st)
            if(.not.is_ok(st)) then
                profile%status=st;return
            end if
            call eos(p,SIDE_REGULAR,s,st)
            if(.not.is_ok(st)) then
                profile%status=st;return
            end if
            if(.not.s%valid.or..not.all(ieee_is_finite([s%energy_j_m3,s%n_b_m3,s%mu_j,s%chi]))) then
                call fail_status(profile%status,NS_INVALID_INPUT,'physics','invalid background EOS fields');return
            end if
            call tov_rhs(eos,profile%r_m(i),[profile%m_kg(i),p,profile%nu(i)],dy,st)
            if(.not.is_ok(st)) then
                profile%status=st;return
            end if
            profile%p_pa(i)=p;profile%energy_j_m3(i)=s%energy_j_m3
            profile%n_b_m3(i)=s%n_b_m3;profile%mu_j(i)=s%mu_j;profile%chi(i)=s%chi
            profile%dnu_dr(i)=dy(3)
            profile%lambda(i)=-0.5_wp*log(1._wp-2._wp*(g_si/c_si**2)*(profile%m_kg(i)/profile%r_m(i)))
            profile%segment(i)=1;profile%side(i)=SIDE_REGULAR
        end do
        if(any(profile%r_m(2:)<=profile%r_m(:n-1)).or.any(profile%m_kg(2:)<profile%m_kg(:n-1))) then
            call fail_status(profile%status,NS_STAGNATION,'numerics','radial or mass progress lost in pressure coordinate');return
        end if
        profile%radius_m=profile%r_m(n);profile%mass_kg=profile%m_kg(n)
        profile%surface_residual_pa=profile%p_pa(n)-options%p_cut_pa
        profile%lapse_residual=profile%nu(n)-0.5_wp*log(1._wp-compact)
        profile%valid=.true.;call clear_status(profile%status)
    contains
        subroutine q_rhs(q,y,dydq,status)
            real(wp),intent(in)::q,y(:)
            real(wp),intent(out)::dydq(:)
            type(status_t),intent(out)::status
            real(wp)::pr,deriv(3),drdq
            dydq=ieee_value(0._wp,ieee_quiet_nan)
            if(y(1)>options%r_limit_m) then
                call fail_status(status,NS_OUT_OF_DOMAIN,'physics','radius exceeds declared integration limit');return
            end if
            call pressure_at_q(q,p0,options%p_cut_pa,qe,pr,status)
            if(.not.is_ok(status)) return
            call tov_rhs(eos,y(1),[y(2),pr,y(3)],deriv,status)
            if(.not.is_ok(status)) return
            if(deriv(2)>=0._wp) then
                call fail_status(status,NS_OUT_OF_DOMAIN,'physics','pressure is not strictly decreasing');return
            end if
            drdq=-pr/deriv(2)
            dydq=[drdq,deriv(1)*drdq,deriv(3)*drdq]
        end subroutine
    end subroutine
end module
