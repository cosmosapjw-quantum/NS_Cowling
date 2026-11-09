! Native Fortran phase-segmented TOV background for one sharp Maxwell jump.
!
! The two EOS callbacks are never mixed in one interpolation. The high-density
! core is integrated to P_t, then the low-density envelope is restarted at the
! same (r,m,nu). Each phase owns a separate monotone-log-pressure background.
module ns_maxwell_tov
    use ns_precision, only: wp
    use ns_status
    use ns_physical_constants, only: g_si,c_si,pi
    use ns_eos_types, only: eos_state_t,SIDE_INNER,SIDE_OUTER
    use ns_numerics_types, only: ode_options_t,ode_result_t
    use ns_adaptive_rk, only: adaptive_integrate
    use ns_background_access, only: background_table_t
    use ns_background_records, only: star_options_t
    use ns_tov, only: eos_pressure_callback
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
    implicit none
    private

    type, public :: maxwell_interface_t
        real(wp) :: pressure_pa=0._wp, radius_m=0._wp, mass_kg=0._wp
        real(wp) :: nu=0._wp, lambda=0._wp, dnu_dr=0._wp
        type(eos_state_t) :: inner,outer
        real(wp) :: delta_energy_j_m3=0._wp,delta_density_m3=0._wp
        real(wp) :: chemical_potential_j=0._wp,chemical_potential_mismatch_j=0._wp
        real(wp) :: thermodynamic_residual_j_m3=0._wp
        real(wp) :: seidov_ratio=0._wp,seidov_critical=0._wp
        logical :: rapid_onset_stable_diagnostic=.false.
    end type

    type, public :: maxwell_star_t
        type(background_table_t) :: core,envelope
        type(maxwell_interface_t) :: jump
        real(wp) :: pc_pa=0._wp,p_cut_pa=0._wp
        real(wp) :: radius_m=0._wp,mass_kg=0._wp,lapse_shift=0._wp
        integer :: accepted=0,rejected=0,evaluations=0
        logical :: valid=.false.
        type(status_t) :: status
    end type

    public :: integrate_maxwell_tov
contains
    subroutine integrate_maxwell_tov(eos_high,eos_low,pc,pt,options,star)
        procedure(eos_pressure_callback) :: eos_high,eos_low
        real(wp),intent(in) :: pc,pt
        type(star_options_t),intent(in) :: options
        type(maxwell_star_t),intent(out) :: star
        type(ode_options_t) :: core_options,outer_options
        type(ode_result_t) :: core_traj,outer_traj
        type(eos_state_t) :: high_t,low_t
        type(status_t) :: st
        real(wp) :: y0(3),yt(3),p_start,q_core,q_outer,compact,shift,mu_t,mu_res,thermo_res,scale
        integer :: ncore,nouter

        star=maxwell_star_t()
        star%pc_pa=pc;star%p_cut_pa=options%p_cut_pa
        call fail_status(star%status,NS_INVALID_INPUT,'physics', &
            'Maxwell TOV requires finite Pc>Pt>Pcut>0 and declared surface policy')
        if(.not.all(ieee_is_finite([pc,pt,options%p_cut_pa,options%r0_m,options%r_limit_m,options%nu_start]))) return
        if(options%surface_policy/='finite_cutoff_schwarzschild_time') return
        if(len_trim(options%eos_id)==0.or.options%eos_id=='unspecified') return
        if(pc<=pt.or.pt<=options%p_cut_pa.or.options%r0_m<=0._wp.or. &
           options%r_limit_m<=options%r0_m) return
        if(.not.allocated(options%ode%atol).or.size(options%ode%atol)/=3) return
        if(.not.all(ieee_is_finite(options%ode%atol)).or.any(options%ode%atol<=0._wp)) return
        if(options%ode%h_initial<=0._wp.or.options%ode%h_min<=0._wp.or. &
           options%ode%h_max<options%ode%h_min.or.options%ode%max_steps<=0) return

        call center_start_phase(eos_high,SIDE_INNER,pc,options%r0_m,options%nu_start,y0,star%status)
        if(.not.is_ok(star%status)) return
        p_start=y0(2)
        if(p_start<=pt) then
            call fail_status(star%status,NS_OUT_OF_DOMAIN,'physics', &
                'central series endpoint must remain above the Maxwell pressure')
            return
        end if
        q_core=log(p_start)-log(pt)
        q_outer=log(pt)-log(options%p_cut_pa)
        if(.not.all(ieee_is_finite([q_core,q_outer])).or.min(q_core,q_outer)<=0._wp) then
            call fail_status(star%status,NS_STAGNATION,'numerics','Maxwell phase pressure spans are unresolved')
            return
        end if

        core_options=options%ode
        core_options%h_initial=min(core_options%h_initial,0.1_wp*(pc-y0(2))/pc,q_core)
        call integrate_segment(eos_high,SIDE_INNER,p_start,pt,q_core,options%r_limit_m, &
            y0,core_options,core_traj,star%status)
        if(.not.is_ok(star%status).or..not.core_traj%valid) return

        ncore=size(core_traj%x)
        if(ncore<2) then
            call fail_status(star%status,NS_STAGNATION,'numerics','Maxwell core segment has fewer than two points')
            return
        end if
        yt=core_traj%y(:,ncore)
        outer_options=options%ode
        outer_options%h_initial=min(outer_options%h_initial,q_outer)
        call integrate_segment(eos_low,SIDE_OUTER,pt,options%p_cut_pa,q_outer,options%r_limit_m, &
            yt,outer_options,outer_traj,star%status)
        if(.not.is_ok(star%status).or..not.outer_traj%valid) return
        nouter=size(outer_traj%x)
        if(nouter<2) then
            call fail_status(star%status,NS_STAGNATION,'numerics','Maxwell envelope segment has fewer than two points')
            return
        end if

        call build_background(eos_high,SIDE_INNER,p_start,pt,q_core,core_traj,star%core,star%status)
        if(.not.is_ok(star%status)) return
        call build_background(eos_low,SIDE_OUTER,pt,options%p_cut_pa,q_outer,outer_traj,star%envelope,star%status)
        if(.not.is_ok(star%status)) return

        star%radius_m=star%envelope%r_m(nouter)
        star%mass_kg=star%envelope%m_kg(nouter)
        compact=2._wp*(g_si/c_si**2)*(star%mass_kg/star%radius_m)
        if(.not.ieee_is_finite(compact).or.compact<=0._wp.or.compact>=1._wp) then
            call fail_status(star%status,NS_OUT_OF_DOMAIN,'physics','invalid Maxwell-star surface compactness')
            return
        end if
        shift=0.5_wp*log(1._wp-compact)-star%envelope%nu(nouter)
        star%core%nu=star%core%nu+shift
        star%envelope%nu=star%envelope%nu+shift
        star%lapse_shift=shift

        call eos_high(pt,SIDE_INNER,high_t,star%status);if(.not.is_ok(star%status)) return
        call eos_low(pt,SIDE_OUTER,low_t,star%status);if(.not.is_ok(star%status)) return
        ! Enforce the defining Maxwell coexistence conditions at the exact
        ! interface: equal pressure, equal baryon chemical potential, and
        ! the zero-temperature identity Delta epsilon = mu_t Delta n_B.
        call fail_status(star%status,NS_MATCHING_INCONSISTENT,'physics', &
            'EOS callbacks do not satisfy thermodynamic Maxwell coexistence')
        if(.not.all(ieee_is_finite([high_t%pressure_pa,low_t%pressure_pa, &
            high_t%energy_j_m3,low_t%energy_j_m3,high_t%n_b_m3,low_t%n_b_m3,high_t%mu_j,low_t%mu_j]))) return
        if(min(high_t%energy_j_m3,low_t%energy_j_m3,high_t%n_b_m3,low_t%n_b_m3, &
            high_t%mu_j,low_t%mu_j)<=0._wp) return
        scale=max(abs(pt),tiny(1._wp))
        if(abs(high_t%pressure_pa-pt)>1e-10_wp*scale.or.abs(low_t%pressure_pa-pt)>1e-10_wp*scale) return
        mu_t=0.5_wp*(high_t%mu_j+low_t%mu_j)
        mu_res=high_t%mu_j-low_t%mu_j
        if(abs(mu_res)>1e-10_wp*max(abs(mu_t),tiny(1._wp))) return
        thermo_res=(high_t%energy_j_m3-low_t%energy_j_m3)- &
                   mu_t*(high_t%n_b_m3-low_t%n_b_m3)
        scale=max(abs(high_t%energy_j_m3-low_t%energy_j_m3), &
                  abs(mu_t*(high_t%n_b_m3-low_t%n_b_m3)),tiny(1._wp))
        if(abs(thermo_res)>1e-9_wp*scale) return
        star%jump%pressure_pa=pt
        star%jump%radius_m=star%core%r_m(ncore)
        star%jump%mass_kg=star%core%m_kg(ncore)
        star%jump%nu=star%core%nu(ncore)
        star%jump%lambda=star%core%lambda(ncore)
        star%jump%dnu_dr=star%core%dnu_dr(ncore)
        star%jump%inner=high_t;star%jump%outer=low_t
        star%jump%delta_energy_j_m3=high_t%energy_j_m3-low_t%energy_j_m3
        star%jump%delta_density_m3=high_t%n_b_m3-low_t%n_b_m3
        star%jump%chemical_potential_j=mu_t
        star%jump%chemical_potential_mismatch_j=mu_res
        star%jump%thermodynamic_residual_j_m3=thermo_res
        if(star%jump%delta_energy_j_m3<=0._wp.or.star%jump%delta_density_m3<=0._wp) then
            call fail_status(star%status,NS_MATCHING_INCONSISTENT,'physics','Maxwell jump ordering is not positive')
            return
        end if
        star%jump%seidov_ratio=star%jump%delta_energy_j_m3/low_t%energy_j_m3
        star%jump%seidov_critical=0.5_wp+1.5_wp*pt/low_t%energy_j_m3
        star%jump%rapid_onset_stable_diagnostic= &
            star%jump%seidov_ratio<star%jump%seidov_critical

        ! The same geometric state is used to restart the outer phase. Verify
        ! the metric and pressure-gradient quantities that must be continuous.
        call fail_status(star%status,NS_MATCHING_INCONSISTENT,'physics','Maxwell TOV interface geometry is discontinuous')
        if(star%envelope%r_m(1)/=star%jump%radius_m) return
        if(abs(star%envelope%m_kg(1)-star%jump%mass_kg)> &
           64._wp*epsilon(1._wp)*max(abs(star%jump%mass_kg),1._wp)) return
        if(abs(star%envelope%nu(1)-star%jump%nu)> &
           64._wp*epsilon(1._wp)*max(abs(star%jump%nu),1._wp)) return
        if(abs(star%envelope%lambda(1)-star%jump%lambda)> &
           64._wp*epsilon(1._wp)*max(abs(star%jump%lambda),1._wp)) return
        if(abs(star%envelope%dnu_dr(1)-star%jump%dnu_dr)> &
           128._wp*epsilon(1._wp)*max(abs(star%jump%dnu_dr),1._wp)) return
        if(star%core%p_pa(ncore)/=pt.or.star%envelope%p_pa(1)/=pt) return

        star%accepted=core_traj%accepted+outer_traj%accepted
        star%rejected=core_traj%rejected+outer_traj%rejected
        star%evaluations=core_traj%evaluations+outer_traj%evaluations
        star%valid=.true.;call clear_status(star%status)
    end subroutine

    subroutine center_start_phase(eos,side,pc,r0,nu0,y,status)
        procedure(eos_pressure_callback) :: eos
        integer,intent(in) :: side
        real(wp),intent(in) :: pc,r0,nu0
        real(wp),intent(out) :: y(3)
        type(status_t),intent(out) :: status
        type(eos_state_t) :: s,starter
        real(wp) :: alpha,drop
        y=ieee_value(0._wp,ieee_quiet_nan)
        call fail_status(status,NS_INVALID_INPUT,'physics','finite central pressure, radius and lapse required')
        if(.not.all(ieee_is_finite([pc,r0,nu0])).or.pc<=0._wp.or.r0<=0._wp) return
        call eos(pc,side,s,status);if(.not.is_ok(status)) return
        if(.not.s%valid.or.s%energy_j_m3<=0._wp) then
            call fail_status(status,NS_INVALID_INPUT,'physics','invalid central EOS state');return
        end if
        alpha=(2._wp*pi*g_si/3._wp)*((s%energy_j_m3+pc)/c_si**2)* &
              ((s%energy_j_m3+3._wp*pc)/c_si**2)
        drop=alpha*r0**2
        if(.not.ieee_is_finite(drop).or.drop<=0._wp.or.drop>=0.01_wp*pc) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','central TOV expansion is not resolved')
            return
        end if
        y=[(4._wp*pi/3._wp)*(s%energy_j_m3/c_si**2)*r0**3,pc-drop,nu0]
        if(y(2)==pc.or..not.all(ieee_is_finite(y))) then
            call fail_status(status,NS_STAGNATION,'numerics','central pressure correction is not representable')
            return
        end if
        call eos(y(2),side,starter,status)
        if(.not.is_ok(status)) return
        if(starter%pressure_pa/=y(2)) then
            call fail_status(status,NS_INVALID_INPUT,'physics','central EOS callback returned a different pressure')
            return
        end if
        call clear_status(status)
    end subroutine

    subroutine integrate_segment(eos,side,p_start,p_end,q_end,r_limit,y_start,options,traj,status)
        procedure(eos_pressure_callback) :: eos
        integer,intent(in) :: side
        real(wp),intent(in) :: p_start,p_end,q_end,r_limit,y_start(3)
        type(ode_options_t),intent(in) :: options
        type(ode_result_t),intent(out) :: traj
        type(status_t),intent(out) :: status
        type(ode_options_t) :: local_options
        call fail_status(status,NS_INVALID_INPUT,'numerics','invalid Maxwell segment integration request')
        if(.not.all(ieee_is_finite([p_start,p_end,q_end,r_limit])).or. &
           p_start<=p_end.or.p_end<=0._wp.or.q_end<=0._wp.or.r_limit<=0._wp) return
        local_options=options
        call adaptive_integrate(rhs_q,0._wp,q_end,y_start,local_options,traj)
        status=traj%status
        if(.not.traj%valid.or..not.is_ok(status)) return
        if(any(traj%y(1,:)>r_limit)) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','Maxwell TOV radius limit exceeded')
            traj%valid=.false.;return
        end if
        call clear_status(status)
    contains
        subroutine rhs_q(q,y,dydq,st)
            real(wp),intent(in) :: q,y(:)
            real(wp),intent(out) :: dydq(:)
            type(status_t),intent(out) :: st
            type(eos_state_t) :: state
            real(wp) :: p,compact,phi,drdq,dm
            dydq=ieee_value(0._wp,ieee_quiet_nan)
            call pressure_at_q(q,p_start,p_end,q_end,p,st);if(.not.is_ok(st)) return
            if(size(y)/=3.or..not.all(ieee_is_finite(y)).or.y(1)<=0._wp.or.y(2)<0._wp) then
                call fail_status(st,NS_INVALID_INPUT,'physics','invalid Maxwell TOV stage state');return
            end if
            call eos(p,side,state,st);if(.not.is_ok(st)) return
            if(.not.state%valid.or..not.all(ieee_is_finite([state%energy_j_m3,state%pressure_pa]))) then
                call fail_status(st,NS_INVALID_INPUT,'physics','invalid one-sided Maxwell EOS state');return
            end if
            if(state%energy_j_m3<=0._wp.or.abs(state%pressure_pa/p-1._wp)>256._wp*epsilon(1._wp)) then
                call fail_status(st,NS_MATCHING_INCONSISTENT,'physics','Maxwell EOS state does not match requested pressure');return
            end if
            compact=2._wp*g_si*y(2)/(y(1)*c_si**2)
            if(.not.ieee_is_finite(compact).or.compact>=1._wp) then
                call fail_status(st,NS_OUT_OF_DOMAIN,'physics','Maxwell TOV trapped-surface denominator nonpositive');return
            end if
            phi=(g_si/c_si**2)*(y(2)/y(1)**2+4._wp*pi*y(1)*p/c_si**2)/(1._wp-compact)
            dm=4._wp*pi*y(1)**2*state%energy_j_m3/c_si**2
            drdq=p/((state%energy_j_m3+p)*phi)
            if(.not.all(ieee_is_finite([phi,dm,drdq])).or.phi<=0._wp.or.drdq<=0._wp) then
                call fail_status(st,NS_NONFINITE_RESULT,'numerics','Maxwell TOV q derivative invalid');return
            end if
            dydq=[drdq,dm*drdq,phi*drdq]
            if(.not.all(ieee_is_finite(dydq))) then
                call fail_status(st,NS_NONFINITE_RESULT,'numerics','Maxwell TOV derivative nonfinite');return
            end if
            call clear_status(st)
        end subroutine
    end subroutine

    subroutine build_background(eos,side,p_start,p_end,q_end,traj,bg,status)
        procedure(eos_pressure_callback) :: eos
        integer,intent(in) :: side
        real(wp),intent(in) :: p_start,p_end,q_end
        type(ode_result_t),intent(in) :: traj
        type(background_table_t),intent(out) :: bg
        type(status_t),intent(out) :: status
        type(eos_state_t) :: s
        real(wp) :: p,compact,phi
        integer :: i,n,astat
        bg=background_table_t();n=size(traj%x)
        call fail_status(status,NS_INVALID_INPUT,'data','invalid Maxwell segment trajectory')
        if(.not.traj%valid.or.n<2.or.size(traj%y,1)/=3.or.size(traj%y,2)/=n) return
        allocate(bg%q(n),bg%r_m(n),bg%m_kg(n),bg%p_pa(n),bg%energy_j_m3(n), &
            bg%n_b_m3(n),bg%mu_j(n),bg%chi(n),bg%nu(n),bg%lambda(n),bg%dnu_dr(n),stat=astat)
        if(astat/=0) then
            call fail_status(status,NS_RANGE_ERROR,'environment','Maxwell background allocation failed');return
        end if
        do i=1,n
            call pressure_at_q(traj%x(i),p_start,p_end,q_end,p,status);if(.not.is_ok(status)) return
            call eos(p,side,s,status);if(.not.is_ok(status)) return
            if(.not.s%valid.or..not.all(ieee_is_finite([s%energy_j_m3,s%n_b_m3,s%mu_j,s%chi]))) then
                call fail_status(status,NS_INVALID_INPUT,'physics','Maxwell background EOS fields invalid');return
            end if
            if(s%energy_j_m3<=0._wp.or.s%n_b_m3<=0._wp.or.s%mu_j<=0._wp.or.s%chi<=0._wp.or.s%chi>1._wp) then
                call fail_status(status,NS_OUT_OF_DOMAIN,'physics','Maxwell Cowling background requires 0<chi<=1');return
            end if
            compact=2._wp*g_si*traj%y(2,i)/(traj%y(1,i)*c_si**2)
            if(compact<=0._wp.or.compact>=1._wp) then
                call fail_status(status,NS_OUT_OF_DOMAIN,'physics','Maxwell background compactness invalid');return
            end if
            phi=(g_si/c_si**2)*(traj%y(2,i)/traj%y(1,i)**2+4._wp*pi*traj%y(1,i)*p/c_si**2)/(1._wp-compact)
            bg%q(i)=traj%x(i);bg%r_m(i)=traj%y(1,i);bg%m_kg(i)=traj%y(2,i);bg%p_pa(i)=p
            bg%energy_j_m3(i)=s%energy_j_m3;bg%n_b_m3(i)=s%n_b_m3;bg%mu_j(i)=s%mu_j
            bg%chi(i)=s%chi;bg%nu(i)=traj%y(3,i);bg%lambda(i)=-0.5_wp*log(1._wp-compact);bg%dnu_dr(i)=phi
        end do
        if(any(bg%q(2:)<=bg%q(:n-1)).or.any(bg%r_m(2:)<=bg%r_m(:n-1)).or. &
           any(bg%m_kg(2:)<bg%m_kg(:n-1)).or.any(bg%p_pa(2:)>=bg%p_pa(:n-1))) then
            call fail_status(status,NS_STAGNATION,'numerics','Maxwell segment lost monotone progress');return
        end if
        bg%n=n;bg%ready=.true.;call clear_status(status)
    end subroutine

    subroutine pressure_at_q(q,p_start,p_end,q_end,p,status)
        real(wp),intent(in) :: q,p_start,p_end,q_end
        real(wp),intent(out) :: p
        type(status_t),intent(out) :: status
        p=ieee_value(0._wp,ieee_quiet_nan)
        call fail_status(status,NS_OUT_OF_DOMAIN,'numerics','Maxwell pressure coordinate outside segment')
        if(.not.ieee_is_finite(q).or.q<0._wp.or.q>q_end) return
        if(q==0._wp) then
            p=p_start
        else if(q==q_end) then
            p=p_end
        else
            p=exp(log(p_start)-q)
        end if
        if(.not.ieee_is_finite(p).or.p<=0._wp) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','Maxwell pressure unrepresentable');return
        end if
        call clear_status(status)
    end subroutine
end module ns_maxwell_tov
