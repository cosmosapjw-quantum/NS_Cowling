! npe beta equilibrium for the parabolic teaching model. No muons, RMF, or superfluid degrees of freedom.
module ns_beta_equilibrium
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use ns_precision, only: wp
    use ns_physical_constants, only: mev_in_j,inverse_fm3_in_si
    use ns_status
    use ns_electron_gas, only: electron_state_t,electron_state
    use ns_nuclear_toy, only: nuclear_model_t,nuclear_state_t,evaluate_nuclear,symmetry_energy_mev
    implicit none
    private
    type, public :: beta_state_t
        logical :: valid=.false.,charge_neutrality_closed=.false.,causal_local=.false.
        real(wp)::n_b_m3=0._wp,proton_fraction=0._wp,alpha=0._wp,n_p_m3=0._wp,n_e_m3=0._wp
        real(wp)::beta_residual_mev=0._wp,energy_j_m3=0._wp,pressure_pa=0._wp,chi_equilibrium=0._wp
        type(nuclear_state_t)::nuclear
        type(electron_state_t)::electron
        integer::iterations=0
    end type
    public :: solve_beta_equilibrium
contains
    subroutine poison(b)
        type(beta_state_t),intent(out)::b;real(wp)::q
        q=ieee_value(0._wp,ieee_quiet_nan)
        b%n_b_m3=q;b%proton_fraction=q;b%alpha=q;b%n_p_m3=q;b%n_e_m3=q;b%beta_residual_mev=q
        b%energy_j_m3=q;b%pressure_pa=q;b%chi_equilibrium=q
        b%valid=.false.;b%charge_neutrality_closed=.false.
        b%causal_local=.false.;b%iterations=0
    end subroutine

    subroutine residual(model,n,x,r,nuc,elec,status)
        type(nuclear_model_t),intent(in)::model;real(wp),intent(in)::n,x;real(wp),intent(out)::r
        type(nuclear_state_t),intent(out)::nuc;type(electron_state_t),intent(out)::elec;type(status_t),intent(out)::status
        real(wp)::alpha,sym
        if(x<0._wp.or.x>0.5_wp) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics',&
                'proton fraction outside [0,1/2]')
            return
        end if
        alpha=1._wp-2._wp*x
        call evaluate_nuclear(model,n,alpha,nuc,status);if(.not.is_ok(status))return
        call electron_state(n*x,elec,status);if(.not.is_ok(status))return
        call symmetry_energy_mev(model,n/(model%n0_fm3*inverse_fm3_in_si),sym,status);if(.not.is_ok(status))return
        r=4._wp*alpha*sym-elec%mu_j/mev_in_j
        if(.not.ieee_is_finite(r)) call fail_status(status,NS_NONFINITE_RESULT,'numerics','beta residual is nonfinite')
    end subroutine

    subroutine solve_basic(model,n,b,status)
        type(nuclear_model_t),intent(in)::model
        real(wp),intent(in)::n
        type(beta_state_t),intent(out)::b
        type(status_t),intent(out)::status
        real(wp)::lo,hi,mid,fl,fh,fm,tol;integer::iter
        type(nuclear_state_t)::nl,nh,nm;type(electron_state_t)::el,eh,em
        call poison(b)
        if(.not.model%ready) then;call fail_status(status,NS_UNINITIALIZED,'implementation','nuclear model not ready');return;end if
        if(.not.ieee_is_finite(n)) then
            call fail_status(status,NS_NONFINITE_INPUT,'data','baryon density must be finite')
            return
        end if
        if(n<=0._wp) then;call fail_status(status,NS_OUT_OF_DOMAIN,'physics','baryon density must be positive');return;end if
        lo=0._wp;hi=0.5_wp
        call residual(model,n,lo,fl,nl,el,status);if(.not.is_ok(status))return
        call residual(model,n,hi,fh,nh,eh,status);if(.not.is_ok(status))return
        if(fl*fh>0._wp) then
            call fail_status(status,NS_NO_BRACKET,'physics',&
                'beta equilibrium has no proton-fraction bracket')
            return
        end if
        tol=1e-10_wp ! MeV residual, deliberately much tighter than teaching-scale energies.
        do iter=1,160
            mid=lo+0.5_wp*(hi-lo)
            call residual(model,n,mid,fm,nm,em,status);if(.not.is_ok(status))return
            if(abs(fm)<=tol.or.abs(hi-lo)<=2e-14_wp) exit
            if((fl<0._wp.and.fm<0._wp).or.(fl>0._wp.and.fm>0._wp)) then;lo=mid;fl=fm;else;hi=mid;fh=fm;end if
        end do
        if(iter>160) then;call fail_status(status,NS_MAX_ITERATIONS,'numerics','beta equilibrium did not converge');return;end if
        b%n_b_m3=n;b%proton_fraction=mid;b%alpha=1._wp-2._wp*mid;b%n_p_m3=n*mid;b%n_e_m3=n*mid
        b%beta_residual_mev=fm;b%nuclear=nm;b%electron=em;b%energy_j_m3=nm%energy_nuclear_j_m3+em%energy_j_m3
        b%pressure_pa=nm%pressure_nuclear_pa+em%pressure_pa;b%iterations=iter
        b%charge_neutrality_closed=(b%n_p_m3==b%n_e_m3);b%valid=.true.;call clear_status(status)
    end subroutine

    subroutine solve_beta_equilibrium(model,n,b,status)
        type(nuclear_model_t),intent(in)::model
        real(wp),intent(in)::n
        type(beta_state_t),intent(out)::b
        type(status_t),intent(out)::status
        type(beta_state_t)::bm,bp;type(status_t)::st;real(wp)::h,nm,np,de,dp
        call solve_basic(model,n,b,status);if(.not.is_ok(status))return
        h=1e-5_wp;nm=n*(1._wp-h);np=n*(1._wp+h)
        call solve_basic(model,nm,bm,st);if(.not.is_ok(st)) then;call poison(b);status=st;return;end if
        call solve_basic(model,np,bp,st);if(.not.is_ok(st)) then;call poison(b);status=st;return;end if
        de=bp%energy_j_m3-bm%energy_j_m3;dp=bp%pressure_pa-bm%pressure_pa
        if(de<=0._wp.or..not.all(ieee_is_finite([de,dp]))) then
            call poison(b)
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics',&
                'equilibrium response is not monotone')
            return
        end if
        b%chi_equilibrium=dp/de;b%causal_local=b%chi_equilibrium>0._wp.and.b%chi_equilibrium<=1._wp
        if(.not.ieee_is_finite(b%chi_equilibrium)) then
            call poison(b)
            call fail_status(status,NS_NONFINITE_RESULT,'numerics',&
                'equilibrium chi nonfinite')
            return
        end if
        call clear_status(status)
    end subroutine
end module ns_beta_equilibrium
