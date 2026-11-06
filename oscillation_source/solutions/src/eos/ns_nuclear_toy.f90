! Phenomenological zero-temperature nucleonic matter for the teaching path.
! Energy per baryon e excludes the common nucleon rest energy unless stated otherwise.
module ns_nuclear_toy
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use ns_precision, only: wp
    use ns_physical_constants, only: hbar_si,c_si,m_n_si,m_p_si,mev_in_j,inverse_fm3_in_si,pi
    use ns_status
    implicit none
    private
    type, public :: nuclear_model_t
        logical :: ready=.false., regulated=.false.
        real(wp) :: n0_fm3=0._wp,binding_mev=0._wp,k0_mev=0._wp,s0_mev=0._wp,regulator_c=0._wp
        real(wp) :: m_nucleon_kg=0._wp,t0_mev=0._wp,a_mev=0._wp,b_mev=0._wp,sigma=0._wp
    end type
    type, public :: nuclear_state_t
        logical :: valid=.false.
        real(wp) :: n_b_m3=0._wp,n_b_fm3=0._wp,u=0._wp,alpha=0._wp,proton_fraction=0._wp
        real(wp) :: symmetric_energy_mev=0._wp,symmetry_mev=0._wp,energy_excess_mev=0._wp
        real(wp) :: exact_kinetic_mev=0._wp,exact_kinetic_excess_mev=0._wp
        real(wp) :: pressure_nuclear_pa=0._wp,energy_nuclear_j_m3=0._wp,chi_nuclear=0._wp
    end type
    public :: init_nuclear_model,evaluate_nuclear,symmetry_energy_mev,exact_kinetic_per_baryon
    public :: scan_nuclear_causality
contains
    pure real(wp) function common_nucleon_mass()
        common_nucleon_mass=0.5_wp*(m_n_si+m_p_si)
    end function

    pure subroutine poison_state(s)
        type(nuclear_state_t),intent(out)::s
        real(wp)::q
        q=ieee_value(0._wp,ieee_quiet_nan)
        s%n_b_m3=q;s%n_b_fm3=q;s%u=q;s%alpha=q;s%proton_fraction=q
        s%symmetric_energy_mev=q;s%symmetry_mev=q;s%energy_excess_mev=q
        s%exact_kinetic_mev=q;s%exact_kinetic_excess_mev=q
        s%pressure_nuclear_pa=q;s%energy_nuclear_j_m3=q;s%chi_nuclear=q;s%valid=.false.
    end subroutine

    pure subroutine shape_values(u,sigma,c,f,fp,fpp,status)
        real(wp),intent(in)::u,sigma,c
        real(wp),intent(out)::f,fp,fpp
        type(status_t),intent(out)::status
        real(wp)::a,d
        f=0._wp;fp=0._wp;fpp=0._wp
        if(.not.ieee_is_finite(u).or.u<=0._wp.or..not.ieee_is_finite(sigma).or.sigma<=1._wp.or.c<0._wp) then
            call fail_status(status,NS_INVALID_INPUT,'physics','invalid regulated-interaction point');return
        end if
        a=u**(sigma-1._wp);d=1._wp+c*a
        f=u*a/d
        fp=a*(sigma+c*a)/(d*d)
        fpp=(sigma-1._wp)*a/u*(sigma+c*a*(2._wp-sigma))/(d*d*d)
        if(.not.all(ieee_is_finite([f,fp,fpp]))) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','interaction shape overflow');return
        end if
        call clear_status(status)
    end subroutine

    pure subroutine calibration_for_sigma(t0,be,k0,c,sigma,a,b,kcalc,status)
        real(wp),intent(in)::t0,be,k0,c,sigma
        real(wp),intent(out)::a,b,kcalc
        type(status_t),intent(out)::status
        real(wp)::f,fp,fpp,d
        call shape_values(1._wp,sigma,c,f,fp,fpp,status);if(.not.is_ok(status))return
        d=t0/3._wp-be
        if(fp<=f.or.d<=0._wp) then
            call fail_status(status,NS_INVALID_INPUT,'physics','calibration denominator is not positive');return
        end if
        b=(sigma+1._wp)*d/(fp-f)
        a=2._wp*(be-t0-b*f/(sigma+1._wp))
        kcalc=-2._wp*t0+9._wp*b*fpp/(sigma+1._wp)
        if(.not.all(ieee_is_finite([a,b,kcalc]))) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','calibration overflow');return
        end if
        call clear_status(status)
    end subroutine

    subroutine init_nuclear_model(model,regulated,regulator_c,n0_fm3,binding_mev,k0_mev,s0_mev,status)
        type(nuclear_model_t),intent(out)::model
        logical,intent(in)::regulated
        real(wp),intent(in)::regulator_c,n0_fm3,binding_mev,k0_mev,s0_mev
        type(status_t),intent(out)::status
        real(wp)::n0_m3,kf,t0,d,sigma,b,a,kcalc,lo,hi,mid,fa,fb,fm
        integer::iter
        model=nuclear_model_t()
        if(.not.all(ieee_is_finite([regulator_c,n0_fm3,binding_mev,k0_mev,s0_mev]))) then
            call fail_status(status,NS_NONFINITE_INPUT,'data','nuclear parameters must be finite');return
        end if
        if(n0_fm3<=0._wp.or.binding_mev>=0._wp.or.k0_mev<=0._wp.or.s0_mev<=0._wp.or.regulator_c<0._wp.or.regulator_c>1._wp) then
            call fail_status(status,NS_INVALID_INPUT,'physics','invalid nuclear calibration parameters');return
        end if
        n0_m3=n0_fm3*inverse_fm3_in_si
        model%m_nucleon_kg=common_nucleon_mass()
        kf=hbar_si*(3._wp*pi*pi*n0_m3/2._wp)**(1._wp/3._wp)
        t0=(3._wp/5._wp)*kf*kf/(2._wp*model%m_nucleon_kg)/mev_in_j
        if(.not.regulated) then
            d=t0/3._wp-binding_mev
            sigma=(k0_mev+2._wp*t0)/(9._wp*d)
            if(sigma<=1._wp) then
                call fail_status(status,NS_OUT_OF_DOMAIN,'physics','simple model requires sigma>1');return
            end if
            b=(sigma+1._wp)/(sigma-1._wp)*d
            a=binding_mev-(5._wp/3._wp)*t0-b
            kcalc=k0_mev
        else
            lo=1.001_wp;hi=8._wp
            call calibration_for_sigma(t0,binding_mev,k0_mev,regulator_c,lo,a,b,kcalc,status);if(.not.is_ok(status))return
            fa=kcalc-k0_mev
            call calibration_for_sigma(t0,binding_mev,k0_mev,regulator_c,hi,a,b,kcalc,status);if(.not.is_ok(status))return
            fb=kcalc-k0_mev
            if(fa*fb>0._wp) then
                call fail_status(status,NS_NO_BRACKET,'physics','regulated sigma calibration has no bracket');return
            end if
            do iter=1,120
                mid=lo+0.5_wp*(hi-lo)
                call calibration_for_sigma(t0,binding_mev,k0_mev,regulator_c,mid,a,b,kcalc,status);if(.not.is_ok(status))return
                fm=kcalc-k0_mev
                if(abs(fm)<=1e-11_wp*max(1._wp,abs(k0_mev)).or.abs(hi-lo)<=1e-13_wp*max(1._wp,abs(mid))) exit
                if((fa<0._wp.and.fm<0._wp).or.(fa>0._wp.and.fm>0._wp)) then
                    lo=mid;fa=fm
                else
                    hi=mid;fb=fm
                end if
            end do
            sigma=mid
        end if
        model%ready=.true.;model%regulated=regulated;model%regulator_c=merge(regulator_c,0._wp,regulated)
        model%n0_fm3=n0_fm3;model%binding_mev=binding_mev;model%k0_mev=k0_mev;model%s0_mev=s0_mev
        model%t0_mev=t0;model%a_mev=a;model%b_mev=b;model%sigma=sigma
        call clear_status(status)
    end subroutine

    pure subroutine symmetric_terms(model,u,e,ep,epp,status)
        type(nuclear_model_t),intent(in)::model
        real(wp),intent(in)::u
        real(wp),intent(out)::e,ep,epp
        type(status_t),intent(out)::status
        real(wp)::f,fp,fpp,c
        if(.not.model%ready.or.u<=0._wp.or..not.ieee_is_finite(u)) then
            call fail_status(status,NS_INVALID_INPUT,'physics','nuclear model/state not ready');return
        end if
        c=merge(model%regulator_c,0._wp,model%regulated)
        call shape_values(u,model%sigma,c,f,fp,fpp,status);if(.not.is_ok(status))return
        e=model%t0_mev*u**(2._wp/3._wp)+0.5_wp*model%a_mev*u+model%b_mev*f/(model%sigma+1._wp)
        ep=(2._wp/3._wp)*model%t0_mev*u**(-1._wp/3._wp)+0.5_wp*model%a_mev+model%b_mev*fp/(model%sigma+1._wp)
        epp=(-2._wp/9._wp)*model%t0_mev*u**(-4._wp/3._wp)+model%b_mev*fpp/(model%sigma+1._wp)
        call clear_status(status)
    end subroutine

    pure subroutine symmetry_terms(model,u,s,sp,spp,status)
        type(nuclear_model_t),intent(in)::model
        real(wp),intent(in)::u
        real(wp),intent(out)::s,sp,spp
        type(status_t),intent(out)::status
        real(wp)::skin0
        if(.not.model%ready.or.u<=0._wp.or..not.ieee_is_finite(u)) then
            call fail_status(status,NS_INVALID_INPUT,'physics','symmetry-energy point invalid');return
        end if
        skin0=(2._wp**(2._wp/3._wp)-1._wp)*model%t0_mev
        s=skin0*(u**(2._wp/3._wp)-u)+model%s0_mev*u
        sp=(2._wp/3._wp)*skin0*u**(-1._wp/3._wp)-skin0+model%s0_mev
        spp=(-2._wp/9._wp)*skin0*u**(-4._wp/3._wp)
        call clear_status(status)
    end subroutine

    pure subroutine symmetry_energy_mev(model,u,s,status)
        type(nuclear_model_t),intent(in)::model;real(wp),intent(in)::u;real(wp),intent(out)::s
        type(status_t),intent(out)::status;real(wp)::sp,spp
        call symmetry_terms(model,u,s,sp,spp,status)
    end subroutine

    pure subroutine exact_kinetic_per_baryon(model,n_b_m3,alpha,t,status)
        type(nuclear_model_t),intent(in)::model;real(wp),intent(in)::n_b_m3,alpha;real(wp),intent(out)::t
        type(status_t),intent(out)::status
        real(wp)::u,factor
        t=ieee_value(0._wp,ieee_quiet_nan)
        if(.not.model%ready.or..not.all(ieee_is_finite([n_b_m3,alpha])).or.n_b_m3<=0._wp.or.abs(alpha)>1._wp) then
            call fail_status(status,NS_INVALID_INPUT,'physics','exact kinetic query outside domain');return
        end if
        u=n_b_m3/(model%n0_fm3*inverse_fm3_in_si)
        factor=0.5_wp*((1._wp+alpha)**(5._wp/3._wp)+(1._wp-alpha)**(5._wp/3._wp))
        t=model%t0_mev*u**(2._wp/3._wp)*factor
        call clear_status(status)
    end subroutine

    subroutine evaluate_nuclear(model,n_b_m3,alpha,s,status)
        type(nuclear_model_t),intent(in)::model;real(wp),intent(in)::n_b_m3,alpha
        type(nuclear_state_t),intent(out)::s;type(status_t),intent(out)::status
        real(wp)::u,e0,e0p,e0pp,sym,sp,spp,e,ep,epp,dpdu,dedu,mrest_mev,tsym,texact
        call poison_state(s)
        if(.not.model%ready) then;call fail_status(status,NS_UNINITIALIZED,'implementation','nuclear model not ready');return;end if
        if(.not.all(ieee_is_finite([n_b_m3,alpha]))) then
            call fail_status(status,NS_NONFINITE_INPUT,'data','nuclear state must be finite')
            return
        end if
        if(n_b_m3<=0._wp.or.abs(alpha)>1._wp) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics',&
                'nuclear state outside n>0 and |alpha|<=1')
            return
        end if
        u=n_b_m3/(model%n0_fm3*inverse_fm3_in_si)
        call symmetric_terms(model,u,e0,e0p,e0pp,status);if(.not.is_ok(status))return
        call symmetry_terms(model,u,sym,sp,spp,status);if(.not.is_ok(status))return
        e=e0+alpha*alpha*sym;ep=e0p+alpha*alpha*sp;epp=e0pp+alpha*alpha*spp
        mrest_mev=model%m_nucleon_kg*c_si*c_si/mev_in_j
        dpdu=model%n0_fm3*(2._wp*u*ep+u*u*epp)
        dedu=model%n0_fm3*(mrest_mev+e+u*ep)
        s%n_b_m3=n_b_m3;s%n_b_fm3=n_b_m3/inverse_fm3_in_si;s%u=u;s%alpha=alpha;s%proton_fraction=0.5_wp*(1._wp-alpha)
        s%symmetric_energy_mev=e0;s%symmetry_mev=sym;s%energy_excess_mev=e
        call exact_kinetic_per_baryon(model,n_b_m3,alpha,texact,status);if(.not.is_ok(status))return
        tsym=model%t0_mev*u**(2._wp/3._wp);s%exact_kinetic_mev=texact;s%exact_kinetic_excess_mev=texact-tsym
        s%pressure_nuclear_pa=(model%n0_fm3*u*u*ep)*(mev_in_j/(1e-15_wp**3))
        s%energy_nuclear_j_m3=n_b_m3*(model%m_nucleon_kg*c_si*c_si+e*mev_in_j)
        s%chi_nuclear=dpdu/dedu
        if(.not.all(ieee_is_finite([s%pressure_nuclear_pa,s%energy_nuclear_j_m3,s%chi_nuclear]))) then
            call poison_state(s);call fail_status(status,NS_RANGE_ERROR,'numerics','nuclear state overflow');return
        end if
        s%valid=.true.;call clear_status(status)
    end subroutine

    subroutine scan_nuclear_causality(model,u_max,max_chi,causal,status)
        type(nuclear_model_t),intent(in)::model;real(wp),intent(in)::u_max;real(wp),intent(out)::max_chi
        logical,intent(out)::causal;type(status_t),intent(out)::status
        integer::i,j;real(wp)::u,alpha,n;type(nuclear_state_t)::s
        max_chi=-huge(1._wp);causal=.false.
        if(.not.ieee_is_finite(u_max).or.u_max<1._wp) then
            call fail_status(status,NS_INVALID_INPUT,'physics','u_max must be >=1')
            return
        end if
        do j=0,1
            alpha=real(j,wp)
            do i=0,200
                u=0.2_wp+(u_max-0.2_wp)*real(i,wp)/200._wp
                n=u*model%n0_fm3*inverse_fm3_in_si
                call evaluate_nuclear(model,n,alpha,s,status);if(.not.is_ok(status))return
                max_chi=max(max_chi,s%chi_nuclear)
            end do
        end do
        causal=max_chi>0._wp.and.max_chi<=1._wp;call clear_status(status)
    end subroutine
end module ns_nuclear_toy
