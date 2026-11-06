! Zero-temperature ideal relativistic electron gas, SI units.
! Provided infrastructure: students use this helper rather than reimplementing cancellation-prone formulas.
module ns_electron_gas
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use ns_precision, only: wp
    use ns_physical_constants, only: hbar_si, c_si, m_e_si, pi
    use ns_status
    implicit none
    private
    type, public :: electron_state_t
        real(wp) :: n_e_m3=0._wp, p_f_kg_m_s=0._wp, relativity_x=0._wp
        real(wp) :: mu_j=0._wp, energy_j_m3=0._wp, pressure_pa=0._wp
        logical :: valid=.false., used_small_x_series=.false.
    end type
    public :: electron_state
contains
    subroutine poison(s)
        type(electron_state_t),intent(out)::s
        real(wp)::nan
        nan=ieee_value(0._wp,ieee_quiet_nan)
        s%n_e_m3=nan;s%p_f_kg_m_s=nan;s%relativity_x=nan;s%mu_j=nan
        s%energy_j_m3=nan;s%pressure_pa=nan;s%valid=.false.;s%used_small_x_series=.false.
    end subroutine

    subroutine electron_state(n_e_m3,s,status)
        real(wp),intent(in)::n_e_m3
        type(electron_state_t),intent(out)::s
        type(status_t),intent(out)::status
        real(wp)::pf,x,root,mec2,pref,bracket,x2,x4,x6,x8
        call poison(s)
        if(.not.ieee_is_finite(n_e_m3)) then
            call fail_status(status,NS_NONFINITE_INPUT,'data','electron number density must be finite');return
        end if
        if(n_e_m3<0._wp) then
            call fail_status(status,NS_INVALID_INPUT,'physics','electron number density must be nonnegative');return
        end if
        mec2=m_e_si*c_si*c_si
        s%n_e_m3=n_e_m3
        if(n_e_m3==0._wp) then
            s%p_f_kg_m_s=0._wp;s%relativity_x=0._wp;s%mu_j=mec2
            s%energy_j_m3=0._wp;s%pressure_pa=0._wp;s%valid=.true.
            call clear_status(status);return
        end if
        pf=hbar_si*(3._wp*pi*pi*n_e_m3)**(1._wp/3._wp)
        x=pf/(m_e_si*c_si)
        if(.not.ieee_is_finite(pf).or..not.ieee_is_finite(x)) then
            call poison(s);call fail_status(status,NS_RANGE_ERROR,'numerics','electron Fermi momentum overflow');return
        end if
        root=sqrt(1._wp+x*x);s%mu_j=mec2*root;s%p_f_kg_m_s=pf;s%relativity_x=x
        if(x<0.05_wp) then
            x2=x*x;x4=x2*x2;x6=x4*x2;x8=x4*x4
            s%energy_j_m3=n_e_m3*mec2*(1._wp+0.3_wp*x2-(3._wp/56._wp)*x4+(1._wp/48._wp)*x6-(15._wp/1408._wp)*x8)
            s%pressure_pa=n_e_m3*mec2*((1._wp/5._wp)*x2-(1._wp/14._wp)*x4+(1._wp/24._wp)*x6-(5._wp/176._wp)*x8)
            s%used_small_x_series=.true.
        else
            pref=m_e_si**4*c_si**5/(8._wp*pi*pi*hbar_si**3)
            bracket=x*(1._wp+2._wp*x*x)*root-asinh(x)
            s%energy_j_m3=pref*bracket
            s%pressure_pa=n_e_m3*s%mu_j-s%energy_j_m3
        end if
        if(.not.ieee_is_finite(s%energy_j_m3).or.&
            .not.ieee_is_finite(s%pressure_pa).or.&
            s%energy_j_m3<0._wp.or.s%pressure_pa<0._wp) then
            call poison(s);call fail_status(status,NS_RANGE_ERROR,'numerics','electron state is not representable');return
        end if
        s%valid=.true.;call clear_status(status)
    end subroutine
end module ns_electron_gas
