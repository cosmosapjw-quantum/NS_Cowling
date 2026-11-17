! Native Cowling shooting through a sharp Maxwell interface.
! The core and envelope are integrated on separate pressure coordinates; no
! EOS/background interpolation crosses the discontinuity.
module ns_maxwell_cowling
    use ns_precision, only: wp
    use ns_status
    use ns_physical_constants, only: c_si
    use ns_numerics_types, only: ode_result_t
    use ns_adaptive_rk, only: adaptive_integrate
    use ns_background_access, only: background_table_t,background_point_t,background_bounds,background_at_q
    use ns_cowling, only: cowling_options_t,cowling_rhs_point
    use ns_center_state, only: cowling_center_state
    use ns_surface_residual, only: cowling_surface_residual
    use ns_eos_types, only: SIDE_INNER,SIDE_OUTER
    use ns_interface_transfer, only: material_interface_transfer,rapid_interface_transfer
    use ns_maxwell_tov, only: maxwell_star_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
    implicit none
    private

    type,public :: maxwell_cowling_shot_t
        real(wp),allocatable :: q(:),r_m(:),w_m3(:),v_m2(:)
        real(wp) :: omega_rad_s=0._wp,residual=0._wp,q_over_r2=0._wp
        real(wp) :: delta_p_inner_pa=0._wp,delta_p_outer_pa=0._wp
        real(wp) :: w_inner_m3=0._wp,v_inner_m2=0._wp,w_outer_m3=0._wp,v_outer_m2=0._wp
        integer :: accepted=0,rejected=0,evaluations=0
        logical :: valid=.false.
        type(status_t) :: status
    end type

    public :: shoot_maxwell_cowling
contains
    subroutine shoot_maxwell_cowling(star,omega,conversion,options,result)
        type(maxwell_star_t),intent(in) :: star
        real(wp),intent(in) :: omega
        character(len=*),intent(in) :: conversion ! 'slow' or 'rapid'
        type(cowling_options_t),intent(in) :: options
        type(maxwell_cowling_shot_t),intent(out) :: result
        type(ode_result_t) :: core_ode,outer_ode
        type(background_point_t) :: inner_point,outer_point,surface_point
        type(status_t) :: st
        real(wp) :: q0,qti,r0,rti,q0o,qto,r0o,rto
        real(wp) :: y0(2),y_interface(2),w0,v0,amp_scale,orientation
        real(wp) :: dpi,dpo,qterm
        integer :: nc,no,n,astat,i

        result=maxwell_cowling_shot_t()
        call fail_status(result%status,NS_INVALID_INPUT,'physics','invalid native Maxwell Cowling shot')
        if(.not.star%valid.or.options%ell/=2.or..not.ieee_is_finite(omega).or.omega<=0._wp) return
        if(options%amplitude==0._wp.or..not.ieee_is_finite(options%amplitude)) return
        if(.not.allocated(options%ode%atol).or.size(options%ode%atol)/=2) return
        if(.not.all(ieee_is_finite(options%ode%atol)).or.any(options%ode%atol<=0._wp)) return
        if(conversion/='slow'.and.conversion/='rapid') return
        call background_bounds(star%core,q0,qti,r0,rti,st);if(.not.is_ok(st)) then;result%status=st;return;end if
        call background_bounds(star%envelope,q0o,qto,r0o,rto,st);if(.not.is_ok(st)) then;result%status=st;return;end if
        if(abs(rti-star%jump%radius_m)>64._wp*epsilon(1._wp)*max(rti,1._wp).or. &
           abs(r0o-star%jump%radius_m)>64._wp*epsilon(1._wp)*max(r0o,1._wp)) then
            call fail_status(result%status,NS_MATCHING_INCONSISTENT,'physics','Cowling background does not meet at Maxwell interface')
            return
        end if

        amp_scale=abs(options%amplitude);orientation=sign(1._wp,options%amplitude)
        call cowling_center_state(2,r0,orientation,w0,v0,st)
        if(.not.is_ok(st)) then;result%status=st;return;end if
        y0=[w0,v0]
        call adaptive_integrate(core_rhs,q0,qti,y0,options%ode,core_ode)
        result%status=core_ode%status
        if(.not.core_ode%valid.or..not.is_ok(result%status)) return
        call background_at_q(star%core,qti,inner_point,st,options%nu_shift)
        if(.not.is_ok(st)) then;result%status=st;return;end if
        call background_at_q(star%envelope,q0o,outer_point,st,options%nu_shift)
        if(.not.is_ok(st)) then;result%status=st;return;end if
        if(conversion=='slow') then
            call material_interface_transfer(star%jump%radius_m,omega,star%jump%pressure_pa, &
                inner_point%energy_j_m3,outer_point%energy_j_m3,inner_point%nu,inner_point%lambda, &
                inner_point%dnu_dr,core_ode%y(1,size(core_ode%x)),core_ode%y(2,size(core_ode%x)), &
                SIDE_INNER,SIDE_OUTER,y_interface(1),y_interface(2),dpi,dpo,st)
        else
            call rapid_interface_transfer(star%jump%radius_m,omega,star%jump%pressure_pa, &
                inner_point%energy_j_m3,outer_point%energy_j_m3,inner_point%nu,inner_point%lambda, &
                inner_point%dnu_dr,core_ode%y(1,size(core_ode%x)),core_ode%y(2,size(core_ode%x)), &
                SIDE_INNER,SIDE_OUTER,y_interface(1),y_interface(2),dpi,dpo,st)
        end if
        if(.not.is_ok(st)) then;result%status=st;return;end if
        result%delta_p_inner_pa=amp_scale*dpi;result%delta_p_outer_pa=amp_scale*dpo
        result%w_inner_m3=amp_scale*core_ode%y(1,size(core_ode%x))
        result%v_inner_m2=amp_scale*core_ode%y(2,size(core_ode%x))
        result%w_outer_m3=amp_scale*y_interface(1)
        result%v_outer_m2=amp_scale*y_interface(2)

        call adaptive_integrate(envelope_rhs,q0o,qto,y_interface,options%ode,outer_ode)
        result%status=outer_ode%status
        if(.not.outer_ode%valid.or..not.is_ok(result%status)) return
        nc=size(core_ode%x);no=size(outer_ode%x);n=nc+no-1
        allocate(result%q(n),result%r_m(n),result%w_m3(n),result%v_m2(n),stat=astat)
        if(astat/=0) then
            call fail_status(result%status,NS_RANGE_ERROR,'environment','Maxwell Cowling result allocation failed')
            return
        end if
        do i=1,nc
            call background_at_q(star%core,core_ode%x(i),surface_point,st,options%nu_shift)
            if(.not.is_ok(st)) then;result%status=st;return;end if
            result%q(i)=core_ode%x(i);result%r_m(i)=surface_point%r_m
            result%w_m3(i)=amp_scale*core_ode%y(1,i)
            result%v_m2(i)=amp_scale*core_ode%y(2,i)
        end do
        do i=2,no
            call background_at_q(star%envelope,outer_ode%x(i),surface_point,st,options%nu_shift)
            if(.not.is_ok(st)) then;result%status=st;return;end if
            result%q(nc+i-1)=qti+outer_ode%x(i);result%r_m(nc+i-1)=surface_point%r_m
            result%w_m3(nc+i-1)=amp_scale*outer_ode%y(1,i)
            result%v_m2(nc+i-1)=amp_scale*outer_ode%y(2,i)
        end do
        call background_at_q(star%envelope,qto,surface_point,st,options%nu_shift)
        if(.not.is_ok(st)) then;result%status=st;return;end if
        call cowling_surface_residual(surface_point%r_m,omega,surface_point%nu,surface_point%lambda, &
            surface_point%dnu_dr,result%w_m3(n),result%v_m2(n),result%residual,result%q_over_r2,st)
        if(.not.is_ok(st)) then;result%status=st;return;end if
        result%accepted=core_ode%accepted+outer_ode%accepted
        result%rejected=core_ode%rejected+outer_ode%rejected
        result%evaluations=core_ode%evaluations+outer_ode%evaluations
        result%omega_rad_s=omega;result%valid=.true.;call clear_status(result%status)
    contains
        subroutine core_rhs(q,y,dydq,status)
            real(wp),intent(in) :: q,y(:)
            real(wp),intent(out) :: dydq(:)
            type(status_t),intent(out) :: status
            call phase_rhs(star%core,q,omega,options%nu_shift,y,dydq,status)
        end subroutine
        subroutine envelope_rhs(q,y,dydq,status)
            real(wp),intent(in) :: q,y(:)
            real(wp),intent(out) :: dydq(:)
            type(status_t),intent(out) :: status
            call phase_rhs(star%envelope,q,omega,options%nu_shift,y,dydq,status)
        end subroutine
    end subroutine

    subroutine phase_rhs(background,q,omega,nu_shift,y,dydq,status)
        type(background_table_t),intent(in) :: background
        real(wp),intent(in) :: q,omega,nu_shift,y(:)
        real(wp),intent(out) :: dydq(:)
        type(status_t),intent(out) :: status
        type(background_point_t) :: p
        real(wp) :: dydr(2),qterm,drdq
        dydq=ieee_value(0._wp,ieee_quiet_nan)
        call background_at_q(background,q,p,status,nu_shift);if(.not.is_ok(status)) return
        call cowling_rhs_point(p%r_m,omega,p%chi,p%nu,p%lambda,p%dnu_dr,y,dydr,qterm,status)
        if(.not.is_ok(status)) return
        drdq=p%p_pa/((p%energy_j_m3+p%p_pa)*p%dnu_dr)
        if(.not.ieee_is_finite(drdq).or.drdq<=0._wp) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','Maxwell Cowling dr/dq invalid');return
        end if
        dydq=dydr*drdq
        if(.not.all(ieee_is_finite(dydq))) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','Maxwell Cowling q derivative nonfinite')
            return
        end if
        call clear_status(status)
    end subroutine
end module ns_maxwell_cowling
