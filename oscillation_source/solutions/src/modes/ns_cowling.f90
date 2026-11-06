module ns_cowling
    use ns_precision, only: wp
    use ns_physical_constants, only: c_si
    use ns_status
    use ns_numerics_types, only: ode_options_t,ode_result_t
    use ns_adaptive_rk, only: adaptive_integrate
    use ns_background_access
    use ns_center_state, only: cowling_center_state
    use ns_surface_residual, only: cowling_surface_residual
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    type,public :: cowling_options_t
        integer :: ell=2
        real(wp) :: amplitude=1._wp,nu_shift=0._wp
        type(ode_options_t) :: ode
    end type
    type,public :: cowling_shot_t
        real(wp),allocatable :: q(:),r_m(:),w_m3(:),v_m2(:)
        real(wp) :: omega_rad_s=0._wp,residual=0._wp,q_over_r2=0._wp
        integer :: accepted=0,rejected=0,evaluations=0
        logical :: valid=.false.
        type(status_t) :: status
    end type
    public :: cowling_rhs_point,shoot_cowling
contains
    subroutine cowling_rhs_point(r,omega,chi,nu,lambda,dnu,y,dydr,qterm,status)
        real(wp),intent(in)::r,omega,chi,nu,lambda,dnu,y(:)
        real(wp),intent(out)::dydr(:),qterm
        type(status_t),intent(out)::status
        real(wp)::k,expo,el,term,L
        dydr=0._wp;qterm=0._wp
        call fail_status(status,NS_INVALID_INPUT,'physics','Cowling RHS requires ell=2 smooth barotropic state')
        if(size(y)/=2.or.size(dydr)/=2) return
        if(.not.all(ieee_is_finite([r,omega,chi,nu,lambda,dnu])).or..not.all(ieee_is_finite(y))) return
        if(r<=0._wp.or.omega<=0._wp.or.chi<=0._wp.or.chi>1._wp) return
        k=omega/c_si;expo=lambda-2._wp*nu
        if(expo>log(huge(1._wp)).or.expo<log(tiny(1._wp))) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','Cowling exponential unrepresentable');return
        end if
        el=exp(lambda);term=(k*r)**2*exp(expo)*y(2)
        if(.not.ieee_is_finite(term)) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','frequency term overflow');return
        end if
        if(abs(dnu)>0._wp.and.abs(y(1))>huge(y(1))/abs(dnu)) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','metric term overflow');return
        end if
        qterm=term+dnu*y(1);L=6._wp
        dydr(1)=qterm/chi-L*el*y(2)
        dydr(2)=2._wp*dnu*y(2)-el*y(1)/(r*r)
        if(.not.all(ieee_is_finite(dydr)).or..not.ieee_is_finite(qterm)) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','Cowling derivative nonfinite');return
        end if
        call clear_status(status)
    end subroutine

    subroutine shoot_cowling(background,omega,options,result)
        type(background_table_t),intent(in)::background
        real(wp),intent(in)::omega
        type(cowling_options_t),intent(in)::options
        type(cowling_shot_t),intent(out)::result
        type(ode_result_t)::ode
        type(status_t)::st
        type(background_point_t)::p
        real(wp)::w0,v0,q0,q1,r0,r1,qterm,drdq,amp_scale,amp_orientation
        real(wp)::y0(2)
        integer::i,n,astat
        result=cowling_shot_t()
        call fail_status(result%status,NS_INVALID_INPUT,'physics','invalid Cowling shot')
        if(.not.background%ready.or.options%ell/=2.or.omega<=0._wp.or.options%amplitude==0._wp) return
        if(.not.ieee_is_finite(omega).or..not.ieee_is_finite(options%amplitude).or..not.ieee_is_finite(options%nu_shift)) return
        if(.not.allocated(options%ode%atol).or.size(options%ode%atol)/=2) return
        call background_bounds(background,q0,q1,r0,r1,st);if(.not.is_ok(st)) then;result%status=st;return;end if
        amp_scale=abs(options%amplitude);amp_orientation=sign(1._wp,options%amplitude)
        call cowling_center_state(2,r0,amp_orientation,w0,v0,st)
        if(.not.is_ok(st)) then;result%status=st;return;end if
        y0=[w0,v0]
        call adaptive_integrate(rhs_q,q0,q1,y0,options%ode,ode)
        result%status=ode%status;result%accepted=ode%accepted;result%rejected=ode%rejected;result%evaluations=ode%evaluations
        if(.not.ode%valid.or..not.is_ok(ode%status)) return
        n=size(ode%x);allocate(result%q(n),result%r_m(n),result%w_m3(n),result%v_m2(n),stat=astat)
        if(astat/=0) then
            call fail_status(result%status,NS_RANGE_ERROR,'environment','Cowling profile allocation failed');return
        end if
        do i=1,n
            call background_at_q(background,ode%x(i),p,st,options%nu_shift);if(.not.is_ok(st)) then;result%status=st;return;end if
            result%q(i)=ode%x(i);result%r_m(i)=p%r_m
            if ((ode%y(1,i)/=0._wp .and. amp_scale>huge(1._wp)/abs(ode%y(1,i))) .or. &
                (ode%y(2,i)/=0._wp .and. amp_scale>huge(1._wp)/abs(ode%y(2,i)))) then
                call fail_status(result%status,NS_RANGE_ERROR,'numerics','requested Cowling output amplitude is unrepresentable')
                return
            end if
            result%w_m3(i)=amp_scale*ode%y(1,i);result%v_m2(i)=amp_scale*ode%y(2,i)
            if (.not.all(ieee_is_finite([result%w_m3(i),result%v_m2(i)]))) then
                call fail_status(result%status,NS_RANGE_ERROR,'numerics','requested Cowling output amplitude is unrepresentable')
                return
            end if
        end do
        call background_at_q(background,q1,p,st,options%nu_shift);if(.not.is_ok(st)) then;result%status=st;return;end if
        call cowling_surface_residual(p%r_m,omega,p%nu,p%lambda,p%dnu_dr, &
             result%w_m3(n),result%v_m2(n),result%residual,result%q_over_r2,st)
        if(.not.is_ok(st)) then;result%status=st;return;end if
        result%omega_rad_s=omega;result%valid=.true.;call clear_status(result%status)
    contains
        subroutine rhs_q(q,y,dydq,status)
            real(wp),intent(in)::q,y(:)
            real(wp),intent(out)::dydq(:)
            type(status_t),intent(out)::status
            real(wp)::dydr(2)
            call background_at_q(background,q,p,status,options%nu_shift);if(.not.is_ok(status)) return
            call cowling_rhs_point(p%r_m,omega,p%chi,p%nu,p%lambda,p%dnu_dr,y,dydr,qterm,status);if(.not.is_ok(status)) return
            drdq=p%p_pa/((p%energy_j_m3+p%p_pa)*p%dnu_dr)
            if(.not.ieee_is_finite(drdq).or.drdq<=0._wp) then
                call fail_status(status,NS_NONFINITE_RESULT,'numerics','dr/dq invalid');return
            end if
            dydq=dydr*drdq
            if(.not.all(ieee_is_finite(dydq))) then
                call fail_status(status,NS_NONFINITE_RESULT,'numerics','Cowling q derivative nonfinite')
            end if
        end subroutine
    end subroutine
end module ns_cowling
