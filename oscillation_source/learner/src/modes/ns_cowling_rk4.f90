! Independently stepped RK4 in each fixed background q-cell.
! Shares radial equations and background with the adaptive path, not an oracle.
module ns_cowling_rk4
    use ns_precision, only: wp
    use ns_status
    use ns_background_access
    use ns_cowling, only: cowling_options_t,cowling_shot_t,cowling_rhs_point
    use ns_center_state, only: cowling_center_state
    use ns_surface_residual, only: cowling_surface_residual
    use ns_rk4, only: rk4_step
    use ns_numerics_types, only: step_result_t,ode_result_t
    use ns_adaptive_rk, only: adaptive_integrate
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: shoot_cowling_rk4,shoot_cowling_aligned
contains
    subroutine shoot_cowling_rk4(bg,omega,options,substeps,result)
        type(background_table_t),intent(in)::bg
        real(wp),intent(in)::omega
        type(cowling_options_t),intent(in)::options
        integer,intent(in)::substeps
        type(cowling_shot_t),intent(out)::result
        call integrate_piecewise(bg,omega,options,substeps,.false.,result)
    end subroutine
    subroutine shoot_cowling_aligned(bg,omega,options,result)
        type(background_table_t),intent(in)::bg
        real(wp),intent(in)::omega
        type(cowling_options_t),intent(in)::options
        type(cowling_shot_t),intent(out)::result
        call integrate_piecewise(bg,omega,options,1,.true.,result)
    end subroutine
    subroutine integrate_piecewise(bg,omega,options,substeps,adaptive,result)
        type(background_table_t),intent(in)::bg
        real(wp),intent(in)::omega
        type(cowling_options_t),intent(in)::options
        integer,intent(in)::substeps
        type(cowling_shot_t),intent(out)::result
        type(background_point_t)::p
        type(step_result_t)::step
        type(ode_result_t)::ode
        logical,intent(in)::adaptive
        type(status_t)::st
        integer::i,j,t,n,astat
        real(wp)::q,nextq,h,y(2),a,q0,q1,r0,r1
        result=cowling_shot_t()
        call fail_status(result%status,NS_INVALID_INPUT,'numerics','invalid RK4 Cowling arguments')
        if(.not.bg%ready.or.options%ell/=2) return
        if(.not.all(ieee_is_finite([omega,options%amplitude,options%nu_shift]))) return
        if(omega<=0._wp.or.options%amplitude==0._wp.or.substeps<1.or.substeps>128) return
        if(bg%n<2.or.bg%n>100000) return
        n=(bg%n-1)*substeps+1
        if(n>2000000.or.n>options%ode%max_steps+1) then
            call fail_status(result%status,NS_MAX_STEPS,'numerics','RK4 cell subdivision exceeds max_steps');return
        end if
        call background_bounds(bg,q0,q1,r0,r1,st)
        if(.not.is_ok(st)) then;result%status=st;return;end if
        call cowling_center_state(2,r0,sign(1._wp,options%amplitude),y(1),y(2),st)
        if(.not.is_ok(st)) then;result%status=st;return;end if
        if(adaptive) then
            call adaptive_integrate(rhs_q,q0,q1,y,options%ode,ode,bg%q(2:bg%n-1))
            result%status=ode%status
            result%accepted=ode%accepted;result%rejected=ode%rejected;result%evaluations=ode%evaluations
            if(.not.ode%valid.or..not.is_ok(ode%status)) return
            n=size(ode%x)
        end if
        allocate(result%q(n),result%r_m(n),result%w_m3(n),result%v_m2(n),stat=astat)
        if(astat/=0) then
            call fail_status(result%status,NS_RANGE_ERROR,'environment','RK4 mode allocation failed');return
        end if
        a=abs(options%amplitude);q=q0;t=1
        call save();if(.not.is_ok(result%status)) return
        if(adaptive) then
            do t=1,n
                q=ode%x(t);y=ode%y(:,t)
                call save();if(.not.is_ok(result%status)) return
            end do
        else
        do i=1,bg%n-1
            do j=1,substeps
                nextq=bg%q(i)+(bg%q(i+1)-bg%q(i))*real(j,wp)/real(substeps,wp)
                if(j==substeps) nextq=bg%q(i+1)
                h=nextq-q
                if(q+h/=nextq) then
                    call fail_status(result%status,NS_STAGNATION,'numerics','RK4 endpoint cannot be represented');return
                end if
                call rk4_step(rhs_q,q,y,h,step)
                result%evaluations=result%evaluations+step%evaluations
                if(.not.is_ok(step%status)) then;result%status=step%status;return;end if
                y=step%y;q=nextq;t=t+1;result%accepted=result%accepted+1
                call save();if(.not.is_ok(result%status)) return
            end do
        end do
        end if
        call background_at_q(bg,q1,p,st,options%nu_shift)
        if(.not.is_ok(st)) then;result%status=st;return;end if
        call cowling_surface_residual(p%r_m,omega,p%nu,p%lambda,p%dnu_dr, &
            result%w_m3(n),result%v_m2(n),result%residual,result%q_over_r2,st)
        if(.not.is_ok(st)) then;result%status=st;return;end if
        result%omega_rad_s=omega;result%valid=.true.;call clear_status(result%status)
    contains
        subroutine rhs_q(x,u,du,status)
            real(wp),intent(in)::x,u(:)
            real(wp),intent(out)::du(:)
            type(status_t),intent(out)::status
            type(background_point_t)::pt
            real(wp)::dr(2),qt,factor
            call background_at_q(bg,x,pt,status,options%nu_shift)
            if(.not.is_ok(status)) return
            call cowling_rhs_point(pt%r_m,omega,pt%chi,pt%nu,pt%lambda,pt%dnu_dr,u,dr,qt,status)
            if(.not.is_ok(status)) return
            factor=pt%p_pa/((pt%energy_j_m3+pt%p_pa)*pt%dnu_dr)
            du=factor*dr
            if(.not.all(ieee_is_finite(du))) call fail_status(status,NS_NONFINITE_RESULT,'numerics','RK4 q RHS nonfinite')
        end subroutine
        subroutine save()
            call background_at_q(bg,q,p,result%status)
            if(.not.is_ok(result%status)) return
            result%q(t)=q;result%r_m(t)=p%r_m
            if(any(abs(y)>huge(1._wp)/max(a,1._wp))) then
                call fail_status(result%status,NS_RANGE_ERROR,'numerics','output amplitude overflow');return
            end if
            result%w_m3(t)=a*y(1);result%v_m2(t)=a*y(2)
        end subroutine
    end subroutine
end module ns_cowling_rk4
