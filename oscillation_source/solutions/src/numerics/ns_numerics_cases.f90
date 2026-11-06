! Benchmark I/O adapter. It calls this workspace's kernels, never another solver.
module ns_numerics_cases
    use ns_precision, only: wp
    use ns_status
    use ns_numerics_types
    use ns_bisection, only: bisect
    use ns_rk4, only: euler_step,rk4_step
    use ns_adaptive_rk, only: adaptive_integrate,scaled_error_norm,step_factor
    use ns_dopri_stages, only: dopri_trial
    use ns_records, only: emit_real,emit_int,emit_text,emit_bool
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
    implicit none
    private
    public :: is_numerics_case,run_numerics_case
contains
    logical function is_numerics_case(name)
        character(len=*),intent(in)::name
        select case(name)
        case('root_sqrt2','root_probe','rk4_exp','adaptive_exp','ode','control','dopri_pair')
            is_numerics_case=.true.
        case default
            is_numerics_case=.false.
        end select
    end function
    subroutine run_numerics_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token
        character(len=40)::method,problem,action
        real(wp)::lo,hi,xtol,rtol,ftol,root_at,scale,x0,x1,y_initial,y2_initial,rate
        real(wp)::h_initial,h_min,h_max,atol1,atol2,event1,event2,error_value
        real(wp)::old1,old2,new1,new2,error1,error2,val,factor
        real(wp)::x,dx,next_x,max_stage_x,min_stage_x
        real(wp),allocatable::y(:)
        integer::max_iter,n_steps,max_steps,u,ios,i,n,accepted,calls,points,j
        logical::record_trials,use_events,previous_rejected,rollback,crossed,accept_ok
        type(root_options_t)::ro
        type(root_result_t)::rr
        type(ode_options_t)::oo
        type(ode_result_t)::res
        type(step_result_t)::st
        type(status_t)::status
        namelist /request/ request_token
        namelist /parameters/ lo,hi,xtol,rtol,ftol,max_iter,root_at,scale,problem,method,x0,x1, &
            y_initial,y2_initial,rate,n_steps,h_initial,h_min,h_max,atol1,atol2,max_steps, &
            record_trials,use_events,event1,event2,action,error_value,previous_rejected, &
            old1,old2,new1,new2,error1,error2
        request_token='';lo=-1._wp;hi=1._wp;xtol=1e-12_wp;rtol=1e-9_wp;ftol=0._wp
        max_iter=100;root_at=0._wp;scale=1._wp;problem='exp';method='adaptive'
        x0=0._wp;x1=1._wp;y_initial=1._wp;y2_initial=0._wp;rate=1._wp;n_steps=32
        h_initial=0.1_wp;h_min=1e-14_wp;h_max=1._wp;atol1=1e-12_wp;atol2=1e-12_wp;max_steps=10000
        record_trials=.true.;use_events=.false.;event1=0.37_wp;event2=0.8_wp
        action='factor';error_value=0._wp;previous_rejected=.false.
        old1=1._wp;old2=1._wp;new1=1._wp;new2=1._wp;error1=0._wp;error2=0._wp
        if (case_id=='root_sqrt2') then
            lo=1._wp;hi=2._wp;problem='sqrt2';rtol=0._wp
        else if (case_id=='root_probe') then
            problem='shift';rtol=0._wp
        else if (case_id=='rk4_exp') then
            method='rk4'
        end if
        call emit_text('case_id',case_id)
        open(newunit=u,file=path,status='old',action='read',iostat=ios)
        if (ios==0) then
            read(u,nml=request,iostat=ios)
            if (ios==0) read(u,nml=parameters,iostat=ios)
            close(u)
        end if
        call emit_text('request_token',trim(request_token))
        if (ios/=0) then
            call emit_text('status','invalid_input');return
        end if
        if (case_id=='root_sqrt2' .or. case_id=='root_probe') then
            if (problem=='adjacent') then
                lo=1._wp;hi=nearest(1._wp,1._wp)
            end if
            ro=root_options_t(xtol,rtol,ftol,max_iter)
            call bisect(root_function,lo,hi,ro,rr)
            call output_status(rr%status)
            call emit_bool('result_valid',rr%valid)
            call emit_int('iterations',int(rr%iterations,int64));call emit_int('evaluations',int(rr%evaluations,int64))
            call emit_real('lo_final',rr%lo);call emit_real('hi_final',rr%hi)
            call emit_text('termination',trim(rr%termination))
            if (rr%valid) then
                call emit_real('root',rr%root);call emit_real('residual',rr%residual)
            end if
            return
        end if
        if (case_id=='control') then
            if (action=='factor') then
                call step_factor(error_value,previous_rejected,factor,status)
                call output_status(status)
                if (is_ok(status)) call emit_real('factor',factor)
            else
                if (action=='bad_size') then
                    call scaled_error_norm([old1,old2],[new1,new2],[error1],[atol1,atol2],rtol,val,status)
                else
                    call scaled_error_norm([old1,old2],[new1,new2],[error1,error2],[atol1,atol2],rtol,val,status)
                end if
                call output_status(status)
                if (is_ok(status)) call emit_real('error_norm',val)
            end if
            return
        end if
        n=1;if (problem=='oscillator' .or. problem=='two_scale') n=2
        allocate(y(n));y(1)=y_initial;if (n==2) y(2)=y2_initial
        max_stage_x=x0;min_stage_x=x0
        if (case_id=='dopri_pair') then
            call dopri_trial(test_rhs,x0,y,h_initial,st)
            call output_status(st%status);call emit_int('rhs_calls',int(st%evaluations,int64))
            if (is_ok(st%status)) then
                call emit_real('y1',st%y(1));call emit_real('error1',st%error(1))
            end if
            return
        end if
        if (method=='adaptive') then
            oo%atol=[atol1];if(n==2) oo%atol=[atol1,atol2]
            oo%rtol=rtol;oo%h_initial=h_initial;oo%h_min=h_min;oo%h_max=h_max
            oo%max_steps=max_steps;oo%record_trials=record_trials
            if (use_events) then
                call adaptive_integrate(test_rhs,x0,x1,y,oo,res,[event1,event2])
            else
                call adaptive_integrate(test_rhs,x0,x1,y,oo,res)
            end if
            call output_status(res%status);call emit_bool('result_valid',res%valid)
            call emit_int('accepted',int(res%accepted,int64));call emit_int('rejected',int(res%rejected,int64))
            call emit_int('attempts',int(res%attempts,int64));call emit_int('rhs_calls',int(res%evaluations,int64))
            call emit_int('stops_hit',int(res%stops_hit,int64))
            x=x0;points=0
            if (allocated(res%x)) then
                points=size(res%x);x=res%x(points)
                call emit_real('last_accepted_y1',res%y(1,points))
            end if
            call emit_real('last_x',x);call emit_int('profile_points',int(points,int64))
            if (res%valid) then
                call emit_real('y1',res%y(1,points))
                if(n==2) call emit_real('y2',res%y(2,points))
            end if
            rollback=.true.;accept_ok=.true.;crossed=.false.
            if (allocated(res%trial_x)) then
                do i=1,size(res%trial_x)
                    if (res%trial_accepted(i) .and. res%trial_error(i)>1._wp) accept_ok=.false.
                    if (.not.res%trial_accepted(i) .and. i<size(res%trial_x)) then
                        if (res%trial_x(i+1)/=res%trial_x(i) .or. &
                            any(res%trial_y(:,i+1)/=res%trial_y(:,i))) rollback=.false.
                    end if
                    if (use_events) then
                        do j=1,2
                            val=event1;if(j==2) val=event2
                            if ((res%trial_x(i)<val .and. res%trial_x(i)+res%trial_h(i)>val) .or. &
                                (res%trial_x(i)>val .and. res%trial_x(i)+res%trial_h(i)<val)) crossed=.true.
                        end do
                    end if
                end do
            end if
            call emit_bool('trial_audit_available',allocated(res%trial_x))
            if (allocated(res%trial_x)) then
                call emit_bool('rejection_rollback_ok',rollback);call emit_bool('all_accepted_errors_le_one',accept_ok)
                call emit_bool('crossed_stop',crossed)
            end if
        else
            accepted=0;calls=0;x=x0
            call fail_status(status,NS_INVALID_INPUT,'numerics','invalid fixed-step method/count')
            if (n_steps>0 .and. (method=='euler' .or. method=='rk4')) then
                call clear_status(status)
                do i=1,n_steps
                    next_x=x0+(x1-x0)*(real(i,wp)/real(n_steps,wp));if(i==n_steps) next_x=x1
                    dx=next_x-x
                    if (method=='euler') then
                        call euler_step(test_rhs,x,y,dx,st)
                    else
                        call rk4_step(test_rhs,x,y,dx,st)
                    end if
                    calls=calls+st%evaluations;status=st%status
                    if (.not.is_ok(status)) exit
                    x=next_x;y=st%y;accepted=accepted+1
                end do
            end if
            call output_status(status);call emit_bool('result_valid',is_ok(status))
            call emit_real('last_x',x);call emit_int('accepted',int(accepted,int64))
            call emit_int('rhs_calls',int(calls,int64))
            if(is_ok(status)) then
                call emit_real('y1',y(1));if(n==2) call emit_real('y2',y(2))
            end if
        end if
        call emit_real('max_stage_x',max_stage_x);call emit_real('min_stage_x',min_stage_x)
    contains
        subroutine root_function(x,value,status)
            real(wp),intent(in)::x
            real(wp),intent(out)::value
            type(status_t),intent(out)::status
            type(root_options_t)::inner_options
            type(root_result_t)::inner_result
            call clear_status(status)
            select case(trim(problem))
            case('nested')
                inner_options%xtol=1e-14_wp
                call bisect(inner_sqrt,1._wp,2._wp,inner_options,inner_result)
                status=inner_result%status
                if (is_ok(status)) value=inner_result%root-x
            case('sqrt2');value=x*x-2._wp
            case('shift');value=scale*(x-root_at)
            case('wide');value=x/1e308_wp-0.1_wp
            case('adjacent');value=-1._wp;if(x>1._wp) value=1._wp
            case('nan');value=ieee_value(0._wp,ieee_quiet_nan)
            case('failure');call fail_status(status,NS_OUT_OF_DOMAIN,'physics','deliberate scalar domain failure')
            case('failure_mid')
                if(x>0.5_wp .and. x<1.5_wp) then
                    call fail_status(status,NS_OUT_OF_DOMAIN,'physics','deliberate midpoint domain failure')
                else
                    value=x-1._wp
                end if
            case default
                call fail_status(status,NS_INVALID_INPUT,'numerics','unknown scalar benchmark')
            end select
        end subroutine
        subroutine test_rhs(t,state,dydt,status)
            real(wp),intent(in)::t,state(:)
            real(wp),intent(out)::dydt(:)
            type(status_t),intent(out)::status
            call clear_status(status)
            max_stage_x=max(max_stage_x,t);min_stage_x=min(min_stage_x,t)
            select case(trim(problem))
            case('failure_after')
                if (t>0.3_wp) then
                    call fail_status(status,NS_OUT_OF_DOMAIN,'physics','intentional late RHS domain failure')
                else
                    dydt=state
                end if
            case('exp');dydt=rate*state
            case('oscillator');dydt(1)=state(2);dydt(2)=-state(1)
            case('two_scale');dydt(1)=state(1);dydt(2)=-3._wp*state(2)
            case('cubic');dydt=t**3
            case('quartic');dydt=t**4
            case('constant');dydt=1._wp
            case('zero');dydt=0._wp
            case('nan');dydt=ieee_value(0._wp,ieee_quiet_nan)
            case('failure');call fail_status(status,NS_OUT_OF_DOMAIN,'physics','deliberate RHS domain failure')
            case default;call fail_status(status,NS_INVALID_INPUT,'numerics','unknown ODE benchmark')
            end select
        end subroutine
    end subroutine
    subroutine inner_sqrt(x,value,status)
        real(wp),intent(in)::x
        real(wp),intent(out)::value
        type(status_t),intent(out)::status
        value=x*x-2._wp
        call clear_status(status)
    end subroutine
    subroutine output_status(status)
        type(status_t),intent(in)::status
        call emit_text('status',trim(status_name(status%code)))
        call emit_text('status_category',trim(status%category));call emit_text('status_detail',trim(status%detail))
    end subroutine
end module ns_numerics_cases
