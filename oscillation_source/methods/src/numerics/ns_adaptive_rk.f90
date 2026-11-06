module ns_adaptive_rk
    use ns_precision, only: wp
    use ns_status
    use ns_numerics_types, only: rhs_callback, ode_options_t, ode_result_t, step_result_t
    use ns_dopri_stages, only: dopri_trial
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: adaptive_integrate, scaled_error_norm, step_factor
contains
    subroutine scaled_error_norm(old,new,error,atol,rtol,value,status)
        real(wp),intent(in)::old(:),new(:),error(:),atol(:),rtol
        real(wp),intent(out)::value
        type(status_t),intent(out)::status
        real(wp)::scale,relative,ratio,ea
        integer::i,n
        value=0._wp;n=size(old)
        call fail_status(status,NS_INVALID_INPUT,'numerics','inconsistent error-norm inputs')
        if (n==0 .or. size(new)/=n .or. size(error)/=n .or. size(atol)/=n) return
        if (.not.ieee_is_finite(rtol)) return
        if (rtol<0._wp .or. any(atol<=0._wp)) return
        if (.not.all(ieee_is_finite(old)) .or. .not.all(ieee_is_finite(new))) return
        if (.not.all(ieee_is_finite(error)) .or. .not.all(ieee_is_finite(atol))) return
        do i=1,n
            relative=rtol*max(abs(old(i)),abs(new(i)))
            if (.not.ieee_is_finite(relative)) then
                call fail_status(status,NS_RANGE_ERROR,'numerics','relative scale overflow');return
            end if
            if (relative>huge(scale)-atol(i)) then
                call fail_status(status,NS_RANGE_ERROR,'numerics','absolute+relative scale overflow');return
            end if
            scale=atol(i)+relative;ea=abs(error(i))
            ! A ratio too large to represent still unambiguously rejects the step.
            if (scale<1._wp) then
                if (ea>huge(ea)*scale) then
                    ratio=huge(ea)
                else
                    ratio=ea/scale
                end if
            else
                ratio=ea/scale
            end if
            value=max(value,ratio)
        end do
        call clear_status(status)
    end subroutine

    subroutine step_factor(error_norm,previous_rejected,factor,status)
        real(wp),intent(in)::error_norm
        logical,intent(in)::previous_rejected
        real(wp),intent(out)::factor
        type(status_t),intent(out)::status
        real(wp)::logfactor
        factor=0._wp
        call fail_status(status,NS_INVALID_INPUT,'numerics','nonnegative finite error norm required')
        if (.not.ieee_is_finite(error_norm)) return
        if (error_norm<0._wp) return
        if (error_norm==0._wp) then
            factor=5._wp
        else
            logfactor=log(0.9_wp)-log(error_norm)/5._wp
            factor=exp(max(log(0.2_wp),min(log(5._wp),logfactor)))
        end if
        if (previous_rejected) factor=min(1._wp,factor)
        call clear_status(status)
    end subroutine

    subroutine adaptive_integrate(rhs,x0,x1,y0,options,result,event_schedule)
        procedure(rhs_callback)::rhs
        real(wp),intent(in)::x0,x1,y0(:)
        type(ode_options_t),intent(in)::options
        type(ode_result_t),intent(out)::result
        real(wp),intent(in),optional::event_schedule(:)
        type(step_result_t)::trial
        type(status_t)::local
        real(wp),allocatable::stops(:),current(:)
        real(wp)::x,h,habs,proposed,target,remaining,dir,next_x,err,factor,span
        integer::n,i,stop_index,stop_count,capacity,trace_count,alloc_status
        logical::at_stop,was_rejected
        result=ode_result_t();trace_count=0
        call fail_status(result%status,NS_INVALID_INPUT,'numerics','invalid adaptive inputs/options')
        n=size(y0)
        if (n==0) return
        if (.not.all(ieee_is_finite(y0)) .or. .not.all(ieee_is_finite([x0,x1]))) return
        if (.not.allocated(options%atol)) return
        if (size(options%atol)/=n) return
        if (.not.all(ieee_is_finite(options%atol)) .or. any(options%atol<=0._wp)) return
        if (.not.all(ieee_is_finite([options%rtol,options%h_initial,options%h_min,options%h_max]))) return
        if (options%rtol<0._wp .or. min(options%h_initial,options%h_min,options%h_max)<=0._wp) return
        if (options%h_min>options%h_max .or. options%max_steps<=0) return
        if (options%max_steps>huge(0)/2-1) return
        span=x1-x0
        if (.not.ieee_is_finite(span)) then
            call fail_status(result%status,NS_RANGE_ERROR,'numerics','unrepresentable integration span');return
        end if
        dir=1._wp;if (x1<x0) dir=-1._wp
        stop_count=1
        if (present(event_schedule)) then
            if (.not.all(ieee_is_finite(event_schedule))) return
            stop_count=size(event_schedule)+1
            do i=1,stop_count-1
                if (dir*(event_schedule(i)-x0)<=0._wp .or. dir*(x1-event_schedule(i))<=0._wp) return
                if (i>1) then
                    if (dir*(event_schedule(i)-event_schedule(i-1))<=0._wp) return
                end if
            end do
        end if
        allocate(stops(stop_count),current(n),stat=alloc_status)
        if (alloc_status/=0) then
            call fail_status(result%status,NS_RANGE_ERROR,'environment','allocation failure');return
        end if
        if (present(event_schedule)) stops(:stop_count-1)=event_schedule
        stops(stop_count)=x1
        capacity=min(64,options%max_steps+1)
        call allocate_history(capacity,local)
        if (.not.is_ok(local)) then
            result%status=local;return
        end if
        x=x0;current=y0;result%x(1)=x;result%y(:,1)=current
        if (x0==x1) then
            result%valid=.true.;call clear_status(result%status);call trim_history();return
        end if
        proposed=min(options%h_initial,options%h_max)
        stop_index=1;was_rejected=.false.
        call clear_status(result%status)
        do while (result%attempts<options%max_steps)
            target=stops(stop_index);remaining=abs(target-x)
            habs=min(proposed,remaining);at_stop=habs==remaining
            if (habs<options%h_min .and. .not.at_stop) then
                call fail_status(result%status,NS_STEP_TOO_SMALL,'numerics','requested retry below h_min');exit
            end if
            h=dir*habs;next_x=x+h
            if (at_stop .and. next_x/=target) then
                call fail_status(result%status,NS_STAGNATION,'numerics','endpoint not reached by x+h; reduce step scale')
                exit
            end if
            ! A rounded sum can hit the target even for a nominally smaller h.
            ! Recompute h only when it is exactly representable as target-x.
            if (next_x==target) then
                h=target-x
                if (x+h/=target) then
                    call fail_status(result%status,NS_STAGNATION,'numerics','endpoint subtraction/addition mismatch');exit
                end if
                habs=abs(h);at_stop=.true.
            end if
            if (next_x==x .or. x+h==x) then
                call fail_status(result%status,NS_STAGNATION,'numerics','no representable integration progress');exit
            end if
            result%attempts=result%attempts+1
            call dopri_trial(rhs,x,current,h,trial)
            result%evaluations=result%evaluations+trial%evaluations
            if (.not.is_ok(trial%status)) then
                result%status=trial%status;exit
            end if
            call scaled_error_norm(current,trial%y,trial%error,options%atol,options%rtol,err,local)
            if (.not.is_ok(local)) then
                result%status=local;exit
            end if
            call grow_if_needed(max(result%accepted+2,merge(trace_count+1,0,options%record_trials)),local)
            if (.not.is_ok(local)) then
                result%status=local;exit
            end if
            if (options%record_trials) then
                trace_count=trace_count+1
                result%trial_x(trace_count)=x;result%trial_h(trace_count)=h
                result%trial_y(:,trace_count)=current
                result%trial_error(trace_count)=err;result%trial_accepted(trace_count)=err<=1._wp
            end if
            if (err<=1._wp) then
                x=next_x;current=trial%y;result%accepted=result%accepted+1
                result%x(result%accepted+1)=x;result%y(:,result%accepted+1)=current
                if (at_stop) then
                    if (stop_index==stop_count) then
                        result%valid=.true.;exit
                    end if
                    result%stops_hit=result%stops_hit+1;stop_index=stop_index+1
                end if
                call step_factor(err,was_rejected,factor,local)
                was_rejected=.false.
            else
                result%rejected=result%rejected+1
                call step_factor(err,.true.,factor,local)
                was_rejected=.true.
                if (habs<=options%h_min) then
                    call fail_status(result%status,NS_STEP_TOO_SMALL,'numerics','error test failed at minimum step');exit
                end if
                ! x/current/profile have deliberately NOT advanced.
            end if
            if (.not.is_ok(local)) then
                result%status=local;exit
            end if
            if (habs>options%h_max/factor) then
                proposed=options%h_max
            else
                proposed=min(options%h_max,habs*factor)
            end if
        end do
        if (result%valid) then
            call clear_status(result%status)
        else if (is_ok(result%status)) then
            call fail_status(result%status,NS_MAX_STEPS,'numerics','attempt budget exhausted; accepted prefix retained')
        end if
        call trim_history()
    contains
        subroutine allocate_history(cap,status)
            integer,intent(in)::cap
            type(status_t),intent(out)::status
            integer::ios
            allocate(result%x(cap),result%y(n,cap),stat=ios)
            if (ios==0 .and. options%record_trials) then
                allocate(result%trial_x(cap),result%trial_h(cap),result%trial_error(cap), &
                    result%trial_y(n,cap),result%trial_accepted(cap),stat=ios)
            end if
            if (ios/=0) then
                call fail_status(status,NS_RANGE_ERROR,'environment','history allocation failed')
            else
                call clear_status(status)
            end if
        end subroutine
        subroutine grow_if_needed(need,status)
            integer,intent(in)::need
            type(status_t),intent(out)::status
            real(wp),allocatable::v(:),m(:,:),tx(:),th(:),te(:),ty(:,:)
            logical,allocatable::ta(:)
            integer::cap,ios,np
            call clear_status(status);if (need<=capacity) return
            cap=min(options%max_steps+1,max(need,2*capacity));np=result%accepted+1
            allocate(v(cap),m(n,cap),stat=ios)
            if (ios/=0) then
                call fail_status(status,NS_RANGE_ERROR,'environment','history growth failed');return
            end if
            v(:np)=result%x(:np);m(:,:np)=result%y(:,:np)
            call move_alloc(v,result%x);call move_alloc(m,result%y)
            if (options%record_trials) then
                allocate(tx(cap),th(cap),te(cap),ty(n,cap),ta(cap),stat=ios)
                if (ios/=0) then
                    call fail_status(status,NS_RANGE_ERROR,'environment','trace growth failed');return
                end if
                tx(:trace_count)=result%trial_x(:trace_count);th(:trace_count)=result%trial_h(:trace_count)
                te(:trace_count)=result%trial_error(:trace_count);ty(:,:trace_count)=result%trial_y(:,:trace_count)
                ta(:trace_count)=result%trial_accepted(:trace_count)
                call move_alloc(tx,result%trial_x);call move_alloc(th,result%trial_h);call move_alloc(te,result%trial_error)
                call move_alloc(ty,result%trial_y);call move_alloc(ta,result%trial_accepted)
            end if
            capacity=cap
        end subroutine
        subroutine trim_history()
            result%x=result%x(:result%accepted+1);result%y=result%y(:,:result%accepted+1)
            if (options%record_trials) then
                result%trial_x=result%trial_x(:trace_count);result%trial_h=result%trial_h(:trace_count)
                result%trial_error=result%trial_error(:trace_count);result%trial_y=result%trial_y(:,:trace_count)
                result%trial_accepted=result%trial_accepted(:trace_count)
            end if
        end subroutine
    end subroutine
end module ns_adaptive_rk
