! Analytic-test BVP. No sin/cos or analytic eigenvalues in the solver.
module ns_shooting_toy
    use ns_precision, only: wp
    use ns_status
    use ns_numerics_types
    use ns_scan_types
    use ns_scan, only: build_brackets
    use ns_bisection, only: bisect
    use ns_rk4, only: rk4_step
    use ns_adaptive_rk, only: adaptive_integrate
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
    implicit none
    private
    public :: shoot_string,scan_string,refine_string
contains
    subroutine shoot_string(length_m,k_per_m,options,result,fault)
        real(wp),intent(in)::length_m,k_per_m
        type(shooting_options_t),intent(in)::options
        type(string_shot_t),intent(out)::result
        character(len=*),intent(in),optional::fault
        type(ode_result_t)::ode
        type(step_result_t)::step
        real(wp)::x,next_x,h,current(2),k2,a
        integer::i,n,istat
        result=string_shot_t();a=options%amplitude
        call fail_status(result%status,NS_INVALID_INPUT,'numerics','invalid length, wavenumber, amplitude or method')
        if (.not.all(ieee_is_finite([length_m,k_per_m,a]))) return
        if(length_m<=0._wp .or. k_per_m<0._wp .or. a==0._wp) return
        if(k_per_m>sqrt(huge(k_per_m))) then
            call fail_status(result%status,NS_RANGE_ERROR,'numerics','k squared not representable');return
        end if
        k2=k_per_m*k_per_m
        if(k_per_m>0._wp .and. k2==0._wp) then
            call fail_status(result%status,NS_RANGE_ERROR,'numerics','k squared underflows');return
        end if
        current=[0._wp,1._wp]
        select case(trim(options%method))
        case('adaptive')
            call adaptive_integrate(string_rhs,0._wp,length_m,current,options%ode,ode)
            result%status=ode%status;result%accepted=ode%accepted
            result%rejected=ode%rejected;result%evaluations=ode%evaluations
            if(.not.ode%valid .or. .not.is_ok(ode%status)) return
        case('rk4')
            n=options%fixed_steps
            if(n<1 .or. n>1000000 .or. n>options%ode%max_steps) return
            allocate(ode%x(n+1),ode%y(2,n+1),stat=istat)
            if(istat/=0) then
                call fail_status(result%status,NS_RANGE_ERROR,'environment','toy history allocation failed');return
            end if
            x=0._wp;ode%x(1)=x;ode%y(:,1)=current
            do i=1,n
                next_x=length_m*(real(i,wp)/real(n,wp))
                if(i==n) next_x=length_m
                h=next_x-x
                if(x+h/=next_x) then
                    call fail_status(result%status,NS_STAGNATION,'numerics','fixed step does not reach coordinate');return
                end if
                call rk4_step(string_rhs,x,current,h,step)
                result%evaluations=result%evaluations+step%evaluations;result%status=step%status
                if(.not.is_ok(step%status)) return
                x=next_x;current=step%y
                ode%x(i+1)=x;ode%y(:,i+1)=current;result%accepted=i
            end do
        case default
            return
        end select
        n=size(ode%x)
        if (ode%x(n)/=length_m .or. .not.all(ieee_is_finite(ode%y))) then
            call fail_status(result%status,NS_NONFINITE_RESULT,'numerics','unfinished/nonfinite shot');return
        end if
        result%residual=sign(1._wp,a)*(ode%y(1,n)/length_m)
        if(.not.ieee_is_finite(result%residual)) then
            call fail_status(result%status,NS_RANGE_ERROR,'numerics','residual overflow');return
        end if
        result%x_m=ode%x;result%y_m=a*ode%y(1,:);result%slope=a*ode%y(2,:)
        if (.not.all(ieee_is_finite(result%y_m)) .or. .not.all(ieee_is_finite(result%slope))) then
            call fail_status(result%status,NS_RANGE_ERROR,'numerics','physical amplitude output overflows');return
        end if
        if(any(result%y_m==0._wp .and. ode%y(1,:)/=0._wp) .or. &
           any(result%slope==0._wp .and. ode%y(2,:)/=0._wp)) then
            call fail_status(result%status,NS_RANGE_ERROR,'numerics','physical amplitude output underflows');return
        end if
        result%k_per_m=k_per_m;result%boundary_y_m=result%y_m(n);result%boundary_slope=result%slope(n)
        result%valid=.true.;call clear_status(result%status)
    contains
        subroutine string_rhs(x,y,dydx,status)
            real(wp),intent(in)::x,y(:)
            real(wp),intent(out)::dydx(:)
            type(status_t),intent(out)::status
            call clear_status(status)
            if(present(fault)) then
                if(fault=='failed_rhs' .and. x>0._wp) then
                    call fail_status(status,NS_OUT_OF_DOMAIN,'numerics','deliberately failed RHS');return
                else if(fault=='nan_rhs') then
                    dydx=ieee_value(0._wp,ieee_quiet_nan);return
                end if
            end if
            dydx(1)=y(2);dydx(2)=-k2*y(1)
        end subroutine
    end subroutine

    subroutine refine_string(length_m,candidate,options,result)
        real(wp),intent(in)::length_m
        type(bracket_candidate_t),intent(in)::candidate
        type(shooting_options_t),intent(in)::options
        type(string_root_t),intent(out)::result
        result=string_root_t()
        call fail_status(result%status,NS_INVALID_INPUT,'numerics','invalid refinement options')
        if(.not.all(ieee_is_finite([options%residual_tol,options%root%xtol,options%root%rtol]))) return
        if(options%residual_tol<=0._wp .or. options%root%xtol<0._wp .or. options%root%rtol<0._wp) return
        if(options%root%xtol==0._wp .and. options%root%rtol==0._wp) return
        if(options%root%ftol/=0._wp .or. options%root%max_iter<=0) return
        call fail_status(result%status,NS_NO_BRACKET,'numerics','sampled near-zero value is not a root bracket')
        if(candidate%kind==CAND_SAMPLE) return
        if(candidate%kind/=CAND_SIGN .and. candidate%kind/=CAND_NEAR_BRACKET) return
        if(candidate%lo>=candidate%hi) return
        call bisect(residual_callback,candidate%lo,candidate%hi,options%root,result%search)
        result%status=result%search%status
        if(.not.result%search%valid .or. .not.is_ok(result%status)) return
        call shoot_string(length_m,result%search%root,options,result%shot)
        result%status=result%shot%status
        if(.not.result%shot%valid .or. .not.is_ok(result%status)) return
        if(abs(result%shot%residual)>options%residual_tol) then
            call fail_status(result%status,NS_STAGNATION,'numerics','root width alone insufficient; residual not closed');return
        end if
        result%valid=.true.;call clear_status(result%status)
    contains
        subroutine residual_callback(k,value,status)
            real(wp),intent(in)::k
            real(wp),intent(out)::value
            type(status_t),intent(out)::status
            type(string_shot_t)::shot
            value=ieee_value(0._wp,ieee_quiet_nan)
            call shoot_string(length_m,k,options,shot);status=shot%status
            if(shot%valid .and. is_ok(status)) value=shot%residual
        end subroutine
    end subroutine

    subroutine scan_string(length_m,grid,options,result,failed_index)
        real(wp),intent(in)::length_m,grid(:)
        type(shooting_options_t),intent(in)::options
        type(string_spectrum_t),intent(out)::result
        integer,intent(in),optional::failed_index
        type(trial_record_t),allocatable::trials(:)
        type(string_root_t),allocatable::roots(:)
        type(string_root_t)::rr
        type(status_t)::first_failure
        type(string_shot_t)::shot
        integer::i,n,c,nroots,ori
        result=string_spectrum_t();n=size(grid);nroots=0
        call fail_status(result%status,NS_INVALID_INPUT,'numerics','invalid spectrum grid/options')
        if(n<2 .or. n>100000) return
        if(.not.all(ieee_is_finite(grid)) .or. any(grid<0._wp)) return
        if(.not.all(ieee_is_finite([options%zero_tol,options%residual_tol,options%root%xtol, &
                                   options%root%rtol,options%amplitude,length_m]))) return
        if(options%zero_tol<0._wp .or. options%residual_tol<=0._wp) return
        if(options%root%xtol<0._wp .or. options%root%rtol<0._wp .or. options%root%max_iter<=0) return
        if(options%root%xtol==0._wp .and. options%root%rtol==0._wp) return
        if(options%root%ftol/=0._wp .or. options%amplitude==0._wp .or. length_m<=0._wp) return
        do i=2,n
            if(grid(i)<=grid(i-1)) return
        end do
        allocate(trials(n));ori=1;if(options%amplitude<0._wp) ori=-1
        do i=1,n
            if(present(failed_index)) then
                if(i==failed_index) then
                    trials(i)%parameter=grid(i);trials(i)%orientation=ori
                    call fail_status(trials(i)%status,NS_OUT_OF_DOMAIN,'numerics','deliberately failed scan shot');cycle
                end if
            end if
            call shoot_string(length_m,grid(i),options,shot)
            trials(i)=trial_record_t(grid(i),shot%residual,ori,shot%valid,shot%status)
        end do
        call build_brackets(trials,options%zero_tol,result%scan)
        result%status=result%scan%status
        if(.not.result%scan%valid) return
        c=size(result%scan%candidates);allocate(roots(c),result%refinements(c))
        call clear_status(first_failure)
        do i=1,c
            call refine_string(length_m,result%scan%candidates(i),options,rr);rr%candidate_index=i
            result%refinements(i)=rr
            if(rr%valid) then
                nroots=nroots+1;roots(nroots)=rr
            else
                result%unresolved_candidates=result%unresolved_candidates+1
                if(is_ok(first_failure)) first_failure=rr%status
            end if
        end do
        result%roots=roots(:nroots);result%roots_found=nroots
        if(.not.result%scan%all_trials_succeeded) then
            call fail_status(result%status,NS_OUT_OF_DOMAIN,'numerics','scan contains failed/nonfinite trial barriers');return
        end if
        if(result%unresolved_candidates>0) then
            result%status=first_failure;return
        end if
        if(nroots==0) then
            call fail_status(result%status,NS_NO_BRACKET,'numerics','no admitted root in the finite scan');return
        end if
        result%valid=.true.;call clear_status(result%status)
    end subroutine
end module ns_shooting_toy
