module ns_target_mass
    use ns_precision, only: wp
    use ns_status
    use ns_background_records
    use ns_numerics_types, only: root_options_t
    use ns_bisection, only: bisect
    use ns_tov, only: eos_pressure_callback,integrate_tov
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: solve_target_mass
contains
    subroutine solve_target_mass(eos,m_target,pc_bracket,star_options,options,result,events)
        procedure(eos_pressure_callback)::eos
        real(wp),intent(in)::m_target,pc_bracket(2)
        type(star_options_t),intent(in)::star_options
        type(mass_options_t),intent(in)::options
        type(star_result_t),intent(out)::result
        real(wp),optional,intent(in)::events(:)
        type(star_profile_t)::tmp
        type(root_options_t)::ro
        real(wp)::span,x,pc,tol
        integer::i,n
        result=star_result_t();result%target_mass_kg=m_target
        call fail_status(result%status,NS_INVALID_INPUT,'physics', &
            'finite positive target, increasing Pc bracket and valid scan required')
        if(.not.all(ieee_is_finite([m_target,pc_bracket,options%mass_atol_kg,options%mass_rtol]))) return
        if(m_target<=0._wp.or.pc_bracket(1)<=0._wp.or.pc_bracket(2)<=pc_bracket(1)) return
        if(options%scan_n<3.or.options%scan_n>65.or.len_trim(options%branch_id)==0) return
        if(min(options%mass_atol_kg,options%mass_rtol)<0._wp) return
        tol=options%mass_atol_kg+options%mass_rtol*m_target
        if(.not.ieee_is_finite(tol).or.tol<=0._wp) return
        span=log(pc_bracket(2))-log(pc_bracket(1))
        if(span<=0._wp) then
            call fail_status(result%status,NS_STAGNATION,'numerics','Pc bracket cannot be resolved');return
        end if
        n=options%scan_n
        allocate(result%scan_pc_pa(n),result%scan_mass_kg(n))
        do i=1,n
            x=span*real(i-1,wp)/real(n-1,wp);pc=pc_value(x)
            call compute(pc,tmp)
            if(.not.tmp%valid.or..not.is_ok(tmp%status)) then
                result%status=tmp%status;return
            end if
            result%scan_pc_pa(i)=pc;result%scan_mass_kg(i)=tmp%mass_kg
        end do
        ! A sampled branch-selection guard, NOT a radial-stability theorem.
        if(any(result%scan_mass_kg(2:)<=result%scan_mass_kg(:n-1))) then
            call fail_status(result%status,NS_OUT_OF_DOMAIN,'physics','Pc bracket is not a sampled increasing-mass branch');return
        end if
        result%scan_monotone=.true.
        if(m_target<result%scan_mass_kg(1).or.m_target>result%scan_mass_kg(n)) then
            call fail_status(result%status,NS_NO_BRACKET,'numerics','target outside computed branch mass interval');return
        end if
        ro=options%root;ro%ftol=tol/m_target
        call bisect(residual,0._wp,span,ro,result%root)
        result%status=result%root%status
        if(.not.result%root%valid.or..not.is_ok(result%root%status)) return
        pc=pc_value(result%root%root)
        call compute(pc,result%profile)
        result%status=result%profile%status
        if(.not.result%profile%valid.or..not.is_ok(result%status)) return
        result%profile%branch_id=options%branch_id
        result%mass_residual_kg=result%profile%mass_kg-m_target
        if(abs(result%mass_residual_kg)>tol) then
            call fail_status(result%status,NS_MAX_ITERATIONS,'numerics', &
            'Pc root closed but requested mass residual did not');return
        end if
        result%valid=.true.;call clear_status(result%status)
    contains
        function pc_value(t) result(p)
            real(wp),intent(in)::t
            real(wp)::p
            if(t==0._wp) then
                p=pc_bracket(1)
            else if(t==span) then
                p=pc_bracket(2)
            else
                p=exp(log(pc_bracket(1))+t)
            end if
        end function
        subroutine compute(p,profile)
            real(wp),intent(in)::p
            type(star_profile_t),intent(out)::profile
            if(present(events)) then
                call integrate_tov(eos,p,star_options,profile,events)
            else
                call integrate_tov(eos,p,star_options,profile)
            end if
            result%star_evaluations=result%star_evaluations+1
        end subroutine
        subroutine residual(t,value,status)
            real(wp),intent(in)::t
            real(wp),intent(out)::value
            type(status_t),intent(out)::status
            type(star_profile_t)::profile
            call compute(pc_value(t),profile);status=profile%status
            value=0._wp ! valid only with success status, never used after failure
            if(is_ok(status).and.profile%valid) value=(profile%mass_kg-m_target)/m_target
        end subroutine
    end subroutine
end module
