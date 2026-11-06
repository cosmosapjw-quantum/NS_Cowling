! Root convergence and smooth-mode identity are independent output states.
module ns_mode_spectrum
    use ns_precision, only: wp
    use ns_physical_constants, only: two_pi
    use ns_status
    use ns_scan_types
    use ns_numerics_types, only: root_options_t,root_result_t
    use ns_scan, only: build_brackets
    use ns_bisection, only: bisect
    use ns_background_access
    use ns_cowling
    use ns_cowling_rk4, only: shoot_cowling_rk4,shoot_cowling_aligned
    use ns_mode_normalization
    use ns_mode_identity
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
    implicit none
    private
    type,public :: mode_options_t
        type(cowling_options_t) :: shot
        type(root_options_t) :: root
        character(len=24) :: method='adaptive'
        integer :: cell_substeps=8
        real(wp)::residual_tol=1e-5_wp,zero_tol=1e-10_wp,node_zero=1e-7_wp
        real(wp)::node_xlo=0.001_wp,node_xhi=0.999_wp
    end type
    type,public :: mode_result_t
        type(root_result_t)::search
        type(cowling_shot_t)::shot
        real(wp),allocatable::x(:),wn(:),vn(:)
        integer::w_nodes=-1,v_nodes=-1
        character(len=16)::label='unresolved'
        type(status_t)::root_status,identity_status
        logical::root_valid=.false.,diagnostics_valid=.false.
    end type
    type,public :: mode_spectrum_t
        type(scan_result_t)::scan
        type(mode_result_t),allocatable::refinements(:)
        integer::roots_found=0,f_index=0,p1_index=0
        logical::valid=.false.
        type(status_t)::status
    end type
    public::refine_mode,scan_modes,evaluate_mode_shot
contains
    subroutine evaluate_mode_shot(bg,omega,options,shot)
        type(background_table_t),intent(in)::bg
        real(wp),intent(in)::omega
        type(mode_options_t),intent(in)::options
        type(cowling_shot_t),intent(out)::shot
        select case(trim(options%method))
        case('adaptive');call shoot_cowling_aligned(bg,omega,options%shot,shot)
        case('adaptive_unaligned');call shoot_cowling(bg,omega,options%shot,shot)
        case('rk4');call shoot_cowling_rk4(bg,omega,options%shot,options%cell_substeps,shot)
        case default
            shot=cowling_shot_t()
            call fail_status(shot%status,NS_INVALID_INPUT,'numerics','unsupported mode integrator')
        end select
    end subroutine

    subroutine refine_mode(bg,candidate,options,result)
        type(background_table_t),intent(in)::bg
        type(bracket_candidate_t),intent(in)::candidate
        type(mode_options_t),intent(in)::options
        type(mode_result_t),intent(out)::result
        type(status_t)::st1,st2
        result=mode_result_t()
        call fail_status(result%identity_status,NS_IDENTITY_UNRESOLVED,'physics','no verified root')
        call fail_status(result%root_status,NS_INVALID_INPUT,'numerics','invalid root settings')
        if(.not.all(ieee_is_finite([options%residual_tol,options%root%xtol,options%root%rtol,options%node_zero]))) return
        if(options%residual_tol<=0._wp.or.options%residual_tol>=1._wp) return
        if(options%root%ftol/=0._wp.or.options%root%xtol<0._wp.or.options%root%rtol<0._wp) return
        if(options%root%xtol==0._wp.and.options%root%rtol==0._wp) return
        if(options%root%max_iter<=0) return
        call fail_status(result%root_status,NS_NO_BRACKET,'numerics','isolated sample is not a root')
        if(candidate%kind/=CAND_SIGN.and.candidate%kind/=CAND_NEAR_BRACKET) return
        if(candidate%lo<=0._wp.or.candidate%lo>=candidate%hi) return
        call bisect(residual,candidate%lo,candidate%hi,options%root,result%search)
        result%root_status=result%search%status
        if(.not.result%search%valid.or..not.is_ok(result%root_status)) return
        call evaluate_mode_shot(bg,result%search%root,options,result%shot)
        result%root_status=result%shot%status
        if(.not.result%shot%valid.or..not.is_ok(result%root_status)) return
        if(abs(result%shot%residual)>options%residual_tol) then
            call fail_status(result%root_status,NS_STAGNATION,'numerics', &
                'coordinate bracket closed but surface residual not resolved')
            return
        end if
        result%root_valid=.true.
        call normalize_mode(result%shot%r_m,result%shot%w_m3,result%shot%v_m2,result%x,result%wn,result%vn,st1)
        if(is_ok(st1)) then
            call count_smooth_nodes(result%x,result%wn,options%node_xlo,options%node_xhi, &
                options%node_zero,result%w_nodes,st1)
            call count_smooth_nodes(result%x,result%vn,options%node_xlo,options%node_xhi, &
                options%node_zero,result%v_nodes,st2)
            result%diagnostics_valid=is_ok(st1).and.is_ok(st2)
        end if
        call identify_smooth_mode(result%root_valid,.true.,result%diagnostics_valid, &
            result%w_nodes,result%v_nodes,result%label,result%identity_status)
    contains
        subroutine residual(omega,value,status)
            real(wp),intent(in)::omega
            real(wp),intent(out)::value
            type(status_t),intent(out)::status
            type(cowling_shot_t)::trial
            value=ieee_value(0._wp,ieee_quiet_nan)
            call evaluate_mode_shot(bg,omega,options,trial);status=trial%status
            if(trial%valid.and.is_ok(status)) value=trial%residual
        end subroutine
    end subroutine

    subroutine scan_modes(bg,grid,options,result)
        type(background_table_t),intent(in)::bg
        real(wp),intent(in)::grid(:)
        type(mode_options_t),intent(in)::options
        type(mode_spectrum_t),intent(out)::result
        type(trial_record_t),allocatable::trials(:)
        type(cowling_shot_t)::shot
        type(status_t)::first_failure
        integer::n,i,j,nf,np
        result=mode_spectrum_t();n=size(grid);nf=0;np=0
        call fail_status(result%status,NS_INVALID_INPUT,'numerics','invalid mode grid/background')
        if(.not.bg%ready.or.n<2.or.n>2000) return
        if(.not.all(ieee_is_finite(grid)).or.any(grid<=0._wp)) return
        do i=2,n
            if(grid(i)<=grid(i-1)) return
        end do
        allocate(trials(n))
        do i=1,n
            call evaluate_mode_shot(bg,grid(i),options,shot)
            trials(i)=trial_record_t(grid(i),shot%residual,int(sign(1._wp,options%shot%amplitude)),shot%valid,shot%status)
        end do
        call build_brackets(trials,options%zero_tol,result%scan)
        result%status=result%scan%status
        if(.not.result%scan%valid) return
        allocate(result%refinements(size(result%scan%candidates)));call clear_status(first_failure)
        do j=1,size(result%refinements)
            call refine_mode(bg,result%scan%candidates(j),options,result%refinements(j))
            if(.not.result%refinements(j)%root_valid) then
                if(is_ok(first_failure)) first_failure=result%refinements(j)%root_status
                cycle
            end if
            result%roots_found=result%roots_found+1
            if(.not.is_ok(result%refinements(j)%identity_status)) cycle
            select case(trim(result%refinements(j)%label))
            case('f');nf=nf+1;result%f_index=j
            case('p1');np=np+1;result%p1_index=j
            end select
        end do
        ! Ambiguous duplicate identities are not resolved by ordering.
        if(nf/=1) result%f_index=0
        if(np/=1) result%p1_index=0
        if(result%f_index>0.and.result%p1_index>0.and.result%scan%all_trials_succeeded.and.is_ok(first_failure)) then
            result%valid=.true.;call clear_status(result%status)
        else if(.not.is_ok(first_failure)) then
            result%status=first_failure
        else
            call fail_status(result%status,NS_IDENTITY_UNRESOLVED,'physics','scan does not establish unique smooth f and p1')
        end if
    end subroutine
end module ns_mode_spectrum
