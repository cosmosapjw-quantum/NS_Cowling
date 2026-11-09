! Phase-segmented Cowling spectrum scan for one sharp Maxwell interface.
!
! A root is a converged zero of the native two-phase shooting residual. Mode
! family labels are deliberately not assigned here: node counting and
! continuation across slow/rapid conversion branches remain separate gates.
module ns_maxwell_mode_spectrum
    use ns_precision, only: wp
    use ns_status
    use ns_scan_types, only: trial_record_t,bracket_candidate_t,scan_result_t, &
        CAND_SIGN,CAND_NEAR_BRACKET,CAND_SAMPLE
    use ns_numerics_types, only: root_result_t
    use ns_scan, only: build_brackets
    use ns_bisection, only: bisect
    use ns_maxwell_tov, only: maxwell_star_t
    use ns_maxwell_cowling, only: maxwell_cowling_shot_t,shoot_maxwell_cowling
    use ns_mode_spectrum, only: mode_options_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
    implicit none
    private

    type, public :: maxwell_node_count_t
        integer :: w_nodes_inner=-1,w_nodes_outer=-1,w_nodes_total=-1
        integer :: v_nodes_inner=-1,v_nodes_outer=-1,v_nodes_total=-1
        integer :: samples_inner=0,samples_outer=0
        logical :: valid=.false.
    end type

    type, public :: maxwell_mode_root_t
        type(bracket_candidate_t) :: candidate
        type(root_result_t) :: search
        type(maxwell_cowling_shot_t) :: shot
        type(maxwell_node_count_t) :: nodes
        integer :: candidate_index=0
        logical :: root_valid=.false.
        character(len=16) :: label='unresolved'
        type(status_t) :: root_status,identity_status,node_status
    end type

    type, public :: maxwell_mode_spectrum_t
        type(scan_result_t) :: scan
        type(maxwell_mode_root_t), allocatable :: refinements(:)
        integer :: roots_found=0,unresolved_candidates=0
        logical :: valid=.false.
        type(status_t) :: status
    end type

    public :: scan_maxwell_modes,count_maxwell_cowling_nodes
contains
    subroutine count_maxwell_cowling_nodes(shot,interface_radius,nodes,status,relative_zero_tol)
        type(maxwell_cowling_shot_t),intent(in) :: shot
        real(wp),intent(in) :: interface_radius
        type(maxwell_node_count_t),intent(out) :: nodes
        type(status_t),intent(out) :: status
        real(wp),intent(in),optional :: relative_zero_tol
        real(wp) :: zero_tol,radius_tol
        integer :: n,n_inner,i

        nodes=maxwell_node_count_t()
        call fail_status(status,NS_INVALID_INPUT,'numerics', &
            'phase-aware node count requires finite piecewise Cowling W profile and interface radius')
        if(.not.ieee_is_finite(interface_radius).or.interface_radius<=0._wp) return
        if(.not.allocated(shot%r_m).or..not.allocated(shot%w_m3).or..not.allocated(shot%v_m2)) return
        n=size(shot%r_m)
        if(n<4.or.size(shot%w_m3)/=n.or.size(shot%v_m2)/=n) return
        if(.not.all(ieee_is_finite(shot%r_m)).or..not.all(ieee_is_finite(shot%w_m3)).or. &
           .not.all(ieee_is_finite(shot%v_m2))) return
        if(any(shot%r_m<=0._wp)) return
        do i=2,n
            if(shot%r_m(i)<=shot%r_m(i-1)) return
        end do
        radius_tol=64._wp*epsilon(1._wp)*max(abs(interface_radius),1._wp)
        if(interface_radius<=shot%r_m(1)+radius_tol.or. &
           interface_radius>=shot%r_m(n)-radius_tol) return
        if(present(relative_zero_tol)) then
            zero_tol=relative_zero_tol
            if(.not.ieee_is_finite(zero_tol).or.zero_tol<0._wp.or.zero_tol>=0.1_wp) return
        else
            zero_tol=128._wp*epsilon(1._wp)
        end if

        ! The interface sample belongs to the inner segment. The first
        ! envelope sample is counted independently; no sign comparison is
        ! ever made between the one-sided interface values.
        n_inner=count(shot%r_m<=interface_radius+radius_tol)
        nodes%samples_inner=n_inner
        nodes%samples_outer=n-n_inner
        if(n_inner<2.or.n-n_inner<2) then
            call fail_status(status,NS_STAGNATION,'numerics', &
                'phase-aware node count is under-resolved in one phase segment')
            return
        end if
        call count_segment_nodes(shot%w_m3(1:n_inner),zero_tol,nodes%w_nodes_inner,status)
        if(.not.is_ok(status)) return
        call count_segment_nodes(shot%w_m3(n_inner+1:n),zero_tol,nodes%w_nodes_outer,status)
        if(.not.is_ok(status)) return
        call count_segment_nodes(shot%v_m2(1:n_inner),zero_tol,nodes%v_nodes_inner,status)
        if(.not.is_ok(status)) return
        call count_segment_nodes(shot%v_m2(n_inner+1:n),zero_tol,nodes%v_nodes_outer,status)
        if(.not.is_ok(status)) return
        nodes%w_nodes_total=nodes%w_nodes_inner+nodes%w_nodes_outer
        nodes%v_nodes_total=nodes%v_nodes_inner+nodes%v_nodes_outer
        nodes%valid=.true.
        call clear_status(status)
    end subroutine count_maxwell_cowling_nodes

    subroutine count_segment_nodes(values,relative_zero_tol,nodes,status)
        real(wp),intent(in) :: values(:),relative_zero_tol
        integer,intent(out) :: nodes
        type(status_t),intent(out) :: status
        real(wp) :: scale,threshold
        integer :: i,last_sign,current_sign

        nodes=-1
        call fail_status(status,NS_INVALID_INPUT,'numerics','node-count segment must be finite and nonempty')
        if(size(values)<1.or..not.all(ieee_is_finite(values))) return
        scale=maxval(abs(values))
        if(.not.ieee_is_finite(scale)) return
        nodes=0
        if(scale<=tiny(1._wp)) then
            call clear_status(status)
            return
        end if
        threshold=relative_zero_tol*scale
        last_sign=0
        do i=1,size(values)
            if(abs(values(i))<=threshold) cycle
            current_sign=merge(1,-1,values(i)>0._wp)
            if(last_sign/=0.and.current_sign/=last_sign) nodes=nodes+1
            last_sign=current_sign
        end do
        call clear_status(status)
    end subroutine count_segment_nodes

    subroutine scan_maxwell_modes(star,grid,conversion,options,result)
        type(maxwell_star_t),intent(in) :: star
        real(wp),intent(in) :: grid(:)
        character(len=*),intent(in) :: conversion
        type(mode_options_t),intent(in) :: options
        type(maxwell_mode_spectrum_t),intent(out) :: result
        type(trial_record_t),allocatable :: trials(:)
        type(maxwell_cowling_shot_t) :: shot
        type(status_t) :: first_failure
        integer :: n,i,j,istat
        logical :: all_roots_resolved

        result=maxwell_mode_spectrum_t()
        n=size(grid)
        call fail_status(result%status,NS_INVALID_INPUT,'numerics', &
            'Maxwell spectrum requires a valid star and a strictly increasing positive frequency grid')
        if(.not.star%valid.or.n<2.or.n>2000) return
        if(conversion/='slow'.and.conversion/='rapid') return
        if(.not.all(ieee_is_finite(grid)).or.any(grid<=0._wp)) return
        do i=2,n
            if(grid(i)<=grid(i-1)) return
        end do
        if(options%shot%ell/=2.or..not.ieee_is_finite(options%shot%amplitude).or. &
           options%shot%amplitude==0._wp) return
        if(options%method/='adaptive') then
            call fail_status(result%status,NS_INVALID_INPUT,'numerics', &
                'Maxwell spectrum currently uses the adaptive phase-segmented integrator')
            return
        end if
        if(.not.all(ieee_is_finite([options%zero_tol,options%residual_tol,options%root%xtol, &
            options%root%rtol,options%root%ftol]))) return
        if(options%zero_tol<0._wp.or.options%residual_tol<=0._wp.or.options%residual_tol>=1._wp) return
        if(min(options%root%xtol,options%root%rtol,options%root%ftol)<0._wp) return
        if(options%root%xtol==0._wp.and.options%root%rtol==0._wp) return
        if(options%root%max_iter<=0) return

        allocate(trials(n),stat=istat)
        if(istat/=0) then
            call fail_status(result%status,NS_RANGE_ERROR,'environment','Maxwell spectrum trial allocation failed')
            return
        end if
        call clear_status(first_failure)
        do i=1,n
            call shoot_maxwell_cowling(star,grid(i),conversion,options%shot,shot)
            trials(i)%parameter=grid(i)
            trials(i)%orientation=merge(1,-1,options%shot%amplitude>0._wp)
            trials(i)%residual=shot%residual
            trials(i)%valid=shot%valid
            trials(i)%status=shot%status
            if(.not.shot%valid.or..not.is_ok(shot%status)) then
                if(is_ok(first_failure)) first_failure=shot%status
            end if
        end do
        call build_brackets(trials,options%zero_tol,result%scan)
        result%status=result%scan%status
        if(.not.result%scan%valid) return
        if(.not.result%scan%all_trials_succeeded) then
            if(is_ok(first_failure)) then
                call fail_status(result%status,NS_STAGNATION,'numerics', &
                    'Maxwell frequency scan contains failed trials and is incomplete')
            else
                result%status=first_failure
            end if
            return
        end if

        allocate(result%refinements(size(result%scan%candidates)),stat=istat)
        if(istat/=0) then
            call fail_status(result%status,NS_RANGE_ERROR,'environment','Maxwell root allocation failed')
            return
        end if
        all_roots_resolved=.true.
        do j=1,size(result%refinements)
            result%refinements(j)%candidate_index=j
            result%refinements(j)%candidate=result%scan%candidates(j)
            call fail_status(result%refinements(j)%identity_status,NS_IDENTITY_UNRESOLVED,'physics', &
                'Maxwell root family is not identified by this scan')
            if(result%scan%candidates(j)%kind==CAND_SAMPLE) then
                call fail_status(result%refinements(j)%root_status,NS_NO_BRACKET,'numerics', &
                    'near-zero sample has no isolating frequency bracket')
                result%unresolved_candidates=result%unresolved_candidates+1
                all_roots_resolved=.false.
                if(is_ok(first_failure)) first_failure=result%refinements(j)%root_status
                cycle
            end if
            if(result%scan%candidates(j)%kind/=CAND_SIGN.and. &
               result%scan%candidates(j)%kind/=CAND_NEAR_BRACKET) then
                call fail_status(result%refinements(j)%root_status,NS_INVALID_INPUT,'numerics', &
                    'unsupported Maxwell frequency candidate')
                all_roots_resolved=.false.
                if(is_ok(first_failure)) first_failure=result%refinements(j)%root_status
                cycle
            end if
            if(result%scan%candidates(j)%lo<=0._wp.or. &
               result%scan%candidates(j)%lo>=result%scan%candidates(j)%hi) then
                call fail_status(result%refinements(j)%root_status,NS_NO_BRACKET,'numerics', &
                    'Maxwell frequency candidate is not an isolating interval')
                all_roots_resolved=.false.
                if(is_ok(first_failure)) first_failure=result%refinements(j)%root_status
                cycle
            end if
            call bisect(residual,result%scan%candidates(j)%lo,result%scan%candidates(j)%hi, &
                options%root,result%refinements(j)%search)
            result%refinements(j)%root_status=result%refinements(j)%search%status
            if(.not.result%refinements(j)%search%valid.or. &
               .not.is_ok(result%refinements(j)%root_status)) then
                all_roots_resolved=.false.
                if(is_ok(first_failure)) first_failure=result%refinements(j)%root_status
                cycle
            end if
            call shoot_maxwell_cowling(star,result%refinements(j)%search%root,conversion, &
                options%shot,result%refinements(j)%shot)
            result%refinements(j)%root_status=result%refinements(j)%shot%status
            if(.not.result%refinements(j)%shot%valid.or. &
               .not.is_ok(result%refinements(j)%root_status)) then
                all_roots_resolved=.false.
                if(is_ok(first_failure)) first_failure=result%refinements(j)%root_status
                cycle
            end if
            if(abs(result%refinements(j)%shot%residual)>options%residual_tol) then
                call fail_status(result%refinements(j)%root_status,NS_STAGNATION,'numerics', &
                    'Maxwell frequency bracket closed but surface residual is unresolved')
                all_roots_resolved=.false.
                if(is_ok(first_failure)) first_failure=result%refinements(j)%root_status
                cycle
            end if
            result%refinements(j)%root_valid=.true.
            call count_maxwell_cowling_nodes(result%refinements(j)%shot,star%jump%radius_m, &
                result%refinements(j)%nodes,result%refinements(j)%node_status)
            result%roots_found=result%roots_found+1
        end do

        if(all_roots_resolved) then
            result%valid=.true.
            call clear_status(result%status)
        else
            if(is_ok(first_failure)) then
                call fail_status(result%status,NS_STAGNATION,'numerics', &
                    'one or more Maxwell frequency candidates remain unresolved')
            else
                result%status=first_failure
            end if
        end if
    contains
        subroutine residual(omega,value,status)
            real(wp),intent(in) :: omega
            real(wp),intent(out) :: value
            type(status_t),intent(out) :: status
            type(maxwell_cowling_shot_t) :: trial
            value=ieee_value(0._wp,ieee_quiet_nan)
            call shoot_maxwell_cowling(star,omega,conversion,options%shot,trial)
            status=trial%status
            if(trial%valid.and.is_ok(status)) value=trial%residual
        end subroutine
    end subroutine
end module ns_maxwell_mode_spectrum
