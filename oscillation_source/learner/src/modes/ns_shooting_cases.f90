! Case I/O and synthetic failure fixtures, not an alternate solver.
module ns_shooting_cases
    use ns_precision, only: wp
    use ns_status
    use ns_scan_types
    use ns_numerics_types, only: root_options_t
    use ns_scan, only: linear_scan_grid,build_brackets
    use ns_shooting_toy, only: shoot_string,scan_string,refine_string
    use ns_records, only: emit_real,emit_int,emit_text,emit_bool
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
    implicit none
    private
    public :: is_shooting_case,run_shooting_case
contains
    logical function is_shooting_case(name)
        character(len=*),intent(in)::name
        is_shooting_case=name=='string_refine' .or. name=='string_shot' .or. name=='string_spectrum' .or. &
            name=='bracket_probe' .or. name=='scan_with_failed_trial'
    end function
    subroutine run_shooting_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token,method,special_input,scenario
        character(len=48)::key
        real(wp)::length_m,k_per_m,k_min,k_max,amplitude,ode_rtol,atol_y_m,atol_slope
        real(wp)::h_initial_m,h_max_m,root_xtol,root_rtol,zero_tol,residual_tol,width,resid
        real(wp),allocatable::grid(:)
        integer::scan_n,n_steps,max_steps,max_iter,unit,ios,i,n,fail_index
        logical::emit_profile
        type(shooting_options_t)::options
        type(string_shot_t)::shot
        type(string_root_t)::refined
        type(bracket_candidate_t)::candidate
        type(string_spectrum_t)::spectrum
        type(scan_result_t)::scan
        type(status_t)::status
        namelist /request/ request_token
        namelist /parameters/ length_m,k_per_m,k_min,k_max,scan_n,amplitude,method,n_steps,ode_rtol, &
            atol_y_m,atol_slope,h_initial_m,h_max_m,max_steps,root_xtol,root_rtol,max_iter, &
            zero_tol,residual_tol,emit_profile,special_input,scenario
        request_token='';length_m=1._wp;k_per_m=2._wp;k_min=0.2_wp;k_max=8._wp;scan_n=200
        amplitude=1._wp;method='adaptive';n_steps=256;ode_rtol=1e-10_wp
        atol_y_m=1e-12_wp;atol_slope=1e-12_wp;h_initial_m=0.01_wp;h_max_m=0.1_wp;max_steps=10000
        root_xtol=1e-11_wp;root_rtol=0._wp;max_iter=100;zero_tol=1e-10_wp;residual_tol=1e-8_wp
        emit_profile=.false.;special_input='none';scenario='gap'
        call emit_text('case_id',case_id)
        open(newunit=unit,file=path,status='old',action='read',iostat=ios)
        if(ios==0) then
            read(unit,nml=request,iostat=ios)
            if(ios==0) read(unit,nml=parameters,iostat=ios)
            close(unit)
        end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0) then
            call emit_text('status','invalid_input');return
        end if
        if(case_id=='bracket_probe' .or. case_id=='scan_with_failed_trial') then
            call probe(trim(scenario),zero_tol,scan)
            call print_status(scan%status);call emit_bool('result_valid',scan%valid)
            if(scan%valid) then
                call print_scan(scan)
                call emit_int('cross_failure_brackets',int(count_crossed_failures(scan),int64))
            end if
            return
        end if
        if(special_input=='nan_k') k_per_m=ieee_value(0._wp,ieee_quiet_nan)
        options%amplitude=amplitude;options%method=method;options%fixed_steps=n_steps
        options%zero_tol=zero_tol;options%residual_tol=residual_tol
        options%root=root_options_t(root_xtol,root_rtol,0._wp,max_iter)
        options%ode%atol=[atol_y_m,atol_slope];options%ode%rtol=ode_rtol
        options%ode%h_initial=h_initial_m;options%ode%h_max=h_max_m;options%ode%h_min=1e-14_wp
        options%ode%max_steps=max_steps
        if(case_id=='string_refine') then
            candidate=bracket_candidate_t(k_min,k_max,1,2,CAND_SIGN)
            if(special_input=='sample_candidate') candidate%kind=CAND_SAMPLE
            if(special_input=='unknown_candidate') candidate%kind=99
            if(special_input=='nan_residual_tol') options%residual_tol=ieee_value(0._wp,ieee_quiet_nan)
            if(special_input=='positive_ftol') options%root%ftol=1._wp
            call refine_string(length_m,candidate,options,refined)
            call print_status(refined%status);call emit_bool('result_valid',refined%valid)
            if(.not.refined%valid) return
            call emit_real('root_per_m',refined%search%root);call emit_real('residual',refined%shot%residual)
            call emit_real('bracket_lo',refined%search%lo);call emit_real('bracket_hi',refined%search%hi)
            if(emit_profile) then
                call emit_int('profile_n',int(size(refined%shot%x_m),int64))
                do i=1,size(refined%shot%x_m)
                    write(key,'(A,I5.5)') 'x_',i;call emit_real(trim(key),refined%shot%x_m(i))
                    write(key,'(A,I5.5)') 'y_',i;call emit_real(trim(key),refined%shot%y_m(i))
                end do
            end if
            return
        end if
        if(case_id=='string_shot') then
            call shoot_string(length_m,k_per_m,options,shot,trim(special_input))
            call print_status(shot%status);call emit_bool('result_valid',shot%valid)
            call emit_int('accepted',int(shot%accepted,int64));call emit_int('rejected',int(shot%rejected,int64))
            call emit_int('rhs_calls',int(shot%evaluations,int64))
            if(.not.shot%valid) return
            n=size(shot%x_m)
            call emit_real('k_per_m',shot%k_per_m);call emit_real('boundary_y_m',shot%boundary_y_m)
            call emit_real('boundary_slope',shot%boundary_slope);call emit_real('residual',shot%residual)
            call emit_real('x_end_m',shot%x_m(n));call emit_text('normalization','unit_amplitude_then_physical_output')
            if(emit_profile) then
                call emit_int('profile_n',int(n,int64))
                do i=1,n
                    write(key,'(A,I5.5)') 'x_',i;call emit_real(trim(key),shot%x_m(i))
                    write(key,'(A,I5.5)') 'y_',i;call emit_real(trim(key),shot%y_m(i))
                    write(key,'(A,I5.5)') 'slope_',i;call emit_real(trim(key),shot%slope(i))
                end do
            end if
        else
            call linear_scan_grid(k_min,k_max,scan_n,grid,status)
            if(.not.is_ok(status)) then
                call print_status(status);call emit_bool('result_valid',.false.);return
            end if
            fail_index=0;if(special_input=='fail_scan_midpoint') fail_index=(scan_n+1)/2
            call scan_string(length_m,grid,options,spectrum,fail_index)
            call print_status(spectrum%status);call emit_bool('result_valid',spectrum%valid)
            call emit_int('roots_found',int(spectrum%roots_found,int64))
            call emit_int('unresolved_candidates',int(spectrum%unresolved_candidates,int64))
            if(spectrum%scan%valid) call print_scan(spectrum%scan)
            if(emit_profile .and. allocated(spectrum%scan%trials)) then
                call emit_int('trial_n',int(size(spectrum%scan%trials),int64))
                do i=1,size(spectrum%scan%trials)
                    write(key,'(A,I5.5)') 'trial_k_',i
                    call emit_real(trim(key),spectrum%scan%trials(i)%parameter)
                    write(key,'(A,I5.5)') 'trial_status_',i
                    call emit_text(trim(key),trim(status_name(spectrum%scan%trials(i)%status%code)))
                    if(.not.spectrum%scan%trials(i)%valid) cycle
                    write(key,'(A,I5.5)') 'trial_residual_',i
                    call emit_real(trim(key),spectrum%scan%trials(i)%residual)
                end do
            end if
            if(allocated(spectrum%refinements)) then
                do i=1,size(spectrum%refinements)
                    write(key,'(A,I5.5)') 'refinement_status_',i
                    call emit_text(trim(key),trim(status_name(spectrum%refinements(i)%status%code)))
                end do
            end if
            if(.not.spectrum%valid) return
            width=0._wp;resid=0._wp
            do i=1,spectrum%roots_found
                write(key,'(A,I0,A)') 'k',i,'_per_m'
                call emit_real(trim(key),spectrum%roots(i)%search%root)
                write(key,'(A,I0)') 'root_residual_',i
                call emit_real(trim(key),spectrum%roots(i)%shot%residual)
                width=max(width,spectrum%roots(i)%search%hi-spectrum%roots(i)%search%lo)
                resid=max(resid,abs(spectrum%roots(i)%shot%residual))
            end do
            call emit_real('max_bracket_width',width);call emit_real('max_root_residual',resid)
            call emit_text('labels','ordered roots inside requested interval; no stellar mode labels')
        end if
    end subroutine

    subroutine print_status(s)
        type(status_t),intent(in)::s
        call emit_text('status',trim(status_name(s%code)));call emit_text('category',trim(s%category))
        call emit_text('detail',trim(s%detail))
    end subroutine
    subroutine print_scan(r)
        type(scan_result_t),intent(in)::r
        integer::i
        character(len=48)::key
        call emit_int('barriers',int(r%barriers,int64));call emit_int('sign_brackets',int(r%sign_brackets,int64))
        call emit_int('sample_candidates',int(r%sample_candidates,int64))
        call emit_int('candidates',int(size(r%candidates),int64))
        call emit_bool('all_trials_succeeded',r%all_trials_succeeded)
        call emit_bool('spectrum_complete_claim',.false.)
        do i=1,size(r%candidates)
            write(key,'(A,I5.5)') 'candidate_lo_',i;call emit_real(trim(key),r%candidates(i)%lo)
            write(key,'(A,I5.5)') 'candidate_hi_',i;call emit_real(trim(key),r%candidates(i)%hi)
            write(key,'(A,I5.5)') 'candidate_kind_',i;call emit_int(trim(key),int(r%candidates(i)%kind,int64))
        end do
    end subroutine

    integer function count_crossed_failures(r) result(n)
        type(scan_result_t),intent(in)::r
        integer::i,j
        n=0
        do i=1,size(r%candidates)
            do j=r%candidates(i)%left_index,r%candidates(i)%right_index
                if(.not.is_ok(r%trials(j)%status) .or. .not.r%trials(j)%valid) then
                    n=n+1;exit
                end if
            end do
        end do
    end function

    subroutine probe(scenario,tol,r)
        character(len=*),intent(in)::scenario
        real(wp),intent(in)::tol
        type(scan_result_t),intent(out)::r
        type(trial_record_t),allocatable::t(:)
        integer::i,n
        n=3;if(scenario=='huge') n=2
        if(scenario=='empty') n=0
        allocate(t(n))
        do i=1,n
            t(i)%parameter=real(i,wp);t(i)%residual=real(2*i-4,wp)
            t(i)%valid=.true.;call clear_status(t(i)%status)
        end do
        select case(scenario)
        case('gap','failed_finite')
            call fail_status(t(2)%status,NS_OUT_OF_DOMAIN,'numerics','synthetic gap')
            if(scenario=='gap') t(2)%valid=.false.
        case('nan');t(2)%residual=ieee_value(0._wp,ieee_quiet_nan)
        case('inf');t(2)%residual=ieee_value(0._wp,ieee_positive_inf)
        case('invalid_record');t(2)%valid=.false.
        case('huge');t(1)%residual=-huge(1._wp);t(2)%residual=huge(1._wp)
        case('two_roots');t%residual=[-1._wp,1._wp,-1._wp]
        case('sample_zero')
        case('nearzero_covered');t%residual=[-1._wp,1e-12_wp,1._wp]
        case('sample_tangent');t%residual=[1._wp,0._wp,1._wp]
        case('zero_gap');t(1)%valid=.false.
        case('all_zero');t%residual=0._wp
        case('unsorted');t%parameter=[1._wp,3._wp,2._wp]
        case('duplicate');t(2)%parameter=t(1)%parameter
        case('nan_parameter');t(2)%parameter=ieee_value(0._wp,ieee_quiet_nan)
        case('orientation_change');t(2)%orientation=-1
        case('zero_orientation');t(2)%orientation=0
        case('empty')
        case default
            r=scan_result_t();call fail_status(r%status,NS_INVALID_INPUT,'numerics','unknown synthetic scan');return
        end select
        call build_brackets(t,tol,r)
    end subroutine
end module ns_shooting_cases
