module ns_matching_cases
    use ns_precision, only: wp
    use ns_status
    use ns_records, only: emit_real,emit_text,emit_bool,emit_int
    use ns_eos_types
    use ns_polytrope
    use ns_two_polytrope
    use ns_direct_join
    use ns_maxwell
    use ns_mu_bridge
    use ns_matching_types
    use ns_numerics_types, only: root_options_t
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_finite
    implicit none
    private
    public :: is_matching_case,run_matching_case
contains
    logical function is_matching_case(name)
        character(len=*),intent(in)::name
        is_matching_case=any(name==[character(len=20)::'two_poly_knot','direct_join','maxwell_match','mu_bridge'])
    end function

    subroutine run_matching_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token
        character(len=48)::action,scenario,special_input
        real(wp)::pressure_pa,mu_left_j,mu_right_j,mu_low_j,mu_high_j,mu_query_j
        real(wp)::root_xtol,root_rtol,residual_tol,gamma_crust,transition_n_m3,crust_n_min_m3
        real(wp)::core_n_anchor_m3,core_p_anchor_pa,core_gamma,core_c_j,core_n_min_m3,core_n_max_m3
        integer::max_iter,u,ios
        type(polytrope_t)::low,high,core
        type(two_polytrope_t)::two
        type(eos_state_t)::a,b,s,ba,bb,before,after
        type(direct_join_report_t)::direct
        type(transition_state_t)::transition
        type(mu_bridge_t)::bridge
        type(root_options_t)::opts
        type(status_t)::status,extra
        real(wp)::dn,den,mu_t,p_t,u0,nlow,nhigh,klow,khigh,clow,chigh
        logical::core_unchanged
        namelist /request/ request_token
        namelist /parameters/ action,scenario,pressure_pa,mu_left_j,mu_right_j,mu_low_j,mu_high_j,mu_query_j, &
            root_xtol,root_rtol,residual_tol,max_iter,gamma_crust,transition_n_m3,crust_n_min_m3, &
            core_n_anchor_m3,core_p_anchor_pa,core_gamma,core_c_j,core_n_min_m3,core_n_max_m3,special_input
        request_token='';action='evaluate';scenario='maxwell_exact';special_input='none'
        pressure_pa=8e32_wp;mu_left_j=1.46e-10_wp;mu_right_j=1.54e-10_wp
        mu_low_j=0.95_wp*1.5e-10_wp;mu_high_j=1.05_wp*1.5e-10_wp;mu_query_j=1.5e-10_wp
        root_xtol=1e-24_wp;root_rtol=0._wp;residual_tol=1e17_wp;max_iter=200
        gamma_crust=1.5_wp;transition_n_m3=1e44_wp;crust_n_min_m3=1e40_wp
        core_n_anchor_m3=1e44_wp;core_p_anchor_pa=1e32_wp;core_gamma=2.5_wp;core_c_j=1.5e-10_wp
        core_n_min_m3=1e44_wp;core_n_max_m3=5e44_wp
        call clear_eos_state(a);call clear_eos_state(b);call clear_eos_state(s)
        call clear_eos_state(ba);call clear_eos_state(bb)
        call emit_text('case_id',case_id)
        open(newunit=u,file=path,status='old',action='read',iostat=ios)
        if(ios==0) then
            read(u,nml=request,iostat=ios);if(ios==0) read(u,nml=parameters,iostat=ios);close(u)
        end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0) then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid T09 matching NAMELIST');call output_common();return
        end if
        select case(special_input)
        case('none')
        case('nan_mu_left');mu_left_j=ieee_value(0._wp,ieee_quiet_nan)
        case('nan_query');mu_query_j=ieee_value(0._wp,ieee_quiet_nan)
        case default
            call fail_status(status,NS_INVALID_INPUT,'data','unknown T09 special input');call output_common();return
        end select
        select case(trim(case_id))
        case('two_poly_knot')
            call init_polytrope(core,core_n_anchor_m3,core_p_anchor_pa,core_gamma,core_c_j, &
                                core_n_min_m3,core_n_max_m3,status)
            core_unchanged=.false.
            if(.not.is_ok(status)) then;call output_two();return;end if
            call eos_at_n(core,core_n_anchor_m3,before,status)
            if(.not.is_ok(status)) then;call output_two();return;end if
            call init_continuous_two_polytrope(two,core,gamma_crust,crust_n_min_m3,transition_n_m3,status)
            call eos_at_n(core,core_n_anchor_m3,after,extra)
            core_unchanged=is_ok(extra) .and. before%pressure_pa==after%pressure_pa .and. &
                           before%energy_j_m3==after%energy_j_m3 .and. &
                           before%mu_j==after%mu_j .and. before%chi==after%chi
            if(is_ok(status)) call two_polytrope_transition(two,a,b,status)
            call output_two()
        case('direct_join')
            call build_scenario(trim(scenario),low,high,status)
            if(.not.is_ok(status)) then;call output_direct();return;end if
            call evaluate_direct_join(low,high,pressure_pa,root_xtol,direct,status);call output_direct()
        case('maxwell_match')
            call build_scenario(trim(scenario),low,high,status)
            if(.not.is_ok(status)) then;call output_maxwell();return;end if
            opts=root_options_t(root_xtol,root_rtol,residual_tol,max_iter)
            call solve_maxwell_match(low,high,mu_left_j,mu_right_j,opts,transition,status);call output_maxwell()
        case('mu_bridge')
            if(trim(scenario)=='maxwell_exact') then
                call build_scenario('maxwell_exact',low,high,status)
                if(is_ok(status)) call construct_branch_bridge(low,high,mu_low_j,mu_high_j,bridge,status)
                if(is_ok(status)) then
                    call evaluate_mu_bridge(bridge,mu_query_j,s,dn,status)
                    call eos_at_mu(low,mu_low_j,ba,den,extra)
                    if(is_ok(extra)) call eos_at_mu(high,mu_high_j,bb,den,extra)
                end if
            else if(trim(scenario)=='negative_compressibility' .or. trim(scenario)=='endpoint_hidden_failure') then
                call construct_bridge_endpoints(1e-10_wp,1e32_wp,1e42_wp,2e-10_wp,2e32_wp,1.1e42_wp,bridge,status)
            else if(trim(scenario)=='acausal') then
                call construct_bridge_endpoints(1e-10_wp,1e32_wp,1e42_wp,2e-10_wp,2.05e32_wp,1.1e42_wp,bridge,status)
            else
                call fail_status(status,NS_INVALID_INPUT,'data','unknown bridge scenario')
            end if
            call output_bridge()
        case default
            call fail_status(status,NS_INVALID_INPUT,'data','unknown matching case');call output_common()
        end select
    contains
        subroutine build_scenario(name,lo,hi,st)
            character(len=*),intent(in)::name
            type(polytrope_t),intent(out)::lo,hi
            type(status_t),intent(out)::st
            type(status_t)::tmp
            u0=1.5e-11_wp;p_t=1e33_wp;mu_t=10._wp*u0
            nlow=p_t/u0;nhigh=2._wp*p_t/u0
            klow=p_t/(nlow*nlow);khigh=p_t/(nhigh*nhigh)
            clow=mu_t-2._wp*p_t/nlow;chigh=mu_t-2._wp*p_t/nhigh
            select case(name)
            case('maxwell_exact')
                call init_polytrope_k(lo,klow,2._wp,clow,nlow,1e43_wp,3e44_wp,st);if(.not.is_ok(st))return
                call init_polytrope_k(hi,khigh,2._wp,chigh,nhigh,1e43_wp,3e44_wp,st)
            case('maxwell_reversed')
                call init_polytrope_k(hi,klow,2._wp,clow,nlow,1e43_wp,3e44_wp,st);if(.not.is_ok(st))return
                call init_polytrope_k(lo,khigh,2._wp,chigh,nhigh,1e43_wp,3e44_wp,st)
            case('no_crossing')
                call init_polytrope_k(lo,klow,2._wp,clow,nlow,1e43_wp,3e44_wp,st);if(.not.is_ok(st))return
                call init_polytrope_k(hi,klow,2._wp,1.25e-10_wp,nlow,1e43_wp,3e44_wp,st)
            case default
                call fail_status(st,NS_INVALID_INPUT,'data','unknown phase scenario')
            end select
        end subroutine
        subroutine output_common()
            call emit_text('status',trim(status_name(status%code)))
            call emit_text('category',trim(status%category));call emit_text('detail',trim(status%detail))
        end subroutine
        subroutine output_two()
            call output_common();call emit_bool('model_ready',two_polytrope_ready(two))
            call emit_bool('core_unchanged',core_unchanged)
            if(is_ok(status) .and. a%valid .and. b%valid) then
                call emit_real('crust_pressure_pa',a%pressure_pa);call emit_real('core_pressure_pa',b%pressure_pa)
                call emit_real('crust_energy_j_m3',a%energy_j_m3);call emit_real('core_energy_j_m3',b%energy_j_m3)
                call emit_real('crust_mu_j',a%mu_j);call emit_real('core_mu_j',b%mu_j)
                call emit_real('crust_n_m3',a%n_b_m3);call emit_real('core_n_m3',b%n_b_m3)
                call emit_real('crust_chi',a%chi);call emit_real('core_chi',b%chi)
                call emit_bool('pressure_closed',a%pressure_pa==b%pressure_pa)
                call emit_bool('energy_closed',a%energy_j_m3==b%energy_j_m3)
                call emit_bool('chemical_closed',a%mu_j==b%mu_j)
                call emit_bool('density_closed',a%n_b_m3==b%n_b_m3)
                call emit_bool('chi_jump_nonzero',a%chi/=b%chi)
            end if
        end subroutine
        subroutine output_direct()
            call output_common();call emit_bool('report_ready',direct%ready)
            if(direct%ready) then
                call emit_real('pressure_pa',direct%pressure_pa)
                call emit_real('pressure_residual_pa',direct%pressure_residual_pa)
                call emit_real('mu_low_j',direct%mu_low_j);call emit_real('mu_high_j',direct%mu_high_j)
                call emit_real('mu_mismatch_j',direct%mu_mismatch_j)
                call emit_real('n_low_m3',direct%n_low_m3);call emit_real('n_high_m3',direct%n_high_m3)
                call emit_bool('pressure_closed',direct%pressure_closed)
                call emit_bool('chemical_closed',direct%chemical_closed)
                call emit_bool('thermodynamically_consistent',direct%thermodynamically_consistent)
            end if
        end subroutine
        subroutine output_maxwell()
            call output_common();call emit_bool('transition_valid',transition%valid)
            call emit_bool('branch_selection_valid',transition%branch_selection_valid)
            if(transition%mu_t_j>0._wp) then
                call emit_real('mu_t_j',transition%mu_t_j);call emit_real('p_t_pa',transition%p_t_pa)
                call emit_real('pressure_residual_pa',transition%pressure_residual_pa)
                call emit_real('mu_mismatch_j',transition%mu_mismatch_j)
                call emit_real('n_low_m3',transition%n_low_m3);call emit_real('n_high_m3',transition%n_high_m3)
                call emit_real('delta_n_m3',transition%delta_n_m3)
                call emit_real('energy_low_j_m3',transition%energy_low_j_m3)
                call emit_real('energy_high_j_m3',transition%energy_high_j_m3)
                call emit_real('delta_e_j_m3',transition%delta_e_j_m3)
                call emit_int('root_iterations',int(transition%root_iterations,int64))
                call emit_bool('thermodynamically_consistent',transition%thermodynamically_consistent)
            end if
        end subroutine
        subroutine output_bridge()
            type(eos_state_t)::el,eh
            real(wp)::d0,d1
            call output_common();call emit_bool('bridge_ready',bridge%ready)
            call emit_bool('interval_extrema_checked',bridge%interval_extrema_checked)
            call emit_bool('stable_interval',bridge%stable_interval)
            call emit_bool('causal_interval',bridge%causal_interval)
            call emit_bool('admitted',bridge%ready .and. bridge%stable_interval .and. bridge%causal_interval)
            if(ba%valid .and. bb%valid) then
                call emit_bool('endpoint_pressure_closed', &
                    abs(bridge%p_low_pa-ba%pressure_pa)<=5e-12_wp*max(abs(ba%pressure_pa),tiny(1._wp)) .and. &
                    abs(bridge%p_high_pa-bb%pressure_pa)<=5e-12_wp*max(abs(bb%pressure_pa),tiny(1._wp)))
                call emit_bool('endpoint_density_closed', &
                    abs(bridge%n_low_m3-ba%n_b_m3)<=5e-12_wp*max(abs(ba%n_b_m3),tiny(1._wp)) .and. &
                    abs(bridge%n_high_m3-bb%n_b_m3)<=5e-12_wp*max(abs(bb%n_b_m3),tiny(1._wp)))
            end if
            if(bridge%mu_low_j>0._wp) then
                call emit_real('bridge_low_pressure_pa',bridge%p_low_pa)
                call emit_real('bridge_high_pressure_pa',bridge%p_high_pa)
                call emit_real('bridge_low_density_m3',bridge%n_low_m3)
                call emit_real('bridge_high_density_m3',bridge%n_high_m3)
            end if
            if(bridge%interval_extrema_checked) then
                call emit_real('min_d2p',bridge%min_d2p);call emit_real('min_causal_margin',bridge%min_causal_margin)
            end if
            if(is_ok(status) .and. s%valid) then
                call emit_real('mu_query_j',mu_query_j);call emit_real('pressure_pa',s%pressure_pa)
                call emit_real('n_b_m3',s%n_b_m3)
                call emit_real('energy_j_m3',s%energy_j_m3);call emit_real('chi',s%chi)
                call emit_real('dn_dmu_m3_per_j',dn)
                call evaluate_mu_bridge(bridge,bridge%mu_low_j,el,d0,extra)
                if(is_ok(extra)) then
                    call emit_real('low_pressure_error_pa',el%pressure_pa-bridge%p_low_pa)
                    call emit_real('low_density_relative_error',(el%n_b_m3-bridge%n_low_m3)/bridge%n_low_m3)
                end if
                call evaluate_mu_bridge(bridge,bridge%mu_high_j,eh,d1,extra)
                if(is_ok(extra)) then
                    call emit_real('high_pressure_error_pa',eh%pressure_pa-bridge%p_high_pa)
                    call emit_real('high_density_relative_error',(eh%n_b_m3-bridge%n_high_m3)/bridge%n_high_m3)
                end if
            end if
        end subroutine
    end subroutine
end module ns_matching_cases
