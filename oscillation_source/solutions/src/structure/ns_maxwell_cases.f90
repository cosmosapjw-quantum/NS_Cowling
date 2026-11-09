! CLI adapter for separate high-/low-density normalized EOS tables.
! The table providers remain distinct across the exact Maxwell pressure.
module ns_maxwell_cases
    use, intrinsic :: iso_fortran_env, only: int64
    use ns_precision, only: wp
    use ns_status
    use ns_physical_constants, only: nominal_solar_mass_kg,two_pi
    use ns_records, only: emit_text,emit_real,emit_int,emit_bool
    use ns_eos_types, only: eos_state_t
    use ns_table_io, only: normalized_table_t,load_normalized_table
    use ns_table_eos, only: table_eos_t,table_audit_t,init_table_eos,eos_table_at_pressure
    use ns_background_records, only: star_options_t
    use ns_maxwell_tov, only: maxwell_star_t
    use ns_production_solvers, only: solve_tov_star,solve_cowling_spectrum
    use ns_mode_spectrum, only: mode_options_t
    use ns_maxwell_mode_spectrum, only: maxwell_mode_spectrum_t
    use ns_scan, only: linear_scan_grid
    use ns_cowling, only: cowling_options_t
    use ns_maxwell_cowling, only: maxwell_cowling_shot_t,shoot_maxwell_cowling
    implicit none
    private
    public :: is_maxwell_case,run_maxwell_case
contains
    logical function is_maxwell_case(name)
        character(len=*),intent(in)::name
        is_maxwell_case=name=='maxwell_table_star'.or.name=='maxwell_table_cowling'.or. &
            name=='maxwell_table_spectrum'
    end function

    subroutine run_maxwell_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token,conversion
        character(len=512)::core_table_path,envelope_table_path
        real(wp)::pc_pa,pt_pa,p_cut_pa,r0_m,r_limit_m,rtol,atol_r_m,atol_mass_kg,atol_nu
        real(wp)::dq_initial,dq_min,dq_max,nu_start,frequency_hz,atol_w,atol_v
        real(wp)::frequency_min_hz,frequency_max_hz,root_xtol,root_rtol,residual_tol,zero_tol
        real(wp),allocatable::frequency_grid(:)
        integer::max_steps,u,ios,scan_n,max_iter,i
        character(len=32)::suffix
        type(normalized_table_t)::core_raw,envelope_raw
        type(table_eos_t)::core_eos,envelope_eos
        type(table_audit_t)::core_audit,envelope_audit
        type(star_options_t)::star_options
        type(maxwell_star_t)::star
        type(cowling_options_t)::cowling_options
        type(maxwell_cowling_shot_t)::shot
        type(mode_options_t)::mode_options
        type(maxwell_mode_spectrum_t)::spectrum
        type(eos_state_t)::inner,outer
        type(status_t)::status

        namelist /request/ request_token
        namelist /parameters/ core_table_path,envelope_table_path,pc_pa,pt_pa,p_cut_pa,r0_m,r_limit_m, &
            rtol,atol_r_m,atol_mass_kg,atol_nu,dq_initial,dq_min,dq_max,nu_start,max_steps, &
            frequency_hz,conversion,atol_w,atol_v,frequency_min_hz,frequency_max_hz,scan_n, &
            root_xtol,root_rtol,residual_tol,zero_tol,max_iter

        request_token='';core_table_path='';envelope_table_path=''
        pc_pa=1.4e32_wp;pt_pa=1e32_wp;p_cut_pa=1e30_wp
        r0_m=1._wp;r_limit_m=1e6_wp;rtol=1e-9_wp
        atol_r_m=1e-5_wp;atol_mass_kg=1e14_wp;atol_nu=1e-12_wp
        dq_initial=1e-5_wp;dq_min=1e-15_wp;dq_max=0.05_wp;nu_start=0._wp;max_steps=100000
        frequency_hz=0._wp;conversion='slow';atol_w=1e-12_wp;atol_v=1e-12_wp
        frequency_min_hz=1799._wp;frequency_max_hz=1801._wp;scan_n=3
        root_xtol=1e-9_wp;root_rtol=0._wp;residual_tol=1e-5_wp;zero_tol=1e-10_wp;max_iter=100
        call emit_text('case_id',case_id)
        open(newunit=u,file=path,status='old',action='read',iostat=ios)
        if(ios==0) then
            read(u,nml=request,iostat=ios)
            if(ios==0) read(u,nml=parameters,iostat=ios)
            close(u)
        end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0) then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid native Maxwell NAMELIST')
            call finish();return
        end if
        call load_normalized_table(trim(core_table_path),core_raw,status)
        if(.not.is_ok(status)) then;call finish();return;end if
        call load_normalized_table(trim(envelope_table_path),envelope_raw,status)
        if(.not.is_ok(status)) then;call finish();return;end if
        call init_table_eos(core_raw,core_eos,core_audit,status)
        if(.not.is_ok(status).or..not.core_audit%fully_admitted) then
            if(is_ok(status)) call fail_status(status,NS_OUT_OF_DOMAIN,'physics','high-density EOS table not admitted')
            call finish();return
        end if
        call init_table_eos(envelope_raw,envelope_eos,envelope_audit,status)
        if(.not.is_ok(status).or..not.envelope_audit%fully_admitted) then
            if(is_ok(status)) call fail_status(status,NS_OUT_OF_DOMAIN,'physics','low-density EOS table not admitted')
            call finish();return
        end if
        star_options%r0_m=r0_m;star_options%r_limit_m=r_limit_m
        star_options%p_cut_pa=p_cut_pa;star_options%nu_start=nu_start
        star_options%surface_policy='finite_cutoff_schwarzschild_time'
        star_options%eos_id='two_table_sharp_maxwell'
        star_options%ode%atol=[atol_r_m,atol_mass_kg,atol_nu]
        star_options%ode%rtol=rtol;star_options%ode%h_initial=dq_initial
        star_options%ode%h_min=dq_min;star_options%ode%h_max=dq_max;star_options%ode%max_steps=max_steps
        call solve_tov_star(high_eos,low_eos,pc_pa,pt_pa,star_options,star)
        status=star%status
        call print_status(status)
        call emit_bool('star_valid',star%valid)
        if(.not.star%valid.or..not.is_ok(status)) return
        call emit_real('mass_kg',star%mass_kg)
        call emit_real('mass_nominal_solar',star%mass_kg/nominal_solar_mass_kg)
        call emit_real('radius_m',star%radius_m)
        call emit_real('central_pressure_pa',pc_pa)
        call emit_real('transition_pressure_pa',pt_pa)
        call emit_real('surface_cutoff_pressure_pa',p_cut_pa)
        call emit_real('interface_radius_m',star%jump%radius_m)
        call emit_real('interface_mass_kg',star%jump%mass_kg)
        call emit_real('energy_jump_j_m3',star%jump%delta_energy_j_m3)
        call emit_real('baryon_density_jump_m3',star%jump%delta_density_m3)
        call emit_real('chemical_potential_j',star%jump%chemical_potential_j)
        call emit_real('chemical_potential_mismatch_j',star%jump%chemical_potential_mismatch_j)
        call emit_real('thermodynamic_jump_residual_j_m3',star%jump%thermodynamic_residual_j_m3)
        call emit_real('seidov_ratio',star%jump%seidov_ratio)
        call emit_real('seidov_critical',star%jump%seidov_critical)
        call emit_bool('rapid_onset_stable_diagnostic',star%jump%rapid_onset_stable_diagnostic)
        call emit_int('accepted_steps',int(star%accepted,int64))
        call emit_int('rejected_steps',int(star%rejected,int64))
        call emit_int('rhs_evaluations',int(star%evaluations,int64))
        call emit_bool('result_valid',.true.)
        if(case_id=='maxwell_table_spectrum') then
            mode_options%shot%ell=2
            mode_options%shot%amplitude=1._wp
            mode_options%shot%nu_shift=0._wp
            mode_options%shot%ode%atol=[atol_w,atol_v]
            mode_options%shot%ode%rtol=rtol
            mode_options%shot%ode%h_initial=dq_initial
            mode_options%shot%ode%h_min=dq_min
            mode_options%shot%ode%h_max=dq_max
            mode_options%shot%ode%max_steps=max_steps
            mode_options%root%xtol=root_xtol
            mode_options%root%rtol=root_rtol
            mode_options%root%ftol=0._wp
            mode_options%root%max_iter=max_iter
            mode_options%residual_tol=residual_tol
            mode_options%zero_tol=zero_tol
            call linear_scan_grid(two_pi*frequency_min_hz,two_pi*frequency_max_hz,scan_n,frequency_grid,status)
            if(.not.is_ok(status)) then
                call print_status(status)
                call emit_bool('spectrum_scan_valid',.false.)
                call emit_bool('spectrum_valid',.false.)
                return
            end if
            call solve_cowling_spectrum(star,frequency_grid,mode_options,trim(conversion),spectrum)
            status=spectrum%status
            call print_status(status)
            call emit_bool('spectrum_scan_valid',spectrum%scan%valid)
            call emit_bool('spectrum_all_trials_succeeded',spectrum%scan%all_trials_succeeded)
            call emit_bool('spectrum_valid',spectrum%valid)
            call emit_int('spectrum_scan_n',int(size(frequency_grid),int64))
            call emit_int('spectrum_barriers',int(spectrum%scan%barriers,int64))
            call emit_int('spectrum_sign_brackets',int(spectrum%scan%sign_brackets,int64))
            call emit_int('spectrum_sample_candidates',int(spectrum%scan%sample_candidates,int64))
            call emit_int('spectrum_roots_found',int(spectrum%roots_found,int64))
            call emit_int('spectrum_unresolved_candidates',int(spectrum%unresolved_candidates,int64))
            call emit_text('phase_conversion',trim(conversion))
            if(allocated(spectrum%scan%trials)) then
                do i=1,size(spectrum%scan%trials)
                    write(suffix,'(I5.5)')i
                    call emit_real('spectrum_trial_hz_'//trim(suffix), &
                        spectrum%scan%trials(i)%parameter/two_pi)
                    call emit_text('spectrum_trial_status_'//trim(suffix), &
                        trim(status_name(spectrum%scan%trials(i)%status%code)))
                    if(spectrum%scan%trials(i)%valid) call emit_real('spectrum_trial_residual_'//trim(suffix), &
                        spectrum%scan%trials(i)%residual)
                end do
            end if
            if(allocated(spectrum%refinements)) then
                do i=1,size(spectrum%refinements)
                    write(suffix,'(I5.5)')i
                    call emit_bool('spectrum_root_valid_'//trim(suffix),spectrum%refinements(i)%root_valid)
                    call emit_text('spectrum_root_status_'//trim(suffix), &
                        trim(status_name(spectrum%refinements(i)%root_status%code)))
                    call emit_text('spectrum_identity_'//trim(suffix),'unresolved')
                    if(spectrum%refinements(i)%search%valid) then
                        call emit_real('spectrum_root_hz_'//trim(suffix), &
                            spectrum%refinements(i)%search%root/two_pi)
                        call emit_real('spectrum_root_residual_'//trim(suffix), &
                            spectrum%refinements(i)%shot%residual)
                    end if
                end do
            end if
        else if(case_id=='maxwell_table_cowling'.and.frequency_hz>0._wp) then
            cowling_options%ell=2;cowling_options%amplitude=1._wp;cowling_options%nu_shift=0._wp
            cowling_options%ode%atol=[atol_w,atol_v]
            cowling_options%ode%rtol=rtol;cowling_options%ode%h_initial=dq_initial
            cowling_options%ode%h_min=dq_min;cowling_options%ode%h_max=dq_max
            cowling_options%ode%max_steps=max_steps
            call shoot_maxwell_cowling(star,two_pi*frequency_hz,trim(conversion),cowling_options,shot)
            call print_status(shot%status)
            call emit_bool('cowling_shot_valid',shot%valid)
            if(shot%valid.and.is_ok(shot%status)) then
                call emit_real('frequency_hz',frequency_hz)
                call emit_real('surface_residual',shot%residual)
                call emit_real('surface_q_over_r2',shot%q_over_r2)
                call emit_real('interface_delta_p_inner_pa',shot%delta_p_inner_pa)
                call emit_real('interface_delta_p_outer_pa',shot%delta_p_outer_pa)
                call emit_real('interface_w_inner_m3',shot%w_inner_m3)
                call emit_real('interface_v_inner_m2',shot%v_inner_m2)
                call emit_real('interface_w_outer_m3',shot%w_outer_m3)
                call emit_real('interface_v_outer_m2',shot%v_outer_m2)
                call emit_int('cowling_accepted_steps',int(shot%accepted,int64))
                call emit_int('cowling_rejected_steps',int(shot%rejected,int64))
                call emit_int('cowling_rhs_evaluations',int(shot%evaluations,int64))
                call emit_int('cowling_profile_points',int(size(shot%q),int64))
                call emit_text('phase_conversion',trim(conversion))
            end if
        else
            call emit_bool('cowling_shot_requested',.false.)
        end if
    contains
        subroutine high_eos(p,side,state,st)
            real(wp),intent(in)::p
            integer,intent(in)::side
            type(eos_state_t),intent(out)::state
            type(status_t),intent(out)::st
            call eos_table_at_pressure(core_eos,p,side,state,st)
        end subroutine
        subroutine low_eos(p,side,state,st)
            real(wp),intent(in)::p
            integer,intent(in)::side
            type(eos_state_t),intent(out)::state
            type(status_t),intent(out)::st
            call eos_table_at_pressure(envelope_eos,p,side,state,st)
        end subroutine
        subroutine print_status(st)
            type(status_t),intent(in)::st
            call emit_text('status',trim(status_name(st%code)))
            call emit_text('category',trim(st%category))
            call emit_text('detail',trim(st%detail))
        end subroutine
        subroutine finish()
            if(.not.is_ok(status)) call print_status(status)
            call emit_bool('result_valid',.false.)
        end subroutine
    end subroutine
end module ns_maxwell_cases
