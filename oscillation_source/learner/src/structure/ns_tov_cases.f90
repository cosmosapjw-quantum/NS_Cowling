! Public I/O adapter, invokes only the selected workspace's structural kernels.
module ns_tov_cases
    use ns_precision, only: wp
    use ns_physical_constants, only: c_si,g_si,pi,nominal_solar_mass_kg
    use ns_status
    use ns_records
    use ns_eos_types
    use ns_polytrope
    use ns_background_records
    use ns_tov
    use ns_target_mass
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_finite
    implicit none
    private
    public :: is_tov_case,run_tov_case
contains
    logical function is_tov_case(name)
        character(len=*),intent(in)::name
        is_tov_case=name=='tov_baseline'.or.name=='tov_rhs'.or.name=='target_mass'
    end function
    subroutine run_tov_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token
        character(len=64)::eos_kind,special_input,surface_policy,branch_id
        real(wp)::n_anchor_m3,p_anchor_pa,gamma,c_j,n_min_m3,n_max_m3,pc_pa,p_cut_pa,r0_m,r_limit_m
        real(wp)::rtol,atol_r_m,atol_mass_kg,atol_nu,dq_initial,dq_min,dq_max,nu_start
        real(wp)::radius_m,mass_kg,pressure_pa,energy_uniform,pc_lo_pa,pc_hi_pa,m_target_kg
        real(wp)::mass_atol_kg,mass_rtol,root_xtol,stop1_pa,stop2_pa
        integer::max_steps,scan_n,root_max_iter,ios,u,n,i
        logical::emit_profile,use_stops
        type(polytrope_t)::model
        type(eos_state_t)::lower,upper,state
        type(star_options_t)::options
        type(mass_options_t)::mass_options
        type(star_result_t)::mass_result
        type(star_profile_t)::profile
        type(status_t)::status
        real(wp)::dydr(3),y(3),stops(2),mcheck,redshift
        character(len=32)::suffix
        namelist /request/ request_token
        namelist /parameters/ n_anchor_m3,p_anchor_pa,gamma,c_j,n_min_m3,n_max_m3,pc_pa,p_cut_pa,r0_m,r_limit_m, &
            rtol,atol_r_m,atol_mass_kg,atol_nu,dq_initial,dq_min,dq_max,nu_start,max_steps,emit_profile, &
            eos_kind,special_input,surface_policy,radius_m,mass_kg,pressure_pa,energy_uniform, &
            pc_lo_pa,pc_hi_pa,m_target_kg,mass_atol_kg,mass_rtol,root_xtol,scan_n,root_max_iter,branch_id, &
            use_stops,stop1_pa,stop2_pa
        request_token='';eos_kind='polytrope';special_input='none'
        surface_policy='finite_cutoff_schwarzschild_time';branch_id='sampled_increasing_pc'
        n_anchor_m3=1e44_wp;p_anchor_pa=1e32_wp;gamma=2._wp;c_j=1.5e-10_wp
        n_min_m3=1e34_wp;n_max_m3=1e46_wp;pc_pa=1e34_wp;p_cut_pa=1e22_wp
        r0_m=1._wp;r_limit_m=1e6_wp;nu_start=0._wp
        rtol=1e-8_wp;atol_r_m=1e-6_wp;atol_mass_kg=1e14_wp;atol_nu=1e-12_wp
        dq_initial=1e-5_wp;dq_min=1e-16_wp;dq_max=0.2_wp;max_steps=20000
        radius_m=1e3_wp;mass_kg=1e27_wp;pressure_pa=1e34_wp;energy_uniform=1e35_wp
        pc_lo_pa=3e32_wp;pc_hi_pa=3e34_wp;m_target_kg=1.2e30_wp
        mass_atol_kg=1e20_wp;mass_rtol=1e-8_wp;root_xtol=1e-12_wp;scan_n=9;root_max_iter=70
        emit_profile=.false.;use_stops=.false.;stop1_pa=1e31_wp;stop2_pa=1e28_wp
        call emit_text('case_id',case_id)
        open(newunit=u,file=path,status='old',action='read',iostat=ios)
        if(ios==0) then
            read(u,nml=request,iostat=ios)
            if(ios==0) read(u,nml=parameters,iostat=ios)
            close(u)
        end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0) then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid structural NAMELIST');call finish();return
        end if
        select case(special_input)
        case('none','fail_callback','wrong_pressure','invalid_state','nan_energy','lifecycle','sharp_event')
        case('nan_pc');pc_pa=ieee_value(0._wp,ieee_quiet_nan)
        case('nan_r');radius_m=ieee_value(0._wp,ieee_quiet_nan)
        case default
            call fail_status(status,NS_INVALID_INPUT,'data','unknown structural special input');call finish();return
        end select
        if(eos_kind=='polytrope') then
            call init_polytrope(model,n_anchor_m3,p_anchor_pa,gamma,c_j,n_min_m3,n_max_m3,status)
            if(.not.is_ok(status)) then
                call finish();return
            end if
        else if(eos_kind=='uniform_benchmark') then
            if(.not.ieee_is_finite(energy_uniform).or.energy_uniform<=0._wp) then
                call fail_status(status,NS_INVALID_INPUT,'physics', &
            'uniform benchmark epsilon must be positive');call finish();return
            end if
        else
            call fail_status(status,NS_INVALID_INPUT,'physics','unknown EOS kind');call finish();return
        end if
        options%r0_m=r0_m;options%r_limit_m=r_limit_m;options%p_cut_pa=p_cut_pa
        options%nu_start=nu_start;options%surface_policy=surface_policy
        options%eos_id=eos_kind
        options%ode%atol=[atol_r_m,atol_mass_kg,atol_nu];options%ode%rtol=rtol
        options%ode%h_initial=dq_initial;options%ode%h_min=dq_min;options%ode%h_max=dq_max
        options%ode%max_steps=max_steps
        if(special_input=='sharp_event') then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','sharp interfaces are not enabled in T05');call finish();return
        end if
        if(case_id=='tov_rhs') then
            y=[mass_kg,pressure_pa,nu_start]
            call tov_rhs(eos,radius_m,y,dydr,status)
            if(is_ok(status)) then
                call emit_real('dm_dr_kg_m',dydr(1));call emit_real('dp_dr_pa_m',dydr(2))
                call emit_real('dnu_dr_per_m',dydr(3))
            end if
        else
            if(case_id=='target_mass') then
                mass_options%scan_n=scan_n;mass_options%mass_atol_kg=mass_atol_kg;mass_options%mass_rtol=mass_rtol
                mass_options%root%xtol=root_xtol;mass_options%root%max_iter=root_max_iter
                mass_options%branch_id=branch_id
                if(use_stops) then
                    call solve_target_mass(eos,m_target_kg,[pc_lo_pa,pc_hi_pa],options,mass_options,mass_result, &
                                           [stop1_pa,stop2_pa])
                else
                    call solve_target_mass(eos,m_target_kg,[pc_lo_pa,pc_hi_pa],options,mass_options,mass_result)
                end if
                profile=mass_result%profile;status=mass_result%status
                call emit_bool('target_converged',mass_result%valid)
                call emit_bool('sampled_branch_increasing',mass_result%scan_monotone)
                call emit_int('star_evaluations',int(mass_result%star_evaluations,int64))
                if(mass_result%valid) then
                    call emit_real('mass_residual_kg',mass_result%mass_residual_kg)
                    call emit_real('target_mass_kg',m_target_kg)
                end if
            else
                if(use_stops) then
                    call integrate_tov(eos,pc_pa,options,profile,[stop1_pa,stop2_pa])
                else
                    call integrate_tov(eos,pc_pa,options,profile)
                end if
                if(special_input=='lifecycle'.and.profile%valid) call integrate_tov(eos,-1._wp,options,profile)
                status=profile%status
            end if
            if(profile%valid.and.is_ok(status)) call output_profile()
        end if
        call finish()
    contains
        subroutine eos(p,side,s,st)
            real(wp),intent(in)::p
            integer,intent(in)::side
            type(eos_state_t),intent(out)::s
            type(status_t),intent(out)::st
            if(special_input=='fail_callback') then
                call clear_eos_state(s);call fail_status(st,NS_OUT_OF_DOMAIN,'physics','deliberate EOS rejection');return
            end if
            if(eos_kind=='polytrope') then
                call eos_at_pressure(model,p,side,s,st)
            else
                call clear_eos_state(s)
                call fail_status(st,NS_OUT_OF_DOMAIN,'physics','formal uniform benchmark requires P>0')
                if(p<=0._wp.or..not.ieee_is_finite(p)) return
                ! Formal incompressible benchmark ONLY. chi=1 is a finite unused placeholder;
                ! no compressional response or Cowling admission is claimed for this test EOS.
                s=eos_state_t(p,energy_uniform,energy_uniform/c_j,(energy_uniform+p)/(energy_uniform/c_j), &
                              1._wp,side,.true.)
                call clear_status(st)
            end if
            if(.not.is_ok(st)) return
            if(special_input=='wrong_pressure') s%pressure_pa=2._wp*p
            if(special_input=='invalid_state') s%valid=.false.
            if(special_input=='nan_energy') s%energy_j_m3=ieee_value(0._wp,ieee_quiet_nan)
        end subroutine
        subroutine finish()
            call emit_text('status',trim(status_name(status%code)))
            call emit_text('category',trim(status%category));call emit_text('detail',trim(status%detail))
            call emit_bool('result_valid',profile%valid.and.is_ok(status))
        end subroutine
        subroutine output_profile()
            n=size(profile%r_m)
            call emit_real('radius_m',profile%radius_m);call emit_real('mass_kg',profile%mass_kg)
            call emit_real('mass_nominal_solar',profile%mass_kg/nominal_solar_mass_kg)
            call emit_real('pc_pa',profile%pc_pa);call emit_real('p_start_pa',profile%p_start_pa)
            call emit_real('p_cut_pa',profile%p_cut_pa);call emit_real('surface_pressure_pa',profile%p_pa(n))
            call emit_real('surface_residual_pa',profile%surface_residual_pa)
            call emit_real('lapse_residual',profile%lapse_residual);call emit_real('lapse_shift',profile%lapse_shift)
            call emit_real('nu_surface',profile%nu(n));call emit_real('r0_m',profile%r0_m)
            call emit_text('surface_policy',profile%surface_policy);call emit_text('coordinate',profile%coordinate)
            call emit_text('eos_kind',eos_kind);call emit_text('branch_id',profile%branch_id)
            call emit_bool('exact_vacuum_surface',.false.)
            call emit_bool('cowling_background_admitted',eos_kind=='polytrope')
            call emit_int('profile_n',int(n,int64));call emit_int('accepted',int(profile%accepted,int64))
            call emit_int('rejected',int(profile%rejected,int64));call emit_int('stops_hit',int(profile%stops_hit,int64))
            call emit_bool('mass_non_decreasing',all(profile%m_kg(2:)>=profile%m_kg(:n-1)))
            call emit_bool('pressure_non_increasing',all(profile%p_pa(2:)<=profile%p_pa(:n-1)))
            call emit_bool('surface_event_closed',profile%surface_residual_pa==0._wp)
            call emit_bool('lapse_contract_closed',abs(profile%lapse_residual)<128._wp*epsilon(1._wp))
            redshift=maxval(abs((profile%mu_j/profile%mu_j(1))*exp(profile%nu-profile%nu(1))-1._wp))
            call emit_real('redshifted_mu_rel_spread',redshift)
            if(emit_profile) then
                do i=1,n
                    write(suffix,'(I5.5)') i
                    call emit_real('q_'//trim(suffix),profile%q(i))
                    call emit_real('r_'//trim(suffix),profile%r_m(i));call emit_real('m_'//trim(suffix),profile%m_kg(i))
                    call emit_real('p_'//trim(suffix),profile%p_pa(i));call emit_real('eps_'//trim(suffix),profile%energy_j_m3(i))
                    call emit_real('n_'//trim(suffix),profile%n_b_m3(i));call emit_real('mu_'//trim(suffix),profile%mu_j(i))
                    call emit_real('chi_'//trim(suffix),profile%chi(i));call emit_real('nu_'//trim(suffix),profile%nu(i))
                    call emit_real('lambda_'//trim(suffix), &
            profile%lambda(i));call emit_real('dnu_'//trim(suffix),profile%dnu_dr(i))
                end do
            end if
        end subroutine
    end subroutine
end module
