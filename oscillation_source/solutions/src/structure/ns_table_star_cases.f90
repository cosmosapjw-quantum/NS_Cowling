! T14 adapter: an admitted normalized table is exposed to the existing SI TOV solver.
! No TOV equation, integrator, or legacy polytrope path is duplicated here.
module ns_table_star_cases
    use ns_precision, only: wp
    use ns_physical_constants, only: nominal_solar_mass_kg
    use ns_status
    use ns_records, only: emit_real,emit_text,emit_bool,emit_int
    use ns_eos_types, only: eos_state_t,SIDE_REGULAR
    use ns_table_io, only: normalized_table_t,load_normalized_table
    use ns_table_eos, only: table_eos_t,table_audit_t,init_table_eos,eos_table_at_pressure,TABLE_INTERPOLATION_ID
    use ns_background_records, only: star_options_t,star_profile_t
    use ns_production_solvers, only: solve_tov_star
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private
    public :: is_table_star_case,run_table_star_case
contains
    logical function is_table_star_case(name)
        character(len=*),intent(in)::name
        is_table_star_case=name=='real_table_star'.or.name=='real_table_star_split'
    end function is_table_star_case

    subroutine run_table_star_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token,suffix
        character(len=512)::table_path
        real(wp)::pc_pa,p_cut_pa,r0_m,r_limit_m,rtol,atol_r_m,atol_mass_kg,atol_nu,split_pressure_pa
        real(wp)::dq_initial,dq_min,dq_max,nu_start,redshift
        integer::max_steps,u,ios,i,n,split_inner_calls,split_outer_calls
        logical::emit_profile,table_admitted
        type(normalized_table_t)::raw
        type(table_eos_t)::model
        type(table_audit_t)::audit
        type(star_options_t)::options
        type(star_profile_t)::profile
        type(status_t)::status
        namelist /request/ request_token
        namelist /parameters/ table_path,pc_pa,p_cut_pa,r0_m,r_limit_m,rtol,atol_r_m,atol_mass_kg,atol_nu, &
            dq_initial,dq_min,dq_max,nu_start,max_steps,emit_profile,split_pressure_pa

        request_token='';table_path='';pc_pa=1.5e34_wp;p_cut_pa=-1._wp
        r0_m=1._wp;r_limit_m=1e6_wp;rtol=1e-8_wp;atol_r_m=1e-6_wp;atol_mass_kg=1e14_wp;atol_nu=1e-12_wp
        dq_initial=1e-5_wp;dq_min=1e-16_wp;dq_max=0.2_wp;nu_start=0._wp;max_steps=30000
        emit_profile=.false.;table_admitted=.false.;profile=star_profile_t();split_pressure_pa=-1._wp
        split_inner_calls=0;split_outer_calls=0
        call emit_text('case_id',case_id)
        open(newunit=u,file=path,status='old',action='read',iostat=ios)
        if(ios==0) then
            read(u,nml=request,iostat=ios);if(ios==0) read(u,nml=parameters,iostat=ios);close(u)
        end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0) then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid T14 table-star NAMELIST');call finish();return
        end if
        call load_normalized_table(trim(table_path),raw,status)
        if(.not.is_ok(status)) then;call finish();return;end if
        call init_table_eos(raw,model,audit,status)
        table_admitted=is_ok(status).and.audit%fully_admitted
        if(.not.table_admitted) then;call finish();return;end if
        if(p_cut_pa<=0._wp) p_cut_pa=raw%pressure_pa(1)
        if(case_id=='real_table_star_split') then
            if(split_pressure_pa<=p_cut_pa.or.split_pressure_pa>=pc_pa) then
                call fail_status(status,NS_OUT_OF_DOMAIN,'physics','split pressure must lie between cutoff and central pressure')
                call finish();return
            end if
        end if
        options%r0_m=r0_m;options%r_limit_m=r_limit_m;options%p_cut_pa=p_cut_pa;options%nu_start=nu_start
        options%surface_policy='finite_cutoff_schwarzschild_time';options%eos_id='direct_P_eps_PCHIP_logP_firstlaw'
        options%ode%atol=[atol_r_m,atol_mass_kg,atol_nu];options%ode%rtol=rtol
        options%ode%h_initial=dq_initial;options%ode%h_min=dq_min;options%ode%h_max=dq_max
        options%ode%max_steps=max_steps
        call solve_tov_star(eos,pc_pa,options,profile)
        status=profile%status
        if(profile%valid.and.is_ok(status)) call output_profile()
        call finish()
    contains
        subroutine eos(p,side,state,st)
            real(wp),intent(in)::p
            integer,intent(in)::side
            type(eos_state_t),intent(out)::state
            type(status_t),intent(out)::st
            if(case_id=='real_table_star_split') then
                if(p>=split_pressure_pa) then
                    split_inner_calls=split_inner_calls+1
                else
                    split_outer_calls=split_outer_calls+1
                end if
            end if
            call eos_table_at_pressure(model,p,side,state,st)
        end subroutine eos

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
            call emit_text('eos_kind','normalized_table');call emit_text('eos_id',trim(profile%eos_id))
            call emit_bool('exact_vacuum_surface',.false.)
            if(case_id=='real_table_star_split') then
                call emit_real('split_pressure_pa',split_pressure_pa)
                call emit_int('split_inner_calls',int(split_inner_calls,int64))
                call emit_int('split_outer_calls',int(split_outer_calls,int64))
                call emit_bool('split_dispatch_same_provider',.true.)
            end if
            call emit_bool('cowling_background_admitted',audit%fully_admitted .and. &
                all(profile%chi>0._wp) .and. all(profile%chi<=1._wp))
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
                    write(suffix,'(I5.5)')i
                    call emit_real('q_'//trim(suffix),profile%q(i))
                    call emit_real('r_'//trim(suffix),profile%r_m(i));call emit_real('m_'//trim(suffix),profile%m_kg(i))
                    call emit_real('p_'//trim(suffix),profile%p_pa(i));call emit_real('eps_'//trim(suffix),profile%energy_j_m3(i))
                    call emit_real('n_'//trim(suffix),profile%n_b_m3(i));call emit_real('mu_'//trim(suffix),profile%mu_j(i))
                    call emit_real('chi_'//trim(suffix),profile%chi(i));call emit_real('nu_'//trim(suffix),profile%nu(i))
                    call emit_real('lambda_'//trim(suffix),profile%lambda(i));call emit_real('dnu_'//trim(suffix),profile%dnu_dr(i))
                end do
            end if
        end subroutine output_profile

        subroutine finish()
            call emit_text('status',trim(status_name(status%code)))
            call emit_text('category',trim(status%category));call emit_text('detail',trim(status%detail))
            call emit_text('interpolation_id',TABLE_INTERPOLATION_ID)
            call emit_bool('table_admitted',table_admitted)
            call emit_bool('result_valid',profile%valid.and.is_ok(status))
        end subroutine finish
    end subroutine run_table_star_case
end module ns_table_star_cases
