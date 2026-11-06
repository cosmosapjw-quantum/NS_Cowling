! T15 controlled low-density table variant through the existing same-mass TOV solver.
module ns_controlled_table_cases
    use ns_precision, only: wp
    use ns_status
    use ns_records, only: emit_real,emit_text,emit_bool,emit_int
    use ns_physical_constants, only: nominal_solar_mass_kg
    use ns_eos_types, only: eos_state_t
    use ns_table_io, only: normalized_table_t,load_normalized_table
    use ns_table_eos, only: table_eos_t,table_audit_t,init_table_eos,eos_table_at_pressure
    use ns_controlled_table_variant, only: controlled_variant_audit_t,build_controlled_variant, &
        CONTROLLED_VARIANT_ID
    use ns_background_records, only: star_options_t,mass_options_t,star_result_t
    use ns_target_mass, only: solve_target_mass
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private
    public :: is_controlled_table_case,run_controlled_table_case
contains
    logical function is_controlled_table_case(name)
        character(len=*),intent(in)::name
        is_controlled_table_case=name=='controlled_table_target_mass'
    end function is_controlled_table_case

    subroutine run_controlled_table_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token,suffix
        character(len=4096)::table_path
        character(len=64)::branch_id
        real(wp)::alpha,m_target_kg,pc_lo_pa,pc_hi_pa,p_cut_pa,r0_m,r_limit_m,rtol
        real(wp)::atol_r_m,atol_mass_kg,atol_nu,dq_initial,dq_min,dq_max,nu_start
        real(wp)::mass_atol_kg,mass_rtol,root_xtol,root_rtol
        integer::split_index,max_steps,root_max_iter,scan_n,u,ios,i,n
        logical::emit_profile,table_admitted
        type(normalized_table_t)::raw,variant
        type(controlled_variant_audit_t)::variant_audit
        type(table_eos_t)::model
        type(table_audit_t)::table_audit
        type(star_options_t)::star_options
        type(mass_options_t)::mass_options
        type(star_result_t)::result
        type(status_t)::status
        namelist /request/ request_token
        namelist /parameters/ table_path,alpha,split_index,m_target_kg,pc_lo_pa,pc_hi_pa,p_cut_pa, &
            r0_m,r_limit_m,rtol,atol_r_m,atol_mass_kg,atol_nu,dq_initial,dq_min,dq_max,nu_start, &
            max_steps,emit_profile,mass_atol_kg,mass_rtol,root_xtol,root_rtol,root_max_iter,scan_n,branch_id
        request_token='';table_path='';alpha=0._wp;split_index=0;m_target_kg=0._wp
        pc_lo_pa=1e34_wp;pc_hi_pa=2e34_wp;p_cut_pa=1e20_wp;r0_m=1._wp;r_limit_m=1e6_wp;rtol=1e-8_wp
        atol_r_m=1e-6_wp;atol_mass_kg=1e14_wp;atol_nu=1e-12_wp;dq_initial=1e-5_wp;dq_min=1e-16_wp
        dq_max=0.2_wp;nu_start=0._wp;max_steps=30000;emit_profile=.false.;mass_atol_kg=1e22_wp
        mass_rtol=0._wp;root_xtol=1e-10_wp;root_rtol=0._wp;root_max_iter=80;scan_n=7
        branch_id='t15_controlled_variant_same_mass';table_admitted=.false.;result=star_result_t()
        call emit_text('case_id',trim(case_id))
        open(newunit=u,file=path,status='old',action='read',iostat=ios)
        if(ios==0) then
            read(u,nml=request,iostat=ios);if(ios==0) read(u,nml=parameters,iostat=ios);close(u)
        end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0) then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid T15 controlled-table NAMELIST')
            call finish();return
        end if
        call load_normalized_table(trim(table_path),raw,status)
        if(.not.is_ok(status)) then;call finish();return;end if
        call build_controlled_variant(raw,alpha,split_index,variant,variant_audit,status)
        if(.not.is_ok(status)) then;call finish();return;end if
        call init_table_eos(variant,model,table_audit,status)
        table_admitted=is_ok(status).and.table_audit%fully_admitted
        if(.not.table_admitted) then;call finish();return;end if
        star_options%r0_m=r0_m;star_options%r_limit_m=r_limit_m;star_options%p_cut_pa=p_cut_pa
        star_options%nu_start=nu_start;star_options%surface_policy='finite_cutoff_schwarzschild_time'
        star_options%eos_id='controlled_table_variant_low_density'
        star_options%ode%atol=[atol_r_m,atol_mass_kg,atol_nu];star_options%ode%rtol=rtol
        star_options%ode%h_initial=dq_initial;star_options%ode%h_min=dq_min;star_options%ode%h_max=dq_max
        star_options%ode%max_steps=max_steps
        mass_options%scan_n=scan_n;mass_options%mass_atol_kg=mass_atol_kg;mass_options%mass_rtol=mass_rtol
        mass_options%root%xtol=root_xtol;mass_options%root%rtol=root_rtol;mass_options%root%ftol=0._wp
        mass_options%root%max_iter=root_max_iter;mass_options%branch_id=branch_id
        call solve_target_mass(eos,m_target_kg,[pc_lo_pa,pc_hi_pa],star_options,mass_options,result)
        status=result%status
        if(result%valid.and.is_ok(status)) call output_result()
        call finish()
    contains
        subroutine eos(p,side,state,st)
            real(wp),intent(in)::p
            integer,intent(in)::side
            type(eos_state_t),intent(out)::state
            type(status_t),intent(out)::st
            call eos_table_at_pressure(model,p,side,state,st)
        end subroutine eos

        subroutine output_result()
            n=size(result%profile%r_m)
            call emit_text('variant_id',CONTROLLED_VARIANT_ID)
            call emit_real('alpha',variant_audit%alpha);call emit_int('split_index',int(variant_audit%split_index,int64))
            call emit_real('target_mass_kg',result%target_mass_kg);call emit_real('mass_kg',result%profile%mass_kg)
            call emit_real('mass_residual_kg',result%mass_residual_kg)
            call emit_real('mass_nominal_solar',result%profile%mass_kg/nominal_solar_mass_kg)
            call emit_real('pc_pa',result%profile%pc_pa);call emit_real('radius_m',result%profile%radius_m)
            call emit_real('p_cut_pa',result%profile%p_cut_pa);call emit_int('profile_n',int(n,int64))
            call emit_int('star_evaluations',int(result%star_evaluations,int64))
            call emit_bool('scan_monotone',result%scan_monotone)
            call emit_bool('cowling_background_admitted',all(result%profile%chi>0._wp).and.all(result%profile%chi<=1._wp))
            if(emit_profile) then
                do i=1,n
                    write(suffix,'(I5.5)')i
                    call emit_real('q_'//trim(suffix),result%profile%q(i))
                    call emit_real('r_'//trim(suffix),result%profile%r_m(i))
                    call emit_real('m_'//trim(suffix),result%profile%m_kg(i))
                    call emit_real('p_'//trim(suffix),result%profile%p_pa(i))
                    call emit_real('eps_'//trim(suffix),result%profile%energy_j_m3(i))
                    call emit_real('n_'//trim(suffix),result%profile%n_b_m3(i))
                    call emit_real('mu_'//trim(suffix),result%profile%mu_j(i))
                    call emit_real('chi_'//trim(suffix),result%profile%chi(i))
                    call emit_real('nu_'//trim(suffix),result%profile%nu(i))
                    call emit_real('lambda_'//trim(suffix),result%profile%lambda(i))
                    call emit_real('dnu_'//trim(suffix),result%profile%dnu_dr(i))
                end do
            end if
        end subroutine output_result

        subroutine finish()
            call emit_text('status',trim(status_name(status%code)))
            call emit_text('category',trim(status%category));call emit_text('detail',trim(status%detail))
            call emit_bool('table_admitted',table_admitted)
            call emit_bool('core_function_unchanged',variant_audit%core_function_unchanged)
            call emit_bool('result_valid',result%valid.and.is_ok(status))
        end subroutine finish
    end subroutine run_controlled_table_case
end module ns_controlled_table_cases
