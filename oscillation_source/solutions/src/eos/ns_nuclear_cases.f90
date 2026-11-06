module ns_nuclear_cases
    use ns_precision, only: wp
    use ns_physical_constants, only: inverse_fm3_in_si,mev_fm3_in_si
    use ns_status
    use ns_records, only: emit_text,emit_real,emit_bool,emit_int
    use ns_nuclear_toy
    use ns_beta_equilibrium
    use ns_liquid_drop
    implicit none
    private
    public :: is_nuclear_case,run_nuclear_case
contains
    logical function is_nuclear_case(name)
        character(len=*),intent(in)::name
        is_nuclear_case=trim(name)=='nuclear_constraints'.or.trim(name)=='nuclear_state'.or.&
            trim(name)=='beta_equilibrium'.or.trim(name)=='liquid_drop_scaling'.or.&
            trim(name)=='liquid_drop_candidate'.or.trim(name)=='liquid_drop_scan'
    end function

    subroutine run_nuclear_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token
        character(len=32)::objective
        integer::unit,ios,a,z,a_min,a_max,a_step,z_min,z_max,z_step
        logical::regulated
        real(wp)::n0_fm3,binding_mev,k0_mev,s0_mev,regulator_c,n_b_fm3,alpha,scan_u_max
        type(nuclear_model_t)::model
        type(nuclear_state_t)::state,sm,sp
        type(beta_state_t)::beta
        type(status_t)::status,st
        real(wp)::h,kfd,p0_mev_fm3,max_chi
        logical::causal
        namelist /request/ request_token
        namelist /parameters/ regulated,n0_fm3,binding_mev,k0_mev,s0_mev,regulator_c,n_b_fm3,alpha,scan_u_max,&
            a,z,a_min,a_max,a_step,z_min,z_max,z_step,objective
        request_token='';regulated=.true.;n0_fm3=0.16_wp;binding_mev=-16._wp;k0_mev=240._wp;s0_mev=30._wp
        regulator_c=0.2_wp;n_b_fm3=0.16_wp;alpha=0._wp;scan_u_max=10._wp
        a=56;z=26;a_min=20;a_max=180;a_step=1;z_min=1;z_max=60;z_step=1;objective='fixed_density_energy'
        call emit_text('case_id',trim(case_id))
        open(newunit=unit,file=path,status='old',action='read',iostat=ios)
        if(ios==0) then
            read(unit,nml=request,iostat=ios)
            if(ios==0) read(unit,nml=parameters,iostat=ios)
            close(unit)
        end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0.or.len_trim(request_token)/=32) then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid nuclear NAMELIST');call emit_status(status);return
        end if
        if(index(trim(case_id),'liquid_drop_')==1) then
            call run_liquid_drop_case(trim(case_id))
            return
        end if
        call init_nuclear_model(model,regulated,regulator_c,n0_fm3,binding_mev,k0_mev,s0_mev,status)
        if(.not.is_ok(status)) then;call emit_status(status);call emit_bool('model_ready',.false.);return;end if
        select case(trim(case_id))
        case('nuclear_constraints')
            call evaluate_nuclear(model,n0_fm3*inverse_fm3_in_si,0._wp,state,status)
            if(is_ok(status)) then
                h=1e-5_wp
                call evaluate_nuclear(model,n0_fm3*(1._wp-h)*inverse_fm3_in_si,0._wp,sm,st)
                if(is_ok(st)) call evaluate_nuclear(model,n0_fm3*(1._wp+h)*inverse_fm3_in_si,0._wp,sp,st)
                if(.not.is_ok(st)) status=st
            end if
            if(is_ok(status)) then
                kfd=9._wp*((sp%pressure_nuclear_pa-sm%pressure_nuclear_pa)/mev_fm3_in_si)/(2._wp*h*n0_fm3)
                p0_mev_fm3=state%pressure_nuclear_pa/mev_fm3_in_si
                call scan_nuclear_causality(model,scan_u_max,max_chi,causal,status)
            end if
            call emit_status(status);call emit_bool('model_ready',model%ready)
            if(.not.is_ok(status))return
            call emit_model(model)
            call emit_real('binding_residual_mev',state%symmetric_energy_mev-binding_mev)
            call emit_real('pressure_residual_mev_fm3',p0_mev_fm3)
            call emit_real('k0_fd_mev',kfd);call emit_real('k0_residual_mev',kfd-k0_mev)
            call emit_bool('regulator_fixed_externally',regulated)
            call emit_real('max_chi_scan',max_chi);call emit_bool('scan_causal',causal)
        case('nuclear_state')
            call evaluate_nuclear(model,n_b_fm3*inverse_fm3_in_si,alpha,state,status)
            call emit_status(status);call emit_bool('model_ready',model%ready);call emit_bool('state_valid',state%valid)
            if(.not.is_ok(status))return
            call emit_model(model);call emit_state(state)
            call emit_bool('kinetic_double_count_detected',.false.)
        case('beta_equilibrium')
            call solve_beta_equilibrium(model,n_b_fm3*inverse_fm3_in_si,beta,status)
            call emit_status(status);call emit_bool('model_ready',model%ready);call emit_bool('state_valid',beta%valid)
            if(.not.is_ok(status))return
            call emit_model(model)
            call emit_real('n_b_m3',beta%n_b_m3);call emit_real('proton_fraction',beta%proton_fraction)
            call emit_real('alpha',beta%alpha);call emit_real('n_p_m3',beta%n_p_m3);call emit_real('n_e_m3',beta%n_e_m3)
            call emit_real('beta_residual_mev',beta%beta_residual_mev);call emit_real('energy_j_m3',beta%energy_j_m3)
            call emit_real('pressure_pa',beta%pressure_pa);call emit_real('chi_equilibrium',beta%chi_equilibrium)
            call emit_bool('charge_neutrality_closed',beta%charge_neutrality_closed)
            call emit_bool('causal_local',beta%causal_local)
            call emit_bool('muons_included',.false.);call emit_int('iterations',int(beta%iterations,8))
        end select
    contains
        subroutine run_liquid_drop_case(which)
            character(len=*),intent(in)::which
            type(liquid_drop_model_t)::ld
            type(liquid_drop_terms_t)::t1,t8,bad
            type(liquid_drop_state_t)::ld_state
            type(liquid_drop_scan_t)::scan
            type(status_t)::ld_status,bad_status
            logical::objective_supported
            call init_liquid_drop_model(ld,ld_status)
            if(.not.is_ok(ld_status)) then;call emit_status(ld_status);return;end if
            select case(trim(which))
            case('liquid_drop_scaling')
                call liquid_drop_nuclear_terms(ld,8,2,t1,ld_status)
                if(is_ok(ld_status)) call liquid_drop_nuclear_terms(ld,64,2,t8,ld_status)
                call liquid_drop_nuclear_terms(ld,10,11,bad,bad_status)
                call emit_status(ld_status)
                if(.not.is_ok(ld_status))return
                call emit_real('surface_a8_over_a',t8%surface_mev/t1%surface_mev)
                call emit_real('coulomb_fixed_z_a8_over_a',t8%coulomb_mev/t1%coulomb_mev)
                call emit_bool('invalid_charge_rejected',.not.is_ok(bad_status))
                call emit_liquid_scope(ld)
            case('liquid_drop_candidate')
                call evaluate_liquid_drop_candidate(ld,n_b_fm3*inverse_fm3_in_si,a,z,ld_state,ld_status)
                call emit_status(ld_status);call emit_bool('candidate_valid',ld_state%valid)
                call emit_liquid_scope(ld)
                if(.not.is_ok(ld_status))return
                call emit_liquid_state(ld_state)
            case('liquid_drop_scan')
                objective_supported=trim(objective)=='fixed_density_energy'
                if(.not.objective_supported) then
                    call fail_status(ld_status,NS_OUT_OF_DOMAIN,'physics',&
                        'T12 does not implement fixed-pressure Gibbs minimization')
                    call emit_status(ld_status);call emit_bool('objective_supported',.false.)
                    call emit_bool('fixed_density_only',.true.);call emit_liquid_scope(ld);return
                end if
                call scan_liquid_drop_candidates(ld,n_b_fm3*inverse_fm3_in_si,a_min,a_max,a_step,&
                    z_min,z_max,z_step,scan,ld_status)
                call emit_status(ld_status);call emit_bool('objective_supported',.true.)
                call emit_bool('fixed_density_only',.true.);call emit_liquid_scope(ld)
                if(.not.is_ok(ld_status))return
                call emit_text('minimum_kind','discrete_candidate_set')
                call emit_bool('stationary_point_claimed',scan%stationary_point_claimed)
                call emit_bool('minimum_on_boundary',scan%minimum_on_boundary)
                call emit_bool('mesh_is_discrete',scan%mesh_is_discrete)
                call emit_int('candidate_count',int(scan%candidate_count,8))
                call emit_int('best_a',int(scan%best_a,8));call emit_int('best_z',int(scan%best_z,8))
                call emit_real('best_total_mev_per_baryon',scan%best%total_mev_per_baryon)
                call emit_bool('charge_neutrality_closed',scan%best%charge_neutrality_closed)
            end select
        end subroutine

        subroutine emit_liquid_scope(ld)
            type(liquid_drop_model_t),intent(in)::ld
            call emit_bool('full_star_provider_enabled',.false.)
            call emit_bool('dripped_neutrons_included',.false.)
            call emit_bool('pairing_included',.false.)
            call emit_bool('lattice_energy_included',.false.)
            call emit_real('n_b_min_fm3',ld%n_b_min_fm3);call emit_real('n_b_max_fm3',ld%n_b_max_fm3)
            call emit_real('a_volume_mev',ld%a_volume_mev);call emit_real('a_surface_mev',ld%a_surface_mev)
            call emit_real('a_coulomb_mev',ld%a_coulomb_mev);call emit_real('a_asymmetry_mev',ld%a_asymmetry_mev)
        end subroutine

        subroutine emit_liquid_state(s)
            type(liquid_drop_state_t),intent(in)::s
            call emit_int('a',int(s%a,8));call emit_int('z',int(s%z,8))
            call emit_real('n_b_m3',s%n_b_m3);call emit_real('n_b_fm3',s%n_b_fm3)
            call emit_real('n_nucleus_m3',s%n_nucleus_m3);call emit_real('n_e_m3',s%n_e_m3)
            call emit_real('volume_mev_per_baryon',s%volume_mev_per_baryon)
            call emit_real('surface_mev_per_baryon',s%surface_mev_per_baryon)
            call emit_real('coulomb_mev_per_baryon',s%coulomb_mev_per_baryon)
            call emit_real('asymmetry_mev_per_baryon',s%asymmetry_mev_per_baryon)
            call emit_real('nuclear_mev_per_baryon',s%nuclear_mev_per_baryon)
            call emit_real('electron_mev_per_baryon',s%electron_mev_per_baryon)
            call emit_real('total_mev_per_baryon',s%total_mev_per_baryon)
            call emit_bool('charge_neutrality_closed',s%charge_neutrality_closed)
        end subroutine

        subroutine emit_status(sta)
            type(status_t),intent(in)::sta
            call emit_text('status',trim(status_name(sta%code)))
            call emit_text('category',trim(sta%category))
            call emit_text('detail',trim(sta%detail))
        end subroutine
        subroutine emit_model(m)
            type(nuclear_model_t),intent(in)::m
            call emit_bool('regulated',m%regulated);call emit_real('n0_fm3',m%n0_fm3);call emit_real('binding_mev',m%binding_mev)
            call emit_real('k0_mev',m%k0_mev);call emit_real('s0_mev',m%s0_mev);call emit_real('regulator_c',m%regulator_c)
            call emit_real('t0_mev',m%t0_mev);call emit_real('a_mev',m%a_mev)
            call emit_real('b_mev',m%b_mev);call emit_real('sigma',m%sigma)
        end subroutine
        subroutine emit_state(s)
            type(nuclear_state_t),intent(in)::s
            call emit_real('n_b_m3',s%n_b_m3);call emit_real('n_b_fm3',s%n_b_fm3)
            call emit_real('u',s%u);call emit_real('alpha',s%alpha)
            call emit_real('proton_fraction',s%proton_fraction);call emit_real('symmetric_energy_mev',s%symmetric_energy_mev)
            call emit_real('symmetry_mev',s%symmetry_mev);call emit_real('energy_excess_mev',s%energy_excess_mev)
            call emit_real('exact_kinetic_mev',s%exact_kinetic_mev)
            call emit_real('exact_kinetic_excess_mev',s%exact_kinetic_excess_mev)
            call emit_real('pressure_nuclear_pa',s%pressure_nuclear_pa);call emit_real('energy_nuclear_j_m3',s%energy_nuclear_j_m3)
            call emit_real('chi_nuclear',s%chi_nuclear)
        end subroutine
    end subroutine
end module ns_nuclear_cases
