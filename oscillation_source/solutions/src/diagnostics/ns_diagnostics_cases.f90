module ns_diagnostics_cases
    use ns_precision, only: wp
    use ns_status
    use ns_records, only: emit_real,emit_text,emit_bool,emit_int
    use ns_background_access, only: background_table_t,load_background_csv
    use ns_table_io, only: normalized_table_t,load_normalized_table
    use ns_table_eos, only: table_audit_t,audit_table_eos
    use ns_controlled_table_variant, only: controlled_variant_audit_t,build_controlled_variant, &
        CONTROLLED_VARIANT_ID
    use ns_structure_diagnostics, only: structure_diagnostics_t,compute_structure_diagnostics,STRUCTURE_DIAGNOSTIC_ID
    use ns_mode_diagnostics, only: mode_difference_t,paired_refinement_t, &
        compare_mode_frequencies,classify_paired_refinement
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private
    public :: is_diagnostics_case,run_diagnostics_case
contains
    logical function is_diagnostics_case(name)
        character(len=*),intent(in)::name
        is_diagnostics_case=name=='structure_diagnostics'.or.name=='mode_difference'.or. &
            name=='paired_refinement'.or.name=='controlled_table_variant'
    end function is_diagnostics_case

    subroutine run_diagnostics_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token
        character(len=4096)::background_path,table_path
        character(len=24)::f_ref_label,f_cmp_label,p1_ref_label,p1_cmp_label
        real(wp)::split_pressure_pa,f_hz,p1_hz,f_ref_hz,f_cmp_hz,p1_ref_hz,p1_cmp_hz,resolution_hz
        real(wp)::delta_f_coarse_hz,delta_f_medium_hz,delta_f_fine_hz
        real(wp)::delta_p1_coarse_hz,delta_p1_medium_hz,delta_p1_fine_hz,absolute_resolution_hz
        real(wp)::mass_residual_ref_kg,mass_residual_cmp_kg,df_dmass_hz_per_kg,dp1_dmass_hz_per_kg,alpha
        integer::split_index
        integer::u,ios
        type(status_t)::status
        type(background_table_t)::bg
        type(structure_diagnostics_t)::sd
        type(mode_difference_t)::md
        type(paired_refinement_t)::pr
        type(normalized_table_t)::raw_table,variant_table
        type(controlled_variant_audit_t)::variant_audit
        type(table_audit_t)::table_audit
        namelist /request/ request_token
        namelist /parameters/ background_path,table_path,split_pressure_pa,f_hz,p1_hz, &
            f_ref_hz,f_cmp_hz,p1_ref_hz,p1_cmp_hz,resolution_hz,delta_f_coarse_hz, &
            delta_f_medium_hz,delta_f_fine_hz,delta_p1_coarse_hz,delta_p1_medium_hz,delta_p1_fine_hz, &
            absolute_resolution_hz,mass_residual_ref_kg,mass_residual_cmp_kg,df_dmass_hz_per_kg, &
            dp1_dmass_hz_per_kg,f_ref_label,f_cmp_label,p1_ref_label,p1_cmp_label,alpha,split_index
        request_token='';background_path='';table_path='';split_pressure_pa=0._wp;f_hz=0._wp;p1_hz=0._wp
        f_ref_hz=0._wp;f_cmp_hz=0._wp;p1_ref_hz=0._wp;p1_cmp_hz=0._wp;resolution_hz=0._wp
        delta_f_coarse_hz=0._wp;delta_f_medium_hz=0._wp;delta_f_fine_hz=0._wp
        delta_p1_coarse_hz=0._wp;delta_p1_medium_hz=0._wp;delta_p1_fine_hz=0._wp;absolute_resolution_hz=0._wp
        mass_residual_ref_kg=0._wp;mass_residual_cmp_kg=0._wp;df_dmass_hz_per_kg=0._wp;dp1_dmass_hz_per_kg=0._wp
        f_ref_label='';f_cmp_label='';p1_ref_label='';p1_cmp_label='';alpha=0._wp;split_index=0
        call emit_text('case_id',trim(case_id))
        open(newunit=u,file=path,status='old',action='read',iostat=ios)
        if(ios==0) then
            read(u,nml=request,iostat=ios);if(ios==0) read(u,nml=parameters,iostat=ios);close(u)
        end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0) then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid T15 diagnostics NAMELIST');call finish(.false.);return
        end if
        select case(trim(case_id))
        case('structure_diagnostics')
            call load_background_csv(trim(background_path),bg,status)
            if(is_ok(status)) call compute_structure_diagnostics(bg,split_pressure_pa,f_hz,p1_hz,sd,status)
            if(is_ok(status)) then
                call emit_real('radius_m',sd%radius_m);call emit_real('mass_kg',sd%mass_kg)
                call emit_real('proper_radius_m',sd%proper_radius_m);call emit_real('acoustic_time_s',sd%acoustic_time_s)
                call emit_real('r_split_m',sd%r_split_m);call emit_real('outer_areal_thickness_m',sd%outer_areal_thickness_m)
                call emit_real('outer_proper_thickness_m',sd%outer_proper_thickness_m)
                call emit_real('outer_acoustic_time_s',sd%outer_acoustic_time_s)
                call emit_real('sound_speed_min_m_s',sd%sound_speed_min_m_s)
                call emit_real('sound_speed_max_m_s',sd%sound_speed_max_m_s)
                call emit_real('f_hat',sd%f_hat);call emit_real('p1_hat',sd%p1_hat)
                call emit_text('diagnostic_id',STRUCTURE_DIAGNOSTIC_ID)
            end if
            call finish(sd%valid)
        case('mode_difference')
            call compare_mode_frequencies(f_ref_hz,f_cmp_hz,p1_ref_hz,p1_cmp_hz,resolution_hz,md,status)
            if(is_ok(status)) then
                call emit_real('delta_f_hz',md%delta_f_hz);call emit_real('delta_p1_hz',md%delta_p1_hz)
                call emit_real('delta_f_relative',md%delta_f_relative);call emit_real('delta_p1_relative',md%delta_p1_relative)
                call emit_text('f_interpretation',trim(md%f_interpretation))
                call emit_text('p1_interpretation',trim(md%p1_interpretation))
                call emit_bool('f_relative_defined',md%f_relative_defined)
                call emit_bool('p1_relative_defined',md%p1_relative_defined)
                call emit_bool('identity_compared',md%identity_compared)
            end if
            call finish(md%valid)
        case('paired_refinement')
            call classify_paired_refinement(delta_f_coarse_hz,delta_f_medium_hz,delta_f_fine_hz, &
                delta_p1_coarse_hz,delta_p1_medium_hz,delta_p1_fine_hz,absolute_resolution_hz, &
                mass_residual_ref_kg,mass_residual_cmp_kg,df_dmass_hz_per_kg,dp1_dmass_hz_per_kg, &
                f_ref_label,f_cmp_label,p1_ref_label,p1_cmp_label,pr,status)
            call emit_bool('mode_identity_consistent',pr%mode_identity_consistent)
            if(is_ok(status)) then
                call emit_real('f_effect_hz',pr%f_effect_hz);call emit_real('p1_effect_hz',pr%p1_effect_hz)
                call emit_real('f_numerical_floor_hz',pr%f_numerical_floor_hz)
                call emit_real('p1_numerical_floor_hz',pr%p1_numerical_floor_hz)
                call emit_real('f_mass_floor_hz',pr%f_mass_floor_hz);call emit_real('p1_mass_floor_hz',pr%p1_mass_floor_hz)
                call emit_real('f_resolution_floor_hz',pr%f_resolution_floor_hz)
                call emit_real('p1_resolution_floor_hz',pr%p1_resolution_floor_hz)
                call emit_bool('f_convergent',pr%f_convergent);call emit_bool('p1_convergent',pr%p1_convergent)
                call emit_text('f_interpretation',trim(pr%f_interpretation))
                call emit_text('p1_interpretation',trim(pr%p1_interpretation))
            end if
            call finish(pr%valid)
        case('controlled_table_variant')
            call load_normalized_table(trim(table_path),raw_table,status)
            if(.not.is_ok(status)) then;call finish(.false.);return;end if
            call build_controlled_variant(raw_table,alpha,split_index,variant_table,variant_audit,status)
            if(.not.is_ok(status)) then
                call emit_bool('core_function_unchanged',variant_audit%core_function_unchanged)
                call finish(.false.);return
            end if
            call audit_table_eos(variant_table,table_audit,status)
            if(is_ok(status).and..not.table_audit%fully_admitted) &
                call fail_status(status,NS_MATCHING_INCONSISTENT,'physics','controlled variant failed native-Hermite admission')
            call emit_text('variant_id',CONTROLLED_VARIANT_ID)
            call emit_int('split_index',int(variant_audit%split_index,int64))
            call emit_real('split_n_b_fm3',variant_audit%split_n_b_m3/1.0e45_wp)
            call emit_real('split_mu_j',variant_audit%split_mu_j)
            call emit_real('max_low_density_pressure_rel_change',variant_audit%max_low_density_pressure_rel_change)
            call emit_real('max_low_density_density_rel_change',variant_audit%max_low_density_density_rel_change)
            call emit_real('perturbation_endpoint_pressure_change_pa',variant_audit%perturbation_endpoint_pressure_change_pa)
            call emit_real('perturbation_endpoint_density_change_m3',variant_audit%perturbation_endpoint_density_change_m3)
            call emit_bool('core_function_unchanged',variant_audit%core_function_unchanged)
            call emit_int('admitted_intervals',int(table_audit%admitted_intervals,int64))
            call emit_int('rejected_intervals',int(table_audit%rejected_intervals,int64))
            call emit_bool('fully_admitted',table_audit%fully_admitted)
            call finish(variant_audit%valid.and.table_audit%fully_admitted.and.is_ok(status))
        case default
            call fail_status(status,NS_INVALID_INPUT,'implementation','unknown T15 diagnostics case');call finish(.false.)
        end select
    contains
        subroutine finish(valid)
            logical,intent(in)::valid
            call emit_text('status',trim(status_name(status%code)))
            call emit_text('category',trim(status%category));call emit_text('detail',trim(status%detail))
            call emit_bool('result_valid',valid.and.is_ok(status))
        end subroutine finish
    end subroutine run_diagnostics_case
end module ns_diagnostics_cases
