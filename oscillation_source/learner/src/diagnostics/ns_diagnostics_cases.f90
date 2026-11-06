module ns_diagnostics_cases
    use ns_precision,only:wp
    use ns_status
    use ns_records,only:emit_real,emit_text,emit_bool
    use ns_background_access,only:background_table_t,load_background_csv
    use ns_structure_diagnostics,only:structure_diagnostics_t,compute_structure_diagnostics,STRUCTURE_DIAGNOSTIC_ID
    use ns_mode_diagnostics,only:mode_difference_t,paired_refinement_t,compare_mode_frequencies,classify_paired_refinement
    implicit none
    private
    public :: is_diagnostics_case,run_diagnostics_case
contains
    logical function is_diagnostics_case(name)
        character(len=*),intent(in)::name
        is_diagnostics_case=name=='structure_diagnostics'.or.name=='mode_difference'.or.name=='paired_refinement'
    end function
    subroutine run_diagnostics_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token
        character(len=4096)::background_path
        character(len=24)::f_ref_label,f_cmp_label,p1_ref_label,p1_cmp_label
        real(wp)::split_pressure_pa,f_hz,p1_hz,f_ref_hz,f_cmp_hz,p1_ref_hz,p1_cmp_hz,resolution_hz
        real(wp)::delta_f_coarse_hz,delta_f_medium_hz,delta_f_fine_hz
        real(wp)::delta_p1_coarse_hz,delta_p1_medium_hz,delta_p1_fine_hz,absolute_resolution_hz
        real(wp)::mass_residual_ref_kg,mass_residual_cmp_kg,df_dmass_hz_per_kg,dp1_dmass_hz_per_kg
        integer::u,ios
        type(status_t)::status
        type(background_table_t)::bg
        type(structure_diagnostics_t)::sd
        type(mode_difference_t)::md
        type(paired_refinement_t)::pr
        namelist /request/ request_token
        namelist /parameters/ background_path,split_pressure_pa,f_hz,p1_hz,f_ref_hz,f_cmp_hz,p1_ref_hz,p1_cmp_hz, &
            resolution_hz,delta_f_coarse_hz,delta_f_medium_hz,delta_f_fine_hz,delta_p1_coarse_hz, &
            delta_p1_medium_hz,delta_p1_fine_hz,absolute_resolution_hz,mass_residual_ref_kg,mass_residual_cmp_kg, &
            df_dmass_hz_per_kg,dp1_dmass_hz_per_kg,f_ref_label,f_cmp_label,p1_ref_label,p1_cmp_label
        request_token='';background_path='';split_pressure_pa=0._wp;f_hz=0._wp;p1_hz=0._wp
        f_ref_hz=0._wp;f_cmp_hz=0._wp;p1_ref_hz=0._wp;p1_cmp_hz=0._wp;resolution_hz=0._wp
        delta_f_coarse_hz=0._wp;delta_f_medium_hz=0._wp;delta_f_fine_hz=0._wp
        delta_p1_coarse_hz=0._wp;delta_p1_medium_hz=0._wp;delta_p1_fine_hz=0._wp;absolute_resolution_hz=0._wp
        mass_residual_ref_kg=0._wp;mass_residual_cmp_kg=0._wp;df_dmass_hz_per_kg=0._wp;dp1_dmass_hz_per_kg=0._wp
        f_ref_label='';f_cmp_label='';p1_ref_label='';p1_cmp_label=''
        call emit_text('case_id',trim(case_id))
        open(newunit=u,file=path,status='old',action='read',iostat=ios)
        if(ios==0)then;read(u,nml=request,iostat=ios);if(ios==0)read(u,nml=parameters,iostat=ios);close(u);end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0)then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid learner diagnostics NAMELIST')
            call finish(.false.);return
        end if
        select case(trim(case_id))
        case('structure_diagnostics')
            call load_background_csv(trim(background_path),bg,status)
            if(is_ok(status))call compute_structure_diagnostics(bg,split_pressure_pa,f_hz,p1_hz,sd,status)
            if(is_ok(status))then
                call emit_real('radius_m',sd%radius_m);call emit_real('mass_kg',sd%mass_kg)
                call emit_real('acoustic_time_s',sd%acoustic_time_s);call emit_text('diagnostic_id',STRUCTURE_DIAGNOSTIC_ID)
            end if
            call finish(sd%valid)
        case('mode_difference')
            call compare_mode_frequencies(f_ref_hz,f_cmp_hz,p1_ref_hz,p1_cmp_hz,resolution_hz,md,status)
            if(is_ok(status))then
                call emit_real('delta_f_hz',md%delta_f_hz);call emit_real('delta_p1_hz',md%delta_p1_hz)
                call emit_text('f_interpretation',trim(md%f_interpretation))
                call emit_text('p1_interpretation',trim(md%p1_interpretation))
            end if
            call finish(md%valid)
        case('paired_refinement')
            call classify_paired_refinement(delta_f_coarse_hz,delta_f_medium_hz,delta_f_fine_hz, &
                delta_p1_coarse_hz,delta_p1_medium_hz,delta_p1_fine_hz,absolute_resolution_hz, &
                mass_residual_ref_kg,mass_residual_cmp_kg,df_dmass_hz_per_kg,dp1_dmass_hz_per_kg, &
                f_ref_label,f_cmp_label,p1_ref_label,p1_cmp_label,pr,status)
            if(is_ok(status))then
                call emit_text('f_interpretation',trim(pr%f_interpretation))
                call emit_text('p1_interpretation',trim(pr%p1_interpretation))
            end if
            call finish(pr%valid)
        case default
            call fail_status(status,NS_INVALID_INPUT,'implementation','unknown learner diagnostics case');call finish(.false.)
        end select
    contains
        subroutine finish(valid)
            logical,intent(in)::valid
            call emit_text('status',trim(status_name(status%code)));call emit_text('category',trim(status%category))
            call emit_text('detail',trim(status%detail));call emit_bool('result_valid',valid.and.is_ok(status))
        end subroutine
    end subroutine
end module ns_diagnostics_cases
