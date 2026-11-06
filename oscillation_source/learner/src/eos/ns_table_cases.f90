module ns_table_cases
    use ns_precision, only: wp
    use ns_status
    use ns_records, only: emit_real,emit_text,emit_bool,emit_int
    use ns_eos_types
    use ns_table_io
    use ns_table_eos
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private
    public :: is_table_case,run_table_case
contains
    logical function is_table_case(name)
        character(len=*),intent(in)::name
        is_table_case=name=='table_diagnostics'.or.name=='table_state'.or.name=='real_table_admission'
    end function is_table_case

    subroutine run_table_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token,query_kind
        character(len=512)::table_path
        real(wp)::query_value,dn
        integer::side,u,ios
        logical::table_admitted,tov_attempted,modes_attempted
        type(normalized_table_t)::raw
        type(table_eos_t)::model
        type(table_audit_t)::audit
        type(eos_state_t)::state
        type(status_t)::status
        namelist /request/ request_token
        namelist /parameters/ table_path,query_kind,query_value,side
        request_token='';table_path='';query_kind='mu';query_value=0._wp;side=SIDE_REGULAR
        table_admitted=.false.;tov_attempted=.false.;modes_attempted=.false.
        call clear_eos_state(state)
        call emit_text('case_id',case_id)
        open(newunit=u,file=path,status='old',action='read',iostat=ios)
        if(ios==0) then
            read(u,nml=request,iostat=ios);if(ios==0) read(u,nml=parameters,iostat=ios);close(u)
        end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0) then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid T14 table NAMELIST');call output_all();return
        end if
        call load_normalized_table(trim(table_path),raw,status)
        if(.not.is_ok(status)) then;call output_all();return;end if
        if(case_id=='table_diagnostics') then
            call audit_table_eos(raw,audit,status)
        else if(case_id=='real_table_admission') then
            call audit_table_eos(raw,audit,status)
            if(is_ok(status)) then
                table_admitted=audit%fully_admitted
                if(.not.table_admitted) call fail_status(status,NS_MATCHING_INCONSISTENT,'physics', &
                    'native cubic Hermite interval admissibility blocks full-table provider')
            end if
        else if(case_id=='table_state') then
            call init_table_eos(raw,model,audit,status)
            if(is_ok(status)) then
                select case(trim(query_kind))
                case('mu');call eos_table_at_mu(model,query_value,side,state,dn,status)
                case('pressure');call eos_table_at_pressure(model,query_value,side,state,status)
                case default;call fail_status(status,NS_INVALID_INPUT,'data','unknown table query kind')
                end select
            end if
        else
            call fail_status(status,NS_INVALID_INPUT,'data','unknown T14 table case')
        end if
        call output_all()
    contains
        subroutine output_all()
            call emit_text('status',trim(status_name(status%code)))
            call emit_text('category',trim(status%category));call emit_text('detail',trim(status%detail))
            call emit_text('interpolation_id',TABLE_INTERPOLATION_ID)
            if(case_id=='real_table_admission') then
                call emit_bool('table_admitted',table_admitted)
                call emit_bool('tov_attempted',tov_attempted)
                call emit_bool('modes_attempted',modes_attempted)
                call emit_text('admission_blocker','native_cubic_hermite_interval_admissibility')
            end if
            call emit_bool('result_valid',is_ok(status) .and. &
                (case_id=='table_diagnostics' .or. case_id=='real_table_admission' .or. state%valid))
            if(audit%row_count>0) then
                call emit_int('row_count',int(audit%row_count,int64))
                call emit_int('interval_count',int(audit%interval_count,int64))
                call emit_int('admitted_intervals',int(audit%admitted_intervals,int64))
                call emit_int('rejected_intervals',int(audit%rejected_intervals,int64))
                call emit_int('first_rejected_interval',int(audit%first_rejected_interval,int64))
                call emit_bool('mu_strict_increasing',audit%mu_strict_increasing)
                call emit_bool('pressure_strict_increasing',audit%pressure_strict_increasing)
                call emit_bool('fully_admitted',audit%fully_admitted)
                call emit_real('raw_energy_closure_max_rel',audit%raw_energy_closure_max_rel)
                call emit_real('raw_trapezoid_first_law_max_rel',audit%raw_trapezoid_first_law_max_rel)
                if(audit%first_rejected_interval>0) then
                    call emit_real('first_bad_min_d2p',audit%first_bad_min_d2p)
                    call emit_real('first_bad_min_causal_margin',audit%first_bad_min_causal_margin)
                end if
            end if
            if(is_ok(status).and.state%valid) then
                call emit_real('pressure_pa',state%pressure_pa);call emit_real('energy_j_m3',state%energy_j_m3)
                call emit_real('n_b_m3',state%n_b_m3);call emit_real('mu_j',state%mu_j);call emit_real('chi',state%chi)
                call emit_int('side',int(state%side,int64))
                if(trim(query_kind)=='mu') call emit_real('dn_dmu_m3_per_j',dn)
            end if
        end subroutine output_all
    end subroutine run_table_case
end module ns_table_cases
