! Provided NAMELIST/TSV adapter. Calls only the selected workspace's EOS implementation.
module ns_polytrope_cases
    use ns_precision, only: wp
    use ns_status
    use ns_eos_types
    use ns_polytrope
    use ns_records, only: emit_real,emit_text,emit_bool,emit_int
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf,ieee_is_nan
    implicit none
    private
    public :: is_polytrope_case,run_polytrope_case
contains
    logical function is_polytrope_case(name)
        character(len=*),intent(in)::name
        is_polytrope_case=name=='poly_state'
    end function
    subroutine run_polytrope_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token
        character(len=40)::query,constructor,action,special_input
        real(wp)::n_anchor_m3,p_anchor_pa,gamma,c_j,n_min_m3,n_max_m3,n_b_m3,value,k_si,dndmu
        integer::side,u,ios
        type(polytrope_t)::model,other
        type(eos_state_t)::state,lower,upper,before,alternate,after,invalid
        type(status_t)::status,extra
        namelist /request/ request_token
        namelist /parameters/ n_anchor_m3,p_anchor_pa,gamma,c_j,n_min_m3,n_max_m3,n_b_m3,value, &
                              k_si,side,query,constructor,action,special_input
        request_token='';query='n';constructor='anchor';action='evaluate';special_input='none'
        n_anchor_m3=1e44_wp;p_anchor_pa=1e32_wp;gamma=2._wp;c_j=1.5e-10_wp
        n_min_m3=1e40_wp;n_max_m3=1e46_wp;n_b_m3=1e44_wp;value=1e32_wp;k_si=1e-56_wp;side=SIDE_REGULAR
        call clear_eos_state(state);call clear_eos_state(lower);call clear_eos_state(upper)
        dndmu=ieee_value(0._wp,ieee_quiet_nan)
        call emit_text('case_id',case_id)
        open(newunit=u,file=path,status='old',action='read',iostat=ios)
        if(ios==0) then
            read(u,nml=request,iostat=ios)
            if(ios==0) read(u,nml=parameters,iostat=ios)
            close(u)
        end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0) then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid EOS NAMELIST');call output_result();return
        end if
        select case(special_input)
        case('none')
        case('nan_n');n_b_m3=ieee_value(0._wp,ieee_quiet_nan)
        case('inf_pressure');query='pressure';value=ieee_value(0._wp,ieee_positive_inf)
        case('nan_mu');query='mu';value=ieee_value(0._wp,ieee_quiet_nan)
        case('nan_gamma');gamma=ieee_value(0._wp,ieee_quiet_nan)
        case('inf_c');c_j=ieee_value(0._wp,ieee_positive_inf)
        case default
            call fail_status(status,NS_INVALID_INPUT,'data','unknown special input');call output_result();return
        end select
        if(constructor=='anchor') then
            call init_polytrope(model,n_anchor_m3,p_anchor_pa,gamma,c_j,n_min_m3,n_max_m3,status)
        else if(constructor=='k') then
            call init_polytrope_k(model,k_si,gamma,c_j,n_anchor_m3,n_min_m3,n_max_m3,status)
        else
            call fail_status(status,NS_INVALID_INPUT,'data','unknown constructor')
        end if
        if(.not.is_ok(status)) then
            call output_result();return
        end if
        call polytrope_domain(model,lower,upper,status)
        if(.not.is_ok(status)) then
            call output_result();return
        end if
        if(action=='lifecycle') then
            call eos_at_n(model,n_anchor_m3,state,status)
            call eos_at_n(model,-1._wp,state,status)
            call init_polytrope(model,n_anchor_m3,-1._wp,gamma,c_j,n_min_m3,n_max_m3,extra)
            call emit_bool('failed_reinit_cleared',.not.polytrope_ready(model))
            call eos_at_n(model,n_anchor_m3,invalid,extra)
            call emit_bool('uninitialized_evaluation_rejected',extra%code==NS_UNINITIALIZED .and. .not.invalid%valid)
        else if(action=='interleaved') then
            call eos_at_n(model,n_anchor_m3,before,status)
            call init_polytrope(other,n_anchor_m3,p_anchor_pa,gamma,2._wp*c_j,n_min_m3,n_max_m3,extra)
            if(is_ok(extra)) call eos_at_n(other,n_anchor_m3,alternate,extra)
            call eos_at_n(model,n_anchor_m3,after,status)
            state=after
            call emit_bool('independent_models',is_ok(extra) .and. alternate%energy_j_m3>before%energy_j_m3 &
                            .and. alternate%pressure_pa==before%pressure_pa)
            call emit_bool('first_model_unchanged',after%energy_j_m3==before%energy_j_m3 &
                            .and. after%mu_j==before%mu_j .and. after%chi==before%chi)
        else if(action=='evaluate') then
            select case(query)
            case('n');call eos_at_n(model,n_b_m3,state,status)
            case('pressure');call eos_at_pressure(model,value,side,state,status)
            case('mu');call eos_at_mu(model,value,state,dndmu,status)
            case default;call fail_status(status,NS_INVALID_INPUT,'data','unknown query')
            end select
        else
            call fail_status(status,NS_INVALID_INPUT,'data','unknown EOS action')
        end if
        call output_result()
    contains
        subroutine output_result()
            call emit_text('status',trim(status_name(status%code)))
            call emit_text('category',trim(status%category));call emit_text('detail',trim(status%detail))
            call emit_bool('model_ready',polytrope_ready(model));call emit_bool('result_valid',state%valid)
            call emit_bool('state_poisoned',all(ieee_is_nan([state%pressure_pa,state%energy_j_m3, &
                                                         state%n_b_m3,state%mu_j,state%chi])))
            if(polytrope_ready(model) .and. lower%valid .and. upper%valid) then
                call emit_real('n_min_m3',lower%n_b_m3);call emit_real('n_max_m3',upper%n_b_m3)
                call emit_real('p_min_pa',lower%pressure_pa);call emit_real('p_max_pa',upper%pressure_pa)
                call emit_real('mu_min_j',lower%mu_j);call emit_real('mu_max_j',upper%mu_j)
                call emit_real('chi_max',upper%chi)
            end if
            if(is_ok(status) .and. state%valid) then
                call emit_real('n_b_m3',state%n_b_m3);call emit_real('pressure_pa',state%pressure_pa)
                call emit_real('energy_j_m3',state%energy_j_m3);call emit_real('mu_j',state%mu_j)
                call emit_real('chi',state%chi);call emit_int('side',int(state%side,int64))
                if(query=='mu' .and. action=='evaluate') call emit_real('dn_dmu_m3_per_j',dndmu)
            end if
        end subroutine
    end subroutine
end module
