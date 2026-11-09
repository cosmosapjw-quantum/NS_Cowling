program case_driver
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use ns_bootstrap_protocol, only: read_bootstrap_input, emit_header, emit_text, &
                                    emit_int, emit_real, emit_bool
    use ns_workspace, only: workspace_id
    use ns_support_cases, only: is_support_case, run_support_case
    use ns_numerics_cases, only: is_numerics_case, run_numerics_case
    use ns_matching_cases, only: is_matching_case,run_matching_case
    use ns_polytrope_cases, only: is_polytrope_case, run_polytrope_case
    use ns_tov_cases, only: is_tov_case, run_tov_case
    use ns_shooting_cases, only: is_shooting_case,run_shooting_case
    use ns_cowling_cases, only: is_cowling_case,run_cowling_case
    use ns_mode_cases, only: is_mode_case,run_mode_case
    use ns_interface_cases, only: is_interface_case,run_interface_case
    use ns_nuclear_cases, only: is_nuclear_case,run_nuclear_case
    use ns_table_cases, only: is_table_case,run_table_case
    use ns_table_star_cases, only: is_table_star_case,run_table_star_case
    use ns_diagnostics_cases, only: is_diagnostics_case,run_diagnostics_case
    use ns_controlled_table_cases, only: is_controlled_table_case,run_controlled_table_case
    use ns_maxwell_cases, only: is_maxwell_case,run_maxwell_case
    implicit none
    character(len=128) :: case_id
    character(len=4096) :: input_path
    character(len=32) :: token
    character(len=256) :: text_value
    real(real64) :: real_value
    integer(int64) :: int_value
    logical :: bool_value, ok
    integer :: arg_status
    case_id = ''; input_path = ''; token = ''
    call emit_header()
    call emit_int('protocol_version', 1_int64)
    call emit_text('workspace', workspace_id)
    if (command_argument_count() /= 2) then
        call emit_text('status', 'invalid_input')
    else
        call get_command_argument(1, case_id, status=arg_status)
        if (arg_status == 0) call get_command_argument(2, input_path, status=arg_status)
        if (arg_status /= 0) then
            call emit_text('status', 'invalid_input')
        else if (is_maxwell_case(trim(case_id))) then
            call run_maxwell_case(trim(case_id),trim(input_path))
        else if (is_controlled_table_case(trim(case_id))) then
            call run_controlled_table_case(trim(case_id),trim(input_path))
        else if (is_diagnostics_case(trim(case_id))) then
            call run_diagnostics_case(trim(case_id),trim(input_path))
        else if (is_table_star_case(trim(case_id))) then
            call run_table_star_case(trim(case_id),trim(input_path))
        else if (is_table_case(trim(case_id))) then
            call run_table_case(trim(case_id),trim(input_path))
        else if (is_nuclear_case(trim(case_id))) then
            call run_nuclear_case(trim(case_id),trim(input_path))
        else if (is_interface_case(trim(case_id))) then
            call run_interface_case(trim(case_id),trim(input_path))
        else if (is_mode_case(trim(case_id))) then
            call run_mode_case(trim(case_id),trim(input_path))
        else if (is_cowling_case(trim(case_id))) then
            call run_cowling_case(trim(case_id),trim(input_path))
        else if (is_shooting_case(trim(case_id))) then
            call run_shooting_case(trim(case_id),trim(input_path))
        else if (is_tov_case(trim(case_id))) then
            call run_tov_case(trim(case_id), trim(input_path))
        else if (is_matching_case(trim(case_id))) then
            call run_matching_case(trim(case_id),trim(input_path))
        else if (is_polytrope_case(trim(case_id))) then
            call run_polytrope_case(trim(case_id), trim(input_path))
        else if (is_numerics_case(trim(case_id))) then
            call run_numerics_case(trim(case_id), trim(input_path))
        else if (is_support_case(trim(case_id))) then
            call run_support_case(trim(case_id), trim(input_path))
        else
            call read_bootstrap_input(trim(input_path), token, real_value, int_value, bool_value, text_value, ok)
            call emit_text('case_id', trim(case_id))
            call emit_text('request_token', token)
            if (.not. ok) then
                call emit_text('status', 'invalid_input')
            else if (trim(case_id) /= 'environment') then
                call emit_text('status', 'unsupported_case')
            else
                call emit_text('status', 'ok')
                call emit_real('probe_real', real_value)
                call emit_int('probe_int', int_value)
                call emit_bool('probe_bool', bool_value)
                call emit_text('probe_text', text_value)
            end if
        end if
    end if
end program case_driver
