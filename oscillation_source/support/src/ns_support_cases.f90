! Public demonstration/verification adapter for PROVIDED modules only.
! It contains no student EOS, integrator, TOV, Cowling or root implementation.
module ns_support_cases
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use ns_precision, only: wp
    use ns_physical_constants
    use ns_status
    use ns_safe_math, only: safe_norm2
    use ns_units
    use ns_electron_gas, only: electron_state_t, electron_state
    use ns_bootstrap_protocol, only: emit_text, emit_int, emit_real, emit_bool
    implicit none
    private
    public :: is_support_case, run_support_case
contains
    logical function is_support_case(case_id)
        character(len=*), intent(in) :: case_id
        select case(trim(case_id))
        case('constants','units','conversion','safe_norm','status_lifecycle','electron_gas')
            is_support_case = .true.
        case default
            is_support_case = .false.
        end select
    end function

    subroutine run_support_case(case_id, path)
        character(len=*), intent(in) :: case_id, path
        character(len=32) :: request_token
        integer :: unit, ios
        type(status_t) :: status
        namelist /request/ request_token
        request_token = ''
        call emit_text('case_id', case_id)
        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios == 0) then
            read(unit, nml=request, iostat=ios)
            if (ios /= 0 .or. len_trim(request_token) /= 32) then
                close(unit)
                ios = 1
            end if
        end if
        call emit_text('request_token', request_token)
        if (ios /= 0) then
            call fail_status(status, NS_INVALID_INPUT, 'data', 'missing or malformed request')
            call emit_status(status)
            return
        end if
        select case(trim(case_id))
        case('constants')
            call emit_constants()
        case('units')
            call emit_units(unit)
        case('conversion')
            call emit_conversion(unit)
        case('safe_norm')
            call emit_norm(unit)
        case('status_lifecycle')
            call emit_lifecycle()
        case('electron_gas')
            call emit_electron(unit)
        end select
        close(unit)
    end subroutine

    subroutine emit_status(status)
        type(status_t), intent(in) :: status
        call emit_text('status', trim(status_name(status%code)))
        call emit_int('status_code', int(status%code,int64))
        call emit_text('category', trim(status%category))
        call emit_text('detail', trim(status%detail))
        call emit_bool('valid_value', is_ok(status))
        call emit_text('constants_id', constants_id)
    end subroutine

    subroutine emit_constants()
        type(status_t) :: status
        call clear_status(status)
        call emit_status(status)
        call emit_real('c_si', c_si)
        call emit_real('h_si', h_si)
        call emit_real('e_si', e_si)
        call emit_real('g_si', g_si)
        call emit_real('g_std_uncertainty', g_std_uncertainty)
        call emit_real('m_e_si', m_e_si)
        call emit_real('m_n_si', m_n_si)
        call emit_real('m_p_si', m_p_si)
        call emit_real('gm_sun_nominal', gm_sun_nominal)
        call emit_real('hbar_si', hbar_si)
        call emit_real('hbar_c_mev_fm', hbar_c_mev_fm)
        call emit_real('two_pi', two_pi)
        call emit_int('precision_bits', int(storage_size(c_si),int64))
        call emit_int('precision_digits', int(precision(c_si),int64))
        call emit_int('precision_radix', int(radix(c_si),int64))
    end subroutine

    subroutine emit_units(unit)
        integer, intent(in) :: unit
        integer :: ios
        real(wp) :: energy_mev_fm3, n_b_fm3, omega_rad_s, q_solar, radius_km
        real(wp) :: energy_j_m3, n_b_m3, frequency_hz, mass_kg, radius_m
        type(status_t) :: status
        namelist /parameters/ energy_mev_fm3, n_b_fm3, omega_rad_s, q_solar, radius_km
        energy_mev_fm3 = 0.0_wp; n_b_fm3 = 0.0_wp; omega_rad_s = 0.0_wp
        q_solar = 0.0_wp; radius_km = 0.0_wp
        read(unit, nml=parameters, iostat=ios)
        if (ios /= 0) then
            call fail_status(status, NS_INVALID_INPUT, 'data', 'invalid unit-conversion parameters')
        else
            call to_si_energy_density(energy_mev_fm3, energy_j_m3, status)
            if (is_ok(status)) call to_si_number_density(n_b_fm3, n_b_m3, status)
            if (is_ok(status)) call omega_to_hz(omega_rad_s, frequency_hz, status)
            if (is_ok(status)) call nominal_solar_mass_to_kg(q_solar, mass_kg, status)
            if (is_ok(status)) call kilometers_to_meters(radius_km, radius_m, status)
        end if
        call emit_status(status)
        if (.not. is_ok(status)) return
        call emit_real('energy_j_m3', energy_j_m3)
        call emit_real('n_b_m3', n_b_m3)
        call emit_real('frequency_hz', frequency_hz)
        call emit_real('mass_kg', mass_kg)
        call emit_real('radius_m', radius_m)
    end subroutine

    subroutine set_special(value, special, status)
        real(wp), intent(inout) :: value
        character(len=*), intent(in) :: special
        type(status_t), intent(out) :: status
        call clear_status(status)
        select case(trim(special))
        case('none')
        case('nan')
            value = ieee_value(0.0_wp, ieee_quiet_nan)
        case('infinity')
            value = ieee_value(0.0_wp, ieee_positive_inf)
        case default
            call fail_status(status, NS_INVALID_INPUT, 'data', 'unrecognized synthetic input')
        end select
    end subroutine

    subroutine emit_conversion(unit)
        integer, intent(in) :: unit
        integer :: ios
        character(len=64) :: operation, special_input
        real(wp) :: value, mass_b_kg, result_value
        type(status_t) :: status
        namelist /parameters/ operation, value, mass_b_kg, special_input
        operation = ''; value = 0.0_wp; mass_b_kg = -1.0_wp; special_input = 'none'
        read(unit, nml=parameters, iostat=ios)
        if (ios /= 0) then
            call fail_status(status, NS_INVALID_INPUT, 'data', 'invalid conversion parameters')
            call emit_status(status)
            return
        end if
        call set_special(value, special_input, status)
        if (.not. is_ok(status)) then
            call emit_status(status)
            return
        end if
        select case(trim(operation))
        case('energy_to_si'); call to_si_energy_density(value, result_value, status)
        case('energy_from_si'); call from_si_energy_density(value, result_value, status)
        case('number_to_si'); call to_si_number_density(value, result_value, status)
        case('number_from_si'); call from_si_number_density(value, result_value, status)
        case('omega_to_hz'); call omega_to_hz(value, result_value, status)
        case('hz_to_omega'); call hz_to_omega(value, result_value, status)
        case('nominal_mass_to_kg'); call nominal_solar_mass_to_kg(value, result_value, status)
        case('kg_to_nominal_mass'); call kg_to_nominal_solar_mass(value, result_value, status)
        case('km_to_m'); call kilometers_to_meters(value, result_value, status)
        case('m_to_km'); call meters_to_kilometers(value, result_value, status)
        case('angular_to_wavenumber'); call angular_to_wavenumber(value, result_value, status)
        case('wavenumber_to_angular'); call wavenumber_to_angular(value, result_value, status)
        case('energy_to_mass_density'); call energy_to_mass_density(value, result_value, status)
        case('number_to_rest_mass_density')
            call number_to_rest_mass_density(value, mass_b_kg, result_value, status)
        case default
            call fail_status(status, NS_INVALID_INPUT, 'data', 'unknown conversion operation')
        end select
        call emit_status(status)
        if (is_ok(status)) call emit_real('value_out', result_value)
    end subroutine

    subroutine emit_norm(unit)
        integer, intent(in) :: unit
        integer :: ios
        real(wp) :: x, y, z, norm, values(3)
        character(len=64) :: special_input
        logical :: empty_vector
        type(status_t) :: status
        namelist /parameters/ x, y, z, special_input, empty_vector
        x = 0.0_wp; y = 0.0_wp; z = 0.0_wp
        special_input = 'none'; empty_vector = .false.
        read(unit, nml=parameters, iostat=ios)
        if (ios /= 0) then
            call fail_status(status, NS_INVALID_INPUT, 'data', 'invalid norm parameters')
            call emit_status(status)
            return
        end if
        call set_special(x, special_input, status)
        if (.not. is_ok(status)) then
            call emit_status(status)
            return
        end if
        values = [x,y,z]
        if (empty_vector) then
            call safe_norm2(values(:0), norm, status)
        else
            call safe_norm2(values, norm, status)
        end if
        call emit_status(status)
        if (is_ok(status)) call emit_real('norm', norm)
    end subroutine

    subroutine emit_electron(unit)
        integer,intent(in)::unit
        integer::ios
        real(wp)::n_e_m3,identity_scale
        type(electron_state_t)::e
        type(status_t)::status
        namelist /parameters/ n_e_m3
        n_e_m3=0._wp
        read(unit,nml=parameters,iostat=ios)
        if(ios/=0) then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid electron-gas parameters')
        else
            call electron_state(n_e_m3,e,status)
        end if
        call emit_status(status)
        if(.not.is_ok(status)) return
        call emit_real('n_e_m3',e%n_e_m3);call emit_real('p_f_kg_m_s',e%p_f_kg_m_s)
        call emit_real('relativity_x',e%relativity_x);call emit_real('mu_e_j',e%mu_j)
        call emit_real('energy_j_m3',e%energy_j_m3);call emit_real('pressure_pa',e%pressure_pa)
        call emit_bool('used_small_x_series',e%used_small_x_series)
        identity_scale=max(abs(e%pressure_pa),abs(e%n_e_m3*e%mu_j),abs(e%energy_j_m3),tiny(1._wp))
        call emit_bool('thermo_identity_closed',abs(e%pressure_pa-(e%n_e_m3*e%mu_j-e%energy_j_m3))<=5e-11_wp*identity_scale)
    end subroutine

    subroutine emit_lifecycle()
        type(status_t) :: status, final_status
        logical :: default_ok, failure_ok, reset_ok
        integer :: reset_code
        character(len=24) :: failure_category, reset_category
        character(len=256) :: failure_detail, reset_detail
        default_ok = is_ok(status)
        call fail_status(status, NS_INVALID_INPUT, 'units', 'synthetic dimension mismatch')
        failure_ok = is_ok(status)
        failure_category = status%category; failure_detail = status%detail
        call clear_status(status)
        reset_ok = is_ok(status); reset_code = status%code
        reset_category = status%category; reset_detail = status%detail
        call clear_status(final_status)
        call emit_status(final_status)
        call emit_bool('default_ok', default_ok)
        call emit_bool('failure_ok', failure_ok)
        call emit_text('failure_category', failure_category)
        call emit_text('failure_detail', failure_detail)
        call emit_bool('reset_ok', reset_ok)
        call emit_int('reset_code', int(reset_code,int64))
        call emit_text('reset_category', reset_category)
        call emit_text('reset_detail', reset_detail)
    end subroutine
end module ns_support_cases
