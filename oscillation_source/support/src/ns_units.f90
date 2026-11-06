! Unit conversion boundaries. Signed scalar conversions do not admit an EOS.
! Failure outputs are NaN and must not be consumed; the driver emits no value.
module ns_units
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use ns_precision, only: wp
    use ns_physical_constants, only: c_si, two_pi, km_in_m, mev_fm3_in_si, &
                                    inverse_fm3_in_si, nominal_solar_mass_kg
    use ns_safe_math, only: checked_scale
    use ns_status, only: status_t, fail_status, NS_INVALID_INPUT, NS_NONFINITE_INPUT
    implicit none
    private
    public :: to_si_energy_density, from_si_energy_density
    public :: to_si_number_density, from_si_number_density
    public :: omega_to_hz, hz_to_omega, nominal_solar_mass_to_kg, kg_to_nominal_solar_mass
    public :: kilometers_to_meters, meters_to_kilometers
    public :: angular_to_wavenumber, wavenumber_to_angular
    public :: energy_to_mass_density, number_to_rest_mass_density
contains
    pure subroutine to_si_energy_density(energy_mev_fm3, energy_j_m3, status)
        real(wp), intent(in) :: energy_mev_fm3
        real(wp), intent(out) :: energy_j_m3
        type(status_t), intent(out) :: status
        call checked_scale(energy_mev_fm3, mev_fm3_in_si, energy_j_m3, status)
    end subroutine
    pure subroutine from_si_energy_density(energy_j_m3, energy_mev_fm3, status)
        real(wp), intent(in) :: energy_j_m3
        real(wp), intent(out) :: energy_mev_fm3
        type(status_t), intent(out) :: status
        call checked_scale(energy_j_m3, 1.0_wp/mev_fm3_in_si, energy_mev_fm3, status)
    end subroutine
    pure subroutine to_si_number_density(n_b_fm3, n_b_m3, status)
        real(wp), intent(in) :: n_b_fm3
        real(wp), intent(out) :: n_b_m3
        type(status_t), intent(out) :: status
        call checked_scale(n_b_fm3, inverse_fm3_in_si, n_b_m3, status)
    end subroutine
    pure subroutine from_si_number_density(n_b_m3, n_b_fm3, status)
        real(wp), intent(in) :: n_b_m3
        real(wp), intent(out) :: n_b_fm3
        type(status_t), intent(out) :: status
        call checked_scale(n_b_m3, 1.0_wp/inverse_fm3_in_si, n_b_fm3, status)
    end subroutine
    pure subroutine omega_to_hz(omega_rad_s, frequency_hz, status)
        real(wp), intent(in) :: omega_rad_s
        real(wp), intent(out) :: frequency_hz
        type(status_t), intent(out) :: status
        call checked_scale(omega_rad_s, 1.0_wp/two_pi, frequency_hz, status)
    end subroutine
    pure subroutine hz_to_omega(frequency_hz, omega_rad_s, status)
        real(wp), intent(in) :: frequency_hz
        real(wp), intent(out) :: omega_rad_s
        type(status_t), intent(out) :: status
        call checked_scale(frequency_hz, two_pi, omega_rad_s, status)
    end subroutine
    pure subroutine nominal_solar_mass_to_kg(q_solar, mass_kg, status)
        real(wp), intent(in) :: q_solar
        real(wp), intent(out) :: mass_kg
        type(status_t), intent(out) :: status
        call checked_scale(q_solar, nominal_solar_mass_kg, mass_kg, status)
    end subroutine
    pure subroutine kg_to_nominal_solar_mass(mass_kg, q_solar, status)
        real(wp), intent(in) :: mass_kg
        real(wp), intent(out) :: q_solar
        type(status_t), intent(out) :: status
        call checked_scale(mass_kg, 1.0_wp/nominal_solar_mass_kg, q_solar, status)
    end subroutine
    pure subroutine kilometers_to_meters(radius_km, radius_m, status)
        real(wp), intent(in) :: radius_km
        real(wp), intent(out) :: radius_m
        type(status_t), intent(out) :: status
        call checked_scale(radius_km, km_in_m, radius_m, status)
    end subroutine
    pure subroutine meters_to_kilometers(radius_m, radius_km, status)
        real(wp), intent(in) :: radius_m
        real(wp), intent(out) :: radius_km
        type(status_t), intent(out) :: status
        call checked_scale(radius_m, 1.0_wp/km_in_m, radius_km, status)
    end subroutine
    pure subroutine angular_to_wavenumber(omega_rad_s, k_per_m, status)
        real(wp), intent(in) :: omega_rad_s
        real(wp), intent(out) :: k_per_m
        type(status_t), intent(out) :: status
        call checked_scale(omega_rad_s, 1.0_wp/c_si, k_per_m, status)
    end subroutine
    pure subroutine wavenumber_to_angular(k_per_m, omega_rad_s, status)
        real(wp), intent(in) :: k_per_m
        real(wp), intent(out) :: omega_rad_s
        type(status_t), intent(out) :: status
        call checked_scale(k_per_m, c_si, omega_rad_s, status)
    end subroutine
    pure subroutine energy_to_mass_density(energy_j_m3, rho_e_kg_m3, status)
        real(wp), intent(in) :: energy_j_m3
        real(wp), intent(out) :: rho_e_kg_m3
        type(status_t), intent(out) :: status
        call checked_scale(energy_j_m3, 1.0_wp/(c_si*c_si), rho_e_kg_m3, status)
    end subroutine
    pure subroutine number_to_rest_mass_density(n_b_m3, mass_b_kg, rho0_kg_m3, status)
        real(wp), intent(in) :: n_b_m3, mass_b_kg
        real(wp), intent(out) :: rho0_kg_m3
        type(status_t), intent(out) :: status
        rho0_kg_m3 = ieee_value(0.0_wp, ieee_quiet_nan)
        if (.not. ieee_is_finite(n_b_m3) .or. .not. ieee_is_finite(mass_b_kg)) then
            call fail_status(status, NS_NONFINITE_INPUT, 'units', 'baryon density and mass must be finite')
            return
        end if
        if (n_b_m3 < 0.0_wp .or. mass_b_kg <= 0.0_wp) then
            call fail_status(status, NS_INVALID_INPUT, 'units', 'n_b must be nonnegative and mass_b positive')
            return
        end if
        call checked_scale(n_b_m3, mass_b_kg, rho0_kg_m3, status)
    end subroutine
end module ns_units
