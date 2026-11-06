! One numerical authority for the course, not a claim of exact floating point.
! SI definitions: BIPM; measured values: CODATA 2022; nominal GM: IAU 2015 B3.
! Metadata and standard uncertainties: fixtures/constants/codata2022.json.
module ns_physical_constants
    use ns_precision, only: wp
    implicit none
    private
    character(len=*), parameter, public :: constants_id = 'SI2019_CODATA2022_IAU2015_v1'
    real(wp), parameter, public :: pi = acos(-1.0_wp)
    real(wp), parameter, public :: two_pi = 2.0_wp*pi
    real(wp), parameter, public :: c_si = 299792458.0_wp                ! m s^-1, exact definition
    real(wp), parameter, public :: h_si = 6.62607015e-34_wp             ! J s, exact definition
    real(wp), parameter, public :: e_si = 1.602176634e-19_wp            ! C, exact definition
    real(wp), parameter, public :: hbar_si = h_si/two_pi                ! J s, derived exact quantity
    real(wp), parameter, public :: g_si = 6.67430e-11_wp                ! m^3 kg^-1 s^-2, measured
    real(wp), parameter, public :: g_std_uncertainty = 1.5e-15_wp       ! same units, NOT a solver tolerance
    real(wp), parameter, public :: m_e_si = 9.1093837139e-31_wp         ! kg, CODATA 2022
    real(wp), parameter, public :: m_n_si = 1.67492750056e-27_wp         ! kg, CODATA 2022
    real(wp), parameter, public :: m_p_si = 1.67262192595e-27_wp         ! kg, CODATA 2022
    real(wp), parameter, public :: gm_sun_nominal = 1.3271244e20_wp     ! m^3 s^-2, nominal exact convention
    real(wp), parameter, public :: fm_in_m = 1.0e-15_wp                 ! exact prefix definition
    real(wp), parameter, public :: km_in_m = 1.0e3_wp                   ! exact prefix definition
    real(wp), parameter, public :: ev_in_j = e_si                       ! e times one volt
    real(wp), parameter, public :: mev_in_j = 1.0e6_wp*ev_in_j
    real(wp), parameter, public :: mev_fm3_in_si = mev_in_j/fm_in_m**3  ! J m^-3 per MeV fm^-3
    real(wp), parameter, public :: inverse_fm3_in_si = 1.0_wp/fm_in_m**3 ! m^-3 per fm^-3
    real(wp), parameter, public :: hbar_c_mev_fm = hbar_si*c_si/(mev_in_j*fm_in_m)
    ! This mass unit uses the pinned G above. It is not an exact measured solar mass.
    real(wp), parameter, public :: nominal_solar_mass_kg = gm_sun_nominal/g_si
end module ns_physical_constants
