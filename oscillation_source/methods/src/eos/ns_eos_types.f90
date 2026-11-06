! Shared EOS wire semantics. This module contains no EOS equation or lookup.
module ns_eos_types
    use ns_precision, only: wp
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    private
    integer, parameter, public :: SIDE_REGULAR=0, SIDE_INNER=1, SIDE_OUTER=-1
    type, public :: eos_state_t
        real(wp) :: pressure_pa=0._wp, energy_j_m3=0._wp, n_b_m3=0._wp, mu_j=0._wp, chi=0._wp
        integer :: side=SIDE_REGULAR
        logical :: valid=.false.
    end type
    public :: clear_eos_state
contains
    pure subroutine clear_eos_state(state)
        type(eos_state_t), intent(out) :: state
        real(wp) :: nan
        nan=ieee_value(0._wp,ieee_quiet_nan)
        state=eos_state_t(nan,nan,nan,nan,nan,SIDE_REGULAR,.false.)
    end subroutine
end module
