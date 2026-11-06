! Shared contracts only: algorithms remain in the selected workspace.
module ns_numerics_types
    use ns_precision, only: wp
    use ns_status
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    private
    type, public :: root_options_t
        real(wp) :: xtol=1e-12_wp, rtol=0._wp, ftol=0._wp
        integer :: max_iter=100
    end type
    type, public :: root_result_t
        real(wp) :: root=0._wp, residual=0._wp, lo=0._wp, hi=0._wp
        real(wp) :: f_lo=0._wp, f_hi=0._wp
        integer :: iterations=0, evaluations=0
        logical :: valid=.false.
        character(len=24) :: termination='none'
        type(status_t) :: status
    end type
    type, public :: step_result_t
        real(wp), allocatable :: y(:), error(:)
        integer :: evaluations=0
        type(status_t) :: status
    end type
    type, public :: ode_options_t
        real(wp), allocatable :: atol(:)
        real(wp) :: rtol=1e-9_wp, h_initial=0.1_wp, h_min=1e-14_wp, h_max=1._wp
        integer :: max_steps=10000
        logical :: record_trials=.false.
    end type
    type, public :: ode_result_t
        ! x/y contain only initial and ACCEPTED states; failure keeps a prefix.
        real(wp), allocatable :: x(:), y(:,:)
        real(wp), allocatable :: trial_x(:), trial_h(:), trial_error(:), trial_y(:,:)
        logical, allocatable :: trial_accepted(:)
        integer :: accepted=0, rejected=0, attempts=0, evaluations=0, stops_hit=0
        logical :: valid=.false.
        type(status_t) :: status
    end type
    abstract interface
        subroutine scalar_callback(x, value, status)
            import wp, status_t
            real(wp), intent(in) :: x
            real(wp), intent(out) :: value
            type(status_t), intent(out) :: status
        end subroutine
        subroutine rhs_callback(x, y, dydx, status)
            import wp, status_t
            real(wp), intent(in) :: x, y(:)
            real(wp), intent(out) :: dydx(:)
            type(status_t), intent(out) :: status
        end subroutine
    end interface
    public :: scalar_callback, rhs_callback, evaluate_rhs, validate_step
contains
    subroutine evaluate_rhs(rhs, x, y, dydx, evaluations, status)
        procedure(rhs_callback) :: rhs
        real(wp), intent(in) :: x, y(:)
        real(wp), intent(out) :: dydx(:)
        integer, intent(inout) :: evaluations
        type(status_t), intent(out) :: status
        dydx=ieee_value(0._wp,ieee_quiet_nan)
        if (.not.ieee_is_finite(x) .or. .not.all(ieee_is_finite(y))) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','nonfinite stage state');return
        end if
        status=status_t()
        call rhs(x,y,dydx,status)
        evaluations=evaluations+1
        if (.not.is_ok(status)) return
        if (.not.all(ieee_is_finite(dydx))) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','RHS returned nonfinite derivative');return
        end if
    end subroutine
    subroutine validate_step(x,y,h,status)
        real(wp), intent(in) :: x,y(:),h
        type(status_t), intent(out) :: status
        call fail_status(status,NS_INVALID_INPUT,'numerics','step requires a finite nonempty state and nonzero h')
        if (size(y)==0) return
        if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(h)) return
        if (.not.all(ieee_is_finite(y)) .or. h==0._wp) return
        if (.not.ieee_is_finite(x+h)) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','x+h overflow');return
        end if
        if (x+h==x) then
            call fail_status(status,NS_STAGNATION,'numerics','x+h equals x in working precision');return
        end if
        call clear_status(status)
    end subroutine
end module ns_numerics_types
