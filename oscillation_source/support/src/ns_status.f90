! Explicit status: outputs of a failed operation are not physical results.
module ns_status
    implicit none
    private
    integer, parameter, public :: NS_OK = 0, NS_UNINITIALIZED = 1, NS_NOT_IMPLEMENTED = 2
    integer, parameter, public :: NS_INVALID_INPUT = 3, NS_NONFINITE_INPUT = 4
    integer, parameter, public :: NS_OUT_OF_DOMAIN = 5, NS_RANGE_ERROR = 6
    integer, parameter, public :: NS_NO_BRACKET = 7, NS_IDENTITY_UNRESOLVED = 8
    integer, parameter, public :: NS_MAX_ITERATIONS = 9, NS_STAGNATION = 10
    integer, parameter, public :: NS_STEP_TOO_SMALL = 11, NS_MAX_STEPS = 12
    integer, parameter, public :: NS_NONFINITE_RESULT = 13
    integer, parameter, public :: NS_MATCHING_INCONSISTENT = 14
    type, public :: status_t
        integer :: code = NS_UNINITIALIZED
        character(len=24) :: category = 'implementation'
        character(len=256) :: detail = 'status not initialized'
    end type status_t
    public :: clear_status, fail_status, is_ok, status_name
contains
    pure subroutine clear_status(status)
        type(status_t), intent(out) :: status
        status%code = NS_OK
        status%category = 'none'
        status%detail = ''
    end subroutine clear_status

    pure subroutine fail_status(status, code, category, detail)
        type(status_t), intent(out) :: status
        integer, intent(in) :: code
        character(len=*), intent(in) :: category, detail
        status%code = code
        status%category = category
        status%detail = detail
        ! A failure must never accidentally carry a successful status.
        if (code == NS_OK) then
            status%code = NS_INVALID_INPUT
            status%category = 'implementation'
            status%detail = 'fail_status called with NS_OK'
        end if
    end subroutine fail_status

    pure logical function is_ok(status)
        type(status_t), intent(in) :: status
        is_ok = status%code == NS_OK
    end function is_ok

    pure function status_name(code) result(name)
        integer, intent(in) :: code
        character(len=32) :: name
        select case(code)
        case(NS_OK); name = 'ok'
        case(NS_UNINITIALIZED); name = 'uninitialized'
        case(NS_NOT_IMPLEMENTED); name = 'not_implemented'
        case(NS_INVALID_INPUT); name = 'invalid_input'
        case(NS_NONFINITE_INPUT); name = 'nonfinite_input'
        case(NS_OUT_OF_DOMAIN); name = 'out_of_domain'
        case(NS_RANGE_ERROR); name = 'range_error'
        case(NS_NO_BRACKET); name = 'no_bracket'
        case(NS_IDENTITY_UNRESOLVED); name = 'identity_unresolved'
        case(NS_MAX_ITERATIONS); name = 'max_iterations'
        case(NS_STAGNATION); name = 'stagnation'
        case(NS_STEP_TOO_SMALL); name = 'step_too_small'
        case(NS_MAX_STEPS); name = 'max_steps'
        case(NS_NONFINITE_RESULT); name = 'nonfinite_result'
        case(NS_MATCHING_INCONSISTENT); name = 'matching_inconsistent'
        case default; name = 'unknown_status'
        end select
    end function status_name
end module ns_status
