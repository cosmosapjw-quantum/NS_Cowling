! Small numerical building blocks, with no EOS or stellar solver hidden here.
module ns_safe_math
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use ns_precision, only: wp
    use ns_status, only: status_t, clear_status, fail_status, NS_INVALID_INPUT, NS_NONFINITE_INPUT, NS_RANGE_ERROR
    implicit none
    private
    public :: safe_norm2, checked_scale
contains
    pure subroutine safe_norm2(vector, norm, status)
        real(wp), intent(in) :: vector(:)
        real(wp), intent(out) :: norm
        type(status_t), intent(out) :: status
        real(wp) :: scale, factor
        norm = ieee_value(0.0_wp, ieee_quiet_nan)
        if (size(vector) == 0) then
            call fail_status(status, NS_INVALID_INPUT, 'numerics', 'safe_norm2 requires a nonempty vector')
            return
        end if
        if (.not. all(ieee_is_finite(vector))) then
            call fail_status(status, NS_NONFINITE_INPUT, 'numerics', 'safe_norm2 input is not finite')
            return
        end if
        scale = maxval(abs(vector))
        if (scale <= 0.0_wp) then
            norm = 0.0_wp
            call clear_status(status)
            return
        end if
        ! At least one normalized component has magnitude one. The sum cannot
        ! underflow; tiny relative components may round to zero harmlessly.
        factor = sqrt(sum((vector/scale)**2))
        if (scale > huge(norm)/factor) then
            call fail_status(status, NS_RANGE_ERROR, 'numerics', 'Euclidean norm exceeds real64 range')
            return
        end if
        norm = scale*factor
        call clear_status(status)
    end subroutine safe_norm2

    pure subroutine checked_scale(value, factor, result_value, status)
        real(wp), intent(in) :: value, factor
        real(wp), intent(out) :: result_value
        type(status_t), intent(out) :: status
        result_value = ieee_value(0.0_wp, ieee_quiet_nan)
        if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(factor)) then
            call fail_status(status, NS_NONFINITE_INPUT, 'numerics', 'conversion input or factor is not finite')
            return
        end if
        if (factor <= 0.0_wp) then
            call fail_status(status, NS_INVALID_INPUT, 'units', 'conversion scale must be positive')
            return
        end if
        if (factor > 1.0_wp) then
            if (abs(value) > huge(value)/factor) then
                call fail_status(status, NS_RANGE_ERROR, 'numerics', 'conversion would overflow real64')
                return
            end if
        end if
        result_value = value*factor
        if (.not. ieee_is_finite(result_value)) then
            result_value = ieee_value(0.0_wp, ieee_quiet_nan)
            call fail_status(status, NS_RANGE_ERROR, 'numerics', 'conversion result is not finite')
            return
        end if
        if (abs(value) > 0.0_wp .and. abs(result_value) <= 0.0_wp) then
            result_value = ieee_value(0.0_wp, ieee_quiet_nan)
            call fail_status(status, NS_RANGE_ERROR, 'numerics', 'nonzero input converted to zero by underflow')
            return
        end if
        call clear_status(status)
    end subroutine checked_scale
end module ns_safe_math
