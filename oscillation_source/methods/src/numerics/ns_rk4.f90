module ns_rk4
    use ns_precision, only: wp
    use ns_status
    use ns_numerics_types, only: rhs_callback, step_result_t, evaluate_rhs, validate_step
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: euler_step, rk4_step
contains
    subroutine euler_step(rhs,x,y,h,result)
        procedure(rhs_callback) :: rhs
        real(wp), intent(in) :: x,y(:),h
        type(step_result_t), intent(out) :: result
        real(wp) :: k(size(y)),trial(size(y))
        result=step_result_t()
        call validate_step(x,y,h,result%status);if (.not.is_ok(result%status)) return
        call evaluate_rhs(rhs,x,y,k,result%evaluations,result%status)
        if (.not.is_ok(result%status)) return
        trial=y+h*k
        if (.not.all(ieee_is_finite(trial))) then
            call fail_status(result%status,NS_NONFINITE_RESULT,'numerics','Euler result nonfinite');return
        end if
        result%y=trial
    end subroutine
    subroutine rk4_step(rhs,x,y,h,result)
        procedure(rhs_callback) :: rhs
        real(wp), intent(in) :: x,y(:),h
        type(step_result_t), intent(out) :: result
        real(wp) :: k(size(y),4),trial(size(y))
        result=step_result_t()
        call validate_step(x,y,h,result%status);if (.not.is_ok(result%status)) return
        call evaluate_rhs(rhs,x,y,k(:,1),result%evaluations,result%status)
        if (.not.is_ok(result%status)) return
        call evaluate_rhs(rhs,x+h/2._wp,y+(h/2._wp)*k(:,1),k(:,2),result%evaluations,result%status)
        if (.not.is_ok(result%status)) return
        call evaluate_rhs(rhs,x+h/2._wp,y+(h/2._wp)*k(:,2),k(:,3),result%evaluations,result%status)
        if (.not.is_ok(result%status)) return
        call evaluate_rhs(rhs,x+h,y+h*k(:,3),k(:,4),result%evaluations,result%status)
        if (.not.is_ok(result%status)) return
        trial=y+h*(k(:,1)/6._wp+k(:,2)/3._wp+k(:,3)/3._wp+k(:,4)/6._wp)
        if (.not.all(ieee_is_finite(trial))) then
            call fail_status(result%status,NS_NONFINITE_RESULT,'numerics','RK4 result nonfinite');return
        end if
        result%y=trial
    end subroutine
end module ns_rk4
