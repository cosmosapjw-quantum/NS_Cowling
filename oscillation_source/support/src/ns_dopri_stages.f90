! Provided Dormand-Prince 5(4) stages. Control/acceptance stays in learner.
! Coefficients: Dormand & Prince (1980); Hairer/Wanner DOPRI5, CDOPRI.
! This educational implementation has no FSAL cache or dense output.
module ns_dopri_stages
    use ns_precision, only: wp
    use ns_status
    use ns_numerics_types, only: rhs_callback, step_result_t, evaluate_rhs, validate_step
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: dopri_trial
contains
    subroutine dopri_trial(rhs,x,y,h,result)
        procedure(rhs_callback) :: rhs
        real(wp), intent(in) :: x,y(:),h
        type(step_result_t), intent(out) :: result
        real(wp) :: a(7,7),c(7),b(7),d(7),k(size(y),7),trial(size(y)),e(size(y)),xs
        integer :: i,j
        result=step_result_t()
        call validate_step(x,y,h,result%status);if (.not.is_ok(result%status)) return
        a=0._wp
        c=[0._wp,1._wp/5,3._wp/10,4._wp/5,8._wp/9,1._wp,1._wp]
        a(2,1)=1._wp/5
        a(3,1:2)=[3._wp/40,9._wp/40]
        a(4,1:3)=[44._wp/45,-56._wp/15,32._wp/9]
        a(5,1:4)=[19372._wp/6561,-25360._wp/2187,64448._wp/6561,-212._wp/729]
        a(6,1:5)=[9017._wp/3168,-355._wp/33,46732._wp/5247,49._wp/176,-5103._wp/18656]
        a(7,1:6)=[35._wp/384,0._wp,500._wp/1113,125._wp/192,-2187._wp/6784,11._wp/84]
        b=[a(7,1:6),0._wp]
        ! Direct embedded difference avoids subtracting two large y0 values.
        d=[71._wp/57600,0._wp,-71._wp/16695,71._wp/1920,-17253._wp/339200,22._wp/525,-1._wp/40]
        do i=1,7
            trial=0._wp
            do j=1,i-1
                trial=trial+a(i,j)*k(:,j)
            end do
            trial=y+h*trial
            xs=x+c(i)*h
            if (i>=6) xs=x+h
            call evaluate_rhs(rhs,xs,trial,k(:,i),result%evaluations,result%status)
            if (.not.is_ok(result%status)) return
        end do
        trial=0._wp;e=0._wp
        do i=1,7
            trial=trial+b(i)*k(:,i);e=e+d(i)*k(:,i)
        end do
        trial=y+h*trial;e=h*e
        if (.not.all(ieee_is_finite(trial)) .or. .not.all(ieee_is_finite(e))) then
            call fail_status(result%status,NS_NONFINITE_RESULT,'numerics','DOPRI result/error nonfinite');return
        end if
        result%y=trial;result%error=e
    end subroutine
end module ns_dopri_stages
