module ns_bisection
    use ns_precision, only: wp
    use ns_status
    use ns_numerics_types, only: scalar_callback, root_options_t, root_result_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    private
    public :: bisect
contains
    ! Continuity is a caller assumption. No multiplication of residual signs.
    recursive subroutine bisect(f,a,b,options,result)
        procedure(scalar_callback) :: f
        real(wp), intent(in) :: a,b
        type(root_options_t), intent(in) :: options
        type(root_result_t), intent(out) :: result
        real(wp) :: lo,hi,flo,fhi,mid,fmid,half,tol,relative
        integer :: i
        result=root_result_t()
        call fail_status(result%status,NS_INVALID_INPUT,'numerics','invalid bisection options or interval')
        if (.not.all(ieee_is_finite([a,b,options%xtol,options%rtol,options%ftol]))) return
        if (a>b .or. min(options%xtol,options%rtol,options%ftol)<0._wp .or. options%max_iter<=0) return
        lo=a;hi=b;result%lo=a;result%hi=b
        call evaluate(lo,flo); if (.not.is_ok(result%status)) return
        result%f_lo=flo
        if (flo==0._wp) then
            call finish(lo,flo,'exact');return
        end if
        call evaluate(hi,fhi); if (.not.is_ok(result%status)) return
        result%f_hi=fhi
        if (fhi==0._wp) then
            call finish(hi,fhi,'exact');return
        end if
        if (same_sign(flo,fhi)) then
            call fail_status(result%status,NS_NO_BRACKET,'numerics','endpoint residuals have the same sign');return
        end if
        if (options%ftol>0._wp) then
            if (abs(flo)<=options%ftol) then
                call finish(lo,flo,'residual');return
            else if (abs(fhi)<=options%ftol) then
                call finish(hi,fhi,'residual');return
            end if
        end if
        do i=1,options%max_iter
            ! Opposite-sign endpoints: splitting each operand avoids b-a overflow.
            if (lo<0._wp .and. hi>0._wp) then
                mid=lo/2._wp+hi/2._wp; half=hi/2._wp-lo/2._wp
            else
                half=(hi-lo)/2._wp;mid=lo+half
            end if
            relative=options%rtol*max(abs(lo),abs(hi))
            if (.not.ieee_is_finite(relative)) relative=huge(relative)
            tol=options%xtol+min(relative,huge(relative)-options%xtol)
            if (mid==lo .or. mid==hi) then
                if (half<=tol .and. tol>0._wp) then
                    if (abs(flo)<=abs(fhi)) then
                        call finish(lo,flo,'coordinate')
                    else
                        call finish(hi,fhi,'coordinate')
                    end if
                else
                    call fail_status(result%status,NS_STAGNATION,'numerics','no representable interior point')
                end if
                return
            end if
            result%iterations=i
            call evaluate(mid,fmid);if (.not.is_ok(result%status)) return
            if (fmid==0._wp) then
                call finish(mid,fmid,'exact');return
            else if (options%ftol>0._wp .and. abs(fmid)<=options%ftol) then
                call finish(mid,fmid,'residual');return
            else if (half<=tol .and. tol>0._wp) then
                call finish(mid,fmid,'coordinate');return
            end if
            if (same_sign(flo,fmid)) then
                lo=mid;flo=fmid
            else
                hi=mid;fhi=fmid
            end if
            result%lo=lo;result%hi=hi;result%f_lo=flo;result%f_hi=fhi
        end do
        call fail_status(result%status,NS_MAX_ITERATIONS,'numerics','bisection iteration limit; bracket retained')
    contains
        recursive subroutine evaluate(x,value)
            real(wp), intent(in) :: x
            real(wp), intent(out) :: value
            value=ieee_value(0._wp,ieee_quiet_nan)
            result%status=status_t()
            call f(x,value,result%status);result%evaluations=result%evaluations+1
            if (.not.is_ok(result%status)) return
            if (.not.ieee_is_finite(value)) &
                call fail_status(result%status,NS_NONFINITE_RESULT,'numerics','nonfinite scalar callback')
        end subroutine
        subroutine finish(x,value,reason)
            real(wp), intent(in) :: x,value
            character(len=*), intent(in) :: reason
            result%root=x;result%residual=value;result%valid=.true.;result%termination=reason
            if (reason=='exact') then
                result%lo=x;result%hi=x;result%f_lo=value;result%f_hi=value
            else
                result%lo=lo;result%hi=hi;result%f_lo=flo;result%f_hi=fhi
            end if
            call clear_status(result%status)
        end subroutine
    end subroutine
    pure logical function same_sign(a,b)
        real(wp),intent(in)::a,b
        same_sign=(a>0._wp .and. b>0._wp) .or. (a<0._wp .and. b<0._wp)
    end function
end module ns_bisection
