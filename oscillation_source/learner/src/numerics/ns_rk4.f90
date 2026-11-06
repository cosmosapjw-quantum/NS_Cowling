module ns_rk4
    use ns_precision, only: wp
    use ns_status
    use ns_numerics_types, only: rhs_callback,step_result_t
    implicit none
contains
    subroutine euler_step(rhs,x,y,h,result)
        procedure(rhs_callback)::rhs
        real(wp),intent(in)::x,y(:),h
        type(step_result_t),intent(out)::result
        result=step_result_t()
        call fail_status(result%status,NS_NOT_IMPLEMENTED,'implementation','implement Euler: lab ch02')
    end subroutine
    subroutine rk4_step(rhs,x,y,h,result)
        procedure(rhs_callback)::rhs
        real(wp),intent(in)::x,y(:),h
        type(step_result_t),intent(out)::result
        result=step_result_t()
        call fail_status(result%status,NS_NOT_IMPLEMENTED,'implementation','implement RK4: lab ch02')
    end subroutine
end module
