module ns_adaptive_rk
    use ns_precision, only: wp
    use ns_status
    use ns_numerics_types, only: rhs_callback,ode_options_t,ode_result_t
    implicit none
contains
    subroutine scaled_error_norm(old,new,error,atol,rtol,value,status)
        real(wp),intent(in)::old(:),new(:),error(:),atol(:),rtol
        real(wp),intent(out)::value
        type(status_t),intent(out)::status
        value=0._wp
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','implement error norm: lab ch02')
    end subroutine
    subroutine step_factor(error_norm,previous_rejected,factor,status)
        real(wp),intent(in)::error_norm
        logical,intent(in)::previous_rejected
        real(wp),intent(out)::factor
        type(status_t),intent(out)::status
        factor=0._wp
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','implement step controller: lab ch02')
    end subroutine
    subroutine adaptive_integrate(rhs,x0,x1,y0,options,result,event_schedule)
        procedure(rhs_callback)::rhs
        real(wp),intent(in)::x0,x1,y0(:)
        type(ode_options_t),intent(in)::options
        type(ode_result_t),intent(out)::result
        real(wp),intent(in),optional::event_schedule(:)
        result=ode_result_t()
        call fail_status(result%status,NS_NOT_IMPLEMENTED,'implementation','implement adaptive accept/reject: lab ch02')
    end subroutine
end module
