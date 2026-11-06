module ns_surface_residual
    use ns_precision, only: wp
    use ns_status
    implicit none
    private
    public :: cowling_surface_residual
contains
    subroutine cowling_surface_residual(r,omega,nu,lambda,dnu,w,v,residual,q_over_r2,status)
        real(wp),intent(in)::r,omega,nu,lambda,dnu,w,v
        real(wp),intent(out)::residual,q_over_r2
        type(status_t),intent(out)::status
        residual=0._wp
        q_over_r2=0._wp
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','T07 learner surface residual not implemented')
    end subroutine cowling_surface_residual
end module ns_surface_residual
