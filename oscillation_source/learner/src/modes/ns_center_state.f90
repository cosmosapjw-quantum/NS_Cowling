module ns_center_state
    use ns_precision, only: wp
    use ns_status
    implicit none
    private
    public :: cowling_center_state
contains
    subroutine cowling_center_state(ell,r0,amplitude,w,v,status)
        integer,intent(in)::ell
        real(wp),intent(in)::r0,amplitude
        real(wp),intent(out)::w,v
        type(status_t),intent(out)::status
        w=0._wp
        v=0._wp
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','T07 learner center state not implemented')
    end subroutine cowling_center_state
end module ns_center_state
