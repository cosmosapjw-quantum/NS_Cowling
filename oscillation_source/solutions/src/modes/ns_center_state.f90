module ns_center_state
    use ns_precision, only: wp
    use ns_status
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: cowling_center_state
contains
    subroutine cowling_center_state(ell,r0,amplitude,w,v,status)
        integer,intent(in)::ell
        real(wp),intent(in)::r0,amplitude
        real(wp),intent(out)::w,v
        type(status_t),intent(out)::status
        integer::i
        real(wp)::rp
        w=0._wp;v=0._wp
        call fail_status(status,NS_INVALID_INPUT,'physics','regular center requires ell>=1, r0>0, nonzero finite amplitude')
        if(ell<1.or.ell>20.or.r0<=0._wp.or.amplitude==0._wp) return
        if(.not.all(ieee_is_finite([r0,amplitude]))) return
        rp=1._wp
        do i=1,ell
            if(abs(rp)>huge(rp)/r0) then;call fail_status(status,NS_RANGE_ERROR,'numerics','center power overflow');return;end if
            rp=rp*r0
        end do
        if(abs(rp)>huge(rp)/r0.or.abs(amplitude)>huge(amplitude)/max(rp,tiny(rp))) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','center state overflow');return
        end if
        v=-amplitude*rp/real(ell,wp);w=amplitude*rp*r0
        if(.not.all(ieee_is_finite([w,v])).or.w==0._wp.or.v==0._wp) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','center state unrepresentable');return
        end if
        call clear_status(status)
    end subroutine
end module ns_center_state
