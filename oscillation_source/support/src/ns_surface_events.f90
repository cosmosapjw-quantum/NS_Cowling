! Monotone pressure-coordinate adapter, not EOS extrapolation or a universal event solver.
module ns_surface_events
    use ns_precision, only: wp
    use ns_status
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: pressure_span,pressure_at_q,pressure_stop_coordinates
contains
    subroutine pressure_span(p_start,p_cut,q_end,status)
        real(wp),intent(in)::p_start,p_cut
        real(wp),intent(out)::q_end
        type(status_t),intent(out)::status
        q_end=0._wp
        call fail_status(status,NS_INVALID_INPUT,'physics','require finite 0<P_cut<P_start')
        if(.not.all(ieee_is_finite([p_start,p_cut]))) return
        if(p_cut<=0._wp.or.p_start<=p_cut) return
        q_end=log(p_start)-log(p_cut)
        if(q_end<=0._wp) then
            call fail_status(status,NS_STAGNATION,'numerics','log pressure endpoints indistinguishable');return
        end if
        call clear_status(status)
    end subroutine
    subroutine pressure_at_q(q,p_start,p_cut,q_end,p,status)
        real(wp),intent(in)::q,p_start,p_cut,q_end
        real(wp),intent(out)::p
        type(status_t),intent(out)::status
        p=0._wp
        call fail_status(status,NS_OUT_OF_DOMAIN,'numerics','pressure coordinate outside closed integration interval')
        if(.not.ieee_is_finite(q)) return
        if(q<0._wp.or.q>q_end) return
        if(q==0._wp) then
            p=p_start
        else if(q==q_end) then
            p=p_cut
        else
            p=exp(log(p_start)-q)
        end if
        if(.not.ieee_is_finite(p).or.p<=0._wp) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','pressure not representable');return
        end if
        call clear_status(status)
    end subroutine
    subroutine pressure_stop_coordinates(events,p_start,p_cut,stops,status)
        real(wp),intent(in)::events(:),p_start,p_cut
        real(wp),allocatable,intent(out)::stops(:)
        type(status_t),intent(out)::status
        integer::i
        allocate(stops(size(events)))
        call fail_status(status,NS_INVALID_INPUT,'physics','smooth pressure stops must strictly decrease inside domain')
        if(.not.all(ieee_is_finite(events))) return
        do i=1,size(events)
            if(events(i)>=p_start.or.events(i)<=p_cut) return
            if(i>1) then
                if(events(i)>=events(i-1)) return
            end if
            stops(i)=log(p_start)-log(events(i))
            if(stops(i)<=0._wp) return
            if(i>1) then
                if(stops(i)<=stops(i-1)) return
            end if
        end do
        call clear_status(status)
    end subroutine
end module
