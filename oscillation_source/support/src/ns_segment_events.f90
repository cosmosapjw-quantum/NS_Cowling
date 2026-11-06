module ns_segment_events
    use ns_precision, only: wp
    use ns_status
    use ns_eos_types, only: SIDE_REGULAR,SIDE_INNER,SIDE_OUTER
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    integer,parameter,public :: EVENT_MATERIAL=1, EVENT_CHI_KNOT=2
    type,public :: segment_event_t
        real(wp) :: pressure_pa=0._wp,q=0._wp
        integer :: kind=0,inner_side=SIDE_REGULAR,outer_side=SIDE_REGULAR
        logical :: ready=.false.
    end type
    public :: make_segment_event,validate_segment_events,event_stops
contains
    subroutine make_segment_event(p_start,p_cut,p_event,kind,inner_side,outer_side,event,status)
        real(wp),intent(in)::p_start,p_cut,p_event
        integer,intent(in)::kind,inner_side,outer_side
        type(segment_event_t),intent(out)::event
        type(status_t),intent(out)::status
        event=segment_event_t();call fail_status(status,NS_INVALID_INPUT,'physics','invalid segment event')
        if(.not.all(ieee_is_finite([p_start,p_cut,p_event]))) return
        if(p_cut<=0._wp.or.p_event<=p_cut.or.p_event>=p_start) return
        select case(kind)
        case(EVENT_MATERIAL)
            if(inner_side/=SIDE_INNER.or.outer_side/=SIDE_OUTER) return
        case(EVENT_CHI_KNOT)
            if(inner_side/=outer_side) return
        case default;return
        end select
        event%pressure_pa=p_event;event%q=log(p_start)-log(p_event)
        event%kind=kind;event%inner_side=inner_side;event%outer_side=outer_side
        if(.not.ieee_is_finite(event%q).or.event%q<=0._wp) return
        event%ready=.true.;call clear_status(status)
    end subroutine
    subroutine validate_segment_events(events,status)
        type(segment_event_t),intent(in)::events(:)
        type(status_t),intent(out)::status
        integer::i
        call fail_status(status,NS_INVALID_INPUT,'physics', &
            'events must be ready and ordered outward by decreasing pressure')
        do i=1,size(events)
            if(.not.events(i)%ready) return
            if(i>1) then
                if(events(i)%pressure_pa>=events(i-1)%pressure_pa.or.events(i)%q<=events(i-1)%q) return
            end if
        end do
        call clear_status(status)
    end subroutine
    subroutine event_stops(events,stops,status)
        type(segment_event_t),intent(in)::events(:)
        real(wp),allocatable,intent(out)::stops(:)
        type(status_t),intent(out)::status
        integer::i
        call validate_segment_events(events,status);if(.not.is_ok(status))then;allocate(stops(0));return;end if
        allocate(stops(size(events)))
        do i=1,size(events);stops(i)=events(i)%q;end do
    end subroutine
end module ns_segment_events
