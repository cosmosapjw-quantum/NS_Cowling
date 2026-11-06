! T14 thermodynamically coupled table representation.
! Each native interval is a cubic Hermite P(mu) with the raw baryon density as dP/dmu.
module ns_table_eos
    use ns_precision, only: wp
    use ns_status
    use ns_eos_types
    use ns_matching_types, only: mu_bridge_t
    use ns_mu_bridge, only: construct_bridge_endpoints,evaluate_mu_bridge
    use ns_table_io, only: normalized_table_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    character(len=*), parameter, public :: TABLE_INTERPOLATION_ID='cubic_Hermite_P_of_mu_raw_n_v1'
    type, public :: table_audit_t
        integer :: row_count=0,interval_count=0,admitted_intervals=0,rejected_intervals=0
        integer :: first_rejected_interval=0
        logical :: mu_strict_increasing=.false.,pressure_strict_increasing=.false.
        logical :: fully_admitted=.false.
        real(wp) :: raw_energy_closure_max_rel=0._wp
        real(wp) :: raw_trapezoid_first_law_max_rel=0._wp
        real(wp) :: first_bad_min_d2p=0._wp,first_bad_min_causal_margin=0._wp
    end type table_audit_t
    type, public :: table_eos_t
        type(normalized_table_t) :: raw
        type(mu_bridge_t), allocatable :: cells(:)
        logical :: ready=.false.
    end type table_eos_t
    public :: audit_table_eos,init_table_eos,eos_table_at_mu,eos_table_at_pressure
contains
    subroutine audit_table_eos(raw,audit,status)
        type(normalized_table_t),intent(in)::raw
        type(table_audit_t),intent(out)::audit
        type(status_t),intent(out)::status
        type(mu_bridge_t)::cell
        type(status_t)::st
        real(wp)::pred,scale,trap,dp
        integer::i
        audit=table_audit_t()
        if(.not.raw%ready .or. raw%n<2) then
            call fail_status(status,NS_UNINITIALIZED,'implementation','normalized table is not loaded');return
        end if
        audit%row_count=raw%n;audit%interval_count=raw%n-1
        audit%mu_strict_increasing=.true.;audit%pressure_strict_increasing=.true.
        do i=1,raw%n
            pred=raw%mu_j(i)*raw%n_b_m3(i)-raw%pressure_pa(i)
            ! Independent raw-row diagnostic: relative mismatch against the supplied
            ! total energy density, not against a fitted/interpolated surrogate.
            scale=max(abs(raw%energy_j_m3(i)),tiny(1._wp))
            audit%raw_energy_closure_max_rel=max(audit%raw_energy_closure_max_rel, &
                                                 abs(pred-raw%energy_j_m3(i))/scale)
        end do
        do i=1,raw%n-1
            if(raw%mu_j(i+1)<=raw%mu_j(i)) audit%mu_strict_increasing=.false.
            if(raw%pressure_pa(i+1)<=raw%pressure_pa(i)) audit%pressure_strict_increasing=.false.
            trap=0.5_wp*(raw%n_b_m3(i)+raw%n_b_m3(i+1))*(raw%mu_j(i+1)-raw%mu_j(i))
            dp=raw%pressure_pa(i+1)-raw%pressure_pa(i)
            scale=max(abs(trap),abs(dp),tiny(1._wp))
            audit%raw_trapezoid_first_law_max_rel=max(audit%raw_trapezoid_first_law_max_rel,abs(trap-dp)/scale)
            cell=mu_bridge_t()
            if(raw%mu_j(i+1)>raw%mu_j(i) .and. raw%pressure_pa(i+1)>raw%pressure_pa(i)) then
                call construct_bridge_endpoints(raw%mu_j(i),raw%pressure_pa(i),raw%n_b_m3(i), &
                                                raw%mu_j(i+1),raw%pressure_pa(i+1),raw%n_b_m3(i+1),cell,st)
            else
                call fail_status(st,NS_MATCHING_INCONSISTENT,'physics','table knot ordering is not strictly increasing')
            end if
            if(is_ok(st)) then
                audit%admitted_intervals=audit%admitted_intervals+1
            else
                audit%rejected_intervals=audit%rejected_intervals+1
                if(audit%first_rejected_interval==0) then
                    audit%first_rejected_interval=i
                    audit%first_bad_min_d2p=cell%min_d2p
                    audit%first_bad_min_causal_margin=cell%min_causal_margin
                end if
            end if
        end do
        audit%fully_admitted=audit%mu_strict_increasing .and. audit%pressure_strict_increasing .and. &
                             audit%rejected_intervals==0
        call clear_status(status)
    end subroutine audit_table_eos

    subroutine init_table_eos(raw,model,audit,status)
        type(normalized_table_t),intent(in)::raw
        type(table_eos_t),intent(out)::model
        type(table_audit_t),intent(out)::audit
        type(status_t),intent(out)::status
        integer::i
        model=table_eos_t()
        call audit_table_eos(raw,audit,status)
        if(.not.is_ok(status)) return
        if(.not.audit%fully_admitted) then
            call fail_status(status,NS_MATCHING_INCONSISTENT,'physics', &
                'full table is not admissible under native-knot cubic Hermite P(mu)')
            return
        end if
        model%raw=raw
        allocate(model%cells(raw%n-1))
        do i=1,raw%n-1
            call construct_bridge_endpoints(raw%mu_j(i),raw%pressure_pa(i),raw%n_b_m3(i), &
                                            raw%mu_j(i+1),raw%pressure_pa(i+1),raw%n_b_m3(i+1), &
                                            model%cells(i),status)
            if(.not.is_ok(status)) then
                model=table_eos_t();return
            end if
        end do
        model%ready=.true.;call clear_status(status)
    end subroutine init_table_eos

    subroutine eos_table_at_mu(model,mu_j,side,state,dn_dmu,status)
        type(table_eos_t),intent(in)::model
        real(wp),intent(in)::mu_j
        integer,intent(in)::side
        type(eos_state_t),intent(out)::state
        real(wp),intent(out)::dn_dmu
        type(status_t),intent(out)::status
        integer::i
        call clear_eos_state(state);dn_dmu=0._wp
        if(.not.model%ready) then
            call fail_status(status,NS_UNINITIALIZED,'implementation','table EOS is not admitted');return
        end if
        if(side/=SIDE_REGULAR .and. side/=SIDE_INNER .and. side/=SIDE_OUTER) then
            call fail_status(status,NS_INVALID_INPUT,'physics','unknown table EOS side');return
        end if
        if(.not.ieee_is_finite(mu_j)) then
            call fail_status(status,NS_NONFINITE_INPUT,'physics','nonfinite table mu query');return
        end if
        if(mu_j<model%raw%mu_j(1) .or. mu_j>model%raw%mu_j(model%raw%n)) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','table mu query outside admitted domain');return
        end if
        i=locate_mu(model,mu_j)
        call evaluate_mu_bridge(model%cells(i),mu_j,state,dn_dmu,status)
        if(is_ok(status)) state%side=side
    end subroutine eos_table_at_mu

    subroutine eos_table_at_pressure(model,p_pa,side,state,status)
        type(table_eos_t),intent(in)::model
        real(wp),intent(in)::p_pa
        integer,intent(in)::side
        type(eos_state_t),intent(out)::state
        type(status_t),intent(out)::status
        type(eos_state_t)::trial,lo_state,hi_state
        real(wp)::lo,hi,mid,next_mu,dummy,dp,tol,best_res,repr_tol
        integer::i,k
        call clear_eos_state(state)
        if(.not.model%ready) then
            call fail_status(status,NS_UNINITIALIZED,'implementation','table EOS is not admitted');return
        end if
        if(side/=SIDE_REGULAR .and. side/=SIDE_INNER .and. side/=SIDE_OUTER) then
            call fail_status(status,NS_INVALID_INPUT,'physics','unknown table EOS side');return
        end if
        if(.not.ieee_is_finite(p_pa)) then
            call fail_status(status,NS_NONFINITE_INPUT,'physics','nonfinite table pressure query');return
        end if
        if(p_pa<model%raw%pressure_pa(1) .or. p_pa>model%raw%pressure_pa(model%raw%n)) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','table pressure query outside admitted domain');return
        end if
        i=locate_pressure(model,p_pa);lo=model%raw%mu_j(i);hi=model%raw%mu_j(i+1)
        if(p_pa==model%raw%pressure_pa(i)) then
            call evaluate_mu_bridge(model%cells(i),lo,state,dummy,status)
        else if(p_pa==model%raw%pressure_pa(i+1)) then
            call evaluate_mu_bridge(model%cells(i),hi,state,dummy,status)
        else
            tol=64._wp*epsilon(1._wp)*max(abs(p_pa),1._wp)
            mid=lo+(p_pa-model%raw%pressure_pa(i))/ &
                (model%raw%pressure_pa(i+1)-model%raw%pressure_pa(i))*(hi-lo)
            do k=1,64
                call evaluate_mu_bridge(model%cells(i),mid,trial,dummy,status)
                if(.not.is_ok(status)) return
                dp=trial%pressure_pa-p_pa
                if(abs(dp)<=tol) then
                    state=trial;exit
                end if
                if(dp<0._wp) then
                    lo=mid
                else
                    hi=mid
                end if
                ! dP/dmu = n.  Newton is used only inside the current monotone bracket;
                ! otherwise fall back to a representable bisection point.
                next_mu=mid-dp/trial%n_b_m3
                if(.not.ieee_is_finite(next_mu).or.next_mu<=lo.or.next_mu>=hi) &
                    next_mu=lo+0.5_wp*(hi-lo)
                if(next_mu==lo.or.next_mu==hi.or.next_mu==mid) exit
                mid=next_mu
            end do
            if(.not.state%valid) then
                ! The pressure root can fall between adjacent representable mu values.
                ! Select the closer endpoint rather than returning an arbitrary midpoint.
                call evaluate_mu_bridge(model%cells(i),lo,lo_state,dummy,status)
                if(.not.is_ok(status)) return
                call evaluate_mu_bridge(model%cells(i),hi,hi_state,dummy,status)
                if(.not.is_ok(status)) return
                best_res=abs(lo_state%pressure_pa-p_pa)
                state=lo_state
                if(abs(hi_state%pressure_pa-p_pa)<best_res) state=hi_state
            end if
        end if
        if(is_ok(status).and.state%valid) then
            ! At very low pressure, P is tiny compared with mu*n.  The nearest
            ! representable mu can therefore imply a relative pressure mismatch
            ! much larger than eps(P), even though the mu error is at one-ULP scale.
            ! Bound the inverse by the physically relevant local resolution dP=n dmu,
            ! then return the caller's independent pressure exactly and recompute eps.
            repr_tol=max(256._wp*epsilon(1._wp)*max(abs(p_pa),1._wp), &
                         4._wp*state%n_b_m3*spacing(state%mu_j))
            if(abs(state%pressure_pa-p_pa)>repr_tol) then
                call clear_eos_state(state)
                call fail_status(status,NS_STAGNATION,'numerics', &
                    'table pressure inverse exceeds local representable dP resolution')
                return
            end if
            state%pressure_pa=p_pa
            state%energy_j_m3=state%mu_j*state%n_b_m3-p_pa
            state%side=side
        end if
    end subroutine eos_table_at_pressure

    integer function locate_mu(model,mu_j) result(idx)
        type(table_eos_t),intent(in)::model
        real(wp),intent(in)::mu_j
        integer::lo,hi,mid
        lo=1;hi=model%raw%n-1
        do while(lo<hi)
            mid=(lo+hi)/2
            if(mu_j<=model%raw%mu_j(mid+1)) then;hi=mid;else;lo=mid+1;end if
        end do
        idx=lo
    end function locate_mu

    integer function locate_pressure(model,p_pa) result(idx)
        type(table_eos_t),intent(in)::model
        real(wp),intent(in)::p_pa
        integer::lo,hi,mid
        lo=1;hi=model%raw%n-1
        do while(lo<hi)
            mid=(lo+hi)/2
            if(p_pa<=model%raw%pressure_pa(mid+1)) then;hi=mid;else;lo=mid+1;end if
        end do
        idx=lo
    end function locate_pressure
end module ns_table_eos
