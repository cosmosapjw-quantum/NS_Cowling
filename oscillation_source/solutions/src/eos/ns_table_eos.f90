! T18: direct pressure-energy table provider.
! Original P and epsilon knots are retained exactly.  A monotone PCHIP of
! g=epsilon-P versus log(P) supplies a causal continuous barotrope; n and mu
! are then recovered from the zero-temperature first law by quadrature.
module ns_table_eos
    use ns_precision, only: wp
    use ns_status
    use ns_eos_types
    use ns_table_io, only: normalized_table_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    character(len=*), parameter, public :: TABLE_INTERPOLATION_ID = &
      'monotone_PCHIP_epsilon_minus_pressure_of_logP_with_firstlaw_density_v1'
    real(wp), parameter :: gl_x(8) = [ &
      0.095012509837637440185_wp, 0.281603550779258913230_wp, &
      0.458016777657227386342_wp, 0.617876244402643748447_wp, &
      0.755404408355003033895_wp, 0.865631202387831743880_wp, &
      0.944575023073232576078_wp, 0.989400934991649932596_wp ]
    real(wp), parameter :: gl_w(8) = [ &
      0.189450610455068496285_wp, 0.182603415044923588867_wp, &
      0.169156519395002538189_wp, 0.149595988816576732081_wp, &
      0.124628971255533872052_wp, 0.095158511682492784810_wp, &
      0.062253523938647892863_wp, 0.027152459411754094852_wp ]
    type, private :: pchip_cell_t
        real(wp) :: u0=0._wp, h=0._wp, a0=0._wp, a1=0._wp, a2=0._wp, a3=0._wp
        real(wp) :: min_dg_du=0._wp
        logical :: ready=.false.
    end type
    type, public :: table_audit_t
        integer :: row_count=0, interval_count=0, admitted_intervals=0, rejected_intervals=0
        integer :: first_rejected_interval=0
        logical :: mu_strict_increasing=.false., pressure_strict_increasing=.false.
        logical :: g_strict_increasing=.false., raw_first_law_closed=.false., fully_admitted=.false.
        real(wp) :: raw_energy_closure_max_rel=0._wp
        real(wp) :: raw_trapezoid_first_law_max_rel=0._wp
        real(wp) :: first_bad_min_d2p=0._wp, first_bad_min_causal_margin=0._wp
        real(wp) :: min_g_u=0._wp, reconstructed_n_max_rel=0._wp
    end type
    type, public :: table_eos_t
        type(normalized_table_t) :: raw
        type(pchip_cell_t), allocatable :: p_cells(:)
        real(wp), allocatable :: logn_knots(:), mu_knots(:)
        logical :: ready=.false.
    end type
    public :: audit_table_eos, init_table_eos, eos_table_at_mu, eos_table_at_pressure
contains
    subroutine audit_table_eos(raw,audit,status)
        type(normalized_table_t),intent(in)::raw
        type(table_audit_t),intent(out)::audit
        type(status_t),intent(out)::status
        type(pchip_cell_t),allocatable::cells(:)
        real(wp),allocatable::slopes(:),logn(:)
        real(wp)::pred,scale,trap,dp,g0,g1,err
        integer::i
        audit=table_audit_t()
        if(.not.raw%ready.or.raw%n<2) then
            call fail_status(status,NS_UNINITIALIZED,'implementation','normalized table is not loaded');return
        end if
        audit%row_count=raw%n; audit%interval_count=raw%n-1
        audit%mu_strict_increasing=.true.; audit%pressure_strict_increasing=.true.
        audit%g_strict_increasing=.true.
        do i=1,raw%n
            pred=raw%mu_j(i)*raw%n_b_m3(i)-raw%pressure_pa(i)
            scale=max(abs(raw%energy_j_m3(i)),tiny(1._wp))
            audit%raw_energy_closure_max_rel=max(audit%raw_energy_closure_max_rel,abs(pred-raw%energy_j_m3(i))/scale)
        end do
        do i=1,raw%n-1
            if(raw%mu_j(i+1)<=raw%mu_j(i)) audit%mu_strict_increasing=.false.
            if(raw%pressure_pa(i+1)<=raw%pressure_pa(i)) audit%pressure_strict_increasing=.false.
            g0=raw%energy_j_m3(i)-raw%pressure_pa(i)
            g1=raw%energy_j_m3(i+1)-raw%pressure_pa(i+1)
            if(g1<=g0) audit%g_strict_increasing=.false.
            trap=0.5_wp*(raw%n_b_m3(i)+raw%n_b_m3(i+1))*(raw%mu_j(i+1)-raw%mu_j(i))
            dp=raw%pressure_pa(i+1)-raw%pressure_pa(i)
            scale=max(abs(trap),abs(dp),tiny(1._wp))
            audit%raw_trapezoid_first_law_max_rel=max(audit%raw_trapezoid_first_law_max_rel,abs(trap-dp)/scale)
        end do
        audit%raw_first_law_closed=audit%raw_energy_closure_max_rel<=2.e-12_wp
        if(.not.audit%raw_first_law_closed.or..not.audit%pressure_strict_increasing.or. &
           .not.audit%g_strict_increasing.or..not.audit%mu_strict_increasing) then
            audit%rejected_intervals=raw%n-1
            audit%first_rejected_interval=1
            audit%fully_admitted=.false.
            audit%min_g_u=-huge(1._wp)
            call clear_status(status)
            return
        end if
        call construct_pchip(raw,cells,slopes,status)
        if(.not.is_ok(status)) return
        audit%admitted_intervals=0; audit%rejected_intervals=0
        audit%min_g_u=huge(1._wp)
        do i=1,size(cells)
            audit%min_g_u=min(audit%min_g_u,cells(i)%min_dg_du)
            if(cells(i)%min_dg_du>=-256._wp*epsilon(1._wp)*max(abs(raw%energy_j_m3(i)),1._wp)) then
                audit%admitted_intervals=audit%admitted_intervals+1
            else
                audit%rejected_intervals=audit%rejected_intervals+1
                if(audit%first_rejected_interval==0) audit%first_rejected_interval=i
            end if
        end do
        audit%fully_admitted=audit%rejected_intervals==0.and.audit%pressure_strict_increasing.and. &
                             audit%g_strict_increasing.and.audit%mu_strict_increasing
        if(audit%fully_admitted) then
            call integrate_all_density(raw,cells,logn,status)
            if(.not.is_ok(status)) return
            do i=1,raw%n
                err=abs(exp(logn(i))/raw%n_b_m3(i)-1._wp)
                audit%reconstructed_n_max_rel=max(audit%reconstructed_n_max_rel,err)
            end do
        end if
        call clear_status(status)
    end subroutine audit_table_eos

    subroutine init_table_eos(raw,model,audit,status)
        type(normalized_table_t),intent(in)::raw
        type(table_eos_t),intent(out)::model
        type(table_audit_t),intent(out)::audit
        type(status_t),intent(out)::status
        real(wp),allocatable::slopes(:)
        integer::i
        model=table_eos_t()
        call audit_table_eos(raw,audit,status)
        if(.not.is_ok(status)) return
        if(.not.audit%fully_admitted) then
            call fail_status(status,NS_MATCHING_INCONSISTENT,'physics', &
              'raw P-epsilon knots are not admissible for monotone causal g(log P) interpolation')
            return
        end if
        call construct_pchip(raw,model%p_cells,slopes,status)
        if(.not.is_ok(status)) return
        call integrate_all_density(raw,model%p_cells,model%logn_knots,status)
        if(.not.is_ok(status)) return
        allocate(model%mu_knots(raw%n))
        do i=1,raw%n
            model%mu_knots(i)=(raw%energy_j_m3(i)+raw%pressure_pa(i))/exp(model%logn_knots(i))
        end do
        if(any(.not.ieee_is_finite(model%mu_knots)).or.any(model%mu_knots<=0._wp)) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','reconstructed chemical potential is invalid')
            model=table_eos_t();return
        end if
        if(any(model%mu_knots(2:)<=model%mu_knots(:raw%n-1))) then
            call fail_status(status,NS_MATCHING_INCONSISTENT,'physics','first-law chemical potential is not increasing')
            model=table_eos_t();return
        end if
        model%raw=raw; model%ready=.true.; call clear_status(status)
    end subroutine init_table_eos

    subroutine construct_pchip(raw,cells,slopes,status)
        type(normalized_table_t),intent(in)::raw
        type(pchip_cell_t),allocatable,intent(out)::cells(:)
        real(wp),allocatable,intent(out)::slopes(:)
        type(status_t),intent(out)::status
        real(wp),allocatable::u(:),g(:),h(:),d(:)
        real(wp)::m0,mn,vertex,vmin,tol
        integer::n,i
        n=raw%n
        if(n<2) then
            call fail_status(status,NS_INVALID_INPUT,'data','PCHIP requires at least two rows');return
        end if
        allocate(u(n),g(n),h(n-1),d(n-1),slopes(n),cells(n-1))
        if(any(raw%pressure_pa<=0._wp)) then
            call fail_status(status,NS_INVALID_INPUT,'data','pressure must be positive for log-pressure interpolation');return
        end if
        u=log(raw%pressure_pa); g=raw%energy_j_m3-raw%pressure_pa
        h=u(2:)-u(:n-1); d=(g(2:)-g(:n-1))/h
        if(any(h<=0._wp).or.any(d<=0._wp)) then
            call fail_status(status,NS_MATCHING_INCONSISTENT,'physics','P and epsilon-P must be strictly increasing');return
        end if
        if(n==2) then
            slopes=d(1)
        else
            m0=((2._wp*h(1)+h(2))*d(1)-h(1)*d(2))/(h(1)+h(2))
            slopes(1)=min(max(m0,0._wp),3._wp*d(1))
            do i=2,n-1
                if(d(i-1)<=0._wp.or.d(i)<=0._wp) then
                    slopes(i)=0._wp
                else
                    m0=2._wp*h(i)+h(i-1)
                    mn=h(i)+2._wp*h(i-1)
                    slopes(i)=(m0+mn)/(m0/d(i-1)+mn/d(i))
                end if
            end do
            m0=((2._wp*h(n-1)+h(n-2))*d(n-1)-h(n-1)*d(n-2))/(h(n-1)+h(n-2))
            slopes(n)=min(max(m0,0._wp),3._wp*d(n-1))
        end if
        do i=1,n-1
            cells(i)%u0=u(i); cells(i)%h=h(i)
            cells(i)%a0=g(i); cells(i)%a1=h(i)*slopes(i)
            cells(i)%a2=3._wp*(g(i+1)-g(i))-h(i)*(2._wp*slopes(i)+slopes(i+1))
            cells(i)%a3=2._wp*(g(i)-g(i+1))+h(i)*(slopes(i)+slopes(i+1))
            cells(i)%min_dg_du=min_derivative(cells(i))
            tol=256._wp*epsilon(1._wp)*max(abs(g(i+1)),1._wp)
            cells(i)%ready=cells(i)%min_dg_du>=-tol
            if(.not.cells(i)%ready) then
                call fail_status(status,NS_MATCHING_INCONSISTENT,'physics','monotone PCHIP derivative audit failed');return
            end if
        end do
        call clear_status(status)
    end subroutine construct_pchip

    pure real(wp) function min_derivative(cell) result(vmin)
        type(pchip_cell_t),intent(in)::cell
        real(wp)::c0,c1,c2,t
        c0=cell%a1/cell%h
        c1=2._wp*cell%a2/cell%h
        c2=3._wp*cell%a3/cell%h
        vmin=min(c0,c0+c1+c2)
        if(c2>0._wp) then
            t=-c1/(2._wp*c2)
            if(t>0._wp.and.t<1._wp) vmin=min(vmin,c0+c1*t+c2*t*t)
        end if
    end function min_derivative

    pure subroutine eval_cell(cell,u,g,gu)
        type(pchip_cell_t),intent(in)::cell
        real(wp),intent(in)::u
        real(wp),intent(out)::g,gu
        real(wp)::t
        t=max(0._wp,min(1._wp,(u-cell%u0)/cell%h))
        g=cell%a0+cell%a1*t+cell%a2*t*t+cell%a3*t*t*t
        gu=(cell%a1+2._wp*cell%a2*t+3._wp*cell%a3*t*t)/cell%h
        if(gu<0._wp.and.gu>-512._wp*epsilon(1._wp)*max(abs(g),1._wp)) gu=0._wp
    end subroutine eval_cell

    subroutine integrate_all_density(raw,cells,logn,status)
        type(normalized_table_t),intent(in)::raw
        type(pchip_cell_t),intent(in)::cells(:)
        real(wp),allocatable,intent(out)::logn(:)
        type(status_t),intent(out)::status
        integer::i
        allocate(logn(raw%n)); logn(1)=log(raw%n_b_m3(1))
        do i=1,raw%n-1
            logn(i+1)=logn(i)+integrate_log_density(cells(i),cells(i)%u0,cells(i)%u0+cells(i)%h)
            if(.not.ieee_is_finite(logn(i+1))) then
                call fail_status(status,NS_NONFINITE_RESULT,'numerics','first-law density quadrature is nonfinite');return
            end if
        end do
        call clear_status(status)
    end subroutine integrate_all_density

    pure real(wp) function integrate_log_density(cell,ua,ub) result(value)
        type(pchip_cell_t),intent(in)::cell
        real(wp),intent(in)::ua,ub
        real(wp)::mid,half,up,um,gp,gm,dgp,dgm,pp,pm,fp,fm
        integer::j
        value=0._wp
        if(ub==ua) return
        mid=0.5_wp*(ua+ub); half=0.5_wp*(ub-ua)
        do j=1,8
            up=mid+half*gl_x(j); um=mid-half*gl_x(j)
            call integrand(up,fp)
            call integrand(um,fm)
            value=value+half*gl_w(j)*(fp+fm)
        end do
    contains
        pure subroutine integrand(u,f)
            real(wp),intent(in)::u
            real(wp),intent(out)::f
            real(wp)::g,gu,p
            call eval_cell(cell,u,g,gu)
            p=exp(u)
            f=(p+gu)/(2._wp*p+g)
        end subroutine integrand
    end function integrate_log_density

    subroutine eos_table_at_pressure(model,p_pa,side,state,status)
        type(table_eos_t),intent(in)::model
        real(wp),intent(in)::p_pa
        integer,intent(in)::side
        type(eos_state_t),intent(out)::state
        type(status_t),intent(out)::status
        real(wp)::u,g,gu,eps,lnn,n,mu,chi
        integer::i
        call clear_eos_state(state)
        if(.not.model%ready) then
            call fail_status(status,NS_UNINITIALIZED,'implementation','direct table EOS is not initialized');return
        end if
        if(side/=SIDE_REGULAR.and.side/=SIDE_INNER.and.side/=SIDE_OUTER) then
            call fail_status(status,NS_INVALID_INPUT,'physics','unknown table EOS side');return
        end if
        if(.not.ieee_is_finite(p_pa)) then
            call fail_status(status,NS_NONFINITE_INPUT,'physics','nonfinite pressure query');return
        end if
        if(p_pa<model%raw%pressure_pa(1).or.p_pa>model%raw%pressure_pa(model%raw%n)) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','pressure query outside original table domain');return
        end if
        i=locate_pressure(model,p_pa); u=log(p_pa)
        call eval_cell(model%p_cells(i),u,g,gu)
        eps=p_pa+g
        lnn=model%logn_knots(i)+integrate_log_density(model%p_cells(i),model%p_cells(i)%u0,u)
        n=exp(lnn); mu=(eps+p_pa)/n
        chi=p_pa/(p_pa+gu)
        if(.not.all(ieee_is_finite([eps,n,mu,chi])).or.min(eps,n,mu,chi)<=0._wp.or.chi>1._wp) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','direct EOS interpolation produced invalid thermodynamics');return
        end if
        state=eos_state_t(p_pa,eps,n,mu,chi,side,.true.)
        call clear_status(status)
    end subroutine eos_table_at_pressure

    subroutine eos_table_at_mu(model,mu_j,side,state,dn_dmu,status)
        type(table_eos_t),intent(in)::model
        real(wp),intent(in)::mu_j
        integer,intent(in)::side
        type(eos_state_t),intent(out)::state
        real(wp),intent(out)::dn_dmu
        type(status_t),intent(out)::status
        type(eos_state_t)::trial
        real(wp)::lo,hi,mid,dummy,tol
        integer::k
        call clear_eos_state(state); dn_dmu=0._wp
        if(.not.model%ready) then
            call fail_status(status,NS_UNINITIALIZED,'implementation','direct table EOS is not initialized');return
        end if
        if(.not.ieee_is_finite(mu_j).or.mu_j<=0._wp) then
            call fail_status(status,NS_NONFINITE_INPUT,'physics','invalid chemical-potential query');return
        end if
        if(mu_j<model%mu_knots(1).or.mu_j>model%mu_knots(model%raw%n)) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','chemical-potential query outside table domain');return
        end if
        lo=model%raw%pressure_pa(1); hi=model%raw%pressure_pa(model%raw%n)
        do k=1,80
            mid=lo+0.5_wp*(hi-lo)
            call eos_table_at_pressure(model,mid,SIDE_REGULAR,trial,status)
            if(.not.is_ok(status)) return
            if(trial%mu_j<mu_j) then
                lo=mid
            else
                hi=mid
            end if
            if(hi-lo<=8._wp*epsilon(1._wp)*max(abs(mid),tiny(1._wp))) exit
        end do
        mid=lo+0.5_wp*(hi-lo)
        call eos_table_at_pressure(model,mid,side,state,status)
        if(.not.is_ok(status)) return
        ! First-law closure gives dmu/dP=1/n and hence dn/dmu=n^2/[chi(epsilon+P)].
        dn_dmu=state%n_b_m3**2/(state%chi*(state%energy_j_m3+state%pressure_pa))
        if(.not.ieee_is_finite(dn_dmu).or.dn_dmu<=0._wp) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','dn/dmu is invalid');return
        end if
        call clear_status(status)
    end subroutine eos_table_at_mu

    integer function locate_pressure(model,p_pa) result(idx)
        type(table_eos_t),intent(in)::model
        real(wp),intent(in)::p_pa
        integer::lo,hi,mid
        lo=1;hi=model%raw%n-1
        do while(lo<hi)
            mid=(lo+hi)/2
            if(p_pa<=model%raw%pressure_pa(mid+1)) then
                hi=mid
            else
                lo=mid+1
            end if
        end do
        idx=lo
    end function locate_pressure
end module ns_table_eos
