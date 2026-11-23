module compose_eos_matching
  use compose_eos_types
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  type, public :: match_report_t
    integer :: i1=0, i2=0, n_bridge=0
    ! A failed preferred construction is not evidence that the source pair
    ! has no matching window. Count the exhaustive finite endpoint search.
    integer :: window_pairs_total=0, window_pairs_examined=0
    integer :: window_candidates_constructed=0, window_candidates_rejected=0
    integer :: window_preferred_candidates=0
    logical :: window_fallback_used=.false., window_search_exhausted=.false.
    real(wp) :: n1=0._wp, n2=0._wp, p1=0._wp, p2=0._wp
    real(wp) :: f1=0._wp, f2=0._wp, mu1=0._wp, mu2=0._wp
    real(wp) :: cs2_source_left=0._wp, cs2_source_right=0._wp
    real(wp) :: cs2_left=0._wp, cs2_right=0._wp, cs2_min=0._wp, cs2_max=0._wp
    real(wp) :: first_law_relative_residual=0._wp, pressure_span=0._wp
    ! Single-speed CSS: a sharp crust->bridge jump and continuous bridge->core join.
    real(wp) :: css_single_sound_speed=0._wp
    real(wp) :: css_transition_nb_low=0._wp, css_transition_nb_high=0._wp
    real(wp) :: css_transition_delta_nb=0._wp, css_transition_delta_eps=0._wp
    real(wp) :: css_transition_epsilon_low=0._wp, css_transition_epsilon_high=0._wp
    real(wp) :: css_transition_pressure=0._wp, css_transition_mu=0._wp
    real(wp) :: seidov_critical_delta_eps=0._wp, seidov_ratio=0._wp
    integer :: css_transition_low_index=-1, css_transition_high_index=-1
    integer :: css_core_bridge_index=-1, bridge_phase_id=0
    logical :: seidov_small_core_stable=.false.
    ! CompOSE normalization preserves physical energy and chemical potentials.
    ! dm_mass_gauge_meV is retained as a legacy report field: the applied
    ! physical energy-per-baryon shift is zero. Changing a normalization
    ! reference only re-encodes Q3--Q7; it cannot add n_B*Delta m to epsilon.
    ! See data/remediation/MASS_NORMALIZATION_CONTRACT.json.
    ! source_euler_residual_meV = max over RETAINED crust rows of
    ! |mu_B,src - (F_src + P/n_B)|: provider tabulation residual repaired by
    ! the output Q3 re-derivation.
    real(wp) :: dm_mass_gauge_meV=0._wp, source_euler_residual_meV=0._wp
    ! source_euler_relative_residual = max over all crust AND core source rows of
    ! |mu_B,src - mu_E|/|mu_E| (mu_E = F + P/n_B, Il=0 charge shift undone).
    ! Gated at the provider Euler tolerance 5e-4 declared by the Maxwell kernel.
    real(wp) :: source_euler_relative_residual=0._wp
    character(len=64) :: method=''
    character(len=96) :: spline_mode=''
    character(len=64) :: construction_kind=''
    logical :: composition_dropped=.false., is_exact_maxwell=.false., seidov_applicable=.false.
  end type match_report_t
  public :: match_crust_core
contains
  subroutine pchip_slopes_local(x,y,m,ok)
    real(wp),intent(in)::x(:),y(:)
    real(wp),intent(out)::m(:)
    logical,intent(out)::ok
    real(wp),allocatable::h(:),d(:)
    real(wp)::w1,w2,edge
    integer::n,i
    n=size(x);ok=.false.;m=0._wp
    if(n<2.or.size(y)/=n.or.size(m)/=n)return
    if(.not.all(ieee_is_finite(x)).or..not.all(ieee_is_finite(y)))return
    if(any(x(2:)<=x(:n-1)))return
    allocate(h(n-1),d(n-1));h=x(2:)-x(:n-1);d=(y(2:)-y(:n-1))/h
    if(n==2)then;m(:)=d(1);ok=all(ieee_is_finite(m));return;endif
    edge=((2._wp*h(1)+h(2))*d(1)-h(1)*d(2))/(h(1)+h(2))
    if(edge*d(1)<=0._wp)then;m(1)=0._wp
    else if(d(1)*d(2)<0._wp.and.abs(edge)>3._wp*abs(d(1)))then;m(1)=3._wp*d(1)
    else;m(1)=edge;endif
    do i=2,n-1
      if(d(i-1)*d(i)<=0._wp)then;m(i)=0._wp
      else
        w1=2._wp*h(i)+h(i-1);w2=h(i)+2._wp*h(i-1)
        m(i)=(w1+w2)/(w1/d(i-1)+w2/d(i))
      endif
    enddo
    edge=((2._wp*h(n-1)+h(n-2))*d(n-1)-h(n-1)*d(n-2))/(h(n-1)+h(n-2))
    if(edge*d(n-1)<=0._wp)then;m(n)=0._wp
    else if(d(n-1)*d(n-2)<0._wp.and.abs(edge)>3._wp*abs(d(n-1)))then;m(n)=3._wp*d(n-1)
    else;m(n)=edge;endif
    ok=all(ieee_is_finite(m))
  end subroutine pchip_slopes_local

  pure subroutine pchip_value_local(x,y,m,xq,value,derivative)
    real(wp),intent(in)::x(:),y(:),m(:),xq
    real(wp),intent(out)::value,derivative
    integer::n,lo,hi,mid,i
    real(wp)::h,t,t2,t3
    n=size(x);lo=1;hi=n
    if(xq<=x(1))then;i=1
    else if(xq>=x(n))then;i=n-1
    else
      do while(hi-lo>1)
        mid=(lo+hi)/2
        if(xq<x(mid))then;hi=mid;else;lo=mid;endif
      enddo
      i=lo
    endif
    h=x(i+1)-x(i);t=(xq-x(i))/h;t2=t*t;t3=t2*t
    value=(2*t3-3*t2+1)*y(i)+(t3-2*t2+t)*h*m(i) &
      +(-2*t3+3*t2)*y(i+1)+(t3-t2)*h*m(i+1)
    derivative=((6*t2-6*t)/h)*y(i)+(3*t2-4*t+1)*m(i) &
      +((-6*t2+6*t)/h)*y(i+1)+(3*t2-2*t)*m(i+1)
  end subroutine pchip_value_local

  subroutine endpoint_cs2(eos,index,side,cs2,ok,message)
    type(eos_data_t),intent(in)::eos
    integer,intent(in)::index
    character(len=*),intent(in)::side
    real(wp),intent(out)::cs2
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(wp),allocatable::x(:),lp(:),le(:),mp(:),me(:)
    real(wp)::pval,dpdx,evalv,dedx,p,eps
    integer::lo,hi,n,j
    ok=.false.;cs2=0._wp;message='invalid endpoint sound-speed request'
    if(side=='left')then
      lo=max(1,index-5);hi=index
    else if(side=='right')then
      lo=index;hi=min(eos%n,index+5)
    else;return;endif
    ! A phase ID is categorical. Endpoint derivatives must use only the
    ! contiguous source phase on the requested side, never the next phase.
    if(side=='left')then
      do j=index-1,lo,-1
        if(eos%iphase(j)/=eos%iphase(index))then;lo=j+1;exit;endif
      enddo
    else
      do j=index+1,hi
        if(eos%iphase(j)/=eos%iphase(index))then;hi=j-1;exit;endif
      enddo
    endif
    n=hi-lo+1
    if(n<3)then;message='not enough points for one-sided endpoint sound speed';return;endif
    allocate(x(n),lp(n),le(n),mp(n),me(n))
    do j=1,n
      p=eos%nb(lo+j-1)*eos%q(1,lo+j-1)
      eps=eos%nb(lo+j-1)*(eos%q(7,lo+j-1)+1._wp)*eos%mn
      if(p<=0._wp.or.eps<=0._wp)then;message='non-positive endpoint P or epsilon';return;endif
      x(j)=log(eos%nb(lo+j-1));lp(j)=log(p);le(j)=log(eps)
    enddo
    call pchip_slopes_local(x,lp,mp,ok);if(.not.ok)then;message='cannot build pressure endpoint PCHIP';return;endif
    call pchip_slopes_local(x,le,me,ok);if(.not.ok)then;message='cannot build energy endpoint PCHIP';return;endif
    if(side=='left')then
      j=n
    else
      j=1
    endif
    call pchip_value_local(x,lp,mp,x(j),pval,dpdx)
    call pchip_value_local(x,le,me,x(j),evalv,dedx)
    if(.not.ieee_is_finite(dpdx+dedx).or.dedx<=0._wp)then
      ok=.false.;message='non-positive/non-finite logarithmic energy derivative';return
    endif
    p=exp(lp(j));eps=exp(le(j));cs2=(p/eps)*(dpdx/dedx)
    ok=ieee_is_finite(cs2).and.cs2>0._wp.and.cs2<1._wp
    if(ok)then;message='ok'
    else;message='endpoint sound speed is outside (0,1)';endif
  end subroutine endpoint_cs2

  pure subroutine cubic_segment(x0,x1,y0,y1,m0,m1,xq,value,derivative,integral)
    real(wp),intent(in)::x0,x1,y0,y1,m0,m1,xq
    real(wp),intent(out)::value,derivative,integral
    real(wp)::h,t,t2,t3,t4
    h=x1-x0;t=max(0._wp,min(1._wp,(xq-x0)/h));t2=t*t;t3=t2*t;t4=t3*t
    value=(2*t3-3*t2+1)*y0+(t3-2*t2+t)*h*m0 &
      +(-2*t3+3*t2)*y1+(t3-t2)*h*m1
    derivative=((6*t2-6*t)/h)*y0+(3*t2-4*t+1)*m0 &
      +((-6*t2+6*t)/h)*y1+(3*t2-2*t)*m1
    integral=h*((0.5_wp*t4-t3+t)*y0 &
      +(0.25_wp*t4-(2._wp/3._wp)*t3+0.5_wp*t2)*h*m0 &
      +(-0.5_wp*t4+t3)*y1+(0.25_wp*t4-(1._wp/3._wp)*t3)*h*m1)
  end subroutine cubic_segment

  subroutine spline_eval(knots,values,slopes,xq,value,derivative,integral)
    real(wp),intent(in)::knots(:),values(:),slopes(:),xq
    real(wp),intent(out)::value,derivative,integral
    real(wp)::v,d,a
    integer::n,k,j
    n=size(knots);integral=0._wp
    if(xq<=knots(1))then
      k=1
    else if(xq>=knots(n))then
      k=n-1
    else
      k=1
      do while(k<n-1)
        if(xq<=knots(k+1))exit
        k=k+1
      enddo
    endif
    do j=1,k-1
      call cubic_segment(knots(j),knots(j+1),values(j),values(j+1),slopes(j),slopes(j+1), &
                         knots(j+1),v,d,a)
      integral=integral+a
    enddo
    call cubic_segment(knots(k),knots(k+1),values(k),values(k+1),slopes(k),slopes(k+1), &
                       xq,value,derivative,a)
    integral=integral+a
  end subroutine spline_eval

  subroutine spline_integral_total(knots,values,slopes,total)
    real(wp),intent(in)::knots(:),values(:),slopes(:)
    real(wp),intent(out)::total
    real(wp)::v,d,a
    integer::k
    total=0._wp
    do k=1,size(knots)-1
      call cubic_segment(knots(k),knots(k+1),values(k),values(k+1),slopes(k),slopes(k+1), &
                         knots(k+1),v,d,a)
      total=total+a
    enddo
  end subroutine spline_integral_total

  subroutine make_c1_knots(ym,um,y0,y1,knots,values,slopes)
    real(wp),intent(in)::ym,um,y0,y1
    real(wp),intent(out)::knots(3),values(3),slopes(3)
    real(wp)::sec1,sec2,sm
    knots=[0._wp,um,1._wp];values=[y0,ym,y1]
    sec1=(ym-y0)/um;sec2=(y1-ym)/(1._wp-um)
    slopes(1)=0.85_wp*sec1;slopes(3)=0.85_wp*sec2
    if(sec1>0._wp.and.sec2>0._wp)then
      sm=2._wp*sec1*sec2/(sec1+sec2)
    else
      sm=0.5_wp*(sec1+sec2)
    endif
    slopes(2)=sm
  end subroutine make_c1_knots

  subroutine make_c2_knots(s,ua,ub,ya,yb,y0,y1,dy0,dy1,knots,values,slopes)
    real(wp),intent(in)::s,ua,ub,ya,yb,y0,y1,dy0,dy1
    real(wp),intent(out)::knots(5),values(5),slopes(5)
    real(wp)::um,ym,secs(4),sl
    integer::k
    um=ua+(0.96_wp-0.92_wp*s)*(ub-ua);ym=ya+s*(yb-ya)
    knots=[0._wp,ua,um,ub,1._wp];values=[y0,ya,ym,yb,y1]
    secs=(values(2:5)-values(1:4))/(knots(2:5)-knots(1:4))
    slopes(1)=dy0;slopes(5)=dy1
    do k=1,3
      if(secs(k)>0._wp.and.secs(k+1)>0._wp)then
        sl=2._wp*secs(k)*secs(k+1)/(secs(k)+secs(k+1))
      else
        sl=0.5_wp*(secs(k)+secs(k+1))
      endif
      slopes(k+1)=sl
    enddo
  end subroutine make_c2_knots

  pure real(wp) function polynomial_value(c,degree,x) result(value)
    real(wp),intent(in)::c(0:5),x
    integer,intent(in)::degree
    integer::j
    value=c(degree)
    do j=degree-1,0,-1
      value=value*x+c(j)
    enddo
  end function polynomial_value

  recursive subroutine polynomial_roots_unit(c,degree,roots,nroots)
    ! Isolate all real roots using derivative roots to partition [0,1]
    ! into monotone intervals. Stationary zeros retain tangencies.
    real(wp),intent(in)::c(0:5)
    integer,intent(in)::degree
    real(wp),intent(out)::roots(32)
    integer,intent(out)::nroots
    real(wp)::dc(0:5),critical(32),bounds(34),lo,hi,mid,flo,fhi,fmid,tol,normalization
    integer::d,ncritical,j,k,iter
    roots=0._wp;nroots=0;d=degree
    do while(d>0)
      if(c(d)/=0._wp)exit
      d=d-1
    enddo
    if(d==0)return
    if(d==1)then
      mid=-c(0)/c(1)
      if(mid>=0._wp.and.mid<=1._wp)then;nroots=1;roots(1)=mid;endif
      return
    endif
    dc=0._wp
    do j=1,d;dc(j-1)=real(j,wp)*c(j);enddo
    normalization=maxval(abs(dc))
    if(normalization>0._wp)dc=dc/normalization
    call polynomial_roots_unit(dc,d-1,critical,ncritical)
    bounds(1)=0._wp;bounds(2:ncritical+1)=critical(:ncritical);bounds(ncritical+2)=1._wp
    tol=32._wp*epsilon(1._wp)*sum(abs(c(:d)))
    do j=1,ncritical+2
      mid=bounds(j);fmid=polynomial_value(c,d,mid)
      if(abs(fmid)<=tol)call append_root(mid)
      if(j==ncritical+2)exit
      lo=mid;hi=bounds(j+1);flo=fmid;fhi=polynomial_value(c,d,hi)
      if(flo==0._wp.or.fhi==0._wp.or.hi<=lo)cycle
      if(sign(1._wp,flo)==sign(1._wp,fhi))cycle
      do iter=1,80
        mid=0.5_wp*(lo+hi);fmid=polynomial_value(c,d,mid)
        if(fmid==0._wp.or.hi-lo<=8._wp*epsilon(1._wp))exit
        if(sign(1._wp,flo)==sign(1._wp,fmid))then
          lo=mid;flo=fmid
        else
          hi=mid
        endif
      enddo
      call append_root(mid)
    enddo
    ! Approximate stationary zeros can be inserted before a nearby root.
    do j=2,nroots
      mid=roots(j);k=j-1
      do while(k>=1)
        if(roots(k)<=mid)exit
        roots(k+1)=roots(k);k=k-1
      enddo
      roots(k+1)=mid
    enddo
  contains
    subroutine append_root(x)
      real(wp),intent(in)::x
      if(nroots>0)then
        if(any(abs(roots(:nroots)-x)<=64._wp*epsilon(1._wp)))return
      endif
      if(nroots<32)then;nroots=nroots+1;roots(nroots)=x;endif
    end subroutine append_root
  end subroutine polynomial_roots_unit

  real(wp) function polynomial_minimum(c) result(value)
    real(wp),intent(in)::c(0:5)
    real(wp)::dc(0:5),roots(32),normalization
    integer::j,nroots
    dc=0._wp
    do j=1,5;dc(j-1)=real(j,wp)*c(j);enddo
    normalization=maxval(abs(dc))
    if(normalization>0._wp)dc=dc/normalization
    call polynomial_roots_unit(dc,4,roots,nroots)
    value=min(polynomial_value(c,5,0._wp),polynomial_value(c,5,1._wp))
    do j=1,nroots;value=min(value,polynomial_value(c,5,roots(j)));enddo
  end function polynomial_minimum

  subroutine certify_profile(coefficients,knots,canonical,axis0,span,potential0,ok,message)
    ! Certify physical inequalities on every polynomial piece, independently
    ! of the diagnostic/output nodes. In x=log n: q=P/n, F_x=q,
    ! D=dP/dn=q+q_x, 0<D<=mu=F+q. In mu: n=P_mu,
    ! n_mu>0 and mu*n_mu-n>=0; epsilon=mu*n-P>0.
    real(wp),intent(in)::coefficients(0:,:),knots(:),axis0,span,potential0
    logical,intent(in)::canonical
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(wp)::y(0:5),dy(0:5),pot(0:5),mu(0:5),stable(0:5),causal(0:5),energy(0:5)
    real(wp)::h,scale,base,axis
    integer::k,j
    ok=.false.;message='matching interval violates positivity, stability or causality'
    base=potential0
    if(.not.all(ieee_is_finite([span,axis0,potential0])).or.span<=0._wp)return
    if(.not.all(ieee_is_finite(coefficients)).or..not.all(ieee_is_finite(knots)))return
    do k=1,size(coefficients,2)
      h=knots(k+1)-knots(k);scale=span*h
      if(scale<=0._wp)return
      y=0._wp;y(:4)=coefficients(:,k);dy=0._wp;pot=0._wp;pot(0)=base
      do j=1,4;dy(j-1)=real(j,wp)*y(j);enddo
      do j=0,4;pot(j+1)=scale*y(j)/real(j+1,wp);enddo
      if(canonical)then
        mu=pot+y;stable=y+dy/scale;causal=mu-stable;energy=pot
      else
        axis=axis0+span*knots(k)
        mu=0._wp;mu(0)=axis;mu(1)=scale
        stable=dy;causal=axis*dy-scale*y
        causal(1:5)=causal(1:5)+scale*dy(:4)
        energy=axis*y-pot;energy(1:5)=energy(1:5)+scale*y(:4)
      endif
      if(.not.all(ieee_is_finite(y)).or..not.all(ieee_is_finite(pot)).or. &
         .not.all(ieee_is_finite(stable)).or..not.all(ieee_is_finite(causal)).or. &
         .not.all(ieee_is_finite(mu)).or..not.all(ieee_is_finite(energy)))return
      if(polynomial_minimum(y)<=0._wp.or.polynomial_minimum(pot)<=0._wp.or. &
         polynomial_minimum(mu)<=0._wp.or.polynomial_minimum(energy)<=0._wp.or. &
         polynomial_minimum(stable)<=0._wp.or.polynomial_minimum(causal)<0._wp)return
      base=polynomial_value(pot,5,1._wp)
    enddo
    ok=.true.;message='ok'
  end subroutine certify_profile

  subroutine certify_spline(knots,values,slopes,canonical,axis0,span,potential0,ok,message)
    real(wp),intent(in)::knots(:),values(:),slopes(:),axis0,span,potential0
    logical,intent(in)::canonical
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(wp)::coefficients(0:4,size(knots)-1),h,change
    integer::k
    coefficients=0._wp
    do k=1,size(knots)-1
      h=knots(k+1)-knots(k);change=values(k+1)-values(k)
      coefficients(0,k)=values(k);coefficients(1,k)=h*slopes(k)
      coefficients(2,k)=3._wp*change-h*(2._wp*slopes(k)+slopes(k+1))
      coefficients(3,k)=-2._wp*change+h*(slopes(k)+slopes(k+1))
    enddo
    call certify_profile(coefficients,knots,canonical,axis0,span,potential0,ok,message)
  end subroutine certify_spline

  subroutine bridge_profile(y0,y1,dy0,dy1,target_avg,match_slopes,u_grid,canonical_dx,has_canonical, &
                            y_grid,dy_grid,int_grid,d0,d1,ok,message,axis0,span,potential0)
    real(wp),intent(in)::axis0,span,potential0
    real(wp),intent(in)::y0,y1,dy0,dy1,target_avg,u_grid(:),canonical_dx
    logical,intent(in)::match_slopes,has_canonical
    real(wp),intent(out)::y_grid(:),dy_grid(:),int_grid(:),d0,d1
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(wp)::coefficients(0:4,1)
    real(wp)::theta,dy_tot,u,yy,ddy,aa,minimum,maximum,dp0,dp1
    real(wp)::c1,c2,a0,a1,a2,c3,c4,c5,yt,dyt,um,ua,ubw,ub,ya,yb
    real(wp)::knots3(3),vals3(3),slopes3(3),knots5(5),vals5(5),slopes5(5)
    real(wp)::lo,hi,mid,flo,fhi,fmid,s,secmin,secmax
    integer::i,j,iter
    logical::candidate
    ok=.false.;message='uninitialized bridge';y_grid=0._wp;dy_grid=0._wp;int_grid=0._wp;d0=0._wp;d1=0._wp
    dy_tot=y1-y0
    if(.not.all(ieee_is_finite([y0,y1,dy0,dy1,target_avg])).or. &
       (.not.has_canonical.and.dy_tot<=0._wp))then
      message='bridge endpoints must be finite and strictly increasing';return
    endif
    if(.not.has_canonical.and.abs(dy_tot)<1.e-14_wp)then
      message='degenerate bridge endpoint difference';return
    endif
    theta=0.5_wp
    if(dy_tot/=0._wp)theta=(target_avg-y0)/dy_tot
    if(.not.match_slopes)then
      if(has_canonical)then
        candidate=theta>=0.28_wp.and.theta<=0.72_wp
      else
        candidate=theta>1._wp/3._wp.and.theta<2._wp/3._wp
      endif
      minimum=huge(1._wp)
      do i=0,400
        u=real(i,wp)/400._wp
        yt=y0*(1._wp-u)*(1._wp-3._wp*u)+y1*u*(3._wp*u-2._wp)+6._wp*target_avg*u*(1._wp-u)
        dyt=y0*(-4._wp+6._wp*u)+y1*(6._wp*u-2._wp)+6._wp*target_avg*(1._wp-2._wp*u)
        if(has_canonical)then
          minimum=min(minimum,yt+dyt/canonical_dx)
        else
          minimum=min(minimum,dyt)
        endif
      enddo
      if(has_canonical)then
        candidate=candidate.and.minimum>=0.15_wp*y0
      else
        candidate=candidate.and.minimum>=0.15_wp*dy_tot
      endif
      if(candidate)then
        do j=1,size(u_grid)
          u=u_grid(j)
          y_grid(j)=y0*(1._wp-u)*(1._wp-3._wp*u)+y1*u*(3._wp*u-2._wp)+6._wp*target_avg*u*(1._wp-u)
          dy_grid(j)=y0*(-4._wp+6._wp*u)+y1*(6._wp*u-2._wp)+6._wp*target_avg*(1._wp-2._wp*u)
          int_grid(j)=y0*u*(1._wp-u)**2+y1*u*u*(u-1._wp)+target_avg*u*u*(3._wp-2._wp*u)
        enddo
        d0=-4._wp*y0-2._wp*y1+6._wp*target_avg;d1=2._wp*y0+4._wp*y1-6._wp*target_avg
        coefficients(:,1)=[y0,d0,3._wp*y0+3._wp*y1-6._wp*target_avg,0._wp,0._wp]
        call certify_profile(coefficients,[0._wp,1._wp],has_canonical,axis0,span,potential0,ok,message)
        if(.not.ok)return
        message='exact_single_cubic';return
      endif
      if(dy_tot<=0._wp)then
        message='unsupported canonical profile: exact polynomial failed interval construction';return
      endif
      um=max(0.18_wp,min(0.82_wp,1._wp-theta))
      lo=y0+1.e-5_wp*dy_tot;hi=y1-1.e-5_wp*dy_tot
      call make_c1_knots(lo,um,y0,y1,knots3,vals3,slopes3)
      call spline_integral_total(knots3,vals3,slopes3,flo);flo=flo-target_avg
      call make_c1_knots(hi,um,y0,y1,knots3,vals3,slopes3)
      call spline_integral_total(knots3,vals3,slopes3,fhi);fhi=fhi-target_avg
      if(flo*fhi>0._wp)then;message='monotone cubic bridge has no integral root';return;endif
      do iter=1,100
        mid=0.5_wp*(lo+hi)
        call make_c1_knots(mid,um,y0,y1,knots3,vals3,slopes3)
        call spline_integral_total(knots3,vals3,slopes3,fmid);fmid=fmid-target_avg
        if(abs(fmid)<2.e-14_wp*max(1._wp,abs(target_avg)))exit
        if(flo*fmid<=0._wp)then;hi=mid;fhi=fmid
        else;lo=mid;flo=fmid;endif
      enddo
      call make_c1_knots(mid,um,y0,y1,knots3,vals3,slopes3)
      d0=slopes3(1);d1=slopes3(3)
      do j=1,size(u_grid)
        call spline_eval(knots3,vals3,slopes3,u_grid(j),y_grid(j),dy_grid(j),int_grid(j))
      enddo
      call certify_spline(knots3,vals3,slopes3,has_canonical,axis0,span,potential0,ok,message)
      if(.not.ok)return
      message='monotone_1knot_cubic';return
    endif

    ! Quartic potential derivative (the interpolated EOS variable) with exact endpoint slopes.
    c1=y0;c2=0.5_wp*dy0;a0=target_avg-c1-c2;a1=y1-c1-2._wp*c2;a2=dy1-2._wp*c2
    c3=10._wp*a0-4._wp*a1+0.5_wp*a2
    c4=-15._wp*a0+7._wp*a1-a2
    c5=6._wp*a0-3._wp*a1+0.5_wp*a2
    minimum=huge(1._wp);maximum=-huge(1._wp);secmin=huge(1._wp);secmax=-huge(1._wp)
    dp0=0._wp;dp1=0._wp
    if(has_canonical)then
      dp0=y0+dy0/canonical_dx;dp1=y1+dy1/canonical_dx
    endif
    do i=0,400
      u=real(i,wp)/400._wp
      yt=c1+2*c2*u+3*c3*u*u+4*c4*u**3+5*c5*u**4
      dyt=2*c2+6*c3*u+12*c4*u*u+20*c5*u**3
      minimum=min(minimum,dyt);maximum=max(maximum,dyt)
      if(has_canonical)then
        secmin=min(secmin,yt+dyt/canonical_dx);secmax=max(secmax,yt+dyt/canonical_dx)
      endif
    enddo
    if(has_canonical)then
      candidate=secmin>=0.35_wp*min(dp0,dp1).and.secmax<=1.8_wp*max(dp0,dp1) &
        .and.(minimum>=-0.1_wp*dy_tot.or.dy_tot<=0._wp)
    else
      candidate=minimum>=0.75_wp*min(dy0,dy1,dy_tot).and.maximum<=2.2_wp*max(dy0,dy1,dy_tot)
    endif
    if(candidate)then
      do j=1,size(u_grid)
        u=u_grid(j)
        y_grid(j)=c1+2*c2*u+3*c3*u*u+4*c4*u**3+5*c5*u**4
        dy_grid(j)=2*c2+6*c3*u+12*c4*u*u+20*c5*u**3
        int_grid(j)=c1*u+c2*u*u+c3*u**3+c4*u**4+c5*u**5
      enddo
      d0=dy0;d1=dy1
      coefficients(:,1)=[c1,2._wp*c2,3._wp*c3,4._wp*c4,5._wp*c5]
      call certify_profile(coefficients,[0._wp,1._wp],has_canonical,axis0,span,potential0,ok,message)
      if(.not.ok)return
      message='exact_single_quintic';return
    endif

    if(dy_tot<=0._wp)then
      message='unsupported canonical profile: exact polynomial failed interval construction';return
    endif
    ua=max(1.e-4_wp,min(0.02_wp,0.025_wp*dy_tot/max(abs(dy0),1.e-12_wp)))
    ubw=max(1.e-4_wp,min(0.02_wp,0.025_wp*dy_tot/max(abs(dy1),1.e-12_wp)));ub=1._wp-ubw
    ya=y0+max(0.002_wp*dy_tot,min(0.025_wp*dy_tot,0.85_wp*max(dy0,0._wp)*ua))
    yb=y1-max(0.002_wp*dy_tot,min(0.025_wp*dy_tot,0.85_wp*max(dy1,0._wp)*ubw))
    lo=1.e-4_wp;hi=1._wp-1.e-4_wp
    call make_c2_knots(lo,ua,ub,ya,yb,y0,y1,dy0,dy1,knots5,vals5,slopes5)
    call spline_integral_total(knots5,vals5,slopes5,flo);flo=flo-target_avg
    call make_c2_knots(hi,ua,ub,ya,yb,y0,y1,dy0,dy1,knots5,vals5,slopes5)
    call spline_integral_total(knots5,vals5,slopes5,fhi);fhi=fhi-target_avg
    if(flo*fhi>0._wp)then;message='adaptive C2 bridge has no integral root';return;endif
    do iter=1,100
      s=0.5_wp*(lo+hi)
      call make_c2_knots(s,ua,ub,ya,yb,y0,y1,dy0,dy1,knots5,vals5,slopes5)
      call spline_integral_total(knots5,vals5,slopes5,fmid);fmid=fmid-target_avg
      if(abs(fmid)<2.e-14_wp*max(1._wp,abs(target_avg)))exit
      if(flo*fmid<=0._wp)then;hi=s;fhi=fmid
      else;lo=s;flo=fmid;endif
    enddo
    call make_c2_knots(s,ua,ub,ya,yb,y0,y1,dy0,dy1,knots5,vals5,slopes5)
    d0=dy0;d1=dy1
    do j=1,size(u_grid)
      call spline_eval(knots5,vals5,slopes5,u_grid(j),y_grid(j),dy_grid(j),int_grid(j))
    enddo
    call certify_spline(knots5,vals5,slopes5,has_canonical,axis0,span,potential0,ok,message)
    if(.not.ok)return
    message='adaptive_2knot_c2_hermite'
  end subroutine bridge_profile

  subroutine make_pressure_spline(pt,p1,p2,dpdu1,dpdu2,knots,values,slopes)
    real(wp),intent(in)::pt,p1,p2,dpdu1,dpdu2
    real(wp),intent(out)::knots(4),values(4),slopes(4)
    real(wp)::ua,ub,dplat,pa,pb,sec1,secp,sec2
    ua=0.22_wp;ub=0.78_wp;dplat=0.015_wp*(p2-p1);pa=pt-dplat;pb=pt+dplat
    sec1=(pa-p1)/ua;secp=(pb-pa)/(ub-ua);sec2=(p2-pb)/(1._wp-ub)
    knots=[0._wp,ua,ub,1._wp];values=[p1,pa,pb,p2]
    slopes(1)=min(dpdu1,2.2_wp*max(1.e-12_wp,sec1))
    slopes(2)=min(1.5_wp*secp,1.5_wp*max(1.e-12_wp,sec1))
    slopes(3)=min(1.5_wp*secp,1.5_wp*max(1.e-12_wp,sec2))
    slopes(4)=min(dpdu2,2.2_wp*max(1.e-12_wp,sec2))
  end subroutine make_pressure_spline

  subroutine gauss_legendre_rule(n,x,w)
    integer,intent(in)::n
    real(wp),intent(out)::x(n),w(n)
    integer::i,j,m,iteration
    real(wp)::z,zold,p1,p2,p3,pp,pi
    pi=acos(-1._wp);m=(n+1)/2
    do i=1,m
      z=cos(pi*(real(i,wp)-0.25_wp)/(real(n,wp)+0.5_wp))
      do iteration=1,100
        p1=1._wp;p2=0._wp
        do j=1,n
          p3=p2;p2=p1
          p1=((2._wp*real(j,wp)-1._wp)*z*p2-(real(j,wp)-1._wp)*p3)/real(j,wp)
        enddo
        pp=real(n,wp)*(z*p1-p2)/(z*z-1._wp)
        zold=z;z=zold-p1/pp
        if(abs(z-zold)<=4._wp*epsilon(1._wp)*max(1._wp,abs(z)))exit
      enddo
      ! Re-evaluate the derivative at the converged root before forming the weight.
      p1=1._wp;p2=0._wp
      do j=1,n
        p3=p2;p2=p1
        p1=((2._wp*real(j,wp)-1._wp)*z*p2-(real(j,wp)-1._wp)*p3)/real(j,wp)
      enddo
      pp=real(n,wp)*(z*p1-p2)/(z*z-1._wp)
      x(i)=-z;x(n+1-i)=z
      w(i)=2._wp/((1._wp-z*z)*pp*pp);w(n+1-i)=w(i)
    enddo
  end subroutine gauss_legendre_rule

  function pressure_work_integral(pt,uend,p1,p2,dpdu1,dpdu2,x1,dx,quadrature_x,quadrature_w) result(val)
    real(wp),intent(in)::pt,uend,p1,p2,dpdu1,dpdu2,x1,dx
    real(wp),intent(in),optional::quadrature_x(16),quadrature_w(16)
    real(wp)::val,kn(4),vv(4),ss(4),u,p,dp,aa,weight
    real(wp)::xg(16),wg(16),lower,upper
    integer::j,k
    val=0._wp
    if(uend<=0._wp)return
    ! Integrate each smooth cubic segment separately. A global quadrature
    ! crosses the C1 knots, and changing its order at u=1 breaks dF/du=dx*P/n.
    if(present(quadrature_x).and.present(quadrature_w))then
      xg=quadrature_x;wg=quadrature_w
    else
      call gauss_legendre_rule(16,xg,wg)
    endif
    call make_pressure_spline(pt,p1,p2,dpdu1,dpdu2,kn,vv,ss)
    do k=1,size(kn)-1
      lower=kn(k);upper=min(uend,kn(k+1))
      if(upper<=lower)exit
      do j=1,size(xg)
        u=lower+0.5_wp*(upper-lower)*(xg(j)+1._wp)
        weight=0.5_wp*(upper-lower)*wg(j)
        call spline_eval(kn,vv,ss,u,p,dp,aa)
        val=val+weight*p/exp(x1+u*dx)
      enddo
    enddo
  end function pressure_work_integral

  subroutine pressure_piece_margin(t,span,energy0,a,xg,wg,energy,margin,gradient)
    ! Stable evaluation of epsilon_t=span*(epsilon+P) with positive
    ! quadrature moments. Avoid subtracting a large exponential primitive.
    real(wp),intent(in)::t,span,energy0,a(0:3),xg(:),wg(:)
    real(wp),intent(out)::energy,margin,gradient
    real(wp)::v,pressure,dp,ddp,integral
    integer::j
    integral=0._wp
    do j=1,size(xg)
      v=0.5_wp*t*(xg(j)+1._wp)
      pressure=((a(3)*v+a(2))*v+a(1))*v+a(0)
      integral=integral+wg(j)*pressure*exp(-span*v)
    enddo
    integral=0.5_wp*t*integral
    energy=exp(span*t)*(energy0+span*integral)
    pressure=((a(3)*t+a(2))*t+a(1))*t+a(0)
    dp=(3._wp*a(3)*t+2._wp*a(2))*t+a(1)
    ddp=6._wp*a(3)*t+2._wp*a(2)
    margin=energy+pressure-dp/span
    gradient=span*margin+2._wp*dp-ddp/span
  end subroutine pressure_piece_margin

  subroutine certify_pressure_spline(knots,values,slopes,x0,dx,f0,ok,message)
    ! For each cubic P(t), epsilon=(epsilon0+A0) exp(s*t)-A(t),
    ! A=P+P_t/s+P_tt/s**2+P_ttt/s**3. The causal margin
    ! G=epsilon+P-P_t/s is exponential minus quadratic, so G'' has at
    ! most one zero. Partition G' by that exact zero, isolate every
    ! stationary point, and test G there and at both endpoints.
    ! Pressure and dP/dx use the polynomial interval-extrema certificate.
    real(wp),intent(in)::knots(:),values(:),slopes(:),x0,dx,f0
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(wp)::a(0:3),poly(0:5),dpoly(0:5),h,span,change,energy0,energy,margin,gradient
    real(wp)::scaled_c,ratio,critical,cuts(3),lo,hi,mid,flo,fhi,fmid,min_margin,dummy,endpoint_energy
    real(wp)::xg(32),wg(32)
    integer::k,j,ncuts,iteration
    ok=.false.;message='first-law pressure interval violates positivity, stability or causality'
    if(.not.all(ieee_is_finite([x0,dx,f0])).or.dx<=0._wp.or.f0<=0._wp)return
    if(size(knots)<2.or.size(values)/=size(knots).or.size(slopes)/=size(knots))return
    if(.not.all(ieee_is_finite(knots)).or..not.all(ieee_is_finite(values)).or. &
       .not.all(ieee_is_finite(slopes)))return
    call gauss_legendre_rule(32,xg,wg)
    energy0=exp(x0)*f0
    if(.not.ieee_is_finite(energy0).or.energy0<=0._wp)return
    do k=1,size(knots)-1
      h=knots(k+1)-knots(k);span=dx*h
      if(span<=0._wp)return
      change=values(k+1)-values(k)
      a(0)=values(k);a(1)=h*slopes(k)
      a(2)=3._wp*change-h*(2._wp*slopes(k)+slopes(k+1))
      a(3)=-2._wp*change+h*(slopes(k)+slopes(k+1))
      poly=0._wp;poly(:3)=a;dpoly=0._wp
      dpoly(:2)=[a(1),2._wp*a(2),3._wp*a(3)]
      if(polynomial_minimum(poly)<=0._wp.or.polynomial_minimum(dpoly)<=0._wp)return
      ! scaled_c=(epsilon0+A0)*s**3 avoids inverse powers of small s.
      scaled_c=(energy0+a(0))*span**3+a(1)*span**2+2._wp*a(2)*span+6._wp*a(3)
      cuts(1)=0._wp;cuts(2)=1._wp;ncuts=2
      if(scaled_c/=0._wp.and.a(3)/=0._wp)then
        ratio=12._wp*a(3)/scaled_c
        if(ieee_is_finite(ratio).and.ratio>0._wp)then
          critical=log(ratio)/span
          if(critical>0._wp.and.critical<1._wp)then
            cuts(2)=critical;cuts(3)=1._wp;ncuts=3
          endif
        endif
      endif
      call pressure_piece_margin(0._wp,span,energy0,a,xg,wg,energy,margin,gradient)
      min_margin=margin
      call pressure_piece_margin(1._wp,span,energy0,a,xg,wg,endpoint_energy,margin,gradient)
      min_margin=min(min_margin,margin)
      do j=1,ncuts-1
        lo=cuts(j);hi=cuts(j+1)
        call pressure_piece_margin(lo,span,energy0,a,xg,wg,energy,margin,flo)
        min_margin=min(min_margin,margin)
        call pressure_piece_margin(hi,span,energy0,a,xg,wg,energy,margin,fhi)
        min_margin=min(min_margin,margin)
        if(.not.all(ieee_is_finite([flo,fhi,min_margin,endpoint_energy])))return
        if(flo==0._wp.or.fhi==0._wp)cycle
        if(sign(1._wp,flo)==sign(1._wp,fhi))cycle
        do iteration=1,80
          mid=0.5_wp*(lo+hi)
          call pressure_piece_margin(mid,span,energy0,a,xg,wg,energy,margin,fmid)
          if(.not.all(ieee_is_finite([margin,fmid])))return
          if(fmid==0._wp.or.hi-lo<=8._wp*epsilon(1._wp))exit
          if(sign(1._wp,flo)==sign(1._wp,fmid))then;lo=mid;flo=fmid
          else;hi=mid;fhi=fmid;endif
        enddo
        min_margin=min(min_margin,margin)
      enddo
      if(min_margin<0._wp.or.endpoint_energy<=0._wp)return
      energy0=endpoint_energy
    enddo
    ok=.true.;message='first_law_pressure_interval_certified'
  end subroutine certify_pressure_spline

  pure function single_css_pressure_span(cs2,n_high,mu_t,mu_core) result(delta_p)
    real(wp),intent(in)::cs2,n_high,mu_t,mu_core
    real(wp)::delta_p,log_ratio,ratio
    delta_p=0._wp
    if(cs2<=0._wp.or.n_high<=0._wp.or.mu_t<=0._wp.or.mu_core<=mu_t)return
    log_ratio=log(mu_core/mu_t)
    ratio=exp(log_ratio)
    delta_p=n_high*mu_t*cs2/(1._wp+cs2)*(ratio-exp(-log_ratio/cs2))
  end function single_css_pressure_span

  subroutine solve_single_css_match(n_crust,n_core,p_transition,p_core,mu_transition,mu_core, &
                                    cs2,nb_bridge,delta_eps,ok,message)
    real(wp),intent(in)::n_crust,n_core,p_transition,p_core,mu_transition,mu_core
    real(wp),intent(out)::cs2,nb_bridge,delta_eps
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(wp)::lo,hi,mid,target,flo,fhi,fmid,log_ratio,delta_nb,epsilon_bridge
    integer::iter
    ok=.false.;cs2=0._wp;nb_bridge=0._wp;delta_eps=0._wp
    message='invalid single-speed CSS matching endpoints'
    if(.not.all(ieee_is_finite([n_crust,n_core,p_transition,p_core,mu_transition,mu_core])))return
    if(n_crust<=0._wp.or.n_core<=n_crust.or.p_transition<0._wp.or.p_core<=p_transition.or. &
       mu_transition<=0._wp.or.mu_core<=mu_transition)then
      message='single-speed CSS endpoints must increase in n_B, P and mu_B';return
    endif
    target=p_core-p_transition
    log_ratio=log(mu_core/mu_transition)
    lo=1.e-10_wp;hi=1._wp
    flo=single_css_pressure_span(lo,n_core,mu_transition,mu_core)-target
    fhi=single_css_pressure_span(hi,n_core,mu_transition,mu_core)-target
    if(.not.ieee_is_finite(flo+fhi).or.flo>0._wp.or.fhi<0._wp)then
      message='no causal single-speed CSS solution (0<cs2<=1) matches the endpoints';return
    endif
    do iter=1,140
      mid=0.5_wp*(lo+hi)
      fmid=single_css_pressure_span(mid,n_core,mu_transition,mu_core)-target
      if(abs(fmid)<=8._wp*epsilon(1._wp)*max(1._wp,abs(target)))exit
      if(fmid<0._wp)then;lo=mid
      else;hi=mid;endif
    enddo
    cs2=0.5_wp*(lo+hi)
    nb_bridge=n_core*exp(-log_ratio/cs2)
    delta_nb=nb_bridge-n_crust
    delta_eps=mu_transition*delta_nb
    epsilon_bridge=mu_transition*nb_bridge-p_transition
    if(.not.all(ieee_is_finite([cs2,nb_bridge,delta_eps,epsilon_bridge])))then
      message='single-speed CSS solution is nonfinite';return
    endif
    if(cs2<=0._wp.or.cs2>1._wp)then
      message='single-speed CSS solution violates 0<cs2<=1';return
    endif
    if(delta_nb<=64._wp*epsilon(1._wp)*max(1._wp,n_crust))then
      message='single-speed CSS solution has no positive crust-to-bridge density jump';return
    endif
    if(epsilon_bridge<=0._wp.or.delta_eps<=0._wp)then
      message='single-speed CSS transition has nonpositive energy density or jump';return
    endif
    ok=.true.;message='single_speed_css_with_crust_density_jump'
  end subroutine solve_single_css_match

  subroutine css_residual(s,n1,n2,mu1,mu2,p1,p2,res,mu_m,n_m,csa,csb)
    real(wp),intent(in)::s,n1,n2,mu1,mu2,p1,p2
    real(wp),intent(out)::res,mu_m,n_m,csa,csb
    real(wp)::ga,gb,ia,ib
    mu_m=mu1*(mu2/mu1)**(1._wp-s);n_m=n1+s*(n2-n1)
    csa=log(mu_m/mu1)/log(n_m/n1);csb=log(mu2/mu_m)/log(n2/n_m)
    if(csa<=0._wp.or.csb<=0._wp)then;res=huge(1._wp);return;endif
    ga=1._wp/csa;gb=1._wp/csb
    ia=(n1*mu1/(ga+1._wp))*((mu_m/mu1)**(ga+1._wp)-1._wp)
    ib=(n_m*mu_m/(gb+1._wp))*((mu2/mu_m)**(gb+1._wp)-1._wp)
    res=ia+ib-(p2-p1)
  end subroutine css_residual

  subroutine construct_bridge(method,u,n1,n2,p1,p2,f1,f2,mu1,mu2,q1_1,q1_2,cs1,cs2, &
                               nb,p,f,mu,q,cs,ok,message)
    character(len=*),intent(in)::method
    real(wp),intent(in)::u(:),n1,n2,p1,p2,f1,f2,mu1,mu2,q1_1,q1_2,cs1,cs2
    real(wp),intent(out)::nb(:),p(:),f(:),mu(:),q(:),cs(:)
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(wp)::dx,dmu,nbar,qbar,dpdu1,dpdu2,lo,hi,mid,fa,fb,fm,pt
    real(wp)::um,nm,csa,csb,ga,gb,pm,nbv,muv,pp,dnbdu,dqdu,d0,d1
    real(wp),allocatable::yy(:),ddy(:),anti(:)
    integer::j,iter
    real(wp)::knots_dummy(4),values_dummy(4),slopes_dummy(4),quadrature_x(16),quadrature_w(16)
    logical::b_ok
    character(len=160)::b_msg
    ok=.false.;message='uninitialized bridge';allocate(yy(size(u)),ddy(size(u)),anti(size(u)))
    nb=0._wp;p=0._wp;f=0._wp;mu=0._wp;q=0._wp;cs=0._wp
    if(size(u)<1)then;message='matching query grid is empty';return;endif
    if(.not.all(ieee_is_finite([n1,n2,p1,p2,f1,f2,mu1,mu2,q1_1,q1_2,cs1,cs2])).or. &
       .not.all(ieee_is_finite(u)))then;message='matching endpoints/query contain nonfinite values';return;endif
    if(min(n1,n2,p1,p2,f1,f2,mu1,mu2)<=0._wp.or.n2<=n1.or.mu2<=mu1.or. &
       any(u<0._wp).or.any(u>1._wp))then
      message='matching endpoints/query outside positive ordered domain';return
    endif
    dx=log(n2/n1);dmu=mu2-mu1
    nbar=(p2-p1)/dmu;qbar=(f2-f1)/dx
    select case(trim(method))
    case('cubic_F','quintic_F')
      dqdu=(mu1*cs1-q1_1)*dx;dpdu2=(mu2*cs2-q1_2)*dx
      call bridge_profile(q1_1,q1_2,dqdu,dpdu2,qbar,trim(method)=='quintic_F',u,dx,.true., &
        yy,ddy,anti,d0,d1,b_ok,b_msg,0._wp,dx,f1)
      if(.not.b_ok)then;message=trim(b_msg);return;endif
      do j=1,size(u)
        nb(j)=exp(log(n1)+u(j)*dx);q(j)=yy(j);p(j)=nb(j)*q(j)
        f(j)=f1+dx*anti(j);mu(j)=f(j)+q(j)
        cs(j)=(q(j)+ddy(j)/dx)/mu(j)
      enddo
      ok=.true.;message=trim(b_msg)
    case('cubic_P','quintic_P')
      dpdu1=(n1/(mu1*cs1))*dmu;dpdu2=(n2/(mu2*cs2))*dmu
      call bridge_profile(n1,n2,dpdu1,dpdu2,nbar,trim(method)=='quintic_P',u,0._wp,.false., &
        yy,ddy,anti,d0,d1,b_ok,b_msg,mu1,dmu,p1)
      if(.not.b_ok)then;message=trim(b_msg);return;endif
      do j=1,size(u)
        nb(j)=yy(j);mu(j)=mu1+u(j)*dmu;p(j)=p1+dmu*anti(j)
        q(j)=p(j)/nb(j);f(j)=mu(j)-q(j)
        cs(j)=(nb(j)*dmu)/(mu(j)*ddy(j))
      enddo
      ok=.true.;message=trim(b_msg)
    case('css')
      ! One constant sound speed over the entire bridge. The low-density
      ! endpoint passed here is the bridge-side state of the crust->bridge
      ! Maxwell jump; u=1 joins continuously to the high-density core.
      if(cs1<=0._wp.or.cs1>1._wp)then
        message='single-speed CSS requires 0<cs2<=1';return
      endif
      if(mu2<=mu1.or.n2<=n1.or.p2<=p1)then
        message='single-speed CSS bridge endpoints are not ordered';return
      endif
      ga=1._wp/cs1
      do j=1,size(u)
        muv=mu1*(mu2/mu1)**u(j)
        nbv=n1*(muv/mu1)**ga
        pp=p1+n1*mu1*cs1/(1._wp+cs1)*((muv/mu1)**(1._wp+ga)-1._wp)
        nb(j)=nbv;p(j)=pp;q(j)=pp/nbv;mu(j)=muv;f(j)=muv-q(j);cs(j)=cs1
      enddo
      ! The diagnostic grid includes u=1 and must hit the selected core
      ! endpoint. The production bridge grid contains interior samples only,
      ! so its last u is n_bridge/(n_bridge+1) and is not an endpoint.
      if(u(size(u))>=1._wp-8._wp*epsilon(1._wp))then
        if(abs(nb(size(u))-n2)>2.e-10_wp*max(1._wp,n2).or. &
           abs(p(size(u))-p2)>2.e-10_wp*max(1._wp,p2).or. &
           abs(mu(size(u))-mu2)>2.e-12_wp*max(1._wp,mu2))then
          message='single-speed CSS failed continuous core endpoint matching';return
        endif
      endif
      ok=.true.;message='single_speed_css_analytic'
    case('first_law_constrained_hermite_pressure_bridge')
      call gauss_legendre_rule(16,quadrature_x,quadrature_w)
      dpdu1=n1*mu1*cs1*dx;dpdu2=n2*mu2*cs2*dx
      lo=p1+1.05_wp*0.015_wp*(p2-p1);hi=p2-1.05_wp*0.015_wp*(p2-p1)
      fa=dx*pressure_work_integral(lo,1._wp,p1,p2,dpdu1,dpdu2,log(n1),dx,quadrature_x,quadrature_w)-(f2-f1)
      fb=dx*pressure_work_integral(hi,1._wp,p1,p2,dpdu1,dpdu2,log(n1),dx,quadrature_x,quadrature_w)-(f2-f1)
      if(fa*fb>0._wp)then;message='near-plateau pressure bridge has no integral root';return;endif
      do iter=1,100
        mid=0.5_wp*(lo+hi)
        fm=dx*pressure_work_integral(mid,1._wp,p1,p2,dpdu1,dpdu2,log(n1),dx,quadrature_x,quadrature_w)-(f2-f1)
        if(abs(fm)<2.e-13_wp*max(1._wp,abs(f2-f1)))exit
        if(fa*fm<=0._wp)then;hi=mid;fb=fm
        else;lo=mid;fa=fm;endif
      enddo
      pt=mid
      call make_pressure_spline(pt,p1,p2,dpdu1,dpdu2,knots_dummy,values_dummy,slopes_dummy)
      call certify_pressure_spline(knots_dummy,values_dummy,slopes_dummy,log(n1),dx,f1,b_ok,b_msg)
      if(.not.b_ok)then;message=trim(b_msg);return;endif
      do j=1,size(u)
        call spline_eval(knots_dummy,values_dummy,slopes_dummy,u(j),pp,dnbdu,anti(j))
        nbv=exp(log(n1)+u(j)*dx)
        nb(j)=nbv;p(j)=pp;f(j)=f1+dx*pressure_work_integral(pt,u(j),p1,p2,dpdu1,dpdu2,log(n1),dx,quadrature_x,quadrature_w)
        q(j)=pp/nbv;mu(j)=f(j)+q(j);cs(j)=(dnbdu/dx)/(nbv*mu(j))
      enddo
      ok=.true.;message='first_law_constrained_hermite_pressure_bridge_no_density_jump'
    case default
      message='unknown smooth matching method; valid methods: first_law_constrained_hermite_pressure_bridge, css, cubic_F, cubic_P, quintic_F, quintic_P'
    end select
  end subroutine construct_bridge

  subroutine find_window(cr,co,i1,i2,nbar,qbar,ok,message,method,report)
    type(eos_data_t),intent(in)::cr,co
    integer,intent(out)::i1,i2
    real(wp),intent(out)::nbar,qbar
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    character(len=*),intent(in)::method
    type(match_report_t),intent(inout)::report
    real(wp)::n1,n2,p1,p2,f1,f2,mu1,mu2,thp,thf,dx
    real(wp)::css_candidate,nb_candidate,de_candidate,flow,qlow
    real(wp)::lower_mean,upper_mean,roundoff,mu_mid,pressure_tol
    real(wp)::u(3),nb(3),p(3),f(3),mu(3),q(3),cs(3)
    real(wp),allocatable::pc(:),pk(:),fc(:),fk(:),ec(:),ek(:),mc(:),mk(:),cl(:),ch(:)
    logical,allocatable::left_valid(:),right_valid(:),left_cached(:),right_cached(:),left_ok(:),right_ok(:)
    integer::j,k,kcore,imin,pass,i
    logical::css_mode,preferred,candidate_ok
    character(len=256)::last_rejection,trial_message
    ok=.false.;i1=0;i2=0;nbar=0._wp;qbar=0._wp
    message='no certified matching window on supplied source endpoint grid'
    last_rejection='no endpoint pair satisfies necessary ordering/integral/source conditions'
    css_mode=trim(method)=='css'
    report%window_pairs_total=max(0,cr%n-2)*max(0,co%n-2)
    allocate(pc(cr%n),fc(cr%n),ec(cr%n),mc(cr%n),cl(cr%n),left_valid(cr%n),left_cached(cr%n),left_ok(cr%n))
    allocate(pk(co%n),fk(co%n),ek(co%n),mk(co%n),ch(co%n),right_valid(co%n),right_cached(co%n),right_ok(co%n))
    pc=cr%nb*cr%q(1,:);fc=(cr%q(6,:)+1._wp)*cr%mn;ec=cr%nb*fc;mc=fc+cr%q(1,:)
    pk=co%nb*co%q(1,:);fk=(co%q(6,:)+1._wp)*co%mn;ek=co%nb*fk;mk=fk+co%q(1,:)
    left_cached=.false.;right_cached=.false.;left_ok=.false.;right_ok=.false.;cl=0._wp;ch=0._wp
    ! Retained thermodynamic rows cannot acquire validity from a bridge.
    ! These are necessary finite-table checks, not an interpolation certificate.
    ! The causal secant bound uses the existing verify C10 tolerance 1e-6:
    ! any causal P(epsilon) connecting the source rows obeys Delta P <= Delta epsilon.
    left_valid=.false.;right_valid=.false.
    left_valid(1)=all(ieee_is_finite([pc(1),fc(1),ec(1),mc(1)])).and.min(pc(1),fc(1),ec(1),mc(1))>0._wp
    pressure_tol=1.e-12_wp*maxval(abs(pc))
    do i=2,cr%n
      roundoff=64._wp*epsilon(1._wp)*max(1._wp,abs(mc(i)),abs(mc(i-1)))
      left_valid(i)=left_valid(i-1).and.all(ieee_is_finite([pc(i),fc(i),ec(i),mc(i)])).and. &
        min(pc(i),fc(i),ec(i),mc(i))>0._wp.and.pc(i)>=pc(i-1)-pressure_tol.and.ec(i)>ec(i-1).and. &
        mc(i)>=mc(i-1)-roundoff.and.pc(i)-pc(i-1)<=(1._wp+1.e-6_wp)*(ec(i)-ec(i-1))
    enddo
    i=co%n
    right_valid(i)=all(ieee_is_finite([pk(i),fk(i),ek(i),mk(i)])).and.min(pk(i),fk(i),ek(i),mk(i))>0._wp
    pressure_tol=1.e-12_wp*maxval(abs(pk))
    do i=co%n-1,1,-1
      roundoff=64._wp*epsilon(1._wp)*max(1._wp,abs(mk(i)),abs(mk(i+1)))
      right_valid(i)=right_valid(i+1).and.all(ieee_is_finite([pk(i),fk(i),ek(i),mk(i)])).and. &
        min(pk(i),fk(i),ek(i),mk(i))>0._wp.and.pk(i+1)>=pk(i)-pressure_tol.and.ek(i+1)>ek(i).and. &
        mk(i+1)>=mk(i)-roundoff.and.pk(i+1)-pk(i)<=(1._wp+1.e-6_wp)*(ek(i+1)-ek(i))
    enddo
    ! Retain historical preferences only as ordering, never as exclusions.
    kcore=0
    if(co%composition_available)then
      if(allocated(co%abar))then
        if(size(co%abar)==co%n)then
          do k=1,co%n
            if(co%nb(k)>=0.03_wp.and.co%abar(k)==0._wp)then;kcore=k;exit;endif
          enddo
        endif
      endif
    endif
    if(kcore==0)then
      do k=1,co%n
        if(co%nb(k)>=0.045_wp)then;kcore=k;exit;endif
      enddo
    endif
    if(kcore==0)kcore=co%n+1
    imin=min(cr%n,max(12,cr%n-148))
    u=[0._wp,0.5_wp,1._wp]
    do pass=1,2
      do j=cr%n,3,-1
        n1=cr%nb(j);p1=pc(j);f1=fc(j);mu1=mc(j)
        do k=1,co%n-2
          n2=co%nb(k);p2=pk(k);f2=fk(k);mu2=mk(k)
          preferred=.false.
          if(n2>1.10_wp*n1.and.p2>1.10_wp*p1.and.mu2>mu1+0.2_wp.and.f2>f1+0.2_wp)then
            nbar=(p2-p1)/(mu2-mu1);dx=log(n2/n1);qbar=(f2-f1)/dx
            thp=(nbar-n1)/(n2-n1)
            if(co%q(1,k)==cr%q(1,j))then
              thf=huge(1._wp)
              if(abs(qbar-cr%q(1,j))<=64._wp*epsilon(1._wp)*max(1._wp,abs(f1),abs(f2))/dx)thf=0.5_wp
            else
              thf=(qbar-cr%q(1,j))/(co%q(1,k)-cr%q(1,j))
            endif
            preferred=j>=imin.and.k>=kcore.and.ieee_is_finite(thp+thf).and. &
              thp>0.12_wp.and.thp<0.88_wp.and.thf>0.12_wp.and.thf<0.88_wp
          endif
          if((pass==1.and..not.preferred).or.(pass==2.and.preferred))cycle
          report%window_pairs_examined=report%window_pairs_examined+1
          if(.not.left_valid(j).or..not.right_valid(k))cycle
          if(n2<=n1.or.p2<=p1.or.mu2<=mu1.or.f2<=f1)cycle
          nbar=(p2-p1)/(mu2-mu1);dx=log(n2/n1);qbar=(f2-f1)/dx
          ! Causality gives d(n/mu)/dmu >= 0. Integrating the endpoint
          ! bounds on n(mu) gives the necessary pressure-average bounds.
          mu_mid=0.5_wp*(mu1+mu2)
          lower_mean=n1*mu_mid/mu1;upper_mean=n2*mu_mid/mu2
          roundoff=64._wp*epsilon(1._wp)*max(1._wp,abs(nbar),abs(lower_mean),abs(upper_mean))
          if(nbar<lower_mean-roundoff.or.nbar>upper_mean+roundoff)cycle
          if(.not.left_cached(j))then
            call endpoint_cs2(cr,j,'left',cl(j),left_ok(j),trial_message)
            left_cached(j)=.true.
          endif
          if(.not.right_cached(k))then
            call endpoint_cs2(co,k,'right',ch(k),right_ok(k),trial_message)
            right_cached(k)=.true.
          endif
          if(.not.left_ok(j).or..not.right_ok(k))then
            last_rejection='endpoint sound speed could not be certified on its retained physical side'
            cycle
          endif
          report%window_candidates_constructed=report%window_candidates_constructed+1
          if(preferred)report%window_preferred_candidates=report%window_preferred_candidates+1
          if(css_mode)then
            call solve_single_css_match(n1,n2,p1,p2,mu1,mu2,css_candidate,nb_candidate,de_candidate,candidate_ok,trial_message)
            if(candidate_ok)then
              flow=mu1-p1/nb_candidate;qlow=p1/nb_candidate
              call construct_bridge(method,u,nb_candidate,n2,p1,p2,flow,f2,mu1,mu2,qlow,co%q(1,k), &
                css_candidate,css_candidate,nb,p,f,mu,q,cs,candidate_ok,trial_message)
            endif
          else
            call construct_bridge(method,u,n1,n2,p1,p2,f1,f2,mu1,mu2,cr%q(1,j),co%q(1,k),cl(j),ch(k), &
              nb,p,f,mu,q,cs,candidate_ok,trial_message)
          endif
          if(candidate_ok)call validate_bridge(nb,p,f,mu,cs,candidate_ok,trial_message)
          if(candidate_ok)then
            if(abs(f(3)-f2)>5.e-9_wp*max(1._wp,abs(f2-f1)))then
              candidate_ok=.false.;trial_message='candidate violates endpoint first-law integral tolerance'
            endif
          endif
          if(.not.candidate_ok)then
            report%window_candidates_rejected=report%window_candidates_rejected+1
            last_rejection=trim(trial_message)
            cycle
          endif
          i1=j;i2=k;ok=.true.;message='ok';report%window_fallback_used=pass==2
          return
        enddo
      enddo
    enddo
    report%window_search_exhausted=.true.
    write(message,'(a,i0,a,i0,a,a)')'no certified window after exhaustive native endpoint search; constructed=', &
      report%window_candidates_constructed,'; rejected=',report%window_candidates_rejected,'; last=',trim(last_rejection)
  end subroutine find_window

  subroutine validate_bridge(nb,p,f,mu,cs,ok,message)
    real(wp),intent(in)::nb(:),p(:),f(:),mu(:),cs(:)
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(wp),allocatable::eps(:)
    ok=.false.;message='matching bridge contains nonfinite values'
    if(.not.all(ieee_is_finite(nb)).or..not.all(ieee_is_finite(p)).or. &
       .not.all(ieee_is_finite(f)).or..not.all(ieee_is_finite(mu)).or..not.all(ieee_is_finite(cs)))return
    if(any(nb<=0._wp).or.any(p<=0._wp).or.any(f<=0._wp).or.any(mu<=0._wp))then
      message='matching bridge violates positive n_B, P, F or mu_B';return
    endif
    if(any(nb(2:)<=nb(:size(nb)-1)))then;message='matching bridge density is not strictly increasing';return;endif
    if(any(p(2:)<=p(:size(p)-1)))then;message='matching bridge pressure is not strictly increasing';return;endif
    allocate(eps(size(nb)));eps=nb*f
    if(any(eps(2:)<=eps(:size(eps)-1)))then;message='matching bridge is not strictly increasing in energy density';return;endif
    if(any(cs<=0._wp).or.any(cs>1._wp))then;message='matching bridge violates 0 < c_s^2 <= 1';return;endif
    ok=.true.;message='ok'
  end subroutine validate_bridge

  subroutine match_crust_core(cr,co,method,n_bridge,drop_composition,out,report,status)
    type(eos_data_t),intent(in)::cr,co
    character(len=*),intent(in)::method
    integer,intent(in)::n_bridge
    logical,intent(in)::drop_composition
    type(eos_data_t),intent(out)::out
    type(match_report_t),intent(out)::report
    character(len=*),intent(out)::status
    real(wp),allocatable::u(:),nb(:),p(:),f(:),mu(:),q(:),cs(:),museq(:)
    real(wp),allocatable::ud(:),nbd(:),pd(:),fd(:),mud(:),qd(:),csd(:)
    real(wp)::nbar,qbar,cs1,cs2,muq1,muq2,su,ff,qq
    real(wp)::css_cs2,css_nb_high,css_delta_eps,css_p_t,css_mu_t,css_f_high
    real(wp)::css_eps_low,css_eps_high,css_critical,css_mu_core,css_f_core
    real(wp)::css_source_cs1,css_source_cs2
    integer::i1,i2,nleft,nright,nout,j,k
    logical::ok1,ok2,ok3,css_mode
    character(len=512)::msg
    status='uninitialized';report=match_report_t()
    css_mode=trim(method)=='css'
    ! The public method="maxwell" remains reserved for the separate sharp
    ! coexistence kernel. The CSS method below is a distinct high-density
    ! constant-sound-speed phase attached by one crust-side Maxwell jump.
    if(trim(method)=='maxwell')then
      status='maxwell is the sharp coexistence method and is not yet integrated in the common composition-aware matcher; require Delta n_B>0'
      return
    endif
    call validate_eos(cr,ok1,msg);if(.not.ok1)then;status='invalid crust EOS: '//trim(msg);return;endif
    call validate_eos(co,ok1,msg);if(.not.ok1)then;status='invalid core EOS: '//trim(msg);return;endif
    if(.not.allocated(cr%iphase).or..not.allocated(co%iphase))then
      status='source phase arrays are not allocated';return
    endif
    if(size(cr%iphase)/=cr%n.or.size(co%iphase)/=co%n)then
      status='source phase arrays have inconsistent lengths';return
    endif
    if(.not.ieee_is_finite(cr%temperature_mev).or..not.ieee_is_finite(co%temperature_mev).or. &
       .not.ieee_is_finite(cr%charge_fraction_grid).or..not.ieee_is_finite(co%charge_fraction_grid))then
      status='native matching requires finite thermodynamic metadata';return
    endif
    if((cr%il/=0.and.cr%il/=1).or.(co%il/=0.and.co%il/=1))then
      status='native matching requires a supported CompOSE lepton-content flag Il=0 or Il=1';return
    endif
    ! Il determines which particles contribute to the thermodynamic table;
    ! matching cannot convert a baryonic table into a lepton-inclusive one.
    ! The physical Yq value identifies the slice; its integer index does not.
    if(cr%il/=co%il)then
      status='native matching requires the same CompOSE lepton-content convention Il';return
    endif
    if(abs(cr%charge_fraction_grid-co%charge_fraction_grid)>1.e-12_wp)then
      status='native matching requires the same physical charge-fraction slice Yq';return
    endif
    if(abs(cr%temperature_mev)>1.e-12_wp.or.abs(co%temperature_mev)>1.e-12_wp)then
      status='native matching currently requires cold T=0 barotropes';return
    endif
    if(.not.cr%beta_equilibrated.or..not.co%beta_equilibrated)then
      status='native matching requires beta-equilibrated cold barotropes';return
    endif
    ! A caller-supplied beta flag cannot override the lepton chemical potential.
    ! For Il=1 the established cold beta-equilibrium gate is |Q5| <= 1e-8.
    if((cr%il==1.and.maxval(abs(cr%q(5,:)))>1.e-8_wp).or. &
       (co%il==1.and.maxval(abs(co%q(5,:)))>1.e-8_wp))then
      status='native matching beta-equilibrium requires vanishing lepton potential Q5';return
    endif
    if(maxval(abs(cr%q(2,:)))>1.e-10_wp.or.maxval(abs(co%q(2,:)))>1.e-10_wp)then
      status='cold T=0 matching requires zero entropy per baryon Q2';return
    endif
    if(maxval(abs(cr%q(6,:)-cr%q(7,:)))>2.e-6_wp.or. &
       maxval(abs(co%q(6,:)-co%q(7,:)))>2.e-6_wp)then
      status='cold matching requires Q6 and Q7 to agree within 2e-6';return
    endif
    if((cr%composition_file_present.or.cr%composition_available.or. &
       co%composition_file_present.or.co%composition_available).and..not.drop_composition)then
      status='composition-bearing inputs require explicit drop_composition=True; native composition bridge is not yet ported';return
    endif
    if(n_bridge<8)then;status='n_bridge must be >= 8';return;endif
    if(trim(method)/='maxwell'.and.trim(method)/='first_law_constrained_hermite_pressure_bridge'.and. &
       trim(method)/='css'.and.trim(method)/='cubic_F'.and.trim(method)/='cubic_P'.and. &
       trim(method)/='quintic_F'.and.trim(method)/='quintic_P')then
      status='unknown matching method';return
    endif

    ! Output mu_B is re-derived from F + P/n_B, so an inconsistent
    ! provider chemical potential (or a misplaced pressure column) would vanish
    ! silently. Gate every source row, before the window search so that
    ! a defective table is not misreported as a physical no-window case, with the provider Euler
    ! tolerance already declared by the sharp-Maxwell kernel (5e-4).
    report%source_euler_relative_residual=max( &
      max_source_euler_relative(cr,1,cr%n),max_source_euler_relative(co,1,co%n))
    if(.not.ieee_is_finite(report%source_euler_relative_residual).or. &
       report%source_euler_relative_residual>5.e-4_wp)then
      status='source rows violate the provider Q3 Euler relation mu_B = F + P/n_B '// &
        '(relative tolerance 5e-4)';return
    endif
    call find_window(cr,co,i1,i2,nbar,qbar,ok1,msg,method,report)
    if(.not.ok1)then;status=trim(msg);return;endif
    call endpoint_cs2(cr,i1,'left',cs1,ok1,msg)
    if(.not.ok1)then;status='invalid crust endpoint sound speed: '//trim(msg);return;endif
    call endpoint_cs2(co,i2,'right',cs2,ok2,msg)
    if(.not.ok2)then;status='invalid core endpoint sound speed: '//trim(msg);return;endif
    css_source_cs1=cs1;css_source_cs2=cs2
    css_f_core=(co%q(6,i2)+1._wp)*co%mn
    css_mu_core=css_f_core+co%q(1,i2)

    css_cs2=0._wp;css_nb_high=0._wp;css_delta_eps=0._wp
    css_p_t=cr%nb(i1)*cr%q(1,i1)
    css_f_high=0._wp;css_mu_t=0._wp;css_eps_low=0._wp
    css_eps_high=0._wp;css_critical=0._wp
    if(css_mode)then
      ! Decode each source with its own reference mass. The physical crust
      ! energy and Euler chemical potential are unchanged by normalization.
      css_mu_t=(cr%q(6,i1)+1._wp)*cr%mn+cr%q(1,i1)
      css_mu_core=css_f_core+co%q(1,i2)
      call solve_single_css_match(cr%nb(i1),co%nb(i2),css_p_t,co%nb(i2)*co%q(1,i2), &
          css_mu_t,css_mu_core,css_cs2,css_nb_high,css_delta_eps,ok1,msg)
      if(.not.ok1)then;status='single-speed CSS candidate rejected: '//trim(msg);return;endif
      css_f_high=css_mu_t-css_p_t/css_nb_high
      css_eps_low=cr%nb(i1)*((cr%q(6,i1)+1._wp)*cr%mn)
      css_eps_high=css_nb_high*css_f_high
      css_critical=0.5_wp*css_eps_low+1.5_wp*css_p_t
      if(css_critical<=0._wp)then;status='single-speed CSS has nonpositive Seidov critical jump';return;endif
      cs1=css_cs2;cs2=css_cs2
      nleft=i1;nright=co%n-i2+1;nout=nleft+1+n_bridge+nright
    else
      nleft=i1;nright=co%n-i2+1;nout=nleft+n_bridge+nright
    endif

    allocate(ud(401),nbd(401),pd(401),fd(401),mud(401),qd(401),csd(401))
    do j=1,401;ud(j)=real(j-1,wp)/400._wp;enddo
    if(css_mode)then
      ! u=0 is the bridge-side transition row; u=1 is the continuous core join.
      call construct_bridge(trim(method),ud,css_nb_high,co%nb(i2),css_p_t, &
        co%nb(i2)*co%q(1,i2),css_f_high,css_f_core,css_mu_t,css_mu_core, &
        css_p_t/css_nb_high,co%q(1,i2),cs1,cs2,nbd,pd,fd,mud,qd,csd,ok2,msg)
    else
      call construct_bridge(trim(method),ud,cr%nb(i1),co%nb(i2),cr%nb(i1)*cr%q(1,i1), &
        co%nb(i2)*co%q(1,i2),(cr%q(6,i1)+1._wp)*cr%mn,css_f_core, &
        ((cr%q(6,i1)+1._wp)*cr%mn)+cr%q(1,i1),css_mu_core, &
        cr%q(1,i1),co%q(1,i2),cs1,cs2,nbd,pd,fd,mud,qd,csd,ok2,msg)
    endif
    if(.not.ok2)then;status='matching construction failed: '//trim(msg);return;endif
    call validate_bridge(nbd,pd,fd,mud,csd,ok2,msg)
    if(.not.ok2)then;status='matching candidate rejected: '//trim(msg);return;endif

    allocate(u(n_bridge),nb(n_bridge),p(n_bridge),f(n_bridge),mu(n_bridge),q(n_bridge),cs(n_bridge),museq(n_bridge))
    do j=1,n_bridge;u(j)=real(j,wp)/real(n_bridge+1,wp);enddo
    if(css_mode)then
      call construct_bridge(trim(method),u,css_nb_high,co%nb(i2),css_p_t,co%nb(i2)*co%q(1,i2), &
        css_f_high,css_f_core,css_mu_t,css_mu_core,css_p_t/css_nb_high,co%q(1,i2), &
        cs1,cs2,nb,p,f,mu,q,cs,ok3,msg)
    else
      call construct_bridge(trim(method),u,cr%nb(i1),co%nb(i2),cr%nb(i1)*cr%q(1,i1), &
        co%nb(i2)*co%q(1,i2),(cr%q(6,i1)+1._wp)*cr%mn,css_f_core, &
        ((cr%q(6,i1)+1._wp)*cr%mn)+cr%q(1,i1),css_mu_core, &
        cr%q(1,i1),co%q(1,i2),cs1,cs2,nb,p,f,mu,q,cs,ok3,msg)
    endif
    if(.not.ok3)then;status='matching construction failed on output grid: '//trim(msg);return;endif

    call allocate_eos(out,nout)
    out%mn=co%mn;out%mp=co%mp;out%il=co%il;out%iT=co%iT;out%iYq=co%iYq
    out%inb_start=1;out%nadd=0;out%temperature_mev=0._wp
    out%charge_fraction_grid=co%charge_fraction_grid
    ! This matcher exports one cold, selected thermodynamic slice. Record that
    ! as an explicit 1x1 (T,Yq) grid rather than leaving nT/nYq at their
    ! allocate_eos defaults of zero. The CompOSE writer serializes one value
    ! on each axis; the readback therefore must recover lengths 1 and 1.
    out%nT=1;out%nYq=1
    out%t_first=out%iT;out%t_last=out%iT
    out%yq_first=out%iYq;out%yq_last=out%iYq
    out%selected_iT=out%iT;out%selected_iYq=out%iYq
    allocate(out%t_axis(1),out%yq_axis(1))
    out%t_axis(1)=out%temperature_mev;out%yq_axis(1)=out%charge_fraction_grid
    ! Iphase lives in eos.compo independently of abundances. Emit phase-only
    ! rows (Npairs=Nquad=0); all presence/closure masks remain false.
    out%composition_file_present=.true.;out%composition_available=.false.;out%beta_equilibrated=.true.
    out%row_status=1
    out%q=0._wp;out%q_extra=0._wp;out%iphase=0
    ! Pick a deterministic bridge label absent from BOTH provider namespaces.
    ! Copy source labels verbatim; equality never assigns physical semantics.
    do j=0,cr%n+co%n
      if(any(cr%iphase==j).or.any(co%iphase==j))cycle
      report%bridge_phase_id=j;exit
    enddo
    out%iphase=report%bridge_phase_id
    out%yn=0._wp;out%yp=0._wp;out%ye=0._wp;out%ymu=0._wp;out%abar=0._wp;out%zbar=0._wp;out%yh=0._wp
    ! Retain crust source rows, then (for CSS) an explicit high-density
    ! transition row, the bridge interior, and finally the continuous core.
    out%nb(1:nleft)=cr%nb(1:i1)
    out%iphase(1:nleft)=cr%iphase(1:i1)
    out%q(1,1:nleft)=cr%q(1,1:i1)
    do j=1,nleft
      ff=(cr%q(6,j)+1._wp)*cr%mn;qq=cr%q(1,j)
      out%q(3,j)=(ff+qq)/out%mn-1._wp;out%q(6,j)=ff/out%mn-1._wp
      out%q(4,j)=cr%q(4,j)*cr%mn/out%mn
    enddo
    muq1=cr%q(4,i1)*cr%mn;muq2=co%q(4,i2)*co%mn
    if(css_mode)then
      ! Two adjacent rows share P_t and mu_t, but differ in n_B and epsilon.
      k=nleft+1;out%nb(k)=css_nb_high;out%q(1,k)=css_p_t/css_nb_high
      out%q(3,k)=css_mu_t/out%mn-1._wp;out%q(6,k)=css_f_high/out%mn-1._wp
      out%q(4,k)=muq1/out%mn
      do j=1,n_bridge
        k=nleft+1+j;out%nb(k)=nb(j);out%q(1,k)=q(j)
        out%q(3,k)=mu(j)/out%mn-1._wp;out%q(6,k)=f(j)/out%mn-1._wp
        su=u(j)**3*(10._wp-15._wp*u(j)+6._wp*u(j)**2)
        museq(j)=muq1+(muq2-muq1)*su;out%q(4,k)=museq(j)/out%mn
      enddo
      report%css_transition_low_index=i1-1
      report%css_transition_high_index=i1
      report%css_core_bridge_index=i1+n_bridge+1
    else
      do j=1,n_bridge
        k=nleft+j;out%nb(k)=nb(j);out%q(1,k)=q(j)
        out%q(3,k)=mu(j)/out%mn-1._wp;out%q(6,k)=f(j)/out%mn-1._wp
        su=u(j)**3*(10._wp-15._wp*u(j)+6._wp*u(j)**2)
        museq(j)=muq1+(muq2-muq1)*su;out%q(4,k)=museq(j)/out%mn
      enddo
    endif
    do j=1,nright
      if(css_mode)then
        k=nleft+n_bridge+1+j
      else
        k=nleft+n_bridge+j
      endif
      out%nb(k)=co%nb(i2+j-1);out%q(1,k)=co%q(1,i2+j-1)
      out%iphase(k)=co%iphase(i2+j-1)
      ff=(co%q(6,i2+j-1)+1._wp)*co%mn;qq=co%q(1,i2+j-1)
      out%q(3,k)=(ff+qq)/out%mn-1._wp;out%q(6,k)=ff/out%mn-1._wp
      out%q(4,k)=co%q(4,i2+j-1)*co%mn/out%mn
    enddo
    ! The construction evolves the Euler combination (epsilon+P)/n_B.  For
    ! Il=0, CompOSE Q3 excludes Yq*mu_q; preserve that source convention at
    ! the retained nonzero charge fraction rather than relabelling Yq=0.
    if(out%il==0)out%q(3,:)=out%q(3,:)-out%charge_fraction_grid*out%q(4,:)
    out%q(7,:)=out%q(6,:);out%q(2,:)=0._wp;out%q(5,:)=0._wp
    if(any(out%nb(2:)<=out%nb(:nout-1)))then;status='assembled EOS density grid is not strictly increasing';return;endif

    report%i1=i1-1;report%i2=i2-1;report%n_bridge=n_bridge
    report%n1=cr%nb(i1);report%n2=co%nb(i2)
    report%p1=css_p_t;report%p2=co%nb(i2)*co%q(1,i2)
    report%f1=(cr%q(6,i1)+1._wp)*cr%mn;report%f2=css_f_core
    report%mu1=report%f1+cr%q(1,i1);report%mu2=css_mu_core
    if(css_mode)then
      report%cs2_source_left=css_source_cs1;report%cs2_source_right=css_source_cs2
      report%is_exact_maxwell=.true.;report%seidov_applicable=.true.
      report%construction_kind='css_single_speed_with_crust_density_jump'
      report%css_single_sound_speed=css_cs2
      report%css_transition_nb_low=cr%nb(i1);report%css_transition_nb_high=css_nb_high
      report%css_transition_delta_nb=css_nb_high-cr%nb(i1)
      report%css_transition_delta_eps=css_delta_eps
      report%css_transition_epsilon_low=css_eps_low;report%css_transition_epsilon_high=css_eps_high
      report%css_transition_pressure=css_p_t;report%css_transition_mu=css_mu_t
      report%seidov_critical_delta_eps=css_critical
      report%seidov_ratio=css_delta_eps/css_critical
      report%seidov_small_core_stable=css_delta_eps<css_critical
    else
      report%cs2_source_left=cs1;report%cs2_source_right=cs2
      report%is_exact_maxwell=.false.;report%seidov_applicable=.false.
    endif
    report%cs2_left=csd(1);report%cs2_right=csd(size(csd))
    report%cs2_min=minval(csd);report%cs2_max=maxval(csd)
    report%pressure_span=maxval(pd)-minval(pd);report%spline_mode=trim(msg)
    report%first_law_relative_residual=abs(fd(size(fd))-report%f2)/max(1._wp,abs(report%f2-report%f1))
    report%method=trim(method)
    if(.not.css_mode)then
      select case(trim(method))
      case('first_law_constrained_hermite_pressure_bridge');report%construction_kind='smooth_near_pressure_plateau'
      case('cubic_F');report%construction_kind='cubic_free_energy_bridge'
      case('quintic_F');report%construction_kind='quintic_C2_free_energy_bridge'
      case('cubic_P');report%construction_kind='cubic_pressure_chemical_potential_bridge'
      case('quintic_P');report%construction_kind='quintic_C2_pressure_chemical_potential_bridge'
      end select
    endif
    ! No physical energy shift is allowed by a normalization-only operation.
    ! The provider Euler residual is disclosed separately from that operation.
    report%dm_mass_gauge_meV=0._wp
    report%source_euler_residual_meV=0._wp
    do j=1,i1
      ff=(cr%q(3,j)+1._wp)*cr%mn
      if(cr%il==0)ff=ff+cr%charge_fraction_grid*cr%q(4,j)*cr%mn
      report%source_euler_residual_meV=max(report%source_euler_residual_meV, &
        abs(ff-((cr%q(6,j)+1._wp)*cr%mn+cr%q(1,j))))
    enddo
    report%composition_dropped=cr%composition_file_present.or.cr%composition_available.or. &
      co%composition_file_present.or.co%composition_available
    if(report%first_law_relative_residual>5.e-9_wp)then
      status='matching candidate violates the endpoint first-law integral tolerance';return
    endif
    status='ok'
  end subroutine match_crust_core

  real(wp) function max_source_euler_relative(eos,first,last) result(worst)
    !! max |mu_B,src - mu_E|/|mu_E| over rows first..last, with
    !! mu_B,src=(Q3+1)m_n (+Yq*mu_q for Il=0) and mu_E=(Q6+1)m_n+Q1.
    type(eos_data_t),intent(in)::eos
    integer,intent(in)::first,last
    integer::j
    real(wp)::mub,mue
    worst=0._wp
    do j=first,last
      mub=(eos%q(3,j)+1._wp)*eos%mn
      if(eos%il==0)mub=mub+eos%charge_fraction_grid*eos%q(4,j)*eos%mn
      mue=(eos%q(6,j)+1._wp)*eos%mn+eos%q(1,j)
      if(.not.ieee_is_finite(mub).or..not.ieee_is_finite(mue).or.mue==0._wp)then
        worst=huge(1._wp);return
      endif
      worst=max(worst,abs(mub-mue)/abs(mue))
    enddo
  end function max_source_euler_relative
end module compose_eos_matching
