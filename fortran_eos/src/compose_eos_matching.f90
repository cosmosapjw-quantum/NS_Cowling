module compose_eos_matching
  use compose_eos_types
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  type, public :: match_report_t
    integer :: i1=0, i2=0, n_bridge=0
    real(wp) :: n1=0._wp, n2=0._wp, p1=0._wp, p2=0._wp
    real(wp) :: f1=0._wp, f2=0._wp, mu1=0._wp, mu2=0._wp
    real(wp) :: cs2_source_left=0._wp, cs2_source_right=0._wp
    real(wp) :: cs2_left=0._wp, cs2_right=0._wp, cs2_min=0._wp, cs2_max=0._wp
    real(wp) :: first_law_relative_residual=0._wp, pressure_span=0._wp
    character(len=32) :: method='', construction_kind='', spline_mode=''
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

  subroutine bridge_profile(y0,y1,dy0,dy1,target_avg,match_slopes,u_grid,canonical_dx,has_canonical, &
                            y_grid,dy_grid,int_grid,d0,d1,ok,message)
    real(wp),intent(in)::y0,y1,dy0,dy1,target_avg,u_grid(:),canonical_dx
    logical,intent(in)::match_slopes,has_canonical
    real(wp),intent(out)::y_grid(:),dy_grid(:),int_grid(:),d0,d1
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(wp)::theta,dy_tot,u,yy,ddy,aa,minimum,maximum,dp0,dp1
    real(wp)::c1,c2,a0,a1,a2,c3,c4,c5,yt,dyt,um,ua,ubw,ub,ya,yb
    real(wp)::knots3(3),vals3(3),slopes3(3),knots5(5),vals5(5),slopes5(5)
    real(wp)::lo,hi,mid,flo,fhi,fmid,s,secmin,secmax
    integer::i,j,iter
    logical::candidate
    ok=.false.;message='uninitialized bridge';y_grid=0._wp;dy_grid=0._wp;int_grid=0._wp;d0=0._wp;d1=0._wp
    dy_tot=y1-y0
    if(.not.all(ieee_is_finite([y0,y1,dy0,dy1,target_avg])).or.dy_tot<=0._wp)then
      message='bridge endpoints must be finite and strictly increasing';return
    endif
    if(abs(dy_tot)<1.e-14_wp)then;message='degenerate bridge endpoint difference';return;endif
    theta=(target_avg-y0)/dy_tot
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
        ok=.true.;message='exact_single_cubic';return
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
      ok=.true.;message='monotone_1knot_cubic';return
    endif

    ! Quartic potential derivative (the interpolated EOS variable) with exact endpoint slopes.
    c1=y0;c2=0.5_wp*dy0;a0=target_avg-c1-c2;a1=y1-c1-2._wp*c2;a2=dy1-2._wp*c2
    c3=10._wp*a0-4._wp*a1+0.5_wp*a2
    c4=-15._wp*a0+7._wp*a1-a2
    c5=6._wp*a0-3._wp*a1+0.5_wp*a2
    minimum=huge(1._wp);maximum=-huge(1._wp);secmin=huge(1._wp);secmax=-huge(1._wp)
    dp0=y0+dy0/canonical_dx;dp1=y1+dy1/canonical_dx
    do i=0,400
      u=real(i,wp)/400._wp
      yt=c1+2*c2*u+3*c3*u*u+4*c4*u**3+5*c5*u**4
      dyt=2*c2+6*c3*u+12*c4*u*u+20*c5*u**3
      minimum=min(minimum,dyt);maximum=max(maximum,dyt)
      secmin=min(secmin,yt+dyt/canonical_dx);secmax=max(secmax,yt+dyt/canonical_dx)
    enddo
    if(has_canonical)then
      candidate=secmin>=0.35_wp*min(dp0,dp1).and.secmax<=1.8_wp*max(dp0,dp1) &
        .and.minimum>=-0.1_wp*dy_tot
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
      d0=dy0;d1=dy1;ok=.true.;message='exact_single_quintic';return
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
    ok=.true.;message='adaptive_2knot_c2_hermite'
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

  function pressure_work_integral(pt,uend,p1,p2,dpdu1,dpdu2,x1,dx) result(val)
    real(wp),intent(in)::pt,uend,p1,p2,dpdu1,dpdu2,x1,dx
    real(wp)::val,kn(4),vv(4),ss(4),u,p,dp,aa,weight
    real(wp)::xg(64),wg(64)
    integer::j,nquad
    val=0._wp
    if(uend<=0._wp)return
    ! Match the historical Python implementation: 64-point Gauss-Legendre
    ! for the full-interval plateau root and 16-point for each partial F(u).
    nquad=16
    if(uend>=1._wp-8._wp*epsilon(1._wp))nquad=64
    call gauss_legendre_rule(nquad,xg(:nquad),wg(:nquad))
    call make_pressure_spline(pt,p1,p2,dpdu1,dpdu2,kn,vv,ss)
    do j=1,nquad
      u=0.5_wp*uend*(xg(j)+1._wp)
      weight=0.5_wp*uend*wg(j)
      call spline_eval(kn,vv,ss,u,p,dp,aa)
      val=val+weight*p/exp(x1+u*dx)
    enddo
  end function pressure_work_integral

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
    real(wp)::knots_dummy(4),values_dummy(4),slopes_dummy(4)
    logical::b_ok
    character(len=160)::b_msg
    ok=.false.;message='uninitialized bridge';allocate(yy(size(u)),ddy(size(u)),anti(size(u)))
    dx=log(n2/n1);dmu=mu2-mu1
    nbar=(p2-p1)/dmu;qbar=(f2-f1)/dx
    if(dx<=0._wp.or.dmu<=0._wp)then;message='matching endpoints are not ordered in density and chemical potential';return;endif
    select case(trim(method))
    case('cubic_F','quintic_F')
      dqdu=(mu1*cs1-q1_1)*dx;dpdu2=(mu2*cs2-q1_2)*dx
      call bridge_profile(q1_1,q1_2,dqdu,dpdu2,qbar,trim(method)=='quintic_F',u,dx,.true., &
        yy,ddy,anti,d0,d1,b_ok,b_msg)
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
        yy,ddy,anti,d0,d1,b_ok,b_msg)
      if(.not.b_ok)then;message=trim(b_msg);return;endif
      do j=1,size(u)
        nb(j)=yy(j);mu(j)=mu1+u(j)*dmu;p(j)=p1+dmu*anti(j)
        q(j)=p(j)/nb(j);f(j)=mu(j)-q(j)
        cs(j)=(nb(j)*dmu)/(mu(j)*ddy(j))
      enddo
      ok=.true.;message=trim(b_msg)
    case('css')
      lo=1.e-6_wp;hi=1._wp-1.e-6_wp
      call css_residual(lo,n1,n2,mu1,mu2,p1,p2,fa,um,nm,csa,csb)
      call css_residual(hi,n1,n2,mu1,mu2,p1,p2,fb,um,nm,csa,csb)
      if(fa*fb>0._wp)then;message='CSS bridge scalar equation has no bracketed root';return;endif
      do iter=1,120
        mid=0.5_wp*(lo+hi)
        call css_residual(mid,n1,n2,mu1,mu2,p1,p2,fm,um,nm,csa,csb)
        if(abs(fm)<2.e-13_wp*max(1._wp,abs(p2-p1)))exit
        if(fa*fm<=0._wp)then;hi=mid;fb=fm
        else;lo=mid;fa=fm;endif
      enddo
      call css_residual(mid,n1,n2,mu1,mu2,p1,p2,fm,um,nm,csa,csb)
      if(csa<=0._wp.or.csb<=0._wp.or.csa>=1._wp.or.csb>=1._wp)then;message='CSS segment sound speed outside (0,1)';return;endif
      ga=1._wp/csa;gb=1._wp/csb
      pm=p1+(n1*mu1/(ga+1._wp))*((um/mu1)**(ga+1._wp)-1._wp)
      do j=1,size(u)
        muv=mu1*(mu2/mu1)**u(j)
        if(muv<=um)then
          nbv=n1*(muv/mu1)**ga
          pp=p1+(n1*mu1/(ga+1._wp))*((muv/mu1)**(ga+1._wp)-1._wp)
          cs(j)=csa
        else
          nbv=nm*(muv/um)**gb
          pp=pm+(nm*um/(gb+1._wp))*((muv/um)**(gb+1._wp)-1._wp)
          cs(j)=csb
        endif
        nb(j)=nbv;p(j)=pp;q(j)=pp/nbv;mu(j)=muv;f(j)=muv-q(j)
      enddo
      ok=.true.;message='two_segment_css'
    case('first_law_constrained_hermite_pressure_bridge')
      dpdu1=n1*mu1*cs1*dx;dpdu2=n2*mu2*cs2*dx
      lo=p1+1.05_wp*0.015_wp*(p2-p1);hi=p2-1.05_wp*0.015_wp*(p2-p1)
      fa=dx*pressure_work_integral(lo,1._wp,p1,p2,dpdu1,dpdu2,log(n1),dx)-(f2-f1)
      fb=dx*pressure_work_integral(hi,1._wp,p1,p2,dpdu1,dpdu2,log(n1),dx)-(f2-f1)
      if(fa*fb>0._wp)then;message='near-plateau pressure bridge has no integral root';return;endif
      do iter=1,100
        mid=0.5_wp*(lo+hi)
        fm=dx*pressure_work_integral(mid,1._wp,p1,p2,dpdu1,dpdu2,log(n1),dx)-(f2-f1)
        if(abs(fm)<2.e-13_wp*max(1._wp,abs(f2-f1)))exit
        if(fa*fm<=0._wp)then;hi=mid;fb=fm
        else;lo=mid;fa=fm;endif
      enddo
      pt=mid
      call make_pressure_spline(pt,p1,p2,dpdu1,dpdu2,knots_dummy,values_dummy,slopes_dummy)
      do j=1,size(u)
        call spline_eval(knots_dummy,values_dummy,slopes_dummy,u(j),pp,dnbdu,anti(j))
        nbv=exp(log(n1)+u(j)*dx)
        nb(j)=nbv;p(j)=pp;f(j)=f1+dx*pressure_work_integral(pt,u(j),p1,p2,dpdu1,dpdu2,log(n1),dx)
        q(j)=pp/nbv;mu(j)=f(j)+q(j);cs(j)=(dnbdu/dx)/(nbv*mu(j))
      enddo
      ok=.true.;message='first_law_constrained_hermite_pressure_bridge_no_density_jump'
    case default
      message='unknown smooth matching method; valid methods: first_law_constrained_hermite_pressure_bridge, css, cubic_F, cubic_P, quintic_F, quintic_P'
    end select
  end subroutine construct_bridge

  subroutine find_window(cr,co,i1,i2,nbar,qbar,ok,message)
    type(eos_data_t),intent(in)::cr,co
    integer,intent(out)::i1,i2
    real(wp),intent(out)::nbar,qbar
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(wp)::n1,n2,p1,p2,f1,f2,mu1,mu2,thp,thf,dmcr
    integer::j,k,kcore,imin
    ok=.false.;i1=0;i2=0;nbar=0._wp;qbar=0._wp;message='no convex matching window found'
    kcore=0
    if(co%composition_available)then
      do k=1,co%n
        if(co%nb(k)>=0.03_wp.and.co%abar(k)==0._wp)then;kcore=k;exit;endif
      enddo
    endif
    if(kcore==0)then
      do k=1,co%n
        if(co%nb(k)>=0.045_wp)then;kcore=k;exit;endif
      enddo
    endif
    if(kcore==0)then;message='core table has no point at or above n_B=0.045 fm^-3';return;endif
    imin=max(12,cr%n-148);if(imin>cr%n)imin=cr%n
    dmcr=co%mn-cr%mn
    do j=cr%n,imin,-1
      n1=cr%nb(j);p1=n1*cr%q(1,j);f1=(cr%q(6,j)+1._wp)*cr%mn+dmcr;mu1=(cr%q(3,j)+1._wp)*cr%mn+dmcr
      do k=kcore,co%n
        n2=co%nb(k);p2=n2*co%q(1,k);f2=(co%q(6,k)+1._wp)*co%mn;mu2=(co%q(3,k)+1._wp)*co%mn
        if(n2<=1.10_wp*n1.or.p2<=1.10_wp*p1.or.mu2<=mu1+0.2_wp.or.f2<=f1+0.2_wp)cycle
        if(mu2<=mu1.or.n2<=n1)cycle
        nbar=(p2-p1)/(mu2-mu1);qbar=(f2-f1)/log(n2/n1)
          thp=(nbar-n1)/(n2-n1);thf=(qbar-cr%q(1,j))/(co%q(1,k)-cr%q(1,j))
        if(ieee_is_finite(thp+thf).and.thp>0.12_wp.and.thp<0.88_wp.and.thf>0.12_wp.and.thf<0.88_wp)then
          i1=j;i2=k;ok=.true.;message='ok';return
        endif
      enddo
    enddo
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
    if(any(cs<=0._wp).or.any(cs>=1._wp))then;message='matching bridge violates 0 < c_s^2 < 1';return;endif
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
    real(wp)::nbar,qbar,cs1,cs2,dmcr,muq1,muq2,su,ff,qq
    integer::i1,i2,nleft,nright,nout,j,k
    logical::ok1,ok2,ok3
    character(len=512)::msg
    status='uninitialized';report=match_report_t()
    ! The smooth pressure-bridge kernel has its own thermodynamic method ID.
    ! The public method="maxwell" is reserved for the physical sharp coexistence
    ! construction and must never enter this continuous-density bridge routine.
    if(trim(method)=='maxwell')then
      status='maxwell requires a sharp coexistence construction with Delta n_B>0; the legacy smooth plateau is disabled'
      return
    endif
    call validate_eos(cr,ok1,msg);if(.not.ok1)then;status='invalid crust EOS: '//trim(msg);return;endif
    call validate_eos(co,ok1,msg);if(.not.ok1)then;status='invalid core EOS: '//trim(msg);return;endif
    if(abs(cr%temperature_mev)>1.e-12_wp.or.abs(co%temperature_mev)>1.e-12_wp)then
      status='native smooth matching currently requires cold T=0 barotropes';return
    endif
    if(.not.cr%beta_equilibrated.or..not.co%beta_equilibrated)then
      status='native smooth matching requires beta-equilibrated cold barotropes';return
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
    if(trim(method)/='maxwell'.and.trim(method)/='first_law_constrained_hermite_pressure_bridge'.and.trim(method)/='css'.and.trim(method)/='cubic_F'.and. &
       trim(method)/='cubic_P'.and.trim(method)/='quintic_F'.and.trim(method)/='quintic_P')then
      status='unknown matching method';return
    endif
    call find_window(cr,co,i1,i2,nbar,qbar,ok1,msg)
    if(.not.ok1)then;status=trim(msg);return;endif
    call endpoint_cs2(cr,i1,'left',cs1,ok1,msg)
    if(.not.ok1)then;status='invalid crust endpoint sound speed: '//trim(msg);return;endif
    call endpoint_cs2(co,i2,'right',cs2,ok2,msg)
    if(.not.ok2)then;status='invalid core endpoint sound speed: '//trim(msg);return;endif
    dmcr=co%mn-cr%mn
    nleft=i1;nright=co%n-i2+1;nout=nleft+n_bridge+nright
    allocate(ud(401),nbd(401),pd(401),fd(401),mud(401),qd(401),csd(401))
    do j=1,401;ud(j)=real(j-1,wp)/400._wp;enddo
    call construct_bridge(trim(method),ud,cr%nb(i1),co%nb(i2),cr%nb(i1)*cr%q(1,i1), &
      co%nb(i2)*co%q(1,i2),(cr%q(6,i1)+1._wp)*cr%mn+dmcr,(co%q(6,i2)+1._wp)*co%mn, &
      ((cr%q(6,i1)+1._wp)*cr%mn+dmcr)+cr%q(1,i1),((co%q(6,i2)+1._wp)*co%mn)+co%q(1,i2), &
      cr%q(1,i1),co%q(1,i2),cs1,cs2,nbd,pd,fd,mud,qd,csd,ok2,msg)
    if(.not.ok2)then;status='matching construction failed: '//trim(msg);return;endif
    call validate_bridge(nbd,pd,fd,mud,csd,ok2,msg)
    if(.not.ok2)then;status='matching candidate rejected: '//trim(msg);return;endif
    allocate(u(n_bridge),nb(n_bridge),p(n_bridge),f(n_bridge),mu(n_bridge),q(n_bridge),cs(n_bridge),museq(n_bridge))
    do j=1,n_bridge;u(j)=real(j,wp)/real(n_bridge+1,wp);enddo
    call construct_bridge(trim(method),u,cr%nb(i1),co%nb(i2),cr%nb(i1)*cr%q(1,i1), &
      co%nb(i2)*co%q(1,i2),(cr%q(6,i1)+1._wp)*cr%mn+dmcr,(co%q(6,i2)+1._wp)*co%mn, &
      ((cr%q(6,i1)+1._wp)*cr%mn+dmcr)+cr%q(1,i1),((co%q(6,i2)+1._wp)*co%mn)+co%q(1,i2), &
      cr%q(1,i1),co%q(1,i2),cs1,cs2,nb,p,f,mu,q,cs,ok3,msg)
    if(.not.ok3)then;status='matching construction failed on output grid: '//trim(msg);return;endif
    call allocate_eos(out,nout)
    out%mn=co%mn;out%mp=co%mp;out%il=co%il;out%iT=co%iT;out%iYq=co%iYq
    out%inb_start=1;out%nadd=0;out%temperature_mev=0._wp;out%charge_fraction_grid=0._wp
    out%composition_file_present=.false.;out%composition_available=.false.;out%beta_equilibrated=.true.
    out%q=0._wp;out%q_extra=0._wp;out%iphase=0
    out%yn=0._wp;out%yp=0._wp;out%ye=0._wp;out%ymu=0._wp;out%abar=0._wp;out%zbar=0._wp;out%yh=0._wp
    ! Retain the selected crust and core source rows, applying only the documented rest-mass gauge alignment.
    out%nb(1:nleft)=cr%nb(1:i1)
    out%q(1,1:nleft)=cr%q(1,1:i1)
    do j=1,nleft
      ff=(cr%q(6,j)+1._wp)*cr%mn+dmcr;qq=cr%q(1,j)
      out%q(3,j)=(ff+qq)/out%mn-1._wp;out%q(6,j)=ff/out%mn-1._wp
      out%q(4,j)=cr%q(4,j)*cr%mn/out%mn
    enddo
    do j=1,n_bridge
      k=nleft+j;out%nb(k)=nb(j);out%q(1,k)=q(j)
      out%q(3,k)=mu(j)/out%mn-1._wp;out%q(6,k)=f(j)/out%mn-1._wp
      su=u(j)**3*(10._wp-15._wp*u(j)+6._wp*u(j)**2)
      muq1=cr%q(4,i1)*cr%mn;muq2=co%q(4,i2)*co%mn
      museq(j)=muq1+(muq2-muq1)*su;out%q(4,k)=museq(j)/out%mn
    enddo
    do j=1,nright
      k=nleft+n_bridge+j;out%nb(k)=co%nb(i2+j-1);out%q(1,k)=co%q(1,i2+j-1)
      ff=(co%q(6,i2+j-1)+1._wp)*co%mn;qq=co%q(1,i2+j-1)
      out%q(3,k)=(ff+qq)/out%mn-1._wp;out%q(6,k)=ff/out%mn-1._wp
      out%q(4,k)=co%q(4,i2+j-1)*co%mn/out%mn
    enddo
    out%q(7,:)=out%q(6,:);out%q(2,:)=0._wp;out%q(5,:)=0._wp
    if(any(out%nb(2:)<=out%nb(:nout-1)))then;status='assembled EOS density grid is not strictly increasing';return;endif
    report%i1=i1-1;report%i2=i2-1;report%n_bridge=n_bridge
    report%n1=cr%nb(i1);report%n2=co%nb(i2)
    report%p1=cr%nb(i1)*cr%q(1,i1);report%p2=co%nb(i2)*co%q(1,i2)
    report%f1=(cr%q(6,i1)+1._wp)*cr%mn+dmcr;report%f2=(co%q(6,i2)+1._wp)*co%mn
    report%mu1=report%f1+cr%q(1,i1);report%mu2=report%f2+co%q(1,i2)
    report%cs2_source_left=cs1;report%cs2_source_right=cs2
    report%cs2_left=csd(1);report%cs2_right=csd(size(csd))
    report%cs2_min=minval(csd);report%cs2_max=maxval(csd)
    report%pressure_span=maxval(pd)-minval(pd);report%spline_mode=trim(msg)
    report%first_law_relative_residual=abs(fd(size(fd))-report%f2)/max(1._wp,abs(report%f2-report%f1))
    report%method=trim(method);report%is_exact_maxwell=.false.;report%seidov_applicable=.false.
    select case(trim(method))
    case('first_law_constrained_hermite_pressure_bridge');report%construction_kind='smooth_near_pressure_plateau'
    case('css');report%construction_kind='two_segment_constant_sound_speed'
    case('cubic_F');report%construction_kind='cubic_free_energy_bridge'
    case('quintic_F');report%construction_kind='quintic_C2_free_energy_bridge'
    case('cubic_P');report%construction_kind='cubic_pressure_chemical_potential_bridge'
    case('quintic_P');report%construction_kind='quintic_C2_pressure_chemical_potential_bridge'
    end select
    report%composition_dropped=cr%composition_file_present.or.cr%composition_available.or. &
      co%composition_file_present.or.co%composition_available
    if(report%first_law_relative_residual>5.e-9_wp)then
      status='matching candidate violates the endpoint first-law integral tolerance';return
    endif
    status='ok'
  end subroutine match_crust_core
end module compose_eos_matching
