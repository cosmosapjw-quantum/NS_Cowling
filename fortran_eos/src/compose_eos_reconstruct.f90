module compose_eos_reconstruct
  use compose_eos_types
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: fill_density_gap, reconstruct_master
contains
  subroutine pchip_slopes(x,y,m,ok,message)
    real(wp),intent(in)::x(:),y(:)
    real(wp),intent(out)::m(:)
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(wp),allocatable::h(:),d(:)
    real(wp)::w1,w2,edge
    integer::n,i
    n=size(x);ok=.false.;message='invalid PCHIP arrays'
    if(n<2.or.size(y)/=n.or.size(m)/=n)return
    if(.not.all(ieee_is_finite(x)).or..not.all(ieee_is_finite(y)))then
      message='PCHIP inputs contain nonfinite values';return
    endif
    if(any(x(2:)<=x(:n-1)))then
      message='PCHIP abscissae must be strictly increasing';return
    endif
    allocate(h(n-1),d(n-1));h=x(2:)-x(:n-1);d=(y(2:)-y(:n-1))/h
    if(n==2)then
      m=d(1);ok=.true.;message='ok';return
    endif
    edge=((2._wp*h(1)+h(2))*d(1)-h(1)*d(2))/(h(1)+h(2))
    if(edge*d(1)<=0._wp)then
      m(1)=0._wp
    else if(d(1)*d(2)<0._wp.and.abs(edge)>3._wp*abs(d(1)))then
      m(1)=3._wp*d(1)
    else
      m(1)=edge
    endif
    do i=2,n-1
      if(d(i-1)*d(i)<=0._wp)then
        m(i)=0._wp
      else
        w1=2._wp*h(i)+h(i-1);w2=h(i)+2._wp*h(i-1)
        m(i)=(w1+w2)/(w1/d(i-1)+w2/d(i))
      endif
    enddo
    edge=((2._wp*h(n-1)+h(n-2))*d(n-1)-h(n-1)*d(n-2))/(h(n-1)+h(n-2))
    if(edge*d(n-1)<=0._wp)then
      m(n)=0._wp
    else if(d(n-1)*d(n-2)<0._wp.and.abs(edge)>3._wp*abs(d(n-1)))then
      m(n)=3._wp*d(n-1)
    else
      m(n)=edge
    endif
    ok=all(ieee_is_finite(m));message='ok'
  end subroutine pchip_slopes

  pure subroutine pchip_value(x,y,m,xq,value,derivative)
    real(wp),intent(in)::x(:),y(:),m(:),xq
    real(wp),intent(out)::value,derivative
    integer::n,lo,hi,mid,i
    real(wp)::h,t,t2,t3
    n=size(x);lo=1;hi=n
    if(xq<=x(1))then
      i=1
    else if(xq>=x(n))then
      i=n-1
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
  end subroutine pchip_value

  function integrate_exp_pchip(x,y,m,a,b) result(value)
    real(wp),intent(in)::x(:),y(:),m(:),a,b
    real(wp)::value,step,s0,s1,s2,xx,lv,dv
    integer,parameter::nseg=2048
    integer::j
    if(a==b)then;value=0._wp;return;endif
    step=(b-a)/real(nseg,wp);s0=0._wp;s1=0._wp;s2=0._wp
    do j=0,nseg
      xx=a+real(j,wp)*step
      call pchip_value(x,y,m,xx,lv,dv)
      if(j==0.or.j==nseg)then
        s0=s0+exp(lv)
      else if(mod(j,2)==1)then
        s1=s1+exp(lv)
      else
        s2=s2+exp(lv)
      endif
    enddo
    value=step*(s0+4._wp*s1+2._wp*s2)/3._wp
  end function integrate_exp_pchip

  pure subroutine hermite_value(x0,x1,y0,y1,d0,d1,xq,value,derivative)
    real(wp),intent(in)::x0,x1,y0,y1,d0,d1,xq
    real(wp),intent(out)::value,derivative
    real(wp)::h,t,t2,t3
    h=x1-x0;t=(xq-x0)/h;t2=t*t;t3=t2*t
    value=(2*t3-3*t2+1)*y0+(t3-2*t2+t)*h*d0 &
         +(-2*t3+3*t2)*y1+(t3-t2)*h*d1
    derivative=((6*t2-6*t)/h)*y0+(3*t2-4*t+1)*d0 &
               +((-6*t2+6*t)/h)*y1+(3*t2-2*t)*d1
  end subroutine hermite_value

  pure subroutine hermite_second_derivative(x0,x1,y0,y1,d0,d1,xq,second_derivative)
    real(wp),intent(in)::x0,x1,y0,y1,d0,d1,xq
    real(wp),intent(out)::second_derivative
    real(wp)::h,t
    h=x1-x0;t=(xq-x0)/h
    second_derivative=((12._wp*t-6._wp)/(h*h))*y0 &
      +((6._wp*t-4._wp)/h)*d0+((-12._wp*t+6._wp)/(h*h))*y1 &
      +((6._wp*t-2._wp)/h)*d1
  end subroutine hermite_second_derivative

  subroutine interp_fit_field(xfit,yfit,xq,yq,ok,message)
    real(wp),intent(in)::xfit(:),yfit(:),xq(:)
    real(wp),intent(out)::yq(:)
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(wp),allocatable::slope(:)
    real(wp)::dy
    integer::i
    allocate(slope(size(xfit)))
    call pchip_slopes(xfit,yfit,slope,ok,message);if(.not.ok)return
    do i=1,size(xq)
      call pchip_value(xfit,yfit,slope,xq(i),yq(i),dy)
    enddo
    ok=all(ieee_is_finite(yq));if(.not.ok)message='nonfinite PCHIP interpolation result'
  end subroutine interp_fit_field

  subroutine fill_density_gap(eos,n_fill,fit_points,gap_phase,has_gap_phase,out,status)
    type(eos_data_t),intent(in)::eos
    integer,intent(in)::n_fill,fit_points,gap_phase
    logical,intent(in)::has_gap_phase
    type(eos_data_t),intent(out)::out
    character(len=*),intent(out)::status
    integer::n,igap,nfit,leftlo,righthi,k,j,ndense,lo,hi,phase_left,phase_right
    real(wp)::ratio,best,x0,x1,dx,f0,f1,q0,q1,ipc,deltaq,u,phi,intphi
    real(wp)::lv,dv,d2v,qbase,dqdx,ff,qq,mu,dpdn,cs2
    real(wp),allocatable::xfit(:),logq(:),slope(:),xfill(:),qfill(:),ffill(:),q4fit(:),q4fill(:)
    real(wp),allocatable::densex(:),densef(:),denseq(:),fitwork(:),fitout(:)
    real(wp),allocatable::ye_fill(:),yp_fill(:),ymu_fill(:),abar_fill(:),zbar_fill(:),yh_fill(:)
    logical::ok,physical,has_nuclei
    character(len=256)::why
    status='uninitialized';n=eos%n
    if(n<3)then;status='EOS requires at least three density points';return;endif
    if(n_fill<1)then;status='n_fill must be >= 1';return;endif
    if(fit_points<2)then;status='fit_points must be >= 2';return;endif
    best=0._wp;igap=1
    do j=1,n-1
      ratio=eos%nb(j+1)/eos%nb(j)
      if(ratio>best)then;best=ratio;igap=j;endif
    enddo
    if(best<5._wp)then;status='no density gap found (maximum adjacent ratio is below 5)';return;endif
    phase_left=eos%iphase(igap);phase_right=eos%iphase(igap+1)
    if(eos%composition_available.and.phase_left/=phase_right.and..not.has_gap_phase)then
      status='density gap crosses phase IDs; pass gap_phase explicitly';return
    endif
    leftlo=max(1,igap-fit_points);righthi=min(n,igap+fit_points)
    nfit=(igap-leftlo+1)+(righthi-igap)
    if(nfit<3)then;status='insufficient fit points around density gap';return;endif
    allocate(xfit(nfit),logq(nfit),slope(nfit),xfill(n_fill),qfill(n_fill),ffill(n_fill),q4fit(nfit),q4fill(n_fill))
    k=0
    do j=leftlo,igap
      k=k+1;xfit(k)=log(eos%nb(j));logq(k)=log(eos%q(1,j));q4fit(k)=eos%q(4,j)
    enddo
    do j=igap+1,righthi
      k=k+1;xfit(k)=log(eos%nb(j));logq(k)=log(eos%q(1,j));q4fit(k)=eos%q(4,j)
    enddo
    if(any(.not.ieee_is_finite(logq)).or.any(eos%q(1,leftlo:righthi)<=0._wp))then
      status='log-Q1 gap interpolation requires finite positive Q1 anchors';return
    endif
    call pchip_slopes(xfit,logq,slope,ok,why);if(.not.ok)then;status=trim(why);return;endif
    x0=log(eos%nb(igap));x1=log(eos%nb(igap+1));dx=x1-x0
    if(dx<=0._wp)then;status='density gap width is not positive';return;endif
    f0=(eos%q(6,igap)+1._wp)*eos%mn;f1=(eos%q(6,igap+1)+1._wp)*eos%mn
    q0=eos%q(1,igap);q1=eos%q(1,igap+1)
    ipc=integrate_exp_pchip(xfit,logq,slope,x0,x1)
    deltaq=(f1-f0-ipc)/dx
    do j=1,n_fill
      u=real(j,wp)/real(n_fill+1,wp);xfill(j)=x0+u*dx
      call pchip_value(xfit,logq,slope,xfill(j),lv,dv)
      qbase=exp(lv);phi=30._wp*u*u*(1._wp-u)**2
      qfill(j)=qbase+deltaq*phi
      intphi=u**3*(10._wp-15._wp*u+6._wp*u*u)
      ffill(j)=f0+integrate_exp_pchip(xfit,logq,slope,x0,xfill(j))+(f1-f0-ipc)*intphi
    enddo
    ndense=max(512,8*(n_fill+2));allocate(densex(ndense),densef(ndense),denseq(ndense))
    physical=.true.
    do j=1,ndense
      u=real(j-1,wp)/real(ndense-1,wp);densex(j)=x0+u*dx
      call pchip_value(xfit,logq,slope,densex(j),lv,dv)
      qbase=exp(lv);dqdx=qbase*dv+deltaq*(60._wp*u*(1._wp-u)*(1._wp-2._wp*u)/dx)
      qq=qbase+deltaq*(30._wp*u*u*(1._wp-u)**2)
      intphi=u**3*(10._wp-15._wp*u+6._wp*u*u)
      ff=f0+integrate_exp_pchip(xfit,logq,slope,x0,densex(j))+(f1-f0-ipc)*intphi
      mu=ff+qq;dpdn=qq+dqdx
      if(mu<=0._wp.or.qq<=0._wp.or.dpdn<=0._wp)then;physical=.false.;exit;endif
      cs2=dpdn/mu
      if(.not.ieee_is_finite(cs2).or.cs2< -1.e-8_wp.or.cs2>=1._wp+1.e-6_wp)then;physical=.false.;exit;endif
    enddo
    if(.not.physical)then
      ! The PCHIP candidate failed; evaluate the Hermite fallback from a fresh
      ! admissibility state rather than inheriting the previous failure.
      physical=.true.
      do j=1,n_fill
        call hermite_value(x0,x1,f0,f1,q0,q1,xfill(j),ffill(j),qfill(j))
      enddo
      do j=1,ndense
        ! Do not reuse densex here: the PCHIP admissibility loop may have
        ! exited early, leaving its remaining elements undefined.
        u=real(j-1,wp)/real(ndense-1,wp);lv=x0+u*dx
        call hermite_value(x0,x1,f0,f1,q0,q1,lv,ff,dv)
        call hermite_second_derivative(x0,x1,f0,f1,q0,q1,lv,d2v)
        ! Q1=F' and dP/dn_B=Q1+dQ1/dln(n_B)=F'+F''.
        mu=ff+dv;dpdn=dv+d2v
        if(.not.ieee_is_finite(ff).or..not.ieee_is_finite(dv).or. &
           .not.ieee_is_finite(d2v).or..not.ieee_is_finite(mu).or. &
           .not.ieee_is_finite(dpdn).or.dv<=0._wp.or.mu<=0._wp.or.dpdn<=0._wp)then
          physical=.false.;exit
        endif
        cs2=dpdn/mu
        if(.not.ieee_is_finite(cs2).or.cs2< -1.e-8_wp.or.cs2>=1._wp+1.e-6_wp)then;physical=.false.;exit;endif
      enddo
      if(.not.physical)then
        status='no positive, mechanically stable and causal First-Law gap interpolant satisfies endpoint constraints';return
      endif
    endif
    call interp_fit_field(xfit,q4fit,xfill,q4fill,ok,why);if(.not.ok)then;status=trim(why);return;endif
    allocate(ye_fill(n_fill),yp_fill(n_fill),ymu_fill(n_fill),abar_fill(n_fill),zbar_fill(n_fill),yh_fill(n_fill))
    allocate(fitwork(nfit),fitout(n_fill))
    fitwork=0._wp
    do j=1,nfit
      if(j<=igap-leftlo+1)then;k=leftlo+j-1;else;k=igap+1+j-(igap-leftlo+1)-1;endif
      fitwork(j)=eos%ye(k)
    enddo
    call interp_fit_field(xfit,fitwork,xfill,ye_fill,ok,why);if(.not.ok)then;status=trim(why);return;endif
    do j=1,nfit
      if(j<=igap-leftlo+1)then;k=leftlo+j-1;else;k=igap+1+j-(igap-leftlo+1)-1;endif
      fitwork(j)=eos%yp(k)
    enddo
    call interp_fit_field(xfit,fitwork,xfill,yp_fill,ok,why);if(.not.ok)then;status=trim(why);return;endif
    do j=1,nfit
      if(j<=igap-leftlo+1)then;k=leftlo+j-1;else;k=igap+1+j-(igap-leftlo+1)-1;endif
      fitwork(j)=eos%ymu(k)
    enddo
    call interp_fit_field(xfit,fitwork,xfill,ymu_fill,ok,why);if(.not.ok)then;status=trim(why);return;endif
    do j=1,nfit
      if(j<=igap-leftlo+1)then;k=leftlo+j-1;else;k=igap+1+j-(igap-leftlo+1)-1;endif
      fitwork(j)=eos%zbar(k)
    enddo
    call interp_fit_field(xfit,fitwork,xfill,zbar_fill,ok,why);if(.not.ok)then;status=trim(why);return;endif
    do j=1,nfit
      if(j<=igap-leftlo+1)then;k=leftlo+j-1;else;k=igap+1+j-(igap-leftlo+1)-1;endif
      fitwork(j)=eos%yh(k)
    enddo
    call interp_fit_field(xfit,fitwork,xfill,yh_fill,ok,why);if(.not.ok)then;status=trim(why);return;endif
    do j=1,nfit
      if(j<=igap-leftlo+1)then;k=leftlo+j-1;else;k=igap+1+j-(igap-leftlo+1)-1;endif
      fitwork(j)=eos%abar(k)
    enddo
    call interp_fit_field(xfit,fitwork,xfill,abar_fill,ok,why);if(.not.ok)then;status=trim(why);return;endif
    has_nuclei=.true.
    do j=1,n_fill
      if(ye_fill(j)+ymu_fill(j)<=1.e-8_wp.or.zbar_fill(j)<=0._wp)has_nuclei=.false.
    enddo
    if(has_nuclei)then
      abar_fill=zbar_fill/(ye_fill+ymu_fill)
    else
      zbar_fill=(ye_fill+ymu_fill)*abar_fill
    endif
    call allocate_eos(out,n+n_fill)
    out%mn=eos%mn;out%mp=eos%mp;out%il=eos%il;out%iT=eos%iT;out%iYq=eos%iYq
    out%inb_start=eos%inb_start;out%nadd=eos%nadd;out%temperature_mev=eos%temperature_mev
    out%charge_fraction_grid=eos%charge_fraction_grid;out%composition_file_present=eos%composition_file_present
    out%composition_available=eos%composition_available;out%beta_equilibrated=eos%beta_equilibrated
    lo=igap+1;hi=igap+n_fill
    out%nb(:igap)=eos%nb(:igap);out%nb(hi+1:)=eos%nb(igap+1:)
    out%q(:, :igap)=eos%q(:, :igap);out%q(:, hi+1:)=eos%q(:, igap+1:)
    out%q(2,:)=0._wp;out%q(5,:)=0._wp;out%q(7,:)=out%q(6,:)
    out%q_extra(:igap)=eos%q_extra(:igap);out%q_extra(hi+1:)=eos%q_extra(igap+1:)
    out%iphase(:igap)=eos%iphase(:igap);out%iphase(hi+1:)=eos%iphase(igap+1:)
    out%yn(:igap)=eos%yn(:igap);out%yn(hi+1:)=eos%yn(igap+1:)
    out%yp(:igap)=eos%yp(:igap);out%yp(hi+1:)=eos%yp(igap+1:)
    out%ye(:igap)=eos%ye(:igap);out%ye(hi+1:)=eos%ye(igap+1:)
    out%ymu(:igap)=eos%ymu(:igap);out%ymu(hi+1:)=eos%ymu(igap+1:)
    out%abar(:igap)=eos%abar(:igap);out%abar(hi+1:)=eos%abar(igap+1:)
    out%zbar(:igap)=eos%zbar(:igap);out%zbar(hi+1:)=eos%zbar(igap+1:)
    out%yh(:igap)=eos%yh(:igap);out%yh(hi+1:)=eos%yh(igap+1:)
    do j=1,n_fill
      k=igap+j;out%nb(k)=exp(xfill(j));out%q(1,k)=qfill(j)
      out%q(2,k)=0._wp;out%q(3,k)=(ffill(j)+qfill(j))/eos%mn-1._wp
      out%q(4,k)=q4fill(j);out%q(5,k)=0._wp;out%q(6,k)=ffill(j)/eos%mn-1._wp
      out%q(7,k)=out%q(6,k)
      if(any(eos%q_extra/=0._wp))then;out%q_extra(k)=1._wp;else;out%q_extra(k)=0._wp;endif
      if(has_gap_phase)then;out%iphase(k)=gap_phase;else;out%iphase(k)=phase_left;endif
      out%yn(k)=0._wp;out%yp(k)=yp_fill(j);out%ye(k)=ye_fill(j);out%ymu(k)=ymu_fill(j)
      out%abar(k)=abar_fill(j);out%zbar(k)=zbar_fill(j);out%yh(k)=yh_fill(j)
    enddo
    status='ok'
  end subroutine fill_density_gap


  integer function pchip_interval(x,xq) result(i)
    real(wp),intent(in)::x(:),xq
    integer::lo,hi,mid,n
    n=size(x)
    if(xq<=x(1))then
      i=1
      return
    endif
    if(xq>=x(n))then
      i=n-1
      return
    endif
    lo=1;hi=n
    do while(hi-lo>1)
      mid=(lo+hi)/2
      if(xq<x(mid))then
        hi=mid
      else
        lo=mid
      endif
    enddo
    i=lo
  end function pchip_interval

  pure subroutine pchip_value_second(x,y,m,i,xq,value,derivative,second)
    real(wp),intent(in)::x(:),y(:),m(:),xq
    integer,intent(in)::i
    real(wp),intent(out)::value,derivative,second
    real(wp)::h,t,t2,t3
    h=x(i+1)-x(i)
    t=(xq-x(i))/h;t2=t*t;t3=t2*t
    value=(2*t3-3*t2+1)*y(i)+(t3-2*t2+t)*h*m(i) &
      +(-2*t3+3*t2)*y(i+1)+(t3-t2)*h*m(i+1)
    derivative=((6*t2-6*t)/h)*y(i)+(3*t2-4*t+1)*m(i) &
      +((-6*t2+6*t)/h)*y(i+1)+(3*t2-2*t)*m(i+1)
    second=((12*t-6)/(h*h))*y(i)+((6*t-4)/h)*m(i) &
      +((-12*t+6)/(h*h))*y(i+1)+((6*t-2)/h)*m(i+1)
  end subroutine pchip_value_second

  pure real(wp) function pchip_segment_primitive(x0,x1,y0,y1,m0,m1,xq) result(value)
    real(wp),intent(in)::x0,x1,y0,y1,m0,m1,xq
    real(wp)::h,t,a00,a10,a01,a11
    h=x1-x0;t=(xq-x0)/h
    a00=0.5_wp*t**4-t**3+t
    a10=0.25_wp*t**4-(2._wp/3._wp)*t**3+0.5_wp*t**2
    a01=-0.5_wp*t**4+t**3
    a11=0.25_wp*t**4-(1._wp/3._wp)*t**3
    value=h*(a00*y0+a10*h*m0+a01*y1+a11*h*m1)
  end function pchip_segment_primitive

  subroutine pchip_primitive(x,y,m,xq,value)
    real(wp),intent(in)::x(:),y(:),m(:),xq
    real(wp),intent(out)::value
    integer::i,k,n
    real(wp)::xcut
    n=size(x);value=0._wp
    if(xq<=x(1))then
      i=1
      value=pchip_segment_primitive(x(i),x(i+1),y(i),y(i+1),m(i),m(i+1),xq)
      return
    endif
    if(xq>=x(n))then
      do k=1,n-1
        value=value+pchip_segment_primitive(x(k),x(k+1),y(k),y(k+1),m(k),m(k+1),x(k+1))
      enddo
      if(xq>x(n)) value=value+(xq-x(n))*y(n)
      return
    endif
    i=pchip_interval(x,xq)
    do k=1,i-1
      value=value+pchip_segment_primitive(x(k),x(k+1),y(k),y(k+1),m(k),m(k+1),x(k+1))
    enddo
    xcut=max(x(i),min(x(i+1),xq))
    value=value+pchip_segment_primitive(x(i),x(i+1),y(i),y(i+1),m(i),m(i+1),xcut)
  end subroutine pchip_primitive

  subroutine pchip_integral(x,y,m,a,b,value)
    real(wp),intent(in)::x(:),y(:),m(:),a,b
    real(wp),intent(out)::value
    real(wp)::pa,pb
    call pchip_primitive(x,y,m,a,pa)
    call pchip_primitive(x,y,m,b,pb)
    value=pb-pa
  end subroutine pchip_integral

  subroutine solve_normal4(matrix,rhs,solution,ok)
    real(wp),intent(in)::matrix(4,4),rhs(4)
    real(wp),intent(out)::solution(4)
    logical,intent(out)::ok
    real(wp)::a(4,5),row(5),pivot,factor,best,tmp
    integer::i,j,k,p
    a(:,1:4)=matrix
    a(:,5)=rhs
    ok=.false.;solution=0._wp
    do k=1,4
      p=k;best=abs(a(k,k))
      do i=k+1,4
        if(abs(a(i,k))>best)then
          p=i;best=abs(a(i,k))
        endif
      enddo
      if(best<=1.e-14_wp*max(1._wp,maxval(abs(a(:,1:4)))))return
      if(p/=k)then
        row=a(k,:);a(k,:)=a(p,:);a(p,:)=row
      endif
      pivot=a(k,k)
      a(k,k:5)=a(k,k:5)/pivot
      do i=1,4
        if(i==k)cycle
        factor=a(i,k)
        a(i,k:5)=a(i,k:5)-factor*a(k,k:5)
      enddo
    enddo
    solution=a(:,5)
    ok=all(ieee_is_finite(solution))
  end subroutine solve_normal4

  subroutine smooth_master_samples(x,values,sigma,smoothed,ok,message)
    real(wp),intent(in)::x(:),values(:),sigma
    real(wp),intent(out)::smoothed(:)
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(wp)::normal(4,4),rhs(4),coeff(4),z,w,powers(0:6),h,range
    integer::n,i,j,jlo,jhi,radius,a,b
    logical::solved
    n=size(x);ok=.false.;message='invalid local master-function smoothing arrays'
    if(n<5.or.size(values)/=n.or.size(smoothed)/=n)return
    if(sigma<=0._wp.or.sigma>0.5_wp.or..not.ieee_is_finite(sigma))then
      message='smooth_sigma must lie in (0, 0.5]'
      return
    endif
    if(any(.not.ieee_is_finite(values)))then
      message='master-function smoothing received nonfinite values'
      return
    endif
    range=x(n)-x(1)
    if(range<=0._wp)then
      message='master-function smoothing requires an increasing log-density grid'
      return
    endif
    radius=max(2,min((n-1)/2,nint(sigma*real(n,wp))))
    do i=1,n
      jlo=max(1,i-radius);jhi=min(n,i+radius)
      if(jhi-jlo+1<4)then
        jlo=max(1,min(i-2,n-3))
        jhi=jlo+3
      endif
      h=max(x(i)-x(jlo),x(jhi)-x(i))
      if(h<=tiny(1._wp))then
        message='degenerate local smoothing bandwidth'
        return
      endif
      normal=0._wp;rhs=0._wp
      do j=jlo,jhi
        z=(x(j)-x(i))/h
        w=exp(-2._wp*z*z)
        powers(0)=1._wp
        do a=1,6
          powers(a)=powers(a-1)*z
        enddo
        do a=0,3
          rhs(a+1)=rhs(a+1)+w*powers(a)*values(j)
          do b=0,3
            normal(a+1,b+1)=normal(a+1,b+1)+w*powers(a+b)
          enddo
        enddo
      enddo
      call solve_normal4(normal,rhs,coeff,solved)
      if(.not.solved)then
        message='local cubic least-squares fit is singular'
        return
      endif
      smoothed(i)=coeff(1)
    enddo
    ok=all(ieee_is_finite(smoothed))
    if(ok)then
      message='ok'
    else
      message='smoothed master function contains nonfinite values'
    endif
  end subroutine smooth_master_samples

  pure real(wp) function lepton_chemical_potential(y,nb,mass) result(mu)
    real(wp),intent(in)::y,nb,mass
    real(wp),parameter::hbarc=197.3269804_wp,pi=acos(-1._wp)
    real(wp)::kf
    mu=mass
    if(y<=0._wp.or.nb<=0._wp)return
    kf=(3._wp*pi*pi*nb*y)**(1._wp/3._wp)
    mu=sqrt(mass*mass+(hbarc*kf)**2)
  end function lepton_chemical_potential

  pure real(wp) function electron_pressure(nb,ye) result(pressure)
    real(wp),intent(in)::nb,ye
    real(wp),parameter::hbarc=197.3269804_wp,me=0.51099895_wp,pi=acos(-1._wp)
    real(wp)::ne,kf,x,bracket,pref
    pressure=0._wp
    if(nb<=0._wp.or.ye<=0._wp)return
    ne=nb*ye
    kf=(3._wp*pi*pi*ne)**(1._wp/3._wp)
    x=hbarc*kf/me
    pref=me**4/(24._wp*pi*pi*hbarc**3)
    if(x<0.05_wp)then
      bracket=x**5*(1.6_wp+x*x*(-4._wp/7._wp+x*x/3._wp))
    else
      bracket=x*(2._wp*x*x-3._wp)*sqrt(1._wp+x*x)+3._wp*asinh(x)
    endif
    pressure=pref*bracket
  end function electron_pressure

  subroutine beta_lepton_fractions(nb,ycharge,ye,ymu,mu_e,ok)
    real(wp),intent(in)::nb,ycharge
    real(wp),intent(out)::ye,ymu,mu_e
    logical,intent(out)::ok
    real(wp),parameter::me=0.51099895_wp,mmu=105.6583745_wp
    real(wp)::lo,hi,mid,flo,fmid
    integer::iter
    ok=.false.;ye=0._wp;ymu=0._wp;mu_e=me
    if(nb<=0._wp.or.ycharge<0._wp.or.ycharge>1._wp)return
    if(ycharge<=tiny(1._wp))then
      ok=.true.
      return
    endif
    mu_e=lepton_chemical_potential(ycharge,nb,me)
    if(mu_e<=mmu)then
      ye=ycharge;ymu=0._wp;ok=.true.;return
    endif
    lo=0._wp;hi=ycharge
    flo=mu_e-mmu
    do iter=1,100
      mid=0.5_wp*(lo+hi)
      fmid=lepton_chemical_potential(ycharge-mid,nb,me) &
        -lepton_chemical_potential(mid,nb,mmu)
      if(fmid>0._wp)then
        lo=mid
      else
        hi=mid
      endif
      if(hi-lo<=2.e-14_wp*max(1._wp,ycharge))exit
    enddo
    ymu=0.5_wp*(lo+hi)
    ye=ycharge-ymu
    mu_e=lepton_chemical_potential(ye,nb,me)
    ok=ieee_is_finite(ye+ymu+mu_e).and.ye>=0._wp.and.ymu>=0._wp
  end subroutine beta_lepton_fractions

  subroutine reconstruct_master(eos,smooth,smooth_sigma,fix_dissolution,dissolution_tol, &
                                recompute_extra_cs2,out,outer_crust_corrected,status)
    type(eos_data_t),intent(in)::eos
    logical,intent(in)::smooth,fix_dissolution,recompute_extra_cs2
    real(wp),intent(in)::smooth_sigma,dissolution_tol
    type(eos_data_t),intent(out)::out
    logical,intent(out)::outer_crust_corrected
    character(len=*),intent(out)::status
    integer,parameter::dense_per_interval=32
    real(wp),parameter::me=0.51099895_wp,alpha=1._wp/137.035999084_wp
    real(wp),parameter::hbarc=197.3269804_wp,pi=acos(-1._wp)
    real(wp),allocatable::x(:),fraw(:),fsm(:),fslope(:),q1(:),dqdx(:),cs2(:)
    real(wp),allocatable::xbridge(:),qbridge(:),bslopes(:),za(:),ycharge(:)
    real(wp)::nq,yeval,ymuval,mueval,t,resid,yecand,yeprev,yenext,yqlo,yqhi,neval
    real(wp)::kf,pel,elat,plat,qval,dqval,fval,muval,dpdn,cs,second,integral
    real(wp)::xq
    integer::n,i,j,k,last_crust,ilo,ihi,idrip,ijoin,nbri,offset_int
    logical::ok,anom_drop,has_comp,has_outer,physical
    character(len=256)::why

    n=eos%n
    status='uninitialized'
    outer_crust_corrected=.false.
    if(n<5)then
      status='reconstruct_master requires at least five source density points'
      return
    endif
    if(.not.allocated(eos%nb).or..not.allocated(eos%q).or..not.allocated(eos%q_extra) &
       .or..not.allocated(eos%iphase))then
      status='EOS primary arrays are not allocated'
      return
    endif
    if(size(eos%nb)/=n.or.any(shape(eos%q)/=[7,n]).or.size(eos%q_extra)/=n &
       .or.size(eos%iphase)/=n)then
      status='EOS primary array shapes are inconsistent'
      return
    endif
    if(.not.all(ieee_is_finite(eos%nb)).or..not.all(ieee_is_finite(eos%q)) &
       .or..not.all(ieee_is_finite(eos%q_extra)))then
      status='EOS source contains nonfinite primary values'
      return
    endif
    if(any(eos%nb<=0._wp).or.any(eos%nb(2:)<=eos%nb(:n-1)))then
      status='reconstruct_master requires a positive strictly increasing density grid'
      return
    endif
    if(abs(eos%temperature_mev)>1.e-12_wp.or..not.ieee_is_finite(eos%temperature_mev))then
      status='reconstruct_master requires a cold T=0 table'
      return
    endif
    if(.not.eos%beta_equilibrated)then
      status='reconstruct_master requires beta-equilibrated source metadata'
      return
    endif
    ! A phase label change is not by itself a thermodynamic discontinuity.
    ! Reject only an equal-pressure phase boundary, which is the signature of
    ! a density plateau/jump that a single-valued F(ln n_B) spline must not cross.
    do i=1,n-1
      if(eos%iphase(i+1)/=eos%iphase(i))then
        if(abs(eos%nb(i+1)*eos%q(1,i+1)-eos%nb(i)*eos%q(1,i)) <= &
          1.e-12_wp*max(1._wp,abs(eos%nb(i)*eos%q(1,i)),abs(eos%nb(i+1)*eos%q(1,i+1))))then
          status='reconstruct_master v1 refuses a pressure plateau across an Iphase boundary'
          return
        endif
      endif
    enddo
    if(eos%nadd<0.or.eos%nadd>1)then
      status='reconstruct_master supports Nadd=0 or 1 only'
      return
    endif
    if(recompute_extra_cs2.and.eos%nadd/=1)then
      status='recomputing q_extra as c_s^2 requires exactly one extra thermo column (Nadd=1)'
      return
    endif
    if(dissolution_tol<0._wp.or..not.ieee_is_finite(dissolution_tol))then
      status='dissolution_tol must be finite and non-negative'
      return
    endif
    if(smooth.and.(smooth_sigma<=0._wp.or.smooth_sigma>0.5_wp.or..not.ieee_is_finite(smooth_sigma)))then
      status='smooth_sigma must lie in (0, 0.5]'
      return
    endif
    has_comp=eos%composition_available
    if(has_comp.and..not.eos%composition_file_present)then
      status='composition_available is true but eos.compo provenance is absent'
      return
    endif
    if(has_comp)then
      if(.not.allocated(eos%yp).or..not.allocated(eos%ye).or..not.allocated(eos%ymu) &
        .or..not.allocated(eos%abar).or..not.allocated(eos%zbar).or..not.allocated(eos%yh) &
        .or..not.allocated(eos%yn))then
        status='composition arrays are not allocated'
        return
      endif
      if(size(eos%yp)/=n.or.size(eos%ye)/=n.or.size(eos%ymu)/=n.or.size(eos%abar)/=n &
        .or.size(eos%zbar)/=n.or.size(eos%yh)/=n.or.size(eos%yn)/=n)then
        status='composition array shapes are inconsistent'
        return
      endif
      if(any(.not.ieee_is_finite(eos%yp)).or.any(.not.ieee_is_finite(eos%ymu)) &
        .or.any(.not.ieee_is_finite(eos%abar)).or.any(.not.ieee_is_finite(eos%zbar)) &
        .or.any(.not.ieee_is_finite(eos%yh)).or.any(.not.ieee_is_finite(eos%ye)))then
        status='composition arrays contain nonfinite values'
        return
      endif
    endif

    allocate(x(n),fraw(n),fsm(n),fslope(n),q1(n),dqdx(n),cs2(n),za(n),ycharge(n))
    x=log(eos%nb)
    fraw=eos_free_energy(eos)
    if(any(.not.ieee_is_finite(fraw)).or.any(fraw<=0._wp))then
      status='master free energy per baryon must be finite and positive'
      return
    endif
    if(smooth)then
      call smooth_master_samples(x,fraw,smooth_sigma,fsm,ok,why)
      if(.not.ok)then
        status=trim(why)
        return
      endif
    else
      fsm=fraw
    endif
    call pchip_slopes(x,fsm,fslope,ok,why)
    if(.not.ok)then
      status='master-function PCHIP failed: '//trim(why)
      return
    endif

    out=eos
    if(has_comp)then
      out%yp=max(0._wp,min(1._wp,eos%yp))
      out%ymu=max(0._wp,min(1._wp,eos%ymu))
      out%yh=max(0._wp,min(1._wp,eos%yh))
      out%abar=eos%abar
      out%zbar=eos%zbar
      out%iphase=eos%iphase
      if(fix_dissolution)then
        last_crust=0
        do i=1,n
          if(out%abar(i)>0._wp)last_crust=i
        enddo
        if(last_crust>0.and.last_crust<n)then
          if(out%abar(last_crust+1)<=0._wp)then
            ilo=last_crust-1;ihi=last_crust+1
            za=0._wp
            do i=1,n
              if(out%abar(i)>0._wp)za(i)=out%zbar(i)/out%abar(i)
            enddo
            resid=abs(out%yp(last_crust)+za(last_crust)*out%yh(last_crust) &
              -(eos%ye(last_crust)+out%ymu(last_crust)))
            anom_drop=.false.
            if(ilo>=1.and.ihi<=n)then
              if(out%abar(ihi)<=0._wp)then
                yecand=out%yp(last_crust)+za(last_crust)*out%yh(last_crust)-out%ymu(last_crust)
                yeprev=out%yp(ilo)+za(ilo)*out%yh(ilo)-out%ymu(ilo)
                yenext=out%yp(ihi)-out%ymu(ihi)
                if(yecand<0.5_wp*min(yeprev,yenext))anom_drop=.true.
              endif
            endif
            if((resid>dissolution_tol.or.anom_drop).and.last_crust+2<=n.and.ilo>=1)then
              if(out%abar(ihi)<=0._wp.and.out%abar(ilo)>0._wp)then
                yqlo=za(ilo)*out%yh(ilo)
                yqhi=out%yp(ihi)
                t=(x(last_crust)-x(ilo))/(x(ihi)-x(ilo))
                nq=yqlo*(1._wp-t)+yqhi*t
                out%iphase(last_crust)=out%iphase(ihi)
                out%yh(last_crust)=0._wp
                out%abar(last_crust)=0._wp
                out%zbar(last_crust)=0._wp
                out%yp(last_crust)=nq
                out%ymu(last_crust)=0._wp
              endif
            endif
          endif
        endif
      endif

      out%yn=max(0._wp,min(1._wp,1._wp-out%yp-out%yh))
      out%yh=max(0._wp,min(1._wp,1._wp-out%yp-out%yn))
      do i=1,n
        if(out%abar(i)>0._wp)then
          za(i)=out%zbar(i)/out%abar(i)
        else
          za(i)=0._wp
        endif
      enddo
      ycharge=max(0._wp,min(1._wp,out%yp+za*out%yh))
      do i=1,n
        call beta_lepton_fractions(eos%nb(i),ycharge(i),yeval,ymuval,mueval,ok)
        if(.not.ok)then
          status='electron-muon beta-equilibrium solve failed at a source row'
          return
        endif
        out%ye(i)=yeval
        out%ymu(i)=ymuval
        out%q(4,i)=-mueval/eos%mn
      enddo
      out%composition_available=.true.
    endif

    has_outer=.false.;idrip=0;ijoin=1
    if(has_comp)then
      do i=1,n
        if(out%abar(i)>0._wp.and.out%yn(i)<=1.e-12_wp)idrip=i
      enddo
      if(idrip>5.and.idrip<n)then
        if(maxval(abs(out%ymu(:idrip)))<=1.e-8_wp)then
          if(all(out%abar(:idrip)>0._wp).and.all(out%yn(:idrip)<=1.e-12_wp))then
            offset_int=max(8,int(0.025_wp*real(n,wp)))
            ijoin=min(n,idrip+offset_int)
            if(ijoin>idrip)has_outer=.true.
          endif
        endif
      endif
    endif

    if(has_outer)then
      outer_crust_corrected=.true.
      nbri=idrip+(n-ijoin+1)
      allocate(xbridge(nbri),qbridge(nbri),bslopes(nbri))
      do i=1,idrip
        xbridge(i)=x(i)
        neval=eos%nb(i)*out%ye(i)
        kf=(3._wp*pi*pi*max(neval,0._wp))**(1._wp/3._wp)
        pel=electron_pressure(eos%nb(i),out%ye(i))
        elat=-0.9_wp*(4._wp/(9._wp*pi))**(1._wp/3._wp)*alpha &
          *max(out%zbar(i),0._wp)**(2._wp/3._wp)*hbarc*kf*neval
        plat=elat/3._wp
        qbridge(i)=(pel+plat)/eos%nb(i)
      enddo
      do i=ijoin,n
        j=idrip+i-ijoin+1
        xbridge(j)=x(i)
        call pchip_value(x,fsm,fslope,x(i),qval,dqval)
        qbridge(j)=dqval
      enddo
      if(any(qbridge<=0._wp).or.any(.not.ieee_is_finite(qbridge)))then
        status='outer-crust electron+lattice bridge has nonpositive or nonfinite Q1'
        return
      endif
      call pchip_slopes(xbridge,qbridge,bslopes,ok,why)
      if(.not.ok)then
        status='outer-crust Q1 bridge PCHIP failed: '//trim(why)
        return
      endif
      do i=1,ijoin
        call pchip_value(xbridge,qbridge,bslopes,x(i),q1(i),dqdx(i))
        call pchip_integral(xbridge,qbridge,bslopes,x(i),x(ijoin),integral)
        fsm(i)=fsm(ijoin)-integral
      enddo
    endif

    if(.not.has_outer)ijoin=1
    do i=ijoin+1,n
      j=pchip_interval(x,x(i))
      call pchip_value_second(x,fsm,fslope,j,x(i),fsm(i),q1(i),dqdx(i))
    enddo
    if(ijoin==1)then
      j=1
      call pchip_value_second(x,fsm,fslope,j,x(1),fsm(1),q1(1),dqdx(1))
    endif
    do i=1,ijoin
      if(has_outer)cycle
      j=pchip_interval(x,x(i))
      call pchip_value_second(x,fsm,fslope,j,x(i),fsm(i),q1(i),dqdx(i))
    enddo
    if(has_outer)then
      call pchip_value(xbridge,qbridge,bslopes,x(ijoin),q1(ijoin),dqdx(ijoin))
    endif

    ! Inspect the actual interpolating polynomials on a dense mesh.  This is
    ! an explicit physical admissibility gate, not a legacy-parity comparison.
    physical=.true.
    do i=1,n-1
      do k=0,dense_per_interval
        xq=x(i)+(x(i+1)-x(i))*real(k,wp)/real(dense_per_interval,wp)
        if(has_outer.and.i<ijoin)then
          call pchip_value(xbridge,qbridge,bslopes,xq,qval,dqval)
          call pchip_integral(xbridge,qbridge,bslopes,xq,x(ijoin),integral)
          fval=fsm(ijoin)-integral
        else
          j=pchip_interval(x,xq)
          call pchip_value_second(x,fsm,fslope,j,xq,fval,qval,second)
          dqval=second
        endif
        muval=fval+qval
        dpdn=qval+dqval
        if(.not.ieee_is_finite(fval+qval+dqval+muval+dpdn).or. &
          fval<=0._wp.or.qval<=0._wp.or.muval<=0._wp.or.dpdn<=0._wp)then
          physical=.false.
          write(status,'(a,i0,a,6(es13.5,1x))') &
            'master reconstruction thermodynamic positivity failure at interval ',i, &
            ' (n_B,F,Q1,mu_B,dPdn,Fpp)= ',exp(xq),fval,qval,muval,dpdn,dqval
          exit
        endif
        cs=dpdn/muval
        if(.not.ieee_is_finite(cs).or.cs<=0._wp.or.cs>1._wp+1.e-6_wp)then
          physical=.false.
          write(status,'(a,i0,a,6(es13.5,1x))') &
            'master reconstruction causality failure at interval ',i, &
            ' (n_B,F,Q1,mu_B,dPdn,cs2)= ',exp(xq),fval,qval,muval,dpdn,cs
          exit
        endif
      enddo
      if(.not.physical)exit
    enddo
    if(.not.physical)then
      if(index(status,'master reconstruction ') /= 1)then
        status='master reconstruction failed its interpolation-mesh physical admissibility gate'
      endif
      return
    endif

    do i=1,n
      if(.not.has_outer.or.i>ijoin)then
        j=pchip_interval(x,x(i))
        call pchip_value_second(x,fsm,fslope,j,x(i),fsm(i),q1(i),dqdx(i))
      endif
      if(has_outer.and.i<=ijoin)then
        call pchip_value(xbridge,qbridge,bslopes,x(i),q1(i),dqdx(i))
      endif
      muval=fsm(i)+q1(i)
      cs2(i)=(q1(i)+dqdx(i))/muval
    enddo

    out%q(1,:)=q1
    out%q(2,:)=0._wp
    out%q(3,:)=(fsm+q1)/eos%mn-1._wp
    out%q(5,:)=0._wp
    out%q(6,:)=fsm/eos%mn-1._wp
    out%q(7,:)=out%q(6,:)
    if(recompute_extra_cs2)then
      out%q_extra=cs2
    endif
    out%beta_equilibrated=.true.
    status='ok'
  end subroutine reconstruct_master

end module compose_eos_reconstruct
