module compose_eos_reconstruct
  use compose_eos_types
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: fill_density_gap
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
        u=real(j-1,wp)/real(ndense-1,wp);lv=densex(j)
        call hermite_value(x0,x1,f0,f1,q0,q1,lv,ff,dv)
        call hermite_second_derivative(x0,x1,f0,f1,q0,q1,lv,d2v)
        ! Q1=F' and dP/dn_B=Q1+dQ1/dln(n_B)=F'+F''.
        mu=ff+dv;dpdn=dv+d2v
        if(j>=ndense-2)write(*,'(a,i0,6(1x,es24.16))') 'DBG_HERMITE ',j,u,dv,d2v,mu,dpdn,ff
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
end module compose_eos_reconstruct
