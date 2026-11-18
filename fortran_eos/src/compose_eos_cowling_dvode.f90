! Equilibrium-barotropic polar Cowling shooting in the smooth-star sector.
! Geometric r in km; omega_geom = 2*pi*f/(c in km/s).
! The reference Python TOV background is tabulated against s=sqrt(h_start-h).
! Monotone PCHIP coefficient/slopes are supplied once. No interpolation across
! a Maxwell interface is permitted in this module.
! DVODE uses COMMON workspace: process-local, non-reentrant, serialized in Python.
module compose_eos_cowling_dvode
  use iso_c_binding, only: c_int, c_double
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  real(c_double), allocatable, save :: knots(:), coeff(:,:), slopes(:,:)
  integer, save :: ntab=0, rhs_failed=0
  public :: cowling_native_setup, cowling_native_shot, cowling_native_segment, cowling_native_interface_transfer
contains
  subroutine cowling_native_setup(n,s,values,derivatives,status) bind(c)
    integer(c_int),value::n
    real(c_double),intent(in)::s(n),values(5,n),derivatives(5,n)
    integer(c_int),intent(out)::status
    integer::astat
    status=1
    ntab=0
    if(allocated(knots))deallocate(knots)
    if(allocated(coeff))deallocate(coeff)
    if(allocated(slopes))deallocate(slopes)
    if(n<4)return
    if(.not.all(ieee_is_finite(s)))return
    if(.not.all(ieee_is_finite(values)).or..not.all(ieee_is_finite(derivatives)))return
    if(s(1)<0._c_double.or.s(n)<=s(1))return
    if(any(s(2:n)<=s(1:n-1)))return
    ! r, dPhi/dr and sound-speed squared must be physically positive.
    if(any(values(1,:)<=0._c_double).or.any(values(4,:)<=0._c_double) &
       .or.any(values(5,:)<=0._c_double))return
    if(any(values(1,2:n)<=values(1,1:n-1)))return
    allocate(knots(n),coeff(5,n),slopes(5,n),stat=astat)
    if(astat/=0)return
    knots=s;coeff=values;slopes=derivatives;ntab=n;status=0
  end subroutine

  subroutine sample_background(s, v, ok)
    real(c_double),intent(in)::s
    real(c_double),intent(out)::v(5)
    logical,intent(out)::ok
    integer::lo,hi,mid
    real(c_double)::dx,z,z2,z3
    ok=.false.;v=0._c_double
    if(ntab<4.or..not.ieee_is_finite(s))return
    if(s<knots(1)-1.e-11_c_double.or.s>knots(ntab)+1.e-11_c_double)return
    if(s<=knots(1))then
      v=coeff(:,1)
    else if(s>=knots(ntab))then
      v=coeff(:,ntab)
    else
      lo=1;hi=ntab
      do while(hi-lo>1)
        mid=(lo+hi)/2
        if(knots(mid)<=s)then
          lo=mid
        else
          hi=mid
        endif
      enddo
      dx=knots(lo+1)-knots(lo)
      z=(s-knots(lo))/dx;z2=z*z;z3=z2*z
      v=(2*z3-3*z2+1)*coeff(:,lo)+(z3-2*z2+z)*dx*slopes(:,lo) &
         +(-2*z3+3*z2)*coeff(:,lo+1)+(z3-z2)*dx*slopes(:,lo+1)
    endif
    ok=all(ieee_is_finite(v)).and.v(1)>0._c_double &
       .and.v(4)>0._c_double.and.v(5)>0._c_double
  end subroutine

  subroutine ode_rhs(neq,s,y,yds,rpar,ipar)
    integer,intent(in)::neq
    real(c_double),intent(in)::s,y(neq),rpar(*)
    integer,intent(inout)::ipar(*)
    real(c_double),intent(out)::yds(neq)
    real(c_double)::v(5),r,phi,lam,gp,cs,el,om2,ell,dr(2)
    logical::ok
    yds=0._c_double
    if(neq/=2.or.rhs_failed/=0)return
    call sample_background(s,v,ok)
    if(.not.ok.or..not.all(ieee_is_finite(y)))then
      rhs_failed=1;return
    endif
    r=v(1);phi=v(2);lam=v(3);gp=v(4);cs=v(5)
    ell=real(ipar(1),c_double);om2=rpar(1)
    if(lam-2*phi>600._c_double.or.lam>600._c_double)then
      rhs_failed=2;return
    endif
    el=exp(lam)
    dr(1)=(om2*el*exp(-2*phi)*r*y(2)+gp*y(1))/cs &
         -(ell+1._c_double)*(ell*el*y(2)+y(1))/r
    dr(2)=(2*gp-ell/r)*y(2)-el*y(1)/r
    yds=2._c_double*s/gp*dr
    if(.not.all(ieee_is_finite(yds)))then
      rhs_failed=2;yds=0._c_double
    endif
  end subroutine

  subroutine ode_jac(neq,s,y,ml,mu,pd,nrowpd,rpar,ipar)
    integer,intent(in)::neq,ml,mu,nrowpd,ipar(*)
    real(c_double),intent(in)::s,y(neq),rpar(*)
    real(c_double),intent(out)::pd(nrowpd,neq)
    ! MF=10 Adams, MF=22 BDF internal finite-difference Jacobian.
    pd=0._c_double
  end subroutine

  subroutine integrate_shot(s0,targets,npts,y0,om2,ell,rtol,atol,hmax,mf,output,nfev,status)
    integer,intent(in)::npts,ell,mf
    real(c_double),intent(in)::s0,targets(npts),y0(2),om2,rtol,atol,hmax
    real(c_double),intent(out)::output(2,npts)
    integer,intent(out)::nfev,status
    integer,parameter::neq=2,lrw=300,liw=100
    integer::itol,itask,istate,iopt,i,iwork(liw),ipar(1)
    real(c_double)::work(lrw),y(neq),t,rpar(1),abs_tol(2),delta
    external::dvode
    status=3;nfev=0;output=0._c_double
    if(npts<2)return
    delta=targets(npts)-s0
    if(abs(delta)<=tiny(delta))return
    if(abs(targets(1)-s0)>1.e-12_c_double)return
    if(any((targets(2:npts)-targets(1:npts-1))*delta<=0._c_double))return
    y=y0;t=s0;work=0._c_double;iwork=0
    rpar(1)=om2;ipar(1)=ell;abs_tol=atol
    itol=2;itask=4;istate=1;iopt=1
    work(1)=targets(npts)
    work(5)=sign(min(1.e-5_c_double,0.1_c_double*hmax),delta)
    work(6)=hmax;iwork(6)=100000
    output(:,1)=y
    rhs_failed=0
    do i=2,npts
      call dvode(ode_rhs,neq,y,t,targets(i),itol,rtol,abs_tol,itask, &
                 istate,iopt,work,lrw,iwork,liw,ode_jac,mf,rpar,ipar)
      if(rhs_failed/=0.or.istate/=2)then
        status=4;return
      endif
      if(.not.all(ieee_is_finite(y)))then
        status=5;return
      endif
      output(:,i)=y
    enddo
    nfev=iwork(12)
    status=0
  end subroutine

  ! Each trajectory includes both ends: inside 0 -> s_match, outside
  ! s_surface -> s_match. The Python root scan uses just the endpoints.
  subroutine cowling_native_shot(omega,ell,rtol,atol,max_ds,mf,ni,si,no,so, &
                                  inner,outer,residual,nfev,status) bind(c)
    real(c_double),value::omega,rtol,atol,max_ds
    integer(c_int),value::ell,mf,ni,no
    real(c_double),intent(in)::si(ni),so(no)
    real(c_double),intent(out)::inner(2,ni),outer(2,no),residual
    integer(c_int),intent(out)::nfev,status
    real(c_double)::v(5),om2,ys(2),ui(2),uo(2),den
    integer::i1,i2,k1,k2
    logical::ok
    status=1;nfev=0;inner=0._c_double;outer=0._c_double;residual=0._c_double
    if(ntab<4.or.ell<2.or.ni<2.or.no<2)return
    if(mf/=10.and.mf/=22)return
    if(.not.all(ieee_is_finite([omega,rtol,atol,max_ds])))return
    if(omega<=0._c_double.or.rtol<1.e-13_c_double.or.rtol>1.e-3_c_double &
       .or.atol<=0._c_double.or.max_ds<=0._c_double)return
    if(abs(si(1))>1.e-12_c_double.or.abs(so(1)-knots(ntab))>1.e-12_c_double)return
    if(si(ni)<=0._c_double.or.abs(si(ni)-so(no))>1.e-12_c_double)return
    if(si(ni)>=knots(ntab))return
    call sample_background(so(1),v,ok)
    if(.not.ok)return
    if(v(3)-2*v(2)>600._c_double)return
    om2=omega*omega
    ys=[1._c_double,-v(4)/(om2*exp(v(3)-2*v(2))*v(1))]
    if(.not.all(ieee_is_finite(ys)))return
    call integrate_shot(0._c_double,si,ni,[1._c_double,-1._c_double/real(ell,c_double)], &
                        om2,ell,rtol,atol,max_ds,mf,inner,i1,k1)
    if(k1/=0)then
      status=10+k1;return
    endif
    call integrate_shot(knots(ntab),so,no,ys,om2,ell,rtol,atol,max_ds,mf,outer,i2,k2)
    if(k2/=0)then
      status=20+k2;return
    endif
    ui=inner(:,ni);uo=outer(:,no)
    den=sqrt(sum(ui*ui))*sqrt(sum(uo*uo))
    if(den<=tiny(den).or..not.ieee_is_finite(den))then
      status=30;return
    endif
    residual=(ui(1)*uo(2)-ui(2)*uo(1))/den
    if(.not.ieee_is_finite(residual))then
      status=31;return
    endif
    nfev=i1+i2;status=0
  end subroutine

  ! One native integration segment on current one-sided background.
  ! Source providers, branch choice, interface transfer remain explicit.
  subroutine cowling_native_segment(omega,ell,rtol,atol,max_ds,mf,npts,targets,initial, &
                                     output,nfev,status) bind(c)
    real(c_double),value :: omega,rtol,atol,max_ds
    integer(c_int),value :: ell,mf,npts
    real(c_double),intent(in) :: targets(npts),initial(2)
    real(c_double),intent(out) :: output(2,npts)
    integer(c_int),intent(out) :: nfev,status
    status=1;nfev=0
    if(ntab<4.or.npts<2.or.ell<2)return
    if(mf/=10.and.mf/=22)return
    if(.not.all(ieee_is_finite([omega,rtol,atol,max_ds])))return
    if(omega<=0._c_double.or.rtol<=0._c_double.or.atol<=0._c_double.or.max_ds<=0._c_double)return
    if(targets(1)<knots(1)-1.e-10_c_double.or.targets(npts)>knots(ntab)+1.e-10_c_double)return
    if(.not.all(ieee_is_finite(initial)))return
    call integrate_shot(targets(1),targets,npts,initial,omega*omega,ell, &
                        rtol,atol,max_ds,mf,output,nfev,status)
  end subroutine

  ! Scaled w = W/r^(ell+1), v = V/r^ell. Traction ∆P is proportional to
  ! H [omega² exp(lambda-2Phi) r v + Phi' w].
  ! Mode 1 (slow): [w]=0; mode 2 (rapid): [v]=0, both with [∆P]=0.
  subroutine cowling_native_interface_transfer(omega,r,phi,lambda,phi_prime, &
      p,epsilon_in,epsilon_out,conversion,y_before,y_after, &
      traction_relative,kinematic_relative,status) bind(c)
    real(c_double),value :: omega,r,phi,lambda,phi_prime,p,epsilon_in,epsilon_out
    integer(c_int),value :: conversion
    real(c_double),intent(in) :: y_before(2)
    real(c_double),intent(out) :: y_after(2),traction_relative,kinematic_relative
    integer(c_int),intent(out) :: status
    real(c_double) :: a,b,h1,h2,q1,q2,check1,check2
    status=1;y_after=0._c_double;traction_relative=0._c_double
    kinematic_relative=0._c_double
    if(conversion/=1.and.conversion/=2)return
    if(.not.all(ieee_is_finite([omega,r,phi,lambda,phi_prime,p,epsilon_in,epsilon_out])))return
    if(.not.all(ieee_is_finite(y_before)))return
    if(r<=0._c_double.or.omega<=0._c_double.or.phi_prime<=0._c_double.or.p<0._c_double)return
    h1=epsilon_in+p;h2=epsilon_out+p
    if(h1<=0._c_double.or.h2<=0._c_double)return
    if(lambda-2._c_double*phi>600._c_double.or.lambda>600._c_double)return
    a=omega*omega*exp(lambda-2._c_double*phi)*r;b=phi_prime
    if(a<=tiny(a))return
    q1=a*y_before(2)+b*y_before(1)
    q2=(h1/h2)*q1
    if(conversion==1)then
      y_after(1)=y_before(1)
      y_after(2)=(q2-b*y_after(1))/a
    else
      y_after(2)=y_before(2)
      y_after(1)=(q2-a*y_after(2))/b
    endif
    if(.not.all(ieee_is_finite(y_after)))return
    check1=h1*q1;check2=h2*(a*y_after(2)+b*y_after(1))
    traction_relative=abs(check1-check2)/max(abs(check1),abs(check2),tiny(1._c_double))
    if(conversion==1)then
      kinematic_relative=abs(y_after(1)-y_before(1))/ &
          max(abs(y_after(1)),abs(y_before(1)),tiny(1._c_double))
    else
      kinematic_relative=abs(y_after(2)-y_before(2))/ &
          max(abs(y_after(2)),abs(y_before(2)),tiny(1._c_double))
    endif
    if(.not.all(ieee_is_finite([traction_relative,kinematic_relative])))return
    if(traction_relative>1.e-9_c_double.or.kinematic_relative>1.e-9_c_double)then
      status=2;return
    endif
    status=0
  end subroutine
end module compose_eos_cowling_dvode
