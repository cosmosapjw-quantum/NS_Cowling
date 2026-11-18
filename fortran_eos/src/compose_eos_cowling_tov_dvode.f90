! COWLING-NATIVE-02: native TOV interior profile in enthalpy-derived s-coordinate.
! s = sqrt(h_start-h), r and m in geometric km, P and epsilon in km^-2.
! Smooth equilibrium-barotropic source only. EOS coefficients are pretabulated
! from a locked contract; this module does all TOV ODE evolution with DVODE.
! Serial / per-process only: Netlib DVODE uses hidden COMMON working state.
module compose_eos_cowling_tov_dvode
  use iso_c_binding, only: c_int, c_double
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  real(c_double),parameter :: pi=3.1415926535897932384626433832795_c_double
  real(c_double),allocatable,save :: ss(:), pvals(:), evals(:), pd(:), ed(:)
  integer,save :: ntable=0, rhs_bad=0
  public :: tov_native_profile
contains
  subroutine interp_eos(s,p,e,ok)
    real(c_double),intent(in)::s
    real(c_double),intent(out)::p,e
    logical,intent(out)::ok
    integer::lo,hi,mid
    real(c_double)::dx,z,z2,z3,ax,ay,bx,by
    p=0._c_double;e=0._c_double;ok=.false.
    if(ntable<4.or..not.ieee_is_finite(s))return
    if(s<ss(1)-1.e-10_c_double.or.s>ss(ntable)+1.e-10_c_double)return
    if(s<=ss(1))then
      p=pvals(1);e=evals(1)
    else if(s>=ss(ntable))then
      p=pvals(ntable);e=evals(ntable)
    else
      lo=1;hi=ntable
      do while(hi-lo>1)
        mid=(lo+hi)/2
        if(ss(mid)<=s)then
          lo=mid
        else
          hi=mid
        endif
      enddo
      dx=ss(lo+1)-ss(lo)
      z=(s-ss(lo))/dx;z2=z*z;z3=z2*z
      ax=2*z3-3*z2+1
      ay=(z3-2*z2+z)*dx
      bx=-2*z3+3*z2
      by=(z3-z2)*dx
      p=ax*pvals(lo)+ay*pd(lo)+bx*pvals(lo+1)+by*pd(lo+1)
      e=ax*evals(lo)+ay*ed(lo)+bx*evals(lo+1)+by*ed(lo+1)
    endif
    ok=ieee_is_finite(p).and.ieee_is_finite(e).and.p>=0._c_double.and.e>=0._c_double
  end subroutine

  subroutine tov_rhs_native(neq,s,y,ydot,rpar,ipar)
    integer,intent(in)::neq
    real(c_double),intent(in)::s,y(neq),rpar(*)
    integer,intent(inout)::ipar(*)
    real(c_double),intent(out)::ydot(neq)
    real(c_double)::r,m,p,e,denominator,ds_factor
    logical::ok
    ydot=0._c_double
    if(rhs_bad/=0)return
    if(neq/=2.or..not.all(ieee_is_finite(y)))then
      rhs_bad=1;return
    endif
    r=y(1);m=y(2)
    if(r<=0._c_double.or.m<=0._c_double.or.2._c_double*m>=r)then
      rhs_bad=2;return
    endif
    call interp_eos(s,p,e,ok)
    if(.not.ok)then
      rhs_bad=3;return
    endif
    denominator=m+4._c_double*pi*r**3*p
    if(denominator<=0._c_double)then
      rhs_bad=4;return
    endif
    ds_factor=2._c_double*s*r*(r-2._c_double*m)/denominator
    ydot(1)=ds_factor
    ydot(2)=4._c_double*pi*r*r*e*ds_factor
    if(.not.all(ieee_is_finite(ydot)))then
      rhs_bad=5;ydot=0._c_double
    endif
  end subroutine

  subroutine dummy_jac(neq,s,y,ml,mu,pdout,nrowpd,rpar,ipar)
    integer,intent(in)::neq,ml,mu,nrowpd,ipar(*)
    real(c_double),intent(in)::s,y(neq),rpar(*)
    real(c_double),intent(out)::pdout(nrowpd,neq)
    ! MF=10 Adams or MF=22 BDF with internal finite differences.
    pdout=0._c_double
  end subroutine

  subroutine tov_native_profile(n,s,p,e,dp,de,r0,m0,rtol,atol_r,atol_m, &
                                ds_max,mf,profile,nfev,status) bind(c)
    integer(c_int),value::n,mf
    real(c_double),intent(in)::s(n),p(n),e(n),dp(n),de(n)
    real(c_double),value::r0,m0,rtol,atol_r,atol_m,ds_max
    real(c_double),intent(out)::profile(2,n)
    integer(c_int),intent(out)::nfev,status
    integer,parameter::neq=2,lrw=300,liw=100
    integer::astat,i,istate,itask,iopt,itol,ipar(1),iwork(liw)
    real(c_double)::t,y(neq),work(lrw),atol(neq),rpar(1)
    external::dvode
    profile=0._c_double;nfev=0;status=1
    if(n<4.or.(mf/=10.and.mf/=22))return
    if(.not.all(ieee_is_finite(s)).or..not.all(ieee_is_finite(p)) &
       .or..not.all(ieee_is_finite(e)).or..not.all(ieee_is_finite(dp)) &
       .or..not.all(ieee_is_finite(de)))return
    if(.not.all(ieee_is_finite([r0,m0,rtol,atol_r,atol_m,ds_max])))return
    if(r0<=0._c_double.or.m0<=0._c_double.or.2*m0>=r0)return
    if(abs(s(1))>1.e-15_c_double.or.any(s(2:n)<=s(1:n-1)))return
    if(any(p<0._c_double).or.any(e<0._c_double))return
    if(any(p(2:n)>=p(1:n-1)).or.any(e(2:n)>e(1:n-1)*(1._c_double+1.e-10_c_double)))return
    if(rtol<1.e-13_c_double.or.rtol>1.e-3_c_double &
       .or.atol_r<=0._c_double.or.atol_m<=0._c_double.or.ds_max<=0._c_double)return
    if(allocated(ss))deallocate(ss,pvals,evals,pd,ed)
    ntable=0;rhs_bad=0
    allocate(ss(n),pvals(n),evals(n),pd(n),ed(n),stat=astat)
    if(astat/=0)then
      status=9;return
    endif
    ss=s;pvals=p;evals=e;pd=dp;ed=de;ntable=n
    t=s(1);y=[r0,m0];profile(:,1)=y
    ! ITASK=4 enforces TCRIT so the internal DVODE step never evaluates
    ! EOS outside the physical surface_h boundary. (ITASK=1 may overshoot.)
    itol=2;itask=4;istate=1;iopt=1
    work=0._c_double;iwork=0
    work(1)=ss(ntable)
    work(5)=min(1.e-7_c_double,0.1_c_double*ds_max)
    work(6)=ds_max
    iwork(6)=100000
    atol=[atol_r,atol_m];rpar=0._c_double;ipar=0
    do i=2,n
      call dvode(tov_rhs_native,neq,y,t,s(i),itol,rtol,atol,itask,istate, &
                 iopt,work,lrw,iwork,liw,dummy_jac,mf,rpar,ipar)
      if(rhs_bad/=0)then
        status=20+rhs_bad;return
      endif
      if(istate/=2)then
        status=40-istate;return
      endif
      if(.not.all(ieee_is_finite(y)))then
        status=60;return
      endif
      ! Near nearly coincident EOS knots, dense interpolation may round two
      ! positive increments to identical machine values. Reject meaningful
      ! nonmonotonicity while respecting floating-point resolution.
      if(y(1)<profile(1,i-1)-3.e-12_c_double*max(1._c_double,y(1)) &
         .or.y(2)<profile(2,i-1)-3.e-12_c_double*max(1._c_double,y(2)) &
         .or.2._c_double*y(2)>=y(1))then
        status=61;return
      endif
      profile(:,i)=y
    enddo
    nfev=iwork(12)
    status=0
  end subroutine
end module compose_eos_cowling_tov_dvode
