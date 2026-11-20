! Equilibrium-barotropic polar Cowling shooting in the smooth-star sector.
! Geometric r in km; omega_geom = 2*pi*f/(c in km/s).
! The reference Python TOV background is tabulated against s=sqrt(h_start-h).
! Monotone PCHIP coefficient/slopes are supplied once. No interpolation across
! a Maxwell interface is permitted in this module.
! DVODE uses COMMON workspace: process-local, non-reentrant, serialized in Python.
module compose_eos_cowling_dvode
  use iso_c_binding, only: c_int, c_double, c_funptr, c_associated, c_f_procpointer
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  real(c_double), allocatable, save :: knots(:), coeff(:,:), slopes(:,:)
  integer, save :: ntab=0, rhs_failed=0, last_bin=1
  ! A source-evaluated callback is opt-in, process-local, and only enabled
  ! during a Python-serialized native segment. The default remains PCHIP.
  abstract interface
    subroutine source_state_callback(s, values, status) bind(c)
      import :: c_int, c_double
      real(c_double), value :: s
      real(c_double), intent(out) :: values(5)
      integer(c_int), intent(out) :: status
    end subroutine
  end interface
  procedure(source_state_callback), pointer, save :: active_source => null()
  ! At a reset knot the RHS derivative has two distinct one-sided limits.
  ! Query the limit inside the currently integrated subinterval rather
  ! than evaluating an arbitrary branch at the exact C0 EOS knot.
  logical,save :: source_restart_active=.false.
  real(c_double),save :: restart_start=0._c_double,restart_end=0._c_double
  public :: cowling_native_setup, cowling_native_shot, cowling_native_segment, &
            cowling_native_interface_transfer, cowling_native_set_source_callback, &
            cowling_native_segment_restarts, cowling_native_scan, &
            cowling_native_root
contains
  subroutine cowling_native_set_source_callback(ptr,status) bind(c)
    type(c_funptr), value :: ptr
    integer(c_int),intent(out) :: status
    status=0
    if(c_associated(ptr))then
      call c_f_procpointer(ptr,active_source)
      if(.not.associated(active_source))status=1
    else
      nullify(active_source)
    endif
  end subroutine

  subroutine cowling_native_setup(n,s,values,derivatives,status) bind(c)
    integer(c_int),value::n
    real(c_double),intent(in)::s(n),values(5,n),derivatives(5,n)
    integer(c_int),intent(out)::status
    integer::astat
    status=1
    ntab=0
    last_bin=1
    ! Disallow persistence of a callback across independent native setups.
    nullify(active_source)
    source_restart_active=.false.
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
    real(c_double)::dx,z,z2,z3,s_query,boundary_epsilon,orientation
    integer(c_int)::cbstatus
    ok=.false.;v=0._c_double
    if(ntab<4.or..not.ieee_is_finite(s))return
    if(s<knots(1)-1.e-11_c_double.or.s>knots(ntab)+1.e-11_c_double)return
    if(associated(active_source))then
      cbstatus=1
      s_query=s
      if(source_restart_active)then
        ! One-sided limiting value for a derivative jump. A 1e-9 shift
        ! in s is always further restricted to <=1% of this subinterval.
        boundary_epsilon=min(1.e-9_c_double, &
                         0.01_c_double*abs(restart_end-restart_start))
        orientation=sign(1._c_double,restart_end-restart_start)
        if(abs(s-restart_start)<=5.e-12_c_double)then
          s_query=restart_start+orientation*boundary_epsilon
        else if(abs(s-restart_end)<=5.e-12_c_double)then
          s_query=restart_end-orientation*boundary_epsilon
        endif
      endif
      call active_source(s_query,v,cbstatus)
      ok=(cbstatus==0).and.all(ieee_is_finite(v)).and. &
         v(1)>0._c_double.and.v(4)>0._c_double.and.v(5)>0._c_double
      return
    endif
    if(s<=knots(1))then
      v=coeff(:,1)
    else if(s>=knots(ntab))then
      v=coeff(:,ntab)
    else
      ! A native shot evaluates nearby s repeatedly. Reuse the last
      ! valid Hermite cell, with a full binary fallback for large DVODE
      ! steps, arbitrary direction changes, or restart jumps.
      ! DVODE COMMON and these module arrays remain process-serial.
      lo=max(1,min(last_bin,ntab-1))
      if(.not.(knots(lo)<=s.and.s<knots(lo+1)))then
        if(lo>1)then
          if(knots(lo-1)<=s.and.s<knots(lo))lo=lo-1
        endif
        if(lo<ntab-1)then
          if(knots(lo+1)<=s.and.s<knots(lo+2))lo=lo+1
        endif
        if(.not.(knots(lo)<=s.and.s<knots(lo+1)))then
          lo=1;hi=ntab
          do while(hi-lo>1)
            mid=(lo+hi)/2
            if(knots(mid)<=s)then
              lo=mid
            else
              hi=mid
            endif
          enddo
        endif
      endif
      last_bin=lo
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

  ! Native bulk smooth-star scan. All ODE integrations and sign-change
  ! decisions occur in Fortran. The caller only formats the results.
  ! Each frequency starts a fresh DVODE history via cowling_native_shot;
  ! DVODE COMMON and SAVEd coefficient arrays require process isolation
  ! for any concurrent scans (not OpenMP threads).
  subroutine cowling_native_scan(nfreq,freq_hz,ell,rtol,atol,max_ds,mf, &
      match_s,surface_s,residuals,bracket_flags,total_nfev,status) bind(c)
    integer(c_int),value :: nfreq,ell,mf
    real(c_double),intent(in) :: freq_hz(nfreq)
    real(c_double),value :: rtol,atol,max_ds,match_s,surface_s
    real(c_double),intent(out) :: residuals(nfreq)
    integer(c_int),intent(out) :: bracket_flags(nfreq-1),total_nfev,status
    real(c_double),parameter :: pi=3.1415926535897932384626433832795_c_double
    real(c_double),parameter :: c_km_s=299792.458_c_double
    real(c_double) :: si(2),so(2),inner(2,2),outer(2,2)
    real(c_double) :: omega,root_residual
    integer(c_int) :: single_nfev,single_status
    integer :: j
    status=1
    total_nfev=0
    residuals=0._c_double
    bracket_flags=0
    if(nfreq<2.or.ntab<4.or.ell<2)return
    ! A live callback belongs to a separate one-sided sharp interface path.
    ! Never silently promote that phase to a smooth whole-star spectrum.
    if(associated(active_source))then
      status=2
      return
    endif
    if(mf/=10.and.mf/=22)return
    if(.not.all(ieee_is_finite(freq_hz)))return
    if(.not.all(ieee_is_finite([rtol,atol,max_ds,match_s,surface_s])))return
    if(any(freq_hz<=0._c_double))return
    if(any(freq_hz(2:nfreq)<=freq_hz(1:nfreq-1)))return
    if(.not.(0._c_double<match_s.and.match_s<surface_s))return
    if(abs(surface_s-knots(ntab))>1.e-12_c_double)return
    if(.not.(rtol>=1.e-13_c_double.and.rtol<=1.e-3_c_double))return
    if(atol<=0._c_double.or.max_ds<=0._c_double)return
    si=[0._c_double,match_s]
    so=[surface_s,match_s]
    do j=1,nfreq
      omega=2._c_double*pi*freq_hz(j)/c_km_s
      call cowling_native_shot(omega,ell,rtol,atol,max_ds,mf, &
             2_c_int,si,2_c_int,so,inner,outer,root_residual, &
             single_nfev,single_status)
      if(single_status/=0.or..not.ieee_is_finite(root_residual))then
        ! Fail closed, even if earlier shots were valid. No partial scan
        ! can be misinterpreted as a complete physical bracket census.
        residuals=0._c_double
        bracket_flags=0
        total_nfev=0
        status=10+single_status
        return
      endif
      residuals(j)=root_residual
      total_nfev=total_nfev+single_nfev
      if(j>1)then
        if(residuals(j-1)*residuals(j)<0._c_double) &
           bracket_flags(j-1)=1
      endif
    enddo
    status=0
  end subroutine

  ! Evaluate one smooth-star normalized determinant in the native kernel.
  ! Every shot starts fresh DVODE history; the background remains unchanged.
  subroutine cowling_root_shot(f_hz,ell,rtol,atol,max_ds,mf,si,so, &
                                residual,nfev,status)
    real(c_double),intent(in) :: f_hz,rtol,atol,max_ds,si(2),so(2)
    integer(c_int),intent(in) :: ell,mf
    real(c_double),intent(out) :: residual
    integer(c_int),intent(out) :: nfev,status
    real(c_double),parameter :: hz_to_omega= &
       2._c_double*3.1415926535897932384626433832795_c_double/299792.458_c_double
    real(c_double) :: inner(2,2),outer(2,2)
    call cowling_native_shot(f_hz*hz_to_omega,ell,rtol,atol,max_ds,mf, &
        2_c_int,si,2_c_int,so,inner,outer,residual,nfev,status)
  end subroutine

  ! Bracketed Brent-Dekker eigenfrequency polishing, smooth equilibrium EOS.
  ! Unlike a Python-brentq callback, every residual and root iteration is
  ! computed inside Fortran. IEEE double, frequency units Hz; xtol_hz is the
  ! requested bracket uncertainty width, not a residual magnitude.
  ! Fails closed on malformed or unbracketed intervals, DVODE failure,
  ! nonfinite function values, or iteration exhaustion. Never crosses a
  ! sharp phase interface; that requires the separate phasewise solver.
  subroutine cowling_native_root(lo_hz,hi_hz,ell,rtol,atol,max_ds,mf, &
     match_s,surface_s,xtol_hz,frequency_hz,residual,iterations, &
     total_nfev,status) bind(c)
    real(c_double),value :: lo_hz,hi_hz,rtol,atol,max_ds
    real(c_double),value :: match_s,surface_s,xtol_hz
    integer(c_int),value :: ell,mf
    real(c_double),intent(out) :: frequency_hz,residual
    integer(c_int),intent(out) :: iterations,total_nfev,status
    integer,parameter :: max_iterations=128
    real(c_double) :: a,b,c,fa,fb,fc,d,e,s,p,q,r,tol,xmid
    real(c_double) :: si(2),so(2)
    integer(c_int) :: nfev,shot_status
    integer :: iteration

    frequency_hz=0._c_double
    residual=0._c_double
    iterations=0
    total_nfev=0
    status=1
    if(ntab<4.or.ell<2)return
    ! Keep a native C-ABI call from silently treating an active one-sided
    ! sharp EOS as a continuous whole-star fluid.
    if(associated(active_source))then
      status=2
      return
    endif
    if(mf/=10.and.mf/=22)return
    if(.not.all(ieee_is_finite([lo_hz,hi_hz,rtol,atol,max_ds, &
        match_s,surface_s,xtol_hz])))return
    if(lo_hz<=0._c_double.or.hi_hz<=lo_hz)return
    if(xtol_hz<=0._c_double.or.xtol_hz>=hi_hz-lo_hz)return
    if(rtol<1.e-13_c_double.or.rtol>1.e-3_c_double)return
    if(atol<=0._c_double.or.max_ds<=0._c_double)return
    if(match_s<=0._c_double.or.match_s>=surface_s)return
    if(abs(surface_s-knots(ntab))>1.e-12_c_double)return
    si=[0._c_double,match_s]
    so=[surface_s,match_s]
    a=lo_hz
    b=hi_hz
    call cowling_root_shot(a,ell,rtol,atol,max_ds,mf,si,so, &
                           fa,nfev,shot_status)
    if(shot_status/=0.or..not.ieee_is_finite(fa))then
      status=10+shot_status
      return
    endif
    total_nfev=total_nfev+nfev
    call cowling_root_shot(b,ell,rtol,atol,max_ds,mf,si,so, &
                           fb,nfev,shot_status)
    if(shot_status/=0.or..not.ieee_is_finite(fb))then
      total_nfev=0
      status=10+shot_status
      return
    endif
    total_nfev=total_nfev+nfev
    if(fa==0._c_double)then
      frequency_hz=a
      residual=fa
      status=0
      return
    endif
    if(fb==0._c_double)then
      frequency_hz=b
      residual=fb
      status=0
      return
    endif
    if(sign(1._c_double,fa)==sign(1._c_double,fb))then
      total_nfev=0
      status=3
      return
    endif
    c=a
    fc=fa
    d=b-a
    e=d
    do iteration=1,max_iterations
      ! Maintain the best endpoint b and opposite-sign endpoint c.
      if(abs(fc)<abs(fb))then
        a=b; b=c; c=a
        fa=fb; fb=fc; fc=fa
      endif
      tol=2._c_double*epsilon(1._c_double)*abs(b)+0.5_c_double*xtol_hz
      xmid=0.5_c_double*(c-b)
      if(abs(xmid)<=tol.or.fb==0._c_double)then
        ! A sign flip of the normalized matching determinant alone is
        ! not sufficient: an underresolved or unphysical source can
        ! jump across zero without ever approaching a true eigenvalue.
        ! Refuse root certification if the computed determinant remains
        ! large at numerical bracket collapse. This is a conservative
        ! root-admission gate, not a new physical perturbation equation.
        if(abs(fb)>1.e-5_c_double)then
          total_nfev=0
          status=33
          return
        endif
        frequency_hz=b
        residual=fb
        iterations=iteration-1
        status=0
        return
      endif
      if(abs(e)>=tol.and.abs(fa)>abs(fb))then
        s=fb/fa
        if(a==c)then
          p=2._c_double*xmid*s
          q=1._c_double-s
        else
          q=fa/fc
          r=fb/fc
          p=s*(2._c_double*xmid*q*(q-r)-(b-a)*(r-1._c_double))
          q=(q-1._c_double)*(r-1._c_double)*(s-1._c_double)
        endif
        if(p>0._c_double)q=-q
        p=abs(p)
        if(q/=0._c_double.and. &
           2._c_double*p<3._c_double*xmid*q-abs(tol*q).and. &
           p<abs(0.5_c_double*e*q))then
          e=d
          d=p/q
        else
          d=xmid
          e=d
        endif
      else
        d=xmid
        e=d
      endif
      a=b
      fa=fb
      if(abs(d)>tol)then
        b=b+d
      else
        b=b+sign(tol,xmid)
      endif
      call cowling_root_shot(b,ell,rtol,atol,max_ds,mf,si,so, &
                             fb,nfev,shot_status)
      if(shot_status/=0.or..not.ieee_is_finite(fb))then
        total_nfev=0
        status=10+shot_status
        return
      endif
      total_nfev=total_nfev+nfev
      if(sign(1._c_double,fb)==sign(1._c_double,fc))then
        c=a
        fc=fa
        d=b-a
        e=d
      endif
    enddo
    total_nfev=0
    status=32
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

  ! Integrate with FRESH DVODE history at every selected source EOS sound-speed
  ! discontinuity. The caller merges the exact breakpoints into the ordered
  ! output grid and passes their one-based indices. This is opt-in and leaves
  ! existing segment integration and all physical jump conditions untouched.
  subroutine cowling_native_segment_restarts(omega,ell,rtol,atol,max_ds,mf, &
      npts,targets,ncuts,cut_indices,initial,output,nfev,status) bind(c)
    real(c_double),value :: omega,rtol,atol,max_ds
    integer(c_int),value :: ell,mf,npts,ncuts
    real(c_double),intent(in) :: targets(npts),initial(2)
    integer(c_int),intent(in) :: cut_indices(ncuts)
    real(c_double),intent(out) :: output(2,npts)
    integer(c_int),intent(out) :: nfev,status
    real(c_double) :: scratch(2,npts),ystart(2)
    integer :: i,nlocal,first,last,calls,local_status
    status=1;nfev=0;output=0._c_double
    if(ntab<4.or.npts<2.or.ell<2.or.ncuts<0.or.ncuts>=npts-1)return
    if(mf/=10.and.mf/=22)return
    if(.not.all(ieee_is_finite([omega,rtol,atol,max_ds])))return
    if(omega<=0._c_double.or.rtol<=0._c_double.or.atol<=0._c_double.or.max_ds<=0._c_double)return
    if(.not.all(ieee_is_finite(initial)).or..not.all(ieee_is_finite(targets)))return
    if(any(cut_indices<2).or.any(cut_indices>npts-1))return
    if(ncuts>1)then
      if(any(cut_indices(2:ncuts)<=cut_indices(1:ncuts-1)))return
    endif
    ystart=initial;first=1
    do i=1,ncuts+1
      if(i<=ncuts)then
        last=cut_indices(i)
      else
        last=npts
      endif
      nlocal=last-first+1
      if(nlocal<2)then
        status=3;return
      endif
      source_restart_active=.true.
      restart_start=targets(first);restart_end=targets(last)
      call integrate_shot(targets(first),targets(first:last),nlocal,ystart, &
                          omega*omega,ell,rtol,atol,max_ds,mf, &
                          scratch(:,1:nlocal),calls,local_status)
      source_restart_active=.false.
      if(local_status/=0)then
        status=10+local_status;return
      endif
      output(:,first:last)=scratch(:,1:nlocal)
      nfev=nfev+calls
      ystart=scratch(:,nlocal)
      first=last
    enddo
    status=0
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
