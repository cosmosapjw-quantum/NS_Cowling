! One-sharp-Maxwell-interface, phase-local PCHIP Cowling root polishing.
! All root iterations, DVODE integrations and slow/rapid transfers are Fortran.
! This is deliberately NOT the exact-source/hybrid-cs2 callback backend.
! Existing DVODE COMMON and module-level coefficient buffers are non-reentrant:
! caller must serialize within a process, and use separate processes in parallel.
module compose_eos_cowling_sharp_root
  use iso_c_binding, only: c_int,c_double
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use compose_eos_cowling_dvode, only: cowling_native_setup, &
        cowling_native_segment, cowling_native_interface_transfer, &
        cowling_native_source_callback_is_active
  implicit none
  private
  real(c_double),parameter :: pi=3.1415926535897932384626433832795_c_double
  real(c_double),parameter :: c_km_s=299792.458_c_double
  public :: cowling_native_sharp_onejump_brent
contains
  ! The geometry is common to both sides of a phase jump, while the
  ! microphysical enthalpy density (epsilon+p) is discontinuous.
  ! geom=(r_km,Phi,lambda,PhiPrime,p,epsilon_low,epsilon_high).
  subroutine sharp_onejump_determinant(f_hz,ell,rtol,atol,max_ds,mf, &
      nh,sh,vh,dh,nl,sl,vl,dl,smatch,geom,conversion,residual,nfev,status)
    real(c_double),intent(in) :: f_hz,rtol,atol,max_ds
    integer(c_int),intent(in) :: ell,mf,nh,nl,conversion
    real(c_double),intent(in) :: sh(nh),vh(5,nh),dh(5,nh)
    real(c_double),intent(in) :: sl(nl),vl(5,nl),dl(5,nl)
    real(c_double),intent(in) :: smatch,geom(7)
    real(c_double),intent(out):: residual
    integer(c_int),intent(out) :: nfev,status
    real(c_double) :: omega,om2,asurf,bsurf
    real(c_double) :: yi(2,2),ylo(2,2),yhi(2,2),y0(2),jstate(2)
    real(c_double) :: tr,kr,den,ti(2),tlo(2),thi(2)
    integer(c_int) :: setup_status,istatus,i1,i2,i3
    residual=0._c_double;nfev=0;status=1
    omega=2._c_double*pi*f_hz/c_km_s
    om2=omega*omega
    ti=[sh(1),smatch]
    tlo=[sl(nl),sl(1)]
    thi=[sh(nh),smatch]
    call cowling_native_setup(nh,sh,vh,dh,setup_status)
    if(setup_status/=0)then
      status=11;return
    endif
    call cowling_native_segment(omega,ell,rtol,atol,max_ds,mf, &
       2_c_int,ti,[1._c_double,-1._c_double/real(ell,c_double)], &
       yi,i1,istatus)
    if(istatus/=0)then
      status=12;return
    endif
    call cowling_native_setup(nl,sl,vl,dl,setup_status)
    if(setup_status/=0)then
      status=13;return
    endif
    if(vl(3,nl)-2._c_double*vl(2,nl)>600._c_double)then
      status=14;return
    endif
    asurf=om2*exp(vl(3,nl)-2._c_double*vl(2,nl))*vl(1,nl)
    bsurf=vl(4,nl)
    if(.not.ieee_is_finite(asurf).or.asurf<=tiny(asurf) &
       .or..not.ieee_is_finite(bsurf).or.bsurf<=0._c_double)then
      status=15;return
    endif
    y0=[1._c_double,-bsurf/asurf]
    if(.not.all(ieee_is_finite(y0)))then
      status=16;return
    endif
    call cowling_native_segment(omega,ell,rtol,atol,max_ds,mf, &
       2_c_int,tlo,y0,ylo,i2,istatus)
    if(istatus/=0)then
      status=17;return
    endif
    ! Transfer from low-density outer phase to high-density inner phase
    ! at *exactly* the prevalidated Maxwell coexistence surface.
    call cowling_native_interface_transfer(omega,geom(1),geom(2),geom(3),geom(4), &
       geom(5),geom(6),geom(7),conversion,ylo(:,2),jstate,tr,kr,istatus)
    if(istatus/=0)then
      status=18;return
    endif
    call cowling_native_setup(nh,sh,vh,dh,setup_status)
    if(setup_status/=0)then
      status=19;return
    endif
    call cowling_native_segment(omega,ell,rtol,atol,max_ds,mf, &
       2_c_int,thi,jstate,yhi,i3,istatus)
    if(istatus/=0)then
      status=20;return
    endif
    den=sqrt(sum(yi(:,2)**2))*sqrt(sum(yhi(:,2)**2))
    if(.not.ieee_is_finite(den).or.den<=tiny(den))then
      status=21;return
    endif
    residual=(yi(1,2)*yhi(2,2)-yi(2,2)*yhi(1,2))/den
    if(.not.ieee_is_finite(residual))then
      residual=0._c_double
      status=22
      return
    endif
    nfev=i1+i2+i3
    status=0
  end subroutine

  ! Strictly bracketed Brent-Dekker polish on the phasewise native determinant.
  ! No numerical EOS callbacks, no Python iterations. Underlying PCHIP tables
  ! are supplied as 2 separate, one-sided arrays with a common interface s.
  subroutine cowling_native_sharp_onejump_brent( &
      f_lo,f_hi,ell,rtol,atol,max_ds,mf, &
      nh,sh,vh,dh,nl,sl,vl,dl,smatch,geom,conversion,xtol_hz, &
      root_hz,root_residual,iterations,total_nfev,status) bind(c)
    real(c_double),value :: f_lo,f_hi,rtol,atol,max_ds,smatch,xtol_hz
    integer(c_int),value :: ell,mf,nh,nl,conversion
    real(c_double),intent(in) :: sh(nh),vh(5,nh),dh(5,nh)
    real(c_double),intent(in) :: sl(nl),vl(5,nl),dl(5,nl),geom(7)
    real(c_double),intent(out) :: root_hz,root_residual
    integer(c_int),intent(out) :: iterations,total_nfev,status
    real(c_double) :: a,b,c,fa,fb,fc,d,e,s,p,q,r,tol,xmid
    integer(c_int) :: nfev,shot_status
    integer :: iteration
    root_hz=0._c_double;root_residual=0._c_double
    iterations=0;total_nfev=0;status=1
    ! The source-exact/hybrid callback has a different physical meaning.
    ! Do not silently erase an active callback by setting up PCHIP tables.
    if(cowling_native_source_callback_is_active())then
      status=2
      return
    endif
    if(nh<4.or.nl<4.or.ell<2)return
    if(mf/=10.and.mf/=22)return
    if(conversion/=1.and.conversion/=2)return
    if(.not.all(ieee_is_finite([f_lo,f_hi,rtol,atol,max_ds, &
        smatch,xtol_hz])))return
    if(.not.all(ieee_is_finite(sh)).or..not.all(ieee_is_finite(sl)))return
    if(.not.all(ieee_is_finite(vh)).or..not.all(ieee_is_finite(vl)))return
    if(.not.all(ieee_is_finite(dh)).or..not.all(ieee_is_finite(dl)))return
    if(.not.all(ieee_is_finite(geom)))return
    if(f_lo<=0._c_double.or.f_hi<=f_lo)return
    if(xtol_hz<=0._c_double.or.xtol_hz>=f_hi-f_lo)return
    if(rtol<1.e-13_c_double.or.rtol>1.e-3_c_double)return
    if(atol<=0._c_double.or.max_ds<=0._c_double)return
    if(abs(sh(1))>1.e-12_c_double)return
    if(sh(nh)<=smatch.or.smatch<=0._c_double)return
    if(abs(sl(1)-sh(nh))>1.e-12_c_double)return
    if(sl(nl)<=sl(1))return
    if(any(sh(2:nh)<=sh(1:nh-1)))return
    if(any(sl(2:nl)<=sl(1:nl-1)))return
    if(any(abs(vh(1:4,nh)-vl(1:4,1))> &
       1.e-8_c_double*max(1._c_double,abs(geom(1:4)))))return
    if(any(abs(vh(1:4,nh)-geom(1:4))> &
       1.e-8_c_double*max(1._c_double,abs(geom(1:4)))))return
    if(geom(1)<=0._c_double.or.geom(4)<=0._c_double.or.geom(5)<0._c_double)return
    if(geom(6)<=0._c_double.or.geom(7)<=geom(6))return
    a=f_lo;b=f_hi
    call sharp_onejump_determinant(a,ell,rtol,atol,max_ds,mf, &
       nh,sh,vh,dh,nl,sl,vl,dl,smatch,geom,conversion,fa,nfev,shot_status)
    if(shot_status/=0)then
      status=10+shot_status
      return
    endif
    total_nfev=total_nfev+nfev
    call sharp_onejump_determinant(b,ell,rtol,atol,max_ds,mf, &
       nh,sh,vh,dh,nl,sl,vl,dl,smatch,geom,conversion,fb,nfev,shot_status)
    if(shot_status/=0)then
      total_nfev=0
      status=10+shot_status
      return
    endif
    total_nfev=total_nfev+nfev
    if(fa==0._c_double.or.fb==0._c_double)then
      if(fa==0._c_double)then
        root_hz=a;root_residual=fa
      else
        root_hz=b;root_residual=fb
      endif
      status=0
      return
    endif
    if(sign(1._c_double,fa)==sign(1._c_double,fb))then
      total_nfev=0
      status=3
      return
    endif
    c=a;fc=fa;d=b-a;e=d
    do iteration=1,128
      if(abs(fc)<abs(fb))then
        a=b;b=c;c=a
        fa=fb;fb=fc;fc=fa
      endif
      tol=2._c_double*epsilon(1._c_double)*abs(b)+0.5_c_double*xtol_hz
      xmid=0.5_c_double*(c-b)
      if(abs(xmid)<=tol.or.fb==0._c_double)then
        ! Refuse false sign-flips: determinant must approach actual zero.
        if(abs(fb)>1.e-5_c_double)then
          total_nfev=0
          status=33
          return
        endif
        root_hz=b;root_residual=fb
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
          q=fa/fc;r=fb/fc
          p=s*(2._c_double*xmid*q*(q-r)-(b-a)*(r-1._c_double))
          q=(q-1._c_double)*(r-1._c_double)*(s-1._c_double)
        endif
        if(p>0._c_double)q=-q
        p=abs(p)
        if(q/=0._c_double.and. &
            2._c_double*p<3._c_double*xmid*q-abs(tol*q).and. &
            p<abs(0.5_c_double*e*q))then
          e=d;d=p/q
        else
          d=xmid;e=d
        endif
      else
        d=xmid;e=d
      endif
      a=b;fa=fb
      if(abs(d)>tol)then
        b=b+d
      else
        b=b+sign(tol,xmid)
      endif
      call sharp_onejump_determinant(b,ell,rtol,atol,max_ds,mf, &
         nh,sh,vh,dh,nl,sl,vl,dl,smatch,geom,conversion,fb,nfev,shot_status)
      if(shot_status/=0.or..not.ieee_is_finite(fb))then
        total_nfev=0
        status=10+shot_status
        return
      endif
      total_nfev=total_nfev+nfev
      if(sign(1._c_double,fb)==sign(1._c_double,fc))then
        c=a;fc=fa;d=b-a;e=d
      endif
    enddo
    total_nfev=0
    status=32
  end subroutine
end module compose_eos_cowling_sharp_root
