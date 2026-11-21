! C16: native source-exact determinant batch frequency scan.
!
! At fixed original one-sided Maxwell star, preselect source
! derivative knots ONE time, then evaluate the original determinant
! entirely inside Fortran for each frequency. All root candidates must
! subsequently pass native Brent residual, independent Maxwell
! matching and nodal mode-identity gates; a sign flip alone is NOT a
! physical mode certificate.
!
! Fortran ODEPACK DVODE has COMMON process-state; this C ABI is
! non-reentrant and must be externally serialized.
module compose_eos_cowling_source_scan
  use iso_c_binding, only:c_int,c_double
  use, intrinsic::ieee_arithmetic, only:ieee_is_finite
  use compose_eos_cowling_sharp_source_root,only:exact_onejump_det
  use compose_eos_cowling_source_exact,only:cowling_native_source_select_restarts
  use compose_eos_cowling_dvode,only:cowling_native_source_callback_is_active
  implicit none
  private
  public :: cowling_native_sharp_source_scan
contains
  subroutine cowling_native_sharp_source_scan( &
      nfreq,frequencies,ell,rtol,atol,max_ds,mf, &
      nh,sh,vh,dh,nl,sl,vl,dl,hstart,smatch,geom,conversion,relative_jump, &
      determinants,nfevs,statuses,source_knots,status) bind(c)
    integer(c_int),value::nfreq,ell,mf,nh,nl,conversion
    real(c_double),intent(in)::frequencies(nfreq)
    real(c_double),value::rtol,atol,max_ds,hstart,smatch,relative_jump
    real(c_double),intent(in)::sh(nh),vh(5,nh),dh(5,nh),sl(nl),vl(5,nl),dl(5,nl),geom(7)
    real(c_double),intent(out)::determinants(nfreq)
    integer(c_int),intent(out)::nfevs(nfreq),statuses(nfreq),source_knots,status
    real(c_double),allocatable::cuts_h(:),cuts_l(:)
    real(c_double)::det
    integer(c_int)::nhc,nlc,shstat,slstat,nfev,shot_stat
    integer::i
    determinants=0._c_double;nfevs=0;statuses=1;source_knots=0;status=1
    if(nfreq<2.or.nfreq>512.or.nh<4.or.nl<4)return
    if(cowling_native_source_callback_is_active())then
      status=2;return
    endif
    if(ell<2.or.(mf/=10.and.mf/=22))return
    if(conversion/=1.and.conversion/=2)return
    if(.not.all(ieee_is_finite(frequencies)).or.any(frequencies<=0._c_double))return
    if(any(frequencies(2:nfreq)<=frequencies(1:nfreq-1)))return
    if(.not.all(ieee_is_finite([rtol,atol,max_ds,hstart,smatch,relative_jump])))return
    if(rtol<1.e-13_c_double.or.rtol>1.e-3_c_double.or.atol<=0._c_double) return
    if(max_ds<=0._c_double.or.relative_jump<=0._c_double.or.relative_jump>=1._c_double)return
    if(hstart<=0._c_double.or.smatch<=sh(1).or.smatch>=sh(nh))return
    if(.not.all(ieee_is_finite(sh)).or..not.all(ieee_is_finite(sl)))return
    if(.not.all(ieee_is_finite(vh)).or..not.all(ieee_is_finite(vl)))return
    if(.not.all(ieee_is_finite(dh)).or..not.all(ieee_is_finite(dl)))return
    if(.not.all(ieee_is_finite(geom)))return
    if(any(sh(2:nh)<=sh(1:nh-1)).or.any(sl(2:nl)<=sl(1:nl-1)))return
    if(abs(sh(1))>1.e-12_c_double.or.abs(sl(1)-sh(nh))>1.e-12_c_double)return
    if(sl(nl)<=sl(1).or.hstart-sl(nl)**2<0._c_double)return
    if(geom(1)<=0._c_double.or.geom(4)<=0._c_double.or.geom(5)<0._c_double)return
    if(geom(6)<=0._c_double.or.geom(7)<=geom(6))return
    if(any(abs(vh(1:4,nh)-geom(1:4))> &
       1.e-8_c_double*max(1._c_double,abs(geom(1:4)))))return
    if(any(abs(vl(1:4,1)-geom(1:4))> &
       1.e-8_c_double*max(1._c_double,abs(geom(1:4)))))return
    allocate(cuts_h(nh),cuts_l(nl))
    call cowling_native_source_select_restarts(nh,hstart, &
        hstart-sh(1)**2,hstart-sh(nh)**2,relative_jump,2_c_int,cuts_h,nhc,shstat)
    if(shstat/=0)then
      status=10+shstat;return
    endif
    call cowling_native_source_select_restarts(nl,hstart, &
        hstart-sl(1)**2,hstart-sl(nl)**2,relative_jump,1_c_int,cuts_l,nlc,slstat)
    if(slstat/=0)then
      status=20+slstat;return
    endif
    source_knots=nhc+nlc
    do i=1,nfreq
      call exact_onejump_det(frequencies(i),ell,rtol,atol,max_ds,mf, &
          nh,sh,vh,dh,nl,sl,vl,dl,hstart,smatch,geom,conversion, &
          cuts_h(1:nhc),cuts_l(1:nlc),det,nfev,shot_stat)
      if(shot_stat/=0)then
        statuses(i)=100+shot_stat
        cycle
      endif
      if(.not.ieee_is_finite(det))then
        statuses(i)=199
        cycle
      endif
      determinants(i)=det
      nfevs(i)=nfev
      statuses(i)=0
    enddo
    status=0
  end subroutine
end module compose_eos_cowling_source_scan
