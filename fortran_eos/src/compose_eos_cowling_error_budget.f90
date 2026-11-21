! C14 native original-EOS eigenfrequency error-budget batch.
! Numerical controls (relative tolerance, maximum step in s, Adams/BDF)
! are varied *without* changing the source EOS, stellar geometry, mode
! bracket, Maxwell conditions, or root stopping criterion.
! Every root, DVODE integration and source-cs2 evaluation is Fortran.
! The untouched VODE COMMON is nonreentrant: host serializes calls,
! independent OS processes are required for parallel star workloads.
module compose_eos_cowling_error_budget
  use iso_c_binding, only:c_double,c_int
  use, intrinsic::ieee_arithmetic,only:ieee_is_finite
  use compose_eos_cowling_sharp_source_root, only: &
    cowling_native_sharp_source_brent
  implicit none
  private
  public :: cowling_native_sharp_source_sweep
contains
  subroutine cowling_native_sharp_source_sweep( &
    nctrl,f_lo,f_hi,ell,atol,hstart,smatch,geom,conversion,xtol_hz,relative_jump, &
    nh,sh,vh,dh,nl,sl,vl,dl,rtols,max_steps,mfs, &
    frequencies,residuals,iterations,nfevs,cut_counts,statuses) bind(c)
    integer(c_int),value :: nctrl,ell,conversion,nh,nl
    real(c_double),value :: f_lo,f_hi,atol,hstart,smatch,xtol_hz,relative_jump
    real(c_double),intent(in) :: geom(7)
    real(c_double),intent(in) :: sh(nh),vh(5,nh),dh(5,nh)
    real(c_double),intent(in) :: sl(nl),vl(5,nl),dl(5,nl)
    real(c_double),intent(in) :: rtols(nctrl),max_steps(nctrl)
    integer(c_int),intent(in) :: mfs(nctrl)
    real(c_double),intent(out)::frequencies(nctrl),residuals(nctrl)
    integer(c_int),intent(out)::iterations(nctrl),nfevs(nctrl),cut_counts(nctrl),statuses(nctrl)
    integer :: j
    real(c_double) :: f,det
    integer(c_int)::niter,ncalls,ncuts,stat
    frequencies=0._c_double
    residuals=0._c_double
    iterations=0;nfevs=0;cut_counts=0;statuses=1
    if(nctrl<1.or.nctrl>10000)return
    if(nh<4.or.nl<4.or.ell<2)return
    if(conversion/=1.and.conversion/=2)return
    if(.not.all(ieee_is_finite([f_lo,f_hi,atol,hstart,smatch,xtol_hz, &
        relative_jump])))return
    if(f_lo<=0._c_double.or.f_hi<=f_lo.or.xtol_hz<=0._c_double)return
    if(.not.all(ieee_is_finite(geom)))return
    if(.not.all(ieee_is_finite(sh)).or..not.all(ieee_is_finite(sl)))return
    if(.not.all(ieee_is_finite(vh)).or..not.all(ieee_is_finite(vl)))return
    if(.not.all(ieee_is_finite(dh)).or..not.all(ieee_is_finite(dl)))return
    do j=1,nctrl
      if(.not.ieee_is_finite(rtols(j)).or..not.ieee_is_finite(max_steps(j)))cycle
      if(rtols(j)<=0._c_double.or.max_steps(j)<=0._c_double)cycle
      if(mfs(j)/=10.and.mfs(j)/=22)cycle
      call cowling_native_sharp_source_brent(f_lo,f_hi,ell, &
        rtols(j),atol,max_steps(j),mfs(j), &
        nh,sh,vh,dh,nl,sl,vl,dl,hstart,smatch,geom, &
        conversion,xtol_hz,relative_jump, &
        f,det,niter,ncalls,ncuts,stat)
      if(stat/=0_c_int)then
        statuses(j)=stat
        cycle
      endif
      if(.not.ieee_is_finite(f).or..not.ieee_is_finite(det))then
        statuses(j)=90
        cycle
      endif
      frequencies(j)=f;residuals(j)=det
      iterations(j)=niter;nfevs(j)=ncalls;cut_counts(j)=ncuts
      statuses(j)=0
    enddo
  end subroutine
end module compose_eos_cowling_error_budget
