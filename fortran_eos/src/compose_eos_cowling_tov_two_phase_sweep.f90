! C15 original-EOS two-phase TOV numerical-control sweep.
! Python is only the source-EOS/geometry-grid frontend; the PCHIP
! derivatives, BOTH phasewise TOV differential equations, independent
! DVODE histories and Maxwell transfer are entirely native Fortran.
! The stellar mass and radius are geometrized in km. Pressure and
! epsilon inputs have been converted by frontend to km**(-2).
! All process-global DVODE state is NONREENTRANT; independent OS workers
! are required for true parallelism. No OpenMP over this routine.
module compose_eos_cowling_tov_two_phase_sweep
  use iso_c_binding, only:c_double,c_int
  use, intrinsic::ieee_arithmetic, only:ieee_is_finite
  use compose_eos_cowling_tov_dvode, only:tov_native_profile
  use compose_eos_cowling_pchip, only:cowling_native_pchip_slopes
  implicit none
  private
  public :: tov_native_two_phase_sweep
contains
  subroutine tov_native_two_phase_sweep( &
    nctrl,nh,sh,ph,eh,nl,sl,pl,el,r0,m0,rtols,atols,dsmaxes,mfs, &
    yh,yl,nfevh,nfevl,statuses) bind(c)
    integer(c_int),value::nctrl,nh,nl
    real(c_double),intent(in)::sh(nh),ph(nh),eh(nh)
    real(c_double),intent(in)::sl(nl),pl(nl),el(nl)
    real(c_double),value::r0,m0
    real(c_double),intent(in)::rtols(nctrl),atols(nctrl),dsmaxes(nctrl)
    integer(c_int),intent(in)::mfs(nctrl)
    real(c_double),intent(out)::yh(2,nh,nctrl),yl(2,nl,nctrl)
    integer(c_int),intent(out)::nfevh(nctrl),nfevl(nctrl),statuses(nctrl)
    real(c_double),allocatable::therm_h(:,:),therm_l(:,:),der_h(:,:),der_l(:,:)
    real(c_double),allocatable::dp_h(:),de_h(:),dp_l(:),de_l(:)
    real(c_double)::mstart,rstart,relative_pressure
    integer(c_int)::st,ih,il
    integer::i,err
    yh=0._c_double;yl=0._c_double
    nfevh=0;nfevl=0;statuses=1
    if(nctrl<1.or.nctrl>32.or.nh<4.or.nl<4)return
    if(.not.all(ieee_is_finite([r0,m0])).or.r0<=0._c_double &
       .or.m0<=0._c_double.or.2._c_double*m0>=r0)return
    if(.not.all(ieee_is_finite(sh)).or..not.all(ieee_is_finite(sl)))return
    if(.not.all(ieee_is_finite(ph)).or..not.all(ieee_is_finite(eh)) &
       .or..not.all(ieee_is_finite(pl)).or..not.all(ieee_is_finite(el)))return
    if(any(sh(2:nh)<=sh(1:nh-1)).or.any(sl(2:nl)<=sl(1:nl-1)))return
    if(abs(sh(1))>1.e-12_c_double.or. &
       abs(sh(nh)-sl(1))>1.e-11_c_double)return
    if(any(ph<0._c_double).or.any(pl<0._c_double).or. &
       any(eh<=0._c_double).or.any(el<=0._c_double))return
    if(any(ph(2:nh)>=ph(1:nh-1)).or. &
       any(pl(2:nl)>=pl(1:nl-1)))return
    if(eh(nh)<=el(1))return
    relative_pressure=abs(ph(nh)-pl(1))/max(abs(ph(nh)),abs(pl(1)),1.e-30_c_double)
    if(relative_pressure>1.e-8_c_double)return
    allocate(therm_h(5,nh),therm_l(5,nl), &
        der_h(5,nh),der_l(5,nl),stat=err)
    if(err/=0)then
      statuses=9
      return
    endif
    therm_h=0._c_double;therm_l=0._c_double
    therm_h(1,:)=ph;therm_h(2,:)=eh
    therm_l(1,:)=pl;therm_l(2,:)=el
    call cowling_native_pchip_slopes(nh,sh,therm_h,der_h,st)
    if(st/=0_c_int)then
      statuses=20+st
      return
    endif
    call cowling_native_pchip_slopes(nl,sl,therm_l,der_l,st)
    if(st/=0_c_int)then
      statuses=30+st
      return
    endif
    ! Four contiguous physical derivative channels are constructed ONCE
    ! for the entire sweep. der_h(1,:) is strided in column-major 5xN;
    ! passing it to TOV created a temporary at each DVODE profile call.
    dp_h=der_h(1,:);de_h=der_h(2,:)
    dp_l=der_l(1,:);de_l=der_l(2,:)
    deallocate(therm_h,therm_l,der_h,der_l)
    do i=1,nctrl
      if(.not.ieee_is_finite(rtols(i)).or..not.ieee_is_finite(atols(i)) &
         .or..not.ieee_is_finite(dsmaxes(i)))cycle
      if(rtols(i)<1.e-13_c_double.or.rtols(i)>1.e-3_c_double &
         .or.atols(i)<=0._c_double.or.dsmaxes(i)<=0._c_double)cycle
      if(mfs(i)/=10.and.mfs(i)/=22)cycle
      call tov_native_profile(nh,sh,ph,eh,dp_h,de_h, &
         r0,m0,rtols(i),atols(i),min(1.e-15_c_double,atols(i)*1.e-3_c_double), &
         dsmaxes(i),mfs(i),yh(:,:,i),ih,st)
      if(st/=0_c_int)then
        statuses(i)=100+st
        yh(:,:,i)=0._c_double
        cycle
      endif
      rstart=yh(1,nh,i);mstart=yh(2,nh,i)
      if(rstart<=2._c_double*mstart.or.mstart<=0._c_double)then
        statuses(i)=199
        yh(:,:,i)=0._c_double
        cycle
      endif
      call tov_native_profile(nl,sl,pl,el,dp_l,de_l, &
         rstart,mstart,rtols(i),atols(i), &
         min(1.e-15_c_double,atols(i)*1.e-3_c_double), &
         dsmaxes(i),mfs(i),yl(:,:,i),il,st)
      if(st/=0_c_int)then
        statuses(i)=200+st
        yh(:,:,i)=0._c_double;yl(:,:,i)=0._c_double
        cycle
      endif
      if(yl(1,nl,i)<=2._c_double*yl(2,nl,i).or.yl(2,nl,i)<=0._c_double &
         .or.yl(1,nl,i)<=r0)then
        statuses(i)=299
        yh(:,:,i)=0._c_double;yl(:,:,i)=0._c_double
        cycle
      endif
      nfevh(i)=ih;nfevl(i)=il
      statuses(i)=0
    enddo
  end subroutine
end module compose_eos_cowling_tov_two_phase_sweep
