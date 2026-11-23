program test_matching_search_exhaustive
  use compose_eos_types
  use compose_eos_matching
  implicit none
  type(eos_data_t)::cr,co,out
  type(match_report_t)::r
  character(len=512)::status
  integer::failures,j
  failures=0
  call make_source(cr,.01_wp,.04_wp,20,100._wp)
  call make_source(co,.06_wp,.15_wp,30,100._wp)
  cr%iphase(cr%n)=7
  call match_crust_core(cr,co,'quintic_P',24,.true.,out,r,status)
  print '(a,a)','LATE_ENDPOINT_RETRY: ',trim(status)
  if(trim(status)/='ok')then
    failures=failures+1
  else if(r%i1>=cr%n-1)then
    failures=failures+1
  endif
  call make_source(cr,.01_wp,.04_wp,220,100._wp)
  do j=51,cr%n;cr%iphase(j)=j;enddo
  call match_crust_core(cr,co,'quintic_P',24,.true.,out,r,status)
  print '(a,a)','BEYOND_LAST149: ',trim(status)
  if(trim(status)/='ok')then
    failures=failures+1
  else if(r%i1>=50)then
    failures=failures+1
  endif
  call make_source(cr,.01_wp,.04_wp,20,1._wp)
  call make_source(co,.06_wp,.15_wp,30,1._wp)
  call match_crust_core(cr,co,'cubic_F',24,.true.,out,r,status)
  print '(a,a)','SUB_POINT2_MEV_SPAN: ',trim(status)
  if(trim(status)/='ok')failures=failures+1
  call make_source(cr,.01_wp,.02_wp,12,100._wp)
  call make_source(co,.025_wp,.04_wp,12,100._wp)
  call match_crust_core(cr,co,'quintic_P',24,.true.,out,r,status)
  print '(a,a)','BELOW_POINT045_CORE: ',trim(status)
  if(trim(status)/='ok')failures=failures+1
  if(.not.r%window_fallback_used)failures=failures+1
  call make_source(cr,.01_wp,.04_wp,20,100._wp)
  call make_source(co,.06_wp,.15_wp,30,100._wp)
  call match_crust_core(cr,co,'first_law_constrained_hermite_pressure_bridge',24,.true.,out,r,status)
  print '(a,a)','FIRSTLAW_METHOD_IDENTITY: ',trim(r%method)
  if(trim(status)/='ok'.or.trim(r%method)/='first_law_constrained_hermite_pressure_bridge')failures=failures+1
  call make_source(cr,.01_wp,.04_wp,20,100._wp)
  cr%q(6,1)=-2._wp;cr%q(7,1)=-2._wp;cr%q(3,1)=-2._wp+cr%q(1,1)/cr%mn
  call match_crust_core(cr,co,'cubic_P',24,.true.,out,r,status)
  print '(a,a)','INVALID_RETAINED_ENERGY: ',trim(status)
  if(trim(status)=='ok'.or..not.r%window_search_exhausted.or. &
    r%window_pairs_examined/=r%window_pairs_total)failures=failures+1
  call make_gamma3_source(cr,.01_wp,.04_wp,30)
  call make_gamma3_source(co,.06_wp,.15_wp,50)
  call match_crust_core(cr,co,'quintic_P',24,.true.,out,r,status)
  print '(a,a)','ACAUSAL_RETAINED_CORE: ',trim(status)
  if(trim(status)=='ok')failures=failures+1
  print '(a,i0)','FAILURES=',failures
  if(failures/=0)stop 1
contains
  subroutine make_gamma3_source(eos,nlo,nhi,n)
    type(eos_data_t),intent(out)::eos
    real(wp),intent(in)::nlo,nhi
    integer,intent(in)::n
    integer::i
    real(wp)::ff,qq
    call allocate_eos(eos,n)
    eos%mn=939._wp;eos%mp=938._wp;eos%beta_equilibrated=.true.
    do i=1,n
      eos%nb(i)=nlo+(nhi-nlo)*real(i-1,wp)/real(n-1,wp)
      qq=50000._wp*eos%nb(i)**2;ff=eos%mn+0.5_wp*qq
      eos%q(1,i)=qq;eos%q(3,i)=(ff+qq)/eos%mn-1._wp
      eos%q(6,i)=ff/eos%mn-1._wp;eos%q(7,i)=eos%q(6,i)
    enddo
  end subroutine
  subroutine make_source(eos,nlo,nhi,n,scale)
    type(eos_data_t),intent(out)::eos
    real(wp),intent(in)::nlo,nhi,scale
    integer,intent(in)::n
    integer::i
    real(wp)::ff,qq
    call allocate_eos(eos,n)
    eos%mn=939._wp;eos%mp=938._wp;eos%beta_equilibrated=.true.
    do i=1,n
      eos%nb(i)=nlo+(nhi-nlo)*real(i-1,wp)/real(n-1,wp)
      qq=scale*eos%nb(i);ff=eos%mn+qq
      eos%q(1,i)=qq;eos%q(3,i)=(ff+qq)/eos%mn-1._wp
      eos%q(6,i)=ff/eos%mn-1._wp;eos%q(7,i)=eos%q(6,i)
    enddo
  end subroutine
end program
