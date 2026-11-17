program test_maxwell_reference_normalization
  use compose_eos_types
  use compose_eos_phase_transition
  implicit none
  type(eos_data_t)::low,high,reencoded,reference,result
  type(maxwell_transition_report_t)::rr,r
  character(len=512)::status,reference_status
  real(wp)::ratio,residual,mn_new,mp_new
  integer::leptons,scenario,failures,j
  failures=0
  do leptons=0,1
    call make_phase(low,.false.,leptons)
    call make_phase(high,.true.,leptons)
    call construct_maxwell_transition(low,high,.true.,reference,rr,reference_status)
    if(reference_status/='ok')then
      print '(a,a)','BASELINE_FAILURE=',trim(reference_status)
      stop 2
    endif
    do scenario=1,3
      reencoded=high
      mn_new=high%mn;mp_new=high%mp
      if(scenario==1.or.scenario==3)mn_new=950._wp
      if(scenario==2.or.scenario==3)mp_new=938.27203_wp
      ratio=high%mn/mn_new
      reencoded%q(3,:)=(1._wp+high%q(3,:))*ratio-1._wp
      reencoded%q(6,:)=(1._wp+high%q(6,:))*ratio-1._wp
      reencoded%q(7,:)=(1._wp+high%q(7,:))*ratio-1._wp
      reencoded%q(4,:)=high%q(4,:)*ratio
      reencoded%q(5,:)=high%q(5,:)*ratio
      reencoded%mn=mn_new;reencoded%mp=mp_new
      call construct_maxwell_transition(low,reencoded,.true.,result,r,status)
      print '(a,2(1x,i0),1x,a)','REFERENCE_CASE',leptons,scenario,trim(status)
      if(status/='ok')then
        failures=failures+1;cycle
      endif
      if(result%n/=reference%n)then
        failures=failures+1;cycle
      endif
      residual=maxval(abs(eos_energy(result)-eos_energy(reference)))
      print '(a,es24.16)','EPSILON_RESIDUAL=',residual
      if(residual>2.e-9_wp)failures=failures+1
      if(maxval(abs(eos_pressure(result)-eos_pressure(reference)))>2.e-10_wp)failures=failures+1
      if(maxval(abs(eos_mu_b(result)-eos_mu_b(reference)))>2.e-9_wp)failures=failures+1
      if(maxval(abs(eos_mu_q(result)-eos_mu_q(reference)))>2.e-10_wp)failures=failures+1
      if(abs(r%mu_transition-rr%mu_transition)>2.e-9_wp)failures=failures+1
      if(abs(r%n_high-rr%n_high)>2.e-11_wp.or.abs(r%n_low-rr%n_low)>2.e-11_wp)failures=failures+1
      if(maxval(abs(eos_energy(reencoded)-eos_energy(high)))>2.e-11_wp)failures=failures+1
    enddo
    ! A physical +100 MeV per baryon shift is not a reference re-encoding.
    ! It moves mu support; the kernel must not fit/remove this energy offset.
    reencoded=high
    reencoded%q(3,:)=high%q(3,:)+100._wp/high%mn
    reencoded%q(6,:)=high%q(6,:)+100._wp/high%mn
    reencoded%q(7,:)=high%q(7,:)+100._wp/high%mn
    call construct_maxwell_transition(low,reencoded,.true.,result,r,status)
    print '(a,1x,i0,1x,a)','PHYSICAL_OFFSET',leptons,trim(status)
    if(status=='ok'.or.maxwell_status_category(status)/='NO_COMMON_MU_SUPPORT')failures=failures+1
  enddo
  print '(a,i0)','FAILURES=',failures
  if(failures/=0)stop 1
contains
  subroutine make_phase(eos,upper,il)
    type(eos_data_t),intent(out)::eos
    logical,intent(in)::upper
    integer,intent(in)::il
    real(wp),parameter::mugrid(6)=[950._wp,960._wp,980._wp,1000._wp,1030._wp,1040._wp]
    real(wp)::x,p,n,mu,eps,mq
    integer::i
    call allocate_eos(eos,6)
    eos%mn=939._wp;eos%mp=939._wp;eos%beta_equilibrated=.true.
    eos%il=il;eos%charge_fraction_grid=.1_wp
    do i=1,6
      mu=mugrid(i);x=mu-990._wp
      p=8._wp+.12_wp*x+.001_wp*x*x;n=.12_wp+.002_wp*x
      if(upper)then;p=p+.01_wp*x;n=n+.01_wp;endif
      eps=mu*n-p;mq=-40._wp
      eos%nb(i)=n;eos%q(1,i)=p/n;eos%q(4,i)=mq/eos%mn
      eos%q(3,i)=mu/eos%mn-1._wp
      if(il==0)eos%q(3,i)=eos%q(3,i)-eos%charge_fraction_grid*eos%q(4,i)
      eos%q(6,i)=eps/(n*eos%mn)-1._wp;eos%q(7,i)=eos%q(6,i)
    enddo
  end subroutine
end program
