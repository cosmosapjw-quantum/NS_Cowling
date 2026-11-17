program test_matching_slice_conventions
  use compose_eos_types
  use compose_eos_matching
  implicit none
  type(eos_data_t) :: crust,core,output,reference
  type(match_report_t) :: report,reference_report
  character(len=512) :: status,reference_status
  character(len=16),parameter :: methods(4)=[character(len=16)::'cubic_F','cubic_P','quintic_F','quintic_P']
  integer :: failures,k,nleft,nright
  real(wp) :: residual
  failures=0
  do k=1,size(methods)
    call source(crust,.01_wp,.04_wp,20)
    call source(core,.06_wp,.15_wp,30)
    call match_crust_core(crust,core,trim(methods(k)),40,.true.,reference,reference_report,reference_status)
    if(trim(reference_status)/='ok')then
      print '(a,1x,a,1x,a)','BASELINE_FAIL',trim(methods(k)),trim(reference_status)
      failures=failures+1;cycle
    endif

    call encode_il0(core,.1_wp)
    call match_crust_core(crust,core,trim(methods(k)),40,.true.,output,report,status)
    print '(a,1x,a,1x,a)','MIXED_IL',trim(methods(k)),trim(status)
    if(trim(status)=='ok')failures=failures+1

    call encode_il0(crust,.1_wp)
    crust%iYq=3;core%iYq=9
    call match_crust_core(crust,core,trim(methods(k)),40,.true.,output,report,status)
    print '(a,1x,a,1x,a)','SAME_IL0',trim(methods(k)),trim(status)
    if(trim(status)/='ok')then
      failures=failures+1;cycle
    endif
    print '(a,1x,f12.8,1x,i0)','OUTPUT_YQ_AND_INDEX',output%charge_fraction_grid,output%iYq
    if(abs(output%charge_fraction_grid-.1_wp)>1.e-14_wp.or.output%iYq/=9)failures=failures+1
    if(output%il/=0.or.abs(output%yq_axis(1)-.1_wp)>1.e-14_wp)failures=failures+1
    if(output%n/=reference%n)then
      failures=failures+1
    else
      ! The normalization retains epsilon, P and Euler mu.  Il=0 Q3 itself
      ! stores the unshifted baryon potential, not this Euler combination.
      residual=maxval(abs(output%q(3,:)+output%charge_fraction_grid*output%q(4,:)-reference%q(3,:)))
      print '(a,1x,es24.16)','EULER_MU_RESIDUAL',residual
      if(residual>2.e-12_wp)failures=failures+1
      if(maxval(abs(output%q(6,:)-reference%q(6,:)))>2.e-12_wp)failures=failures+1
      nleft=report%i1+1;nright=core%n-report%i2
      residual=maxval(abs(output%q(3,:nleft)-crust%q(3,:nleft)))
      residual=max(residual,maxval(abs(output%q(3,output%n-nright+1:)-core%q(3,report%i2+1:))))
      print '(a,1x,es24.16)','RETAINED_Q3_RESIDUAL',residual
      if(residual>2.e-12_wp)failures=failures+1
      if(report%source_euler_residual_meV>2.e-10_wp)failures=failures+1
    endif

    call source(core,.06_wp,.15_wp,30)
    call encode_il0(core,.2_wp)
    call match_crust_core(crust,core,trim(methods(k)),40,.true.,output,report,status)
    print '(a,1x,a,1x,a)','DIFFERENT_YQ',trim(methods(k)),trim(status)
    if(trim(status)=='ok')failures=failures+1
  enddo
  print '(a,i0)','FAILURES=',failures
  if(failures/=0)stop 1
contains
  subroutine source(eos,nlo,nhi,n)
    type(eos_data_t),intent(out)::eos
    real(wp),intent(in)::nlo,nhi
    integer,intent(in)::n
    integer::j
    real(wp)::f,q
    call allocate_eos(eos,n)
    eos%mn=939._wp;eos%mp=938._wp;eos%beta_equilibrated=.true.
    do j=1,n
      eos%nb(j)=nlo+(nhi-nlo)*real(j-1,wp)/real(n-1,wp)
      q=100._wp*eos%nb(j);f=eos%mn+q
      eos%q(1,j)=q;eos%q(3,j)=(f+q)/eos%mn-1._wp
      eos%q(6,j)=f/eos%mn-1._wp;eos%q(7,j)=eos%q(6,j)
    enddo
  end subroutine source
  subroutine encode_il0(eos,yq)
    type(eos_data_t),intent(inout)::eos
    real(wp),intent(in)::yq
    eos%il=0;eos%charge_fraction_grid=yq;eos%q(4,:)=-.01_wp
    eos%q(3,:)=eos%q(3,:)-yq*eos%q(4,:)
  end subroutine encode_il0
end program test_matching_slice_conventions
