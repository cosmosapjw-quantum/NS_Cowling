program test_matching_thermodynamic_adversarial
  use compose_eos_types
  use compose_eos_matching
  implicit none
  type(eos_data_t)::cr,co,out,reference
  type(match_report_t)::r,rr
  character(len=512)::status,rs
  character(len=64)::method
  real(wp)::u(3),nb(3),p(3),f(3),mu(3),q(3),cs(3),res,dx
  real(wp)::kn(4),vv(4),ss(4),v,d,a
  integer::j,k,failures
  logical::ok
  character(len=64),parameter::methods(5)=[character(len=64):: &
      'cubic_F','cubic_P','quintic_F','quintic_P','first_law_constrained_hermite_pressure_bridge']
  failures=0
  call make_source(cr,0.01_wp,0.04_wp,20)
  call make_source(co,0.06_wp,0.15_wp,30)
  do k=1,4
    method=methods(k)
    call match_crust_core(cr,co,trim(method),40,.true.,reference,rr,rs)
    co%q(3,:)=(co%q(3,:)+1._wp)*(1._wp+1.e-4_wp)-1._wp
    call match_crust_core(cr,co,trim(method),40,.true.,out,r,status)
    print '(a,1x,a,1x,a)', 'Q3_SHIFT',trim(method),trim(status)
    if(trim(rs)/='ok'.or.trim(status)/='ok')then
      failures=failures+1
    else if(out%n/=reference%n)then
      failures=failures+1
    else
      res=maxval(abs(out%q-reference%q))
      print '(a,1x,es24.16)', 'Q3_OUTPUT_RESIDUAL',res
      if(res>1.e-12_wp.or.abs(r%mu2-rr%mu2)>1.e-11_wp)failures=failures+1
    endif
    call make_source(co,0.06_wp,0.15_wp,30)
  enddo
  call make_source(cr,.01_wp,.04_wp,20)
  call make_source(co,.06_wp,.15_wp,30)
  call match_crust_core(cr,co,'quintic_F',40,.true.,reference,rr,rs)
  do j=1,cr%n
    cr%q(3,j)=(1._wp+cr%q(3,j))*cr%mn/950._wp-1._wp
    cr%q(6,j)=(1._wp+cr%q(6,j))*cr%mn/950._wp-1._wp
    cr%q(7,j)=(1._wp+cr%q(7,j))*cr%mn/950._wp-1._wp
  enddo
  cr%mn=950._wp
  call match_crust_core(cr,co,'quintic_F',40,.true.,out,r,status)
  print '(a,1x,a)','REST_MASS_REFERENCE',trim(status)
  if(trim(status)/='ok'.or.trim(rs)/='ok')then
    failures=failures+1
  else if(maxval(abs(reference%q-out%q))>1.e-11_wp)then
    failures=failures+1
  endif
  co%il=0;co%charge_fraction_grid=.1_wp;co%q(4,:)=-.01_wp
  co%q(3,:)=co%q(3,:)-co%charge_fraction_grid*co%q(4,:)
  call match_crust_core(cr,co,'quintic_P',40,.true.,out,r,status)
  print '(a,1x,a)','MIXED_LEPTON_CONTENT_REJECT',trim(status)
  if(trim(status)=='ok')failures=failures+1
  call make_soft_source(cr,.005_wp,.04_wp,41)
  call make_soft_source(co,.045_wp,.12_wp,61)
  do k=1,4
    call match_crust_core(cr,co,trim(methods(k)),40,.true.,out,r,status)
    print '(a,1x,a,1x,a)','DECREASING_Q',trim(methods(k)),trim(status)
    if(trim(status)/='ok')failures=failures+1
  enddo
  call make_soft_source(cr,.005_wp,.04_wp,41)
  call make_soft_source(co,.045_wp,.12_wp,61)
  cr%q(1,:)=10._wp;co%q(1,:)=10._wp
  cr%q(6,:)=(940._wp+10._wp*log(cr%nb/.04_wp))/cr%mn-1._wp
  co%q(6,:)=(940._wp+10._wp*log(co%nb/.04_wp))/co%mn-1._wp
  cr%q(7,:)=cr%q(6,:);co%q(7,:)=co%q(6,:)
  cr%q(3,:)=cr%q(6,:)+10._wp/cr%mn;co%q(3,:)=co%q(6,:)+10._wp/co%mn
  do k=1,4
    call match_crust_core(cr,co,trim(methods(k)),40,.true.,out,r,status)
    print '(a,1x,a,1x,a)','CONSTANT_Q',trim(methods(k)),trim(status)
    if(trim(status)/='ok')failures=failures+1
  enddo
  ! A smooth pressure spline must obey dF/du = Delta log(n) P/n
  ! at interior points, not only at the endpoint used by root finding.
  dx=log(2._wp)
  call make_pressure_spline(2._wp,1._wp,4._wp,4._wp,8._wp,kn,vv,ss)
  do j=1,7
    u(2)=real(j,wp)/10._wp;u(1)=u(2)-1.e-5_wp;u(3)=u(2)+1.e-5_wp
    call spline_eval(kn,vv,ss,u(2),v,d,a)
    res=(pressure_work_integral(2._wp,u(3),1._wp,4._wp,4._wp,8._wp,0._wp,dx)- &
         pressure_work_integral(2._wp,u(1),1._wp,4._wp,4._wp,8._wp,0._wp,dx))/(2.e-5_wp)
    res=abs(res-v/exp(u(2)*dx))
    print '(a,2(1x,es24.16))','PRESSURE_WORK_DERIVATIVE',u(2),res
    if(res>2.e-8_wp)failures=failures+1
  enddo
  call check_pressure_work_knots(failures)
  print '(a,i0)','FAILURES=',failures
  if(failures/=0)stop 1
contains
  subroutine check_pressure_work_knots(errors)
    ! Independent composite GL4 integration tests both sides of the C1 knots,
    ! tiny and wide log-density spans, and queries in nonmonotone order.
    integer,intent(inout)::errors
    real(wp),parameter::nodes(4)=[-.8611363115940526_wp,-.3399810435848563_wp,.3399810435848563_wp,.8611363115940526_wp]
    real(wp),parameter::weights(4)=[.3478548451374539_wp,.6521451548625461_wp,.6521451548625461_wp,.3478548451374539_wp]
    real(wp),parameter::spans(2)=[.002_wp,17.3_wp]
    real(wp)::knots(4),values(4),slopes(4),query(10),upper,width,point,actual,expected,pressure,deriv,anti
    real(wp)::worst,x1
    integer::i,k,panel,q,j
    call make_pressure_spline(2._wp,1._wp,4._wp,4._wp,8._wp,knots,values,slopes)
    query=[1._wp,knots(3),nearest(knots(2),-1._wp),0._wp,nearest(knots(3),1._wp), &
           .1_wp,knots(2),nearest(knots(3),-1._wp),.5_wp,nearest(knots(2),1._wp)]
    worst=0._wp
    x1=log(1.e-9_wp)
    do j=1,size(spans)
      do i=1,size(query)
        expected=0._wp
        do k=1,3
          upper=min(query(i),knots(k+1))
          if(upper<=knots(k))exit
          width=(upper-knots(k))/32._wp
          do panel=0,31
            do q=1,4
              point=knots(k)+(panel+.5_wp*(1._wp+nodes(q)))*width
              call spline_eval(knots,values,slopes,point,pressure,deriv,anti)
              expected=expected+.5_wp*width*weights(q)*pressure/exp(x1+point*spans(j))
            enddo
          enddo
        enddo
        actual=pressure_work_integral(2._wp,query(i),1._wp,4._wp,4._wp,8._wp,x1,spans(j))
        worst=max(worst,abs(actual-expected)/max(1._wp,abs(expected)))
      enddo
    enddo
    print '(a,1x,es24.16)','PRESSURE_WORK_KNOTS_MAX_RELATIVE',worst
    if(worst>2.e-12_wp)errors=errors+1
  end subroutine
  subroutine make_soft_source(eos,nlo,nhi,n)
    type(eos_data_t),intent(out)::eos
    real(wp),intent(in)::nlo,nhi
    integer,intent(in)::n
    integer::i
    real(wp)::x,ff,qq
    call allocate_eos(eos,n)
    eos%mn=939._wp;eos%mp=938._wp;eos%beta_equilibrated=.true.
    do i=1,n
      eos%nb(i)=exp(log(nlo)+log(nhi/nlo)*real(i-1,wp)/real(n-1,wp))
      x=log(eos%nb(i)/.04_wp);qq=10._wp-x;ff=940._wp+10._wp*x-.5_wp*x*x
      eos%q(1,i)=qq;eos%q(3,i)=(ff+qq)/eos%mn-1._wp
      eos%q(6,i)=ff/eos%mn-1._wp;eos%q(7,i)=eos%q(6,i)
    enddo
  end subroutine
  subroutine make_source(eos,nlo,nhi,n)
    type(eos_data_t),intent(out)::eos
    real(wp),intent(in)::nlo,nhi
    integer,intent(in)::n
    integer::i
    real(wp)::ff,qq
    call allocate_eos(eos,n)
    eos%mn=939._wp;eos%mp=938._wp;eos%beta_equilibrated=.true.
    do i=1,n
      eos%nb(i)=nlo+(nhi-nlo)*real(i-1,wp)/real(n-1,wp)
      qq=100._wp*eos%nb(i);ff=eos%mn+qq
      eos%q(1,i)=qq;eos%q(3,i)=(ff+qq)/eos%mn-1._wp
      eos%q(6,i)=ff/eos%mn-1._wp;eos%q(7,i)=eos%q(6,i)
    enddo
  end subroutine
end program
