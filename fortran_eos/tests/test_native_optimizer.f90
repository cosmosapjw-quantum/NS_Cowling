program test_native_optimizer
  use compose_eos_types, only: wp
  use compose_eos_optimize_native, only: minimize_slsqp_native, native_optimizer_result_t
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  real(wp)::x(2),lo(2),hi(2),aeq(1,2),beq(1),empty_a(0,2),empty_b(0)
  type(native_optimizer_result_t)::r
  integer::callback_kind
  lo=-2;hi=3;aeq(1,:)=[1._wp,1._wp];beq=1; x=[.7_wp,.3_wp];callback_kind=1
  call minimize_slsqp_native(x,lo,hi,aeq,beq,evaluate,1e-12_wp,100,r)
  call check(r%success.and.maxval(abs(x-[0._wp,1._wp]))<1e-8_wp,'equality optimum')
  call check(r%max_constraint_error<1e-10_wp,'equality residual')
  lo=0;hi=2;x=[1._wp,1._wp];callback_kind=2
  call minimize_slsqp_native(x,lo,hi,empty_a,empty_b,evaluate,1e-12_wp,100,r)
  call check(r%success.and.maxval(abs(x))<1e-8_wp,'active lower-bound optimum')
  lo=-2;hi=2;x=[-1.2_wp,1._wp];callback_kind=3
  call minimize_slsqp_native(x,lo,hi,empty_a,empty_b,evaluate,1e-12_wp,200,r)
  call check(r%success.and.maxval(abs(x-1._wp))<2e-6_wp,'Rosenbrock analytic gradient')
  x=[-1.2_wp,1._wp]
  call minimize_slsqp_native(x,lo,hi,empty_a,empty_b,evaluate,1e-12_wp,1,r)
  call check(.not.r%success.and.r%mode==9,'iteration limit retained')
  x=0;callback_kind=4
  call minimize_slsqp_native(x,lo,hi,empty_a,empty_b,evaluate,1e-12_wp,100,r)
  call check(.not.r%success.and.r%mode==-102,'nonfinite callback rejected')
  x=.5_wp;callback_kind=5
  call minimize_slsqp_native(x,lo,hi,empty_a,empty_b,evaluate,1e-12_wp,100,r)
  call check(.not.r%success.and.r%mode==-101.and.all(x==.5_wp),'failed callback state preserved')
  lo=0;hi=1;x=.5_wp;aeq(1,:)=[1._wp,0._wp];beq=2;callback_kind=1
  call minimize_slsqp_native(x,lo,hi,aeq,beq,evaluate,1e-12_wp,30,r)
  call check(.not.r%success.and.r%max_constraint_error>.9_wp,'infeasible equality not successful')
  lo=1;hi=0;x=.5_wp
  call minimize_slsqp_native(x,lo,hi,empty_a,empty_b,evaluate,1e-12_wp,100,r)
  call check(.not.r%success.and.r%mode==-100,'invalid bounds rejected')
  lo=0;hi=10;x=0;callback_kind=6
  call minimize_slsqp_native(x,lo,hi,empty_a,empty_b,evaluate,1e-12_wp,100,r)
  call check(.not.r%success.and.r%mode==-101.and.any(x>1._wp),'failed later trial preserved')
  call check(.not.r%objective_valid,'failed objective marked unavailable')
  call large_box_equality()
  print '(a)','PASS: 10 native optimizer cases'
contains
  subroutine check(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    if(.not.ok) then
      print *,'FAIL ',label,' mode=',r%mode,' x=',x,' msg=',trim(r%message)
      stop 1
    endif
  end subroutine
  subroutine large_box_equality()
    integer,parameter::nn=129
    integer::j
    real(wp)::q(nn),ql(nn),qu(nn),aa(1,nn),bb(1),expected(nn)
    do j=1,nn
      expected(j)=min(.4_wp,max(0._wp,.3_wp*sin(real(j,wp)*.17_wp)+.2_wp-.07_wp))
    enddo
    q=.2_wp;ql=0;qu=.4_wp;aa=1._wp/nn;bb(1)=sum(expected)/nn
    call minimize_slsqp_native(q,ql,qu,aa,bb,large_objective,1e-12_wp,100,r)
    call check(r%success.and.maxval(abs(q-expected))<2e-7_wp,'129-variable active-bound equality solution')
  end subroutine
  subroutine large_objective(z,f,g,ok)
    real(wp),intent(in)::z(:)
    real(wp),intent(out)::f,g(:)
    logical,intent(out)::ok
    integer::j
    real(wp)::target
    f=0;ok=.true.
    do j=1,size(z)
      target=.3_wp*sin(real(j,wp)*.17_wp)+.2_wp
      f=f+(z(j)-target)**2;g(j)=2*(z(j)-target)
    enddo
  end subroutine
  subroutine evaluate(z,f,g,ok)
    real(wp),intent(in)::z(:)
    real(wp),intent(out)::f,g(:)
    logical,intent(out)::ok
    ok=.true.
    select case(callback_kind)
    case(1)
      f=(z(1)-1)**2+(z(2)-2)**2;g=2*(z-[1._wp,2._wp])
    case(2)
      f=sum((z+1)**2);g=2*(z+1)
    case(3)
      f=(1-z(1))**2+100*(z(2)-z(1)**2)**2
      g(1)=2*(z(1)-1)-400*z(1)*(z(2)-z(1)**2);g(2)=200*(z(2)-z(1)**2)
    case(4)
      f=ieee_value(0._wp,ieee_quiet_nan);g=0
    case(5)
      ok=.false.
    case(6)
      if(any(z>1._wp)) then
        ok=.false.;return
      endif
      f=sum((z-5._wp)**2);g=2*(z-5._wp)
    end select
  end subroutine
end program
