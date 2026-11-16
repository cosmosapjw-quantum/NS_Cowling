module compose_eos_optimize_native
  !! Native SLSQP reverse-communication adapter. Physics lives in callbacks;
  !! this module neither replaces failed objectives nor certifies physics.
  use compose_eos_types, only: wp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: minimize_slsqp_native, native_optimizer_result_t, native_objective_callback
  type :: native_optimizer_result_t
    logical :: success=.false., objective_valid=.false.
    integer :: mode=-100, iterations=0, evaluations=0
    real(wp) :: objective=huge(1._wp), max_constraint_error=huge(1._wp)
    real(wp) :: max_bound_violation=huge(1._wp)
    character(len=160) :: message='not executed'
  end type
  abstract interface
    subroutine native_objective_callback(x,f,g,ok)
      import wp
      real(wp), intent(in) :: x(:)
      real(wp), intent(out) :: f,g(:)
      logical, intent(out) :: ok
    end subroutine
  end interface
  interface
    subroutine slsqp(m,meq,la,n,x,xl,xu,f,c,g,a,acc,iter,mode,w,l_w,jw,l_jw, &
                     alpha,f0,gs,h1,h2,h3,h4,t,t0,tol,iexact,incons,ireset,itermx,line,n1,n2,n3)
      import wp
      integer :: m,meq,la,n,iter,mode,l_w,l_jw,jw(l_jw)
      integer :: iexact,incons,ireset,itermx,line,n1,n2,n3
      real(wp) :: x(n),xl(n),xu(n),f,c(la),g(n+1),a(la,n+1),acc,w(l_w)
      real(wp) :: alpha,f0,gs,h1,h2,h3,h4,t,t0,tol
    end subroutine
  end interface
contains
  subroutine minimize_slsqp_native(x,lower,upper,aeq,beq,evaluate,ftol,maxiter,result)
    real(wp), intent(inout) :: x(:)
    real(wp), intent(in) :: lower(:),upper(:),aeq(:,:),beq(:),ftol
    procedure(native_objective_callback) :: evaluate
    integer, intent(in) :: maxiter
    type(native_optimizer_result_t), intent(out) :: result
    integer :: n,m,la,n1,mineq,lw,ljw,iter,mode,iexact,incons,ireset,itermx,line,n2,n3,ios
    integer :: requests,maxrequests
    integer(int64) :: n64,m64,n164,mineq64,lw64
    integer, allocatable :: jw(:)
    real(wp), allocatable :: w(:),c(:),g(:),a(:,:),xl(:),xu(:),xeval(:),last_x(:),last_g(:)
    real(wp) :: f,acc,alpha,f0,gs,h1,h2,h3,h4,t,t0,tol,last_f
    logical :: ok,cached
    result=native_optimizer_result_t()
    n=size(x);m=size(beq)
    if(n<1.or.size(lower)/=n.or.size(upper)/=n.or.size(aeq,1)/=m.or.size(aeq,2)/=n) then
      result%message='invalid optimizer array shapes'; return
    endif
    if(m>n) then
      result%mode=2;result%message='more equality constraints than independent variables';return
    endif
    if(.not.all(ieee_is_finite(x)).or..not.all(ieee_is_finite(lower)).or. &
       .not.all(ieee_is_finite(upper)).or..not.all(ieee_is_finite(aeq)).or..not.all(ieee_is_finite(beq))) then
      result%message='nonfinite optimizer input';return
    endif
    if(any(lower>upper).or..not.ieee_is_finite(ftol).or.ftol<=0._wp.or.maxiter<1) then
      result%message='invalid bounds, accuracy or iteration limit';return
    endif
    n64=int(n,int64);m64=int(m,int64);n164=n64+1;mineq64=2*n164
    lw64=(3*n164+m64)*(n164+1)+(n164-m64+1)*(mineq64+2)+2*mineq64+ &
         (n164+mineq64)*(n164-m64)+2*m64+n164+n64*n164/2+2*m64+3*n64+3*n164+1
    if(lw64>huge(lw).or.lw64<1) then
      result%message='SLSQP workspace exceeds supported integer index range';return
    endif
    n1=n+1;mineq=2*n1;lw=int(lw64);ljw=mineq;la=max(1,m)
    allocate(w(lw),jw(ljw),c(la),g(n+1),a(la,n+1),xl(n),xu(n),xeval(n),last_x(n),last_g(n),stat=ios)
    if(ios/=0) then
      result%message='SLSQP workspace allocation failed';return
    endif
    w=0;jw=0;c=0;g=0;a=0;xl=lower;xu=upper
    if(m>0)a(1:m,1:n)=aeq
    ! The same initial box projection as SciPy's legacy SLSQP adapter.
    x=min(max(x,lower),upper);cached=.false.
    iter=maxiter-1;mode=0;acc=ftol
    alpha=0;f0=0;gs=0;h1=0;h2=0;h3=0;h4=0;t=0;t0=0;tol=0
    iexact=0;incons=0;ireset=0;itermx=0;line=0;n1=0;n2=0;n3=0
    call get_values(ok)
    if(.not.ok)return
    requests=0
    maxrequests=max(100,100*min(maxiter,huge(maxrequests)/100))
    do
      call slsqp(m,m,la,n,x,xl,xu,f,c,g,a,acc,iter,mode,w,lw,jw,ljw, &
                 alpha,f0,gs,h1,h2,h3,h4,t,t0,tol,iexact,incons,ireset,itermx,line,n1,n2,n3)
      requests=requests+1
      if(abs(mode)/=1)exit
      if(requests>maxrequests)then
        result%mode=-103;result%message='reverse-communication request limit exceeded'
        result%iterations=iter;call diagnostics();return
      endif
      call get_values(ok)
      if(.not.ok)then
        result%iterations=iter;return
      endif
    enddo
    ! Evaluate/report the returned point, rather than any cached preceding point.
    call get_values(ok)
    if(.not.ok)then
      result%iterations=iter;return
    endif
    x=xeval
    result%mode=mode;result%iterations=iter;result%success=(mode==0)
    call diagnostics()
    select case(mode)
    case(0);result%message='SLSQP optimization terminated successfully; physical checks remain caller-owned'
    case(2);result%message='more equality constraints than independent variables'
    case(3);result%message='more than 3*n iterations in least-squares subproblem'
    case(4);result%message='inequality constraints incompatible'
    case(5);result%message='singular matrix E in least-squares subproblem'
    case(6);result%message='singular matrix C in least-squares subproblem'
    case(7);result%message='rank-deficient equality-constraint subproblem'
    case(8);result%message='positive directional derivative in line search'
    case(9);result%message='iteration limit reached'
    case default;result%message='SLSQP workspace or unrecognized termination code'
    end select
  contains
    subroutine get_values(valid)
      logical,intent(out)::valid
      ! Roundoff-size bound excursions are clipped for objective evaluation,
      ! as in the upstream adapter. No finite penalty substitutes a failed solve.
      xeval=min(max(x,lower),upper)
      valid=.true.
      if(cached)then
        if(all(xeval==last_x))then
          f=last_f;g(1:n)=last_g;g(n+1)=0
          if(m>0)c(1:m)=matmul(aeq,x)-beq
          return
        endif
      endif
      result%objective_valid=.false.
      call evaluate(xeval,f,g(1:n),valid)
      result%evaluations=result%evaluations+1;g(n+1)=0
      if(.not.valid)then
        f=huge(1._wp);g=0._wp
        result%mode=-101;result%message='objective callback reported failure at retained trial point'
        call diagnostics();return
      endif
      if(.not.ieee_is_finite(f).or..not.all(ieee_is_finite(g)))then
        f=huge(1._wp);g=0._wp
        valid=.false.;result%mode=-102;result%message='objective callback returned nonfinite objective or gradient'
        call diagnostics();return
      endif
      if(m>0)c(1:m)=matmul(aeq,x)-beq
      result%objective_valid=.true.
      last_x=xeval;last_f=f;last_g=g(1:n);cached=.true.
    end subroutine
    subroutine diagnostics()
      result%objective=f
      result%max_constraint_error=0._wp
      if(m>0)result%max_constraint_error=maxval(abs(matmul(aeq,x)-beq))
      result%max_bound_violation=max(0._wp,maxval(lower-x),maxval(x-upper))
    end subroutine
  end subroutine
end module
