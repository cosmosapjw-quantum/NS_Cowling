module compose_eos_causality
  use compose_eos_types
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: causal_completion_report_t, causal_complete

  type :: causal_completion_report_t
    integer :: crossing_left_index = 0
    real(wp) :: n_match = 0._wp
    real(wp) :: pressure_match = 0._wp
    real(wp) :: energy_match = 0._wp
    real(wp) :: mu_match = 0._wp
    real(wp) :: source_mu_match = 0._wp
    real(wp) :: mu_match_relative_residual = 0._wp
    real(wp) :: source_max_cs2 = 0._wp
    real(wp) :: cs2_extension = 1._wp
    real(wp) :: n_max = 0._wp
    integer :: n_tail_points = 0
    logical :: composition_omitted = .false.
    logical :: sound_speed_continuous = .true.
  end type causal_completion_report_t

contains

  subroutine pchip_slopes(x,y,m,ok,message)
    real(wp),intent(in) :: x(:),y(:)
    real(wp),intent(out) :: m(:)
    logical,intent(out) :: ok
    character(len=*),intent(out) :: message
    real(wp),allocatable :: h(:),d(:)
    real(wp) :: edge,w1,w2
    integer :: n,i

    n=size(x)
    ok=.false.
    message='invalid PCHIP arrays'
    if(n<2.or.size(y)/=n.or.size(m)/=n) return
    if(any(.not.ieee_is_finite(x)).or.any(.not.ieee_is_finite(y))) then
      message='PCHIP inputs contain nonfinite values'
      return
    endif
    if(any(x(2:)<=x(:n-1))) then
      message='PCHIP abscissae must be strictly increasing'
      return
    endif

    allocate(h(n-1),d(n-1))
    h=x(2:)-x(:n-1)
    d=(y(2:)-y(:n-1))/h
    if(n==2) then
      m=d(1)
      ok=.true.
      message='ok'
      return
    endif

    edge=((2._wp*h(1)+h(2))*d(1)-h(1)*d(2))/(h(1)+h(2))
    if(edge*d(1)<=0._wp) then
      m(1)=0._wp
    else if(d(1)*d(2)<0._wp.and.abs(edge)>3._wp*abs(d(1))) then
      m(1)=3._wp*d(1)
    else
      m(1)=edge
    endif
    do i=2,n-1
      if(d(i-1)*d(i)<=0._wp) then
        m(i)=0._wp
      else
        w1=2._wp*h(i)+h(i-1)
        w2=h(i)+2._wp*h(i-1)
        m(i)=(w1+w2)/(w1/d(i-1)+w2/d(i))
      endif
    enddo
    edge=((2._wp*h(n-1)+h(n-2))*d(n-1)-h(n-1)*d(n-2))/(h(n-1)+h(n-2))
    if(edge*d(n-1)<=0._wp) then
      m(n)=0._wp
    else if(d(n-1)*d(n-2)<0._wp.and.abs(edge)>3._wp*abs(d(n-1))) then
      m(n)=3._wp*d(n-1)
    else
      m(n)=edge
    endif
    ok=all(ieee_is_finite(m))
    if(ok) then
      message='ok'
    else
      message='PCHIP slopes are nonfinite'
    endif
  end subroutine pchip_slopes

  pure subroutine pchip_segment_value(x,y,m,i,xq,value,derivative)
    real(wp),intent(in) :: x(:),y(:),m(:),xq
    integer,intent(in) :: i
    real(wp),intent(out) :: value,derivative
    real(wp) :: h,t,t2,t3

    h=x(i+1)-x(i)
    t=(xq-x(i))/h
    t2=t*t
    t3=t2*t
    value=(2._wp*t3-3._wp*t2+1._wp)*y(i) &
      +(t3-2._wp*t2+t)*h*m(i) &
      +(-2._wp*t3+3._wp*t2)*y(i+1) &
      +(t3-t2)*h*m(i+1)
    derivative=((6._wp*t2-6._wp*t)/h)*y(i) &
      +(3._wp*t2-4._wp*t+1._wp)*m(i) &
      +((-6._wp*t2+6._wp*t)/h)*y(i+1) &
      +(3._wp*t2-2._wp*t)*m(i+1)
  end subroutine pchip_segment_value

  subroutine evaluate_cs2(i,xq,x,logp,loge,slopep,slopee,cs2,ok)
    integer,intent(in) :: i
    real(wp),intent(in) :: xq,x(:),logp(:),loge(:),slopep(:),slopee(:)
    real(wp),intent(out) :: cs2
    logical,intent(out) :: ok
    real(wp) :: lp,le,dp,de

    call pchip_segment_value(x,logp,slopep,i,xq,lp,dp)
    call pchip_segment_value(x,loge,slopee,i,xq,le,de)
    ok=.false.
    cs2=0._wp
    if(.not.ieee_is_finite(lp).or..not.ieee_is_finite(le).or. &
       .not.ieee_is_finite(dp).or..not.ieee_is_finite(de)) return
    if(de<=0._wp) return
    cs2=exp(lp-le)*dp/de
    ok=ieee_is_finite(cs2)
  end subroutine evaluate_cs2

  subroutine causal_complete(eos,cs2_extension,n_max_requested,n_points,out,report,status)
    type(eos_data_t),intent(in) :: eos
    real(wp),intent(in) :: cs2_extension,n_max_requested
    integer,intent(in) :: n_points
    type(eos_data_t),intent(out) :: out
    type(causal_completion_report_t),intent(out) :: report
    character(len=*),intent(out) :: status

    real(wp),allocatable :: x(:),p(:),eps(:),mu_source(:),logp(:),loge(:)
    real(wp),allocatable :: slopep(:),slopee(:),slopemu(:),cs2_knot(:)
    logical,allocatable :: cs2_valid(:)
    real(wp),allocatable :: nb_out(:),p_out(:),eps_out(:),mu_out(:)
    real(wp) :: nmax,fa,fb,fm,a,b,mid,cs2,nt,pt,et,mu_t,mu_interp
    real(wp) :: mismatch,w_t,q,power,free_per_baryon
    real(wp) :: lv,dv,dmu
    integer :: n,i,ilo,ihi,nkeep,j,iter
    logical :: ok,found
    character(len=256) :: why

    status='uninitialized'
    report=causal_completion_report_t()
    n=eos%n
    if(n<4) then
      status='causal completion requires at least four source points'
      return
    endif
    if(.not.ieee_is_finite(cs2_extension).or.cs2_extension<=0._wp.or.cs2_extension>1._wp) then
      status='cs2_extension must be finite and in (0,1]'
      return
    endif
    if(.not.ieee_is_finite(n_max_requested)) then
      status='n_max must be finite; use zero to select the source maximum'
      return
    endif
    if(n_max_requested<0._wp) then
      status='n_max must be non-negative; use zero to select the source maximum'
      return
    endif
    if(n_points<4) then
      status='n_points must be >= 4'
      return
    endif
    if(.not.allocated(eos%nb).or..not.allocated(eos%q).or..not.allocated(eos%iphase)) then
      status='EOS primary arrays are not allocated'
      return
    endif
    if(size(eos%nb)/=n.or.any(shape(eos%q)/=[7,n]).or.size(eos%iphase)/=n) then
      status='EOS primary array shapes are inconsistent'
      return
    endif
    if(.not.all(ieee_is_finite(eos%nb)).or..not.all(ieee_is_finite(eos%q))) then
      status='EOS source contains nonfinite primary values'
      return
    endif

    allocate(x(n),p(n),eps(n),mu_source(n),logp(n),loge(n))
    allocate(slopep(n),slopee(n),slopemu(n),cs2_knot(n),cs2_valid(n))
    if(any(eos%nb<=0._wp)) then
      status='source baryon density must be positive'
      return
    endif
    x=log(eos%nb)
    p=eos_pressure(eos)
    eps=eos_energy(eos)
    mu_source=eos_mu_b(eos)
    if(any(p<=0._wp).or.any(eps<=0._wp).or.any(mu_source<=0._wp)) then
      status='causal completion requires positive pressure, energy density and mu_B'
      return
    endif
    if(any(eos%nb(2:)<=eos%nb(:n-1))) then
      status='source baryon density must be strictly increasing'
      return
    endif
    if(any(p(2:)<=p(:n-1)).or.any(eps(2:)<=eps(:n-1))) then
      status='causal completion requires strictly increasing pressure and energy density; phase plateaus/jumps are not bridged'
      return
    endif

    logp=log(p)
    loge=log(eps)
    call pchip_slopes(x,logp,slopep,ok,why)
    if(.not.ok) then
      status='pressure PCHIP failed: '//trim(why)
      return
    endif
    call pchip_slopes(x,loge,slopee,ok,why)
    if(.not.ok) then
      status='energy-density PCHIP failed: '//trim(why)
      return
    endif
    call pchip_slopes(x,mu_source,slopemu,ok,why)
    if(.not.ok) then
      status='chemical-potential PCHIP failed: '//trim(why)
      return
    endif

    cs2_knot=0._wp
    cs2_valid=.false.
    do i=1,n
      if(slopee(i)>0._wp) then
        cs2_knot(i)=exp(logp(i)-loge(i))*slopep(i)/slopee(i)
        cs2_valid(i)=ieee_is_finite(cs2_knot(i))
      endif
    enddo
    if(.not.any(cs2_valid)) then
      status='no finite energy-density derivative is available for the causal scan'
      return
    endif
    report%source_max_cs2=maxval(pack(cs2_knot,cs2_valid))

    found=.false.
    do i=1,n-1
      if(.not.cs2_valid(i).or..not.cs2_valid(i+1)) cycle
      if(cs2_knot(i)< -1.e-8_wp) then
        status='source EOS is mechanically unstable below its first causal-limit crossing'
        return
      endif
      if(cs2_knot(i)>1._wp.or.cs2_knot(i+1)<=1._wp) cycle
      ! A local PCHIP slope uses neighbouring knots.  Keep the match stencil
      ! inside one phase so a derivative is never blended across a phase edge.
      ilo=max(1,i-1)
      ihi=min(n,i+2)
      if(any(eos%iphase(ilo:ihi)/=eos%iphase(i))) cycle
      if(x(i+1)-x(i)>0.25_wp) cycle
      fa=cs2_knot(i)-1._wp
      fb=cs2_knot(i+1)-1._wp
      if(fa>0._wp.or.fb<=0._wp) cycle

      if(abs(fa)<=1.e-13_wp) then
        mid=x(i)
      else
        a=x(i)
        b=x(i+1)
        mid=0.5_wp*(a+b)
        ok=.true.
        do iter=1,100
          mid=0.5_wp*(a+b)
          call evaluate_cs2(i,mid,x,logp,loge,slopep,slopee,cs2,ok)
          if(.not.ok) exit
          fm=cs2-1._wp
          if(abs(fm)<=1.e-13_wp.or.(b-a)<=1.e-13_wp*max(1._wp,abs(mid))) exit
          if(fm<=0._wp) then
            a=mid
          else
            b=mid
          endif
        enddo
        if(.not.ok) cycle
        call evaluate_cs2(i,mid,x,logp,loge,slopep,slopee,cs2,ok)
        if(.not.ok.or.abs(cs2-1._wp)>1.e-8_wp) cycle
      endif

      nt=exp(mid)
      call pchip_segment_value(x,mu_source,slopemu,i,mid,mu_interp,dmu)
      call pchip_segment_value(x,logp,slopep,i,mid,lv,dv)
      pt=exp(lv)
      call pchip_segment_value(x,loge,slopee,i,mid,lv,dv)
      et=exp(lv)
      mu_t=(et+pt)/nt
      mismatch=abs(mu_t-mu_interp)/max(abs(mu_t),1._wp)
      if(.not.ieee_is_finite(mismatch).or.mismatch>5.e-3_wp) then
        status='first-law chemical-potential mismatch at causal match exceeds 5e-3'
        return
      endif
      found=.true.
      report%crossing_left_index=i
      exit
    enddo

    if(.not.found) then
      status='no resolved first c_s^2=1 crossing found in a single-phase interval with delta ln(n_B) <= 0.25'
      return
    endif

    nmax=eos%nb(n)
    if(n_max_requested>0._wp) nmax=n_max_requested
    if(nmax<=nt) then
      status='n_max must exceed the first causal-limit density'
      return
    endif
    w_t=et+pt
    nkeep=count(eos%nb<nt)
    if(nkeep<3) then
      status='fewer than three source points remain below the causal matching density'
      return
    endif

    allocate(nb_out(nkeep+n_points),p_out(nkeep+n_points),eps_out(nkeep+n_points),mu_out(nkeep+n_points))
    nb_out(:nkeep)=eos%nb(:nkeep)
    p_out(:nkeep)=p(:nkeep)
    eps_out(:nkeep)=eps(:nkeep)
    mu_out(:nkeep)=(eps_out(:nkeep)+p_out(:nkeep))/nb_out(:nkeep)

    do j=1,n_points
      q=real(j-1,wp)/real(n_points-1,wp)
      nb_out(nkeep+j)=nt*(nmax/nt)**q
      q=nb_out(nkeep+j)/nt
      power=q**(1._wp+cs2_extension)
      p_out(nkeep+j)=pt+(cs2_extension/(1._wp+cs2_extension))*w_t*(power-1._wp)
      mu_out(nkeep+j)=mu_t*q**cs2_extension
      eps_out(nkeep+j)=mu_out(nkeep+j)*nb_out(nkeep+j)-p_out(nkeep+j)
    enddo

    if(any(.not.ieee_is_finite(nb_out)).or.any(.not.ieee_is_finite(p_out)).or. &
       any(.not.ieee_is_finite(eps_out)).or.any(.not.ieee_is_finite(mu_out))) then
      status='causal extension generated nonfinite thermodynamic values'
      return
    endif
    if(any(nb_out<=0._wp).or.any(p_out<=0._wp).or.any(eps_out<=0._wp).or.any(mu_out<=0._wp)) then
      status='causal extension generated nonpositive thermodynamic values'
      return
    endif
    if(any(nb_out(2:)<=nb_out(:size(nb_out)-1)).or.any(p_out(2:)<=p_out(:size(p_out)-1)).or. &
       any(eps_out(2:)<=eps_out(:size(eps_out)-1))) then
      status='causal completion is not strictly monotone in n_B, pressure and energy density'
      return
    endif

    call allocate_eos(out,size(nb_out))
    out%mn=eos%mn
    out%mp=eos%mp
    out%il=eos%il
    out%iT=eos%iT
    out%iYq=eos%iYq
    out%inb_start=1
    out%nadd=0
    out%temperature_mev=0._wp
    out%charge_fraction_grid=0._wp
    out%composition_file_present=.false.
    out%composition_available=.false.
    out%beta_equilibrated=.true.
    out%nb=nb_out
    out%q(1,:)=p_out/nb_out
    out%q(2,:)=0._wp
    out%q(3,:)=mu_out/out%mn-1._wp
    out%q(4,:)=0._wp
    out%q(5,:)=0._wp
    free_per_baryon=eps_out/nb_out
    out%q(6,:)=free_per_baryon/out%mn-1._wp
    out%q(7,:)=out%q(6,:)
    out%q_extra=0._wp
    out%iphase=0
    out%yn=0._wp
    out%yp=0._wp
    out%ye=0._wp
    out%ymu=0._wp
    out%abar=0._wp
    out%zbar=0._wp
    out%yh=0._wp

    report%n_match=nt
    report%pressure_match=pt
    report%energy_match=et
    report%mu_match=mu_t
    report%source_mu_match=mu_interp
    report%mu_match_relative_residual=mismatch
    report%cs2_extension=cs2_extension
    report%n_max=nmax
    report%n_tail_points=n_points
    report%composition_omitted=eos%composition_available.or.eos%composition_file_present
    report%sound_speed_continuous=abs(cs2_extension-1._wp)<=1.e-12_wp
    status='ok'
  end subroutine causal_complete

end module compose_eos_causality
