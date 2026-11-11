module compose_eos_phase_transition
  use compose_eos_types
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: maxwell_transition_report_t, construct_maxwell_transition

  real(wp), parameter :: euler_tolerance = 5.e-4_wp
  real(wp), parameter :: algebraic_tolerance = 2.e-6_wp
  real(wp), parameter :: entropy_tolerance = 2.e-8_wp
  real(wp), parameter :: lepton_potential_tolerance = 1.e-8_wp
  real(wp), parameter :: causal_tolerance = 1.e-8_wp
  real(wp), parameter :: pressure_relative_tolerance = 1.e-10_wp

  type :: maxwell_transition_report_t
    integer :: pressure_crossings = 0
    integer :: low_phase_rows = 0
    integer :: high_phase_rows = 0
    real(wp) :: mu_transition = 0._wp
    real(wp) :: pressure_transition = 0._wp
    real(wp) :: n_low = 0._wp
    real(wp) :: n_high = 0._wp
    real(wp) :: epsilon_low = 0._wp
    real(wp) :: epsilon_high = 0._wp
    real(wp) :: density_jump = 0._wp
    real(wp) :: latent_energy_density = 0._wp
    real(wp) :: maxwell_identity_residual = 0._wp
    real(wp) :: max_euler_residual_low = 0._wp
    real(wp) :: max_euler_residual_high = 0._wp
    real(wp) :: min_cs2_low = huge(1._wp)
    real(wp) :: max_cs2_low = 0._wp
    real(wp) :: min_cs2_high = huge(1._wp)
    real(wp) :: max_cs2_high = 0._wp
  end type maxwell_transition_report_t

contains

  pure subroutine hermite_pressure(mu0,mu1,p0,p1,n0,n1,muq,p,n,dndmu)
    real(wp),intent(in) :: mu0,mu1,p0,p1,n0,n1,muq
    real(wp),intent(out) :: p,n,dndmu
    real(wp) :: h,t,t2,t3
    h=mu1-mu0
    t=(muq-mu0)/h
    t2=t*t
    t3=t2*t
    p=(2._wp*t3-3._wp*t2+1._wp)*p0 &
      +(t3-2._wp*t2+t)*h*n0 &
      +(-2._wp*t3+3._wp*t2)*p1 &
      +(t3-t2)*h*n1
    n=((6._wp*t2-6._wp*t)/h)*p0 &
      +(3._wp*t2-4._wp*t+1._wp)*n0 &
      +((-6._wp*t2+6._wp*t)/h)*p1 &
      +(3._wp*t2-2._wp*t)*n1
    dndmu=((12._wp*t-6._wp)/(h*h))*p0 &
      +((6._wp*t-4._wp)/h)*n0 &
      +((-12._wp*t+6._wp)/(h*h))*p1 &
      +((6._wp*t-2._wp)/h)*n1
  end subroutine hermite_pressure

  integer function interval_index(x,xq) result(i)
    real(wp),intent(in) :: x(:),xq
    integer :: lo,hi,mid,n
    n=size(x)
    if(xq<=x(1)) then
      i=1
      return
    endif
    if(xq>=x(n)) then
      i=n-1
      return
    endif
    lo=1
    hi=n
    do while(hi-lo>1)
      mid=(lo+hi)/2
      if(xq<x(mid)) then
        hi=mid
      else
        lo=mid
      endif
    enddo
    i=lo
  end function interval_index

  subroutine evaluate_phase(eos,mu_grid,p_grid,muq,p,n,dndmu,q4,q5,qextra,ok)
    type(eos_data_t),intent(in) :: eos
    real(wp),intent(in) :: mu_grid(:),p_grid(:),muq
    real(wp),intent(out) :: p,n,dndmu,q4,q5,qextra
    logical,intent(out) :: ok
    integer :: i
    real(wp) :: h,t
    ok=.false.
    p=0._wp
    n=0._wp
    dndmu=0._wp
    q4=0._wp
    q5=0._wp
    qextra=0._wp
    if(muq<mu_grid(1).or.muq>mu_grid(size(mu_grid))) return
    i=interval_index(mu_grid,muq)
    call hermite_pressure(mu_grid(i),mu_grid(i+1),p_grid(i),p_grid(i+1), &
      eos%nb(i),eos%nb(i+1),muq,p,n,dndmu)
    h=mu_grid(i+1)-mu_grid(i)
    t=(muq-mu_grid(i))/h
    q4=(1._wp-t)*eos%q(4,i)+t*eos%q(4,i+1)
    q5=(1._wp-t)*eos%q(5,i)+t*eos%q(5,i+1)
    qextra=(1._wp-t)*eos%q_extra(i)+t*eos%q_extra(i+1)
    ok=ieee_is_finite(p).and.ieee_is_finite(n).and.ieee_is_finite(dndmu) &
      .and.ieee_is_finite(q4).and.ieee_is_finite(q5).and.ieee_is_finite(qextra)
  end subroutine evaluate_phase

  pure subroutine quadratic_interval_minimum(y0,ymid,y1,ymin)
    real(wp),intent(in) :: y0,ymid,y1
    real(wp),intent(out) :: ymin
    real(wp) :: a,b,tv,yv
    a=2._wp*(y1+y0-2._wp*ymid)
    b=y1-y0-a
    ymin=min(y0,y1)
    if(abs(a)>1.e-15_wp*max(1._wp,abs(a),abs(b),abs(y0))) then
      tv=-b/(2._wp*a)
      if(tv>0._wp.and.tv<1._wp) then
        yv=a*tv*tv+b*tv+y0
        ymin=min(ymin,yv)
      endif
    endif
  end subroutine quadratic_interval_minimum

  subroutine validate_phase(eos,mu_grid,p_grid,eps_grid,max_euler,min_cs2,max_cs2,status)
    type(eos_data_t),intent(in) :: eos
    real(wp),allocatable,intent(out) :: mu_grid(:),p_grid(:),eps_grid(:)
    real(wp),intent(out) :: max_euler,min_cs2,max_cs2
    character(len=*),intent(out) :: status
    integer :: npts,i,j
    real(wp) :: mu_source,residual,muq,p,nbval,dndmu,eps,cs2
    real(wp) :: mu0,mu1,mumid,p0,pmid,p1,n0,nmid,n1,d0,dmid,d1
    real(wp) :: h0,hmid,h1,min_n,min_h,min_dn,scale

    status='uninitialized'
    max_euler=0._wp
    min_cs2=huge(1._wp)
    max_cs2=0._wp
    npts=eos%n
    if(npts<3) then
      status='each Maxwell phase requires at least three source points'
      return
    endif
    if(.not.allocated(eos%nb).or..not.allocated(eos%q).or. &
       .not.allocated(eos%q_extra).or..not.allocated(eos%iphase)) then
      status='phase EOS arrays are not allocated'
      return
    endif
    if(size(eos%nb)/=npts.or.any(shape(eos%q)/=[7,npts]).or. &
       size(eos%q_extra)/=npts.or.size(eos%iphase)/=npts) then
      status='phase EOS array shapes are inconsistent'
      return
    endif
    if(.not.all(ieee_is_finite(eos%nb)).or..not.all(ieee_is_finite(eos%q)).or. &
       .not.all(ieee_is_finite(eos%q_extra))) then
      status='phase EOS contains nonfinite source values'
      return
    endif
    if(eos%composition_file_present.or.eos%composition_available) then
      status='v1 Maxwell construction is thermo-only; composition-bearing inputs are rejected, not discarded'
      return
    endif
    if(.not.eos%beta_equilibrated) then
      status='Maxwell construction requires a cold beta-equilibrated branch (Q5 approximately zero)'
      return
    endif
    if(abs(eos%temperature_mev)>1.e-12_wp) then
      status='Maxwell construction requires T=0'
      return
    endif
    if(maxval(abs(eos%q(2,:)))>entropy_tolerance) then
      status='Maxwell construction requires zero entropy per baryon Q2'
      return
    endif
    if(maxval(abs(eos%q(5,:)))>lepton_potential_tolerance) then
      status='Maxwell construction requires beta-equilibrium lepton potential Q5 approximately zero'
      return
    endif
    if(maxval(abs(eos%q(6,:)-eos%q(7,:)))>algebraic_tolerance) then
      status='Maxwell construction requires cold equality Q6=Q7'
      return
    endif
    if(any(eos%iphase(2:)/=eos%iphase(:npts-1))) then
      status='each input must be a single smooth phase; an input Iphase boundary was found'
      return
    endif
    if(eos%nadd/=0) then
      status='v1 Maxwell construction rejects provider-specific Nadd because its phase semantics are not defined generically'
      return
    endif
    if(eos%mn<=0._wp.or.eos%mp<=0._wp.or..not.ieee_is_finite(eos%mn).or. &
       .not.ieee_is_finite(eos%mp).or..not.ieee_is_finite(eos%temperature_mev).or. &
       .not.ieee_is_finite(eos%charge_fraction_grid)) then
      status='invalid nucleon masses or thermodynamic metadata in phase EOS'
      return
    endif
    if(eos%il/=0.and.eos%il/=1) then
      status='unsupported CompOSE lepton flag in phase EOS'
      return
    endif

    allocate(mu_grid(npts),p_grid(npts),eps_grid(npts))
    p_grid=eos%nb*eos%q(1,:)
    eps_grid=eos_energy(eos)
    if(.not.all(ieee_is_finite(p_grid)).or..not.all(ieee_is_finite(eps_grid))) then
      status='phase pressure or energy density overflowed or became nonfinite'
      return
    endif
    if(any(eos%nb<=0._wp).or.any(p_grid<=0._wp).or.any(eps_grid<=0._wp)) then
      status='Maxwell construction requires positive n_B, pressure and energy density'
      return
    endif
    if(any(eos%nb(2:)<=eos%nb(:npts-1)).or.any(p_grid(2:)<=p_grid(:npts-1)).or. &
       any(eps_grid(2:)<=eps_grid(:npts-1))) then
      status='each phase must have strictly increasing n_B, P and epsilon'
      return
    endif

    mu_grid=(eps_grid+p_grid)/eos%nb
    if(any(.not.ieee_is_finite(mu_grid)).or.any(mu_grid<=0._wp)) then
      status='Euler-derived baryon chemical potential must be finite and positive'
      return
    endif
    if(any(mu_grid(2:)<=mu_grid(:npts-1))) then
      status='each phase must have strictly increasing Euler-derived mu_B'
      return
    endif
    max_euler=0._wp
    do i=1,npts
      mu_source=eos%q(3,i)*eos%mn+eos%mn
      ! CompOSE Il=0 stores the charge-shifted baryon chemical potential:
      ! mu_Q3 = mu_Euler - Yq*mu_q. Recover Euler mu before comparing.
      ! Il=1 with Q5~0 uses mu_Q3=mu_Euler.
      if(eos%il==0)mu_source=mu_source+eos%charge_fraction_grid*eos%q(4,i)*eos%mn
      residual=abs(mu_source-mu_grid(i))/max(abs(mu_grid(i)),1._wp)
      max_euler=max(max_euler,residual)
    enddo
    if(max_euler>euler_tolerance) then
      status='provider Q3 chemical potential disagrees with the Euler value (relative tolerance 5e-4)'
      return
    endif

    ! Within one Hermite interval, n(mu) is quadratic and dn/dmu is linear.
    ! Therefore n>0, dn/dmu>0 and mu*dn/dmu-n>=0 can be checked over the
    ! entire interval by testing endpoint values and the exact quadratic vertex.
    ! The last inequality is precisely c_s^2/c^2 = n/(mu*dn/dmu) <= 1.
    do i=1,npts-1
      mu0=mu_grid(i)
      mu1=mu_grid(i+1)
      mumid=0.5_wp*(mu0+mu1)
      call hermite_pressure(mu0,mu1,p_grid(i),p_grid(i+1),eos%nb(i),eos%nb(i+1), &
        mu0,p0,n0,d0)
      call hermite_pressure(mu0,mu1,p_grid(i),p_grid(i+1),eos%nb(i),eos%nb(i+1), &
        mumid,pmid,nmid,dmid)
      call hermite_pressure(mu0,mu1,p_grid(i),p_grid(i+1),eos%nb(i),eos%nb(i+1), &
        mu1,p1,n1,d1)
      h0=mu0*d0-n0
      hmid=mumid*dmid-nmid
      h1=mu1*d1-n1
      call quadratic_interval_minimum(n0,nmid,n1,min_n)
      call quadratic_interval_minimum(h0,hmid,h1,min_h)
      min_dn=min(d0,d1)
      scale=max(tiny(1._wp),abs(n0),abs(nmid),abs(n1),abs(mu0*d0),abs(mu1*d1))
      if(min_n<=0._wp.or.min_dn<=0._wp) then
        status='P(mu_B) interpolation violates positive density or mechanical stability'
        return
      endif
      if(min_h < -causal_tolerance*scale) then
        status='P(mu_B) interpolation violates the causal bound c_s^2/c^2<=1'
        return
      endif
      if(min(p0,p1)<=0._wp.or.min(eps_grid(i),eps_grid(i+1))<=0._wp) then
        status='P(mu_B) interpolation has nonpositive pressure or energy density'
        return
      endif
      ! The following samples only summarize the sound-speed range; acceptance
      ! above is based on whole-interval polynomial extrema, not these samples.
      do j=0,32
        muq=mu0+(mu1-mu0)*real(j,wp)/32._wp
        call hermite_pressure(mu0,mu1,p_grid(i),p_grid(i+1),eos%nb(i),eos%nb(i+1), &
          muq,p,nbval,dndmu)
        cs2=nbval/(muq*dndmu)
        min_cs2=min(min_cs2,cs2)
        max_cs2=max(max_cs2,cs2)
      enddo
    enddo
    status='ok'
  end subroutine validate_phase

  recursive subroutine quick_sort(values,left,right)
    real(wp),intent(inout) :: values(:)
    integer,intent(in) :: left,right
    integer :: i,j
    real(wp) :: pivot,tmp
    if(left>=right) return
    i=left
    j=right
    pivot=values((left+right)/2)
    do
      do while(values(i)<pivot)
        i=i+1
      enddo
      do while(values(j)>pivot)
        j=j-1
      enddo
      if(i<=j) then
        tmp=values(i)
        values(i)=values(j)
        values(j)=tmp
        i=i+1
        j=j-1
      endif
      if(i>j) exit
    enddo
    if(left<j) call quick_sort(values,left,j)
    if(i<right) call quick_sort(values,i,right)
  end subroutine quick_sort

  subroutine build_union_grid(mu_a,mu_b,mu_min,mu_max,grid)
    real(wp),intent(in) :: mu_a(:),mu_b(:),mu_min,mu_max
    real(wp),allocatable,intent(out) :: grid(:)
    real(wp),allocatable :: raw(:),unique(:)
    integer :: i,k,nraw,nu
    real(wp) :: tol
    nraw=size(mu_a)+size(mu_b)+2
    allocate(raw(nraw))
    k=1
    raw(k)=mu_min
    k=k+1
    raw(k)=mu_max
    k=k+1
    do i=1,size(mu_a)
      if(mu_a(i)>mu_min.and.mu_a(i)<mu_max) then
        raw(k)=mu_a(i)
        k=k+1
      endif
    enddo
    do i=1,size(mu_b)
      if(mu_b(i)>mu_min.and.mu_b(i)<mu_max) then
        raw(k)=mu_b(i)
        k=k+1
      endif
    enddo
    call quick_sort(raw,1,k-1)
    allocate(unique(k-1))
    nu=0
    do i=1,k-1
      tol=1.e-12_wp*max(1._wp,abs(raw(i)))
      if(nu==0) then
        nu=1
        unique(nu)=raw(i)
      else if(abs(raw(i)-unique(nu))>tol) then
        nu=nu+1
        unique(nu)=raw(i)
      endif
    enddo
    if(nu>=2) then
      allocate(grid(nu))
      grid=unique(:nu)
    else
      allocate(grid(0))
    endif
  end subroutine build_union_grid

  subroutine evaluate_difference(low,high,mu_low,p_low,mu_high,p_high,muq,delta_p,delta_n,pressure_scale,ok)
    type(eos_data_t),intent(in) :: low,high
    real(wp),intent(in) :: mu_low(:),p_low(:),mu_high(:),p_high(:),muq
    real(wp),intent(out) :: delta_p,delta_n,pressure_scale
    logical,intent(out) :: ok
    real(wp) :: pa,na,ca,pb,nb,cb,q4,q5,qe
    logical :: oka,okb
    call evaluate_phase(low,mu_low,p_low,muq,pa,na,ca,q4,q5,qe,oka)
    call evaluate_phase(high,mu_high,p_high,muq,pb,nb,cb,q4,q5,qe,okb)
    ok=oka.and.okb
    delta_p=pa-pb
    delta_n=na-nb
    pressure_scale=max(tiny(1._wp),abs(pa),abs(pb))
  end subroutine evaluate_difference

  subroutine add_unique_root(root,roots,nroots)
    real(wp),intent(in) :: root
    real(wp),intent(inout) :: roots(:)
    integer,intent(inout) :: nroots
    integer :: i
    real(wp) :: tol
    tol=1.e-10_wp*max(1._wp,abs(root))
    do i=1,nroots
      if(abs(roots(i)-root)<=tol) return
    enddo
    if(nroots<size(roots)) then
      nroots=nroots+1
      roots(nroots)=root
    endif
  end subroutine add_unique_root

  subroutine bisect_pressure_root(low,high,mu_low,p_low,mu_high,p_high,a0,b0,fa0,root,ok)
    type(eos_data_t),intent(in) :: low,high
    real(wp),intent(in) :: mu_low(:),p_low(:),mu_high(:),p_high(:)
    real(wp),intent(in) :: a0,b0,fa0
    real(wp),intent(out) :: root
    logical,intent(out) :: ok
    real(wp) :: a,b,fa,mid,fmid,dn,scale
    integer :: iter
    logical :: valid
    a=a0
    b=b0
    fa=fa0
    root=0._wp
    ok=.false.
    do iter=1,120
      mid=0.5_wp*(a+b)
      call evaluate_difference(low,high,mu_low,p_low,mu_high,p_high,mid,fmid,dn,scale,valid)
      if(.not.valid) return
      if((b-a)<=1.e-13_wp*max(1._wp,abs(mid))) then
        root=mid
        ok=abs(fmid)<=10._wp*pressure_relative_tolerance*scale
        return
      endif
      if((fa<0._wp.and.fmid>0._wp).or.(fa>0._wp.and.fmid<0._wp)) then
        b=mid
      else
        a=mid
        fa=fmid
      endif
    enddo
    root=0.5_wp*(a+b)
    call evaluate_difference(low,high,mu_low,p_low,mu_high,p_high,root,fmid,dn,scale,valid)
    ok=valid.and.abs(fmid)<=10._wp*pressure_relative_tolerance*scale
  end subroutine bisect_pressure_root

  subroutine find_coexistence(low,high,mu_low,p_low,mu_high,p_high,mu_t,p_t,root_count,status)
    type(eos_data_t),intent(in) :: low,high
    real(wp),intent(in) :: mu_low(:),p_low(:),mu_high(:),p_high(:)
    real(wp),intent(out) :: mu_t,p_t
    integer,intent(out) :: root_count
    character(len=*),intent(out) :: status
    real(wp),allocatable :: grid(:),roots(:)
    real(wp) :: mu_min,mu_max,mu0,mu1,mumid,dp0,dpmid,dp1,dn
    real(wp) :: g0,gm,g1,qa,qb,qc,disc,r1,r2
    real(wp) :: parts(4),tmp,fa,fb,scalea,scaleb,scale0,scalemid,scale1,root
    real(wp) :: p,n,dndmu,q4,q5,qe
    integer :: i,j,nparts,ncritical,nroots,k
    logical :: ok,root_ok

    status='uninitialized'
    mu_t=0._wp
    p_t=0._wp
    root_count=0
    mu_min=max(mu_low(1),mu_high(1))
    mu_max=min(mu_low(size(mu_low)),mu_high(size(mu_high)))
    if(mu_max<=mu_min) then
      status='the two phases have no overlapping baryon-chemical-potential interval'
      return
    endif
    call build_union_grid(mu_low,mu_high,mu_min,mu_max,grid)
    if(size(grid)<2) then
      status='insufficient common chemical-potential support for coexistence search'
      return
    endif
    allocate(roots(3*(size(grid)-1)))
    nroots=0

    ! On each union interval both Hermite segments are fixed.  Their pressure
    ! difference is cubic; its derivative n_low-n_high is exactly quadratic.
    ! Solving that quadratic isolates all stationary points, making the
    ! pressure-difference root search monotone on each subinterval.
    do i=1,size(grid)-1
      mu0=grid(i)
      mu1=grid(i+1)
      mumid=0.5_wp*(mu0+mu1)
      call evaluate_difference(low,high,mu_low,p_low,mu_high,p_high,mu0,dp0,g0,scale0,ok)
      if(.not.ok) then
        status='failed to evaluate coexistence pressure at a union-grid endpoint'
        return
      endif
      call evaluate_difference(low,high,mu_low,p_low,mu_high,p_high,mumid,dpmid,gm,scalemid,ok)
      if(.not.ok) then
        status='failed to evaluate coexistence pressure at a union-grid midpoint'
        return
      endif
      call evaluate_difference(low,high,mu_low,p_low,mu_high,p_high,mu1,dp1,g1,scale1,ok)
      if(.not.ok) then
        status='failed to evaluate coexistence pressure at a union-grid endpoint'
        return
      endif

      ! Fit g(u)=n_low-n_high=a*u^2+b*u+c from u=0, 1/2, 1.
      qc=g0
      qa=2._wp*(g1+g0-2._wp*gm)
      qb=g1-g0-qa
      parts=0._wp
      parts(1)=0._wp
      nparts=1
      ncritical=0
      if(abs(qa)<=1.e-13_wp*max(1._wp,abs(qb),abs(qc))) then
        if(abs(qb)>1.e-13_wp*max(1._wp,abs(qc))) then
          r1=-qc/qb
          if(r1>1.e-12_wp.and.r1<1._wp-1.e-12_wp) then
            nparts=nparts+1
            parts(nparts)=r1
            ncritical=1
          endif
        endif
      else
        disc=qb*qb-4._wp*qa*qc
        if(disc>=-1.e-13_wp*max(1._wp,qb*qb,abs(4._wp*qa*qc))) then
          disc=max(0._wp,disc)
          ! Stable quadratic formula: avoid cancellation when |qb|~sqrt(disc).
          if(qb>=0._wp) then
            tmp=-0.5_wp*(qb+sqrt(disc))
          else
            tmp=-0.5_wp*(qb-sqrt(disc))
          endif
          if(abs(tmp)>tiny(1._wp)) then
            r1=tmp/qa
            r2=qc/tmp
          else
            r1=(-qb-sqrt(disc))/(2._wp*qa)
            r2=(-qb+sqrt(disc))/(2._wp*qa)
          endif
          if(r1>1.e-12_wp.and.r1<1._wp-1.e-12_wp) then
            nparts=nparts+1
            parts(nparts)=r1
            ncritical=ncritical+1
          endif
          if(r2>1.e-12_wp.and.r2<1._wp-1.e-12_wp.and.abs(r2-r1)>1.e-10_wp) then
            nparts=nparts+1
            parts(nparts)=r2
            ncritical=ncritical+1
          endif
        endif
      endif
      nparts=nparts+1
      parts(nparts)=1._wp
      do j=2,nparts
        tmp=parts(j)
        k=j-1
        do while(k>=1)
          if(parts(k)<=tmp) exit
          parts(k+1)=parts(k)
          k=k-1
        enddo
        parts(k+1)=tmp
      enddo

      do j=1,nparts-1
        mu0=grid(i)+(grid(i+1)-grid(i))*parts(j)
        mu1=grid(i)+(grid(i+1)-grid(i))*parts(j+1)
        call evaluate_difference(low,high,mu_low,p_low,mu_high,p_high,mu0,fa,dn,scalea,ok)
        if(.not.ok) then
          status='failed to evaluate pressure difference on a monotone subinterval'
          return
        endif
        call evaluate_difference(low,high,mu_low,p_low,mu_high,p_high,mu1,fb,dn,scaleb,ok)
        if(.not.ok) then
          status='failed to evaluate pressure difference on a monotone subinterval'
          return
        endif
        if(abs(fa)<=pressure_relative_tolerance*scalea) call add_unique_root(mu0,roots,nroots)
        if(abs(fb)<=pressure_relative_tolerance*scaleb) call add_unique_root(mu1,roots,nroots)
        if(fa*fb<0._wp) then
          call bisect_pressure_root(low,high,mu_low,p_low,mu_high,p_high, &
            mu0,mu1,fa,root,root_ok)
          if(.not.root_ok) then
            status='pressure equality root did not converge to the requested tolerance'
            return
          endif
          call add_unique_root(root,roots,nroots)
        endif
      enddo
    enddo

    root_count=nroots
    if(nroots==0) then
      status='no pressure equality was found in the common mu_B interval'
      return
    endif
    if(nroots/=1) then
      status='multiple pressure-equality roots found; v1 accepts exactly one isolated transition'
      return
    endif
    mu_t=roots(1)
    call evaluate_phase(low,mu_low,p_low,mu_t,p,n,dndmu,q4,q5,qe,ok)
    if(.not.ok) then
      status='failed to evaluate low-density phase at the coexistence point'
      return
    endif
    p_t=p
    call evaluate_phase(high,mu_high,p_high,mu_t,p,n,dndmu,q4,q5,qe,ok)
    if(.not.ok) then
      status='failed to evaluate high-density phase at the coexistence point'
      return
    endif
    if(abs(p-p_t)>10._wp*pressure_relative_tolerance*max(tiny(1._wp),abs(p),abs(p_t))) then
      status='coexistence pressure equality failed its final residual check'
      return
    endif
    status='ok'
  end subroutine find_coexistence

  subroutine set_transition_row(out,k,eos,mu_grid,p_grid,mu_t,pressure,n_transition,phase_id)
    type(eos_data_t),intent(inout) :: out
    integer,intent(in) :: k,phase_id
    type(eos_data_t),intent(in) :: eos
    real(wp),intent(in) :: mu_grid(:),p_grid(:),mu_t,pressure,n_transition
    real(wp) :: p,n,dndmu,q4,q5,qextra,eps
    logical :: ok
    call evaluate_phase(eos,mu_grid,p_grid,mu_t,p,n,dndmu,q4,q5,qextra,ok)
    out%nb(k)=n_transition
    out%q(1,k)=pressure/n_transition
    out%q(2,k)=0._wp
    out%q(3,k)=mu_t/out%mn-1._wp
    out%q(4,k)=q4
    if(eos%il==0)out%q(3,k)=out%q(3,k)-eos%charge_fraction_grid*out%q(4,k)
    out%q(5,k)=0._wp
    eps=mu_t*n_transition-pressure
    out%q(6,k)=eps/(n_transition*out%mn)-1._wp
    out%q(7,k)=out%q(6,k)
    out%q_extra(k)=0._wp
    out%iphase(k)=phase_id
  end subroutine set_transition_row

  subroutine copy_source_row(out,k,eos,i,mu_euler,phase_id)
    type(eos_data_t),intent(inout) :: out
    integer,intent(in) :: k,i,phase_id
    type(eos_data_t),intent(in) :: eos
    real(wp),intent(in) :: mu_euler
    out%nb(k)=eos%nb(i)
    out%q(:,k)=eos%q(:,i)
    out%q(2,k)=0._wp
    out%q(5,k)=0._wp
    out%q(3,k)=mu_euler/out%mn-1._wp
    if(eos%il==0)out%q(3,k)=out%q(3,k)-eos%charge_fraction_grid*out%q(4,k)
    out%q(6,k)=eos%q(7,i)
    out%q(7,k)=eos%q(7,i)
    out%q_extra(k)=0._wp
    out%iphase(k)=phase_id
  end subroutine copy_source_row

  subroutine construct_maxwell_transition(low,high,local_neutrality_confirmed,out,report,status)
    type(eos_data_t),intent(in) :: low,high
    logical,intent(in) :: local_neutrality_confirmed
    type(eos_data_t),intent(out) :: out
    type(maxwell_transition_report_t),intent(out) :: report
    character(len=*),intent(out) :: status
    real(wp),allocatable :: mu_low(:),p_low(:),eps_low(:),mu_high(:),p_high(:),eps_high(:)
    real(wp) :: euler_low,euler_high,mincs_low,maxcs_low,mincs_high,maxcs_high
    real(wp) :: mu_t,p_t,pa,na,ca,pb,nb,cb,q4,q5,qe,epsa,epsb
    real(wp) :: identity_residual,mu_tol
    real(wp),allocatable :: p_out(:),eps_out(:)
    integer :: roots,nlow,nhigh,i,k
    logical :: ok

    report=maxwell_transition_report_t()
    status='uninitialized'
    if(.not.local_neutrality_confirmed) then
      status='caller must explicitly confirm that both input branches are locally charge neutral'
      return
    endif
    if(abs(low%mn-high%mn)>1.e-12_wp*max(low%mn,high%mn).or. &
       abs(low%mp-high%mp)>1.e-12_wp*max(low%mp,high%mp)) then
      status='phase tables must use the same neutron and proton rest-mass conventions'
      return
    endif
    if(low%il/=high%il.or.low%iyq/=high%iyq.or. &
       abs(low%charge_fraction_grid-high%charge_fraction_grid)>1.e-12_wp) then
      status='phase tables must use the same CompOSE lepton-content and charge-fraction conventions'
      return
    endif

    call validate_phase(low,mu_low,p_low,eps_low,euler_low,mincs_low,maxcs_low,status)
    if(status/='ok') then
      status='low-density phase: '//trim(status)
      return
    endif
    call validate_phase(high,mu_high,p_high,eps_high,euler_high,mincs_high,maxcs_high,status)
    if(status/='ok') then
      status='high-density phase: '//trim(status)
      return
    endif

    call find_coexistence(low,high,mu_low,p_low,mu_high,p_high,mu_t,p_t,roots,status)
    report%pressure_crossings=roots
    if(status/='ok') return

    call evaluate_phase(low,mu_low,p_low,mu_t,pa,na,ca,q4,q5,qe,ok)
    if(.not.ok) then
      status='failed to evaluate low-density phase at coexistence'
      return
    endif
    call evaluate_phase(high,mu_high,p_high,mu_t,pb,nb,cb,q4,q5,qe,ok)
    if(.not.ok) then
      status='failed to evaluate high-density phase at coexistence'
      return
    endif
    if(nb<=na+1.e-8_wp*max(na,nb)) then
      status='pressure equality does not connect a distinct denser high-density phase'
      return
    endif
    epsa=mu_t*na-p_t
    epsb=mu_t*nb-p_t
    if(epsa<=0._wp.or.epsb<=epsa) then
      status='coexistence energy densities are nonpositive or not ordered'
      return
    endif

    mu_tol=1.e-10_wp*max(1._wp,abs(mu_t))
    nlow=count(mu_low<mu_t-mu_tol)
    nhigh=count(mu_high>mu_t+mu_tol)
    if(nlow<2.or.nhigh<2) then
      status='fewer than two retained source knots remain on one side of the coexistence point'
      return
    endif
    call allocate_eos(out,nlow+2+nhigh)
    out%mn=low%mn
    out%mp=low%mp
    out%il=low%il
    out%iT=low%iT
    out%iYq=low%iYq
    out%inb_start=1
    out%nadd=0
    out%temperature_mev=0._wp
    out%charge_fraction_grid=low%charge_fraction_grid
    out%composition_file_present=.false.
    out%composition_available=.false.
    out%beta_equilibrated=.true.

    k=0
    do i=1,low%n
      if(mu_low(i)<mu_t-mu_tol) then
        k=k+1
        call copy_source_row(out,k,low,i,mu_low(i),1)
      endif
    enddo
    k=k+1
    call set_transition_row(out,k,low,mu_low,p_low,mu_t,p_t,na,1)
    k=k+1
    call set_transition_row(out,k,high,mu_high,p_high,mu_t,p_t,nb,2)
    do i=1,high%n
      if(mu_high(i)>mu_t+mu_tol) then
        k=k+1
        call copy_source_row(out,k,high,i,mu_high(i),2)
      endif
    enddo
    if(k/=out%n) then
      status='internal Maxwell output row count mismatch'
      return
    endif

    out%yn=0._wp
    out%yp=0._wp
    out%ye=0._wp
    out%ymu=0._wp
    out%abar=0._wp
    out%zbar=0._wp
    out%yh=0._wp

    if(any(.not.ieee_is_finite(out%nb)).or.any(.not.ieee_is_finite(out%q))) then
      status='Maxwell output contains nonfinite thermodynamic rows'
      return
    endif
    allocate(p_out(out%n),eps_out(out%n))
    p_out=eos_pressure(out)
    eps_out=eos_energy(out)
    if(any(out%nb(2:)<=out%nb(:out%n-1)).or.any(p_out(2:)<p_out(:out%n-1))) then
      status='Maxwell output is not ordered by strictly increasing density and nondecreasing pressure'
      return
    endif
    if(any(eps_out(2:)<=eps_out(:out%n-1))) then
      status='Maxwell output energy density is not strictly increasing'
      return
    endif

    identity_residual=abs((epsb-epsa)-mu_t*(nb-na))/max(abs(epsb-epsa),1._wp)
    report%pressure_crossings=roots
    report%low_phase_rows=nlow+1
    report%high_phase_rows=nhigh+1
    report%mu_transition=mu_t
    report%pressure_transition=p_t
    report%n_low=na
    report%n_high=nb
    report%epsilon_low=epsa
    report%epsilon_high=epsb
    report%density_jump=nb-na
    report%latent_energy_density=epsb-epsa
    report%maxwell_identity_residual=identity_residual
    report%max_euler_residual_low=euler_low
    report%max_euler_residual_high=euler_high
    report%min_cs2_low=mincs_low
    report%max_cs2_low=maxcs_low
    report%min_cs2_high=mincs_high
    report%max_cs2_high=maxcs_high
    if(identity_residual>1.e-10_wp) then
      status='Maxwell latent-energy identity failed its numerical residual check'
      return
    endif
    status='ok'
  end subroutine construct_maxwell_transition

end module compose_eos_phase_transition
