module compose_eos_bbp_original
  !! BBP (1971), Nucl. Phys. A175, 225--271, DOI 10.1016/0375-9474(71)90281-8.
  !! Original bulk, surface and finite-thickness functional, spherical cells.
  !! Published closure omits ONLY the W_g derivative in w_s (p.251).
  !! Variational closure retains it. Neither path integrates a pressure table.
  use compose_eos_types, only: wp,eos_data_t,allocate_eos,validate_eos
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: generate_bbp_original_published_native,generate_bbp_original_variational_native
  public :: generate_bbp_variational_maxwell_native,bbp_original_maxwell
  public :: bbp_original_bulk,bbp_original_point,bbp_original_solve,bbp_original_state_t
  real(wp),parameter :: pi=acos(-1._wp),hbarc=197.3269804_wp,e2=1.43996454784255_wp
  real(wp),parameter :: mn=939.56542052_wp,mp=938.27208816_wp
  real(wp),parameter :: w0=16.5_wp,k0=1.43_wp,compress=143._wp,symmetry=33._wp,sigma=21._wp
  real(wp),parameter :: tncoef=19.74_wp,nb_lo=2.67e-4_wp,nb_hi=0.08_wp
  integer,parameter :: nd=4,nparams=21
  type :: dual_t
    real(wp) :: v=0._wp,d(nd)=0._wp
  end type
  type :: bbp_original_state_t
    real(wp) :: nb=0,k=0,x=0,kg=0,gas_fraction=0,u=0,a=0,z=0
    real(wp) :: energy=0,pressure=0,mu=0,mue=0,ws=0,wc=0,wthick=0,binding=0
    real(wp) :: gradient(3)=0,true_gradient(3)=0,hessian_min=0,uniform_delta=0
    real(wp) :: firstlaw_local=0
    integer :: iterations=0
    logical :: gas_free=.false.,converged=.false.
  end type
  interface operator(+)
    module procedure add_dd,add_dr,add_rd
  end interface
  interface operator(-)
    module procedure sub_dd,sub_dr,sub_rd,neg_d
  end interface
  interface operator(*)
    module procedure mul_dd,mul_dr,mul_rd
  end interface
  interface operator(/)
    module procedure div_dd,div_dr,div_rd
  end interface
  interface operator(**)
    module procedure pow_dr,pow_di
  end interface
contains
  pure function seed(v,i) result(r)
    real(wp),intent(in)::v
    integer,intent(in)::i
    type(dual_t)::r
    r%v=v;r%d=0
    if(i>0)r%d(i)=1
  end function
  pure function add_dd(a,b) result(r)
    type(dual_t),intent(in)::a,b
    type(dual_t)::r
    r%v=a%v+b%v;r%d=a%d+b%d
  end function
  pure function add_dr(a,b) result(r)
    type(dual_t),intent(in)::a
    real(wp),intent(in)::b
    type(dual_t)::r
    r=a;r%v=a%v+b
  end function
  pure function add_rd(a,b) result(r)
    real(wp),intent(in)::a
    type(dual_t),intent(in)::b
    type(dual_t)::r
    r=b+a
  end function
  pure function sub_dd(a,b) result(r)
    type(dual_t),intent(in)::a,b
    type(dual_t)::r
    r%v=a%v-b%v;r%d=a%d-b%d
  end function
  pure function sub_dr(a,b) result(r)
    type(dual_t),intent(in)::a
    real(wp),intent(in)::b
    type(dual_t)::r
    r=a;r%v=a%v-b
  end function
  pure function sub_rd(a,b) result(r)
    real(wp),intent(in)::a
    type(dual_t),intent(in)::b
    type(dual_t)::r
    r%v=a-b%v;r%d=-b%d
  end function
  pure function neg_d(a) result(r)
    type(dual_t),intent(in)::a
    type(dual_t)::r
    r%v=-a%v;r%d=-a%d
  end function
  pure function mul_dd(a,b) result(r)
    type(dual_t),intent(in)::a,b
    type(dual_t)::r
    r%v=a%v*b%v;r%d=a%d*b%v+a%v*b%d
  end function
  pure function mul_dr(a,b) result(r)
    type(dual_t),intent(in)::a
    real(wp),intent(in)::b
    type(dual_t)::r
    r%v=a%v*b;r%d=a%d*b
  end function
  pure function mul_rd(a,b) result(r)
    real(wp),intent(in)::a
    type(dual_t),intent(in)::b
    type(dual_t)::r
    r=b*a
  end function
  pure function div_dd(a,b) result(r)
    type(dual_t),intent(in)::a,b
    type(dual_t)::r
    r%v=a%v/b%v;r%d=(a%d-r%v*b%d)/b%v
  end function
  pure function div_dr(a,b) result(r)
    type(dual_t),intent(in)::a
    real(wp),intent(in)::b
    type(dual_t)::r
    r%v=a%v/b;r%d=a%d/b
  end function
  pure function div_rd(a,b) result(r)
    real(wp),intent(in)::a
    type(dual_t),intent(in)::b
    type(dual_t)::r
    r%v=a/b%v;r%d=-r%v*b%d/b%v
  end function
  pure function pow_dr(a,b) result(r)
    type(dual_t),intent(in)::a
    real(wp),intent(in)::b
    type(dual_t)::r
    r%v=a%v**b;r%d=b*a%v**(b-1._wp)*a%d
  end function
  pure function pow_di(a,b) result(r)
    type(dual_t),intent(in)::a
    integer,intent(in)::b
    type(dual_t)::r
    r%v=a%v**b
    if(b==0)then
      r%d=0
    else
      r%d=real(b,wp)*a%v**(b-1)*a%d
    endif
  end function

  pure function neutron_w(k) result(w)
    type(dual_t),intent(in)::k
    type(dual_t)::w
    ! Eq.3.8: denominator is 1+2.545*k, not k**2.
    w=tncoef*k**2-k**3*(40.4_wp-1.088_wp*k**3)/(1._wp+2.545_wp*k)
  end function
  pure function neutron_dw(k) result(w)
    type(dual_t),intent(in)::k
    type(dual_t)::w,num,dnum,den
    num=40.4_wp*k**3-1.088_wp*k**6
    dnum=121.2_wp*k**2-6.528_wp*k**5
    den=1._wp+2.545_wp*k
    w=2._wp*tncoef*k-(dnum*den-2.545_wp*num)/den**2
  end function
  pure function bulk_w(k,x) result(w)
    type(dual_t),intent(in)::k,x
    type(dual_t)::w,y,ts,tn,ws,wn,alpha,mp0,mn0,kin
    y=k/k0;tn=tncoef*k**2;ts=tn/2._wp**(2._wp/3._wp)
    ws=ts*(1._wp-y)**3-w0*y**3*(1._wp+(1._wp-y)*(9._wp-6._wp*y)) &
       +0.5_wp*compress*(1._wp-y)**2*y**3
    wn=neutron_w(k);mn0=wn+k*neutron_dw(k)/3._wp
    mp0=-k**3*(218._wp+277._wp*k)/(1._wp+8.57_wp*k**2)
    alpha=1._wp-2._wp*x
    kin=tn*(x**(5._wp/3._wp)+(1._wp-x)**(5._wp/3._wp))
    w=(ws-ts)*(1._wp-3._wp*alpha**4+2._wp*alpha**6) &
      +(symmetry*y**2-5._wp*ts/9._wp)*alpha**2*(1._wp-alpha**2)**2 &
      +(wn-tn)*(3._wp*alpha**4-2._wp*alpha**6) &
      +0.25_wp*(mp0-mn0+5._wp*tn/3._wp)*(alpha**4-alpha**6)+kin
  end function

  subroutine bbp_original_bulk(k,x,w,dw_dk,dw_dx,wn,mun0,mup0,ok)
    real(wp),intent(in)::k,x
    real(wp),intent(out)::w,dw_dk,dw_dx,wn,mun0,mup0
    logical,intent(out)::ok
    type(dual_t)::kd,xd,v,vn,dn
    w=0;dw_dk=0;dw_dx=0;wn=0;mun0=0;mup0=0;ok=.false.
    if(.not.ieee_is_finite(k).or..not.ieee_is_finite(x))return
    if(k<=0.or.k>1.5_wp.or.x<0.or.x>1) return
    kd=seed(k,1);xd=seed(x,2)
    ! The fractional kinetic power has a regular first derivative at x=0,1.
    v=bulk_w(kd,xd);vn=neutron_w(kd);dn=neutron_dw(kd)
    w=v%v;dw_dk=v%d(1);dw_dx=v%d(2);wn=vn%v
    mun0=vn%v+k*dn%v/3._wp
    mup0=-k**3*(218._wp+277._wp*k)/(1._wp+8.57_wp*k**2)
    ok=all(ieee_is_finite([w,dw_dk,dw_dx,wn,mun0,mup0]))
  end subroutine

  subroutine bbp_original_point(nb,z,published,state,f,g,ok)
    real(wp),intent(in)::nb,z(3)
    logical,intent(in)::published
    type(bbp_original_state_t),intent(out)::state
    real(wp),intent(out)::f,g(4)
    logical,intent(out)::ok
    type(dual_t)::k,x,t,b,ni,ng,kg,u,wi,wg,wgs,ws,wc,wth,a,ne,ke,ee,enet,ff,d,kc
    state=bbp_original_state_t();f=huge(1._wp);g=0;ok=.false.
    if(.not.ieee_is_finite(nb).or..not.all(ieee_is_finite(z)))return
    if(nb<nb_lo.or.nb>nb_hi.or.z(1)<=0.or.z(1)>1.5_wp.or.z(2)<=0.or.z(2)>=1.or.z(3)<0.or.z(3)>=1) return
    if(z(3)>0.and.z(3)<1.e-9_wp)return
    k=seed(z(1),1);x=seed(z(2),2);t=seed(z(3),3);b=seed(nb,4)
    ni=k**3/(1.5_wp*pi**2);ng=b*t
    if(ni%v<=nb.or.ni%v<=ng%v)return
    u=(b-ng)/(ni-ng)
    if(u%v<=0.or.u%v>=1)return
    wi=bulk_w(k,x)
    if(z(3)==0)then
      kg=seed(0._wp,0);wg=seed(0._wp,0)
      ! The exact gas-free state is evaluated on its two-dimensional face;
      ! its singular one-sided surface derivative is not replaced by zero.
      ng=seed(0._wp,0);u=b/ni
    else
      kg=(1.5_wp*pi**2*ng)**(1._wp/3._wp);wg=neutron_w(kg)
    endif
    wgs=wg
    if(published)wgs%d=0._wp ! BBP p.251: omit dW_g/dn_g ONLY in surface term.
    ws=(sigma/w0)*(wgs-wi)*(1._wp-ng/ni)**(2._wp/3._wp)
    wc=(3._wp/5._wp)*(8._wp/(9._wp*pi))**(1._wp/3._wp)*e2*k &
       *(1._wp-1.5_wp*u**(1._wp/3._wp)+0.5_wp*u)
    if(ws%v<=0.or.wc%v<=0)return
    a=ws/(2._wp*wc*x**2) ! Eq.6.4; exact minimization in A at fixed u.
    kc=(k**3-kg**3)**(1._wp/3._wp);d=0.74_wp/kc
    wth=-(4._wp*pi/9._wp)*e2*d**2*k**3*x**2 ! Eq.5.3 divided by A.
    ne=u*ni*x;ke=(3._wp*pi**2*ne)**(1._wp/3._wp)
    ee=0.75_wp*hbarc*ke*ne ! Eq.5.13, ultrarelativistic electrons.
    ! Exchange and screening omitted as authorized on pp.246--248.
    enet=u*ni*(wi+wth+(mp-mn)*x+1.5_wp*ws/a**(1._wp/3._wp)) &
         +(1._wp-u)*ng*wg+ee
    ff=enet/b;f=ff%v;g=ff%d
    state%nb=nb;state%k=k%v;state%x=x%v;state%kg=kg%v
    state%gas_fraction=z(3);state%u=u%v;state%a=a%v;state%z=a%v*x%v
    state%energy=nb*(mn+f);state%pressure=nb**2*g(4);state%binding=f
    state%mu=mn+f+nb*g(4);state%mue=hbarc*ke%v
    state%ws=ws%v;state%wc=wc%v;state%wthick=wth%v
    state%gradient=g(:3);state%gas_free=z(3)==0
    ok=all(ieee_is_finite([f,g,state%a,state%energy,state%pressure,state%mu]))
  end subroutine

  subroutine bounds(nb,lower,upper)
    real(wp),intent(in)::nb
    real(wp),intent(out)::lower(3),upper(3)
    lower=[max(1.05_wp,(1.5_wp*pi**2*nb*1.000001_wp)**(1._wp/3._wp)),0.025_wp,1.e-9_wp]
    upper=[1.5_wp,0.45_wp,0.999999_wp]
  end subroutine

  subroutine gradient_jacobian(nb,z,published,active,jac,ok)
    real(wp),intent(in)::nb,z(3)
    logical,intent(in)::published
    integer,intent(in)::active
    real(wp),intent(out)::jac(3,3)
    logical,intent(out)::ok
    real(wp)::zz(3),gp(4),gm(4),f,h,lo(3),hi(3)
    type(bbp_original_state_t)::st
    logical::op,om
    integer::j
    call bounds(nb,lo,hi);jac=0;ok=.true.
    do j=1,active
      h=2.e-5_wp*max(abs(z(j)),0.02_wp)
      h=min(h,0.2_wp*min(z(j)-lo(j),hi(j)-z(j)))
      if(h<=1.e-13_wp)then
        ok=.false.;return
      endif
      zz=z;zz(j)=z(j)+h;call bbp_original_point(nb,zz,published,st,f,gp,op)
      zz=z;zz(j)=z(j)-h;call bbp_original_point(nb,zz,published,st,f,gm,om)
      if(.not.op.or..not.om)then
        ok=.false.;return
      endif
      jac(:active,j)=(gp(:active)-gm(:active))/(2*h)
    enddo
  end subroutine

  subroutine solve_linear(a,b,n,x,ok)
    real(wp),intent(in)::a(3,3),b(3)
    integer,intent(in)::n
    real(wp),intent(out)::x(3)
    logical,intent(out)::ok
    real(wp)::m(3,4),r(4),v
    integer::i,j,p
    m=0;m(:n,:n)=a(:n,:n);m(:n,n+1)=b(:n);x=0;ok=.false.
    do i=1,n
      p=i-1+maxloc(abs(m(i:n,i)),dim=1)
      if(abs(m(p,i))<1.e-14_wp)return
      if(p/=i)then
        r=m(i,:);m(i,:)=m(p,:);m(p,:)=r
      endif
      v=m(i,i);m(i,i:n+1)=m(i,i:n+1)/v
      do j=1,n
        if(j/=i)then
          v=m(j,i);m(j,i:n+1)=m(j,i:n+1)-v*m(i,i:n+1)
        endif
      enddo
    enddo
    x(:n)=m(:n,n+1);ok=all(ieee_is_finite(x))
  end subroutine

  subroutine polish(nb,z,published,active,iterations,ok)
    real(wp),intent(in)::nb
    real(wp),intent(inout)::z(3)
    logical,intent(in)::published
    integer,intent(in)::active
    integer,intent(out)::iterations
    logical,intent(out)::ok
    real(wp)::f,g(4),gt(4),jac(3,3),step(3),rhs(3),zz(3),lo(3),hi(3),fac,norm,trial,alpha
    type(bbp_original_state_t)::st
    logical::valid
    integer::it,ls
    call bounds(nb,lo,hi);ok=.false.;iterations=0
    do it=1,70
      iterations=it
      call bbp_original_point(nb,z,published,st,f,g,valid)
      if(.not.valid)return
      norm=maxval(abs(g(:active)))
      if(norm<2.e-9_wp)then
        ok=.true.;return
      endif
      call gradient_jacobian(nb,z,published,active,jac,valid)
      if(.not.valid)return
      rhs=-g(:3);call solve_linear(jac,rhs,active,step,valid)
      if(.not.valid)return
      fac=max(1._wp,maxval(abs(step)/[0.10_wp,0.05_wp,0.2_wp]));step=step/fac
      alpha=1
      do ls=1,30
        zz=z+alpha*step
        if(all(zz(:active)>lo(:active)).and.all(zz(:active)<hi(:active)))then
          call bbp_original_point(nb,zz,published,st,f,gt,valid)
          if(valid)then
            trial=maxval(abs(gt(:active)))
            if(trial<norm*(1._wp-1.e-4_wp*alpha))exit
          endif
        endif
        alpha=alpha*0.5_wp
      enddo
      if(ls>30)return
      z=zz
    enddo
  end subroutine

  pure function minimum_eigenvalue(a,n) result(ev)
    real(wp),intent(in)::a(3,3)
    integer,intent(in)::n
    real(wp)::ev,m(3,3),c,s,t,tau,ap,aq,b
    integer::i,p,q,j,it
    m=0.5_wp*(a+transpose(a))
    do it=1,30
      b=0;p=1;q=2
      do i=1,n
        do j=i+1,n
          if(abs(m(i,j))>b)then
            b=abs(m(i,j));p=i;q=j
          endif
        enddo
      enddo
      if(b<1.e-11_wp)exit
      tau=(m(q,q)-m(p,p))/(2*m(p,q));t=sign(1._wp,tau)/(abs(tau)+sqrt(1+tau*tau))
      c=1/sqrt(1+t*t);s=t*c
      ap=m(p,p);aq=m(q,q);b=m(p,q)
      m(p,p)=ap-t*b;m(q,q)=aq+t*b;m(p,q)=0;m(q,p)=0
      do j=1,n
        if(j==p.or.j==q)cycle
        ap=m(j,p);aq=m(j,q)
        m(j,p)=c*ap-s*aq;m(p,j)=m(j,p)
        m(j,q)=s*ap+c*aq;m(q,j)=m(j,q)
      enddo
    enddo
    ev=m(1,1)
    do i=2,n
      ev=min(ev,m(i,i))
    enddo
  end function

  function uniform_energy_per_baryon(nb) result(f)
    real(wp),intent(in)::nb
    real(wp)::f,k,l,r,x,w,dwk,dwx,wn,mug,mup,der
    logical::ok
    integer::i
    k=(1.5_wp*pi**2*nb)**(1._wp/3._wp);l=1.e-10_wp;r=0.499999_wp
    do i=1,65
      x=(l+r)*0.5_wp
      call bbp_original_bulk(k,x,w,dwk,dwx,wn,mug,mup,ok)
      der=dwx+(mp-mn)+hbarc*(3*pi**2*nb*x)**(1._wp/3._wp)
      if(der>0)then
        r=x
      else
        l=x
      endif
    enddo
    f=w+(mp-mn)*x+0.75_wp*hbarc*(3*pi**2*nb*x)**(1._wp/3._wp)*x
  end function

  subroutine bbp_original_solve(nb,published,state,status,guess,phase_mode)
    real(wp),intent(in)::nb
    logical,intent(in)::published
    type(bbp_original_state_t),intent(out)::state
    character(len=*),intent(out)::status
    real(wp),intent(in),optional::guess(3)
    integer,intent(in),optional::phase_mode ! 0: lower energy, 1: gas-free, 2: gas-bearing.
    real(wp)::z(3),bestz(3),lo(3),hi(3),f,g(4),best,jac(3,3)
    real(wp)::gasfree_f,gasfree_z(3)
    type(bbp_original_state_t)::st,truth
    logical::ok,good,have,gasfree_ok
    integer::attempt,its,totalits,active,mode,last_attempt
    state=bbp_original_state_t();status='invalid BBP baryon density'
    if(.not.ieee_is_finite(nb))return
    if(nb<nb_lo*(1._wp-1.e-12_wp).or.nb>nb_hi*(1._wp+1.e-12_wp))return
    call bounds(nb,lo,hi);best=huge(1._wp);have=.false.;totalits=0
    mode=0
    if(present(phase_mode))mode=phase_mode
    if(mode<0.or.mode>2)then
      status='invalid BBP phase selection';return
    endif
    ! Solve the gas-free face separately. The exact functional has a
    ! nonanalytic +n_g**(2/3) surface increment, so it is a legitimate candidate.
    gasfree_z=[1.32_wp,0.31_wp,0._wp]
    call polish(nb,gasfree_z,.false.,2,its,gasfree_ok)
    call bbp_original_point(nb,gasfree_z,.false.,st,gasfree_f,g,ok)
    gasfree_ok=gasfree_ok.and.ok
    if(gasfree_ok)then
      call gradient_jacobian(nb,gasfree_z,.false.,2,jac,ok)
      gasfree_ok=ok
      if(ok)gasfree_ok=minimum_eigenvalue(jac,2)>=-1.e-7_wp
    endif
    if(.not.published.and.gasfree_ok.and.mode/=2)then
      best=gasfree_f;bestz=gasfree_z;have=.true.
    endif
    last_attempt=5
    if(mode==1)last_attempt=0
    do attempt=1,last_attempt
      z=[1.30_wp,0.28_wp,min(0.90_wp,max(0.05_wp,1._wp-2.5e-4_wp/nb))]
      select case(attempt)
      case(1)
        if(present(guess))z=guess
        z(3)=max(z(3),0.005_wp)
      case(2)
        z=[1.25_wp,0.20_wp,0.85_wp]
      case(3)
        z=[1.31_wp,0.31_wp,0.25_wp]
      case(4)
        z=[1.20_wp,0.12_wp,0.98_wp]
      case(5)
        z=[1.32_wp,0.32_wp,0.02_wp]
      end select
      z=min(max(z,lo+1.e-6_wp),hi-1.e-6_wp)
      call polish(nb,z,published,3,its,good);totalits=totalits+its
      if(.not.good)cycle
      call bbp_original_point(nb,z,published,st,f,g,ok)
      if(.not.ok)cycle
      call gradient_jacobian(nb,z,published,3,jac,ok)
      if(.not.ok)cycle
      ! A published-closure fixed point need not minimize the exact energy;
      ! retain it for comparison, but a variational candidate must be stable.
      if(.not.published.and.minimum_eigenvalue(jac,3)<-1.e-7_wp)cycle
      if(f<best)then
        best=f;bestz=z;have=.true.
      endif
      if(published.and.have)exit
    enddo
    if(.not.have)then
      status='BBP equilibrium solve failed; no certified stationary candidate';return
    endif
    call bbp_original_point(nb,bestz,published,state,f,g,ok)
    call bbp_original_point(nb,bestz,.false.,truth,f,g,good)
    state%true_gradient=g(:3)
    state%firstlaw_local=nb*(truth%mu-state%mu)/max(abs(state%pressure),tiny(1._wp))
    active=3
    if(state%gas_free)active=2
    call gradient_jacobian(nb,bestz,.false.,active,jac,good)
    if(.not.good)then
      status='BBP final curvature evaluation failed';return
    endif
    state%hessian_min=minimum_eigenvalue(jac,active)
    state%uniform_delta=uniform_energy_per_baryon(nb)-f
    state%iterations=totalits;state%converged=ok
    if(.not.ok.or.state%pressure<=0.or.state%a<=0.or.state%uniform_delta<=0)then
      status='BBP stationary state violates positive pressure or uniform-phase comparison';return
    endif
    status='ok'
  end subroutine

  subroutine bbp_original_maxwell(left,right,status)
    !! Two locally neutral spherical phases; common tangent in epsilon(nb).
    !! This construction is a derived thermodynamic completion, not a BBP
    !! claim about a physical first-order neutron-drip transition.
    type(bbp_original_state_t),intent(out)::left,right
    character(len=*),intent(out)::status
    type(bbp_original_state_t)::sa,sb,sp,sm,ta,tb
    real(wp)::y(2),yt(2),rr(3),rt(3),jac(3,3),step(3),da,db,h,alpha,norm
    integer::it,ls
    logical::ok
    character(len=160)::msg
    y=log([2.73e-4_wp,3.0e-4_wp]);status='BBP Maxwell common-tangent solve failed'
    left=bbp_original_state_t();right=bbp_original_state_t()
    do it=1,35
      call tangent_residual(y,sa,sb,rr,ok)
      if(.not.ok)return
      norm=maxval(abs(rr(:2)))
      if(norm<5.e-10_wp)then
        left=sa;right=sb;status='ok';return
      endif
      h=1.e-4_wp
      call bbp_original_solve(exp(y(1)+h),.false.,sp,msg,phase_mode=1)
      if(trim(msg)/='ok')return
      call bbp_original_solve(exp(y(1)-h),.false.,sm,msg,phase_mode=1)
      if(trim(msg)/='ok')return
      da=(sp%mu-sm%mu)/(2*h)
      call bbp_original_solve(exp(y(2)+h),.false.,sp,msg,phase_mode=2)
      if(trim(msg)/='ok')return
      call bbp_original_solve(exp(y(2)-h),.false.,sm,msg,phase_mode=2)
      if(trim(msg)/='ok')return
      db=(sp%mu-sm%mu)/(2*h)
      jac=0;jac(1,:2)=[da,-db]
      jac(2,:2)=[exp(y(1))*da,-exp(y(2))*db]/2.8e-4_wp
      call solve_linear(jac,-rr,2,step,ok)
      if(.not.ok)return
      step=step/max(1._wp,maxval(abs(step(:2)))/0.08_wp);alpha=1
      do ls=1,25
        yt=y+alpha*step(:2)
        call tangent_residual(yt,ta,tb,rt,ok)
        if(ok)then
          if(maxval(abs(rt(:2)))<norm)exit
        endif
        alpha=alpha*0.5_wp
      enddo
      if(ls>25)return
      y=yt
    enddo
  contains
    subroutine tangent_residual(yy,a,b,r,valid)
      real(wp),intent(in)::yy(2)
      type(bbp_original_state_t),intent(out)::a,b
      real(wp),intent(out)::r(3)
      logical,intent(out)::valid
      character(len=160)::message
      valid=.false.;r=huge(1._wp)
      if(yy(1)>=yy(2).or.minval(yy)<log(nb_lo).or.maxval(yy)>log(5.e-4_wp))return
      call bbp_original_solve(exp(yy(1)),.false.,a,message,phase_mode=1)
      if(trim(message)/='ok')return
      call bbp_original_solve(exp(yy(2)),.false.,b,message,phase_mode=2)
      if(trim(message)/='ok')return
      r=[a%mu-b%mu,(a%pressure-b%pressure)/2.8e-4_wp,0._wp]
      valid=.true.
    end subroutine
  end subroutine

  subroutine same_phase_firstlaw(st,published,res,ok)
    type(bbp_original_state_t),intent(in)::st
    logical,intent(in)::published
    real(wp),intent(out)::res
    logical,intent(out)::ok
    type(bbp_original_state_t)::s1,s2
    real(wp)::h,n1,n2,der,guess(3)
    integer::mode
    character(len=160)::status
    h=st%nb*1.e-4_wp;mode=2
    if(st%gas_free)mode=1
    guess=[st%k,st%x,st%gas_fraction];res=0;ok=.false.
    n1=st%nb-h;n2=st%nb+h
    if(n1<nb_lo)then
      n1=st%nb+h;n2=st%nb+2*h
    else if(n2>nb_hi)then
      n1=st%nb-h;n2=st%nb-2*h
    endif
    call bbp_original_solve(n1,published,s1,status,guess,mode)
    if(trim(status)/='ok')return
    call bbp_original_solve(n2,published,s2,status,guess,mode)
    if(trim(status)/='ok')return
    if(n1<st%nb.and.n2>st%nb)then
      der=(s2%binding-s1%binding)/(2*h)
    else if(n1>st%nb)then
      der=(-3*st%binding+4*s1%binding-s2%binding)/(2*h)
    else
      der=(3*st%binding-4*s1%binding+s2%binding)/(2*h)
    endif
    res=(st%nb**2*der-st%pressure)/st%pressure;ok=ieee_is_finite(res)
  end subroutine

  subroutine generate_bbp_original_published_native(source_path,n_points,eos,output_dir,status)
    character(len=*),intent(in)::source_path,output_dir
    integer,intent(in)::n_points
    type(eos_data_t),intent(out)::eos
    character(len=*),intent(out)::status
    call generate(source_path,n_points,.true.,eos,output_dir,status)
  end subroutine
  subroutine generate_bbp_original_variational_native(source_path,n_points,eos,output_dir,status)
    character(len=*),intent(in)::source_path,output_dir
    integer,intent(in)::n_points
    type(eos_data_t),intent(out)::eos
    character(len=*),intent(out)::status
    call generate(source_path,n_points,.false.,eos,output_dir,status)
  end subroutine
  subroutine generate_bbp_variational_maxwell_native(source_path,n_points,eos,output_dir,status)
    character(len=*),intent(in)::source_path,output_dir
    integer,intent(in)::n_points
    type(eos_data_t),intent(out)::eos
    character(len=*),intent(out)::status
    call generate(source_path,n_points,.false.,eos,output_dir,status,.true.)
  end subroutine

  subroutine read_source(path,status)
    character(len=*),intent(in)::path
    character(len=*),intent(out)::status
    character(len=1024)::line
    character(len=64)::name,extra,num
    character(len=40),parameter::names(nparams)=[character(len=40):: &
      'w0_MeV','k0_fm^-1','K_MeV','s_MeV','sigma_MeV','tn_coefficient_MeV_fm2', &
      'neutron_a','neutron_b','neutron_c','proton_a','proton_b','proton_c', &
      'thickness_coefficient','hbarc_MeV_fm','e2_MeV_fm','mn_c2_MeV','mp_c2_MeV', &
      'nb_min_fm^-3','nb_max_fm^-3','include_exchange','include_screening']
    real(wp),parameter::values(nparams)=[w0,k0,compress,symmetry,sigma,tncoef, &
      40.4_wp,1.088_wp,2.545_wp,218._wp,277._wp,8.57_wp,0.74_wp,hbarc,e2,mn,mp,nb_lo,nb_hi,0._wp,0._wp]
    real(wp)::v
    integer::unit,ios,i,check,p,ep,power
    status='BBP original parameter source open failed'
    open(newunit=unit,file=trim(path),status='old',action='read',iostat=ios)
    if(ios/=0)return
    read(unit,'(a)',iostat=ios)line
    if(ios/=0.or.trim(line)/='BBP1971_ORIGINAL_PARAMETERS_V1')then
      status='BBP original parameter source has wrong identity';close(unit);return
    endif
    do i=1,nparams
      read(unit,'(a)',iostat=ios)line
      if(ios/=0)then
        status='BBP original parameter source is truncated';close(unit);return
      endif
      name='';extra='';v=0
      line=adjustl(line)
      p=scan(trim(line),' '//achar(9))
      if(p==0.or.p>64)then
        status='BBP original parameter source requires two whitespace tokens';close(unit);return
      endif
      name=line(:p-1);num=adjustl(line(p:))
      if(len_trim(line(p:))>64.or.len_trim(num)==0)then
        status='BBP original parameter source has invalid numeric token';close(unit);return
      endif
      if(verify(trim(num),'0123456789.+-eE')/=0)then
        status='BBP original parameter source has invalid numeric token';close(unit);return
      endif
      ! Reject exponents outside the fixed model's parameter scale BEFORE a
      ! real read, which otherwise can raise SIGFPE under overflow trapping.
      ep=scan(trim(num),'eE')
      if(ep>0)then
        if(len_trim(num)-ep>4.or.len_trim(num)==ep)then
          status='BBP original parameter exponent outside supported range';close(unit);return
        endif
        power=0;read(num(ep+1:),*,iostat=ios)power
        if(ios/=0.or.abs(power)>32)then
          status='BBP original parameter exponent outside supported range';close(unit);return
        endif
      endif
      read(line,*,iostat=ios)name,v
      read(line,*,iostat=check)name,v,extra
      if(ios/=0.or.check>=0.or..not.ieee_is_finite(v))then
        status='BBP original parameter source malformed or nonfinite';close(unit);return
      endif
      if(trim(name)/=trim(names(i)).or.abs(v-values(i))>1.e-13_wp*max(abs(values(i)),1._wp))then
        status='BBP original parameter source differs from declared model';close(unit);return
      endif
    enddo
    read(unit,'(a)',iostat=ios)line
    close(unit)
    if(ios>=0)then
      status='BBP original parameter source has extra records';return
    endif
    status='ok'
  end subroutine

  subroutine generate(source_path,n_points,published,eos,output_dir,status,with_maxwell)
    character(len=*),intent(in)::source_path,output_dir
    integer,intent(in)::n_points
    logical,intent(in)::published
    type(eos_data_t),intent(out)::eos
    character(len=*),intent(out)::status
    logical,intent(in),optional::with_maxwell
    type(bbp_original_state_t),allocatable::states(:)
    type(bbp_original_state_t)::st,left,right
    real(wp),allocatable::grid(:)
    real(wp)::nb,guess(3),ni,ng,res,maxres,maxgrad,maxtrue,minh,mingap,pc,mc,vtemp
    integer::i,j,unit,ios,closeios,gasfree_count,rows,removed,inserted,phase_count,mode
    logical::ok,convex
    character(len=80)::model
    status='BBP generation requires 3..100001 density points'
    if(n_points<3.or.n_points>100001)return
    call read_source(source_path,status)
    if(trim(status)/='ok')return
    convex=.false.
    if(present(with_maxwell))convex=with_maxwell
    allocate(grid(n_points+6));rows=0;removed=0;inserted=0
    if(convex)then
      call bbp_original_maxwell(left,right,status)
      if(trim(status)/='ok')return
      pc=0.5_wp*(left%pressure+right%pressure);mc=0.5_wp*(left%mu+right%mu)
    endif
    do i=1,n_points
      nb=exp(log(nb_lo)+(log(nb_hi)-log(nb_lo))*real(i-1,wp)/real(n_points-1,wp))
      if(i==1)nb=nb_lo
      if(i==n_points)nb=nb_hi
      if(convex)then
        if(nb>left%nb.and.nb<right%nb)then
          removed=removed+1;cycle
        endif
      endif
      rows=rows+1;grid(rows)=nb
    enddo
    if(convex)then
      call insert_grid(left%nb);call insert_grid(right%nb)
      ! At least three one-sided samples per phase, even for sparse requests.
      if(count(grid(:rows)<=left%nb)<3)call insert_grid((2*nb_lo+left%nb)/3)
      if(count(grid(:rows)<=left%nb)<3)call insert_grid((nb_lo+2*left%nb)/3)
      if(count(grid(:rows)>=right%nb)<3)call insert_grid((2*right%nb+nb_hi)/3)
      if(count(grid(:rows)>=right%nb)<3)call insert_grid((right%nb+2*nb_hi)/3)
    endif
    do i=2,rows
      vtemp=grid(i);j=i-1
      do while(j>=1)
        if(grid(j)<=vtemp)exit
        grid(j+1)=grid(j);j=j-1
      enddo
      grid(j+1)=vtemp
    enddo
    allocate(states(rows));call allocate_eos(eos,rows)
    eos%mn=mn;eos%mp=mp;eos%composition_file_present=.true.;eos%composition_available=.true.
    eos%beta_equilibrated=.true.
    model='BBP1971_original_variational'
    if(published)model='BBP1971_original_published'
    if(convex)model='BBP1971_variational_maxwell'
    maxgrad=0;maxtrue=0;minh=huge(1._wp);mingap=huge(1._wp);gasfree_count=0
    guess=[1.32_wp,0.31_wp,0.02_wp]
    do i=1,rows
      nb=grid(i);mode=0
      if(convex)then
        mode=2
        if(nb<=left%nb)mode=1
      endif
      call bbp_original_solve(nb,published,st,status,guess,mode)
      if(trim(status)/='ok')return
      if(convex)then
        if(nb==left%nb.or.nb==right%nb)then
          ! Common tangent represented exactly at its numerically solved
          ! endpoints. Unadjusted values are retained in transitions.csv.
          st%pressure=pc;st%mu=mc;st%energy=mc*nb-pc
        endif
      endif
      states(i)=st;guess=[st%k,st%x,st%gas_fraction]
      eos%nb(i)=nb;eos%q(1,i)=st%pressure/nb;eos%q(3,i)=st%mu/mn-1
      eos%q(4,i)=-st%mue/mn
      ! CompOSE q5 is the conserved lepton-number potential, not the electron
      ! chemical potential. Cold neutrino-free beta equilibrium has mu_L=0.
      eos%q(5,i)=0._wp
      eos%q(6,i)=st%energy/(nb*mn)-1;eos%q(7,i)=eos%q(6,i)
      ni=st%k**3/(1.5_wp*pi**2);ng=nb*st%gas_fraction
      eos%yn(i)=(1-st%u)*ng/nb;eos%yp(i)=0;eos%ye(i)=st%u*ni*st%x/nb;eos%ymu(i)=0
      eos%yh(i)=st%u*ni/nb;eos%abar(i)=st%a;eos%zbar(i)=st%z
      eos%row_status(i)=2;eos%iphase(i)=101
      if(.not.st%gas_free)eos%iphase(i)=102
      maxgrad=max(maxgrad,maxval(abs(st%gradient)))
      maxtrue=max(maxtrue,maxval(abs(st%true_gradient)))
      minh=min(minh,st%hessian_min);mingap=min(mingap,st%uniform_delta)
      if(st%gas_free)gasfree_count=gasfree_count+1
    enddo
    eos%yn_present=.true.;eos%yp_present=.true.;eos%ye_present=.true.;eos%ymu_present=.true.
    eos%abar_present=.true.;eos%zbar_present=.true.;eos%yh_present=.true.
    eos%c3_applicable=.true.;eos%c4_applicable=.true.
    call validate_eos(eos,ok,status)
    if(.not.ok)return
    phase_count=1
    if(gasfree_count>0.and.gasfree_count<rows)phase_count=2
    open(newunit=unit,file=trim(output_dir)//'/bbp_original_states.csv',status='replace',action='write',iostat=ios)
    if(ios/=0)then
      status='BBP state diagnostics open failed';return
    endif
    write(unit,'(a)',iostat=ios)'nb_fm-3,epsilon_MeV_fm-3,P_MeV_fm-3,mu_MeV,k_fm-1,x,kg_fm-1,'// &
      'gas_density_over_nb,u,A,Z,mu_e_MeV,ws_MeV,wc_MeV,wthick_MeV,'// &
      'closure_g_k,closure_g_x,closure_g_t,true_g_k,true_g_x,true_g_t,hessian_min,'// &
      'uniform_minus_crust_MeV,same_phase_firstlaw_relative,iterations,gas_free'
    maxres=0
    do i=1,rows
      call same_phase_firstlaw(states(i),published,res,ok)
      if(.not.ok)then
        close(unit);status='BBP independent same-phase derivative calculation failed';return
      endif
      maxres=max(maxres,abs(res));st=states(i)
      write(unit,'(24(es24.16e3,","),i0,",",l1)',iostat=ios) &
        st%nb,st%energy,st%pressure,st%mu,st%k,st%x,st%kg,st%gas_fraction,st%u,st%a,st%z, &
        st%mue,st%ws,st%wc,st%wthick,st%gradient,st%true_gradient,st%hessian_min,st%uniform_delta,res, &
        st%iterations,st%gas_free
      if(ios/=0)exit
    enddo
    close(unit,iostat=closeios)
    if(ios/=0.or.closeios/=0)then
      status='BBP state diagnostics write failed';return
    endif
    open(newunit=unit,file=trim(output_dir)//'/bbp_original_transitions.csv', &
      status='replace',action='write',iostat=ios)
    if(ios/=0)then
      status='BBP transition diagnostics open failed';return
    endif
    write(unit,'(a)',iostat=ios)'nb_low,nb_high,P_common,mu_common,epsilon_low_raw,epsilon_high_raw,'// &
      'P_low_raw,P_high_raw,mu_low_raw,mu_high_raw'
    if(convex)write(unit,'(10(es24.16e3,:,","))',iostat=ios)left%nb,right%nb,pc,mc, &
      left%energy,right%energy,left%pressure,right%pressure,left%mu,right%mu
    close(unit,iostat=closeios)
    if(ios/=0.or.closeios/=0)then
      status='BBP transition diagnostics write failed';return
    endif
    open(newunit=unit,file=trim(output_dir)//'/bbp_original_diagnostics.json', &
      status='replace',action='write',iostat=ios)
    if(ios/=0)then
      status='BBP original diagnostics open failed';return
    endif
    write(unit,'(a)',iostat=ios)'{'
    write(unit,'(a)',iostat=ios)'  "model": "'//trim(model)//'",'
    write(unit,'(a)',iostat=ios)'  "source_doi": "10.1016/0375-9474(71)90281-8",'
    write(unit,'(a)',iostat=ios)'  "source_pdf_sha256": "f5d15dbc97c8dae24557b55844c2bfc550984aaa50e514e8d9412e11a05ca495",'
    write(unit,'(a)',iostat=ios)'  "equations": "3.8,3.9,3.18,3.19,4.30,5.3,5.9,5.13,6.4; p251 derivative convention",'
    write(unit,'(a)',iostat=ios)'  "representation": "original spherical compressible liquid-drop energy, native stationary solve",'
    write(unit,'(a)',iostat=ios)'  "exchange_and_screening": "omitted as described on pp246-248",'
    write(unit,'(a)',iostat=ios)'  "constants": "modern hbarc,e2,nucleon masses; printed kinetic coefficient19.74 retained",'
    write(unit,'(a)',iostat=ios)'  "bounds": "spherical inner crust 0.000267 <= nb <= 0.08 fm^-3; no pasta extrapolation",'
    write(unit,'(a)',iostat=ios)'  "phase_comparison": "finite stationary candidates plus uniform npe trial; no global spherical-branch or convex-hull proof",'
    write(unit,'(a)',iostat=ios)'  "stationary_search": "one gas-free Newton seed and five gas-bearing seeds; published closure takes first converged root",'
    write(unit,'(a)',iostat=ios)'  "coordinate_box": "k: max(1.05,(1.5*pi^2*nb*1.000001)^(1/3))..1.5; x:0.025..0.45; t=ng/nb:1e-9..0.999999; gas-free t=0 separately",'
    write(unit,'(a)',iostat=ios)'  "maxwell_search": "selected gas-free/gas-bearing branches; seed densities(2.73e-4,3.0e-4); bracket nb_min<=n_low<n_high<=5e-4; not global hull",'
    write(unit,'(a)',iostat=ios)'  "gradient_semantics": "three active coordinates inside; gas-free face has two active coordinates and true_g_t=0 is an inactive placeholder, not its divergent one-sided derivative",'
    write(unit,'(a)',iostat=ios)'  "table_reproduction": "nine Table1/2 anchors retained in benchmark; last nb=0.0789 gives P=0.40107766 vs printed0.384 (+4.4473%), A=1027.55 vs990, x=0.118809 vs0.121; unresolved beyond-rounding mismatch",'
    if(published)then
      write(unit,'(a)',iostat=ios)'  "closure": "published p251 omission of dWg/dng only in surface derivative",'
    else if(convex)then
      write(unit,'(a)',iostat=ios)'  "closure": "full original energy derivatives with explicit common-tangent phase endpoints",'
    else
      write(unit,'(a)',iostat=ios)'  "closure": "all derivatives of original energy retained; pointwise branch minimum",'
    endif
    write(unit,'(a)',iostat=ios)'  "pressure_basis": "nb^2 times selected constrained derivative of epsilon/nb; published gas-surface omission propagated",'
    write(unit,'(a)',iostat=ios)'  "chemical_potential_basis": "Euler (epsilon+P)/nb under selected derivative closure",'
    write(unit,'(a,i0,a)',iostat=ios)'  "requested_points": ',n_points,','
    write(unit,'(a,i0,a)',iostat=ios)'  "rows": ',rows,','
    write(unit,'(a,i0,a)',iostat=ios)'  "phase_count": ',phase_count,','
    write(unit,'(a,i0,a)',iostat=ios)'  "removed_grid_points": ',removed,','
    write(unit,'(a,i0,a)',iostat=ios)'  "inserted_grid_points": ',inserted,','
    write(unit,'(a,i0,a)',iostat=ios)'  "gas_free_rows": ',gasfree_count,','
    write(unit,'(a,es24.16e3,a)',iostat=ios)'  "max_closure_gradient": ',maxgrad,','
    write(unit,'(a,es24.16e3,a)',iostat=ios)'  "max_true_energy_gradient": ',maxtrue,','
    write(unit,'(a,es24.16e3,a)',iostat=ios)'  "minimum_true_hessian_eigenvalue": ',minh,','
    write(unit,'(a,es24.16e3,a)',iostat=ios)'  "minimum_uniform_energy_gap_MeV": ',mingap,','
    write(unit,'(a,es24.16e3,a)',iostat=ios)'  "max_same_phase_firstlaw_relative": ',maxres,','
    write(unit,'(a)',iostat=ios)'  "firstlaw_note": "independently re-solved same-phase finite difference, delta_nb/nb=1e-4; phase gaps are explicit",'
    write(unit,'(a)',iostat=ios)'  "success_meaning": "native states generated; matching admission remains independent"'
    write(unit,'(a)',iostat=ios)'}'
    close(unit,iostat=closeios)
    if(ios/=0.or.closeios/=0)then
      status='BBP original diagnostics write failed';return
    endif
    status='ok'
  contains
    subroutine insert_grid(value)
      real(wp),intent(in)::value
      if(any(abs(grid(:rows)-value)<=32*epsilon(1._wp)*value))return
      rows=rows+1;grid(rows)=value;inserted=inserted+1
    end subroutine
  end subroutine
end module compose_eos_bbp_original
