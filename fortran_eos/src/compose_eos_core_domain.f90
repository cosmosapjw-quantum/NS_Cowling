! Local causal-domain diagnostics; does not alter or extend the source functional.
! n: fm^-3; epsilon/P: MeV fm^-3; mu: MeV; c_s^2/c^2 dimensionless.
module compose_eos_core_domain
 use compose_eos_types, only: wp
 use compose_eos_core_native, only: core_state_t,solve_core_state,core_model_bounds
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 real(wp),parameter::hc=197.3269804_wp,pi=acos(-1._wp)
 public::core_sound_speed_fd,core_sound_speed_analytic,core_causal_bracket,core_sound_speed_gibbs
contains
 subroutine check_density_support(model,n,lo,hi,status)
  character(len=*),intent(in)::model
  real(wp),intent(in)::n
  real(wp),intent(out)::lo,hi
  character(len=*),intent(out)::status
  real(wp)::mn,mp
  lo=0;hi=0;status='density must be finite'
  if(.not.ieee_is_finite(n))return
  call core_model_bounds(model,lo,hi,mn,mp,status);if(status/='ok')return
  status='density outside declared native core model bounds'
  if(n<lo.or.n>hi)return
  status='ok'
 end subroutine
 subroutine core_sound_speed_fd(model,n,rel,cs,cs_mu,firstlaw,status)
  character(len=*),intent(in)::model
  real(wp),intent(in)::n,rel
  real(wp),intent(out)::cs,cs_mu,firstlaw
  character(len=*),intent(out)::status
  type(core_state_t)::s(-2:2)
  real(wp),parameter::bbb_knots(16)=[.08_wp,.16_wp,.2_wp,.3_wp,.4_wp,.5_wp,.6_wp,.7_wp,.8_wp,.9_wp,1._wp,1.1_wp,1.2_wp,1.3_wp,1.4_wp,1.5_wp]
  real(wp)::h,de,dp,dm,nlo,nhi
  integer::i
  cs=0;cs_mu=0;firstlaw=0
  status='invalid density or finite-difference step'
  if(.not.all(ieee_is_finite([n,rel])))return
  if(n<=0.or.rel<=0.or.rel>=.25_wp)return
  call check_density_support(model,n,nlo,nhi,status);if(status/='ok')return
  h=n*rel
  status='finite-difference stencil outside declared native core model bounds'
  if(n-2*h<nlo.or.n+2*h>nhi)return
  if(model=='BBB_BHF2')then
   status='stencil crosses BBB source-energy knot; curvature is branch-sided'
   if(any(bbb_knots>=n-2*h.and.bbb_knots<=n+2*h))return
  endif
  do i=-2,2
   call solve_core_state(model,n+i*h,s(i),status)
   if(status/='ok')return
  enddo
  status='stencil crosses a phase boundary; use branch-sided derivative'
  if(any(s%phase/=s(0)%phase))return
  status='stencil crosses muon onset; reduce step or use branch-sided derivative'
  if(any((s%ymu>0).neqv.(s(0)%ymu>0)))return
  de=(s(-2)%epsilon-8*s(-1)%epsilon+8*s(1)%epsilon-s(2)%epsilon)/(12*h)
  dp=(s(-2)%pressure-8*s(-1)%pressure+8*s(1)%pressure-s(2)%pressure)/(12*h)
  dm=(s(-2)%mu_n-8*s(-1)%mu_n+8*s(1)%mu_n-s(2)%mu_n)/(12*h)
  status='nonpositive energy derivative';if(de<=0.or.s(0)%mu_n<=0)return
  cs=dp/de;cs_mu=n*dm/s(0)%mu_n;firstlaw=abs(de/s(0)%mu_n-1)
  status='nonfinite sound-speed result'
  if(.not.all(ieee_is_finite([cs,cs_mu,firstlaw])))return
  status='ok'
 end subroutine
 subroutine shape(q,y,f)
  real(wp),intent(in)::q,y
  real(wp),intent(out)::f(3)
  f=[(1-y)**q+y**q,q*(y**(q-1)-(1-y)**(q-1)), &
      q*(q-1)*(y**(q-2)+(1-y)**(q-2))]
 end subroutine
 subroutine add(v,r,f)
  real(wp),intent(inout)::v(6)
  real(wp),intent(in)::r(3),f(3)
  v=v+[r(1)*f(1),r(2)*f(1),r(1)*f(2),r(3)*f(1),r(2)*f(2),r(1)*f(3)]
 end subroutine
 function power(n,q) result(r)
  real(wp),intent(in)::n,q
  real(wp)::r(3)
  r=[n**q,q*n**(q-1),q*(q-1)*n**(q-2)]
 end function
 subroutine nuclear_hessian(model,n,y,branch,v,status)
  ! Independent density/composition representation of published energy functionals.
  character(len=*),intent(in)::model
  real(wp),intent(in)::n,y
  integer,intent(in)::branch
  real(wp),intent(out)::v(6)
  character(len=*),intent(out)::status
  real(wp),parameter::a(21)=[337.2_wp,-382._wp,89.8_wp,.457_wp,-59._wp,-19.1_wp,214.6_wp, &
    -384._wp,6.4_wp,69._wp,-33._wp,.35_wp,0._wp,0._wp,287._wp,-1.54_wp,175._wp,-1.45_wp,.32_wp,.195_wp,0._wp]
  real(wp)::t0,t1,t2,t3,x0,x1,x2,x3,alpha,c0,c1,k,f5(3),f8(3),f2(3),f(3),r(3),g1(3),g2(3)
  real(wp)::ex,dx,ddx,d,poly,dpoly,ddpoly,lambda,b(3),mn,mp
  integer::j
  v=0;status='analytic Hessian supports APR/SLy4/KDE0v1 single phases'
  if(model/='APR_APR'.and.model/='RG_SLY4'.and.model/='RG_KDE0v1')return
  status='analytic Hessian requires interior composition'
  if(y<=0.or.y>=1.or.n<=0)return
  k=.6_wp*(3*pi*pi)**(2._wp/3)
  call shape(5._wp/3,y,f5);call shape(8._wp/3,y,f8);call shape(2._wp,y,f2)
  mn=939.5654_wp;mp=938.2720_wp
  if(model=='APR_APR')then
   mn=939.5654133_wp;mp=938.2720813_wp
   if(branch/=1.and.branch/=2)return
  endif
  call add(v,power(n,1._wp),[mn*(1-y)+mp*y,mp-mn,0._wp])
  if(model/='APR_APR')then
   if(model=='RG_KDE0v1')then
    t0=-2553.08_wp;t1=411.70_wp;t2=-419.87_wp;t3=14603.6_wp
    x0=.6483_wp;x1=-.3472_wp;x2=-.9268_wp;x3=.9475_wp;alpha=.1673_wp
   else
    t0=-2488.91_wp;t1=486.82_wp;t2=-546.39_wp;t3=13777._wp
    x0=.834_wp;x1=-.344_wp;x2=-1._wp;x3=1.354_wp;alpha=1._wp/6
   endif
   c0=(t1*(2+x1)+t2*(2+x2))/8
   c1=-(t1*(2*x1+1)-t2*(2*x2+1))/8
   call add(v,20.73553_wp*k*power(n,5._wp/3),f5)
   f=-(2*x0+1)*f2;f(1)=f(1)+(2+x0)
   call add(v,t0/4*power(n,2._wp),f)
   f=-(2*x3+1)*f2;f(1)=f(1)+(2+x3)
   call add(v,t3/24*power(n,2+alpha),f)
   call add(v,k*power(n,8._wp/3),c0*f5+c1*f8)
  else
   f=[hc**2/(2*mn)*(1-y)**(5._wp/3)+hc**2/(2*mp)*y**(5._wp/3), &
     (5._wp/3)*(hc**2/(2*mp)*y**(2._wp/3)-hc**2/(2*mn)*(1-y)**(2._wp/3)), &
     (10._wp/9)*(hc**2/(2*mp)*y**(-1._wp/3)+hc**2/(2*mn)*(1-y)**(-1._wp/3))]
   call add(v,k*power(n,5._wp/3),f)
   b=power(n,8._wp/3);ex=exp(-a(4)*n)
   r=ex*[b(1),b(2)-a(4)*b(1),b(3)-2*a(4)*b(2)+a(4)**2*b(1)]
   call add(v,k*r,a(3)*f5+a(5)*f8)
   ex=exp(-(a(9)*n)**2);dx=-2*a(9)**2*n*ex
   ddx=(-2*a(9)**2+4*a(9)**4*n*n)*ex
   g1=[-a(1)-a(2)*n-a(6)*n*n-(a(10)+a(11)*n)*ex, &
       -a(2)-2*a(6)*n-a(11)*ex-(a(10)+a(11)*n)*dx, &
       -2*a(6)-2*a(11)*dx-(a(10)+a(11)*n)*ddx]
   g2=[-a(12)/n-a(7)-a(8)*n-a(13)*ex,a(12)/n**2-a(8)-a(13)*dx,-2*a(12)/n**3-a(13)*ddx]
   if(branch==2)then
    do j=1,2
     if(j==1)then
      d=n-a(19);lambda=a(18);poly=a(17)*d+a(21)*d*d;dpoly=a(17)+2*a(21)*d;ddpoly=2*a(21)
     else
      d=n-a(20);lambda=a(16);poly=a(15)*d+a(14)*d*d;dpoly=a(15)+2*a(14)*d;ddpoly=2*a(14)
     endif
     r=exp(lambda*d)*[poly,dpoly+lambda*poly,ddpoly+2*lambda*dpoly+lambda**2*poly]
     if(j==1)g1=g1-r
     if(j==2)g2=g2-r
    enddo
   endif
   r=[n*n*g1(1),2*n*g1(1)+n*n*g1(2),2*g1(1)+4*n*g1(2)+n*n*g1(3)]
   call add(v,r,[4*y*(1-y),4*(1-2*y),-8._wp])
   r=[n*n*g2(1),2*n*g2(1)+n*n*g2(2),2*g2(1)+4*n*g2(2)+n*n*g2(3)]
   call add(v,r,[(1-2*y)**2,-4*(1-2*y),8._wp])
  endif
  status='ok'
 end subroutine
 subroutine core_sound_speed_analytic(model,n,cs,frozen,hyy,beta_residual,status)
  character(len=*),intent(in)::model
  real(wp),intent(in)::n
  real(wp),intent(out)::cs,frozen,hyy,beta_residual
  character(len=*),intent(out)::status
  type(core_state_t)::s
  real(wp)::v(6),sus,mu,y,hny,hnn,m,nlo,nhi
  integer::branch,j
  cs=0;frozen=0;hyy=0;beta_residual=0
  call check_density_support(model,n,nlo,nhi,status);if(status/='ok')return
  call solve_core_state(model,n,s,status);if(status/='ok')return
  status='mixed-phase analytic Hessian requires constrained Gibbs response'
  if(s%phase==8)return
  branch=1;if(s%phase==7)branch=2
  y=s%yp;mu=s%mu_e
  call nuclear_hessian(model,n,y,branch,v,status);if(status/='ok')return
  sus=0
  do j=1,2
   m=.51099895_wp;if(j==2)m=105.6583755_wp
   if(mu>m)sus=sus+mu*sqrt(mu*mu-m*m)/(pi*pi*hc**3)
  enddo
  status='nonpositive lepton susceptibility';if(sus<=0)return
  hnn=v(4)+y*y/sus;hny=v(5)+mu+n*y/sus;hyy=v(6)+n*n/sus
  beta_residual=v(3)+n*mu
  status='nonpositive composition Hessian';if(hyy<=0)return
  frozen=n*hnn/s%mu_n
  cs=n*(hnn-hny*hny/hyy)/s%mu_n
  status='ok'
 end subroutine
 subroutine core_sound_speed_gibbs(n,cs,status)
  ! Two conserved charges, global neutrality, homogeneous leptons.
  ! Coexistence tangent v=(Delta q,-Delta n); dmu_B/dn=v_B^2/(v.K.v).
  real(wp),intent(in)::n
  real(wp),intent(out)::cs
  character(len=*),intent(out)::status
  type(core_state_t)::s
  real(wp)::v(6),kmat(2,2),h(2,2),nn,y,weight,det,dn,dq,mu,sus,m,denom,nlo,nhi
  integer::i
  cs=0
  call check_density_support('APR_APR',n,nlo,nhi,status);if(status/='ok')return
  call solve_core_state('APR_APR',n,s,status);if(status/='ok')return
  status='Gibbs response requires an interior APR mixed-phase state';if(s%phase/=8)return
  kmat=0
  do i=1,2
   nn=s%n_low;y=s%yp_low;weight=1-s%chi
   if(i==2)then
    nn=s%n_high;y=s%yp_high;weight=s%chi
   endif
   call nuclear_hessian('APR_APR',nn,y,i,v,status);if(status/='ok')return
   h(1,1)=v(4)-2*y*v(5)/nn+y*y*v(6)/(nn*nn)+2*y*v(3)/(nn*nn)
   h(1,2)=v(5)/nn-(v(3)+y*v(6))/(nn*nn);h(2,1)=h(1,2)
   h(2,2)=v(6)/(nn*nn);det=h(1,1)*h(2,2)-h(1,2)**2
   status='nonpositive APR susceptibility Hessian'
   if(det<=0.or.min(h(1,1),h(2,2))<=0)return
   kmat(1,1)=kmat(1,1)+weight*h(2,2)/det
   kmat(1,2)=kmat(1,2)-weight*h(1,2)/det
   kmat(2,2)=kmat(2,2)+weight*h(1,1)/det
  enddo
  mu=s%mu_e;sus=0
  do i=1,2
   m=.51099895_wp;if(i==2)m=105.6583755_wp
   if(mu>m)sus=sus+mu*sqrt(mu*mu-m*m)/(pi*pi*hc**3)
  enddo
  kmat(2,2)=kmat(2,2)+sus
  dn=s%n_high-s%n_low;dq=s%n_high*s%yp_high-s%n_low*s%yp_low
  denom=dq*dq*kmat(1,1)-2*dq*dn*kmat(1,2)+dn*dn*kmat(2,2)
  status='nonpositive Gibbs tangent susceptibility';if(denom<=0)return
  cs=n/s%mu_n*dq*dq/denom;status='ok'
 end subroutine
 subroutine core_causal_bracket(model,left,right,rel,tol,lo,hi,status)
  character(len=*),intent(in)::model
  real(wp),intent(in)::left,right,rel,tol
  real(wp),intent(out)::lo,hi
  character(len=*),intent(out)::status
  real(wp)::a,b,c,m,u,f
  integer::i
  lo=left;hi=right;status='invalid crossing bracket/tolerance'
  if(.not.all(ieee_is_finite([left,right,rel,tol])))return
  if(left<=0.or.right<=left.or.tol<=0)return
  call core_sound_speed_fd(model,lo,rel,a,u,f,status);if(status/='ok')return
  call core_sound_speed_fd(model,hi,rel,b,u,f,status);if(status/='ok')return
  status='causal crossing is not bracketed';if(a>=1.or.b<=1)return
  do i=1,80
   if(hi-lo<=tol)exit
   m=(lo+hi)/2
   call core_sound_speed_fd(model,m,rel,c,u,f,status);if(status/='ok')return
   if(c<1)then
    lo=m
   else
    hi=m
   endif
  enddo
  status='ok';if(hi-lo>tol)status='causal bisection tolerance not reached'
 end subroutine
end module
