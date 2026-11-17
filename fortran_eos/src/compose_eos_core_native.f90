! Native, cold composition-aware core generators.
! Units: n fm^-3; chemical potentials MeV; epsilon and pressure MeV fm^-3.
! hbar*c is retained explicitly. No final CompOSE EOS table is read.
! APR: Akmal et al. PRC58,1804; Schneider et al. PRC100,025803 (21 parameters).
! BBB: Baldo/Bombaci/Burgio1997 Tables1/3 published SNM/PNM energy reconstruction,
! not a new microscopic BHF solve. Source calculations stop at .8 fm^-3;
! published rounded/extrapolated source-energy data extend to 1.5 fm^-3.
! BBB input JSON SHA256 b0d646a7b334d824d52af4c8bdd255f937f348f390c41e7ec2ba76736982234c.
! GM1: Glendenning/Moszkowski1991 TableII (rounded original couplings), npe.
! KDE0v1: nucl-th/0505071 Table3; SLy4: astro-ph/0006135 Table1.
! Published parameter precision, kinetic/rest-mass conventions match the
! separately retained Python research implementations; Python is a test oracle.
module compose_eos_core_native
 use compose_eos_types, only: wp,eos_data_t,allocate_eos
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 real(wp),parameter::hc=197.3269804_wp,pi=acos(-1._wp)
 real(wp),parameter::me=.51099895_wp,mmu=105.6583755_wp
 type,public::core_state_t
  real(wp)::n_b=0,epsilon=0,pressure=0,mu_n=0,mu_p=0,mu_e=0
  real(wp)::yn=0,yp=0,ye=0,ymu=0,chi=0,n_low=0,n_high=0,yp_low=0,yp_high=0
  real(wp)::gibbs_residual=0,charge_residual=0,euler_residual=0
  integer::phase=6
 end type
 type::nuclear_t
  real(wp)::e=0,p=0,mn=0,mp=0
 end type
 public::core_model_bounds,solve_core_state,generate_core_native,solve_apr_gibbs_fraction
 public::core_nuclear_state,core_free_fermion
contains
 subroutine core_model_bounds(model,lo,hi,mn,mp,status)
  character(len=*),intent(in)::model
  real(wp),intent(out)::lo,hi,mn,mp
  character(len=*),intent(out)::status
  lo=0;hi=0;mn=939;mp=939;status='ok'
  select case(trim(model))
  case('APR_APR');lo=.085_wp;hi=1.34_wp;mn=939.5654133_wp;mp=938.2720813_wp
  case('BBB_BHF2');lo=.08_wp;hi=1.5_wp
  ! These three analytic functionals also define a homogeneous
  ! beta-equilibrated branch below the customary .08 default sampling limit.
  ! Computational support is not a claim that uniform matter is the global
  ! ground state: mechanical stability, causality and crust coexistence are
  ! checked independently. No table extrapolation or pressure clipping occurs.
  case('GM_GM1');lo=.02_wp;hi=1._wp
  case('RG_KDE0v1');lo=.02_wp;hi=1.8_wp;mn=939.5654_wp;mp=938.2720_wp
  case('RG_SLY4');lo=.02_wp;hi=1.4_wp;mn=939.5654_wp;mp=938.2720_wp
  case default;status='unknown native core model: '//trim(model)
  end select
 end subroutine
 subroutine core_free_fermion(n,m,e,p,mu)
  real(wp),intent(in)::n,m
  real(wp),intent(out)::e,p,mu
  real(wp)::x,x2,c,power,se,sp,root
  integer::k
  e=0;p=0;mu=m
  if(n<=0)return
  x=hc*(3*pi*pi*n)**(1._wp/3)/m;x2=x*x;root=sqrt(1+x2);mu=m*root
  if(x<.2_wp)then
   ! Integrated binomial series avoids threshold cancellation; at x<.2
   ! 20 terms are well below binary64 rounding, including tiny n.
   c=1;power=1;se=1;sp=0
   do k=1,20
    c=c*(1.5_wp-real(k,wp))/real(k,wp);power=power*x2
    se=se+c*power*3/(2*k+3);sp=sp+c*power*(2*k)/(2*k+3)
   enddo
   e=m*n*se;p=m*n*sp
  else
   e=m**4/(8*pi*pi*hc**3)*(x*root*(1+2*x2)-asinh(x))
   p=m**4/(24*pi*pi*hc**3)*(x*root*(2*x2-3)+3*asinh(x))
  endif
 end subroutine
 subroutine lepton(mu,m,n,e,p)
  real(wp),intent(in)::mu,m
  real(wp),intent(out)::n,e,p
  real(wp)::unused,k
  n=0;e=0;p=0
  if(mu<=m)return
  k=sqrt((mu-m)*(mu+m));n=(k/hc)**3/(3*pi*pi)
  call core_free_fermion(n,m,e,p,unused)
 end subroutine
 subroutine leptons(mu,muons,ne,nm,e,p)
  real(wp),intent(in)::mu
  logical,intent(in)::muons
  real(wp),intent(out)::ne,nm,e,p
  real(wp)::em,pm
  call lepton(mu,me,ne,e,p);nm=0
  if(muons)then
   call lepton(mu,mmu,nm,em,pm);e=e+em;p=p+pm
  endif
 end subroutine
 subroutine apr_nuclear(nn,np,branch,s,status)
  real(wp),intent(in)::nn,np
  integer,intent(in)::branch
  type(nuclear_t),intent(out)::s
  character(len=*),intent(out)::status
  real(wp),parameter::a(21)=[337.2_wp,-382._wp,89.8_wp,.457_wp,-59._wp,-19.1_wp,214.6_wp,&
   -384._wp,6.4_wp,69._wp,-33._wp,.35_wp,0._wp,0._wp,287._wp,-1.54_wp,175._wp,-1.45_wp,.32_wp,.195_wp,0._wp]
  real(wp),parameter::mn=939.5654133_wp,mp=938.2720813_wp
  real(wp)::n,ex,dex,g1,g2,dg1,dg2,d1,d2,e1,e2,z,kn,kp,tn,tp,cn,cp,delta,cu
  s=nuclear_t();status='APR density outside .01..2 fm^-3'
  n=nn+np
  if(n<.01_wp.or.n>2._wp)return
  status='APR branch invalid';if(branch/=1.and.branch/=2)return
  ex=exp(-(a(9)*n)**2);dex=-2*a(9)**2*n*ex
  g1=-(a(1)+a(2)*n+a(6)*n*n+(a(10)+a(11)*n)*ex)
  g2=-(a(12)/n+a(7)+a(8)*n+a(13)*ex)
  dg1=-(a(2)+2*a(6)*n+a(11)*ex+(a(10)+a(11)*n)*dex)
  dg2=a(12)/n**2-a(8)-a(13)*dex
  if(branch==2)then
   d1=n-a(19);d2=n-a(20);e1=exp(a(18)*d1);e2=exp(a(16)*d2)
   g1=g1-(a(17)*d1+a(21)*d1*d1)*e1
   g2=g2-(a(15)*d2+a(14)*d2*d2)*e2
   dg1=dg1-(a(17)+2*a(21)*d1+a(18)*(a(17)*d1+a(21)*d1*d1))*e1
   dg2=dg2-(a(15)+2*a(14)*d2+a(16)*(a(15)*d2+a(14)*d2*d2))*e2
  endif
  z=exp(-a(4)*n);kn=hc**2/(2*mn)+(a(3)*n+a(5)*nn)*z
  kp=hc**2/(2*mp)+(a(3)*n+a(5)*np)*z
  status='APR inverse effective mass not positive';if(min(kn,kp)<=0)return
  tn=.6_wp*(3*pi*pi)**(2._wp/3)*nn**(5._wp/3)
  tp=.6_wp*(3*pi*pi)**(2._wp/3)*np**(5._wp/3)
  cn=z*(a(3)-a(4)*(a(3)*n+a(5)*nn));cp=z*(a(3)-a(4)*(a(3)*n+a(5)*np))
  delta=nn-np;cu=4*dg1*nn*np+dg2*delta**2
  s%mn=mn+kn*(3*pi*pi*nn)**(2._wp/3)+(cn+a(5)*z)*tn+cp*tp+cu+4*g1*np+2*g2*delta
  s%mp=mp+kp*(3*pi*pi*np)**(2._wp/3)+cn*tn+(cp+a(5)*z)*tp+cu+4*g1*nn-2*g2*delta
  s%e=mn*nn+mp*np+kn*tn+kp*tp+4*g1*nn*np+g2*delta**2
  s%p=nn*s%mn+np*s%mp-s%e;status='ok'
 end subroutine
 subroutine skyrme_nuclear(model,nn,np,s,status)
  character(len=*),intent(in)::model
  real(wp),intent(in)::nn,np
  type(nuclear_t),intent(out)::s
  character(len=*),intent(out)::status
  real(wp)::t0,t1,t2,t3,x0,x1,x2,x3,alpha,n,tn,tp,tau,ff,a0,a3,c0,c1,e,dsq
  real(wp)::den(2),tq(2),dt(2),mass(2),chem(2)
  integer::j
  s=nuclear_t();status='ok'
  if(model=='RG_KDE0v1')then
   t0=-2553.08_wp;t1=411.70_wp;t2=-419.87_wp;t3=14603.6_wp
   x0=.6483_wp;x1=-.3472_wp;x2=-.9268_wp;x3=.9475_wp;alpha=.1673_wp
  else
   t0=-2488.91_wp;t1=486.82_wp;t2=-546.39_wp;t3=13777._wp
   x0=.834_wp;x1=-.344_wp;x2=-1._wp;x3=1.354_wp;alpha=1._wp/6
  endif
  n=nn+np;ff=(3*pi*pi)**(2._wp/3);tn=.6_wp*ff*nn**(5._wp/3);tp=.6_wp*ff*np**(5._wp/3)
  tau=tn+tp;dsq=nn*nn+np*np
  a0=(2+x0)*n*n-(2*x0+1)*dsq;a3=(2+x3)*n*n-(2*x3+1)*dsq
  c0=(t1*(2+x1)+t2*(2+x2))/8;c1=-(t1*(2*x1+1)-t2*(2*x2+1))/8
  e=20.73553_wp*tau+t0*a0/4+t3*n**alpha*a3/24+c0*n*tau+c1*(nn*tn+np*tp)
  den=[nn,np];tq=[tn,tp];dt=ff*den**(2._wp/3);mass=[939.5654_wp,938.2720_wp]
  do j=1,2
   chem(j)=mass(j)+20.73553_wp*dt(j)+t0*((2+x0)*n-(2*x0+1)*den(j))/2+&
    t3*(alpha*n**(alpha-1)*a3+2*n**alpha*((2+x3)*n-(2*x3+1)*den(j)))/24+&
    c0*(tau+n*dt(j))+c1*(tq(j)+den(j)*dt(j))
  enddo
  s%e=dot_product(mass,den)+e;s%mn=chem(1);s%mp=chem(2);s%p=dot_product(den,chem)-s%e
 end subroutine
 subroutine pchip_source(x,y,n,v,dv)
  ! Shape-preserving cubic Hermite interpolation of published ENERGY samples.
  ! Fritsch-Butland interior harmonic slopes and one-sided endpoint limits.
  real(wp),intent(in)::x(:),y(:),n
  real(wp),intent(out)::v,dv
  real(wp)::h(size(x)-1),d(size(x)-1),m(size(x)),w1,w2,t,hh
  integer::i,k,l
  l=size(x);h=x(2:)-x(:l-1);d=(y(2:)-y(:l-1))/h
  do i=2,l-1
   if(d(i-1)*d(i)<=0)then
    m(i)=0
   else
    w1=2*h(i)+h(i-1);w2=h(i)+2*h(i-1)
    m(i)=(w1+w2)/(w1/d(i-1)+w2/d(i))
   endif
  enddo
  m(1)=((2*h(1)+h(2))*d(1)-h(1)*d(2))/(h(1)+h(2))
  if(m(1)*d(1)<=0)then
   m(1)=0
  else if(d(1)*d(2)<0.and.abs(m(1))>3*abs(d(1)))then
   m(1)=3*d(1)
  endif
  m(l)=((2*h(l-1)+h(l-2))*d(l-1)-h(l-1)*d(l-2))/(h(l-1)+h(l-2))
  if(m(l)*d(l-1)<=0)then
   m(l)=0
  else if(d(l-1)*d(l-2)<0.and.abs(m(l))>3*abs(d(l-1)))then
   m(l)=3*d(l-1)
  endif
  k=l-1
  do i=1,l-1
   if(n<=x(i+1))then
    k=i;exit
   endif
  enddo
  hh=h(k);t=(n-x(k))/hh
  v=(2*t**3-3*t*t+1)*y(k)+(t**3-2*t*t+t)*hh*m(k)+&
    (-2*t**3+3*t*t)*y(k+1)+(t**3-t*t)*hh*m(k+1)
  dv=(6*t*t-6*t)*(y(k)-y(k+1))/hh+(3*t*t-4*t+1)*m(k)+(3*t*t-2*t)*m(k+1)
 end subroutine
 subroutine bbb_nuclear(nn,np,s,status)
  real(wp),intent(in)::nn,np
  type(nuclear_t),intent(out)::s
  character(len=*),intent(out)::status
  call bbb_density(nn+np,np/(nn+np),s,status)
 end subroutine
 subroutine bbb_density(n,yp,s,status)
  ! Preserve the caller's independent total density exactly at source endpoints.
  ! Reconstructing n from two rounded species densities can leave the support.
  real(wp),intent(in)::n,yp
  type(nuclear_t),intent(out)::s
  character(len=*),intent(out)::status
  real(wp),parameter::x(16)=[.08_wp,.16_wp,.2_wp,.3_wp,.4_wp,.5_wp,.6_wp,.7_wp,.8_wp,.9_wp,1._wp,1.1_wp,1.2_wp,1.3_wp,1.4_wp,1.5_wp]
  real(wp),parameter::es(16)=[-10.8_wp,-15.88_wp,-15.74_wp,-9.4_wp,3.2_wp,20.46_wp,41.78_wp,67.17_wp,&
    96.7_wp,130.22_wp,167.4_wp,207.76_wp,250.87_wp,296.27_wp,343.6_wp,392.52_wp]
  real(wp),parameter::en(16)=[9.47_wp,17._wp,22.1_wp,38.8_wp,61.4_wp,90.3_wp,125.55_wp,166.7_wp,&
    212.66_wp,262.3_wp,314.6_wp,368.78_wp,424.12_wp,480.3_wp,537._wp,594.2_wp]
  real(wp)::e0,ep,d0,dp,sym,dsym,delta,e,d,common,diff
  s=nuclear_t();status='BBB published-energy support is .08..1.5 fm^-3'
  if(n<.08_wp.or.n>1.5_wp)return
  call pchip_source(x,es,n,e0,d0);call pchip_source(x,en,n,ep,dp)
  sym=ep-e0;dsym=dp-d0;delta=1-2*yp;e=e0+sym*delta**2;d=d0+dsym*delta**2
  common=939+e+n*d;diff=4*sym*delta;s%mn=common+yp*diff;s%mp=common-(1-yp)*diff
  s%e=n*(939+e);s%p=n*n*d;status='ok'
 end subroutine
 real(wp) function scalar_density(n,m)result(rs)
  real(wp),intent(in)::n,m
  real(wp)::x,c,power,total
  integer::k
  rs=0;if(n<=0)return
  x=(3*pi*pi*n)**(1._wp/3)/m
  if(x<.1_wp)then
   c=1;power=1;total=1
   do k=1,16
    c=c*(.5_wp-real(k,wp))/real(k,wp);power=power*x*x
    total=total+c*power*3/(2*k+3)
   enddo
   rs=n*total
  else
   rs=m**3/(2*pi*pi)*(x*sqrt(1+x*x)-asinh(x))
  endif
 end function
 real(wp) function gm_residual(v,nn,np)result(r)
  real(wp),intent(in)::v,nn,np
  real(wp)::m
  m=939/hc
  r=v/11.79_wp+.002947_wp*m*v*v-.001070_wp*v**3-&
    scalar_density(nn,m-v)-scalar_density(np,m-v)
 end function
 subroutine gm_nuclear(nn,np,s,status)
  real(wp),intent(in)::nn,np
  type(nuclear_t),intent(out)::s
  character(len=*),intent(out)::status
  real(wp)::m,n,l,r,fl,fr,a,b,fa,c,fc,eff,ee,pp,mn,ep,mp,e,field,best,v
  integer::j,k
  s=nuclear_t();m=939/hc;n=nn+np;best=huge(1._wp)
  l=0;fl=gm_residual(l,nn,np);status='no positive-effective-mass GM1 scalar root'
  do j=1,64
   r=m*(1-1.e-10_wp)*real(j,wp)/64;fr=gm_residual(r,nn,np)
   if(fl*fr<=0)then
    a=l;b=r;fa=fl
    do k=1,80
     c=.5_wp*(a+b);fc=gm_residual(c,nn,np)
     if(fa*fc<=0)then
      b=c
     else
      a=c;fa=fc
     endif
     if(b-a<2.e-14_wp)exit
    enddo
    v=.5_wp*(a+b);eff=(m-v)*hc
    call core_free_fermion(nn,eff,ee,pp,mn)
    call core_free_fermion(np,eff,ep,pp,mp)
    field=v*v/(2*11.79_wp)+.002947_wp*m*v**3/3-.001070_wp*v**4/4+&
     7.149_wp*n*n/2+4.411_wp*(np-nn)**2/8
    e=ee+ep+hc*field
    if(e<best)then
     best=e;s%e=e
     s%mn=mn+hc*(7.149_wp*n-4.411_wp*(np-nn)/4)
     s%mp=mp+hc*(7.149_wp*n+4.411_wp*(np-nn)/4)
     s%p=nn*s%mn+np*s%mp-e;status='ok'
    endif
   endif
   l=r;fl=fr
  enddo
 end subroutine
 subroutine nuclear(model,nn,np,branch,s,status)
  character(len=*),intent(in)::model
  real(wp),intent(in)::nn,np
  integer,intent(in)::branch
  type(nuclear_t),intent(out)::s
  character(len=*),intent(out)::status
  s=nuclear_t();status='invalid species densities'
  if(.not.all(ieee_is_finite([nn,np])))return
  if(min(nn,np)<0.or.nn+np<=0)return
  select case(trim(model))
  case('APR_APR');call apr_nuclear(nn,np,branch,s,status)
  case('BBB_BHF2');call bbb_nuclear(nn,np,s,status)
  case('GM_GM1');call gm_nuclear(nn,np,s,status)
  case('RG_KDE0v1','RG_SLY4');call skyrme_nuclear(model,nn,np,s,status)
  case default;status='unknown native core model: '//trim(model)
  end select
  if(status=='ok'.and..not.all(ieee_is_finite([s%e,s%p,s%mn,s%mp])))status='nonfinite nuclear state'
 end subroutine
 subroutine core_nuclear_state(model,nn,np,branch,epsilon,pressure,mun,mup,status)
  character(len=*),intent(in)::model
  real(wp),intent(in)::nn,np
  integer,intent(in)::branch
  real(wp),intent(out)::epsilon,pressure,mun,mup
  character(len=*),intent(out)::status
  type(nuclear_t)::s
  call nuclear(model,nn,np,branch,s,status)
  epsilon=s%e;pressure=s%p;mun=s%mn;mup=s%mp
 end subroutine
 subroutine beta_eval(model,n,yp,branch,s,r,status)
  character(len=*),intent(in)::model
  real(wp),intent(in)::n,yp
  integer,intent(in)::branch
  type(core_state_t),intent(out)::s
  real(wp),intent(out)::r
  character(len=*),intent(out)::status
  type(nuclear_t)::b
  real(wp)::ne,nm,e,p
  s=core_state_t();r=huge(1._wp)
  if(model=='BBB_BHF2')then
   call bbb_density(n,yp,b,status)
  else
   call nuclear(model,n*(1-yp),n*yp,branch,b,status)
  endif
  if(status/='ok')return
  s%mu_e=b%mn-b%mp
  call leptons(s%mu_e,model/='GM_GM1',ne,nm,e,p)
  r=n*yp-ne-nm;s%n_b=n;s%epsilon=b%e+e;s%pressure=b%p+p;s%mu_n=b%mn;s%mu_p=b%mp
  s%yn=1-yp;s%yp=yp;s%ye=ne/n;s%ymu=nm/n;s%charge_residual=r
  s%euler_residual=s%epsilon+s%pressure-n*s%mu_n
  if(model=='APR_APR'.and.branch==2)then
   s%phase=7;s%chi=1
  endif
 end subroutine
 subroutine solve_beta(model,n,branch,s,status)
  character(len=*),intent(in)::model
  real(wp),intent(in)::n
  integer,intent(in)::branch
  type(core_state_t),intent(out)::s
  character(len=*),intent(out)::status
  type(core_state_t)::trial,tmp
  real(wp)::l,r,fl,fr,a,b,fa,c,fc,best
  integer::j,k,roots
  s=core_state_t();best=huge(1._wp);roots=0
  l=0;call beta_eval(model,n,l,branch,tmp,fl,status);if(status/='ok')return
  do j=1,64
   r=.5_wp*real(j,wp)/64
   call beta_eval(model,n,r,branch,tmp,fr,status);if(status/='ok')return
   if(fl*fr<=0)then
    a=l;b=r;fa=fl
    do k=1,80
     c=.5_wp*(a+b);call beta_eval(model,n,c,branch,trial,fc,status)
     if(status/='ok')return
     if(fa*fc<=0)then
      b=c
     else
      a=c;fa=fc
     endif
     if(b-a<4.e-16_wp)exit
    enddo
    call beta_eval(model,n,.5_wp*(a+b),branch,trial,fc,status)
    if(status/='ok')return
    roots=roots+1
    if(trial%epsilon<best)then
     best=trial%epsilon;s=trial
    endif
   endif
   l=r;fl=fr
  enddo
  status='ok';if(roots==0)status='no neutral beta-equilibrium root on 0<=Yp<=1/2'
 end subroutine
 subroutine gibbs_residual(z,chi,r,low,high,ne,nm,el,pl,status)
  real(wp),intent(in)::z(5),chi
  real(wp),intent(out)::r(5),ne,nm,el,pl
  type(nuclear_t),intent(out)::low,high
  character(len=*),intent(out)::status
  call apr_nuclear(z(1)*(1-z(2)),z(1)*z(2),1,low,status);if(status/='ok')return
  call apr_nuclear(z(3)*(1-z(4)),z(3)*z(4),2,high,status);if(status/='ok')return
  call leptons(z(5),.true.,ne,nm,el,pl)
  r=[low%mn-high%mn,low%mp-high%mp,(low%p-high%p)/.2_wp,low%mn-low%mp-z(5),&
    ((1-chi)*z(1)*z(2)+chi*z(3)*z(4)-ne-nm)*5000]
 end subroutine
 subroutine linear5(a,b,x,ok)
  real(wp),intent(in)::a(5,5),b(5)
  real(wp),intent(out)::x(5)
  logical,intent(out)::ok
  real(wp)::w(5,6),row(6),factor
  integer::i,j,k,p
  w(:,:5)=a;w(:,6)=b;ok=.false.;x=0
  do k=1,5
   p=k
   do i=k+1,5
    if(abs(w(i,k))>abs(w(p,k)))p=i
   enddo
   if(abs(w(p,k))<1.e-14_wp)return
   if(p/=k)then
    row=w(k,:);w(k,:)=w(p,:);w(p,:)=row
   endif
   do i=k+1,5
    factor=w(i,k)/w(k,k);w(i,k:)=w(i,k:)-factor*w(k,k:)
   enddo
  enddo
  do i=5,1,-1
   x(i)=(w(i,6)-sum(w(i,i+1:5)*x(i+1:5)))/w(i,i)
  enddo
  ok=all(ieee_is_finite(x))
 end subroutine
 subroutine gibbs_root(chi,seed,s,status)
  real(wp),intent(in)::chi,seed(5)
  type(core_state_t),intent(out)::s
  character(len=*),intent(out)::status
  real(wp),parameter::scale(5)=[.2_wp,.1_wp,.25_wp,.1_wp,140._wp]
  real(wp),parameter::lo(5)=[.10_wp,.001_wp,.20_wp,.001_wp,100._wp],hi(5)=[.25_wp,.30_wp,.36_wp,.30_wp,180._wp]
  real(wp)::z(5),zp(5),zm(5),trial(5),r(5),rp(5),rm(5),rt(5),jac(5,5),dx(5),h,alpha
  real(wp)::ne,nm,el,pl,npt,oldnorm
  type(nuclear_t)::low,high,t1,t2
  integer::iter,j,line
  logical::ok,accepted
  s=core_state_t();z=seed;status='APR Gibbs root did not converge'
  do iter=1,80
   call gibbs_residual(z,chi,r,low,high,ne,nm,el,pl,status);if(status/='ok')return
   if(maxval(abs(r))<2.e-10_wp)exit
   do j=1,5
    h=2.e-6_wp*scale(j);zp=z;zm=z;zp(j)=z(j)+h;zm(j)=z(j)-h
    call gibbs_residual(zp,chi,rp,t1,t2,ne,nm,el,pl,status);if(status/='ok')return
    call gibbs_residual(zm,chi,rm,t1,t2,ne,nm,el,pl,status);if(status/='ok')return
    jac(:,j)=(rp-rm)/(2*h)*scale(j)
   enddo
   call linear5(jac,-r,dx,ok)
   if(.not.ok)then
    status='singular APR Gibbs Jacobian';return
   endif
   alpha=1;accepted=.false.;oldnorm=sum(r*r)
   do line=1,40
    trial=z+alpha*dx*scale
    if(all(trial>lo).and.all(trial<hi))then
     call gibbs_residual(trial,chi,rt,t1,t2,ne,nm,el,pl,status)
     if(status=='ok')then
      if(sum(rt*rt)<oldnorm)then
       accepted=.true.;exit
      endif
     endif
    endif
    alpha=alpha*.5_wp
   enddo
   if(.not.accepted)then
    status='APR Gibbs Newton line search exhausted';return
   endif
   z=trial
  enddo
  call gibbs_residual(z,chi,r,low,high,ne,nm,el,pl,status);if(status/='ok')return
  if(maxval(abs(r))>1.e-7_wp.or.z(3)<=z(1))then
   status='APR Gibbs residual or density ordering failed';return
  endif
  s%n_b=(1-chi)*z(1)+chi*z(3);npt=(1-chi)*z(1)*z(2)+chi*z(3)*z(4)
  s%epsilon=(1-chi)*low%e+chi*high%e+el;s%pressure=.5_wp*(low%p+high%p)+pl
  s%mu_n=.5_wp*(low%mn+high%mn);s%mu_e=z(5);s%mu_p=s%mu_n-s%mu_e
  s%yn=1-npt/s%n_b;s%yp=npt/s%n_b;s%ye=ne/s%n_b;s%ymu=nm/s%n_b
  s%chi=chi;s%n_low=z(1);s%n_high=z(3);s%yp_low=z(2);s%yp_high=z(4);s%phase=8
  s%gibbs_residual=maxval(abs(r));s%charge_residual=npt-ne-nm
  s%euler_residual=s%epsilon+s%pressure-s%n_b*s%mu_n;status='ok'
 end subroutine
 subroutine solve_apr_gibbs_fraction(chi,s,status)
  real(wp),intent(in)::chi
  type(core_state_t),intent(out)::s
  character(len=*),intent(out)::status
  s=core_state_t();status='APR phase fraction must be in [0,1]'
  if(.not.ieee_is_finite(chi))return
  if(chi<0.or.chi>1)return
  call gibbs_root(chi,[.205_wp,.0735_wp,.237_wp,.0588_wp,138._wp],s,status)
 end subroutine
 subroutine apr_endpoints(left,right,status)
  type(core_state_t),intent(out)::left,right
  character(len=*),intent(out)::status
  type(core_state_t)::check,s
  real(wp)::chi
  integer::i
  do i=0,1
   chi=real(i,wp)
   call solve_apr_gibbs_fraction(chi,s,status);if(status/='ok')return
   call gibbs_root(chi,[.19_wp,.06_wp,.26_wp,.07_wp,140._wp],check,status);if(status/='ok')return
   if(maxval(abs([s%n_low-check%n_low,s%n_high-check%n_high,s%yp_low-check%yp_low,s%yp_high-check%yp_high]))>1.e-8_wp)then
    status='APR Gibbs multistart roots disagree';return
   endif
   call gibbs_root(chi,[.22_wp,.10_wp,.28_wp,.09_wp,145._wp],check,status);if(status/='ok')return
   if(maxval(abs([s%n_low-check%n_low,s%n_high-check%n_high,s%yp_low-check%yp_low,s%yp_high-check%yp_high]))>1.e-8_wp)then
    status='APR Gibbs multistart roots disagree';return
   endif
   if(i==0)then
    left=s
   else
    right=s
   endif
  enddo
  status='ok'
  if(right%n_b<=left%n_b)status='APR mean coexistence densities not ordered'
 end subroutine
 subroutine solve_with_context(model,n,s,status,left,right)
  character(len=*),intent(in)::model
  real(wp),intent(in)::n
  type(core_state_t),intent(out)::s
  character(len=*),intent(out)::status
  type(core_state_t),intent(in)::left,right
  real(wp)::a,b,chi
  integer::k
  s=core_state_t();status='density must be finite and positive'
  if(.not.ieee_is_finite(n))return
  if(n<=0)return
  if(model/='APR_APR')then
   call solve_beta(model,n,1,s,status)
  else if(n<=left%n_b)then
   call solve_beta(model,n,1,s,status)
  else if(n>=right%n_b)then
   call solve_beta(model,n,2,s,status)
  else
   a=0;b=1
   do k=1,60
    chi=.5_wp*(a+b);call solve_apr_gibbs_fraction(chi,s,status);if(status/='ok')return
    if(abs(s%n_b-n)<1.e-14_wp)exit
    if(s%n_b>n)then
     b=chi
    else
     a=chi
    endif
   enddo
   if(abs(s%n_b-n)>1.e-11_wp)then
    status='APR Gibbs mean-density inversion failed';return
   endif
  endif
  if(status/='ok')return
  if(abs(s%yp-s%ye-s%ymu)>1.e-10_wp.or.abs(s%epsilon+s%pressure-s%n_b*s%mu_n)>1.e-10_wp*s%n_b*939)then
   status='native core conservation/Euler acceptance failed'
  endif
 end subroutine
 subroutine solve_core_state(model,n,s,status)
  character(len=*),intent(in)::model
  real(wp),intent(in)::n
  type(core_state_t),intent(out)::s
  character(len=*),intent(out)::status
  type(core_state_t)::left,right
  real(wp)::lo,hi,mn,mp
  s=core_state_t()
  call core_model_bounds(model,lo,hi,mn,mp,status);if(status/='ok')return
  left=core_state_t();right=core_state_t()
  if(model=='APR_APR')then
   call apr_endpoints(left,right,status);if(status/='ok')return
  endif
  call solve_with_context(model,n,s,status,left,right)
 end subroutine
 subroutine generate_core_native(model,densities,eos,states,status)
  character(len=*),intent(in)::model
  real(wp),intent(in)::densities(:)
  type(eos_data_t),intent(out)::eos
  type(core_state_t),allocatable,intent(out)::states(:)
  character(len=*),intent(out)::status
  type(core_state_t)::left,right
  real(wp),allocatable::grid(:)
  real(wp)::lo,hi,mn,mp,v,edges(2)
  integer::n,i,j,k
  call core_model_bounds(model,lo,hi,mn,mp,status);if(status/='ok')return
  status='at least three finite positive increasing densities required'
  n=size(densities);if(n<3)return
  if(.not.all(ieee_is_finite(densities)))return
  if(any(densities<=0).or.any(densities(2:)<=densities(:n-1)))return
  allocate(grid(n+2));grid(:n)=densities
  left=core_state_t();right=core_state_t()
  if(model=='APR_APR')then
   call apr_endpoints(left,right,status);if(status/='ok')return
   edges=[left%n_b,right%n_b]
   do k=1,2
    if(edges(k)>densities(1).and.edges(k)<densities(size(densities)))then
     if(.not.any(grid(:n)==edges(k)))then
      n=n+1;grid(n)=edges(k)
     endif
    endif
   enddo
  endif
  do i=2,n
   v=grid(i);j=i-1
   do while(j>=1)
    if(grid(j)<=v)exit
    grid(j+1)=grid(j);j=j-1
   enddo
   grid(j+1)=v
  enddo
  call allocate_eos(eos,n);allocate(states(n));eos%mn=mn;eos%mp=mp
  do i=1,n
   call solve_with_context(model,grid(i),states(i),status,left,right)
   if(status/='ok')return
   eos%nb(i)=states(i)%n_b;eos%q(1,i)=states(i)%pressure/eos%nb(i)
   eos%q(3,i)=states(i)%mu_n/mn-1;eos%q(4,i)=-states(i)%mu_e/mn
   eos%q(6,i)=states(i)%epsilon/(eos%nb(i)*mn)-1;eos%q(7,i)=eos%q(6,i)
   eos%yn(i)=states(i)%yn;eos%yp(i)=states(i)%yp;eos%ye(i)=states(i)%ye;eos%ymu(i)=states(i)%ymu
   eos%iphase(i)=states(i)%phase
  enddo
  eos%abar=1;eos%zbar=0;eos%yh=0
  eos%composition_file_present=.true.;eos%composition_available=.true.;eos%beta_equilibrated=.true.
  eos%yn_present=.true.;eos%yp_present=.true.;eos%ye_present=.true.;eos%ymu_present=.true.
  eos%abar_present=.true.;eos%zbar_present=.true.;eos%yh_present=.true.
  eos%c3_applicable=.true.;eos%c4_applicable=.true.;eos%row_status=2
  status='ok'
 end subroutine
end module
