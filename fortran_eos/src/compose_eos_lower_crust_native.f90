! Native lower/outer-crust construction. Units MeV, fm; cold locally neutral matter.
! Physics sources: Pearson et al., MNRAS481 2994 (2018), Sec2.1/AppB;
! corrected MNRAS486 768 (2019), AppC; author bskfit18.f (2023-02-13).
! BSk API is an analytic representation, NOT a new ETF/Strutinsky minimization.
module compose_eos_lower_crust_native
 use compose_eos_types, only: wp,eos_data_t,allocate_eos
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 real(wp),parameter :: PI=acos(-1._wp),HBARC=197.3269804_wp,ME=.51099895_wp
 real(wp),parameter :: MN=939.56542052_wp,MP=938.27208816_wp,AU=931.49410242_wp
 real(wp),parameter :: ALPHA=1._wp/137.035999084_wp,E2=ALPHA*HBARC,BCC=.895929255682_wp,U1=.5113875_wp
 real(wp),parameter,public :: BSK_DRIP=2.69403e-4_wp,BSK_LOWER_MAX=.002_wp,BSK_EXT_MAX=.01_wp
 real(wp),parameter :: BSK_CORE=.0716068_wp,BSK_MN=939.5654_wp,BSK_E0=-9.1536_wp,RHO_FACTOR=1.782662e12_wp
 real(wp),parameter :: EP(14)=[7.02e8_wp,1.133e11_wp,6.19e7_wp,4.54e6_wp,5.46e5_wp,15.24_wp,.0683_wp,8.86_wp,4611._wp,48.07_wp,2.697_wp,81.7_wp,7.05_wp,1.50_wp]
 real(wp),parameter :: PP(23)=[6.682_wp,5.651_wp,.00459_wp,.14359_wp,2.681_wp,11.972_wp,13.993_wp,1.2904_wp,2.665_wp,-27.787_wp,2.0140_wp,4.09_wp,14.135_wp,28.03_wp,-1.921_wp,1.08_wp,14.89_wp,.098_wp,11.67_wp,4.75_wp,-.037_wp,14.10_wp,11.9_wp]
 type,public :: mass_catalogue_t
  integer :: n=0
  real(wp),allocatable :: z(:),a(:),excess(:),radius(:)
  character(len=64),allocatable :: source(:)
 end type
 type,public :: lower_crust_state_t
  real(wp) :: nb=0,ne=0,epsilon=0,pressure=0,mu_b=0,mu_e=0
  real(wp) :: yn=0,yp=0,ye=0,ymu=0,yh=0,aav=0,zav=0
  real(wp) :: pressure_residual=0,second_best_gap=0
  real(wp) :: source_epsilon=0,source_pressure=0,energy_derivative_pressure=0
  real(wp) :: mu_composition=0,composition_mu_residual=0,source_thermo_residual=0
  real(wp) :: sound2=0,gamma=0,quadrature_error=0
  integer :: mass_index=0,phase=0
  character(len=64) :: representation='',mass_source=''
 end type
 ! Independent source branches. Differences are lower-first minus outer-last.
 type,public :: lower_source_seam_t
  real(wp) :: outer_nb=0,lower_nb=0,delta_nb=0,delta_pressure=0,delta_mu_b=0,delta_epsilon=0
  logical :: density_overlap=.false.,pressure_continuous=.false.,mu_continuous=.false.
 end type
 type,public :: bsk_anchor_t
  real(wp) :: outer_nb=0,outer_p=0,outer_mu=0,outer_epsilon=0,n0=0,e0=0,epsabs=1e-10_wp
  logical :: valid=.false.
 end type
 public :: load_mass_catalogue,outer_isotope_state,outer_state,outer_neutron_drip,outer_sequence
 public :: bsk22_state,bsk22_make_anchor,bsk22_anchored_state,bsk22_generate,lower_states_to_eos
 public :: prepare_lower_source_native
contains
 subroutine load_mass_catalogue(path,cat,status)
  character(len=*),intent(in)::path
  type(mass_catalogue_t),intent(out)::cat
  character(len=*),intent(out)::status
  character(len=1024)::line
  integer::u,ios,k,j,count
  status='invalid nuclear mass catalogue'
  open(newunit=u,file=path,status='old',action='read',iostat=ios)
  if(ios/=0)then;status='cannot open nuclear mass CSV';return;endif
  read(u,'(a)',iostat=ios)line
  if(ios/=0)then;close(u);return;endif
  count=0
  do
   read(u,'(a)',iostat=ios)line
   if(ios<0)exit
   if(ios/=0)then;close(u);return;endif
   if(len_trim(line)>0)count=count+1
  enddo
  if(count<1)then;close(u);return;endif
  rewind(u);read(u,'(a)')line
  allocate(cat%z(count),cat%a(count),cat%excess(count),cat%radius(count),cat%source(count))
  k=0
  do
   read(u,'(a)',iostat=ios)line
   if(ios<0)exit
   if(len_trim(line)==0)cycle
   k=k+1
   read(line,*,iostat=ios)cat%z(k),cat%a(k),cat%excess(k),cat%radius(k),cat%source(k)
   if(ios/=0)then;close(u);return;endif
  enddo
  close(u)
  if(.not.all(ieee_is_finite(cat%z)).or..not.all(ieee_is_finite(cat%a)))return
  if(.not.all(ieee_is_finite(cat%excess)).or..not.all(ieee_is_finite(cat%radius)))return
  if(any(cat%z<=0).or.any(cat%a<cat%z).or.any(cat%radius<0))return
  if(any(cat%z/=anint(cat%z)).or.any(cat%a/=anint(cat%a)))return
  do k=1,count
   do j=1,k-1
    if(cat%z(k)==cat%z(j).and.cat%a(k)==cat%a(j))then
     status='duplicate nuclear mass entry';return
    endif
   enddo
  enddo
  cat%n=count;status='ok'
 end subroutine

 pure subroutine isotope_ep(ne,z,a,excess,radius,energy,pressure,mu_e)
  real(wp),intent(in)::ne,z,a,excess,radius
  real(wp),intent(out)::energy,pressure,mu_e
  real(wp)::x,gam,iv,fac,ee,pe,aa,k,ex,px,cp,corr,pcorr,atom,rc,hw,pl,s,b1,b2,b3,b4,finf,sa,bb,base,screen,df,da,ds,pscreen,lat,fs,zp,eel,pel
  x=HBARC*(3*PI**2*ne)**(1._wp/3)/ME;gam=sqrt(1+x*x);iv=asinh(x)
  fac=ME**4/(8*PI**2*HBARC**3)
  ee=fac*(x*(1+2*x*x)*gam-iv);pe=fac/3*(x*(2*x*x-3)*gam+3*iv)
  aa=x*gam-iv;k=ALPHA*ME**4/(4*PI**3*HBARC**3)
  ex=k*(1.5_wp*aa*aa-x**4);px=k*(2*aa*x**3/gam-x**4/3-1.5_wp*aa*aa)
  cp=ALPHA**2*ME**4/(12*PI**4*HBARC**3);corr=cp*x**4*(-12.51_wp+log10(x))
  pcorr=corr/3+cp*x**4/(3*log(10._wp))
  atom=a*AU+excess+1.44381e-5_wp*z**2.39_wp+1.55468e-12_wp*z**5.35_wp
  rc=(3*z/(4*PI*ne))**(1._wp/3);hw=HBARC*sqrt(4*PI*z*E2*ne/atom);pl=z*z*E2/(rc*hw)
  s=1/(1+.01_wp*log(z)**1.5_wp+.097_wp/z**2)
  b1=1-1.1866_wp/z**.267_wp+.27_wp/z
  b2=1+2.25_wp/z**(1._wp/3)*(1+.684_wp*z**5+.222_wp*z**6)/(1+.222_wp*z**6)
  b3=41.5_wp/(1+log(z));b4=.395_wp*log(z)+.347_wp/z**1.5_wp
  finf=54._wp/175*(12/PI)**(1._wp/3)*ALPHA*z**(2._wp/3)*b1*sqrt(1+b2/x**2)
  sa=(b3+17.9_wp*x*x)/(1+b4*x*x);bb=(.205_wp/pl)**s
  base=(4*PI/3)**(1._wp/3)*E2*z**(2._wp/3)*ne**(4._wp/3)
  screen=-base*finf*(1+sa*bb);df=-finf*b2/(3*(x*x+b2))
  da=2*x*x*(17.9_wp-b4*b3)/(3*(1+b4*x*x)**2)
  ds=df*(1+sa*bb)+finf*bb*(da+sa*s/6);pscreen=screen/3-base*ds
  lat=-BCC*base;fs=-lat*5*radius**2/(9*rc**2);zp=ne/z*1.5_wp*U1*hw
  eel=ee+ex+corr+screen;pel=pe+px+pcorr+pscreen
  energy=ne/z*(atom-z*ME)+eel+lat+fs+zp
  pressure=pel+lat/3+fs+zp/2;mu_e=(eel+pel)/ne
 end subroutine

 subroutine outer_isotope_state(ne,z,a,excess,radius,state,status)
  real(wp),intent(in)::ne,z,a,excess,radius
  type(lower_crust_state_t),intent(out)::state
  character(len=*),intent(out)::status
  real(wp)::x
  state=lower_crust_state_t();status='invalid isotope inputs'
  if(.not.all(ieee_is_finite([ne,z,a,excess,radius])))return
  if(ne<=0.or.z<=0.or.a<z.or.radius<0.or.z/=anint(z).or.a/=anint(a))return
  x=HBARC*(3*PI**2*ne)**(1._wp/3)/ME
  if(x<.1_wp.or.x>100)then;status='outer isotope outside 0.1<=x_e<=100 support';return;endif
  call isotope_ep(ne,z,a,excess,radius,state%epsilon,state%pressure,state%mu_e)
  state%nb=ne*a/z;state%ne=ne;state%mu_b=(state%epsilon+state%pressure)/state%nb
  state%ye=z/a;state%yh=1;state%aav=a;state%zav=z
  state%representation='discrete-nucleus fixed-pressure Gibbs minimum'
  status='ok'
 end subroutine

 subroutine isotope_pressure_root(cat,i,p,ne,e,pressure,mu,g,status)
  type(mass_catalogue_t),intent(in)::cat
  integer,intent(in)::i
  real(wp),intent(in)::p
  real(wp),intent(out)::ne,e,pressure,mu,g
  character(len=*),intent(out)::status
  real(wp)::lo,hi,mid,pm,phi
  integer::j
  lo=1e-11_wp;hi=2e-3_wp;status='ok'
  call isotope_ep(lo,cat%z(i),cat%a(i),cat%excess(i),cat%radius(i),e,pm,mu)
  call isotope_ep(hi,cat%z(i),cat%a(i),cat%excess(i),cat%radius(i),e,phi,mu)
  if(pm>=p.or.phi<=p)then;status='outer pressure bracket does not cover isotope';return;endif
  do j=1,44
   mid=sqrt(lo*hi)
   call isotope_ep(mid,cat%z(i),cat%a(i),cat%excess(i),cat%radius(i),e,pm,mu)
   if(pm<p)then;lo=mid;else;hi=mid;endif
  enddo
  ne=sqrt(lo*hi)
  call isotope_ep(ne,cat%z(i),cat%a(i),cat%excess(i),cat%radius(i),e,pressure,mu)
  g=(e+pressure)/(ne*cat%a(i)/cat%z(i))
 end subroutine

 subroutine outer_candidates(cat,p,nes,energies,pressures,mus,gs,status)
  type(mass_catalogue_t),intent(in)::cat
  real(wp),intent(in)::p
  real(wp),intent(out)::nes(:),energies(:),pressures(:),mus(:),gs(:)
  character(len=*),intent(out)::status
  integer::i
  status='invalid outer pressure or mass catalogue'
  if(.not.ieee_is_finite(p).or.p<1e-13_wp.or.p>1e-2_wp.or.cat%n<1)return
  do i=1,cat%n
   call isotope_pressure_root(cat,i,p,nes(i),energies(i),pressures(i),mus(i),gs(i),status)
   if(status/='ok')return
  enddo
  if(.not.all(ieee_is_finite(gs)))then;status='nonfinite isotope Gibbs energy';return;endif
  status='ok'
 end subroutine

 subroutine outer_state(cat,p,state,status,allow_above_drip)
  type(mass_catalogue_t),intent(in)::cat
  real(wp),intent(in)::p
  type(lower_crust_state_t),intent(out)::state
  character(len=*),intent(out)::status
  logical,intent(in),optional::allow_above_drip
  real(wp),allocatable::nes(:),es(:),ps(:),mus(:),gs(:)
  integer::i
  logical::above
  state=lower_crust_state_t();above=.false.
  if(present(allow_above_drip))above=allow_above_drip
  allocate(nes(cat%n),es(cat%n),ps(cat%n),mus(cat%n),gs(cat%n))
  call outer_candidates(cat,p,nes,es,ps,mus,gs,status)
  if(status/='ok')return
  i=minloc(gs,dim=1)
  if(cat%radius(i)==0)then;status='selected isotope lacks source charge radius';return;endif
  if(gs(i)>MN+1e-9_wp.and..not.above)then;status='pressure is above neutron drip';return;endif
  call outer_isotope_state(nes(i),cat%z(i),cat%a(i),cat%excess(i),cat%radius(i),state,status)
  if(status/='ok')return
  state%pressure=p;state%mu_b=(state%epsilon+p)/state%nb
  state%pressure_residual=abs(ps(i)/p-1);state%mass_index=i;state%mass_source=cat%source(i)
  gs(i)=huge(1._wp);state%second_best_gap=minval(gs)-(es(i)+ps(i))/(nes(i)*cat%a(i)/cat%z(i))
 end subroutine

 subroutine outer_neutron_drip(cat,state,status)
  type(mass_catalogue_t),intent(in)::cat
  type(lower_crust_state_t),intent(out)::state
  character(len=*),intent(out)::status
  real(wp)::lo,hi,mid,flo,fhi
  integer::i
  lo=log(1e-5_wp);hi=log(2e-3_wp)
  call outer_state(cat,exp(lo),state,status,.true.);if(status/='ok')return
  flo=state%mu_b-MN
  call outer_state(cat,exp(hi),state,status,.true.);if(status/='ok')return
  fhi=state%mu_b-MN
  if(flo>=0.or.fhi<=0)then;status='no neutron drip bracket';return;endif
  do i=1,64
   mid=(lo+hi)/2
   call outer_state(cat,exp(mid),state,status,.true.);if(status/='ok')return
   if(state%mu_b<MN)then;lo=mid;else;hi=mid;endif
   if(hi-lo<2e-13_wp)exit
  enddo
  call outer_state(cat,exp((lo+hi)/2),state,status,.true.)
 end subroutine

 subroutine outer_sequence(cat,points,eos,states,status,p_min,refinements)
  type(mass_catalogue_t),intent(in)::cat
  integer,intent(in)::points
  type(eos_data_t),intent(out)::eos
  type(lower_crust_state_t),allocatable,intent(out)::states(:)
  character(len=*),intent(out)::status
  real(wp),intent(in),optional::p_min
  integer,intent(out),optional::refinements
  type(lower_crust_state_t),allocatable::rows(:),sampled(:)
  type(lower_crust_state_t)::drip,tmp
  real(wp)::p0,lp0,lp1
  integer::i,j,k,nr,nclean,nref
  status='outer sequence requires at least 16 samples'
  if(points<16)return
  call outer_neutron_drip(cat,drip,status);if(status/='ok')return
  if(present(p_min))then
   p0=p_min
  else
   k=0
   do i=1,cat%n
    if(cat%z(i)==26.and.cat%a(i)==56)k=i
   enddo
   if(k==0)then;status='default lower boundary requires Fe56';return;endif
   call outer_isotope_state(1e-9_wp*26/56,26._wp,56._wp,cat%excess(k),cat%radius(k),tmp,status)
   if(status/='ok')return
   p0=tmp%pressure
  endif
  if(.not.ieee_is_finite(p0).or.p0<1e-13_wp.or.p0>=drip%pressure)then;status='invalid outer lower pressure';return;endif
  allocate(rows(points+2*cat%n),sampled(points));nr=points;nref=0
  lp0=log(p0);lp1=log(drip%pressure)
  do i=1,points
   call outer_state(cat,exp(lp0+(lp1-lp0)*real(i-1,wp)/(points-1)),rows(i),status,.true.)
   if(status/='ok')return
  enddo
  sampled=rows(:points)
  do i=1,points-1
   call refine(sampled(i),sampled(i+1),0)
   if(status/='ok')return
  enddo
  ! Stable insertion sort retains both one-sided phase limits.
  do i=2,nr
   tmp=rows(i);j=i-1
   do while(j>=1)
    if(rows(j)%nb<=tmp%nb)exit
    rows(j+1)=rows(j);j=j-1
   enddo
   rows(j+1)=tmp
  enddo
  nclean=1
  do i=2,nr
   if(rows(i)%nb>rows(nclean)%nb*(1+1e-13_wp))then
    nclean=nclean+1;rows(nclean)=rows(i)
   endif
  enddo
  allocate(states(nclean));states=rows(:nclean)
  do i=1,nclean
   states(i)%phase=10
   do j=1,i-1
    if(states(j)%mass_index==states(i)%mass_index)then
     states(i)%phase=states(j)%phase;exit
    endif
    states(i)%phase=max(states(i)%phase,states(j)%phase+1)
   enddo
  enddo
  call lower_states_to_eos(states,MN,MP,eos,status)
  if(present(refinements))refinements=nref
 contains
  recursive subroutine refine(left,right,depth)
   type(lower_crust_state_t),intent(in)::left,right
   integer,intent(in)::depth
   type(lower_crust_state_t)::middle,pair
   real(wp),allocatable::nes(:),es(:),ps(:),mus(:),gs(:)
   real(wp)::lo,hi,mid,p,delta,ni,nj,ei,ej,pii,pjj,mui,muj,gi,gj
   integer::i,j,k,iteration,winner,q
   i=left%mass_index;j=right%mass_index
   if(i==j)return
   if(depth>32)then;status='nuclear transition depth exceeded';return;endif
   allocate(nes(cat%n),es(cat%n),ps(cat%n),mus(cat%n),gs(cat%n))
   lo=log(left%pressure);hi=log(right%pressure)
   do iteration=1,64
    mid=(lo+hi)/2;p=exp(mid)
    call isotope_pressure_root(cat,i,p,ni,ei,pii,mui,gi,status);if(status/='ok')return
    call isotope_pressure_root(cat,j,p,nj,ej,pjj,muj,gj,status);if(status/='ok')return
    delta=gi-gj
    if(delta<0)then;lo=mid;else;hi=mid;endif
    if(hi-lo<2e-12_wp)exit
   enddo
   p=exp((lo+hi)/2)
   call outer_candidates(cat,p,nes,es,ps,mus,gs,status);if(status/='ok')return
   winner=minloc(gs,dim=1)
   if(gs(winner)<min(gs(i),gs(j))-1e-8_wp)then
    nref=nref+1
    call outer_state(cat,p,middle,status,.true.);if(status/='ok')return
    call refine(left,middle,depth+1);if(status/='ok')return
    call refine(middle,right,depth+1);return
   endif
   do q=1,2
    k=i;if(q==2)k=j
    call outer_isotope_state(nes(k),cat%z(k),cat%a(k),cat%excess(k),cat%radius(k),pair,status)
    if(status/='ok')return
    pair%pressure=p;pair%mu_b=(pair%epsilon+p)/pair%nb
    pair%mass_index=k;pair%mass_source=cat%source(k);pair%pressure_residual=abs(ps(k)/p-1)
    if(nr>=size(rows))then;status='nuclear sequence storage exhausted';return;endif
    nr=nr+1;rows(nr)=pair
   enddo
  end subroutine
 end subroutine

 pure subroutine bsk_energy(n,e,de,dde)
  real(wp),intent(in)::n
  real(wp),intent(out)::e,de,dde
  real(wp)::u2,u3,u4,u5,low,dl,mid,dm,high,dh,w1,w2,t(3),ddl,ddm,ddh,l1,l2,l3,c1,c2
  u2=sqrt(EP(2)*n);u3=sqrt(EP(3)*n);u4=sqrt(EP(4)*n);u5=sqrt(EP(5)*n)
  low=(EP(1)*n)**(7._wp/6)*(1+u4)/((1+u2)*(1+u3)*(1+u5))
  dl=7._wp/6+.5_wp*u4/(1+u4)-.5_wp*(u2/(1+u2)+u3/(1+u3)+u5/(1+u5))
  mid=EP(6)*n**EP(7)*(1+EP(8)*n);dm=EP(7)+EP(8)*n/(1+EP(8)*n)
  high=(EP(10)*n)**EP(11)/(1+EP(12)*n);dh=EP(11)-EP(12)*n/(1+EP(12)*n)
  w1=1/(1+EP(9)*n);w2=1/(1+(EP(13)*n)**EP(14))
  t=[low*w1,mid*(1-w1)*w2,high*(1-w2)]
  e=BSK_E0+sum(t)
  de=t(1)*(dl-(1-w1))+t(2)*(dm+w1-EP(14)*(1-w2))+t(3)*(dh+EP(14)*w2)
  ddl=.25_wp*(u4/(1+u4)**2-u2/(1+u2)**2-u3/(1+u3)**2-u5/(1+u5)**2)
  ddm=EP(8)*n/(1+EP(8)*n)**2;ddh=-EP(12)*n/(1+EP(12)*n)**2
  l1=dl-(1-w1);l2=dm+w1-EP(14)*(1-w2);l3=dh+EP(14)*w2
  c1=w1*(1-w1);c2=EP(14)**2*w2*(1-w2)
  dde=t(1)*(l1*l1+ddl-c1)+t(2)*(l2*l2+ddm-c1-c2)+t(3)*(l3*l3+ddh-c2)
 end subroutine

 pure real(wp) function logistic_negative(x) result(f)
  real(wp),intent(in)::x
  real(wp)::z
  if(x>=0)then;z=exp(-x);f=z/(1+z);else;f=1/(1+exp(x));endif
 end function

 pure subroutine bsk_pressure(rho,p,chi)
  real(wp),intent(in)::rho
  real(wp),intent(out)::p,chi
  real(wp)::x,den,a,da,f1,f2,f3,f4,b,c,d,g,h,t1,t2,lp
  x=log10(rho);den=1+PP(4)*x
  a=(PP(1)+PP(2)*x+PP(3)*x**3)/den;da=(PP(2)+3*PP(3)*x*x-PP(4)*a)/den
  f1=logistic_negative(PP(5)*(x-PP(6)));f2=logistic_negative(PP(9)*(PP(6)-x))
  f3=logistic_negative(PP(12)*(PP(13)-x));f4=logistic_negative(PP(16)*(PP(17)-x))
  b=PP(7)+PP(8)*x;c=PP(10)+PP(11)*x;d=PP(14)+PP(15)*x
  g=PP(20)*(x-PP(19));h=PP(23)*(x-PP(22))
  t1=PP(18)/(1+g*g);t2=PP(21)/(1+h*h)
  lp=a*f1+b*f2+c*f3+d*f4+t1+t2-33.2047_wp
  chi=da*f1-a*PP(5)*f1*(1-f1)+PP(8)*f2+b*PP(9)*f2*(1-f2)+PP(11)*f3+c*PP(12)*f3*(1-f3)+PP(15)*f4+d*PP(16)*f4*(1-f4)-2*t1*g*PP(20)/(1+g*g)-2*t2*h*PP(23)/(1+h*h)
  p=10._wp**lp
 end subroutine

 subroutine bsk22_state(n,mode,state,status,extended)
  real(wp),intent(in)::n
  character(len=*),intent(in)::mode
  type(lower_crust_state_t),intent(out)::state
  character(len=*),intent(out)::status
  logical,intent(in),optional::extended
  real(wp)::upper,e,de,dde,pfit,chi,penergy,x,acell,nf
  state=lower_crust_state_t();status='unsupported BSk22 density';upper=BSK_LOWER_MAX
  if(present(extended))then;if(extended)upper=BSK_EXT_MAX;endif
  if(.not.ieee_is_finite(n).or.n<BSK_DRIP.or.n>upper)return
  if(mode/='source_fit'.and.mode/='energy_derivative')then;status='invalid BSk22 thermodynamic mode';return;endif
  call bsk_energy(n,e,de,dde)
  state%nb=n;state%epsilon=n*(BSK_MN+e)
  call bsk_pressure(state%epsilon*RHO_FACTOR,pfit,chi);penergy=n*de
  state%source_epsilon=state%epsilon;state%source_pressure=pfit;state%energy_derivative_pressure=penergy
  state%source_thermo_residual=(pfit-penergy)/pfit
  if(mode=='source_fit')then
   state%pressure=pfit;state%sound2=chi*pfit/state%epsilon
   state%gamma=chi*(1+penergy/state%epsilon)
  else
   state%pressure=penergy;state%sound2=n*(de+dde)/(state%epsilon+penergy)
   state%gamma=1+dde/de
  endif
  state%ye=(6.08e-4_wp*n**(-.75_wp)-3.42e-3_wp)*(1+133.6_wp*n)*(1+(15.89_wp*n)**4)
  x=(n-BSK_DRIP)/(BSK_CORE-BSK_DRIP)
  state%yn=(173.2_wp*x+61.4_wp*x**1.86_wp)/(1+204.8_wp*x+57.3_wp*x*x)+.080_wp*x*exp(18*(x-1))
  state%yp=0;state%yh=1-state%yn;state%zav=40
  acell=40/state%ye;state%aav=acell*state%yh
  if(state%ye<=0.or.state%ye>=1.or.state%yn<0.or.state%yn>=1.or.state%aav<=40)then
   status='BSk22 analytic composition outside physical support';return
  endif
  state%ne=n*state%ye
  state%mu_e=.51099893_wp*sqrt(1+(386.159268_wp*(3*PI**2*n*state%ye)**(1._wp/3))**2)
  nf=n*state%yn
  state%mu_composition=BSK_MN+.211_wp*BSK_MN*nf**(2._wp/3)/(1+9.41_wp*sqrt(nf)+7.62_wp*nf-259*nf*nf)
  state%mu_b=(state%epsilon+state%pressure)/n
  state%composition_mu_residual=state%mu_composition-state%mu_b
  state%phase=2;state%representation='published BSk22 analytic ETFSI fit'
  status='ok'
 end subroutine

 subroutine bsk22_make_anchor(outer_density,outer_pressure,outer_mu_b,anchor,status,outer_epsilon,epsabs)
  real(wp),intent(in)::outer_density,outer_pressure,outer_mu_b
  type(bsk_anchor_t),intent(out)::anchor
  character(len=*),intent(out)::status
  real(wp),intent(in),optional::outer_epsilon,epsabs
  type(lower_crust_state_t)::s
  real(wp)::lo,hi,mid,el
  integer::i
  anchor=bsk_anchor_t();status='invalid BSk22 outer anchor'
  if(.not.all(ieee_is_finite([outer_density,outer_pressure,outer_mu_b])))return
  if(min(outer_density,outer_pressure,outer_mu_b)<=0.or.outer_density>=BSK_LOWER_MAX)return
  el=outer_density*outer_mu_b-outer_pressure
  if(present(outer_epsilon))el=outer_epsilon
  if(.not.ieee_is_finite(el).or.el<=0)return
  if(abs(el+outer_pressure-outer_density*outer_mu_b)>1e-9_wp)then
   status='outer anchor violates cold Euler identity';return
  endif
  if(present(epsabs))then
   if(.not.ieee_is_finite(epsabs).or.epsabs<=0)return
   anchor%epsabs=epsabs
  endif
  lo=max(BSK_DRIP,outer_density);hi=BSK_LOWER_MAX
  call bsk22_state(lo,'source_fit',s,status);if(status/='ok')return
  if(s%pressure>=outer_pressure)then;status='no positive-density-jump pressure anchor';return;endif
  call bsk22_state(hi,'source_fit',s,status);if(status/='ok')return
  if(s%pressure<=outer_pressure)then;status='no BSk pressure match in declared support';return;endif
  do i=1,64
   mid=(lo+hi)/2
   call bsk22_state(mid,'source_fit',s,status)
   if(s%pressure<outer_pressure)then;lo=mid;else;hi=mid;endif
   if(hi-lo<1e-17_wp)exit
  enddo
  anchor%n0=(lo+hi)/2
  if(anchor%n0<=outer_density*(1+1e-12_wp))then;status='nonpositive anchor density jump';return;endif
  call bsk22_state(anchor%n0,'source_fit',s,status)
  if(s%gamma<=0)then;status='nonpositive BSk anchor compressibility';return;endif
  anchor%outer_nb=outer_density;anchor%outer_p=outer_pressure;anchor%outer_mu=outer_mu_b;anchor%outer_epsilon=el
  anchor%e0=outer_mu_b-outer_pressure/anchor%n0;anchor%valid=.true.;status='ok'
 end subroutine

 subroutine bsk22_anchored_state(anchor,n,state,status,extended)
  type(bsk_anchor_t),intent(in)::anchor
  real(wp),intent(in)::n
  type(lower_crust_state_t),intent(out)::state
  character(len=*),intent(out)::status
  logical,intent(in),optional::extended
  real(wp)::upper,l,fa,fb,fm,whole,integral,error,e
  logical::ext
  state=lower_crust_state_t();status='density outside anchored BSk22 support';upper=BSK_LOWER_MAX;ext=.false.
  if(present(extended))ext=extended
  if(ext)upper=BSK_EXT_MAX
  if(.not.anchor%valid.or..not.ieee_is_finite(n).or.n<anchor%n0.or.n>upper)return
  call bsk22_state(n,'source_fit',state,status,ext);if(status/='ok')return
  l=log(n/anchor%n0)
  integral=0;error=0
  if(l>0)then
   fa=integrand(0._wp);fm=integrand(l/2);fb=integrand(l)
   whole=l*(fa+4*fm+fb)/6
   call adapt(0._wp,l,fa,fm,fb,whole,anchor%epsabs,0,integral,error)
   if(status/='ok')return
  endif
  e=anchor%e0+integral
  state%epsilon=n*e;state%mu_b=e+state%pressure/n
  state%energy_derivative_pressure=state%pressure;state%quadrature_error=error
  state%composition_mu_residual=state%mu_composition-state%mu_b
  state%sound2=state%gamma*state%pressure/(state%epsilon+state%pressure)
  state%representation='outer-anchored pressure integral; analytic composition'
 contains
  real(wp) function integrand(x) result(v)
   real(wp),intent(in)::x
   real(wp)::density,en,den,dden,pressure,chi
   density=anchor%n0*exp(x)
   call bsk_energy(density,en,den,dden)
   call bsk_pressure(density*(BSK_MN+en)*RHO_FACTOR,pressure,chi)
   v=pressure/density
  end function
  recursive subroutine adapt(a,b,fa,fm,fb,old,tol,depth,value,error)
   real(wp),intent(in)::a,b,fa,fm,fb,old,tol
   integer,intent(in)::depth
   real(wp),intent(out)::value,error
   real(wp)::m,lm,rm,fl,fr,left,right,delta,v1,v2,e1,e2
   m=(a+b)/2;lm=(a+m)/2;rm=(m+b)/2;fl=integrand(lm);fr=integrand(rm)
   left=(m-a)*(fa+4*fl+fm)/6;right=(b-m)*(fm+4*fr+fb)/6;delta=left+right-old
   if(abs(delta)<=15*tol)then
    value=left+right+delta/15;error=abs(delta)/15;return
   endif
   if(depth>=24)then;status='BSk pressure integration failed tolerance';value=left+right;error=abs(delta)/15;return;endif
   call adapt(a,m,fa,fl,fm,left,tol/2,depth+1,v1,e1)
   call adapt(m,b,fm,fr,fb,right,tol/2,depth+1,v2,e2)
   value=v1+v2;error=e1+e2
  end subroutine
 end subroutine

 subroutine lower_states_to_eos(states,mn_header,mp_header,eos,status)
  type(lower_crust_state_t),intent(in)::states(:)
  real(wp),intent(in)::mn_header,mp_header
  type(eos_data_t),intent(out)::eos
  character(len=*),intent(out)::status
  integer::i
  status='invalid native lower-crust rows'
  if(size(states)<1.or.mn_header<=0.or.mp_header<=0)return
  call allocate_eos(eos,size(states))
  eos%mn=mn_header;eos%mp=mp_header;eos%composition_file_present=.true.
  eos%composition_available=.true.;eos%beta_equilibrated=.true.
  do i=1,size(states)
   if(.not.all(ieee_is_finite([states(i)%nb,states(i)%pressure,states(i)%epsilon,states(i)%mu_b])))return
   if(states(i)%nb<=0.or.states(i)%pressure<=0)return
   if(i>1)then;if(states(i)%nb<=states(i-1)%nb)return;endif
   eos%nb(i)=states(i)%nb;eos%q(1,i)=states(i)%pressure/states(i)%nb
   eos%q(3,i)=states(i)%mu_b/mn_header-1;eos%q(4,i)=-states(i)%mu_e/mn_header
   eos%q(6,i)=states(i)%epsilon/(states(i)%nb*mn_header)-1;eos%q(7,i)=eos%q(6,i)
   eos%iphase(i)=states(i)%phase
   eos%yn(i)=states(i)%yn;eos%yp(i)=states(i)%yp;eos%ye(i)=states(i)%ye;eos%ymu(i)=states(i)%ymu
   eos%yh(i)=states(i)%yh;eos%abar(i)=states(i)%aav;eos%zbar(i)=states(i)%zav
  enddo
  eos%row_status=2;eos%yn_present=.true.;eos%yp_present=.true.;eos%ye_present=.true.;eos%ymu_present=.true.
  eos%abar_present=.true.;eos%zbar_present=.true.;eos%yh_present=.true.;eos%c3_applicable=.true.;eos%c4_applicable=.true.
  status='ok'
 end subroutine

 subroutine bsk22_generate(densities,mode,eos,states,status,anchor,extended)
  real(wp),intent(in)::densities(:)
  character(len=*),intent(in)::mode
  type(eos_data_t),intent(out)::eos
  type(lower_crust_state_t),allocatable,intent(out)::states(:)
  character(len=*),intent(out)::status
  type(bsk_anchor_t),intent(in),optional::anchor
  logical,intent(in),optional::extended
  logical::ext
  integer::i
  status='BSk generation requires at least three increasing finite densities'
  if(size(densities)<3.or..not.all(ieee_is_finite(densities)))return
  if(any(densities(2:)<=densities(:size(densities)-1)))return
  ext=.false.;if(present(extended))ext=extended
  allocate(states(size(densities)))
  do i=1,size(densities)
   if(mode=='outer_anchored')then
    if(.not.present(anchor))then;status='outer-anchored BSk requires explicit physical anchor';return;endif
    call bsk22_anchored_state(anchor,densities(i),states(i),status,ext)
   else
    call bsk22_state(densities(i),mode,states(i),status,ext)
   endif
   if(status/='ok')return
  enddo
  call lower_states_to_eos(states,BSK_MN,BSK_MN-1.2933_wp,eos,status)
 end subroutine
 ! Regenerate two distinct representations without inventing a thermodynamic seam.
 ! The BSk22 analytic source drip density lies slightly below regenerated outer drip;
 ! return separate monotone segments and explicit seam residuals, never reorder/offset.
 subroutine prepare_lower_source_native(cat,outer_points,lower_points,outer_eos,lower_eos,outer_states,lower_states,seam,status)
  type(mass_catalogue_t),intent(in)::cat
  integer,intent(in)::outer_points,lower_points
  type(eos_data_t),intent(out)::outer_eos,lower_eos
  type(lower_crust_state_t),allocatable,intent(out)::outer_states(:),lower_states(:)
  type(lower_source_seam_t),intent(out)::seam
  character(len=*),intent(out)::status
  real(wp),allocatable::densities(:)
  integer::i,k
  seam=lower_source_seam_t();status='source lower segment requires at least three points'
  if(lower_points<3)return
  call outer_sequence(cat,outer_points,outer_eos,outer_states,status)
  if(status/='ok')return
  allocate(densities(lower_points))
  do i=1,lower_points
   densities(i)=BSK_DRIP*(BSK_LOWER_MAX/BSK_DRIP)**(real(i-1,wp)/real(lower_points-1,wp))
  enddo
  densities(1)=BSK_DRIP;densities(lower_points)=BSK_LOWER_MAX
  call bsk22_generate(densities,'source_fit',lower_eos,lower_states,status)
  if(status/='ok')return
  k=size(outer_states)
  seam%outer_nb=outer_states(k)%nb;seam%lower_nb=lower_states(1)%nb
  seam%delta_nb=seam%lower_nb-seam%outer_nb
  seam%delta_pressure=lower_states(1)%pressure-outer_states(k)%pressure
  seam%delta_mu_b=lower_states(1)%mu_b-outer_states(k)%mu_b
  seam%delta_epsilon=lower_states(1)%epsilon-outer_states(k)%epsilon
  seam%density_overlap=seam%delta_nb<0
  seam%pressure_continuous=abs(seam%delta_pressure)<1e-12_wp
  seam%mu_continuous=abs(seam%delta_mu_b)<1e-9_wp
  status='ok'
 end subroutine
end module compose_eos_lower_crust_native
