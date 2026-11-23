! Myers & Swiatecki, Nuclear Physics 81 (1966) 1-60,
! DOI 10.1016/S0029-5582(66)80001-9, Eqs. (5),(6),(9)-(11), pp.13,21-26,34.
! This is a NEW formula-generated mass model, not the UCRL-11980 Appendix D table.
! Nuclear local minima of the cubic deformation expansion are used; the polynomial
! is unbounded at large theta, so an unconstrained global minimum is meaningless.
module compose_eos_myers1966
 use compose_eos_types, only: wp,eos_data_t,allocate_eos
 use compose_eos_bps_primary, only: bps_state_t,bps_isotope_state,BPS_MN,BPS_MP,BPS_ME,BPS_HBARC
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 real(wp),parameter :: A1=15.677_wp,A2=18.56_wp,C3=.717_wp,C4=1.21129_wp,KAPPA=1.79_wp
 real(wp),parameter :: SHELL_C=5.8_wp,SHELL_c0=.26_wp,AR0=.27_wp,DELTA_N=8.07144_wp
 real(wp),parameter :: DELTA_H=7.28899_wp,PAIR=11._wp,PI=acos(-1._wp)
 real(wp),parameter :: E2=BPS_HBARC/137.035999084_wp
 real(wp),parameter :: MU_CERT_TOL=5e-11_wp
 integer,parameter :: MAX_PHASES=256
 type,public :: myers_mass_t
  integer :: n=0,z=0
  logical :: valid=.false.
  real(wp) :: excess=0,binding=0,mass=0,theta=0,alpha=0,shell=0,e=0,f=0
  real(wp) :: curvature=0,stationarity=0
 end type
 type,public :: myers_catalog_t
  integer :: count=0,magic=0,nmax=0,zmax=0,shapes=0,particle_rejected=0
  integer,allocatable :: n(:),z(:),a(:)
  real(wp),allocatable :: excess(:),binding(:),mass(:),theta(:),alpha(:),sn(:),sp(:)
 end type
 type :: pressure_point_t
  real(wp) :: p=0
  integer :: winner=0
  real(wp),allocatable :: ne(:),mu(:),nb(:)
 end type
 public :: myers_shell_f,myers_deformation_minimum,myers_mass,load_myers_source
 public :: build_myers_catalog,myers_state_at_pressure,myers_isotope_at_pressure
 public :: generate_myers_bps_native
contains
 pure real(wp) function myers_shell_f(x,magic,neutron) result(v)
  integer,intent(in)::x,magic
  logical,intent(in)::neutron
  integer :: m(10),i,top
  real(wp) :: l,r,q
  m=[0,2,8,magic,28,50,82,126,184,258]
  top=9;if(neutron)top=10
  v=huge(1._wp)
  if(x<0.or.x>m(top).or.(magic/=14.and.magic/=20))return
  do i=2,top
   if(x>m(i))cycle
   l=real(m(i-1),wp);r=real(m(i),wp)
   q=.6_wp*(r**(5._wp/3)-l**(5._wp/3))/(r-l)
   v=q*(x-l)-.6_wp*(real(x,wp)**(5._wp/3)-l**(5._wp/3))
   if(x==m(i).or.x==m(i-1))v=0
   return
  enddo
 end function

 pure subroutine myers_deformation_minimum(e,f,s,theta,v,curvature,residual,valid)
  real(wp),intent(in)::e,f,s
  real(wp),intent(out)::theta,v,curvature,residual
  logical,intent(out)::valid
  real(wp)::target,lo,hi,mid,ta,tb,fa,fb,t,vt,ct
  integer::j
  theta=0;v=s;curvature=2*(e-s);residual=0
  valid=curvature>0
  if(f<0.or..not.all(ieee_is_finite([e,f,s])))then;valid=.false.;return;endif
  if(f==0)then
   if(e==s.and.e>0)valid=.true.
   if(e>0.and.s>e)then
    theta=sqrt(log(s/e));v=e*(1+log(s/e));curvature=4*e*theta**2;valid=.true.
   endif
   return
  endif
  if(s<=0)return
  target=.75_wp*f/s
  if(target>=1/sqrt(2*exp(1._wp)))return
  lo=0;hi=1/sqrt(2._wp)
  do j=1,70
   mid=(lo+hi)/2
   if(mid*exp(-mid*mid)<target)then;lo=mid;else;hi=mid;endif
  enddo
  ta=(lo+hi)/2;lo=1/sqrt(2._wp);hi=2
  do while(hi*exp(-hi*hi)>target)
   hi=hi*1.5_wp
   if(hi>26)exit
  enddo
  do j=1,70
   mid=(lo+hi)/2
   if(mid*exp(-mid*mid)>target)then;lo=mid;else;hi=mid;endif
  enddo
  tb=(lo+hi)/2
  fa=e-1.5_wp*f*ta-s*exp(-ta*ta)
  fb=e-1.5_wp*f*tb-s*exp(-tb*tb)
  if(fa>=0.or.fb<=0)return
  lo=ta;hi=tb
  do j=1,75
   mid=(lo+hi)/2
   if(e-1.5_wp*f*mid-s*exp(-mid*mid)<0)then;lo=mid;else;hi=mid;endif
  enddo
  t=(lo+hi)/2;vt=e*t*t-f*t**3+s*exp(-t*t)
  ct=2*e-6*f*t+s*(4*t*t-2)*exp(-t*t)
  if(ct<=0)return
  if(.not.valid.or.vt<v)then
   theta=t;v=vt;curvature=ct
   residual=abs(2*t*(e-1.5_wp*f*t-s*exp(-t*t)))
   valid=.true.
  endif
 end subroutine

 pure subroutine myers_mass(n,z,magic,m)
  integer,intent(in)::n,z,magic
  type(myers_mass_t),intent(out)::m
  real(wp)::a,a13,asym,c1n,c2n,alpha0,delta,v,cohesion
  m=myers_mass_t();m%n=n;m%z=z
  if(n<1.or.z<1.or.n>258.or.z>184.or.(magic/=14.and.magic/=20))return
  a=real(n+z,wp);a13=a**(1._wp/3)
  asym=real(n-z,wp)/a;cohesion=1-KAPPA*asym**2
  if(cohesion<=0)return
  c1n=A1*cohesion;c2n=A2*cohesion;alpha0=sqrt(5._wp)*AR0/a13
  m%shell=SHELL_C*((myers_shell_f(n,magic,.true.)+myers_shell_f(z,magic,.false.))/ &
    (a/2)**(2._wp/3)-SHELL_c0*a13)
  m%e=(.4_wp*c2n*a13**2-.2_wp*C3*real(z,wp)**2/a13)*alpha0**2
  m%f=4._wp/105*(c2n*a13**2+C3*real(z,wp)**2/a13)*alpha0**3
  call myers_deformation_minimum(m%e,m%f,m%shell,m%theta,v,m%curvature,m%stationarity,m%valid)
  if(.not.m%valid)return
  delta=0
  if(mod(n,2)==0.and.mod(z,2)==0)delta=-PAIR/sqrt(a)
  if(mod(n,2)==1.and.mod(z,2)==1)delta=PAIR/sqrt(a)
  m%alpha=alpha0*m%theta
  m%excess=DELTA_N*n+DELTA_H*z-c1n*a+c2n*a13**2+C3*real(z,wp)**2/a13-C4*real(z,wp)**2/a+delta+v
  ! Retain published mass-excess differences, anchor the neutron to BPS_MN,
  ! and subtract electron rest mass. No extra atomic binding correction is added.
  m%mass=a*(BPS_MN-DELTA_N)+m%excess-z*BPS_ME
  m%binding=(n*BPS_MN+z*BPS_MP-m%mass)/a
  m%valid=all(ieee_is_finite([m%excess,m%binding,m%mass,m%theta,m%curvature]))
 end subroutine

 subroutine load_myers_source(path,status)
  character(*),intent(in)::path
  character(*),intent(out)::status
  character(32),parameter :: keys(11)=[character(32)::'a1_MeV','a2_MeV','c3_MeV','c4_MeV','kappa', &
     'shell_C_MeV','shell_c','a_over_r0','neutron_excess_MeV','hydrogen_excess_MeV','pairing_MeV']
  character(20),parameter :: values(11)=[character(20)::'15.677','18.56','0.717','1.21129','1.79', &
     '5.8','0.26','0.27','8.07144','7.28899','11']
  character(512)::line
  integer::u,ios,i,k,c
  status='Myers cannot open source constants'
  open(newunit=u,file=path,status='old',action='read',iostat=ios);if(ios/=0)return
  i=0
  do
   read(u,'(a)',iostat=ios)line
   if(ios<0)exit
   if(ios/=0.or.len_trim(line)==len(line))then
    status='Myers unreadable or overlong source line';close(u);return
   endif
   line=adjustl(line)
   if(len_trim(line)==0)cycle
   if(line(1:1)=='#')cycle
   if(i==0)then
    if(trim(line)/='parameter,value')then;status='Myers source header mismatch';close(u);return;endif
    i=1;cycle
   endif
   if(i>11)then;status='Myers source requires exactly 11 constants';close(u);return;endif
   k=index(line,',');c=scan(line(k+1:),',')
   if(k<2.or.c>0)then;status='Myers source requires two fields';close(u);return;endif
   ! Canonical decimal text is intentional: this source-locked model is not a fit
   ! interface, and hostile NaN/overflow tokens never enter a real conversion.
   if(trim(line(:k-1))/=trim(keys(i)).or.trim(adjustl(line(k+1:)))/=trim(values(i)))then
    status='Myers source constant or canonical decimal differs from published model';close(u);return
   endif
   i=i+1
  enddo
  close(u)
  if(i/=12)then;status='Myers source requires exactly 11 constants';return;endif
  status='ok'
 end subroutine

 subroutine build_myers_catalog(magic,cat,status,nmax,zmax)
  integer,intent(in)::magic
  type(myers_catalog_t),intent(out)::cat
  character(*),intent(out)::status
  integer,intent(in),optional::nmax,zmax
  type(myers_mass_t),allocatable::grid(:,:)
  real(wp)::sn,sp
  integer::nn,zz,n,z,k,pass
  nn=258;zz=184
  if(present(nmax))nn=nmax
  if(present(zmax))zz=zmax
  status='Myers invalid finite catalogue domain'
  if(nn<4.or.nn>258.or.zz<4.or.zz>184.or.(magic/=14.and.magic/=20))return
  cat%magic=magic;cat%nmax=nn;cat%zmax=zz
  allocate(grid(nn,zz))
  !$omp parallel do collapse(2) default(none) shared(grid,nn,zz,magic) private(n,z) schedule(static)
  do z=1,zz
   do n=1,nn
    call myers_mass(n,z,magic,grid(n,z))
   enddo
  enddo
  !$omp end parallel do
  k=0
  do pass=1,2
   k=0
   do z=2,zz
    do n=2,nn
     if(.not.grid(n,z)%valid)cycle
     if(pass==1)cat%shapes=cat%shapes+1
     if(.not.grid(n-1,z)%valid.or..not.grid(n,z-1)%valid)cycle
     sn=grid(n-1,z)%excess+DELTA_N-grid(n,z)%excess
     sp=grid(n,z-1)%excess+DELTA_H-grid(n,z)%excess
     ! Published Fig.19: neutron bound; proton no more than 3 MeV unbound.
     if(sn<0.or.sp< -3._wp.or.grid(n,z)%binding<=0.or.grid(n,z)%binding>=20)cycle
     k=k+1
     if(pass==1)cycle
     cat%n(k)=n;cat%z(k)=z;cat%a(k)=n+z
     cat%excess(k)=grid(n,z)%excess;cat%binding(k)=grid(n,z)%binding;cat%mass(k)=grid(n,z)%mass
     cat%theta(k)=grid(n,z)%theta;cat%alpha(k)=grid(n,z)%alpha;cat%sn(k)=sn;cat%sp(k)=sp
    enddo
   enddo
   if(pass==1)then
    cat%count=k
    if(k<2)then;status='Myers finite catalogue has fewer than two candidates';return;endif
    allocate(cat%n(k),cat%z(k),cat%a(k),cat%excess(k),cat%binding(k),cat%mass(k),cat%theta(k), &
      cat%alpha(k),cat%sn(k),cat%sp(k))
   endif
  enddo
  cat%particle_rejected=cat%shapes-cat%count
  status='ok'
 end subroutine

 subroutine electron_at_pressure(p,z,ne,hz,status)
  real(wp),intent(in)::p
  integer,intent(in)::z
  real(wp),intent(out)::ne,hz
  character(*),intent(out)::status
  real(wp)::lo,hi,x,next,deriv,wl
  type(bps_state_t)::s
  integer::j
  status='Myers invalid electron pressure'
  ne=0;hz=0
  if(.not.ieee_is_finite(p).or.p<=0.or.p>1e-2_wp.or.z<1.or.z>184)return
  lo=log(1e-24_wp);hi=log(1e-2_wp)
  x=log(max(1e-20_wp,min(1e-3_wp,(12*PI**2*p/BPS_HBARC)**.75_wp/(3*PI**2))))
  do j=1,96
   call bps_isotope_state(exp(x),z,z,1._wp,s,status);if(status/='ok')return
   if(s%pressure<p)then;lo=x;else;hi=x;endif
   if(abs(s%pressure/p-1)<2e-13_wp.or.hi-lo<5e-15_wp)exit
   deriv=s%gamma*s%pressure
   next=(lo+hi)/2
   if(deriv>0)then
    next=x-(s%pressure-p)/deriv
    if(next<=lo.or.next>=hi)next=(lo+hi)/2
   endif
   x=next
  enddo
  ne=exp(x)
  if(s%gamma<=0)then;status='Myers pressure branch has nonpositive compressibility';return;endif
  if(abs(s%pressure/p-1)>1e-10_wp)then;status='Myers electron pressure inversion failed';return;endif
  wl=-1.819620_wp*real(z,wp)**2*E2*(ne/(2*z))**(1._wp/3)
  hz=s%mu_e+4*wl/(3*z)
  status='ok'
 end subroutine

 subroutine make_pressure_point(cat,p,point,status)
  type(myers_catalog_t),intent(in)::cat
  real(wp),intent(in)::p
  type(pressure_point_t),intent(out)::point
  character(*),intent(out)::status
  real(wp)::hz,ne,cacheh(cat%zmax)
  logical::needed(cat%zmax)
  integer::i,z
  point%p=p
  allocate(point%ne(cat%zmax),point%mu(cat%count),point%nb(cat%count))
  needed=.false.;cacheh=0;point%ne=0
  do i=1,cat%count
   needed(cat%z(i))=.true.
  enddo
  do z=1,cat%zmax
   if(.not.needed(z))cycle
   call electron_at_pressure(p,z,ne,hz,status);if(status/='ok')return
   point%ne(z)=ne;cacheh(z)=hz
  enddo
  do i=1,cat%count
   z=cat%z(i)
   point%mu(i)=(cat%mass(i)+z*cacheh(z))/cat%a(i)
   point%nb(i)=point%ne(z)*cat%a(i)/z
  enddo
  point%winner=minloc(point%mu,dim=1);status='ok'
 end subroutine

 subroutine myers_isotope_at_pressure(cat,i,p,s,status)
  type(myers_catalog_t),intent(in)::cat
  integer,intent(in)::i
  real(wp),intent(in)::p
  type(bps_state_t),intent(out)::s
  character(*),intent(out)::status
  real(wp)::ne,hz
  s=bps_state_t();status='Myers invalid isotope index'
  if(i<1.or.i>cat%count)return
  call electron_at_pressure(p,cat%z(i),ne,hz,status);if(status/='ok')return
  call bps_isotope_state(ne*cat%a(i)/cat%z(i),cat%z(i),cat%a(i),cat%binding(i),s,status)
  if(status/='ok')return
  s%isotope=i;s%pressure_residual=abs(s%pressure/p-1);s%pressure=p
 end subroutine

 subroutine myers_state_at_pressure(cat,p,s,status)
  type(myers_catalog_t),intent(in)::cat
  real(wp),intent(in)::p
  type(bps_state_t),intent(out)::s
  character(*),intent(out)::status
  type(pressure_point_t)::point
  call make_pressure_point(cat,p,point,status);if(status/='ok')return
  call myers_isotope_at_pressure(cat,point%winner,p,s,status)
 end subroutine

 recursive subroutine walk_phase_interval(cat,l,r,depth,ids,nph,bounds,left,right,ncert,nquery,status)
  type(myers_catalog_t),intent(in)::cat
  type(pressure_point_t),intent(in)::l,r
  integer,intent(in)::depth
  integer,intent(inout)::ids(:),nph,ncert,nquery
  real(wp),intent(inout)::bounds(0:)
  type(bps_state_t),intent(inout)::left(:),right(:)
  character(*),intent(out)::status
  type(pressure_point_t)::mid,atroot
  type(bps_state_t)::sl,sr
  real(wp)::dp,den,t,up,secant,p,lo,hi,delta
  integer::a,b,k,j
  logical::certified,ordered
  a=l%winner;b=r%winner;dp=r%p-l%p
  ! Each fixed-nucleus Gibbs curve is concave: dmu/dP=1/nb,
  ! d2mu/dP2=-(dnb/dP)/nb**2<0. Secants are lower bounds, tangents
  ! upper bounds. Their piecewise-linear comparison excludes hidden phases.
  den=1/l%nb(a)-1/r%nb(b)
  if(den/=0)then
   t=(r%mu(b)-l%mu(a)-dp/r%nb(b))/den
  else
   t=0
  endif
  t=max(0._wp,min(dp,t))
  up=min(l%mu(a)+t/l%nb(a),r%mu(b)+(t-dp)/r%nb(b))
  certified=.true.
  do k=1,cat%count
   if(k==a.or.k==b)cycle
   secant=l%mu(k)+(r%mu(k)-l%mu(k))*(t/dp)
   if(secant<up-MU_CERT_TOL)then;certified=.false.;exit;endif
  enddo
  ordered=a==b
  if(a/=b)ordered=l%nb(b)>r%nb(a)
  if(certified.and.ordered)then
   ncert=ncert+1
   if(a==b)then;status='ok';return;endif
   if(nph>=size(ids).or.ids(nph)/=a)then;status='Myers phase chain overflow or inconsistency';return;endif
   lo=log(l%p);hi=log(r%p)
   do j=1,70
    p=exp((lo+hi)/2)
    call myers_isotope_at_pressure(cat,a,p,sl,status);if(status/='ok')return
    call myers_isotope_at_pressure(cat,b,p,sr,status);if(status/='ok')return
    delta=sl%mu_b-sr%mu_b
    if(delta<0)then;lo=log(p);else;hi=log(p);endif
    if(hi-lo<2e-13_wp)exit
   enddo
   p=exp((lo+hi)/2)
   call myers_isotope_at_pressure(cat,a,p,sl,status);if(status/='ok')return
   call myers_isotope_at_pressure(cat,b,p,sr,status);if(status/='ok')return
   call make_pressure_point(cat,p,atroot,status);if(status/='ok')return
   nquery=nquery+1
   if(abs(sl%mu_b-sr%mu_b)>1e-9_wp.or.minval(atroot%mu)<min(sl%mu_b,sr%mu_b)-1e-9_wp)then
    status='Myers phase crossing failed global Gibbs check';return
   endif
   if(sr%nb<=sl%nb.or.sr%epsilon<=sl%epsilon)then
    status='Myers nonpositive coexistence density or energy jump';return
   endif
   bounds(nph)=p;left(nph)=sl;right(nph)=sr
   nph=nph+1;ids(nph)=b
   status='ok';return
  endif
  if(depth>=48.or.dp<2e-13_wp*sqrt(l%p*r%p))then
   status='Myers unresolved phase interval at certificate resolution';return
  endif
  call make_pressure_point(cat,sqrt(l%p*r%p),mid,status);if(status/='ok')return
  nquery=nquery+1
  call walk_phase_interval(cat,l,mid,depth+1,ids,nph,bounds,left,right,ncert,nquery,status)
  if(status/='ok')return
  call walk_phase_interval(cat,mid,r,depth+1,ids,nph,bounds,left,right,ncert,nquery,status)
 end subroutine

 subroutine generate_myers_bps_native(source_path,n_points,eos,output_dir,status,magic_choice)
  character(*),intent(in)::source_path,output_dir
  integer,intent(in)::n_points,magic_choice
  type(eos_data_t),intent(out)::eos
  character(*),intent(out)::status
  type(myers_catalog_t)::cat
  type(pressure_point_t)::pstart,pend
  type(bps_state_t),allocatable::rows(:)
  type(bps_state_t)::s,tmp,drip,left(MAX_PHASES),right(MAX_PHASES),s1,s2,sm1,sm2
  real(wp)::p0,pdrip,lo,hi,p,bounds(0:MAX_PHASES),maxpr,maxeuler,maxjump,fdres,h
  integer::ids(MAX_PHASES),i,j,k,nph,nr,nkeep,ncert,nquery,u,ios,closeios,active
  logical::seen(MAX_PHASES),boundary
  character(2)::magic_text
  status='Myers requires 16 <= n_points <= 100000'
  if(n_points<16.or.n_points>100000)return
  if(magic_choice/=14.and.magic_choice/=20)then;status='Myers magic choice must be 14 or 20';return;endif
  call load_myers_source(source_path,status);if(status/='ok')return
  call build_myers_catalog(magic_choice,cat,status);if(status/='ok')return
  ! Same low-pressure control as historical BPS generator, not a density claim.
  call bps_isotope_state(6.295e-12_wp,26,56,8.79_wp,s,status);if(status/='ok')return
  p0=s%pressure;lo=log(1e-5_wp);hi=log(1e-2_wp);nquery=0
  call myers_state_at_pressure(cat,exp(lo),s,status);if(status/='ok')return
  if(s%mu_b>=BPS_MN)then;status='Myers drip lacks low-pressure bracket';return;endif
  call myers_state_at_pressure(cat,1e-2_wp,s,status);if(status/='ok')return
  if(s%mu_b<=BPS_MN)then;status='Myers drip lacks high-pressure bracket';return;endif
  do i=1,64
   p=exp((lo+hi)/2)
   call myers_state_at_pressure(cat,p,s,status);if(status/='ok')return
   if(s%mu_b<BPS_MN)then;lo=log(p);else;hi=log(p);endif
   if(hi-lo<2e-13_wp)exit
  enddo
  pdrip=exp((lo+hi)/2)
  call myers_state_at_pressure(cat,pdrip,drip,status);if(status/='ok')return
  if(abs(drip%mu_b-BPS_MN)>1e-9_wp)then;status='Myers neutron drip equality failed';return;endif
  call make_pressure_point(cat,p0,pstart,status);if(status/='ok')return
  call make_pressure_point(cat,pdrip,pend,status);if(status/='ok')return
  nph=1;ids=0;ids(1)=pstart%winner;ncert=0
  bounds=0;bounds(0)=p0
  call walk_phase_interval(cat,pstart,pend,0,ids,nph,bounds,left,right,ncert,nquery,status)
  if(status/='ok')return
  bounds(nph)=pdrip
  if(ids(nph)/=pend%winner)then;status='Myers terminal phase chain mismatch';return;endif
  do i=1,nph
   if(bounds(i)<=bounds(i-1))then;status='Myers phase pressure intervals are not increasing';return;endif
  enddo
  allocate(rows(n_points+3*nph));nr=0;active=1
  do i=1,n_points
   p=exp(log(p0)+log(pdrip/p0)*real(i-1,wp)/(n_points-1))
   if(i==1)p=p0
   if(i==n_points)p=pdrip
   do while(active<nph)
    if(p<=bounds(active))exit
    active=active+1
   enddo
   nr=nr+1
   call myers_isotope_at_pressure(cat,ids(active),p,rows(nr),status);if(status/='ok')return
  enddo
  do i=1,nph-1
   nr=nr+1;rows(nr)=left(i)
   nr=nr+1;rows(nr)=right(i)
  enddo
  fdres=0
  do i=1,nph
   p=sqrt(bounds(i-1)*bounds(i))
   call myers_state_at_pressure(cat,p,s,status);if(status/='ok')return
   if(s%isotope/=ids(i))then;status='Myers phase interior failed global minimum check';return;endif
   nr=nr+1;rows(nr)=s;k=ids(i);h=1e-3_wp*s%nb
   call bps_isotope_state(s%nb+h,cat%z(k),cat%a(k),cat%binding(k),s1,status);if(status/='ok')return
   call bps_isotope_state(s%nb+2*h,cat%z(k),cat%a(k),cat%binding(k),s2,status);if(status/='ok')return
   call bps_isotope_state(s%nb-h,cat%z(k),cat%a(k),cat%binding(k),sm1,status);if(status/='ok')return
   call bps_isotope_state(s%nb-2*h,cat%z(k),cat%a(k),cat%binding(k),sm2,status);if(status/='ok')return
   fdres=max(fdres,abs((s%nb*(-s2%kinetic_lattice+8*s1%kinetic_lattice- &
     8*sm1%kinetic_lattice+sm2%kinetic_lattice)/(12*h)-s%kinetic_lattice)/s%pressure-1))
  enddo
  if(fdres>3e-7_wp)then;status='Myers finite-difference first-law check failed';return;endif
  ! Only O(number_of_phases*N) insertions: the first requested_points rows
  ! are already ordered; no quadratic full-catalogue sort or per-row Gibbs scan.
  do i=2,nr
   tmp=rows(i);j=i-1
   do while(j>=1)
    if(rows(j)%nb<=tmp%nb)exit
    rows(j+1)=rows(j);j=j-1
   enddo
   rows(j+1)=tmp
  enddo
  nkeep=1
  do i=2,nr
   if(rows(i)%nb>rows(nkeep)%nb*(1+1e-13_wp))then
    nkeep=nkeep+1;rows(nkeep)=rows(i)
   elseif(rows(i)%isotope/=rows(nkeep)%isotope)then
    status='Myers coexistence phases have unresolved density separation';return
   endif
  enddo
  seen=.false.;maxpr=0;maxeuler=0;maxjump=0;boundary=.false.
  do i=1,nkeep
   active=0
   do j=1,nph
    if(ids(j)==rows(i)%isotope.and.rows(i)%pressure>=bounds(j-1)*(1-3e-13_wp).and. &
       rows(i)%pressure<=bounds(j)*(1+3e-13_wp))active=j
   enddo
   if(active==0)then;status='Myers generated row absent from phase plan';return;endif
   seen(active)=.true.;k=rows(i)%isotope
   boundary=boundary.or.cat%n(k)==2.or.cat%n(k)==cat%nmax.or.cat%z(k)==2.or.cat%z(k)==cat%zmax
   maxpr=max(maxpr,rows(i)%pressure_residual)
   maxeuler=max(maxeuler,abs((rows(i)%epsilon+rows(i)%pressure)/rows(i)%nb-rows(i)%mu_b))
   if(.not.all(ieee_is_finite([rows(i)%epsilon,rows(i)%pressure,rows(i)%sound2])))then
    status='Myers generated nonfinite state';return
   endif
   if(rows(i)%pressure<=0.or.rows(i)%sound2<=0.or.rows(i)%sound2>=1)then
    status='Myers generated unstable or acausal branch';return
   endif
  enddo
  if(.not.all(seen(1:nph)))then;status='Myers generator omitted a certified phase';return;endif
  if(maxeuler>1e-9_wp.or.maxpr>1e-10_wp)then;status='Myers thermodynamic identity failed';return;endif
  if(boundary)then;status='Myers selected phase touches finite catalogue boundary';return;endif
  do i=1,nph-1
   maxjump=max(maxjump,abs(left(i)%mu_b-right(i)%mu_b))
  enddo
  call allocate_eos(eos,nkeep);eos%mn=BPS_MN;eos%mp=BPS_MP
  eos%composition_available=.true.;eos%composition_file_present=.true.;eos%beta_equilibrated=.true.
  do i=1,nkeep
   s=rows(i);active=0
   do j=1,nph
    if(ids(j)==s%isotope.and.s%pressure>=bounds(j-1)*(1-3e-13_wp).and. &
       s%pressure<=bounds(j)*(1+3e-13_wp))active=j
   enddo
   eos%nb(i)=s%nb;eos%q(1,i)=s%pressure/s%nb
   eos%q(3,i)=s%mu_b/eos%mn-1;eos%q(4,i)=-s%mu_e/eos%mn
   eos%q(6,i)=s%epsilon/(s%nb*eos%mn)-1;eos%q(7,i)=eos%q(6,i)
   eos%iphase(i)=200+active
   eos%yn(i)=0;eos%yp(i)=0;eos%ye(i)=real(s%z,wp)/s%a;eos%ymu(i)=0
   eos%yh(i)=1;eos%abar(i)=s%a;eos%zbar(i)=s%z
  enddo
  eos%row_status=2;eos%yn_present=.true.;eos%yp_present=.true.;eos%ye_present=.true.;eos%ymu_present=.true.
  eos%yh_present=.true.;eos%abar_present=.true.;eos%zbar_present=.true.
  eos%c3_applicable=.true.;eos%c4_applicable=.true.
  open(newunit=u,file=trim(output_dir)//'/myers_mass_catalog.csv',status='replace',iostat=ios)
  if(ios/=0)then;status='Myers cannot write mass catalogue';return;endif
  write(u,'(a)',iostat=ios)'N,Z,A,atomic_mass_excess_MeV,nuclear_mass_MeV,binding_MeV_per_baryon,theta,alpha,Sn_MeV,Sp_MeV'
  do i=1,cat%count
   if(ios/=0)exit
   write(u,'(3(i0,","),6(es25.16e3,","),es25.16e3)',iostat=ios) &
    cat%n(i),cat%z(i),cat%a(i),cat%excess(i),cat%mass(i),cat%binding(i),cat%theta(i),cat%alpha(i),cat%sn(i),cat%sp(i)
  enddo
  close(u,iostat=closeios)
  if(ios/=0.or.closeios/=0)then;status='Myers mass catalogue write failed';return;endif
  open(newunit=u,file=trim(output_dir)//'/myers_bps_states.csv',status='replace',iostat=ios)
  if(ios/=0)then;status='Myers cannot write states diagnostic';return;endif
  write(u,'(a)',iostat=ios)'n_b,epsilon,pressure,mu_b,mu_e,Z,A,phase,sound2,Gamma,pressure_residual'
  do i=1,nkeep
   if(ios/=0)exit
   s=rows(i)
   write(u,'(5(es25.16e3,","),3(i0,","),2(es25.16e3,","),es25.16e3)',iostat=ios) &
    s%nb,s%epsilon,s%pressure,s%mu_b,s%mu_e,s%z,s%a,eos%iphase(i),s%sound2,s%gamma,s%pressure_residual
  enddo
  close(u,iostat=closeios)
  if(ios/=0.or.closeios/=0)then;status='Myers states diagnostic write failed';return;endif
  open(newunit=u,file=trim(output_dir)//'/myers_bps_transitions.csv',status='replace',iostat=ios)
  if(ios/=0)then;status='Myers cannot write transitions diagnostic';return;endif
  write(u,'(a)',iostat=ios)'phase_low,phase_high,pressure,n_low,n_high,epsilon_low,epsilon_high,mu_low,mu_high'
  do i=1,nph-1
   if(ios/=0)exit
   write(u,'(2(i0,","),6(es25.16e3,","),es25.16e3)',iostat=ios)i+200,i+201,bounds(i),left(i)%nb, &
    right(i)%nb,left(i)%epsilon,right(i)%epsilon,left(i)%mu_b,right(i)%mu_b
  enddo
  close(u,iostat=closeios)
  if(ios/=0.or.closeios/=0)then;status='Myers transitions diagnostic write failed';return;endif
  open(newunit=u,file=trim(output_dir)//'/MYERS_BPS_GENERATION.json',status='replace',iostat=ios)
  if(ios/=0)then;status='Myers cannot write validation diagnostic';return;endif
  write(magic_text,'(i2)')magic_choice
  write(u,'(a)',iostat=ios)'{"model":"BPS_Myers1966_magic'//magic_text//'","doi":"10.1016/S0029-5582(66)80001-9",'
  if(ios==0)write(u,'(a)',iostat=ios)' "BPS_doi":"10.1086/151216","full_UCRL11980_catalogue_reproduced":false,'
  if(ios==0)write(u,'(a)',iostat=ios)' "mass_convention":"published atomic excess; physical neutron anchor; electron rest subtraction",'
  if(ios==0)write(u,'(a)',iostat=ios)' "bound_atomic_electron_corrections":false,"Wigner_term":false,'
  if(ios==0)write(u,'(a)',iostat=ios)' "deformation":"stable local minima of cubic Eq10, prolate and spherical compared",'
  if(ios==0)write(u,'(a)',iostat=ios)' "selection":"Sn>=0, Sp>=-3 MeV, positive binding and stable shape",'
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)' "magic_choice":',magic_choice,','
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)' "requested_points":',n_points,','
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)' "rows":',nkeep,','
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)' "phase_count":',nph,','
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)' "coexistence_pairs":',nph-1,','
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)' "candidate_count":',cat%count,','
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)' "stable_shape_count":',cat%shapes,','
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)' "particle_or_binding_rejected":',cat%particle_rejected,','
  if(ios==0)write(u,'(a)',iostat=ios)' "N_min":2,"N_max":258,"Z_min":2,"Z_max":184,"selected_boundary_hit":false,'
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)' "certified_pressure_intervals":',ncert,','
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)' "phase_search_global_queries":',nquery+2,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "Gibbs_certificate_tolerance_MeV":',MU_CERT_TOL,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "n_min":',rows(1)%nb,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "n_max":',drip%nb,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "pressure_min":',p0,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "pressure_drip":',pdrip,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "neutron_drip_mu_e_MeV":',drip%mu_e,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "neutron_drip_delta_mu_MeV":',abs(drip%mu_b-BPS_MN),','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "first_law_fd_relative_max":',fdres,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "pressure_inversion_relative_max":',maxpr,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "Euler_identity_MeV_max":',maxeuler,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "transition_delta_mu_MeV_max":',maxjump,'}'
  close(u,iostat=closeios)
  if(ios/=0.or.closeios/=0)then;status='Myers validation diagnostic write failed';return;endif
  status='ok'
 end subroutine
end module compose_eos_myers1966
