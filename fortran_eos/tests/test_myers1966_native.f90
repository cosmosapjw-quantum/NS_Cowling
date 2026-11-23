! Independent validation: Eq.(9) in real128 with a scan/golden-section
! deformation search, versus production Eq.(10) root classification in real64.
program test_myers1966_native
 use, intrinsic :: iso_fortran_env, only: real128,int64
 use compose_eos_types, only: wp,eos_data_t
 use compose_eos_bps_primary, only: bps_state_t,bps_isotope_state,BPS_MN
 use compose_eos_myers1966
 implicit none
 integer,parameter::qp=real128,nn(12)=[2,8,20,28,30,36,50,82,102,126,146,184]
 integer,parameter::zz(12)=[2,8,20,24,26,28,36,36,72,82,92,126]
 type(myers_mass_t)::m
 type(myers_catalog_t)::cat,small
 type(bps_state_t)::best,brute,s,ss
 real(wp)::theta,v,curv,resid,excess128,theta128,maxmass,maxtheta,maxmu,maxn,bmin,p,lo,hi,mid
 real(wp)::cached_seconds=0,brute_seconds=0
 integer(int64)::tick0,tick1,rate
 real(qp)::qa13,qc2,qshell,qalpha0,qz
 integer::magic,i,j,k,ib,iv,checks,magics(2)=[14,20]
 logical::valid
 character(256)::status
 maxmass=0;maxtheta=0;maxmu=0;maxn=0;checks=0
 call system_clock(count_rate=rate)
 do iv=1,2
  magic=magics(iv)
  do i=1,12
   call myers_mass(nn(i),zz(i),magic,m)
   call require(m%valid,'chosen source mass has no stable local shape')
   call independent_mass(nn(i),zz(i),magic,excess128,theta128)
   maxmass=max(maxmass,abs(m%excess-excess128))
   maxtheta=max(maxtheta,abs(m%theta-theta128))
   call require(abs(m%excess-excess128)<2e-9_wp,'real128 Eq9 mass disagrees')
   call require(abs(m%theta-theta128)<2e-7_wp,'independent deformation minimum disagrees')
   call require(m%curvature>0.and.m%stationarity<2e-10_wp,'native shape is not stable stationary')
   checks=checks+1
  enddo
  call require(myers_shell_f(28,magic,.true.)==0,'shell does not vanish at magic number')
  call require(myers_shell_f(magic,magic,.false.)==0,'variant shell closure not exact')
  call require(myers_shell_f(184,magic,.false.)==0,'proton domain endpoint is not closure')
  call require(myers_shell_f(258,magic,.true.)==0,'neutron illustrative closure not exact')
  call build_myers_catalog(magic,cat,status);call require(status=='ok','catalogue construction failed')
  call build_myers_catalog(magic,small,status,184,126);call require(status=='ok','inner catalogue failed')
  call require(all(cat%sn>=0).and.all(cat%sp>= -3._wp),'published particle filters violated')
  do i=1,17
   p=exp(log(7e-15_wp)+log(4.8e-4_wp/7e-15_wp)*real(i-1,wp)/16)
   call system_clock(tick0)
   call myers_state_at_pressure(cat,p,best,status);call require(status=='ok','cached Gibbs search failed')
   call system_clock(tick1);cached_seconds=cached_seconds+real(tick1-tick0,wp)/real(rate,wp)
   call myers_state_at_pressure(small,p,ss,status);call require(status=='ok','inner Gibbs search failed')
   call require(best%z==ss%z.and.best%a==ss%a,'winner changed upon finite-domain extension')
   call require(abs(best%mu_b-ss%mu_b)<1e-10_wp,'Gibbs changed upon finite-domain extension')
   ! Independent fixed-P inversion separately for EVERY candidate: no Z cache.
   call system_clock(tick0)
   bmin=huge(1._wp);ib=0
   do k=1,cat%count
    lo=log(1e-20_wp);hi=log(1e-2_wp)
    do j=1,64
     mid=(lo+hi)/2
     call bps_isotope_state(exp(mid),cat%z(k),cat%a(k),cat%binding(k),s,status)
     call require(status=='ok','brute force primitive failed')
     if(s%pressure<p)then;lo=mid;else;hi=mid;endif
    enddo
    call bps_isotope_state(exp((lo+hi)/2),cat%z(k),cat%a(k),cat%binding(k),s,status)
    if(s%mu_b<bmin)then;bmin=s%mu_b;ib=k;brute=s;endif
   enddo
   call system_clock(tick1);brute_seconds=brute_seconds+real(tick1-tick0,wp)/real(rate,wp)
   call require(ib==best%isotope,'Z cache changed Gibbs winner')
   maxmu=max(maxmu,abs(best%mu_b-brute%mu_b))
   maxn=max(maxn,abs(best%nb/brute%nb-1))
   call require(abs(best%mu_b-brute%mu_b)<2e-10_wp,'cached Gibbs differs from brute force')
   call require(abs(best%nb/brute%nb-1)<2e-10_wp,'cached density differs from brute force')
   checks=checks+1
  enddo
 enddo
 call myers_deformation_minimum(2._wp,0._wp,5._wp,theta,v,curv,resid,valid)
 call require(valid.and.abs(theta-sqrt(log(2.5_wp)))<1e-14_wp,'quadratic analytic limit failed')
 call require(abs(v-2*(1+log(2.5_wp)))<1e-14_wp,'quadratic energy analytic limit failed')
 ! Cubic runaway is not a nuclear ground state: choose the locally bound sphere.
 call myers_deformation_minimum(1._wp,.05_wp,-2._wp,theta,v,curv,resid,valid)
 call require(valid.and.theta==0.and.v== -2._wp,'cubic runaway incorrectly selected')
 ! First-order coexistence: sphere remains locally stable but deformed well lower.
 call myers_deformation_minimum(1._wp,.1_wp,1._wp-0.004_wp,theta,v,curv,resid,valid)
 call require(valid.and.theta>0.and.v<.996_wp,'metastable sphere hid lower deformed local minimum')
 call myers_deformation_minimum(-1._wp,.05_wp,2._wp,theta,v,curv,resid,valid)
 call require(.not.valid,'unbound deformation accepted')
 call myers_mass(184,126,14,m)
 call require(abs(m%shell+10.2_wp)<.05_wp,'published p50 doubly magic shell depression mismatch')
 call myers_mass(259,126,14,m);call require(.not.m%valid,'unsupported neutron magic interval accepted')
 call myers_mass(82,185,14,m);call require(.not.m%valid,'unsupported proton magic interval accepted')
 print '(a,i0)', 'MYERS_NATIVE_INDEPENDENT_CHECKS=',checks+8
 print '(a,es24.16)', 'maximum_mass_difference_MeV=',maxmass
 print '(a,es24.16)', 'maximum_theta_difference=',maxtheta
 print '(a,es24.16)', 'cache_vs_brute_mu_difference_MeV=',maxmu
 print '(a,es24.16)', 'cache_vs_brute_density_relative=',maxn
 print '(a,es24.16)', 'cached_newton_seconds_34_queries=',cached_seconds
 print '(a,es24.16)', 'brute_bisection_seconds_34_queries=',brute_seconds
 print '(a)', 'MYERS1966_NATIVE_VERIFIED'
contains
 subroutine require(ok,message)
  logical,intent(in)::ok
  character(*),intent(in)::message
  if(.not.ok)then
   print '(a)',trim(message)
   stop 2
  endif
 end subroutine
 pure real(qp) function fq(x,magic,neutron) result(value)
  integer,intent(in)::x,magic
  logical,intent(in)::neutron
  integer::magic_list(10),k,top
  real(qp)::low,high
  magic_list=[0,2,8,magic,28,50,82,126,184,258]
  top=9;if(neutron)top=10
  value=0
  do k=2,top
   if(x>magic_list(k))cycle
   low=real(magic_list(k-1),qp);high=real(magic_list(k),qp)
   value=3._qp/5*((high**(5._qp/3)-low**(5._qp/3))*(x-low)/(high-low)- &
     real(x,qp)**(5._qp/3)+low**(5._qp/3))
   return
  enddo
 end function
 subroutine independent_mass(n,z,magic,mass,theta)
  integer,intent(in)::n,z,magic
  real(wp),intent(out)::mass,theta
  real(qp)::a,a13,c1,c2,shell,alpha0,pairing,lowest,chosen,sl,s0,sr,t,l,r,u,w,fu,fw,gold
  integer::k,it
  a=real(n+z,qp);a13=a**(1._qp/3)
  c1=15.677_qp*(1-1.79_qp*((n-z)/a)**2)
  c2=18.56_qp*(1-1.79_qp*((n-z)/a)**2)
  shell=5.8_qp*((fq(n,magic,.true.)+fq(z,magic,.false.))/(a/2)**(2._qp/3)-.26_qp*a13)
  alpha0=sqrt(5._qp)*.27_qp/a13
  qa13=a13;qc2=c2;qshell=shell;qalpha0=alpha0;qz=real(z,qp)
  lowest=huge(1._qp);chosen=0
  if(potential(0._qp)<=potential(.000001_qp))lowest=potential(0._qp)
  gold=(sqrt(5._qp)-1)/2
  ! Search for LOCAL minima only: the large-theta cubic runaway is excluded.
  sl=potential(0._qp);s0=potential(.001_qp)
  do k=1,5999
   t=real(k,qp)*.001_qp;sr=potential(t+.001_qp)
   if(s0<sl.and.s0<sr)then
    l=t-.001_qp;r=t+.001_qp
    u=r-gold*(r-l);w=l+gold*(r-l);fu=potential(u);fw=potential(w)
    do it=1,150
     if(fu<fw)then
      r=w;w=u;fw=fu;u=r-gold*(r-l);fu=potential(u)
     else
      l=u;u=w;fu=fw;w=l+gold*(r-l);fw=potential(w)
     endif
    enddo
    if(potential((l+r)/2)<lowest)then;chosen=(l+r)/2;lowest=potential(chosen);endif
   endif
   sl=s0;s0=sr
  enddo
  pairing=0
  if(mod(n,2)==0.and.mod(z,2)==0)pairing=-11._qp/sqrt(a)
  if(mod(n,2)==1.and.mod(z,2)==1)pairing=11._qp/sqrt(a)
  mass=real(8.07144_qp*n+7.28899_qp*z-c1*a-1.21129_qp*z*z/a+pairing+lowest,wp)
  theta=real(chosen,wp)
 end subroutine
  real(qp) function potential(t) result(value)
   real(qp),intent(in)::t
   real(qp)::alpha
   alpha=qalpha0*t
   value=qc2*qa13*qa13*(1+2._qp/5*alpha*alpha-4._qp/105*alpha**3)+ &
     .717_qp*qz*qz/qa13*(1-alpha*alpha/5-4._qp/105*alpha**3)+qshell*exp(-t*t)
  end function
end program
