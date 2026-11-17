! BBB source representation diagnostic; native physics, no final EOS interpolation.
! Inputs are original Table 1 rows and printed-digit halfwidths.
! method=1 PCHIP; method=2 local four-node Lagrange cubic (sensitivity control).
! Truncating x,es,en to the first nine entries gives microscopic-source-only input.
module bbb_source_model
 use compose_eos_types, only: wp
 use compose_eos_core_native, only: core_free_fermion
 implicit none
 private
 real(wp),parameter::hc=197.3269804_wp,pi=acos(-1._wp)
 public::source_state,interpolate,stencil
contains
 subroutine stencil(x,n,k,lo,hi)
  real(wp),intent(in)::x(:),n
  integer,intent(out)::k,lo,hi
  integer::i
  k=size(x)-1
  do i=1,size(x)-1
   if(n<=x(i+1))then
    k=i;exit
   endif
  enddo
  lo=max(1,k-1);hi=min(size(x),k+2)
 end subroutine
 subroutine interpolate(x,y,n,method,v,dv)
  real(wp),intent(in)::x(:),y(:),n
  integer,intent(in)::method
  real(wp),intent(out)::v,dv
  real(wp)::h(size(x)-1),d(size(x)-1),m(size(x)),t,hh,w1,w2,prod,dp
  integer::k,lo,hi,i,j,a,q,l
  l=size(x)
  call stencil(x,n,k,lo,hi)
  if(method==2)then
   lo=max(1,min(k-1,l-3));hi=lo+3;v=0;dv=0
   do i=lo,hi
    prod=1
    do j=lo,hi
     if(j/=i)prod=prod*(n-x(j))/(x(i)-x(j))
    enddo
    dp=0
    do a=lo,hi
     if(a==i)cycle
     w1=1/(x(i)-x(a))
     do q=lo,hi
      if(q/=i.and.q/=a)w1=w1*(n-x(q))/(x(i)-x(q))
     enddo
     dp=dp+w1
    enddo
    v=v+y(i)*prod;dv=dv+y(i)*dp
   enddo
   return
  endif
  h=x(2:)-x(:l-1);d=(y(2:)-y(:l-1))/h
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
  hh=h(k);t=(n-x(k))/hh
  v=(2*t**3-3*t*t+1)*y(k)+(t**3-2*t*t+t)*hh*m(k)+&
   (-2*t**3+3*t*t)*y(k+1)+(t**3-t*t)*hh*m(k+1)
  dv=(6*t*t-6*t)*(y(k)-y(k+1))/hh+(3*t*t-4*t+1)*m(k)+(3*t*t-2*t)*m(k+1)
 end subroutine
 subroutine lepton(mu,m,n,e,p,sus)
  real(wp),intent(in)::mu,m
  real(wp),intent(out)::n,e,p,sus
  real(wp)::kf,unused
  n=0;e=0;p=0;sus=0
  if(mu<=m)return
  kf=sqrt((mu-m)*(mu+m));n=(kf/hc)**3/(3*pi*pi)
  call core_free_fermion(n,m,e,p,unused)
  sus=mu*kf/(pi*pi*hc**3)
 end subroutine
 subroutine source_state(n,x,es,en,method,eps,p,y,mu,pathp,charge,status)
  real(wp),intent(in)::n,x(:),es(:),en(:)
  integer,intent(in)::method
  real(wp),intent(out)::eps,p,y,mu,pathp,charge
  character(len=*),intent(out)::status
  real(wp)::e0,d0,e1,d1,s,ds,yl,yr,ne,nm,ee,em,pe,pm,se,sm,dydn
  integer::i
  eps=0;p=0;y=0;mu=0;pathp=0;charge=0
  status='out_of_source_domain'
  if(n<x(1).or.n>x(size(x)))return
  call interpolate(x,es,n,method,e0,d0);call interpolate(x,en,n,method,e1,d1)
  s=e1-e0;ds=d1-d0
  status='nonpositive_symmetry';if(s<=0)return
  yl=0;yr=.5_wp
  do i=1,80
   y=(yl+yr)/2;mu=4*s*(1-2*y)
   call lepton(mu,.51099895_wp,ne,ee,pe,se)
   call lepton(mu,105.6583755_wp,nm,em,pm,sm)
   if(n*y-ne-nm>0)then
    yr=y
   else
    yl=y
   endif
  enddo
  eps=n*(939+e0+s*(1-2*y)**2)+ee+em
  p=n*n*(d0+ds*(1-2*y)**2)+pe+pm
  dydn=(4*(se+sm)*ds*(1-2*y)-y)/(n+8*s*(se+sm))
  pathp=p-n*n*mu*dydn;charge=n*y-ne-nm;status='ok'
 end subroutine
end module

program bbb_source_diagnostics
 use compose_eos_types, only: wp
 use compose_eos_core_native, only: core_state_t,solve_core_state
 use bbb_source_model
 use, intrinsic::ieee_arithmetic, only: ieee_value,ieee_quiet_nan
 implicit none
 real(wp)::x(16),es(16),en(16),hs(16),hn(16),py(16),pp(16),rho(16)
 real(wp)::grid(35),vals(23),a(16),b(16),n,ep,p,y,mu,pathp,q,ce,cp,cy,cm,cpath,cq
 real(wp)::te,tp,ty,tm,tpath,tq,ee,pm,ym,mum,ppm,qm,er,pr,yr,mur,ppr,qr,fd,fd2
 real(wp)::mini(3),maxi(3),agree(3),v,dv,errnode,maxcharge,maxfirst,maxfirst2,dh
 integer::i,j,k,lo,hi,l,mask,ncorner,u,out,checks
 character(len=512)::input,output,status
 type(core_state_t)::native
 call get_command_argument(1,input);call get_command_argument(2,output)
 if(len_trim(input)==0.or.len_trim(output)==0)error stop 'Usage: bbb_source_diagnostics input output'
 open(newunit=u,file=trim(input),status='old')
 do i=1,16
  read(u,*)x(i),es(i),en(i),hs(i),hn(i),py(i),pp(i),rho(i)
 enddo
 close(u)
 do i=1,16
  grid(i)=x(i)
 enddo
 do i=1,15
  grid(16+i)=(x(i)+x(i+1))/2
 enddo
 grid(32:35)=[.71_wp,.725_wp,.775_wp,.79_wp]
 open(newunit=out,file=trim(output),status='replace')
 write(out,'(a)')'n,stencil_max,epsilon,p,yp,p_local_cubic,yp_local_cubic,epsilon_truncated,p_truncated,yp_truncated,'//&
  'epsilon_corner_min,epsilon_corner_max,p_corner_min,p_corner_max,yp_corner_min,yp_corner_max,'//&
  'p_baryon_path,firstlaw_fd_residual,p_table3,yp_table3,native_p_difference,native_epsilon_difference,firstlaw_fd_residual_halfstep,corner_count'
 agree=0;maxcharge=0;maxfirst=0;maxfirst2=0;errnode=0;checks=0
 do i=1,16
  call interpolate(x,es,x(i),1,v,dv);errnode=max(errnode,abs(v-es(i)))
  call interpolate(x,en,x(i),1,v,dv);errnode=max(errnode,abs(v-en(i)))
 enddo
 do i=1,size(grid)
  n=grid(i)
  call source_state(n,x,es,en,1,ep,p,y,mu,pathp,q,status)
  if(status/='ok')error stop 'base evaluator failed'
  call solve_core_state('BBB_BHF2',n,native,status)
  if(status/='ok')error stop 'production core failed'
  agree=max(agree,[abs(ep-native%epsilon),abs(p-native%pressure),abs(y-native%yp)])
  maxcharge=max(maxcharge,abs(q))
  call source_state(n,x,es,en,2,ce,cp,cy,cm,cpath,cq,status)
  if(status/='ok')error stop 'cubic evaluator failed'
  te=ieee_value(0._wp,ieee_quiet_nan);tp=te;ty=te
  if(n<=.8_wp)then
   call source_state(n,x(:9),es(:9),en(:9),1,te,tp,ty,tm,tpath,tq,status)
   if(status/='ok')error stop 'truncated evaluator failed'
   if(n<=.7_wp.and.max(abs(ep-te),abs(p-tp),abs(y-ty))>1.e-12_wp)error stop 'unexpected truncation contamination'
  endif
  call stencil(x,n,k,lo,hi)
  l=hi-lo+1;ncorner=2**(2*l)
  mini=huge(1._wp);maxi=-huge(1._wp)
  do mask=0,ncorner-1
   a=es;b=en
   do j=0,l-1
    a(lo+j)=es(lo+j)+merge(hs(lo+j),-hs(lo+j),btest(mask,j))
    b(lo+j)=en(lo+j)+merge(hn(lo+j),-hn(lo+j),btest(mask,j+l))
   enddo
   call source_state(n,x,a,b,1,ee,pm,ym,mum,ppm,qm,status)
   if(status/='ok')error stop 'corner evaluator failed'
   mini=min(mini,[ee,pm,ym]);maxi=max(maxi,[ee,pm,ym])
  enddo
  fd=ieee_value(0._wp,ieee_quiet_nan);fd2=fd
  ! Central derivative at midpoints only, away from C1 knot breaks.
  if(i>16)then
   dh=min(n-x(k),x(k+1)-n)*1.e-3_wp
   call source_state(n-dh,x,es,en,1,ee,pm,ym,mum,ppm,qm,status)
   call source_state(n+dh,x,es,en,1,er,pr,yr,mur,ppr,qr,status)
   fd=(n*(er-ee)/(2*dh)-ep-p)/max(1._wp,abs(p))
   maxfirst=max(maxfirst,abs(fd));checks=checks+1
   dh=dh/2
   call source_state(n-dh,x,es,en,1,ee,pm,ym,mum,ppm,qm,status)
   call source_state(n+dh,x,es,en,1,er,pr,yr,mur,ppr,qr,status)
   fd2=(n*(er-ee)/(2*dh)-ep-p)/max(1._wp,abs(p))
   maxfirst2=max(maxfirst2,abs(fd2))
  endif
  vals=[n,x(hi),ep,p,y,cp,cy,te,tp,ty,mini(1),maxi(1),mini(2),maxi(2),mini(3),maxi(3),&
   pathp,fd,ieee_value(0._wp,ieee_quiet_nan),ieee_value(0._wp,ieee_quiet_nan),p-native%pressure,ep-native%epsilon,fd2]
  if(i<=16)then
   vals(19)=pp(i)*6.241509074460763_wp;vals(20)=py(i)
  endif
  write(out,'(*(es24.16,:,","))')vals,real(ncorner,wp)
 enddo
 close(out)
 call source_state(.800001_wp,x(:9),es(:9),en(:9),1,te,tp,ty,tm,tpath,tq,status)
 if(status/='out_of_source_domain')error stop 'truncated upper domain not enforced'
 call source_state(.079999_wp,x(:9),es(:9),en(:9),1,te,tp,ty,tm,tpath,tq,status)
 if(status/='out_of_source_domain')error stop 'truncated lower domain not enforced'
 print '(a,es24.16)','source_node_energy_abs_error_MeV=',errnode
 print '(a,3es24.16)','production_max_abs_epsilon_pressure_yp_difference=',agree
 print '(a,es24.16)','max_charge_residual_fm3=',maxcharge
 print '(a,i0,a,es24.16)','firstlaw_midpoint_checks=',checks,', max_scaled_residual=',maxfirst
 print '(a,es24.16)','firstlaw_halfstep_max_scaled_residual=',maxfirst2
 if(maxfirst2>maxfirst/2)error stop 'finite-difference convergence failed'
 if(errnode>1.e-12_wp.or.agree(1)>1.e-9_wp.or.agree(2)>1.e-9_wp.or.agree(3)>1.e-12_wp)error stop 'native source mismatch'
 if(maxcharge>1.e-12_wp.or.maxfirst>1.e-6_wp)error stop 'thermodynamic check failed'
 print '(a)','PASS: native source parity, printed energy nodes, charge, midpoint first law, truncated domain'
end program
