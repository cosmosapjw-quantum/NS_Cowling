! Parameter-injected, cold uniform RMF evaluation using unchanged TF kernels.
! Units: fields/masses MeV; densities fm^-3; epsilon/P MeV fm^-3.
! rho vertex gr*tau3/2, tau3(p)=+1; M_i*=M_i-gs*phi.
program rmf_candidate_driver
 use compose_eos_types, only: wp
 use compose_eos_tf_models
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 type(tf_model_t) :: p
 type :: state_t
  real(wp) :: n=0,yp=0,ye=0,ymu=0,eps=0,pres=0,mun=0,mup=0,mue=0,mub=0
  real(wp) :: phi=0,w=0,r=0,mn=0,mp=0,sr=0,ns=0,fs=0,fw=0,fr=0,curv=0
  real(wp) :: beta=0,charge=0,euler=0,epsb=0,pb=0
  integer :: roots=0,branch=0
 end type
 type(state_t) :: s
 character(2048) :: parfile,infile,outfile,mode,arg
 integer :: iu,ou,ios,k,rows
 real(wp) :: n,yp
 namelist /candidate/ p
 if(command_argument_count()<4)call fail('usage: driver params.nml density.txt states.csv beta|baryons [Yp]')
 call get_command_argument(1,parfile)
 call get_command_argument(2,infile)
 call get_command_argument(3,outfile)
 call get_command_argument(4,mode)
 yp=.5_wp
 if(trim(mode)=='baryons')then
  if(command_argument_count()>=5)then
   call get_command_argument(5,arg);read(arg,*,iostat=ios)yp
   if(ios/=0)call fail('invalid Yp')
  endif
 else if(trim(mode)/='beta')then
  call fail('mode must be beta or baryons')
 endif
 if(.not.ieee_is_finite(yp).or.yp<0.or.yp>1)call fail('invalid Yp')
 open(newunit=iu,file=trim(parfile),status='old',iostat=ios)
 if(ios/=0)call fail('cannot open parameters')
 read(iu,nml=candidate,iostat=ios);close(iu)
 if(ios/=0)call fail('invalid parameter namelist')
 if(any(.not.ieee_is_finite([p%mn,p%mp,p%ms,p%mw,p%mr,p%gs,p%gw,p%gr,p%a3,p%a4,p%b4,p%cross])))&
  call fail('nonfinite parameters')
 if(min(p%mn,p%mp,p%ms,p%mw,p%mr,p%gs,p%gw,p%gr)<=0)call fail('masses and couplings must be positive')
 if(p%b4<0.or.p%cross/=0)call fail('driver supports nonnegative omega quartic and zero cross coupling only')
 if(p%density_dependent)then
  if(.not.ieee_is_finite(p%n0).or.p%n0<=0)call fail('DD model requires positive n0')
 endif
 open(newunit=iu,file=trim(infile),status='old',iostat=ios)
 if(ios/=0)call fail('cannot open densities')
 open(newunit=ou,file=trim(outfile),status='replace',iostat=ios)
 if(ios/=0)call fail('cannot open output')
 write(ou,'(a)') 'n_b,Y_p,Y_e,Y_mu,epsilon,pressure,mu_B,mu_n,mu_p,mu_e,phi,omega,rho,'//&
  'effective_mass_n,effective_mass_p,rearrangement,scalar_density,scalar_residual_fm3,'//&
  'omega_residual_fm3,rho_residual_fm3,scalar_curvature,beta_residual_mev,charge_residual_fm3,'//&
  'euler_residual_mevfm3,epsilon_baryons,pressure_baryons,scalar_root_count,composition_branch'
 rows=0
 do
  read(iu,*,iostat=ios)n
  if(ios<0)exit
  if(ios/=0)call fail('invalid density input')
  if(.not.ieee_is_finite(n).or.n<=0)call fail('density must be finite and positive')
  if(trim(mode)=='beta')then
   call beta_state(n,s)
  else
   call baryons(n,yp,s);s%mub=(1-yp)*s%mun+yp*s%mup
   s%euler=s%eps+s%pres-n*s%mub
  endif
  write(ou,'(*(es25.16e3,:,","))') s%n,s%yp,s%ye,s%ymu,s%eps,s%pres,s%mub,s%mun,s%mup,s%mue,&
   s%phi,s%w,s%r,s%mn,s%mp,s%sr,s%ns,s%fs,s%fw,s%fr,s%curv,s%beta,s%charge,s%euler,&
   s%epsb,s%pb,real(s%roots,wp),real(s%branch,wp)
  rows=rows+1
 enddo
 close(iu);close(ou)
 if(rows==0)call fail('empty density input')
 write(*,'(a,i0,a)') '{"status":"ok","rows":',rows,'}'
contains
 subroutine fail(message)
 character(*),intent(in)::message
 write(error_unit,'(a)')trim(message)
 error stop 1
 end subroutine

 subroutine scalar_terms(n,y,gs,phi,f,curv,energy,ns)
 real(wp),intent(in)::n,y,gs,phi
 real(wp),intent(out)::f,curv,energy,ns
 real(wp)::en,ep,sn,sp,dn,dp,ef,u,du,ddu
 call tf_fermion(n*(1-y),p%mn-gs*phi,en,sn,dn,ef)
 call tf_fermion(n*y,p%mp-gs*phi,ep,sp,dp,ef)
 call tf_scalar_potential(p,phi,u,du,ddu)
 ns=sn+sp;f=du/tf_hc3-gs*ns;curv=ddu/tf_hc3+gs*gs*(dn+dp)
 energy=en+ep+u/tf_hc3
 end subroutine

 subroutine baryons(n,y,s)
 real(wp),intent(in)::n,y
 type(state_t),intent(out)::s
 real(wp)::gs,gw,gr,dgs,dgw,dgr,nn,np,n3,cap,left,right,fl,fr,mid,f,curv,energy,ns,best
 real(wp)::phi,en,ep,sn,sp,dn,dp,efn,efp,u,du,ddu,v,dw,dr,dww,drr,dwr
 real(wp)::wl,wr,wm
 integer::j,it
 s=state_t();s%n=n;s%yp=y
 call tf_couplings(p,n,gs,gw,gr,dgs,dgw,dgr)
 if(.not.all(ieee_is_finite([gs,gw,gr,dgs,dgw,dgr])).or.min(gs,gw,gr)<=0)call fail('invalid density-dependent couplings')
 nn=n*(1-y);np=n*y;n3=np-nn
 cap=min(p%mn,p%mp)/gs*(1-1.e-10_wp)
 left=0;call scalar_terms(n,y,gs,left,fl,curv,energy,ns)
 best=huge(1._wp);phi=0
 ! Isolate every upward root in the positive-mass interval and select the
 ! lowest scalar+kinetic energy among locally stable stationary points.
 do j=1,256
  right=cap*real(j,wp)/256
  call scalar_terms(n,y,gs,right,fr,curv,energy,ns)
  if(fl<=0.and.fr>=0)then
   wl=left;wr=right
   do it=1,90
    mid=(wl+wr)/2
    call scalar_terms(n,y,gs,mid,f,curv,energy,ns)
    if(f>0)then
     wr=mid
    else
     wl=mid
    endif
    if(wr-wl<2.e-14_wp*max(1._wp,abs(mid)))exit
   enddo
   mid=(wl+wr)/2;call scalar_terms(n,y,gs,mid,f,curv,energy,ns)
   if(curv>0.and.ieee_is_finite(energy))then
    s%roots=s%roots+1
    if(energy<best)then
     phi=mid;best=energy
    endif
   endif
  endif
  left=right;fl=fr
 enddo
 if(s%roots==0)call fail('no locally stable positive-mass scalar root')
 s%phi=phi;s%mn=p%mn-gs*phi;s%mp=p%mp-gs*phi
 call scalar_terms(n,y,gs,phi,s%fs,s%curv,energy,s%ns)
 wl=0;wr=gw*n*tf_hc3/p%mw**2
 if(p%b4==0)then
  s%w=wr
 else
  do it=1,90
   wm=(wl+wr)/2
   if(p%mw**2*wm+p%b4*wm**3>gw*n*tf_hc3)then
    wr=wm
   else
    wl=wm
   endif
   if(wr-wl<2.e-14_wp*max(1._wp,abs(wm)))exit
  enddo
  s%w=(wl+wr)/2
 endif
 s%r=gr*n3*tf_hc3/(2*p%mr**2)
 call tf_fermion(nn,s%mn,en,sn,dn,efn)
 call tf_fermion(np,s%mp,ep,sp,dp,efp)
 call tf_scalar_potential(p,phi,u,du,ddu)
 call tf_vector_potential(p,s%w,s%r,v,dw,dr,dww,drr,dwr)
 s%fw=dw/tf_hc3-gw*n;s%fr=dr/tf_hc3-gr*n3/2
 s%sr=dgw*s%w*n+dgr*s%r*n3/2-dgs*phi*(sn+sp)
 s%mun=efn+gw*s%w-gr*s%r/2+s%sr
 s%mup=efp+gw*s%w+gr*s%r/2+s%sr
 s%eps=en+ep+(u-v)/tf_hc3+gw*s%w*n+gr*s%r*n3/2
 s%pres=nn*efn+np*efp-en-ep+(v-u)/tf_hc3+n*s%sr
 s%epsb=s%eps;s%pb=s%pres
 end subroutine

 pure function number_at_mu(mu,mass) result(number)
 real(wp),intent(in)::mu,mass
 real(wp)::number
 number=max(0._wp,(mu-mass)*(mu+mass))**1.5_wp/(3*tf_pi**2*tf_hc3)
 end function

 subroutine lepton_at_charge(charge,mu,ne,nmu,energy,pressure)
 real(wp),intent(in)::charge
 real(wp),intent(out)::mu,ne,nmu,energy,pressure
 real(wp)::lo,hi,mid,ee,em,scalar,ds,ef
 integer::it
 lo=tf_me;hi=sqrt((tf_hbarc*(3*tf_pi**2*charge)**(1._wp/3))**2+tf_me**2)
 if(hi>tf_mmu)then
  do it=1,90
   mid=(lo+hi)/2
   if(number_at_mu(mid,tf_me)+number_at_mu(mid,tf_mmu)>charge)then
    hi=mid
   else
    lo=mid
   endif
   if(hi-lo<2.e-13_wp)exit
  enddo
  mu=(lo+hi)/2
 else
  mu=hi
 endif
 ne=number_at_mu(mu,tf_me);nmu=number_at_mu(mu,tf_mmu)
 call tf_fermion(ne,tf_me,ee,scalar,ds,ef)
 call tf_fermion(nmu,tf_mmu,em,scalar,ds,ef)
 energy=ee+em;pressure=mu*(ne+nmu)-energy
 end subroutine

 subroutine beta_probe(n,y,s,f)
 real(wp),intent(in)::n,y
 type(state_t),intent(out)::s
 real(wp),intent(out)::f
 real(wp)::ne,nmu,el,pl
 call baryons(n,y,s)
 call lepton_at_charge(n*y,s%mue,ne,nmu,el,pl)
 s%ye=ne/n;s%ymu=nmu/n
 s%eps=s%eps+el;s%pres=s%pres+pl
 f=s%mun-s%mup-s%mue;s%beta=f;s%charge=n*y-ne-nmu
 end subroutine

 subroutine beta_state(n,s)
 real(wp),intent(in)::n
 type(state_t),intent(out)::s
 type(state_t)::s0,s1
 real(wp)::f0,f1,lo,hi,y,f
 integer::it,branch
 call beta_probe(n,0._wp,s0,f0);call beta_probe(n,1._wp,s1,f1)
 if(f0<=0)then
  s=s0;branch=-1
 else if(f1>=0)then
  s=s1;branch=1
 else
  lo=0;hi=1
  do it=1,90
   y=(lo+hi)/2;call beta_probe(n,y,s,f)
   if(f>0)then
    lo=y
   else
    hi=y
   endif
   if(hi-lo<2.e-15_wp)exit
  enddo
  call beta_probe(n,(lo+hi)/2,s,f);branch=0
 endif
 s%branch=branch;s%mub=s%mun
 if(branch==1)s%mub=s%mup+s%mue
 s%euler=s%eps+s%pres-n*s%mub
 if(.not.all(ieee_is_finite([s%eps,s%pres,s%mub,s%euler])))call fail('nonfinite state')
 end subroutine
end program
