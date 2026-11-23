! Orbital Hartree-Fock solver for the NV1973 spherical Wigner-Seitz cell.
! The DME functional is evaluated in compose_eos_nv_dme. Numerics are Fortran.
! Radial finite volumes operate on R(r)=u(r)/r, with symmetric volume weighting.
module compose_eos_nv_hf
 use compose_eos_types,only:wp
 use compose_eos_nv_original,only:nv_original_source_t
 use compose_eos_nv_dme,only:nv_dme_local
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 real(wp),parameter::pi=acos(-1._wp),hc=197.3269804_wp,alpha=1._wp/137.035999084_wp
 real(wp),parameter::mn=939.56542052_wp,mp=938.27208816_wp,me=.51099895_wp
 type,public::nv_hf_config_t
  integer::neutrons=140,protons=40,mesh=360,lmax=0,radial_states=0,max_iterations=600,seed_kind=0,log_unit=0
  integer::coefficient_mode=0
  real(wp)::nb=2.79e-4_wp,mixing=.12_wp,density_tolerance=2e-6_wp,energy_tolerance=2e-7_wp
  real(wp)::cluster_radius=0._wp,cluster_neutron_ratio=1.9_wp
  logical::outer_smoothing=.false.
 end type
 type,public::nv_hf_result_t
  logical::converged=.false.,spectral_complete=.false.,outer_smoothing=.false.,energy_valid=.false.
  integer::iterations=0,mesh=0,lmax_used=0,radial_states_used=0,coefficient_mode=0
  real(wp)::radius=0,volume=0,energy=0,energy_per_baryon=0,pressure=0
  real(wp)::nuclear_energy=0,coulomb_direct=0,coulomb_exchange=0,electron_energy=0
  real(wp)::mu_n=0,mu_p=0,mu_e=0,beta_residual=0,density_residual=huge(1._wp)
  real(wp)::mu_n_plus=0,mu_p_plus=0,beta_interval_lower=0,beta_interval_upper=0
  real(wp)::energy_residual=huge(1._wp),particle_error=0
  real(wp)::tau_residual=huge(1._wp),spin_residual=huge(1._wp),selfconsistent_residual=huge(1._wp)
  real(wp)::boundary_neutron_density=0,minimum_effective_mass_coefficient=0
  real(wp)::dilation_pressure_error=0,smoothing_particle_shift(2)=0
  real(wp)::maximum_raw_density(2)=0,last_admissible_energy_per_baryon=huge(1._wp)
  real(wp),allocatable::r(:),weight(:),rho(:,:),tau(:,:),jp(:)
 end type
 public::nv_hf_solve,nv_hf_radial_spectrum,nv_hf_gradient,nv_hf_gradient_adjoint
 public::nv_hf_cell_energy,nv_hf_electrons,nv_hf_fields
contains
 subroutine nv_hf_gradient(x,h,g)
  real(wp),intent(in)::x(:),h
  real(wp),intent(out)::g(:)
  integer::i,n
  n=size(x)
  g(1)=(x(2)-x(1))/(2*h)
  do i=2,n-1
   g(i)=(x(i+1)-x(i-1))/(2*h)
  enddo
  g(n)=(x(n)-x(n-1))/(2*h)
 end subroutine
 subroutine nv_hf_gradient_adjoint(q,w,h,a)
  ! a=W^-1 D^T W q: exact derivative of sum_i w_i H(rho_i,D rho_i).
  real(wp),intent(in)::q(:),w(:),h
  real(wp),intent(out)::a(:)
  integer::i,n
  n=size(q);a=0
  a(1)=a(1)-w(1)*q(1)/(2*h);a(2)=a(2)+w(1)*q(1)/(2*h)
  do i=2,n-1
   a(i-1)=a(i-1)-w(i)*q(i)/(2*h)
   a(i+1)=a(i+1)+w(i)*q(i)/(2*h)
  enddo
  a(n-1)=a(n-1)-w(n)*q(n)/(2*h);a(n)=a(n)+w(n)*q(n)/(2*h)
  a=a/w
 end subroutine
 subroutine make_mesh(radius,n,r,w,h)
  real(wp),intent(in)::radius
  integer,intent(in)::n
  real(wp),intent(out)::r(n),w(n),h
  integer::i
  h=radius/n
  do i=1,n
   r(i)=(real(i,wp)-.5_wp)*h
   w(i)=((real(i,wp)*h)**3-((real(i,wp)-1)*h)**3)/3
  enddo
 end subroutine
 subroutine nv_hf_radial_spectrum(r,w,h,bmass,u,l,ls,wspin,k,eval,radial,status)
  real(wp),intent(in)::r(:),w(:),h,bmass(:),u(:),ls,wspin(:)
  integer,intent(in)::l,k
  real(wp),intent(out)::eval(k),radial(size(r),k)
  character(*),intent(out)::status
  real(wp)::d(size(r)),e(size(r)),z(size(r),k),work(5*size(r)),face,abstol
  integer::i,n,m,info,iw(5*size(r)),ifail(size(r))
  external::dstevx
  n=size(r);status='HF invalid radial eigenproblem'
  eval=0;radial=0
  if(n<8.or.k<1.or.k>n.or.l<0.or.h<=0)return
  if(any(bmass<=0).or..not.all(ieee_is_finite(bmass)).or..not.all(ieee_is_finite(u)))return
  d=u+bmass*real(l*(l+1),wp)*h/w+wspin*ls/r;e=0
  do i=1,n-1
   face=(real(i,wp)*h)**2*(bmass(i)+bmass(i+1))/(2*h)
   d(i)=d(i)+face/w(i);d(i+1)=d(i+1)+face/w(i+1)
   e(i)=-face/sqrt(w(i)*w(i+1))
  enddo
  ! Even l: R(Rcell)=0. Odd l: R'(Rcell)=0 (u'=u/R).
  if(mod(l,2)==0)d(n)=d(n)+2*(real(n,wp)*h)**2*bmass(n)/(h*w(n))
  abstol=2*tiny(1._wp)
  ! DSTEVX does not perform DSTEVR/ILAENV's deliberate IEEE divide-by-zero probe;
  ! therefore invalid/zero/overflow traps remain meaningful in the caller.
  call dstevx('V','I',n,d,e,0._wp,0._wp,1,k,abstol,m,eval,z,n,work,iw,ifail,info)
  if(info/=0.or.m/=k)then
   write(status,'(a,i0,a,i0)')'HF LAPACK dstevx failure info=',info,' eigenvalues=',m;return
  endif
  do i=1,k
   radial(:,i)=z(:,i)/sqrt(w)
  enddo
  status='ok'
 end subroutine
 subroutine species_density(r,w,h,bmass,u,wspin,number,lmax,nrad,proton,rho,tau,jp,mu,muplus,complete,status)
  real(wp),intent(in)::r(:),w(:),h,bmass(:),u(:),wspin(:)
  integer,intent(in)::number,lmax,nrad
  logical,intent(in)::proton
  real(wp),intent(out)::rho(:),tau(:),jp(:),mu,muplus
  logical,intent(out)::complete
  character(*),intent(out)::status
  real(wp),allocatable::rad(:,:),values(:),spins(:),degen(:),ev(:),rv(:,:)
  integer,allocatable::ell(:),irad(:),order(:)
  real(wp)::ls,occ,remain,fac,edge,tmpv
  integer::n,nch,ns,l,c,k,j,i,idx,nmax,tmpidx,nspin
  n=size(r);nch=2*lmax+1
  if(.not.proton)nch=lmax+1
  nmax=nch*nrad
  allocate(rad(n,nmax),values(nmax),spins(nmax),degen(nmax),ell(nmax),irad(nmax),order(nmax))
  allocate(ev(nrad),rv(n,nrad));ns=0
  rho=0;tau=0;jp=0;mu=0;muplus=0;complete=.false.
  do l=0,lmax
   nspin=1
   if(proton.and.l>0)nspin=2
   do c=1,nspin
    ls=0
    if(proton)then
     ls=real(l,wp)
     if(c==2)ls=-real(l+1,wp)
    endif
    call nv_hf_radial_spectrum(r,w,h,bmass,u,l,ls,wspin,nrad,ev,rv,status)
    if(status/='ok')return
    do k=1,nrad
     ns=ns+1;values(ns)=ev(k);rad(:,ns)=rv(:,k);spins(ns)=ls;ell(ns)=l;irad(ns)=k
     degen(ns)=real(2*(2*l+1),wp)
     if(proton)then
      degen(ns)=real(2*l+2,wp)
      if(c==2)degen(ns)=real(2*l,wp)
     endif
     order(ns)=ns
    enddo
   enddo
  enddo
  ! Stable insertion sort of the few hundred radial channels.
  do i=2,ns
   tmpidx=order(i);tmpv=values(tmpidx);j=i-1
   do while(j>=1)
    if(values(order(j))<=tmpv)exit
    order(j+1)=order(j);j=j-1
   enddo
   order(j+1)=tmpidx
  enddo
  remain=real(number,wp);complete=.true.
  do j=1,ns
   if(remain<=0)exit
   idx=order(j);occ=min(remain,degen(idx));remain=remain-occ;fac=occ/(4*pi)
   rho=rho+fac*rad(:,idx)**2
   tau=tau+fac*real(ell(idx)*(ell(idx)+1),wp)*h*rad(:,idx)**2/w
   do i=1,n-1
    edge=fac*(real(i,wp)*h)**2*(rad(i+1,idx)-rad(i,idx))**2/h
    tau(i)=tau(i)+edge/(2*w(i));tau(i+1)=tau(i+1)+edge/(2*w(i+1))
   enddo
   if(mod(ell(idx),2)==0)tau(n)=tau(n)+fac*2*(real(n,wp)*h)**2*rad(n,idx)**2/(h*w(n))
   if(proton)jp=jp+fac*spins(idx)*rad(:,idx)**2/r
   mu=values(idx)
   if(remain<=0)then
    muplus=mu
    if(occ>=degen(idx)-1e-12_wp.and.j<ns)muplus=values(order(j+1))
   endif
   if(ell(idx)>=lmax.or.irad(idx)>=nrad)complete=.false.
  enddo
  if(remain>1e-8_wp)then;status='HF radial basis cannot hold requested particles';return;endif
  status='ok'
 end subroutine
 pure subroutine nv_hf_electrons(ne,energy,pressure,mu)
  real(wp),intent(in)::ne
  real(wp),intent(out)::energy,pressure,mu
  real(wp)::x,pf,a
  if(ne<=0)then;energy=0;pressure=0;mu=me;return;endif
  pf=hc*(3*pi*pi*ne)**(1._wp/3);x=pf/me;mu=sqrt(pf*pf+me*me)
  if(x<.05_wp)then
   energy=me*ne*(1+3*x*x/10-3*x**4/56+x**6/48-15*x**8/1408)
  else
   a=me**4/(8*pi*pi*hc**3)
   energy=a*(x*sqrt(1+x*x)*(1+2*x*x)-asinh(x))
  endif
  pressure=ne*mu-energy
 end subroutine
 subroutine coulomb(r,w,rhop,ne,v,edirect,eexchange,vexchange)
  real(wp),intent(in)::r(:),w(:),rhop(:),ne
  real(wp),intent(out)::v(:),edirect,eexchange,vexchange(:)
  real(wp)::q(size(r)),inside,outside,ex
  integer::i
  q=4*pi*w*(rhop-ne);inside=0;outside=sum(q/r)
  do i=1,size(r)
   inside=inside+q(i);outside=outside-q(i)/r(i)
   v(i)=alpha*hc*(inside/r(i)+outside)
  enddo
  edirect=.5_wp*sum(q*v)
  ex=alpha*hc*(3/pi)**(1._wp/3)
  vexchange=-ex*max(rhop,0._wp)**(1._wp/3)
  eexchange=-.75_wp*ex*sum(4*pi*w*max(rhop,0._wp)**(4._wp/3))
 end subroutine
 subroutine nv_hf_fields(src,r,w,h,rho,tau,jp,ne,uu,bb,ss,energy,edirect,eexchange,pressure,muec,status,coefficient_mode)
  type(nv_original_source_t),intent(in)::src
  real(wp),intent(in)::r(:),w(:),h,rho(:,:),tau(:,:),jp(:),ne
  real(wp),intent(out)::uu(:,:),bb(:,:),ss(:),energy,edirect,eexchange,pressure,muec
  character(*),intent(out)::status
  integer,intent(in),optional::coefficient_mode
  real(wp)::gr(size(r),2),qq(size(r),2),adj(size(r)),v(size(r)),vx(size(r))
  real(wp)::en,ul(2),bm(2),qg(2),ws,vol,virial,ee,pe,mu,tlocal(2),glocal(2)
  integer::i,s
  do s=1,2
   call nv_hf_gradient(rho(:,s),h,gr(:,s))
  enddo
  energy=0;virial=0
  do i=1,size(r)
   tlocal=tau(i,:);glocal=gr(i,:)
   call nv_dme_local(src,rho(i,1),rho(i,2),tlocal,glocal,jp(i),en,ul,bm,qg,ws,status,coefficient_mode)
   if(status/='ok')return
   uu(i,:)=ul;bb(i,:)=bm;qq(i,:)=qg;ss(i)=ws
   energy=energy+4*pi*w(i)*en
   virial=virial+4*pi*w(i)*(3*sum(rho(i,:)*ul)+5*sum(tau(i,:)*bm)+4*sum(gr(i,:)*qg)+4*jp(i)*ws-3*en)
  enddo
  do s=1,2
   call nv_hf_gradient_adjoint(qq(:,s),w,h,adj)
   uu(:,s)=uu(:,s)+adj
  enddo
  call coulomb(r,w,rho(:,2),ne,v,edirect,eexchange,vx)
  uu(:,2)=uu(:,2)+v+vx
  vol=sum(4*pi*w);muec=-sum(4*pi*w*v)/vol
  call nv_hf_electrons(ne,ee,pe,mu)
  pressure=(virial+edirect+eexchange)/(3*vol)+pe
  if(.not.all(ieee_is_finite(uu)).or..not.all(ieee_is_finite(bb)).or.any(bb<=0))then
   status='HF nonfinite field or nonpositive effective mass';return
  endif
  status='ok'
 end subroutine
 subroutine nv_hf_cell_energy(src,r,w,h,rho,tau,jp,protons,energy,pressure,status,coefficient_mode)
  type(nv_original_source_t),intent(in)::src
  real(wp),intent(in)::r(:),w(:),h,rho(:,:),tau(:,:),jp(:)
  integer,intent(in)::protons
  real(wp),intent(out)::energy,pressure
  character(*),intent(out)::status
  integer,intent(in),optional::coefficient_mode
  real(wp)::u(size(r),2),b(size(r),2),ss(size(r)),ec,ex,ee,pe,mu,muec,vol
  vol=sum(4*pi*w)
  call nv_hf_fields(src,r,w,h,rho,tau,jp,real(protons,wp)/vol,u,b,ss,energy,ec,ex,pressure,muec,status,coefficient_mode)
  if(status/='ok')return
  call nv_hf_electrons(real(protons,wp)/vol,ee,pe,mu)
  energy=energy+ec+ex+vol*ee+real(protons,wp)*(mp-mn)
 end subroutine
 subroutine smooth_outer(r,w,h,rho,tau,number,shift)
  ! Optional implementation of the authors' nonvariational boundary prescription.
  ! Gas interval is explicit: r in [max(0.6R,R-6/kF),R-1/kF].
  ! Renormalization is recorded; raw orbital fields are always retained separately.
  real(wp),intent(in)::r(:),w(:),h
  real(wp),intent(inout)::rho(:,:),tau(:,:)
  integer,intent(in)::number(2)
  real(wp),intent(out)::shift(2)
  real(wp)::rad,kf,delta,low,avgr,avgt,tot,fac
  logical::mask(size(r)),outer(size(r))
  integer::s
  rad=r(size(r))+.5_wp*h;kf=(3*pi*pi*max(rho(size(r),1),1e-12_wp))**(1._wp/3)
  delta=min(1/kf,.2_wp*rad);low=max(.6_wp*rad,rad-6*delta)
  mask=r>=low.and.r<rad-delta;outer=r>=rad-delta;shift=0
  if(count(mask)<3)return
  do s=1,2
   avgr=sum(w*rho(:,s),mask)/sum(w,mask);avgt=sum(w*tau(:,s),mask)/sum(w,mask)
   where(outer)rho(:,s)=avgr
   where(outer)tau(:,s)=avgt
   tot=sum(4*pi*w*rho(:,s));shift(s)=tot-real(number(s),wp)
   if(tot>0)then
    fac=real(number(s),wp)/tot;rho(:,s)=rho(:,s)*fac;tau(:,s)=tau(:,s)*fac
   endif
  enddo
 end subroutine
 subroutine nv_hf_solve(src,cfg,res,status,initial)
  type(nv_original_source_t),intent(in)::src
  type(nv_hf_config_t),intent(in)::cfg
  type(nv_hf_result_t),intent(out)::res
  character(*),intent(out)::status
  type(nv_hf_result_t),intent(in),optional::initial
  real(wp),allocatable::r(:),w(:),rho(:,:),tau(:,:),jp(:),rr(:,:),tt(:,:),jj(:),dummy(:)
  real(wp),allocatable::uu(:,:),bb(:,:),ss(:),shape(:),rawrho(:,:),rawtau(:,:),rawjp(:)
  real(wp)::radius,vol,h,kf,rc,bound,fac,mu(2),muplus(2),en,ed,ex,p,muec,ee,pe,mue
  real(wp)::prev_energy,epart,dr,dt,dj,denom,alpha_mix,oldres,shift(2),epsd,eplus,eminus,ptest,pdilate
  integer::n,lmax,nrad,it,s,numbers(2),a
  logical::complete(2)
  res=nv_hf_result_t();status='HF invalid configuration'
  if(.not.src%valid)return
  if(cfg%neutrons<1.or.cfg%protons<1.or.cfg%neutrons>100000.or.cfg%protons>100000)return
  if(cfg%mesh<32.or.cfg%mesh>4096.or.cfg%max_iterations<1)return
  if(cfg%lmax<0.or.cfg%radial_states<0)return
  if(.not.all(ieee_is_finite([cfg%nb,cfg%mixing,cfg%density_tolerance,cfg%energy_tolerance, &
       cfg%cluster_radius,cfg%cluster_neutron_ratio])))return
  if(cfg%nb<=0.or.cfg%nb>.1_wp.or.cfg%mixing<=0.or.cfg%mixing>1)return
  if(cfg%density_tolerance<=0.or.cfg%energy_tolerance<=0)return
  if(cfg%cluster_radius<0.or.cfg%cluster_neutron_ratio<0.or.cfg%seed_kind<0.or.cfg%seed_kind>2)return
  if(cfg%coefficient_mode<0.or.cfg%coefficient_mode>1)return
  n=cfg%mesh;a=cfg%neutrons+cfg%protons;numbers=[cfg%neutrons,cfg%protons]
  vol=real(a,wp)/cfg%nb;radius=(3*vol/(4*pi))**(1._wp/3)
  kf=(3*pi*pi*cfg%nb)**(1._wp/3)
  lmax=cfg%lmax
  if(lmax<=0)lmax=max(10,ceiling(kf*radius)+5)
  nrad=cfg%radial_states
  if(nrad<=0)nrad=min(n,max(8,ceiling(kf*radius/pi)+7))
  if(lmax>100.or.nrad>n.or.nrad<2)return
  allocate(r(n),w(n),rho(n,2),tau(n,2),jp(n),rr(n,2),tt(n,2),jj(n),dummy(n))
  allocate(uu(n,2),bb(n,2),ss(n),shape(n),rawrho(n,2),rawtau(n,2),rawjp(n))
  call make_mesh(radius,n,r,w,h)
  bound=min(real(cfg%neutrons,wp),cfg%cluster_neutron_ratio*cfg%protons)
  rc=cfg%cluster_radius
  if(rc<=0)rc=(3*(bound+cfg%protons)/(4*pi*.14_wp))**(1._wp/3)
  shape=1/(1+exp(min(700._wp,(r-rc)/.65_wp)));fac=sum(4*pi*w*shape)
  if(maxval(shape)/fac>1/vol)bound=min(bound,.98_wp*(.145_wp-real(cfg%neutrons,wp)/vol)/(maxval(shape)/fac-1/vol))
  bound=max(0._wp,bound)
  rho(:,2)=real(cfg%protons,wp)*shape/fac
  rho(:,1)=bound*shape/fac+(cfg%neutrons-bound)/vol
  if(cfg%seed_kind==1)then
   rho(:,1)=real(cfg%neutrons,wp)/vol;rho(:,2)=real(cfg%protons,wp)/vol
  else if(cfg%seed_kind==2)then
   shape=1/(1+exp(min(700._wp,(r-1.35_wp*rc)/1.2_wp)));fac=sum(4*pi*w*shape)
   rho(:,2)=real(cfg%protons,wp)*shape/fac
   rho(:,1)=bound*shape/fac+(cfg%neutrons-bound)/vol
  endif
  do s=1,2
   tau(:,s)=.6_wp*(3*pi*pi)**(2._wp/3)*rho(:,s)**(5._wp/3)
  enddo
  jp=0
  if(present(initial))then
   if(allocated(initial%rho).and.allocated(initial%tau).and.allocated(initial%jp))then
    if(size(initial%rho,1)==n.and.size(initial%rho,2)==2.and.size(initial%tau,1)==n.and. &
       size(initial%tau,2)==2.and.size(initial%jp)==n.and.ieee_is_finite(initial%radius).and.initial%radius>0)then
     fac=radius/initial%radius
     rho=initial%rho/fac**3;tau=initial%tau/fac**5;jp=initial%jp/fac**4
     do s=1,2
      fac=real(numbers(s),wp)/sum(4*pi*w*rho(:,s));rho(:,s)=rho(:,s)*fac;tau(:,s)=tau(:,s)*fac
      if(s==2)jp=jp*fac
     enddo
    endif
   endif
  endif
  alpha_mix=cfg%mixing;prev_energy=huge(1._wp);oldres=huge(1._wp);shift=0
  res%radius=radius;res%volume=vol;res%mesh=n;res%lmax_used=lmax;res%radial_states_used=nrad
  res%outer_smoothing=cfg%outer_smoothing
  res%coefficient_mode=cfg%coefficient_mode
  do it=1,cfg%max_iterations
   call nv_hf_fields(src,r,w,h,rho,tau,jp,real(cfg%protons,wp)/vol,uu,bb,ss,en,ed,ex,p,muec,status,cfg%coefficient_mode)
   if(status/='ok')return
   call species_density(r,w,h,bb(:,1),uu(:,1),0*ss,cfg%neutrons,lmax,nrad,.false.,rr(:,1),tt(:,1),dummy,mu(1),muplus(1),complete(1),status)
   if(status/='ok')return
   call species_density(r,w,h,bb(:,2),uu(:,2),ss,cfg%protons,lmax,nrad,.true.,rr(:,2),tt(:,2),jj,mu(2),muplus(2),complete(2),status)
   if(status/='ok')return
   rawrho=rr;rawtau=tt;rawjp=jj
   res%maximum_raw_density=[maxval(rr(:,1)),maxval(rr(:,2))]
   if(cfg%outer_smoothing)call smooth_outer(r,w,h,rr,tt,numbers,shift)
   dr=0;dt=0
   do s=1,2
    denom=real(numbers(s),wp)
    dr=max(dr,sum(4*pi*w*abs(rr(:,s)-rho(:,s)))/denom)
    dt=max(dt,sum(4*pi*w*abs(tt(:,s)-tau(:,s)))/max(sum(4*pi*w*tt(:,s)),1e-12_wp))
   enddo
   dj=sum(4*pi*w*abs(jj-jp))/max(sum(4*pi*w*abs(jj)),cfg%protons/radius)
   epart=huge(1._wp)/4
   if(maxval(rawrho)<=.15_wp)then
    call nv_hf_cell_energy(src,r,w,h,rawrho,rawtau,rawjp,cfg%protons,epart,ptest,status,cfg%coefficient_mode)
    if(status/='ok')return
    epart=epart/real(a,wp)
    res%last_admissible_energy_per_baryon=epart
   endif
   if(cfg%log_unit>0.and.(it<10.or.mod(it,20)==0))write(cfg%log_unit,'(a,i0,7(1x,es12.4))') &
    'HFITER ',it,epart,dr,dt,dj,maxval(rawrho(:,1)),maxval(rawrho(:,2)),alpha_mix
   res%energy_residual=abs(epart-prev_energy);res%density_residual=dr;res%iterations=it
   res%tau_residual=dt;res%spin_residual=dj
   if(it>4.and.max(dr,dt,dj)<cfg%density_tolerance.and.res%energy_residual<cfg%energy_tolerance.and.all(complete))then
    res%converged=.true.;exit
   endif
   ! Reduce mixing only after observed oscillatory growth; never relax tolerances.
   if(it>8.and.dr>1.25_wp*oldres)alpha_mix=max(.002_wp,.8_wp*alpha_mix)
   ! Stay inside the original coefficient fit domain while damping trial orbitals.
   do while(maxval((1-alpha_mix)*rho+alpha_mix*rr)>.1499_wp)
    alpha_mix=.5_wp*alpha_mix
    if(alpha_mix<1e-6_wp)then;status='HF_SCF_DOMAIN_STEP_TOO_SMALL';return;endif
   enddo
   rho=(1-alpha_mix)*rho+alpha_mix*rr;tau=(1-alpha_mix)*tau+alpha_mix*tt;jp=(1-alpha_mix)*jp+alpha_mix*jj
   prev_energy=epart;oldres=dr
  enddo
  res%spectral_complete=all(complete);res%smoothing_particle_shift=shift
  res%rho=rawrho;res%tau=rawtau;res%jp=rawjp;res%r=r;res%weight=w
  call nv_hf_fields(src,r,w,h,rawrho,rawtau,rawjp,real(cfg%protons,wp)/vol,uu,bb,ss,en,ed,ex,p,muec,status,cfg%coefficient_mode)
  if(status/='ok')then;res%converged=.false.;return;endif
  call nv_hf_electrons(real(cfg%protons,wp)/vol,ee,pe,mue)
  res%nuclear_energy=en;res%coulomb_direct=ed;res%coulomb_exchange=ex;res%electron_energy=ee*vol
  res%energy=en+ed+ex+ee*vol+cfg%protons*(mp-mn);res%energy_per_baryon=res%energy/a
  res%energy_valid=.true.
  res%pressure=p;res%mu_n=mu(1);res%mu_p=mu(2);res%mu_e=mue+muec
  res%beta_residual=mu(2)-mu(1)+res%mu_e+mp-mn
  res%mu_n_plus=muplus(1);res%mu_p_plus=muplus(2)
  res%beta_interval_lower=mu(2)-muplus(1)+res%mu_e+mp-mn
  res%beta_interval_upper=muplus(2)-mu(1)+res%mu_e+mp-mn
  res%particle_error=maxval(abs([sum(4*pi*w*rawrho(:,1)),sum(4*pi*w*rawrho(:,2))]-real(numbers,wp)))
  res%boundary_neutron_density=rawrho(n,1);res%minimum_effective_mass_coefficient=minval(bb)
  if(res%converged)then
   call species_density(r,w,h,bb(:,1),uu(:,1),0*ss,cfg%neutrons,lmax,nrad,.false.,rr(:,1),tt(:,1),dummy,mu(1),muplus(1),complete(1),status)
   if(status/='ok')then;res%converged=.false.;return;endif
   call species_density(r,w,h,bb(:,2),uu(:,2),ss,cfg%protons,lmax,nrad,.true.,rr(:,2),tt(:,2),jj,mu(2),muplus(2),complete(2),status)
   if(status/='ok')then;res%converged=.false.;return;endif
   dr=0;dt=0
   do s=1,2
    dr=max(dr,sum(4*pi*w*abs(rr(:,s)-rawrho(:,s)))/real(numbers(s),wp))
    dt=max(dt,sum(4*pi*w*abs(tt(:,s)-rawtau(:,s)))/max(sum(4*pi*w*rawtau(:,s)),1e-12_wp))
   enddo
   dj=sum(4*pi*w*abs(jj-rawjp))/max(sum(4*pi*w*abs(rawjp)),cfg%protons/radius)
   res%selfconsistent_residual=max(dr,dt,dj)
   if(res%selfconsistent_residual>cfg%density_tolerance)res%converged=.false.
   res%spectral_complete=all(complete)
   if(.not.res%spectral_complete)res%converged=.false.
   ! Report frontier levels from the same final raw fields used by the replay.
   res%mu_n=mu(1);res%mu_p=mu(2);res%mu_n_plus=muplus(1);res%mu_p_plus=muplus(2)
   res%beta_residual=mu(2)-mu(1)+res%mu_e+mp-mn
   res%beta_interval_lower=mu(2)-muplus(1)+res%mu_e+mp-mn
   res%beta_interval_upper=muplus(2)-mu(1)+res%mu_e+mp-mn
  endif
  epsd=1e-5_wp
  fac=exp(epsd)
  call nv_hf_cell_energy(src,r*fac,w*fac**3,h*fac,rawrho/fac**3,rawtau/fac**5,rawjp/fac**4,cfg%protons,eplus,ptest,status,cfg%coefficient_mode)
  if(status/='ok')then;res%converged=.false.;return;endif
  fac=exp(-epsd)
  call nv_hf_cell_energy(src,r*fac,w*fac**3,h*fac,rawrho/fac**3,rawtau/fac**5,rawjp/fac**4,cfg%protons,eminus,ptest,status,cfg%coefficient_mode)
  if(status/='ok')then;res%converged=.false.;return;endif
  pdilate=-(eplus-eminus)/(2*epsd*3*vol)
  res%dilation_pressure_error=abs(pdilate-p)
  if(res%converged)then;status='ok';else;status='HF_ITERATION_NOT_CONVERGED';endif
 end subroutine
end module
