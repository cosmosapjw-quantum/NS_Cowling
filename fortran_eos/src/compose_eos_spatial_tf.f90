! Spatial finite-range RMF TF with native screened electrons and SLSQP.
! Physical equations and stopping tests preserve the declared Python reference.
module compose_eos_spatial_tf
 use compose_eos_types, only: wp
 use compose_eos_tf_models
 use compose_eos_optimize_native, only: native_optimizer_result_t,minimize_slsqp_native
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 integer,parameter::hp=selected_real_kind(18,300)
 type,public::tf_grid_t
  integer::points=0,dimension=0
  real(wp),allocatable::x(:),weights(:),sqrtw(:),lap(:,:),face(:)
 end type
 type,public::tf_state_t
  character(32)::model_id='',seed='',numerical_status='NOT_RUN',geometry=''
  character(160)::evaluation_status='not evaluated'
  integer::points=0,dimension=0,iterations=0,evaluations=0,optimizer_mode=0
  logical::optimizer_success=.false.,radius_at_search_bound=.false.,density_at_search_bound=.false.,muon_threshold_crossed=.false.
  real(wp)::n_b=0,radius_fm=0,epsilon=0,pressure=0,mu_B=0,mu_Q_lepton=0,energy_per_baryon=0
  real(wp)::Y_n=0,Y_p=0,Y_e=0,Y_mu=0
  real(wp)::scalar_pde_residual=huge(1._wp),vector_pde_residual=huge(1._wp),electron_pde_residual=huge(1._wp)
  real(wp)::mean_baryon_density=0,mean_charge_density=0,baryon_constraint_error=huge(1._wp),charge_constraint_error=huge(1._wp)
  real(wp)::d_epsilon_d_log_radius=0,max_chemical_KKT_residual_MeV=huge(1._wp),radius_stationarity_MeV=huge(1._wp)
  real(wp)::baryon_density_contrast=0,gradient_sigma=0,gradient_omega=0,gradient_rho=0,coulomb_energy=0
  real(wp),allocatable::nn(:),np(:),ne(:),phi(:),omega(:),rho(:),coulomb(:),mu_n(:),mu_p(:),mu_e(:)
 end type
 public::make_tf_grid,tf_poisson_native,screened_electrons_native,evaluate_spatial_tf_profile
 public::evaluate_screened_tf_profile,solve_spatial_tf,initial_tf_profile,tf_geometry_label,select_tf_candidates_native
 public::solve_tf_family_native
contains
 subroutine make_tf_grid(n,d,g,status)
 integer,intent(in)::n,d
 type(tf_grid_t),intent(out)::g
 character(*),intent(out)::status
 integer::i
 real(wp)::edge
 status='grid requires N>=8 and dimension 1,2,3'
 if(n<8.or.d<1.or.d>3)return
 g%points=n;g%dimension=d
 allocate(g%x(n),g%weights(n),g%sqrtw(n),g%lap(n,n),g%face(n-1));g%lap=0
 do i=1,n
  g%x(i)=(i-.5_wp)/n
  g%weights(i)=(real(i,wp)/n)**d-(real(i-1,wp)/n)**d
 enddo
 g%sqrtw=sqrt(g%weights)
 do i=1,n-1
  edge=real(i,wp)/n;g%face(i)=d*edge**(d-1)*n
  g%lap(i,i)=g%lap(i,i)+g%face(i)/g%weights(i)
  g%lap(i,i+1)=-g%face(i)/g%weights(i)
  g%lap(i+1,i+1)=g%lap(i+1,i+1)+g%face(i)/g%weights(i+1)
  g%lap(i+1,i)=-g%face(i)/g%weights(i+1)
 enddo
 status='ok'
 end subroutine

 subroutine tf_poisson_native(g,charge,radius,v)
 type(tf_grid_t),intent(in)::g
 real(wp),intent(in)::charge(:),radius
 real(wp),intent(out)::v(:)
 real(wp)::q(g%points),flux(g%points-1),enclosed
 integer::i
 q=charge-sum(g%weights*charge);enclosed=0
 do i=1,g%points-1
  enclosed=enclosed+g%weights(i)*q(i)
  flux(i)=4*tf_pi*tf_e2*radius**2*enclosed
 enddo
 v(g%points)=0
 do i=g%points-1,1,-1
  v(i)=v(i+1)+flux(i)/g%face(i)
 enddo
 v=v-sum(g%weights*v)
 end subroutine

 subroutine band_solve(a,b,width,x,ok)
 real(wp),intent(in)::a(:,:),b(:,:)
 integer,intent(in)::width
 real(wp),intent(out)::x(:,:)
 logical,intent(out)::ok
 real(wp),allocatable::ab(:,:)
 integer,allocatable::ipiv(:)
 integer::n,nrhs,ldab,i,j,info
 n=size(a,1);nrhs=size(b,2);ldab=3*width+1
 allocate(ab(ldab,n),ipiv(n));ab=0
 do j=1,n
  do i=max(1,j-width),min(n,j+width)
   ab(2*width+1+i-j,j)=a(i,j)
  enddo
 enddo
 x=b
 call dgbsv(n,width,width,nrhs,ab,ldab,ipiv,x,n,info)
 ok=info==0.and.all(ieee_is_finite(x))
 end subroutine

 subroutine screened_electrons_native(g,protons,radius,ne,v,mu,error,status,initial_v)
 type(tf_grid_t),intent(in)::g
 real(wp),intent(in)::protons(:),radius
 real(wp),intent(out)::ne(:),v(:),mu,error
 character(*),intent(out)::status
 real(wp),optional,intent(in)::initial_v(:)
 real(hp)::vh(g%points),neh(g%points),chi(g%points),res(g%points),trial(g%points),tr(g%points)
 real(hp)::mu_h,mu_trial,denom,dmu,before,after,target,avg,weighted(g%points-1),dv(g%points),gauge
 real(wp)::lap(g%points,g%points),jac(g%points-1,g%points-1),rhs(g%points-1,2),sol(g%points-1,2)
 real(wp)::cold(g%points)
 integer::n,k,i,power
 logical::ok,accepted
 n=g%points;status='invalid electron profile/radius';error=huge(1._wp);ne=0;v=0;mu=tf_me
 if(size(protons)/=n.or.size(ne)/=n.or.size(v)/=n)return
 if(.not.all(ieee_is_finite(protons)).or.any(protons<0).or..not.ieee_is_finite(radius).or.radius<=0)return
 avg=sum(real(g%weights,hp)*real(protons,hp))
 if(avg==0)then
  status='ok';error=0;return
 endif
 lap=g%lap/radius**2
 if(present(initial_v))then
  if(size(initial_v)/=n)return
  if(.not.all(ieee_is_finite(initial_v)))return
  vh=real(initial_v,hp)
 else
  call tf_poisson_native(g,protons-real(avg,wp),radius,cold);vh=real(cold,hp)
 endif
 vh=vh-vh(n);call neutral(vh,mu_h)
 target=min(2e-12_hp,real(4*tf_pi*tf_e2,hp)*avg*1e-6_hp)
 do k=1,60
  call terms(vh,mu_h,neh,chi,res)
  if(maxval(abs(res))<target)exit
  jac=lap(:n-1,:n-1)
  do i=1,n-1
   jac(i,i)=jac(i,i)+real(4*tf_pi*tf_e2*chi(i),wp)
  enddo
  rhs(:,1)=-real(res(:n-1),wp);rhs(:,2)=real(4*tf_pi*tf_e2*chi(:n-1),wp)
  call band_solve(jac,rhs,1,sol,ok)
  status='screened electron linear solve failed'
  if(.not.ok)return
  weighted=real(g%weights(:n-1),hp)*chi(:n-1)
  denom=sum(weighted*(1-real(sol(:,2),hp)))+g%weights(n)*chi(n)
  status='screened electron susceptibility outside numerical support'
  if(.not.ieee_is_finite(denom).or.denom<=0)return
  dmu=(-sum(real(g%weights,hp)*(neh-real(protons,hp)))-sum(weighted*real(sol(:,1),hp)))/denom
  dv(:n-1)=real(sol(:,1),hp)-real(sol(:,2),hp)*dmu;dv(n)=0
  before=sqrt(sum(res*res));accepted=.false.
  do power=0,31
   trial=vh+2._hp**(-power)*dv
   call neutral(trial,mu_trial);call terms(trial,mu_trial,neh,chi,tr)
   after=sqrt(sum(tr*tr))
   if(after<before)then
    vh=trial;mu_h=mu_trial;accepted=.true.;exit
   endif
  enddo
  status='screened electron Newton line search failed'
  if(.not.accepted)return
 enddo
 gauge=sum(real(g%weights,hp)*vh);v=real(vh-gauge,wp);mu=real(mu_h+gauge,wp)
 ! Public rounded state is independently re-evaluated in extended arithmetic.
 call terms(real(v,hp),real(mu,hp),neh,chi,res);ne=real(neh,wp);error=real(maxval(abs(res)),wp)
 status='screened electron PDE residual exceeds tolerance'
 if(error>2e-10_wp)return
 status='screened electron neutrality failure'
 if(abs(sum(g%weights*(protons-ne)))>2e-11_wp)return
 status='ok'
 contains
  subroutine density(pot,q,d,dd)
   real(hp),intent(in)::pot(:),q
   real(hp),intent(out)::d(:),dd(:)
   real(hp)::chem(n),mom(n),kin(n)
   chem=q+pot;kin=max(chem-real(tf_me,hp),0._hp)
   mom=sqrt(kin*max(chem+real(tf_me,hp),0._hp))
   d=(mom/real(tf_hbarc,hp))**3/(3*real(tf_pi,hp)**2)
   dd=0
   where(chem>real(tf_me,hp))dd=chem*mom/(real(tf_pi,hp)**2*real(tf_hc3,hp))
  end subroutine
  subroutine neutral(pot,q)
   real(hp),intent(in)::pot(:)
   real(hp),intent(out)::q
   real(hp)::lo,hi,mid,d(n),dd(n),f
   integer::j
   lo=real(tf_me,hp)-maxval(pot)
   hi=sqrt((real(tf_hbarc,hp)*(3*real(tf_pi,hp)**2*avg)**(1._hp/3))**2+real(tf_me,hp)**2)-minval(pot)
   hi=hi+max(1e-10_hp,1e-12_hp*abs(hi))
   do j=1,160
    mid=(lo+hi)/2
    call density(pot,mid,d,dd);f=sum(real(g%weights,hp)*d)-avg
    if(f>0)then
     hi=mid
    else
     lo=mid
    endif
    if(hi-lo<=2*epsilon(1._hp)*max(abs(mid),1._hp))exit
   enddo
   q=(lo+hi)/2
  end subroutine
  subroutine terms(pot,q,d,dd,resid)
   real(hp),intent(in)::pot(:),q
   real(hp),intent(out)::d(:),dd(:),resid(:)
   real(hp)::flux(n-1),action(n)
   call density(pot,q,d,dd)
   flux=real(g%face,hp)*(pot(:n-1)-pot(2:))/real(radius**2,hp)
   action=0;action(:n-1)=flux;action(2:)=action(2:)-flux
   resid=action/real(g%weights,hp)-real(4*tf_pi*tf_e2,hp)*(real(protons,hp)-d)
  end subroutine
 end subroutine

 subroutine meson_fields(p,g,nn,np,radius,s,gs,gw,gr,dgs,dgw,dgr,status)
 type(tf_model_t),intent(in)::p
 type(tf_grid_t),intent(in)::g
 real(wp),intent(in)::nn(:),np(:),radius
 type(tf_state_t),intent(inout)::s
 real(wp),intent(out)::gs(:),gw(:),gr(:),dgs(:),dgw(:),dgr(:)
 character(*),intent(out)::status
 real(wp)::lap(g%points,g%points),phi(g%points),w(g%points),r(g%points),n(g%points)
 real(wp)::res(g%points),jac(g%points,g%points),rhs(g%points,1),sol(g%points,1),trial(g%points),rtest(g%points)
 real(wp)::rw(g%points),rr(g%points),tw(g%points),tr(g%points),rwtest(g%points),rrtest(g%points)
 real(wp)::jvec(2*g%points,2*g%points),bvec(2*g%points,1),dvec(2*g%points,1)
 real(wp)::dw(g%points),dr(g%points),scales(8),before,u,du,ddu,pot,ww,vr,dww,drr,dwr
 integer::i,j,k,ng
 logical::ok,accepted
 ng=g%points;n=nn+np;lap=tf_hbarc**2*g%lap/radius**2
 do i=1,ng
  call tf_couplings(p,n(i),gs(i),gw(i),gr(i),dgs(i),dgw(i),dgr(i))
 enddo
 if(allocated(s%phi))then
  phi=s%phi;w=s%omega;r=s%rho
 else
  phi=gs*n*tf_hc3/p%ms**2;w=gw*n*tf_hc3/p%mw**2;r=gr*(np-nn)*tf_hc3/(2*p%mr**2)
 endif
 phi=max(0._wp,min(phi,.95_wp*min(p%mn,p%mp)/gs))
 scales=[1._wp,.5_wp,.25_wp,.125_wp,.0625_wp,.03125_wp,.015625_wp,.0078125_wp]
 do k=1,45
  call scalar_res(phi,res,jac)
  s%scalar_pde_residual=maxval(abs(res))/tf_hc3
  if(s%scalar_pde_residual<2e-10_wp)exit
  rhs(:,1)=-res;call band_solve(jac,rhs,1,sol,ok)
  status='scalar field linear solve failed';if(.not.ok)return
  before=sqrt(sum(res*res));accepted=.false.
  do j=1,8
   trial=phi+scales(j)*sol(:,1)
   if(any(trial<0).or.any(gs*trial>=.999_wp*min(p%mn,p%mp)))cycle
   call scalar_res(trial,rtest,jac)
   if(sqrt(sum(rtest*rtest))<before)then
    phi=trial;accepted=.true.;exit
   endif
  enddo
  status='scalar field Newton line search failed';if(.not.accepted)return
 enddo
 call scalar_res(phi,res,jac);s%scalar_pde_residual=maxval(abs(res))/tf_hc3
 status='scalar PDE residual exceeds tolerance'
 if(s%scalar_pde_residual>2e-8_wp)return
 do k=1,30
  call vector_res(w,r,rw,rr)
  s%vector_pde_residual=max(maxval(abs(rw)),maxval(abs(rr)))/tf_hc3
  if(s%vector_pde_residual<2e-10_wp)exit
  jvec=0
  do i=1,ng
   call tf_vector_potential(p,w(i),r(i),pot,ww,vr,dww,drr,dwr)
   jvec(2*i-1,2*i-1)=lap(i,i)+dww;jvec(2*i,2*i)=lap(i,i)+drr
   jvec(2*i-1,2*i)=dwr;jvec(2*i,2*i-1)=dwr
   if(i<ng)then
    jvec(2*i-1,2*i+1)=lap(i,i+1);jvec(2*i,2*i+2)=lap(i,i+1)
    jvec(2*i+1,2*i-1)=lap(i+1,i);jvec(2*i+2,2*i)=lap(i+1,i)
   endif
   bvec(2*i-1,1)=-rw(i);bvec(2*i,1)=-rr(i)
  enddo
  call band_solve(jvec,bvec,2,dvec,ok)
  status='vector field linear solve failed';if(.not.ok)return
  dw=dvec(1:2*ng:2,1);dr=dvec(2:2*ng:2,1)
  before=sqrt(sum(rw*rw)+sum(rr*rr));accepted=.false.
  do j=1,6
   tw=w+scales(j)*dw;tr=r+scales(j)*dr;call vector_res(tw,tr,rwtest,rrtest)
   if(sqrt(sum(rwtest*rwtest)+sum(rrtest*rrtest))<before)then
    w=tw;r=tr;accepted=.true.;exit
   endif
  enddo
  status='vector field Newton line search failed';if(.not.accepted)return
 enddo
 call vector_res(w,r,rw,rr)
 s%vector_pde_residual=max(maxval(abs(rw)),maxval(abs(rr)))/tf_hc3
 status='vector PDE residual exceeds tolerance'
 if(s%vector_pde_residual>2e-8_wp)return
 if(.not.allocated(s%phi))allocate(s%phi(ng),s%omega(ng),s%rho(ng))
 s%phi=phi;s%omega=w;s%rho=r;status='ok'
 contains
  subroutine scalar_res(v,res,jac)
   real(wp),intent(in)::v(:)
   real(wp),intent(out)::res(:),jac(:,:)
   real(wp)::en,ep,sn,sp,dn,dp,ef
   integer::i
   res=matmul(lap,v);jac=lap
   do i=1,ng
    call tf_fermion(nn(i),p%mn-gs(i)*v(i),en,sn,dn,ef)
    call tf_fermion(np(i),p%mp-gs(i)*v(i),ep,sp,dp,ef)
    call tf_scalar_potential(p,v(i),u,du,ddu)
    res(i)=res(i)+du-gs(i)*(sn+sp)*tf_hc3
    jac(i,i)=jac(i,i)+ddu+gs(i)**2*(dn+dp)*tf_hc3
   enddo
  end subroutine
  subroutine vector_res(w,r,rw,rr)
   real(wp),intent(in)::w(:),r(:)
   real(wp),intent(out)::rw(:),rr(:)
   integer::i
   rw=matmul(lap,w);rr=matmul(lap,r)
   do i=1,ng
    call tf_vector_potential(p,w(i),r(i),pot,ww,vr,dww,drr,dwr)
    rw(i)=rw(i)+ww-gw(i)*n(i)*tf_hc3
    rr(i)=rr(i)+vr-gr(i)*(np(i)-nn(i))*tf_hc3/2
   enddo
  end subroutine
 end subroutine

 subroutine evaluate_spatial_tf_profile(model_id,nbar,g,nn,np,ne,radius,s,status)
 character(*),intent(in)::model_id
 real(wp),intent(in)::nbar,nn(:),np(:),ne(:),radius
 type(tf_grid_t),intent(in)::g
 type(tf_state_t),intent(inout)::s
 character(*),intent(out)::status
 type(tf_model_t)::p
 real(wp)::gs(g%points),gw(g%points),gr(g%points),dgs(g%points),dgw(g%points),dgr(g%points)
 real(wp)::en,ep,ee,sn,sp,ss,ds,efn,efp,efe,u,du,ddu,vp,dw,dr,dww,drr,dwr,rearr,local(g%points),n,n3
 integer::i,ng
 ng=g%points
 status='invalid TF profile or density domain'
 if(size(nn)/=ng.or.size(np)/=ng.or.size(ne)/=ng)return
 if(.not.all(ieee_is_finite(nn+np+ne)).or.any(nn<0).or.any(np<0).or.any(ne<0))return
 if(.not.ieee_is_finite(radius).or.radius<=0.or..not.ieee_is_finite(nbar).or.nbar<=0.or.nbar>=.3_wp)return
 call load_tf_model(model_id,p,status);if(status/='ok')return
 if(allocated(s%phi))then
  if(size(s%phi)/=ng)then
   s=tf_state_t()
  endif
 endif
 call meson_fields(p,g,nn,np,radius,s,gs,gw,gr,dgs,dgw,dgr,status);if(status/='ok')return
 if(.not.allocated(s%nn))allocate(s%nn(ng),s%np(ng),s%ne(ng),s%coulomb(ng),s%mu_n(ng),s%mu_p(ng),s%mu_e(ng))
 s%nn=nn;s%np=np;s%ne=ne;s%model_id=p%id;s%n_b=nbar;s%points=ng;s%dimension=g%dimension;s%radius_fm=radius
 call tf_poisson_native(g,np-ne,radius,s%coulomb)
 do i=1,ng
  n=nn(i)+np(i);n3=np(i)-nn(i)
  call tf_fermion(nn(i),p%mn-gs(i)*s%phi(i),en,sn,ds,efn)
  call tf_fermion(np(i),p%mp-gs(i)*s%phi(i),ep,sp,ds,efp)
  call tf_fermion(ne(i),tf_me,ee,ss,ds,efe)
  call tf_scalar_potential(p,s%phi(i),u,du,ddu)
  call tf_vector_potential(p,s%omega(i),s%rho(i),vp,dw,dr,dww,drr,dwr)
  local(i)=en+ep+ee+u/tf_hc3+gw(i)*s%omega(i)*n+gr(i)*s%rho(i)*n3/2-vp/tf_hc3
  rearr=dgw(i)*s%omega(i)*n+dgr(i)*s%rho(i)*n3/2-dgs(i)*s%phi(i)*(sn+sp)
  s%mu_n(i)=efn+gw(i)*s%omega(i)-gr(i)*s%rho(i)/2+rearr
  s%mu_p(i)=efp+gw(i)*s%omega(i)+gr(i)*s%rho(i)/2+rearr+s%coulomb(i)
  s%mu_e(i)=efe-s%coulomb(i)
 enddo
 s%gradient_sigma=sum(g%weights*s%phi*matmul(g%lap,s%phi))/(2*tf_hbarc*radius**2)
 s%gradient_omega=sum(g%weights*s%omega*matmul(g%lap,s%omega))/(2*tf_hbarc*radius**2)
 s%gradient_rho=sum(g%weights*s%rho*matmul(g%lap,s%rho))/(2*tf_hbarc*radius**2)
 s%coulomb_energy=.5_wp*sum(g%weights*(np-ne)*s%coulomb)
 s%epsilon=sum(g%weights*local)+s%gradient_sigma-s%gradient_omega-s%gradient_rho+s%coulomb_energy
 s%d_epsilon_d_log_radius=2*(-s%gradient_sigma+s%gradient_omega+s%gradient_rho+s%coulomb_energy)
 s%mean_baryon_density=sum(g%weights*(nn+np));s%mean_charge_density=sum(g%weights*(np-ne))
 s%mu_B=sum(g%weights*s%mu_n);s%mu_Q_lepton=sum(g%weights*s%mu_e)
 s%pressure=nbar*s%mu_B-s%epsilon;s%energy_per_baryon=s%epsilon/nbar
 s%Y_n=sum(g%weights*nn)/nbar;s%Y_p=sum(g%weights*np)/nbar;s%Y_e=sum(g%weights*ne)/nbar;s%Y_mu=0
 s%muon_threshold_crossed=maxval(s%mu_e+s%coulomb)>tf_mmu
 s%baryon_constraint_error=abs(s%mean_baryon_density-nbar);s%charge_constraint_error=abs(s%mean_charge_density)
 s%baryon_density_contrast=(maxval(nn+np)-minval(nn+np))/nbar
 status='ok'
 end subroutine

 subroutine evaluate_screened_tf_profile(model_id,nbar,g,nn,np,radius,s,status)
 character(*),intent(in)::model_id
 real(wp),intent(in)::nbar,nn(:),np(:),radius
 type(tf_grid_t),intent(in)::g
 type(tf_state_t),intent(inout)::s
 character(*),intent(out)::status
 real(wp)::ne(g%points),v(g%points),mu,err
 if(allocated(s%coulomb))then
  call screened_electrons_native(g,np,radius,ne,v,mu,err,status,s%coulomb)
 else
  call screened_electrons_native(g,np,radius,ne,v,mu,err,status)
 endif
 if(status/='ok')return
 call evaluate_spatial_tf_profile(model_id,nbar,g,nn,np,ne,radius,s,status)
 s%mu_Q_lepton=mu;s%electron_pde_residual=err
 end subroutine

 subroutine initial_tf_profile(g,nbar,seed,radius,nn,np,status)
 type(tf_grid_t),intent(in)::g
 real(wp),intent(in)::nbar,radius
 character(*),intent(in)::seed
 real(wp),intent(out)::nn(:),np(:)
 character(*),intent(out)::status
 real(wp)::f,x0,shape(g%points),avg
 status='ok'
 if(seed=='uniform')then
  nn=nbar*(1-.035_wp);np=nbar*.035_wp
 else
  f=min(.65_wp,max((2._wp/g%points)**g%dimension,.45_wp*nbar/.13_wp));x0=f**(1._wp/g%dimension)
  shape=1/(1+exp(max(-60._wp,min(60._wp,(g%x-x0)*radius/.7_wp))))
  if(seed=='hole')then
   shape=1-shape
  elseif(seed/='droplet')then
   status='seed must be uniform, droplet or hole';return
  endif
  avg=max(sum(g%weights*shape),1e-12_wp)
  np=.10_wp*nbar*shape/avg
  nn=.55_wp*nbar+(.45_wp*nbar-sum(g%weights*np))*shape/avg
  nn=nn*(nbar-sum(g%weights*np))/sum(g%weights*nn)
 endif
 end subroutine

 subroutine solve_spatial_tf(model_id,nbar,npoints,dimension,seed,radius,maxiter,ftol,s,status,initial_nn,initial_np,radius_bounds)
 character(*),intent(in)::model_id,seed
 real(wp),intent(in)::nbar,radius,ftol
 integer,intent(in)::npoints,dimension,maxiter
 type(tf_state_t),intent(out)::s
 character(*),intent(out)::status
 real(wp),optional,intent(in)::initial_nn(:),initial_np(:),radius_bounds(2)
 type(tf_grid_t)::grid
 type(tf_model_t)::model
 type(tf_state_t)::work
 type(native_optimizer_result_t)::op
 real(wp)::x(2*npoints+1),lower(2*npoints+1),upper(2*npoints+1),aeq(1,2*npoints+1),beq(1)
 real(wp)::nn(npoints),np(npoints),rb(2),rad,threshold,denom,mub,kkt,eq,fvalue,grad(2*npoints+1)
 integer::i,j,ng
 logical::ok,active
 character(160)::evaluation_status
 s=tf_state_t();ng=npoints;status='invalid spatial TF solve input'
 if(npoints<8.or.dimension<1.or.dimension>3.or.nbar<=0.or.nbar>=.3_wp.or.radius<=0.or.maxiter<1.or.ftol<=0)return
 if(.not.all(ieee_is_finite([nbar,radius,ftol])))return
 rb=[5._wp,180._wp];if(present(radius_bounds))rb=radius_bounds
 if(any(rb<=0).or.rb(2)<=rb(1).or.radius<rb(1).or.radius>rb(2))return
 call load_tf_model(model_id,model,status);if(status/='ok')return
 call make_tf_grid(npoints,dimension,grid,status);if(status/='ok')return
 if(present(initial_nn).and.present(initial_np))then
  if(size(initial_nn)/=ng.or.size(initial_np)/=ng)return
  nn=initial_nn;np=initial_np
 else
  call initial_tf_profile(grid,nbar,seed,radius,nn,np,status);if(status/='ok')return
 endif
 if(any(nn<0).or.any(np<0).or..not.all(ieee_is_finite(nn+np)))then
  status='invalid nonnegative initial densities';return
 endif
 x(:ng)=nn*grid%sqrtw/nbar;x(ng+1:2*ng)=np*grid%sqrtw/nbar;x(2*ng+1)=log(radius)
 lower=0;upper(:ng)=.4_wp*grid%sqrtw/nbar;upper(ng+1:2*ng)=upper(:ng)
 lower(2*ng+1)=log(rb(1));upper(2*ng+1)=log(rb(2))
 aeq(1,:ng)=grid%sqrtw;aeq(1,ng+1:2*ng)=grid%sqrtw;aeq(1,2*ng+1)=0;beq=1
 evaluation_status='ok'
 call minimize_slsqp_native(x,lower,upper,aeq,beq,objective,ftol,maxiter,op)
 nn=max(0._wp,x(:ng))*nbar/grid%sqrtw;np=max(0._wp,x(ng+1:2*ng))*nbar/grid%sqrtw;rad=exp(x(2*ng+1))
 call evaluate_screened_tf_profile(model_id,nbar,grid,nn,np,rad,work,status)
 if(status/='ok')return
 s=work;s%seed=seed;s%optimizer_success=op%success;s%optimizer_mode=op%mode;s%iterations=op%iterations;s%evaluations=op%evaluations
 threshold=max(1e-10_wp,nbar*1e-7_wp);denom=0;mub=0
 do i=1,ng
  if(nn(i)>threshold)then
   denom=denom+grid%weights(i);mub=mub+grid%weights(i)*s%mu_n(i)
  endif
  if(np(i)>threshold)then
   denom=denom+grid%weights(i);mub=mub+grid%weights(i)*(s%mu_p(i)+s%mu_Q_lepton)
  endif
 enddo
 if(denom<=0)then
  status='no active nuclear density support';return
 endif
 mub=mub/denom;kkt=maxval(abs(s%mu_e-s%mu_Q_lepton))
 do i=1,ng
  if(nn(i)>threshold)then
   kkt=max(kkt,abs(s%mu_n(i)-mub))
  else
   kkt=max(kkt,mub-s%mu_n(i))
  endif
  if(np(i)>threshold)then
   kkt=max(kkt,abs(s%mu_p(i)+s%mu_Q_lepton-mub))
  else
   kkt=max(kkt,mub-s%mu_p(i)-s%mu_Q_lepton)
  endif
 enddo
 s%mu_B=mub;s%pressure=nbar*mub-s%epsilon;s%max_chemical_KKT_residual_MeV=kkt
 s%radius_stationarity_MeV=abs(s%d_epsilon_d_log_radius/nbar)
 s%radius_at_search_bound=min(abs(log(rad/rb(1))),abs(log(rad/rb(2))))<1e-5_wp
 s%density_at_search_bound=max(maxval(nn),maxval(np))>.4_wp*(1-1e-7_wp)
 s%numerical_status='REQUIRES_REFINEMENT'
 if(op%success.and.kkt<.02_wp.and.s%radius_stationarity_MeV<2e-4_wp.and..not.s%radius_at_search_bound.and..not.s%density_at_search_bound) &
  s%numerical_status='STATIONARY'
 if(s%baryon_density_contrast<1e-4_wp.and.kkt<.02_wp.and..not.s%density_at_search_bound)s%numerical_status='UNIFORM_STATIONARY'
 s%geometry=tf_geometry_label(s);status='ok'
 contains
  subroutine objective(z,f,gradient,ok)
   real(wp),intent(in)::z(:)
   real(wp),intent(out)::f,gradient(:)
   logical,intent(out)::ok
   real(wp)::lnn(ng),lnp(ng),lr
   lnn=z(:ng)*nbar/grid%sqrtw;lnp=z(ng+1:2*ng)*nbar/grid%sqrtw;lr=exp(z(2*ng+1))
   ! SLSQP may evaluate a tiny negative roundoff outside a bound. Reject
   ! actual unsupported trials; a numerical clamp here matches explicit bounds.
   lnn=max(lnn,0._wp);lnp=max(lnp,0._wp)
   call evaluate_screened_tf_profile(model_id,nbar,grid,lnn,lnp,lr,work,evaluation_status)
   ok=evaluation_status=='ok';f=0;gradient=0
   if(.not.ok)return
   f=(work%epsilon-model%mn*work%mean_baryon_density)/nbar
   gradient(:ng)=grid%sqrtw*(work%mu_n-model%mn)
   gradient(ng+1:2*ng)=grid%sqrtw*(work%mu_p+work%mu_Q_lepton-model%mn)
   gradient(2*ng+1)=work%d_epsilon_d_log_radius/nbar
  end subroutine
 end subroutine

 subroutine solve_tf_family_native(model_id,nbar,npoints,dimension,seed,selected,status,attempts)
 ! Finite, deterministic initial-radius search. All failed attempts remain in
 ! attempts; a family is resolved only by an independently accepted native state.
 character(*),intent(in)::model_id,seed
 real(wp),intent(in)::nbar
 integer,intent(in)::npoints,dimension
 type(tf_state_t),intent(out)::selected
 character(*),intent(out)::status
 type(tf_state_t),allocatable,intent(out)::attempts(:)
 type(tf_state_t)::trial,refined,uniform_seed
 type(tf_grid_t)::grid
 type(tf_state_t),allocatable::saved(:)
 real(wp)::radii(3),initial_radius,shape(npoints),nn(npoints),np(npoints),amplitude
 integer::j,k,chosen,unresolved,lower,level
 integer::dims(1)
 character(16)::seeds(1)
 character(160)::trial_status,selection_status
 allocate(saved(19));k=0;selected=tf_state_t();dims=[dimension];seeds=[seed]
 radii=[30._wp,45._wp,70._wp]
 if(nbar<.003_wp)radii=[45._wp,70._wp,100._wp]
 do j=1,3
  initial_radius=radii(j)
  call solve_spatial_tf(model_id,nbar,npoints,dimension,seed,initial_radius,1200,1e-11_wp,trial,trial_status)
  call save(trial,trial_status,initial_radius)
  call select_tf_candidates_native(saved(:k),dims,seeds,nbar,chosen,selection_status,unresolved,lower)
  if(chosen>0)exit
  if(trial_status=='ok'.and.allocated(trial%nn).and.allocated(trial%np))then
   call solve_spatial_tf(model_id,nbar,npoints,dimension,seed,trial%radius_fm,1800,1e-13_wp,refined,trial_status, &
     initial_nn=trial%nn,initial_np=trial%np)
   call save(refined,trial_status,trial%radius_fm)
   call select_tf_candidates_native(saved(:k),dims,seeds,nbar,chosen,selection_status,unresolved,lower)
   if(chosen>0)exit
  endif
 enddo
 if(chosen==0.and.seed/='uniform')then
  ! A modest, explicit cluster perturbation is a further initial guess, never
  ! a replacement of failed clustered physics by an imposed uniform answer.
  call solve_spatial_tf(model_id,nbar,npoints,dimension,'uniform',30._wp,1200,1e-13_wp,uniform_seed,trial_status)
  call save(uniform_seed,trial_status,30._wp,'uniform')
  if(trial_status=='ok'.and.uniform_seed%numerical_status=='UNIFORM_STATIONARY')then
   call make_tf_grid(npoints,dimension,grid,trial_status)
   shape=exp(-(grid%x/.45_wp)**2);shape=shape-sum(grid%weights*shape)
   radii=[15._wp,30._wp,60._wp]
   weak:do level=1,2
    amplitude=.25_wp;if(level==2)amplitude=.08_wp
    if(seed=='hole')amplitude=-amplitude
    nn=nbar*uniform_seed%Y_n*(1+.3_wp*amplitude*shape)
    np=nbar*uniform_seed%Y_p*(1+amplitude*shape)
    do j=1,3
     call solve_spatial_tf(model_id,nbar,npoints,dimension,seed,radii(j),1200,1e-13_wp,trial,trial_status, &
       initial_nn=nn,initial_np=np)
     call save(trial,trial_status,radii(j))
     call select_tf_candidates_native(saved(:k),dims,seeds,nbar,chosen,selection_status,unresolved,lower)
     if(chosen>0)exit weak
     if(trial_status=='ok'.and.allocated(trial%nn).and.allocated(trial%np))then
      call solve_spatial_tf(model_id,nbar,npoints,dimension,seed,trial%radius_fm,1800,1e-14_wp,refined,trial_status, &
        initial_nn=trial%nn,initial_np=trial%np)
      call save(refined,trial_status,trial%radius_fm)
      call select_tf_candidates_native(saved(:k),dims,seeds,nbar,chosen,selection_status,unresolved,lower)
      if(chosen>0)exit weak
     endif
    enddo
   enddo weak
  endif
 endif
 allocate(attempts(k));attempts=saved(:k)
 call select_tf_candidates_native(attempts,dims,seeds,nbar,chosen,status,unresolved,lower)
 if(chosen>0)then
  selected=attempts(chosen);status='ok'
 endif
 contains
  subroutine save(state,ss,radius,seed_override)
   type(tf_state_t),intent(in)::state
   character(*),intent(in)::ss
   real(wp),intent(in)::radius
   character(*),optional,intent(in)::seed_override
   k=k+1;saved(k)=state
   saved(k)%model_id=model_id;saved(k)%n_b=nbar;saved(k)%points=npoints;saved(k)%dimension=dimension;saved(k)%seed=seed
   if(present(seed_override))saved(k)%seed=seed_override
   saved(k)%evaluation_status=ss
   if(ss/='ok')then
    saved(k)%numerical_status='EVALUATION_FAILED'
    if(saved(k)%radius_fm<=0)saved(k)%radius_fm=radius
   endif
  end subroutine
 end subroutine

 function tf_geometry_label(s)result(label)
 type(tf_state_t),intent(in)::s
 character(40)::label
 real(wp)::b(s%points),q(s%points),scale,contrast,avg,w(s%points)
 logical::uniform,bi,qi,bm,qm
 integer::i,d,n
 n=s%points;d=s%dimension;label='unavailable'
 if(.not.allocated(s%nn).or.n<2)return
 b=s%nn+s%np;q=s%np;scale=s%mean_baryon_density
 do i=1,n
  w(i)=(real(i,wp)/n)**d-(real(i-1,wp)/n)**d
 enddo
 if((maxval(b)-minval(b))/scale<1e-4_wp)then
  uniform=species_uniform(s%nn).and.species_uniform(s%np).and.species_uniform(s%ne)
  label='near_uniform_profile';if(uniform)label='uniform';return
 endif
 call occupancy(b,bi,bm);call occupancy(q,qi,qm)
 if(d==1)then
  label='slab'
 elseif(bi.and.qi)then
  label='noncanonical_spherical_profile';if(d==2)label='noncanonical_cylindrical_profile'
 elseif((bm.and.qm).or.(bi.neqv.qi).or.(bm.neqv.qm))then
  label='noncanonical_radial_profile'
 elseif(d==2)then
  label='tube';if(b(1)>b(n))label='rod'
 else
  label='bubble';if(b(1)>b(n))label='droplet'
 endif
 contains
 logical function species_uniform(a)
 real(wp),intent(in)::a(:)
 real(wp)::delta,mean
 delta=maxval(a)-minval(a);mean=sum(w*a)
 species_uniform=delta<=max(1e-18_wp,1e-12_wp*scale).or.delta<=1e-4_wp*abs(mean)
 end function
 subroutine occupancy(a,interior,multiple)
 real(wp),intent(in)::a(:)
 logical,intent(out)::interior,multiple
 logical::mask(n),last
 integer::j,count
 if(maxval(a)-minval(a)<=0)then
  interior=.false.;multiple=.false.;return
 endif
 mask=(a-minval(a))/(maxval(a)-minval(a))>.5_wp;last=.false.;count=0
 do j=1,n
  if(mask(j).and..not.last)count=count+1
  last=mask(j)
 enddo
 interior=count==1.and..not.mask(1).and..not.mask(n);multiple=count>1
 end subroutine
 end function

 subroutine select_tf_candidates_native(states,expected_dims,expected_seeds,nbar,index,status,unresolved_count,lower_count)
 type(tf_state_t),intent(in)::states(:)
 integer,intent(in)::expected_dims(:)
 character(*),intent(in)::expected_seeds(:)
 real(wp),intent(in)::nbar
 integer,intent(out)::index,unresolved_count,lower_count
 character(*),intent(out)::status
 integer::i,j,best
 logical::valid(size(states)),feasible(size(states)),resolved
 real(wp)::tol,beste
 index=0;best=0;unresolved_count=0;lower_count=0;tol=max(1e-12_wp,1e-8_wp*nbar);beste=huge(1._wp)
 do i=1,size(states)
  feasible(i)=abs(states(i)%n_b-nbar)<=tol.and.states(i)%baryon_constraint_error<=tol.and.states(i)%charge_constraint_error<=tol.and. &
    states(i)%scalar_pde_residual<=1e-8_wp.and.states(i)%vector_pde_residual<=1e-8_wp.and.states(i)%electron_pde_residual<=2e-10_wp
  valid(i)=feasible(i).and.(states(i)%numerical_status=='STATIONARY'.or.states(i)%numerical_status=='UNIFORM_STATIONARY').and. &
    .not.states(i)%muon_threshold_crossed.and.ieee_is_finite(states(i)%energy_per_baryon).and.ieee_is_finite(states(i)%epsilon)
  if(valid(i).and.states(i)%energy_per_baryon<beste)then
   best=i;beste=states(i)%energy_per_baryon
  endif
 enddo
 do j=1,size(expected_dims)
  resolved=.false.
  do i=1,size(states)
   if(valid(i).and.states(i)%dimension==expected_dims(j).and.states(i)%seed==expected_seeds(j))resolved=.true.
  enddo
  if(.not.resolved)unresolved_count=unresolved_count+1
 enddo
 if(best>0)then
  do i=1,size(states)
   if(.not.valid(i).and.feasible(i).and.ieee_is_finite(states(i)%energy_per_baryon))then
    if(states(i)%energy_per_baryon<beste-1e-6_wp)lower_count=lower_count+1
   endif
  enddo
 endif
 status='NO_STATIONARY_CANDIDATE'
 if(best==0)return
 status='HOLD_UNRESOLVED_GEOMETRY';if(unresolved_count>0)return
 status='LOWER_UNCONVERGED_ATTEMPT';if(lower_count>0)return
 index=best;status='LOWEST_EVALUATED_STATIONARY_CANDIDATE'
 end subroutine
end module
