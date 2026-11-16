! Spherical SPG candidate CLD physics in Fortran.
! The sigma0,b prescription is explicit INPUT, never an inferred author calibration.
! Scurto Eqs21-25; uniform electrons, no muons, sharp interface, baryon/charge
! constraints eliminated exactly; radius eliminated by E_surface=2 E_Coulomb.
module compose_eos_spg_native
 use compose_eos_types, only: wp
 use compose_eos_gdfm_native
 use compose_eos_optimize_native, only: minimize_slsqp_native,native_optimizer_result_t
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 real(wp),parameter::pi=acos(-1._wp),e2=1.43996448_wp,me=.51099895_wp,mmu=105.6583755_wp
 real(wp),parameter::scales(4)=[1._wp,.1_wp,.1_wp,1._wp]
 type,public::spg_native_result_t
  real(wp)::nb=0,epsilon=0,pressure=0,mu_b=0,mu_e=0
  real(wp)::yn=0,yp=0,ye=0,ymu=0,yh=0,aav=0,zav=0
  real(wp)::z(4)=0,gradient(4)=0,max_projected_gradient=huge(1._wp),min_scaled_curvature=0
  real(wp)::gas_proton_margin=0,muon_margin=0,radius=0,volume_fraction=0,n_gas=0
  real(wp)::energy_per_baryon=huge(1._wp),epsilon_surface=0,epsilon_coulomb=0
  integer::iterations=0,attempt_count=0,optimizer_mode=0,optimizer_evaluations=0
  integer::uniform_limit_attempts=0,unresolved_lower_attempts=0
  logical::stationary=.false.,evaluated=.false.,optimizer_success=.false.
  character(160)::status='not evaluated'
 end type
 public::spg_native_cell,spg_native_solve,spg_native_drip,spg_native_endpoint,spg_native_crossing,spg_native_transition
contains
 pure subroutine surface(y,sigma0,b,sigma,derivative)
  real(wp),intent(in)::y,sigma0,b
  real(wp),intent(out)::sigma,derivative
  real(wp)::q,k,scale,d
  q=y*(1-y);k=b+16;scale=max(k,1._wp)
  d=(1-2*y)**2*(1+q+4*q*q)/scale+(k/scale)*q**3
  sigma=sigma0*(k/scale)*q**3/d
  derivative=sigma0*(k/scale)*(3*q*q*(1-2*y)*(1-2*q)/scale)/d/d
 end subroutine

 subroutine spg_native_cell(model,sigma0,b,nb,z,energy,grad,result,status)
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::sigma0,b,nb,z(4)
  real(wp),intent(out)::energy,grad(4)
  type(spg_native_result_t),intent(out)::result
  character(*),intent(out)::status
  real(wp)::t,ni,yi,v,f,ng,yg,contrast,x,phi,dphi,sigma,ds,r
  real(wp)::ei,eg,muni,mupi,mung,mupg,press,ne,ee,pe,mue,coul,surf
  real(wp)::mi,mg,g(5),eps,partialnb,aa
  result=spg_native_result_t();energy=huge(1._wp);grad=0
  status='invalid SPG coordinates or surface parameters'
  if(.not.all(ieee_is_finite([sigma0,b,nb,z])))return
  t=z(1);ni=z(2);yi=z(3);v=z(4)
  if(.not.(nb>0.and.ni>nb.and.t>0.and.t<=1.and.yi>=1e-8_wp.and.yi<.5_wp.and.v>=0.and.v<1))return
  if(sigma0<=0.or.b<=-16)return
  f=t*nb/ni;ng=nb*(1-t)/(1-f);yg=v*yi;contrast=ni*yi-ng*yg
  if(contrast<=0)return
  x=f**(1._wp/3)
  phi=((1-f)/(1+x+x*x))**2*(x+2)/5
  dphi=(1-f**(-2._wp/3))/5
  call surface(yi,sigma0,b,sigma,ds)
  r=(3*sigma/(4*pi*e2*phi*contrast**2))**(1._wp/3)
  call gdfm_native_baryons(model,ni,yi,ei,muni,mupi,press,status)
  if(status/='ok')return
  call gdfm_native_baryons(model,ng,yg,eg,mung,mupg,press,status)
  if(status/='ok')return
  ne=f*ni*yi+(1-f)*ng*yg
  call native_fermion(ne,me,ee,pe,mue)
  coul=2*pi*e2*f*phi*r*r*contrast**2;surf=3*f*sigma/r
  mi=muni+(mupi-muni)*yi;mg=mung+(mupg-mung)*yg
  g(1)=ei-eg+mg*(ng-ni)+mue*ni*(yi-yg) &
   +2*pi*e2*r*r*(contrast**2*(phi+f*dphi)+2*f*phi*contrast*yg*(ni-ng)/(1-f))+3*sigma/r
  g(2)=(2*coul-surf)/r
  g(3)=f*(mi-mg+mue*(yi-yg))+2*coul/contrast*(yi+f*yg/(1-f))
  g(4)=f*ni*(mupi-muni+mue)+2*coul*ni/contrast+3*f*ds/r
  g(5)=(1-f)*ng*(mupg-mung+mue)-2*coul*ng/contrast
  eps=f*ei+(1-f)*eg+ee+coul+surf
  partialnb=mg+mue*yg-2*coul*yg/(contrast*(1-f))
  grad=[g(1)/ni,(g(3)-g(1)*f/ni)/nb,(g(4)+v*g(5))/nb,yi*g(5)/nb]
  energy=(eps-gdfm_mass*nb)/nb
  result%nb=nb;result%epsilon=eps;result%energy_per_baryon=energy
  result%z=z;result%gradient=grad;result%mu_b=partialnb+t*grad(1)
  result%pressure=nb*result%mu_b-eps;result%mu_e=mue
  result%volume_fraction=f;result%radius=r;result%n_gas=ng
  result%gas_proton_margin=mupg-mung+mue-2*coul/(contrast*(1-f))
  result%muon_margin=mmu-mue;result%epsilon_surface=surf;result%epsilon_coulomb=coul
  aa=4*pi*r**3*ni/3
  result%yn=(1-t)*(1-yg);result%yp=(1-t)*yg;result%ye=t*yi+(1-t)*yg
  result%yh=t;result%aav=aa;result%zav=aa*yi;result%ymu=0
  result%evaluated=.true.;result%status='evaluated';status='ok'
 end subroutine

 subroutine solve_linear(a,b,x,n,ok)
  integer,intent(in)::n
  real(wp),intent(in)::a(n,n),b(n)
  real(wp),intent(out)::x(n)
  logical,intent(out)::ok
  real(wp)::m(n,n+1),row(n+1),fac
  integer::i,j,k
  m(:,:n)=a;m(:,n+1)=b;ok=.false.;x=0
  do i=1,n
   k=maxloc(abs(m(i:n,i)),dim=1)+i-1
   if(abs(m(k,i))<1e-18_wp)return
   if(k/=i)then
    row=m(i,:);m(i,:)=m(k,:);m(k,:)=row
   endif
   do j=i+1,n
    fac=m(j,i)/m(i,i);m(j,i:n+1)=m(j,i:n+1)-fac*m(i,i:n+1)
   enddo
  enddo
  do i=n,1,-1
   x(i)=(m(i,n+1)-dot_product(m(i,i+1:n),x(i+1:n)))/m(i,i)
  enddo
  ok=all(ieee_is_finite(x))
 end subroutine

 function smallest_eigen(a,n) result(value)
  integer,intent(in)::n
  real(wp),intent(in)::a(n,n)
  real(wp)::value,m(n,n),off,t,c,s,tau,ap,aq
  integer::i,j,k,p,q,it
  m=(a+transpose(a))/2
  do it=1,100
   off=0;p=1;q=1
   do i=1,n
    do j=i+1,n
     if(abs(m(i,j))>off)then
      off=abs(m(i,j));p=i;q=j
     endif
    enddo
   enddo
   if(off<1e-12_wp)exit
   tau=(m(q,q)-m(p,p))/(2*m(p,q))
   t=sign(1._wp,tau)/(abs(tau)+sqrt(1+tau*tau));c=1/sqrt(1+t*t);s=t*c
   ap=m(p,p);aq=m(q,q)
   m(p,p)=ap-t*m(p,q);m(q,q)=aq+t*m(p,q);m(p,q)=0;m(q,p)=0
   do k=1,n
    if(k==p.or.k==q)cycle
    ap=m(k,p);aq=m(k,q);m(k,p)=c*ap-s*aq;m(p,k)=m(k,p)
    m(k,q)=s*ap+c*aq;m(q,k)=m(k,q)
   enddo
  enddo
  value=minval([(m(i,i),i=1,n)])
 end function

 subroutine polish_candidate(model,sigma0,b,nb,z,ndim,result)
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::sigma0,b,nb
  real(wp),intent(inout)::z(4)
  integer,intent(in)::ndim
  type(spg_native_result_t),intent(out)::result
  real(wp)::lo(4),hi(4),g(4),ga(4),gb(4),e,h(4,4),rhs(4),delta(4)
  real(wp)::za(4),zb(4),trial(4),step,alpha,norm0,norm1,pg(4),hf(4,4),q(4)
  integer::free(4),nf,i,j,k,it,ls
  logical::ok
  character(160)::status
  type(spg_native_result_t)::tmp
  lo=[1e-8_wp,max(.025_wp,nb*(1+1e-6_wp)),.0001_wp,0._wp]
  hi=[1._wp,.25_wp,.499999_wp,.999999_wp]
  do it=1,80
   call spg_native_cell(model,sigma0,b,nb,z,e,g,tmp,status)
   if(status/='ok')exit
   if(z(1)>1-1e-7_wp.and.g(1)<0)z(1)=1
   if(ndim==4.and.z(4)<1e-7_wp.and.g(4)>0)z(4)=0
   call spg_native_cell(model,sigma0,b,nb,z,e,g,tmp,status)
   nf=0
   do i=1,ndim
    if(i==1.and.z(1)>=1-1e-12_wp.and.g(1)<0)cycle
    if(i==4.and.z(4)<=1e-12_wp.and.g(4)>0)cycle
    nf=nf+1;free(nf)=i
   enddo
   if(nf==0)exit
   norm0=maxval(abs(g(free(:nf))*scales(free(:nf))))
   if(maxval(abs(g(free(:nf))))<2e-8_wp)exit
   h=0
   do j=1,nf
    k=free(j);step=2e-6_wp*scales(k)
    za=z;zb=z;za(k)=min(hi(k),z(k)+step);zb(k)=max(lo(k),z(k)-step)
    call spg_native_cell(model,sigma0,b,nb,za,e,ga,tmp,status)
    call spg_native_cell(model,sigma0,b,nb,zb,e,gb,tmp,status)
    h(:nf,j)=(ga(free(:nf))-gb(free(:nf)))*scales(free(:nf))*scales(k)/(za(k)-zb(k))
   enddo
   rhs(:nf)=-g(free(:nf))*scales(free(:nf))
   call solve_linear(h(:nf,:nf),rhs(:nf),delta(:nf),nf,ok)
   if(.not.ok)exit
   alpha=1
   do j=1,nf
    k=free(j)
    if(delta(j)>0)alpha=min(alpha,.99_wp*(hi(k)-z(k))/(delta(j)*scales(k)))
    if(delta(j)<0)alpha=min(alpha,.99_wp*(lo(k)-z(k))/(delta(j)*scales(k)))
   enddo
   if(alpha<1e-14_wp)exit
   do ls=1,25
    trial=z
    trial(free(:nf))=z(free(:nf))+alpha*delta(:nf)*scales(free(:nf))
    call spg_native_cell(model,sigma0,b,nb,trial,e,ga,tmp,status)
    norm1=maxval(abs(ga(free(:nf))*scales(free(:nf))))
    if(status=='ok'.and.norm1<norm0)exit
    alpha=alpha/2
   enddo
   if(ls>25)exit
   z=trial
  enddo
  call spg_native_cell(model,sigma0,b,nb,z,e,g,result,status)
  if(status/='ok')then
   result%status=status;return
  endif
  pg=g
  if(ndim==3)pg(4)=0
  if(z(1)>=1-1e-8_wp)pg(1)=max(g(1),0._wp)
  if(ndim==4.and.z(4)<=1e-8_wp)pg(4)=min(g(4),0._wp)
  result%max_projected_gradient=maxval(abs(pg))
  result%stationary=.false.;result%status='STATIONARITY_RESIDUAL'
  if(result%max_projected_gradient>1e-5_wp)return
  result%status='SEARCH_BOUND_ACTIVE'
  do i=1,3
   if(z(i)-lo(i)<1e-6_wp*(hi(i)-lo(i)))return
   if(i>1.and.hi(i)-z(i)<1e-6_wp*(hi(i)-lo(i)))return
  enddo
  if(ndim==4.and.hi(4)-z(4)<1e-6_wp)return
  result%status='MUON_INSERTION_UNSTABLE'
  if(result%muon_margin<0)return
  result%status='GAS_PROTON_INSERTION_UNSTABLE'
  if(ndim==3.and.result%gas_proton_margin<0)return
  ! Curvature is checked in feasible tangent coordinates, not inferred from SLSQP.
  nf=0
  do i=1,ndim
   if(i==1.and.z(1)>=1-1e-8_wp.and.g(1)<-1e-7_wp)cycle
   if(i==4.and.z(4)<=1e-8_wp.and.g(4)>1e-7_wp)cycle
   nf=nf+1;free(nf)=i
  enddo
  hf=0
  do j=1,nf
   k=free(j);step=2e-6_wp*scales(k)
   za=z;zb=z;za(k)=min(hi(k),z(k)+step);zb(k)=max(lo(k),z(k)-step)
   call spg_native_cell(model,sigma0,b,nb,za,e,ga,tmp,status)
   call spg_native_cell(model,sigma0,b,nb,zb,e,gb,tmp,status)
   hf(:nf,j)=(ga(free(:nf))-gb(free(:nf)))*scales(free(:nf))*scales(k)/(za(k)-zb(k))
  enddo
  if(nf>0)result%min_scaled_curvature=smallest_eigen(hf(:nf,:nf),nf)
  result%status='NEGATIVE_LOCAL_CURVATURE'
  if(result%min_scaled_curvature< -1e-6_wp)return
  result%stationary=.true.;result%status='STATIONARY_CANDIDATE'
 end subroutine

 subroutine optimize_attempt(model,sigma0,b,nb,seed,ndim,result)
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::sigma0,b,nb,seed(4)
  integer,intent(in)::ndim
  type(spg_native_result_t),intent(out)::result
  real(wp)::x(ndim),lo(ndim),hi(ndim),z(4),aeq(0,ndim),beq(0)
  type(native_optimizer_result_t)::opt
  lo(1)=1e-8_wp;lo(2)=max(.025_wp,nb*(1+1e-6_wp));lo(3)=.0001_wp
  hi(1)=1;hi(2)=.25_wp;hi(3)=.499999_wp
  if(ndim==4)then
   lo(4)=0;hi(4)=.999999_wp
  endif
  lo=lo/scales(:ndim);hi=hi/scales(:ndim)
  x=max(lo+1e-10_wp,min(hi-1e-10_wp,seed(:ndim)/scales(:ndim)))
  call minimize_slsqp_native(x,lo,hi,aeq,beq,objective,1e-14_wp,1000,opt)
  z=0;z(:ndim)=x*scales(:ndim)
  call polish_candidate(model,sigma0,b,nb,z,ndim,result)
  result%iterations=opt%iterations;result%optimizer_mode=opt%mode
  result%optimizer_evaluations=opt%evaluations;result%optimizer_success=opt%success
 contains
  subroutine objective(x,f,g,ok)
   real(wp),intent(in)::x(:)
   real(wp),intent(out)::f,g(:)
   logical,intent(out)::ok
   real(wp)::zz(4),gg(4)
   type(spg_native_result_t)::rr
   character(160)::ss
   zz=0;zz(:ndim)=x*scales(:ndim)
   call spg_native_cell(model,sigma0,b,nb,zz,f,gg,rr,ss)
   g=gg(:ndim)*scales(:ndim);ok=ss=='ok'
  end subroutine
 end subroutine

 subroutine spg_native_solve(model,sigma0,b,nb,result,status,start,attempts)
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::sigma0,b,nb
  type(spg_native_result_t),intent(out)::result
  character(*),intent(out)::status
  real(wp),intent(in),optional::start(4)
  type(spg_native_result_t),allocatable,intent(out),optional::attempts(:)
  type(spg_native_result_t)::rows(7),candidate
  real(wp)::seed(4),fixed(3,3)
  integer::i,k,chosen,restricted
  result=spg_native_result_t();status='unsupported SPG density'
  if(.not.ieee_is_finite(nb).or.nb<1e-7_wp.or.nb>.12_wp)return
  fixed(:,1)=[.9_wp,.16_wp,.3_wp];fixed(:,2)=[.5_wp,.15_wp,.14_wp];fixed(:,3)=[.2_wp,.14_wp,.08_wp]
  k=0
  if(present(start))then
   k=k+1;seed=start;seed(4)=0
   call optimize_attempt(model,sigma0,b,nb,seed,3,rows(k))
  endif
  do i=1,3
   k=k+1;seed=0;seed(:3)=fixed(:,i)
   call optimize_attempt(model,sigma0,b,nb,seed,3,rows(k))
  enddo
  chosen=minloc(rows(:k)%energy_per_baryon,dim=1);restricted=chosen
  if(rows(chosen)%gas_proton_margin<0)then
   do i=1,3
    seed=rows(restricted)%z
    select case(i)
    case(1);seed(4)=.001_wp
    case(2);seed(4)=.05_wp
    case(3);seed(4)=.3_wp
    end select
    k=k+1
    call optimize_attempt(model,sigma0,b,nb,seed,4,rows(k))
   enddo
   chosen=minloc(rows(:k)%energy_per_baryon,dim=1)
  endif
  result=rows(chosen);result%attempt_count=k
  status=result%status
  if(result%stationary)status='ok'
  if(present(attempts))then
   allocate(attempts(k));attempts=rows(:k)
  endif
 end subroutine

 subroutine vacuum_cell(model,sigma0,b,nb,result,status)
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::sigma0,b,nb
  type(spg_native_result_t),intent(out)::result
  character(*),intent(out)::status
  real(wp)::x(2),lo(2),hi(2),aeq(0,2),beq(0),z(4),g(4),ga(4),gb(4),e
  real(wp)::h(2,2),delta(2),za(4),zb(4),trial(4),alpha,norm0
  logical::ok
  integer::j,it,ls
  type(native_optimizer_result_t)::opt
  type(spg_native_result_t)::tmp
  x=[1.6_wp,3.3_wp];lo=[.5_wp,.1_wp];hi=[2.5_wp,4.99_wp]
  call minimize_slsqp_native(x,lo,hi,aeq,beq,objective,1e-14_wp,500,opt)
  z=[1._wp,.1_wp*x(1),.1_wp*x(2),0._wp]
  do it=1,40
   call spg_native_cell(model,sigma0,b,nb,z,e,g,result,status)
   if(status/='ok')return
   norm0=maxval(abs(g(2:3)))
   if(norm0<2e-8_wp)exit
   do j=1,2
    za=z;zb=z;za(j+1)=za(j+1)+2e-7_wp;zb(j+1)=zb(j+1)-2e-7_wp
    call spg_native_cell(model,sigma0,b,nb,za,e,ga,tmp,status)
    call spg_native_cell(model,sigma0,b,nb,zb,e,gb,tmp,status)
    h(:,j)=(ga(2:3)-gb(2:3))/4e-7_wp
   enddo
   call solve_linear(h,-g(2:3),delta,2,ok)
   if(.not.ok)exit
   alpha=1
   do ls=1,25
    trial=z;trial(2:3)=z(2:3)+alpha*delta
    call spg_native_cell(model,sigma0,b,nb,trial,e,ga,tmp,status)
    if(status=='ok'.and.maxval(abs(ga(2:3)))<norm0)exit
    alpha=alpha/2
   enddo
   if(ls>25)exit
   z=trial
  enddo
  call spg_native_cell(model,sigma0,b,nb,z,e,g,result,status)
  result%max_projected_gradient=maxval(abs(g(2:3)))
  if(result%max_projected_gradient>1e-5_wp)then
   status='vacuum tangent stationarity failed';return
  endif
  result%stationary=.true.;result%status='VACUUM_TANGENT_STATIONARY';status='ok'
 contains
  subroutine objective(x,f,g,ok)
   real(wp),intent(in)::x(:)
   real(wp),intent(out)::f,g(:)
   logical,intent(out)::ok
   real(wp)::zz(4),gg(4)
   type(spg_native_result_t)::rr
   character(160)::ss
   zz=[1._wp,.1_wp*x(1),.1_wp*x(2),0._wp]
   call spg_native_cell(model,sigma0,b,nb,zz,f,gg,rr,ss)
   g=.1_wp*gg(2:3);ok=ss=='ok'
  end subroutine
 end subroutine

 subroutine spg_native_drip(model,sigma0,b,nb,result,status)
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::sigma0,b
  real(wp),intent(out)::nb
  type(spg_native_result_t),intent(out)::result
  character(*),intent(out)::status
  real(wp)::lo,hi,fl,fh,fm
  integer::it
  lo=1e-6_wp;hi=.002_wp;nb=0
  call vacuum_cell(model,sigma0,b,lo,result,status)
  if(status/='ok')return
  fl=result%gradient(1)
  call vacuum_cell(model,sigma0,b,hi,result,status)
  if(status/='ok')return
  fh=result%gradient(1)
  if(fl*fh>=0)then
   status='neutron drip not bracketed';return
  endif
  do it=1,60
   nb=(lo+hi)/2
   call vacuum_cell(model,sigma0,b,nb,result,status)
   if(status/='ok')return
   fm=result%gradient(1)
   if(fl*fm<=0)then
    hi=nb;fh=fm
   else
    lo=nb;fl=fm
   endif
   if(hi-lo<2e-14_wp)exit
  enddo
  result%max_projected_gradient=maxval(abs(result%gradient(:3)))
  result%stationary=result%max_projected_gradient<1e-5_wp
  result%status='NEUTRON_DRIP_KKT';status='ok'
 end subroutine

 subroutine solve_cluster_branch(model,sigma0,b,nb,result,status,start)
  ! The source endpoint compares a metastable clustered branch with uniform
  ! matter. A trial at the artificial ni=nb floor with u->1, or t->0, belongs to the
  ! analytically known uniform limit, not to that finite-cluster branch.
  ! Keep it in attempt provenance, but do not let it hide an independently
  ! stationary finite cluster. Other lower unconverged trials remain a HOLD.
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::sigma0,b,nb
  type(spg_native_result_t),intent(out)::result
  character(*),intent(out)::status
  real(wp),intent(in),optional::start(4)
  type(spg_native_result_t),allocatable::attempts(:)
  type(spg_native_result_t)::initial
  integer::i,best,nlimit,nlower
  real(wp)::energy
  logical::uniform_limit
  call spg_native_solve(model,sigma0,b,nb,initial,status,start,attempts)
  result=initial;best=0;energy=huge(1._wp);nlimit=0;nlower=0
  do i=1,size(attempts)
   if(.not.attempts(i)%stationary)cycle
   if(attempts(i)%energy_per_baryon<energy)then
    best=i;energy=attempts(i)%energy_per_baryon
   endif
  enddo
  if(best==0)return
  do i=1,size(attempts)
   if(attempts(i)%stationary.or..not.attempts(i)%evaluated)cycle
   uniform_limit=(attempts(i)%z(2)<=nb*(1+2e-6_wp).and.attempts(i)%volume_fraction>1-1e-4_wp).or. &
    (attempts(i)%z(1)<=1.0001e-8_wp.and.attempts(i)%volume_fraction<1e-7_wp)
   if(uniform_limit)then
    nlimit=nlimit+1
   elseif(attempts(i)%energy_per_baryon<energy-1e-9_wp)then
    nlower=nlower+1
   endif
  enddo
  result=attempts(best);result%attempt_count=size(attempts)
  result%uniform_limit_attempts=nlimit;result%unresolved_lower_attempts=nlower
  status='CLUSTER_BRANCH_LOWER_UNCONVERGED_ATTEMPT'
  if(nlower>0)return
  result%status='STATIONARY_CLUSTER_BRANCH_FOR_FIXED_DENSITY_COMPARISON';status='ok'
 end subroutine

 subroutine spg_native_crossing(model,sigma0,b,bracket,cross,cluster,uniform,status)
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::sigma0,b,bracket(2)
  real(wp),intent(out)::cross
  type(spg_native_result_t),intent(out)::cluster
  type(gdfm_native_beta_t),intent(out)::uniform
  character(*),intent(out)::status
  real(wp)::lo,hi,fl,fh,fm
  integer::it
  lo=bracket(1);hi=bracket(2);cross=0
  call delta_at(lo,fl,status)
  if(status/='ok')return
  call delta_at(hi,fh,status)
  if(status/='ok')return
  if(fl*fh>=0)then
   status='fixed-density crossing not bracketed';return
  endif
  do it=1,50
   cross=(lo+hi)/2
   call delta_at(cross,fm,status)
   if(status/='ok')return
   if(fl*fm<=0)then
    hi=cross;fh=fm
   else
    lo=cross;fl=fm
   endif
   if(hi-lo<2e-10_wp)exit
  enddo
 contains
  subroutine delta_at(n,delta,ss)
   real(wp),intent(in)::n
   real(wp),intent(out)::delta
   character(*),intent(out)::ss
   call solve_cluster_branch(model,sigma0,b,n,cluster,ss)
   if(ss/='ok')return
   call gdfm_native_beta(model,n,uniform,ss)
   delta=cluster%epsilon-uniform%epsilon
  end subroutine
 end subroutine

 subroutine spg_native_transition(model,sigma0,b,cross,cluster,uniform,status,bracket_out,evaluations)
  ! Find the source fixed-density energy crossing from independent native cell
  ! evaluations. No tabulated transition density or CompOSE target is an input.
  ! A lost/nonstationary branch reduces the forward step; it is never assigned
  ! zero energy difference or silently interpreted as a physical transition.
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::sigma0,b
  real(wp),intent(out)::cross
  type(spg_native_result_t),intent(out)::cluster
  type(gdfm_native_beta_t),intent(out)::uniform
  character(*),intent(out)::status
  real(wp),intent(out),optional::bracket_out(2)
  integer,intent(out),optional::evaluations
  type(spg_native_result_t)::previous,trial
  type(gdfm_native_beta_t)::beta
  real(wp)::left,right,delta,step,bracket(2)
  integer::i,count
  character(160)::cell_status
  cross=0;count=0;cluster=spg_native_result_t();uniform=gdfm_native_beta_t()
  if(present(bracket_out))bracket_out=0
  if(present(evaluations))evaluations=0
  status='invalid SPG surface parameters'
  if(.not.all(ieee_is_finite([sigma0,b])).or.sigma0<=0.or.b<=-16)return
  left=.002_wp;step=.00025_wp
  call solve_cluster_branch(model,sigma0,b,left,previous,status);count=count+1
  if(status/='ok')return
  call gdfm_native_beta(model,left,beta,status)
  if(status/='ok')return
  status='native transition scan lacks clustered lower endpoint'
  if(previous%epsilon>=beta%epsilon)return
  do i=1,1500
   right=min(.12_wp,left+step)
   call solve_cluster_branch(model,sigma0,b,right,trial,cell_status,start=previous%z)
   count=count+1
   if(present(evaluations))evaluations=count
   if(cell_status/='ok')then
    step=step/2
    if(step<1e-8_wp)then
     cluster=trial;status='native crossing unresolved after branch loss: '//trim(cell_status);return
    endif
    cycle
   endif
   call gdfm_native_beta(model,right,beta,status)
   if(status/='ok')return
   delta=trial%epsilon-beta%epsilon
   if(delta>=0)then
    bracket=[left,right]
    if(present(bracket_out))bracket_out=bracket
    call spg_native_crossing(model,sigma0,b,bracket,cross,cluster,uniform,status)
    return
   endif
   left=right;previous=trial
   if(right>=.12_wp)exit
  enddo
  cluster=previous;uniform=beta
  status='native fixed-density crossing not bracketed within support'
 end subroutine

 subroutine spg_native_endpoint(model,sigma0,b,bracket,n_cluster,n_uniform,cross,status)
  ! Explicit alternative common-tangent diagnostic, distinct from source crossing.
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::sigma0,b,bracket(2)
  real(wp),intent(out)::n_cluster,n_uniform,cross
  character(*),intent(out)::status
  real(wp)::x(2),f(2),fa(2),fb(2),jac(2,2),delta(2),xa(2),xb(2),trial(2),alpha,norm0
  integer::it,j,ls
  logical::ok
  type(spg_native_result_t)::cl
  type(gdfm_native_beta_t)::un
  call spg_native_crossing(model,sigma0,b,bracket,cross,cl,un,status)
  if(status/='ok')return
  x=[cross*.995_wp,cross*1.005_wp]
  n_cluster=0;n_uniform=0
  do it=1,30
   call residual(x,f,status)
   if(status/='ok')return
   norm0=maxval(abs(f))
   if(norm0<1e-7_wp)exit
   do j=1,2
    xa=x;xb=x;xa(j)=xa(j)+2e-7_wp;xb(j)=xb(j)-2e-7_wp
    call residual(xa,fa,status)
    if(status/='ok')return
    call residual(xb,fb,status)
    if(status/='ok')return
    jac(:,j)=(fa-fb)/4e-7_wp
   enddo
   call solve_linear(jac,-f,delta,2,ok)
   if(.not.ok)then
    status='singular common tangent Jacobian';return
   endif
   alpha=1
   do ls=1,20
    trial=x+alpha*delta
    if(minval(trial)>bracket(1)*.98_wp.and.maxval(trial)<bracket(2)*1.02_wp)then
     call residual(trial,fa,status)
     if(status=='ok'.and.maxval(abs(fa))<norm0)exit
    endif
    alpha=alpha/2
   enddo
   if(ls>20)then
    status='common tangent line search failed';return
   endif
   x=trial
  enddo
  if(maxval(abs(f))>1e-6_wp.or.x(2)<=x(1))then
   status='common tangent unresolved';return
  endif
  n_cluster=x(1);n_uniform=x(2);status='ok'
 contains
  subroutine residual(v,r,ss)
   real(wp),intent(in)::v(2)
   real(wp),intent(out)::r(2)
   character(*),intent(out)::ss
   call spg_native_solve(model,sigma0,b,v(1),cl,ss)
   if(ss/='ok')return
   call gdfm_native_beta(model,v(2),un,ss)
   r=[cl%mu_b-un%mu_n,10*(cl%pressure-un%pressure)]
  end subroutine
 end subroutine
end module
