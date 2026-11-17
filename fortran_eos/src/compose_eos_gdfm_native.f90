! Author-parameterized GDFM bulk, Scurto 2024 Eqs.1-20.
! Exact parameter rows: data/research_generator_20261116/spg_author_candidates.json.
! Density fm^-3; energies MeV; hbar*c retained. Gamma_rho is the archived vertex,
! equivalent to paper g_rho/2; no EOS table is loaded or interpolated.
module compose_eos_gdfm_native
 use compose_eos_types, only: wp
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 real(wp),parameter,public :: gdfm_hbarc=197.3269804_wp, gdfm_mass=938.9_wp
 real(wp),parameter :: pi=acos(-1._wp), me=.51099895_wp, mmu=105.6583755_wp
 type,public :: gdfm_native_model_t
  character(16) :: id=''
  real(wp) :: n0=0._wp, masses(3)=0._wp, coefficients(4,3)=0._wp
 end type
 type,public :: gdfm_native_beta_t
  real(wp):: nb=0, yp=0, ye=0, ymu=0, epsilon=0, pressure=0, mu_n=0,mu_p=0,mu_e=0
  real(wp):: charge_residual=0, scalar_residual=0
 end type
 public :: load_gdfm_native,gdfm_native_baryons,gdfm_native_beta,native_fermion
contains
 subroutine load_gdfm_native(id,model,status)
  character(*),intent(in)::id
  type(gdfm_native_model_t),intent(out)::model
  character(*),intent(out)::status
  model=gdfm_native_model_t(); model%id=id
  model%masses=[560._wp,782.6_wp,769._wp]
  select case(trim(id))
  case('SPG_M1')
   model%n0=.1572_wp
   model%coefficients(:,1)=[9.8447_wp,2.0646_wp,2.633_wp,3.8992_wp]
   model%coefficients(:,2)=[11.531_wp,2.3992_wp,2.546_wp,8.4805_wp]
   model%coefficients(:,3)=[1.969310528889322_wp,3.9153_wp,.5512318936407548_wp,-.6028_wp]
  case('SPG_M2')
   model%n0=.1551_wp
   model%coefficients(:,1)=[9.4639_wp,2.3656_wp,2.5004_wp,3.544_wp]
   model%coefficients(:,2)=[11.0713_wp,2.3436_wp,2.3324_wp,6.2542_wp]
   model%coefficients(:,3)=[2.1072879846410077_wp,4.8314_wp,1.004343531768306_wp,-.9447_wp]
  case('SPG_M3')
   model%n0=.1612_wp
   model%coefficients(:,1)=[9.9562_wp,1.9862_wp,2.949_wp,4.0312_wp]
   model%coefficients(:,2)=[11.5983_wp,2.3833_wp,2.6646_wp,8.9181_wp]
   model%coefficients(:,3)=[.49636167129873143_wp,5.2825_wp,.3765226791962973_wp,-.5474_wp]
  case('SPG_M4')
   model%n0=.1621_wp
   model%coefficients(:,1)=[9.2457_wp,2.2078_wp,2.6503_wp,4.0604_wp]
   model%coefficients(:,2)=[10.7675_wp,2.0646_wp,2.3546_wp,6.1144_wp]
   model%coefficients(:,3)=[1.377335089991927_wp,5.1015_wp,.7495296034394823_wp,-.6212_wp]
  case('SPG_M5')
   model%n0=.1483_wp
   model%coefficients(:,1)=[10.6344_wp,1.853_wp,3.0812_wp,3.6036_wp]
   model%coefficients(:,2)=[12.5265_wp,2.2047_wp,2.7491_wp,9.0012_wp]
   model%coefficients(:,3)=[1.510450686691207_wp,4.6623_wp,.6383881919207285_wp,-.5763_wp]
  case default
   status='unknown GDFM model';return
  end select
  status='ok'
 end subroutine

 pure subroutine native_fermion(n,m,energy,pressure,chemical,scalar)
  real(wp),intent(in)::n,m
  real(wp),intent(out)::energy,pressure,chemical
  real(wp),intent(out),optional::scalar
  real(wp)::k,x,x2,root,a,ce,cs,pe,ps,se,ss,te,ts,pref
  integer::j
  energy=0;pressure=0;chemical=m
  if(present(scalar))scalar=0
  if(n<=0 .or. m<=0)return
  k=gdfm_hbarc*(3*pi*pi*n)**(1._wp/3)
  x=k/m; x2=x*x; root=sqrt(1+x2);chemical=m*root
  if(x<.15_wp)then
   ! Uniformly convergent low-momentum expansion avoids threshold cancellation.
   ce=1;cs=1;pe=1;ps=1;se=1;ss=1
   do j=1,20
    ce=ce*(1.5_wp-j)/j;cs=cs*(.5_wp-j)/j
    pe=pe*x2;ps=ps*x2
    te=ce*pe*3/(2*j+3);ts=cs*ps*3/(2*j+3)
    se=se+te;ss=ss+ts
    if(max(abs(te),abs(ts))<1e-18_wp)exit
   end do
   energy=m*n*se
   if(present(scalar))scalar=n*ss
  else
   a=asinh(x);pref=m**4/(8*pi*pi*gdfm_hbarc**3)
   energy=pref*(x*(2*x2+1)*root-a)
   if(present(scalar))scalar=m**3/(2*pi*pi*gdfm_hbarc**3)*(x*root-a)
  end if
  pressure=n*chemical-energy
 end subroutine

 subroutine gdfm_native_baryons(model,n,yp,eps,mun,mup,p,status,scalar_residual)
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::n,yp
  real(wp),intent(out)::eps,mun,mup,p
  character(*),intent(out)::status
  real(wp),intent(out),optional::scalar_residual
  real(wp)::x,ga(3),dg(3),co(3),dc(3),a,b,c,d,ex
  real(wp)::nn,np,n3,lo,hi,s,res,mstar,en,ep,pn,pp,kn,kp,sn,sp,rearr
  integer::i,it
  eps=0;p=0;mun=gdfm_mass;mup=gdfm_mass
  if(present(scalar_residual))scalar_residual=0
  status='invalid GDFM input'
  if(.not.all(ieee_is_finite([n,yp,model%n0])))return
  if(n<0.or.yp<0.or.yp>1.or.model%n0<=0)return
  if(n==0)then
   status='ok';return
  end if
  x=n/model%n0
  do i=1,3
   a=model%coefficients(1,i);b=model%coefficients(2,i)
   c=model%coefficients(3,i);d=model%coefficients(4,i)
   ex=exp(-c*x);ga(i)=a+(b+d*x**3)*ex
   dg(i)=(3*d*x*x-c*(b+d*x**3))*ex/model%n0
  end do
  co=(gdfm_hbarc/model%masses)**2*ga*ga
  dc=2*(gdfm_hbarc/model%masses)**2*ga*dg
  nn=n*(1-yp);np=n*yp;n3=np-nn
  lo=0;hi=gdfm_mass/gdfm_hbarc*(1-1e-12_wp)
  do it=1,100
   s=(lo+hi)/2
   mstar=gdfm_mass-s*gdfm_hbarc
   call native_fermion(nn,mstar,en,pn,kn,sn)
   call native_fermion(np,mstar,ep,pp,kp,sp)
   res=s/co(1)-sn-sp
   if(res>0)then
    hi=s
   else
    lo=s
   endif
   if(hi-lo<2e-15_wp)exit
  enddo
  eps=en+ep+gdfm_hbarc*(s*s/co(1)+co(2)*n*n+co(3)*n3*n3)/2
  rearr=gdfm_hbarc*(dc(2)*n*n+dc(3)*n3*n3-dc(1)*(s/co(1))**2)/2
  mun=kn+gdfm_hbarc*(co(2)*n-co(3)*n3)+rearr
  mup=kp+gdfm_hbarc*(co(2)*n+co(3)*n3)+rearr
  p=nn*mun+np*mup-eps
  if(present(scalar_residual))scalar_residual=res
  status='ok'
 end subroutine

 subroutine gdfm_native_beta(model,n,state,status,include_muons)
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::n
  type(gdfm_native_beta_t),intent(out)::state
  character(*),intent(out)::status
  logical,intent(in),optional::include_muons
  real(wp)::lo,hi,y,ep,pp,mun,mup,chem,ne,nmu,ee,em,pe,pm,tmp,res
  logical::muons
  integer::it
  state=gdfm_native_beta_t();muons=.true.
  if(present(include_muons))muons=include_muons
  status='invalid beta density'
  if(.not.ieee_is_finite(n).or.n<=0)return
  lo=1e-12_wp;hi=.5_wp
  do it=1,90
   y=(lo+hi)/2
   call gdfm_native_baryons(model,n,y,ep,mun,mup,pp,status,res)
   if(status/='ok')return
   chem=mun-mup
   ne=0;nmu=0
   if(chem>me)ne=(sqrt((chem-me)*(chem+me))/gdfm_hbarc)**3/(3*pi*pi)
   if(muons.and.chem>mmu)nmu=(sqrt((chem-mmu)*(chem+mmu))/gdfm_hbarc)**3/(3*pi*pi)
   if(ne+nmu>n*y)then
    lo=y
   else
    hi=y
   endif
   if(hi-lo<2e-14_wp)exit
  enddo
  call native_fermion(ne,me,ee,pe,tmp)
  call native_fermion(nmu,mmu,em,pm,tmp)
  state%nb=n;state%yp=y;state%ye=ne/n;state%ymu=nmu/n
  state%epsilon=ep+ee+em;state%pressure=pp+pe+pm
  state%mu_n=mun;state%mu_p=mup;state%mu_e=chem
  state%charge_residual=n*y-ne-nmu;state%scalar_residual=res
  status='ok'
 end subroutine
end module
