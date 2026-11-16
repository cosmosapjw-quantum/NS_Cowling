! Source-defined GPPVA spatial RMF kernels. MeV, fm; rho vertex g_rho*tau3/2.
! Formula/parameter mapping: python_frontend/compose_eos/gppva_rmf.py at 4de4943.
module compose_eos_tf_models
 use compose_eos_types, only: wp
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 real(wp), parameter, public :: tf_hbarc=197.3269804_wp,tf_me=.51099895_wp,tf_mmu=105.6583755_wp
 real(wp), parameter, public :: tf_pi=acos(-1._wp),tf_e2=1.43996448_wp,tf_hc3=tf_hbarc**3
 type, public :: tf_model_t
  character(32) :: id=''
  logical :: density_dependent=.false.
  real(wp)::mn=939._wp,mp=939._wp,n0=0._wp,ms=0._wp,mw=0._wp,mr=0._wp
  real(wp)::gs=0._wp,gw=0._wp,gr=0._wp
  real(wp)::sa=0._wp,sb=0._wp,sc=0._wp,sd=0._wp,wa=0._wp,wb=0._wp,wc=0._wp,wd=0._wp,ra=0._wp
  real(wp)::a3=0._wp,a4=0._wp,b4=0._wp,cross=0._wp
 end type
 public :: load_tf_model,tf_couplings,tf_scalar_potential,tf_vector_potential,tf_fermion
contains
 subroutine load_tf_model(name,p,status)
 character(*),intent(in)::name
 type(tf_model_t),intent(out)::p
 character(*),intent(out)::status
 character(32)::id
 real(wp)::kap,lam,zeta,cross
 p=tf_model_t();id=trim(name);if(index(id,'GPPVA_')==1)id=id(7:)
 p%id='GPPVA_'//trim(id);status='ok'
 select case(trim(id))
 case('DD2')
  p%mn=939.56536_wp;p%mp=938.27203_wp;p%n0=.149065_wp;p%ms=546.212459_wp;p%mw=783;p%mr=763
  p%gs=10.686681_wp;p%gw=13.342362_wp;p%gr=2*3.626940_wp
  p%sa=1.357630_wp;p%sb=.634442_wp;p%sc=1.005358_wp;p%sd=.575810_wp
  p%wa=1.369718_wp;p%wb=.496475_wp;p%wc=.817753_wp;p%wd=.638452_wp;p%ra=.518903_wp
 case('DDME2')
  p%n0=.152_wp;p%ms=550.1238_wp;p%mw=783;p%mr=763
  p%gs=10.5396_wp;p%gw=13.0189_wp;p%gr=2*3.6836_wp
  p%sa=1.3881_wp;p%sb=1.0943_wp;p%sc=1.7057_wp;p%sd=.4421_wp
  p%wa=1.3892_wp;p%wb=.9240_wp;p%wc=1.4620_wp;p%wd=.4775_wp;p%ra=.5647_wp
 case('TW')
  p%n0=.153_wp;p%ms=550;p%mw=783;p%mr=763
  p%gs=10.72854_wp;p%gw=13.29015_wp;p%gr=2*3.66098_wp
  p%sa=1.365469_wp;p%sb=.226061_wp;p%sc=.409704_wp;p%sd=.901995_wp
  p%wa=1.402488_wp;p%wb=.172577_wp;p%wc=.344293_wp;p%wd=.983955_wp;p%ra=.515_wp
 case('FSU2','FSU2R','FSU2H')
  p%ms=497.479_wp;p%mw=782.5_wp;p%mr=763
  select case(trim(id))
  case('FSU2')
   p%gs=sqrt(108.0943_wp);p%gw=sqrt(183.7893_wp);p%gr=sqrt(80.4656_wp)
   kap=3.0029_wp;lam=-.000533_wp;zeta=.0256_wp;cross=.000823_wp
  case('FSU2R')
   p%gs=sqrt(107.5751_wp);p%gw=sqrt(182.3949_wp);p%gr=sqrt(206.4260_wp)
   kap=3.0911_wp;lam=-.001680_wp;zeta=.024_wp;cross=.045_wp
  case default
   p%gs=sqrt(102.7200_wp);p%gw=sqrt(169.5315_wp);p%gr=sqrt(197.2692_wp)
   kap=4.0014_wp;lam=-.013298_wp;zeta=.008_wp;cross=.045_wp
  end select
  p%a3=kap*p%gs**3/2;p%a4=lam*p%gs**4/6;p%b4=zeta*p%gw**4/6;p%cross=cross*p%gw**2*p%gr**2
 case('NL3wrL55')
  p%ms=508.194_wp;p%mw=782.501_wp;p%mr=763;p%gs=10.217_wp;p%gw=12.868_wp;p%gr=11.2766_wp
  p%a3=4.384_wp*939/2;p%a4=-173.31_wp/6;p%cross=.03_wp*p%gw**2*p%gr**2
 case('TM1e')
  p%mn=938;p%mp=938;p%ms=511.198_wp;p%mw=783;p%mr=770
  p%gs=10.0289_wp;p%gw=12.6139_wp;p%gr=13.9714_wp
  p%a3=7.2325_wp*tf_hbarc;p%a4=.6183_wp;p%b4=71.3075_wp;p%cross=.0429_wp*p%gw**2*p%gr**2
 case default
  status='unknown GPPVA spatial TF model'
 end select
 p%density_dependent=p%n0>0
 end subroutine

 pure subroutine tf_couplings(p,n,gs,gw,gr,dgs,dgw,dgr)
 type(tf_model_t),intent(in)::p
 real(wp),intent(in)::n
 real(wp),intent(out)::gs,gw,gr,dgs,dgw,dgr
 real(wp)::x,u
 gs=p%gs;gw=p%gw;gr=p%gr;dgs=0;dgw=0;dgr=0
 if(.not.p%density_dependent)return
 x=n/p%n0;u=x+p%sd
 gs=p%gs*p%sa*(1+p%sb*u*u)/(1+p%sc*u*u)
 dgs=p%gs*2*p%sa*(p%sb-p%sc)*u/(1+p%sc*u*u)**2/p%n0
 u=x+p%wd
 gw=p%gw*p%wa*(1+p%wb*u*u)/(1+p%wc*u*u)
 dgw=p%gw*2*p%wa*(p%wb-p%wc)*u/(1+p%wc*u*u)**2/p%n0
 gr=p%gr*exp(-p%ra*(x-1));dgr=-p%ra*gr/p%n0
 end subroutine

 pure subroutine tf_scalar_potential(p,s,u,du,ddu)
 type(tf_model_t),intent(in)::p
 real(wp),intent(in)::s
 real(wp),intent(out)::u,du,ddu
 u=p%ms**2*s*s/2+p%a3*s**3/3+p%a4*s**4/4
 du=p%ms**2*s+p%a3*s*s+p%a4*s**3
 ddu=p%ms**2+2*p%a3*s+3*p%a4*s*s
 end subroutine
 pure subroutine tf_vector_potential(p,w,r,v,dw,dr,dww,drr,dwr)
 type(tf_model_t),intent(in)::p
 real(wp),intent(in)::w,r
 real(wp),intent(out)::v,dw,dr,dww,drr,dwr
 v=p%mw**2*w*w/2+p%mr**2*r*r/2+p%b4*w**4/4+p%cross*w*w*r*r
 dw=p%mw**2*w+p%b4*w**3+2*p%cross*w*r*r
 dr=p%mr**2*r+2*p%cross*w*w*r
 dww=p%mw**2+3*p%b4*w*w+2*p%cross*r*r
 drr=p%mr**2+2*p%cross*w*w;dwr=4*p%cross*w*r
 end subroutine

 pure subroutine tf_fermion(n,mass,energy,scalar,ds,ef)
 real(wp),intent(in)::n,mass
 real(wp),intent(out)::energy,scalar,ds,ef
 real(wp)::k,x,a
 k=tf_hbarc*(3*tf_pi**2*max(n,0._wp))**(1._wp/3);x=k/mass;ef=sqrt(k*k+mass*mass)
 if(x<.08_wp)then
  scalar=n*(1-3*x*x/10+9*x**4/56-5*x**6/48+105*x**8/1408)
  ds=n/mass*(3*x*x/5-9*x**4/14+5*x**6/8-105*x**8/176)
  energy=n*mass*(1+3*x*x/10-3*x**4/56+x**6/48-15*x**8/1408)
 else
  a=asinh(x)
  scalar=mass*(k*ef-mass**2*a)/(2*tf_pi**2*tf_hc3)
  ds=(k*ef+2*mass**2*k/ef-3*mass**2*a)/(2*tf_pi**2*tf_hc3)
  energy=(k*ef*(2*k*k+mass*mass)-mass**4*a)/(8*tf_pi**2*tf_hc3)
 endif
 end subroutine
end module
