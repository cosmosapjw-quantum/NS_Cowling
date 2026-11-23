! Attach the finite native BPS outer-crust and NV1973 equilibrium-fit tables.
! The intervening cold, scalar barotrope is an interpolation assumption.
program bps_nv1973_bridge
 use compose_eos_types, only: wp, eos_data_t, allocate_eos
 use compose_eos_io, only: read_compose, write_compose
 implicit none
 integer,parameter :: panels=64
 type(eos_data_t) :: outer,inner,combined,readback
 character(1024) :: outer_path,inner_path,target,status
 integer :: i,j,r,k,m,u
 real(wp) :: n0,n1,e0,e1,mu0,mu1,h,secant,a,b,curv0,curv1,csbound, &
             x,n,e,mu,p,upper_pressure
 if(command_argument_count()/=3)error stop 'usage: bps_nv1973_bridge BPS_DIR NV_DIR EXISTING_NEW_OUTPUT_DIR'
 call get_command_argument(1,outer_path)
 call get_command_argument(2,inner_path)
 call get_command_argument(3,target)
 call read_compose(trim(outer_path),outer,status);call check(status)
 call read_compose(trim(inner_path),inner,status);call check(status)
 if(outer%n/=623.or.inner%n/=257)error stop 'unexpected frozen source grid count'
 if(outer%mn/=inner%mn.or.outer%mp/=inner%mp.or..not.outer%composition_available.or.inner%composition_available) &
  error stop 'unexpected source mass or composition convention'
 if(.not.all(outer%iphase(2:)>=outer%iphase(:outer%n-1)).or. &
    outer%iphase(1)/=201.or.outer%iphase(outer%n)/=212.or.any(inner%iphase/=0))error stop 'source phase chain mismatch'
 if(.not.outer%beta_equilibrated.or..not.inner%beta_equilibrated)error stop 'requires two cold equilibrium sources'
 m=outer%n
 n0=outer%nb(m);e0=n0*(1+outer%q(7,m))*outer%mn
 mu0=(1+outer%q(3,m))*outer%mn
 if(abs(mu0-outer%mn)>1e-8_wp)error stop 'outer source does not terminate at drip'
 open(newunit=u,file=trim(target)//'/join_candidates.csv',status='new',action='write')
 write(u,'(a)')'NV_row,n,epsilon,mu,P,mu_minus_secant,curvature_left_MeV_fm3,curvature_right_MeV_fm3,cs2_bound,admissible'
 do r=1,inner%n
  n1=inner%nb(r);e1=n1*(1+inner%q(7,r))*inner%mn
  mu1=(1+inner%q(3,r))*inner%mn
  h=n1-n0
  if(h<=0)error stop 'NV row outside positive-width bridge'
  secant=(e1-e0)/h
  ! epsilon=e0+h*(mu0*x+a*x**2+b*x**3), x=(n-n0)/h.
  a=3*(secant-mu0)-(mu1-mu0)
  b=(mu1-mu0)-2*(secant-mu0)
  curv0=2*a/h;curv1=(2*a+6*b)/h
  csbound=n1*max(curv0,curv1)/mu0
  write(u,'(i0,8(",",es24.15e3),",",l1)')r,n1,e1,mu1,n1*inner%q(1,r),mu1-secant,curv0,curv1,csbound, &
    min(curv0,curv1)>0.and.csbound<1
  if(min(curv0,curv1)>0.and.csbound<1)exit
 enddo
 close(u)
 if(r>inner%n)error stop 'no stable subluminal Hermite bridge'
 if(abs(e0+h*(mu0+a+b)-e1)>1e-12_wp.or.abs(mu0+2*a+3*b-mu1)>1e-8_wp)error stop 'C1 attachment failed'
 upper_pressure=n1*mu1-e1
 if(abs(upper_pressure-n1*inner%q(1,r))>1e-12_wp.or.abs(n0*mu0-e0-n0*outer%q(1,m))>1e-12_wp) &
  error stop 'pressure endpoint mismatch'
 call allocate_eos(combined,m+panels-1+inner%n-r+1)
 combined%mn=outer%mn;combined%mp=outer%mp
 combined%composition_file_present=.true.;combined%composition_available=.true.;combined%beta_equilibrated=.true.
 combined%nb(:m)=outer%nb;combined%q(:,:m)=outer%q;combined%iphase(:m)=outer%iphase
 combined%yn(:m)=outer%yn;combined%yp(:m)=outer%yp;combined%ymu(:m)=outer%ymu;combined%ye(:m)=outer%ye
 combined%abar(:m)=outer%abar;combined%zbar(:m)=outer%zbar;combined%yh(:m)=outer%yh
 combined%yn_present(:m)=outer%yn_present;combined%yp_present(:m)=outer%yp_present
 combined%ye_present(:m)=outer%ye_present;combined%ymu_present(:m)=outer%ymu_present
 combined%abar_present(:m)=outer%abar_present;combined%zbar_present(:m)=outer%zbar_present
 combined%yh_present(:m)=outer%yh_present;combined%row_status(:m)=outer%row_status
 combined%c3_applicable(:m)=outer%c3_applicable;combined%c4_applicable(:m)=outer%c4_applicable
 do i=1,panels-1
  x=real(i,wp)/panels;n=n0+h*x
  e=e0+h*x*(mu0+x*(a+b*x));mu=mu0+x*(2*a+3*b*x)
  p=n*mu-e
  k=m+i
  combined%nb(k)=n;combined%iphase(k)=213
  combined%q(1,k)=p/n;combined%q(3,k)=mu/combined%mn-1
  combined%q(6,k)=e/n/combined%mn-1;combined%q(7,k)=combined%q(6,k)
 enddo
 j=m+panels
 combined%nb(j:)=inner%nb(r:);combined%q(:,j:)=inner%q(:,r:);combined%iphase(j:)=214
 call write_compose(trim(target),combined,status);call check(status)
 call read_compose(trim(target),readback,status);call check(status)
 if(readback%n/=combined%n.or.any(readback%iphase/=combined%iphase).or. &
   any(abs(readback%q-combined%q)>1e-14_wp).or.any(readback%nb/=combined%nb))error stop 'native readback mismatch'
 write(*,'(a,i0,a,i0,a,i0)')'rows=',combined%n,' first_NV_row=',r,' bridge_interior=',panels-1
 write(*,'(a,es24.16)')'global upper bound on c_s^2/c^2 over bridge=',csbound
contains
 subroutine check(s)
  character(*),intent(in)::s
  if(trim(s)/='ok')then
   print *,trim(s)
   error stop 'native I/O failed'
  endif
 end subroutine
end program
