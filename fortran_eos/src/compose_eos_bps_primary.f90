! BPS1971 Eqs. (1),(3)-(5),(14)-(19), DOI 10.1086/151216.
! Published Table 2 sequence and rounded binding energies, NOT the entire 1965
! Myers-Swiatecki candidate mass catalogue. No modern atomic/mass corrections.
! Units: fm, MeV, MeV/fm**3; c and hbar are retained in dimensional constants.
module compose_eos_bps_primary
 use compose_eos_types, only: wp,eos_data_t,allocate_eos
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 real(wp),parameter,public :: BPS_MN=939.56542052_wp,BPS_MP=938.27208816_wp
 real(wp),parameter,public :: BPS_ME=.51099895_wp,BPS_HBARC=197.3269804_wp
 real(wp),parameter,public :: BPS_RHO_PER_MEVFM3=1.7826619216278976e12_wp
 real(wp),parameter,public :: BPS_P_CGS_PER_MEVFM3=1.602176634e33_wp
 real(wp),parameter :: PI=acos(-1._wp),E2=BPS_HBARC/137.035999084_wp
 real(wp),parameter :: NB_MIN=6.295e-12_wp
 integer,parameter :: NC=12
 integer,parameter :: EXPECT_Z(NC)=[26,28,28,34,32,30,28,26,42,40,38,36]
 integer,parameter :: EXPECT_A(NC)=[56,62,64,84,82,80,78,76,124,122,120,118]
 type,public :: bps_catalogue_t
  integer :: n=0
  integer :: z(NC)=0,a(NC)=0
  real(wp) :: binding(NC)=0,rho_max(NC)=0,mu_e_max(NC)=0,jump_percent(NC)=0
 end type
 type,public :: bps_state_t
  real(wp) :: nb=0,ne=0,epsilon=0,pressure=0,mu_b=0,mu_e=0
  real(wp) :: kinetic_lattice=0,sound2=0,gamma=0,pressure_residual=0
  integer :: isotope=0,z=0,a=0
 end type
 public :: load_bps_primary_catalogue,bps_isotope_state,bps_state_at_pressure
 public :: generate_bps_primary_native
contains
 subroutine load_bps_primary_catalogue(path,cat,status)
  character(*),intent(in)::path
  type(bps_catalogue_t),intent(out)::cat
  character(*),intent(out)::status
  character(4096)::line
  character(256)::token
  integer::u,ios,i,j,k,start,cols,physical_line
  real(wp)::v(6)
  cat=bps_catalogue_t();status='BPS cannot open primary Table 2 CSV'
  open(newunit=u,file=path,status='old',action='read',iostat=ios)
  if(ios/=0)return
  physical_line=0;i=0
  do
   read(u,'(a)',iostat=ios)line
   if(ios<0)exit
   physical_line=physical_line+1
   if(ios/=0.or.len_trim(line)==len(line))then
    status='BPS primary CSV unreadable or overlong line';close(u);return
   endif
   line=adjustl(line)
   if(len_trim(line)==0)cycle
   if(line(1:1)=='#')cycle
   if(i==0)then
    if(trim(line)/='Z,A,b_MeV_per_baryon,rho_max_g_cm3,mu_e_max_MeV,delta_rho_percent')then
     status='BPS primary Table 2 CSV header mismatch';close(u);return
    endif
    i=1;cycle
   endif
   if(i>NC)then
    status='BPS Table 2 requires exactly 12 isotope rows';close(u);return
   endif
   cols=1
   do j=1,len_trim(line)
    if(line(j:j)==',')cols=cols+1
   enddo
   if(cols/=6)then
    status='BPS primary CSV requires exactly six numeric fields';close(u);return
   endif
   start=1
   do j=1,6
    k=index(line(start:),',')
    if(j==6)k=len_trim(line)-start+2
    if(k<=1.or.k>len(token))then
     status='BPS primary CSV empty or overlong numeric field';close(u);return
    endif
    token=adjustl(line(start:start+k-2))
    if(len_trim(token)==0.or.verify(trim(token),'0123456789+-.eEdD')/=0)then
     status='BPS primary CSV invalid numeric token';close(u);return
    endif
    read(token,*,iostat=ios)v(j)
    if(ios/=0)then;status='BPS primary CSV invalid number';close(u);return;endif
    start=start+k
   enddo
   if(.not.all(ieee_is_finite(v)))then
    status='BPS primary CSV nonfinite input';close(u);return
   endif
   if(v(1)/=real(EXPECT_Z(i),wp).or.v(2)/=real(EXPECT_A(i),wp))then
    status='BPS primary CSV isotope sequence differs from Table 2';close(u);return
   endif
   if(v(3)<=0.or.v(3)>=20.or.v(4)<=0.or.v(5)<=0)then
    status='BPS primary CSV nonphysical isotope input';close(u);return
   endif
   if((i<NC.and.v(6)<=0).or.(i==NC.and.v(6)/=-1._wp))then
    status='BPS primary CSV jump field invalid (-1 marks drip only)';close(u);return
   endif
   cat%z(i)=nint(v(1));cat%a(i)=nint(v(2));cat%binding(i)=v(3)
   cat%rho_max(i)=v(4);cat%mu_e_max(i)=v(5);cat%jump_percent(i)=v(6)
   i=i+1
  enddo
  close(u)
  if(i/=NC+1)then;status='BPS Table 2 requires exactly 12 isotope rows';return;endif
  cat%n=NC;status='ok'
 end subroutine

 pure subroutine isotope_evaluate(nb,z,a,b,state)
  real(wp),intent(in)::nb,b
  integer,intent(in)::z,a
  type(bps_state_t),intent(out)::state
  real(wp)::ne,ni,x,x2,g,fac,ee,pe,ke,wl,enuc,lat,dpdnlog
  state=bps_state_t()
  ne=nb*real(z,wp)/a;ni=nb/a
  x=BPS_HBARC*(3*PI**2*ne)**(1._wp/3)/BPS_ME
  x2=x*x;g=sqrt(1+x2);fac=BPS_ME**4/(PI**2*BPS_HBARC**3)
  if(x<.1_wp)then
   ! Convergent nonrelativistic series avoids subtraction of O(x) terms.
   ke=fac*x**5*(.1_wp+x2*(-1._wp/56+x2*(1._wp/144+x2*(-5._wp/1408+ &
      x2*(7._wp/3328-x2*21._wp/15360)))))
   pe=fac*x**5/3*(.2_wp+x2*(-1._wp/14+x2*(1._wp/24+x2*(-5._wp/176+ &
      x2*(35._wp/1664-x2*63._wp/3840)))))
   ee=ne*BPS_ME+ke
  else
   ee=fac/8*(x*(1+2*x2)*g-asinh(x))
   pe=fac/24*(x*(2*x2-3)*g+3*asinh(x))
   ke=ee-ne*BPS_ME
  endif
  wl=-1.819620_wp*real(z,wp)**2*E2*(ni/2)**(1._wp/3)
  lat=ni*wl
  enuc=(a-z)*BPS_MN+z*BPS_MP-a*b
  state%nb=nb;state%ne=ne;state%z=z;state%a=a
  state%epsilon=ni*enuc+ee+lat
  state%kinetic_lattice=ke+lat
  state%pressure=pe+lat/3
  state%mu_e=BPS_ME*g
  state%mu_b=(enuc+4*wl/3+z*state%mu_e)/a
  dpdnlog=ne*BPS_ME*x2/(3*g)+4*lat/9
  state%sound2=dpdnlog/(state%epsilon+state%pressure)
  if(state%pressure/=0)state%gamma=dpdnlog/state%pressure
 end subroutine

 subroutine bps_isotope_state(nb,z,a,b,state,status)
  real(wp),intent(in)::nb,b
  integer,intent(in)::z,a
  type(bps_state_t),intent(out)::state
  character(*),intent(out)::status
  state=bps_state_t();status='BPS invalid isotope density or binding energy'
  if(.not.all(ieee_is_finite([nb,b])))return
  if(nb<=0.or.nb>1._wp.or.b<=0.or.b>=20.or.z<=0.or.a<z)return
  call isotope_evaluate(nb,z,a,b,state)
  if(.not.all(ieee_is_finite([state%epsilon,state%pressure,state%mu_b,state%sound2])))return
  status='ok'
 end subroutine

 subroutine isotope_pressure(cat,i,p,state,status)
  type(bps_catalogue_t),intent(in)::cat
  integer,intent(in)::i
  real(wp),intent(in)::p
  type(bps_state_t),intent(out)::state
  character(*),intent(out)::status
  real(wp)::lo,hi,mid
  type(bps_state_t)::s
  integer::j
  state=bps_state_t();status='BPS invalid fixed-pressure request'
  if(cat%n/=NC.or.i<1.or.i>NC.or..not.ieee_is_finite(p))return
  if(p<=0.or.p>1e-2_wp)return
  lo=log(1e-16_wp);hi=log(1e-2_wp)
  call isotope_evaluate(exp(lo),cat%z(i),cat%a(i),cat%binding(i),s)
  if(s%pressure>=p)then;status='BPS pressure below numerical bracket';return;endif
  call isotope_evaluate(exp(hi),cat%z(i),cat%a(i),cat%binding(i),s)
  if(s%pressure<=p)then;status='BPS pressure above numerical bracket';return;endif
  do j=1,64
   mid=(lo+hi)/2
   call isotope_evaluate(exp(mid),cat%z(i),cat%a(i),cat%binding(i),s)
   if(s%pressure<p)then;lo=mid;else;hi=mid;endif
   if(hi-lo<2e-14_wp)exit
  enddo
  call isotope_evaluate(exp((lo+hi)/2),cat%z(i),cat%a(i),cat%binding(i),state)
  state%isotope=i;state%pressure_residual=abs(state%pressure/p-1)
  ! Preserve the computed energy. The common pressure is the solved control.
  state%pressure=p
  if(state%pressure_residual>1e-10_wp)then;status='BPS pressure inversion failed';return;endif
  status='ok'
 end subroutine

 subroutine bps_state_at_pressure(cat,p,state,status)
  type(bps_catalogue_t),intent(in)::cat
  real(wp),intent(in)::p
  type(bps_state_t),intent(out)::state
  character(*),intent(out)::status
  type(bps_state_t)::s
  integer::i
  state=bps_state_t();state%mu_b=huge(1._wp)
  do i=1,NC
   call isotope_pressure(cat,i,p,s,status)
   if(status/='ok')return
   if(s%mu_b<state%mu_b)state=s
  enddo
  status='ok'
 end subroutine

 subroutine transition_pressure(cat,i,pmin,pmax,p,left,right,status)
  type(bps_catalogue_t),intent(in)::cat
  integer,intent(in)::i
  real(wp),intent(in)::pmin,pmax
  real(wp),intent(out)::p
  type(bps_state_t),intent(out)::left,right
  character(*),intent(out)::status
  type(bps_state_t)::global
  real(wp)::lo,hi,mid,delta
  integer::j
  lo=log(pmin);hi=log(pmax)
  call isotope_pressure(cat,i,pmin,left,status);if(status/='ok')return
  call isotope_pressure(cat,i+1,pmin,right,status);if(status/='ok')return
  if(left%mu_b>=right%mu_b)then;status='BPS transition lacks low-P Gibbs bracket';return;endif
  call isotope_pressure(cat,i,pmax,left,status);if(status/='ok')return
  call isotope_pressure(cat,i+1,pmax,right,status);if(status/='ok')return
  if(left%mu_b<=right%mu_b)then;status='BPS transition lacks high-P Gibbs bracket';return;endif
  do j=1,64
   mid=(lo+hi)/2;p=exp(mid)
   call isotope_pressure(cat,i,p,left,status);if(status/='ok')return
   call isotope_pressure(cat,i+1,p,right,status);if(status/='ok')return
   delta=left%mu_b-right%mu_b
   if(delta<0)then;lo=mid;else;hi=mid;endif
   if(hi-lo<2e-13_wp)exit
  enddo
  p=exp((lo+hi)/2)
  call isotope_pressure(cat,i,p,left,status);if(status/='ok')return
  call isotope_pressure(cat,i+1,p,right,status);if(status/='ok')return
  call bps_state_at_pressure(cat,p,global,status);if(status/='ok')return
  if(abs(left%mu_b-right%mu_b)>1e-10_wp)then
   status='BPS transition chemical potentials disagree';return
  endif
  if(global%mu_b<min(left%mu_b,right%mu_b)-1e-10_wp)then
   status='BPS Table 2 phase preempted by another supplied isotope';return
  endif
  if(right%nb<=left%nb.or.right%epsilon<=left%epsilon)then
   status='BPS transition has nonpositive density or energy jump';return
  endif
  status='ok'
 end subroutine

 subroutine generate_bps_primary_native(source_path,n_points,eos,output_dir,status)
  character(*),intent(in)::source_path,output_dir
  integer,intent(in)::n_points
  type(eos_data_t),intent(out)::eos
  character(*),intent(out)::status
  type(bps_catalogue_t)::cat
  type(bps_state_t),allocatable::rows(:)
  type(bps_state_t)::s,tmp,left(NC-1),right(NC-1),drip
  real(wp)::p0,pdrip,lo,hi,p,bounds(0:NC),trans(NC-1),maxpr,maxeuler,maxjump,fdres,h
  type(bps_state_t)::s1,s2,sm1,sm2
  integer::i,j,nr,nkeep,u,ios,closeios
  logical::seen(NC)
  status='BPS requires 16 <= n_points <= 100000'
  if(n_points<16.or.n_points>100000)return
  call load_bps_primary_catalogue(source_path,cat,status);if(status/='ok')return
  call isotope_evaluate(NB_MIN,cat%z(1),cat%a(1),cat%binding(1),s)
  p0=s%pressure
  if(p0<=0)then;status='BPS lower boundary has nonpositive pressure';return;endif
  lo=log(1e-5_wp);hi=log(1e-2_wp)
  call bps_state_at_pressure(cat,exp(lo),s,status);if(status/='ok')return
  if(s%mu_b>=BPS_MN)then;status='BPS neutron drip lacks low-P bracket';return;endif
  call bps_state_at_pressure(cat,1e-2_wp,s,status);if(status/='ok')return
  if(s%mu_b<=BPS_MN)then;status='BPS neutron drip lacks high-P bracket';return;endif
  do i=1,64
   p=exp((lo+hi)/2)
   call bps_state_at_pressure(cat,p,s,status);if(status/='ok')return
   if(s%mu_b<BPS_MN)then;lo=log(p);else;hi=log(p);endif
   if(hi-lo<2e-13_wp)exit
  enddo
  pdrip=exp((lo+hi)/2)
  call bps_state_at_pressure(cat,pdrip,drip,status);if(status/='ok')return
  if(drip%isotope/=NC.or.abs(drip%mu_b-BPS_MN)>1e-9_wp)then
   status='BPS drip nucleus or chemical-potential condition failed';return
  endif
  bounds(0)=p0;bounds(NC)=pdrip
  do i=1,NC-1
   call transition_pressure(cat,i,p0,pdrip,trans(i),left(i),right(i),status)
   if(status/='ok')return
   bounds(i)=trans(i)
   if(bounds(i)<=bounds(i-1))then;status='BPS Table 2 phase ordering failed';return;endif
  enddo
  if(bounds(NC)<=bounds(NC-1))then;status='BPS final phase vanishes before drip';return;endif
  allocate(rows(n_points+3*NC));nr=0
  do i=1,n_points
   p=exp(log(p0)+log(pdrip/p0)*real(i-1,wp)/real(n_points-1,wp))
   if(i==1)p=p0
   if(i==n_points)p=pdrip
   nr=nr+1
   call bps_state_at_pressure(cat,p,rows(nr),status);if(status/='ok')return
  enddo
  do i=1,NC-1
   nr=nr+1;rows(nr)=left(i)
   nr=nr+1;rows(nr)=right(i)
  enddo
  ! At least one interior sample per phase; a global logarithmic grid can miss
  ! short-lived phases even though both transition limits are available.
  fdres=0
  do i=1,NC
   p=sqrt(bounds(i-1)*bounds(i))
   call bps_state_at_pressure(cat,p,s,status);if(status/='ok')return
   if(s%isotope/=i)then;status='BPS phase interior is not supplied-sequence minimum';return;endif
   nr=nr+1;rows(nr)=s
   h=1e-3_wp*s%nb
   call isotope_evaluate(s%nb+h,cat%z(i),cat%a(i),cat%binding(i),s1)
   call isotope_evaluate(s%nb+2*h,cat%z(i),cat%a(i),cat%binding(i),s2)
   call isotope_evaluate(s%nb-h,cat%z(i),cat%a(i),cat%binding(i),sm1)
   call isotope_evaluate(s%nb-2*h,cat%z(i),cat%a(i),cat%binding(i),sm2)
   ! Differentiate only kinetic+lattice energy to avoid cancellation of nuclear
   ! and electron rest energies, whose linear contributions cancel identically.
   fdres=max(fdres,abs((s%nb*(-s2%kinetic_lattice+8*s1%kinetic_lattice- &
     8*sm1%kinetic_lattice+sm2%kinetic_lattice)/(12*h)-s%kinetic_lattice)/s%pressure-1))
  enddo
  if(fdres>3e-7_wp)then;status='BPS finite-difference first-law check failed';return;endif
  ! Insertion sort is O(N^2) worst case; appended limits/interiors are O(NC*N).
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
    status='BPS distinct phases have unresolved baryon-density separation';return
   endif
  enddo
  seen=.false.;maxpr=0;maxeuler=0;maxjump=0
  do i=1,nkeep
   seen(rows(i)%isotope)=.true.
   maxpr=max(maxpr,rows(i)%pressure_residual)
   maxeuler=max(maxeuler,abs((rows(i)%epsilon+rows(i)%pressure)/rows(i)%nb-rows(i)%mu_b))
   if(.not.all(ieee_is_finite([rows(i)%epsilon,rows(i)%pressure,rows(i)%sound2])))then
    status='BPS generated nonfinite state';return
   endif
   if(rows(i)%pressure<=0.or.rows(i)%sound2<=0.or.rows(i)%sound2>=1)then
    status='BPS generated unstable or acausal branch';return
   endif
  enddo
  if(.not.all(seen))then;status='BPS generator omitted a Table 2 phase';return;endif
  if(maxeuler>1e-9_wp)then;status='BPS Euler identity failed';return;endif
  do i=1,NC-1
   maxjump=max(maxjump,abs(left(i)%mu_b-right(i)%mu_b))
  enddo
  call allocate_eos(eos,nkeep);eos%mn=BPS_MN;eos%mp=BPS_MP
  eos%composition_available=.true.;eos%composition_file_present=.true.;eos%beta_equilibrated=.true.
  do i=1,nkeep
   s=rows(i)
   eos%nb(i)=s%nb;eos%q(1,i)=s%pressure/s%nb
   eos%q(3,i)=s%mu_b/eos%mn-1;eos%q(4,i)=-s%mu_e/eos%mn
   eos%q(6,i)=s%epsilon/(s%nb*eos%mn)-1;eos%q(7,i)=eos%q(6,i)
   eos%iphase(i)=100+s%isotope
   eos%yn(i)=0;eos%yp(i)=0;eos%ye(i)=real(s%z,wp)/s%a;eos%ymu(i)=0
   eos%yh(i)=1;eos%abar(i)=s%a;eos%zbar(i)=s%z
  enddo
  eos%row_status=2;eos%yn_present=.true.;eos%yp_present=.true.;eos%ye_present=.true.;eos%ymu_present=.true.
  eos%yh_present=.true.;eos%abar_present=.true.;eos%zbar_present=.true.
  eos%c3_applicable=.true.;eos%c4_applicable=.true.
  open(newunit=u,file=trim(output_dir)//'/bps_states.csv',status='replace',iostat=ios)
  if(ios/=0)then;status='BPS cannot write states diagnostic';return;endif
  write(u,'(a)',iostat=ios)'n_b,epsilon,pressure,mu_b,mu_e,Z,A,phase,sound2,Gamma,pressure_residual'
  do i=1,nkeep
   if(ios/=0)exit
   s=rows(i)
   write(u,'(5(es25.16e3,","),3(i0,","),2(es25.16e3,","),es25.16e3)',iostat=ios) &
    s%nb,s%epsilon,s%pressure,s%mu_b,s%mu_e,s%z,s%a,100+s%isotope,s%sound2,s%gamma,s%pressure_residual
  enddo
  close(u,iostat=closeios)
  if(ios/=0.or.closeios/=0)then;status='BPS states diagnostic write failed';return;endif
  open(newunit=u,file=trim(output_dir)//'/bps_transitions.csv',status='replace',iostat=ios)
  if(ios/=0)then;status='BPS cannot write transitions diagnostic';return;endif
  write(u,'(a)',iostat=ios)'phase_low,phase_high,pressure,n_low,n_high,epsilon_low,epsilon_high,mu_low,mu_high'
  do i=1,NC-1
   if(ios/=0)exit
   write(u,'(2(i0,","),6(es25.16e3,","),es25.16e3)',iostat=ios)i+100,i+101,trans(i),left(i)%nb, &
    right(i)%nb,left(i)%epsilon,right(i)%epsilon,left(i)%mu_b,right(i)%mu_b
  enddo
  close(u,iostat=closeios)
  if(ios/=0.or.closeios/=0)then;status='BPS transitions diagnostic write failed';return;endif
  open(newunit=u,file=trim(output_dir)//'/primary_source_diagnostics.json',status='replace',iostat=ios)
  if(ios/=0)then;status='BPS cannot write validation diagnostic';return;endif
  write(u,'(a)',iostat=ios)'{"model":"BPS1971_sequence","doi":"10.1086/151216",'
  if(ios==0)write(u,'(a)',iostat=ios)' "provenance":"primary_equations_with_Table2_restricted_candidate_set",'
  if(ios==0)write(u,'(a)',iostat=ios)' "binding_energy_printed_resolution_MeV_per_baryon":0.0001,'
  if(ios==0)write(u,'(a)',iostat=ios)' "constants_convention":"CODATA dimensional constants, printed1971 binding data",'
  if(ios==0)write(u,'(a)',iostat=ios)' "full_UCRL11980_catalogue_reproduced":false,"modern_atomic_corrections":false,'
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)' "rows":',nkeep,',"phases":12,"coexistence_pairs":11,'
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "first_law_fd_relative_max":',fdres,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "pressure_inversion_relative_max":',maxpr,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "Euler_identity_MeV_max":',maxeuler,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "transition_delta_mu_MeV_max":',maxjump,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "neutron_drip_n_fm3":',drip%nb,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "neutron_drip_rho_g_cm3":',drip%epsilon*BPS_RHO_PER_MEVFM3,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)' "neutron_drip_mu_e_MeV":',drip%mu_e,'}'
  close(u,iostat=closeios)
  if(ios/=0.or.closeios/=0)then;status='BPS validation diagnostic write failed';return;endif
  status='ok'
 end subroutine
end module compose_eos_bps_primary
