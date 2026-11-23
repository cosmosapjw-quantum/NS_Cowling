! Negele (1974), IAU Symposium 53, pp.1-25, DOI 10.1017/S0074180900099885.
! Table III p.23 and the unnumbered equation p.24. This is the author's
! analytic representation of the NV calculation, not a new Hartree-Fock solve.
! Units: n fm^-3, energy MeV, epsilon/P MeV fm^-3; sound2 = c_s^2/c^2.
module compose_eos_nv_primary
 use compose_eos_types, only: wp,eos_data_t,allocate_eos
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 real(wp),parameter :: source_coeff(0:7)=[-4._wp,.28822899_wp,.59150523_wp,.090185940_wp, &
  -.11025614_wp,.029377479_wp,-.0032618465_wp,.00013543555_wp]
 real(wp),parameter :: source_nb(11)=[2.79e35_wp,4.00e35_wp,6.00e35_wp,8.79e35_wp, &
  1.59e36_wp,3.73e36_wp,5.77e36_wp,8.91e36_wp,2.04e37_wp,4.75e37_wp,7.89e37_wp]
 real(wp),parameter :: source_table(9,11)=reshape([ &
  2.79e35_wp,140._wp,40._wp,.2_wp,-26.8_wp,4e-5_wp,.53_wp,-1.425_wp,.436_wp, &
  4e35_wp,160._wp,40._wp,.3_wp,-29.4_wp,9.7e-5_wp,.53_wp,-.962_wp,.543_wp, &
  6e35_wp,210._wp,40._wp,.6_wp,-29.5_wp,2.6e-4_wp,.53_wp,-.462_wp,.692_wp, &
  8.79e35_wp,280._wp,40._wp,1._wp,-28.5_wp,4.8e-4_wp,.53_wp,-.050_wp,.865_wp, &
  1.59e36_wp,460._wp,40._wp,1.4_wp,-29.4_wp,.0012_wp,.52_wp,.541_wp,1.214_wp, &
  3.73e36_wp,900._wp,50._wp,2.6_wp,-33.6_wp,.003_wp,.46_wp,1.465_wp,1.926_wp, &
  5.77e36_wp,1050._wp,50._wp,3.3_wp,-34.5_wp,.0047_wp,.45_wp,1.966_wp,2.408_wp, &
  8.91e36_wp,1300._wp,50._wp,4.2_wp,-35.8_wp,.0078_wp,.44_wp,2.610_wp,2.981_wp, &
  2.04e37_wp,1750._wp,50._wp,6.5_wp,-43.6_wp,.0184_wp,.35_wp,4.097_wp,4.422_wp, &
  4.75e37_wp,1460._wp,40._wp,10.9_wp,-54._wp,.0436_wp,.28_wp,6.428_wp,6.660_wp, &
  7.89e37_wp,950._wp,32._wp,15._wp,-68.3_wp,.0737_wp,.16_wp,8.611_wp,8.657_wp],[9,11])
 real(wp),parameter :: domain_lo=2.79e-4_wp,domain_hi=7.89e-2_wp,log_reference=1e-4_wp
 real(wp),parameter :: canonical_mn=939.56542052_wp,canonical_mp=938.27208816_wp
 type,public :: nv_primary_fit_t
  real(wp) :: coeff(0:7)=0,mn=0,mp=0,nlo=0,nhi=0,nref=0,table(9,11)=0
  logical :: valid=.false.
 end type
 type,public :: nv_primary_state_t
  real(wp) :: nb=0,e_offset=0,epsilon=0,pressure=0,mu_b=0,sound2=0,gamma=0
  real(wp) :: x=0,poly=0,poly1=0,poly2=0,cancellation=0
 end type
 public :: generate_nv_primary_native,load_nv_primary_source,nv_primary_state
contains
 subroutine next_tokens(u,tokens,nt,eof,status)
  integer,intent(in)::u
  character(*),intent(out)::tokens(:)
  integer,intent(out)::nt
  logical,intent(out)::eof
  character(*),intent(out)::status
  character(4096)::line
  integer::ios,i,j,ll
  tokens='';nt=0;eof=.false.;status='ok'
  do
   read(u,'(a)',iostat=ios)line
   if(ios<0)then;eof=.true.;return;endif
   if(ios/=0)then;status='NV source read error';return;endif
   if(len_trim(line)==len(line))then;status='NV source line exceeds supported length';return;endif
   line=adjustl(line);ll=len_trim(line)
   if(ll==0)cycle
   if(line(1:1)=='#'.or.line(1:1)=='!')cycle
   exit
  enddo
  i=1
  do while(i<=ll)
   if(index(' '//achar(9)//achar(13),line(i:i))/=0)then;i=i+1;cycle;endif
   j=i
   do while(j<=ll)
    if(index(' '//achar(9)//achar(13),line(j:j))/=0)exit
    j=j+1
   enddo
   nt=nt+1
   if(nt>size(tokens))then;status='NV source has extra columns';return;endif
   if(j-i>len(tokens(1)))then;status='NV source token exceeds supported length';return;endif
   tokens(nt)=line(i:j-1);i=j
  enddo
 end subroutine

 subroutine numeric_token(token,value,status)
  character(*),intent(in)::token
  real(wp),intent(out)::value
  character(*),intent(out)::status
  integer::k,ios,iexp,exponent_value,exponent_length
  value=0;status='NV source invalid numeric token'
  if(len_trim(token)==0)return
  do k=1,len_trim(token)
   if(index('0123456789+-.eEdD',token(k:k))==0)return
  enddo
  ! Source facts need exponents only -6..37. Bound the lexical exponent
  ! before conversion so hostile overflow cannot trigger the caller's FPE trap.
  ! A <=128-digit mantissa with |exponent|<=100 remains representable in real64.
  iexp=scan(token,'eEdD')
  if(iexp>0)then
   exponent_length=len_trim(token)-iexp
   if(exponent_length<1.or.exponent_length>6)return
   read(token(iexp+1:),*,iostat=ios)exponent_value
   if(ios/=0)return
   if(abs(exponent_value)>100)return
  endif
  read(token,*,iostat=ios)value
  if(ios/=0)return
  if(.not.ieee_is_finite(value))return
  status='ok'
 end subroutine

 subroutine load_nv_primary_source(path,fit,status)
  character(*),intent(in)::path
  type(nv_primary_fit_t),intent(out)::fit
  character(*),intent(out)::status
  character(128)::tokens(16)
  real(wp)::v(9)
  logical::eof
  integer::u,ios,nt,row,j
  fit=nv_primary_fit_t();status='NV source cannot be opened'
  open(newunit=u,file=path,status='old',action='read',iostat=ios)
  if(ios/=0)return
  do row=1,17
   call next_tokens(u,tokens,nt,eof,status)
   if(status/='ok')exit
   if(eof)then;status='NV source is truncated';exit;endif
   select case(row)
   case(1)
    if(nt/=1.or.tokens(1)/='NV1974_AUTHOR_FIT_V1')then;status='NV source schema mismatch';exit;endif
   case(2)
    if(nt/=3.or.tokens(1)/='MASS_MEV')then;status='NV source mass row requires two values';exit;endif
    do j=1,2
     call numeric_token(tokens(j+1),v(j),status)
     if(status/='ok')exit
    enddo
    if(status/='ok')exit
    fit%mn=v(1);fit%mp=v(2)
   case(3)
    if(nt/=3.or.tokens(1)/='DENSITY_DOMAIN_FM3')then;status='NV source domain row requires two values';exit;endif
    do j=1,2
     call numeric_token(tokens(j+1),v(j),status)
     if(status/='ok')exit
    enddo
    if(status/='ok')exit
    fit%nlo=v(1);fit%nhi=v(2)
   case(4)
    if(nt/=2.or.tokens(1)/='LOG_REFERENCE_FM3')then;status='NV source log reference row malformed';exit;endif
    call numeric_token(tokens(2),fit%nref,status)
    if(status/='ok')exit
   case(5)
    if(nt/=9.or.tokens(1)/='COEFFICIENTS')then;status='NV source coefficients require eight values';exit;endif
    do j=0,7
     call numeric_token(tokens(j+2),fit%coeff(j),status)
     if(status/='ok')exit
    enddo
    if(status/='ok')exit
   case(6)
    if(nt/=1.or.tokens(1)/='TABLE_II_9COL_11')then;status='NV source Table II marker mismatch';exit;endif
   case default
    if(nt/=9)then;status='NV source Table II row must have nine columns';exit;endif
    do j=1,9
     call numeric_token(tokens(j),fit%table(j,row-6),status)
     if(status/='ok')exit
    enddo
    if(status/='ok')exit
   end select
  enddo
  if(status=='ok')then
   call next_tokens(u,tokens,nt,eof,status)
   if(status=='ok'.and..not.eof)status='NV source has trailing data'
  endif
  close(u,iostat=ios)
  if(status/='ok')return
  if(ios/=0)then;status='NV source close failed';return;endif
  ! Model identity: this entry point implements the published coefficients only.
  ! Decimal facts parse to the same binary values as these compile-time literals.
  if(any(fit%coeff/=source_coeff))then;status='NV coefficients differ from primary Table III';return;endif
  if(fit%nlo/=domain_lo.or.fit%nhi/=domain_hi.or.fit%nref/=log_reference)then
   status='NV domain or logarithm convention differs from declared source contract';return
  endif
  if(fit%mn/=canonical_mn.or.fit%mp/=canonical_mp)then
   status='NV mass convention differs from declared project convention';return
  endif
  if(any(fit%table(1,:)/=source_nb))then;status='NV Table II density identity mismatch';return;endif
  if(any(fit%table(2:3,:)<=0).or.any(fit%table(2:3,:)/=anint(fit%table(2:3,:))))then
   status='NV Table II nucleon counts must be positive integers';return
  endif
  if(any(fit%table(6:7,:)<=0).or.any(fit%table(7,:)>=1))then
   status='NV Table II approximate gas/composition values invalid';return
  endif
  if(fit%table(8,7)/=1.966_wp)then
   status='NV1974 Table II row 7 is 1.966; retain the distinct 1973/NUDA value separately';return
  endif
  if(any(fit%table/=source_table))then;status='NV Table II values differ from primary source';return;endif
  fit%valid=.true.;status='ok'
 end subroutine

 subroutine nv_primary_state(fit,nb,state,status)
  type(nv_primary_fit_t),intent(in)::fit
  real(wp),intent(in)::nb
  type(nv_primary_state_t),intent(out)::state
  character(*),intent(out)::status
  real(wp)::x,f,fp,fpp,absolute_sum,ex,first,second,dpdn
  integer::k
  state=nv_primary_state_t();status='NV state requires a validated primary source'
  if(.not.fit%valid)return
  status='NV state density must be finite and inside the source-table domain'
  if(.not.ieee_is_finite(nb))return
  if(nb<fit%nlo.or.nb>fit%nhi)return
  x=log(nb/fit%nref)
  ! Coupled Horner evaluation computes Q,Q',Q'' without differencing or powers.
  f=fit%coeff(7);fp=0;fpp=0;absolute_sum=abs(f)
  do k=6,1,-1
   fpp=fpp*x+2*fp;fp=fp*x+f;f=f*x+fit%coeff(k)
   absolute_sum=absolute_sum*abs(x)+abs(fit%coeff(k))
  enddo
  ex=exp(f);first=ex*fp;second=ex*(fp*fp+fpp)
  state%nb=nb;state%e_offset=fit%coeff(0)+ex
  state%epsilon=nb*(fit%mn+state%e_offset)
  state%pressure=nb*first
  state%mu_b=fit%mn+state%e_offset+first
  dpdn=first+second
  status='NV analytic state is nonfinite or not mechanically stable'
  if(.not.all(ieee_is_finite([state%e_offset,state%epsilon,state%pressure,state%mu_b,dpdn])))return
  if(state%epsilon<=0.or.state%pressure<=0.or.state%mu_b<=0.or.dpdn<=0)return
  state%sound2=dpdn/state%mu_b;state%gamma=dpdn/first
  if(state%sound2>1)then;status='NV analytic state violates causal sound-speed bound';return;endif
  state%x=x;state%poly=f;state%poly1=fp;state%poly2=fpp
  state%cancellation=absolute_sum/max(abs(f),tiny(1._wp))
  status='ok'
 end subroutine

 subroutine generate_nv_primary_native(source_path,n_points,eos,output_dir,status)
  character(*),intent(in)::source_path,output_dir
  integer,intent(in)::n_points
  type(eos_data_t),intent(out)::eos
  character(*),intent(out)::status
  type(nv_primary_fit_t)::fit
  type(nv_primary_state_t)::s
  type(nv_primary_state_t),allocatable::states(:)
  character(64),parameter::files(3)=[character(64)::'nv_primary_states.csv', &
   'nv_primary_source_comparison.csv','NV_PRIMARY_GENERATION.json']
  real(wp)::nb,p_euler,p_round,density_round,e_nuda,max_e,max_nuda,max_mu,max_p_relative,max_cond
  real(wp)::mu_fit_offset,euler,elerr,min_mu_offset
  integer::i,u,ios,close_ios
  logical::exists
  eos=eos_data_t();status='NV grid requires 3 to 100000 points'
  if(n_points<3.or.n_points>100000)return
  call load_nv_primary_source(source_path,fit,status)
  if(status/='ok')return
  do i=1,3
   inquire(file=trim(output_dir)//'/'//trim(files(i)),exist=exists)
   if(exists)then;status='NV diagnostic output already exists';return;endif
  enddo
  allocate(states(n_points))
  do i=1,n_points
   nb=exp(log(fit%nlo)+real(i-1,wp)/real(n_points-1,wp)*log(fit%nhi/fit%nlo))
   if(i==1)nb=fit%nlo
   if(i==n_points)nb=fit%nhi
   call nv_primary_state(fit,nb,states(i),status)
   if(status/='ok')return
  enddo
  call allocate_eos(eos,n_points)
  eos%mn=fit%mn;eos%mp=fit%mp;eos%il=1;eos%beta_equilibrated=.true.
  ! The fit is a total one-charge barotrope. Do not invent interpolated species.
  eos%composition_available=.false.;eos%composition_file_present=.false.
  do i=1,n_points
   eos%nb(i)=states(i)%nb
   eos%q(1,i)=states(i)%pressure/states(i)%nb
   eos%q(3,i)=(states(i)%e_offset+eos%q(1,i))/fit%mn
   eos%q(6,i)=states(i)%e_offset/fit%mn;eos%q(7,i)=eos%q(6,i)
  enddo
  open(newunit=u,file=trim(output_dir)//'/'//trim(files(1)),status='new',action='write',iostat=ios)
  if(ios/=0)then;status='NV cannot create state diagnostics';return;endif
  write(u,'(a)',iostat=ios)'n_b_fm3,e_offset_mev,energy_mev_fm3,pressure_mev_fm3,mu_b_mev,cs2_over_c2,'// &
   'adiabatic_index,x,log_polynomial,log_polynomial_prime,log_polynomial_second,polynomial_cancellation_ratio'
  do i=1,n_points
   if(ios/=0)exit
   s=states(i)
   write(u,'(*(es25.16e3,:,","))',iostat=ios)s%nb,s%e_offset,s%epsilon,s%pressure,s%mu_b, &
    s%sound2,s%gamma,s%x,s%poly,s%poly1,s%poly2,s%cancellation
  enddo
  close(u,iostat=close_ios)
  if(ios==0)ios=close_ios
  if(ios/=0)then;status='NV state diagnostics write failed';return;endif
  open(newunit=u,file=trim(output_dir)//'/'//trim(files(2)),status='new',action='write',iostat=ios)
  if(ios/=0)then;status='NV cannot create source comparison';return;endif
  write(u,'(a)',iostat=ios)'n_b_fm3,e_table_1974_mev,e_historical_nuda1973_mev,e_fit_mev,'// &
   'e_fit_minus_1974_mev,e_fit_minus_nuda1973_mev,mu_n_table_kinetic_mev,mu_fit_offset_mev,'// &
   'mu_fit_minus_table_mev,pressure_fit_mev_fm3,pressure_euler_table_mev_fm3,'// &
   'pressure_euler_rounding_only_mev_fm3,pressure_relative_difference'
  max_e=0;max_nuda=0;max_mu=0;max_p_relative=0;max_cond=maxval(states%cancellation)
  do i=1,11
   if(ios/=0)exit
   nb=fit%table(1,i)*1e-39_wp
   if(i==1)nb=fit%nlo
   if(i==11)nb=fit%nhi
   call nv_primary_state(fit,nb,s,status)
   if(status/='ok')then;close(u,iostat=close_ios);return;endif
   e_nuda=fit%table(8,i)
   if(i==7)e_nuda=1.996_wp
   mu_fit_offset=s%e_offset+s%pressure/nb
   p_euler=nb*(fit%table(4,i)-fit%table(8,i))
   ! Rounding only: half-last digits for 0.1 MeV mu, 0.001 MeV E,
   ! and three-significant-digit n. This excludes HF/fit systematic errors.
   density_round=.005_wp*10._wp**floor(log10(fit%table(1,i)))*1e-39_wp
   p_round=nb*.0505_wp+density_round*(abs(fit%table(4,i)-fit%table(8,i))+.0505_wp)
   elerr=(s%pressure-p_euler)/s%pressure
   write(u,'(*(es25.16e3,:,","))',iostat=ios)nb,fit%table(8,i),e_nuda,s%e_offset, &
    s%e_offset-fit%table(8,i),s%e_offset-e_nuda,fit%table(4,i),mu_fit_offset, &
    mu_fit_offset-fit%table(4,i),s%pressure,p_euler,p_round,elerr
   max_e=max(max_e,abs(s%e_offset-fit%table(8,i)))
   max_nuda=max(max_nuda,abs(s%e_offset-e_nuda))
   max_mu=max(max_mu,abs(mu_fit_offset-fit%table(4,i)))
   max_p_relative=max(max_p_relative,abs(elerr));max_cond=max(max_cond,s%cancellation)
  enddo
  close(u,iostat=close_ios)
  if(ios==0)ios=close_ios
  if(ios/=0)then;status='NV source comparison write failed';return;endif
  euler=maxval(abs((states%epsilon+states%pressure)/states%nb-states%mu_b))
  min_mu_offset=minval(states%mu_b-fit%mn)
  open(newunit=u,file=trim(output_dir)//'/'//trim(files(3)),status='new',action='write',iostat=ios)
  if(ios/=0)then;status='NV cannot create generation receipt';return;endif
  write(u,'(a)',iostat=ios)'{"schema":"ns_cowling.nv_primary_generation.v1",'
  if(ios==0)write(u,'(a)',iostat=ios)'"model":"NV1974_author_fit","status":"NATIVE_AUTHOR_ANALYTIC_FIT_GENERATED",'
  if(ios==0)write(u,'(a)',iostat=ios)'"doi":"10.1017/S0074180900099885","source_pages":[20,23,24],'
  if(ios==0)write(u,'(a)',iostat=ios)'"underlying_hf_reference_doi":"10.1016/0375-9474(73)90349-7",'
  if(ios==0)write(u,'(a)',iostat=ios)'"representation":"published_author_energy_fit_and_analytic_derivatives",'
  if(ios==0)write(u,'(a)',iostat=ios)'"full_hartree_fock_recomputed":false,"source_values_modified":false,'
  if(ios==0)write(u,'(a)',iostat=ios)'"electron_energy_included":true,"composition_available":false,"extrapolation":false,'
  if(ios==0)write(u,'(a)',iostat=ios)'"mass_convention":"project_neutron_rest_energy_explicit_not_1974_numeric_determination",'
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"neutron_rest_energy_mev":',fit%mn,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"proton_header_mass_mev":',fit%mp,','
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)'"rows":',n_points,','
  if(ios==0)write(u,'(a,es25.16e3,a,es25.16e3,a)',iostat=ios)'"density_domain_fm3":[',fit%nlo,',',fit%nhi,'],'
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"max_energy_fit_minus_1974_table_mev":',max_e,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"max_energy_fit_minus_historical_nuda1973_mev":',max_nuda,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"max_mu_fit_minus_hf_eigenvalue_mev":',max_mu,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"max_pressure_relative_euler_diagnostic":',max_p_relative,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"max_euler_roundoff_mev":',euler,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"max_polynomial_cancellation_ratio":',max_cond,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"minimum_fit_mu_minus_neutron_rest_mev":',min_mu_offset,','
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)'"sampled_fit_mu_below_neutron_rest_rows":',count(states%mu_b<fit%mn),','
  if(ios==0)write(u,'(a)',iostat=ios)'"low_density_fit_limitation":"published_fit_mu_falls_below_free_neutron_threshold_near_drip_not_repaired",'
  if(ios==0)write(u,'(a)',iostat=ios)'"table_mu_check":"rounded_HF_last_occupied_eigenvalue_not_an_exact_fit_constraint",'
  if(ios==0)write(u,'(a)',iostat=ios)'"table_discrepancy":"1974_row7_energy_1.966_vs_historical_nuda1973_1.996_retained",'
  if(ios==0)write(u,'(a)',iostat=ios)'"physical_matching_acceptance_asserted":false}'
  close(u,iostat=close_ios)
  if(ios==0)ios=close_ios
  if(ios/=0)then;status='NV generation receipt write failed';return;endif
  status='ok'
 end subroutine
end module compose_eos_nv_primary
