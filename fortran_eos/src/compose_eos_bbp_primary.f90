! BBP primary-author Table 5 (BPS1971 p309) import and distinct pressure integral.
! Source DOI 10.1086/151216; underlying BBP DOI 10.1016/0375-9474(71)90281-8.
! These paths do not implement the original compressible liquid-drop functional.
module compose_eos_bbp_primary
 use compose_eos_types,only:wp,eos_data_t,allocate_eos,validate_eos
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 integer,parameter::raw_count=35,inner_count=26
 real(wp),parameter::c_si=299792458._wp,mev_fm3_to_j_m3=1.602176634e32_wp
 real(wp),parameter::rscale=1000._wp*c_si*c_si/mev_fm3_to_j_m3
 real(wp),parameter::pscale=.1_wp/mev_fm3_to_j_m3,nscale=1e-39_wp
 type::source_t
  real(wp)::rho(raw_count),p(raw_count),nb(raw_count),eps(raw_count),gamma(raw_count)
  real(wp)::er(raw_count),pr(raw_count),nr(raw_count)
  integer::z(raw_count),a(raw_count)
 end type
 public::generate_bbp_primary_native,generate_bbp_pressure_primary_native
contains
 subroutine generate_bbp_primary_native(source_path,n_points,eos,output_dir,status)
  character(*),intent(in)::source_path,output_dir
  integer,intent(in)::n_points
  type(eos_data_t),intent(out)::eos
  character(*),intent(out)::status
  call generate(source_path,n_points,.false.,eos,output_dir,status)
 end subroutine

 subroutine generate_bbp_pressure_primary_native(source_path,n_points,eos,output_dir,status)
  character(*),intent(in)::source_path,output_dir
  integer,intent(in)::n_points
  type(eos_data_t),intent(out)::eos
  character(*),intent(out)::status
  call generate(source_path,n_points,.true.,eos,output_dir,status)
 end subroutine

 subroutine generate(source_path,n_points,derived,eos,output_dir,status)
  character(*),intent(in)::source_path,output_dir
  integer,intent(in)::n_points
  logical,intent(in)::derived
  type(eos_data_t),intent(out)::eos
  character(*),intent(out)::status
  type(source_t)::s
  real(wp)::integral(inner_count),power(inner_count-1),energy,pressure,mu
  real(wp),allocatable::nodes(:),grid(:)
  integer::i,j,k,m
  logical::ok
  eos=eos_data_t()
  status='BBP n_points must be between 3 and 1000000'
  if(n_points<3.or.n_points>1000000)return
  call load_source(source_path,s,status)
  if(status/='ok')return
  call build_integrals(s,power,integral,status)
  if(status/='ok')return
  if(derived)then
   allocate(grid(n_points),nodes(n_points+inner_count))
   do i=1,n_points
    grid(i)=s%nb(1)*exp(log(s%nb(inner_count)/s%nb(1))*real(i-1,wp)/real(n_points-1,wp))
   enddo
   grid(1)=s%nb(1);grid(n_points)=s%nb(inner_count)
   i=1;j=1;m=0
   do while(i<=n_points.or.j<=inner_count)
    m=m+1
    if(i>n_points)then
     nodes(m)=s%nb(j);j=j+1
    else if(j>inner_count)then
     nodes(m)=grid(i);i=i+1
    else if(abs(grid(i)-s%nb(j))<=64*epsilon(1._wp)*max(grid(i),s%nb(j)))then
     nodes(m)=s%nb(j);i=i+1;j=j+1
    else if(grid(i)<s%nb(j))then
     nodes(m)=grid(i);i=i+1
    else
     nodes(m)=s%nb(j);j=j+1
    endif
   enddo
  else
   m=inner_count;allocate(nodes(m));nodes=s%nb(:inner_count)
  endif
  call allocate_eos(eos,m)
  ! Header masses are encoding constants only. Source epsilon already includes rest energy.
  eos%mn=939._wp;eos%mp=939._wp;eos%il=1;eos%beta_equilibrated=.true.
  ! Z/A do not identify the neutron-gas fraction. Composition masks remain unavailable.
  do k=1,m
   if(derived)then
    call curve_state(s,power,integral,nodes(k),energy,pressure,mu)
   else
    energy=s%eps(k);pressure=s%p(k);mu=(energy+pressure)/nodes(k)
   endif
   eos%nb(k)=nodes(k);eos%q(1,k)=pressure/nodes(k)
   eos%q(3,k)=mu/eos%mn-1
   eos%q(6,k)=energy/(nodes(k)*eos%mn)-1;eos%q(7,k)=eos%q(6,k)
  enddo
  call validate_eos(eos,ok,status)
  if(.not.ok)return
  call write_diagnostics(output_dir,s,power,integral,derived,n_points,m,status)
 end subroutine

 subroutine load_source(path,s,status)
  character(*),intent(in)::path
  type(source_t),intent(out)::s
  character(*),intent(out)::status
  character(2048)::line
  character(64)::tokens(6)
  integer::u,ios,rows,j,nt
  real(wp)::raw(3),half(3)
  status='cannot open BBP primary table'
  open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
  if(ios/=0)return
  rows=0
  do
   read(u,'(a)',iostat=ios)line
   if(ios<0)exit
   if(ios/=0)then
    status='BBP primary table read failed';close(u);return
   endif
   if(len_trim(line)==len(line))then
    status='BBP primary table line exceeds 2047 characters';close(u);return
   endif
   line=adjustl(line)
   if(len_trim(line)==0.or.line(1:1)=='#')cycle
   call tokenize(line,tokens,nt,status)
   if(status/='ok')then;close(u);return;endif
   rows=rows+1
   if(rows>raw_count)then
    status='BBP primary table must contain exactly 35 numeric rows';close(u);return
   endif
   do j=1,3
    call decimal_token(tokens(j),raw(j),half(j),status)
    if(status/='ok')then;close(u);return;endif
   enddo
   read(tokens(4),*,iostat=ios)s%z(rows)
   if(ios/=0.or.verify(trim(tokens(4)),'+-0123456789')/=0)then
    status='BBP Z token is not an integer';close(u);return
   endif
   read(tokens(5),*,iostat=ios)s%a(rows)
   if(ios/=0.or.verify(trim(tokens(5)),'+-0123456789')/=0)then
    status='BBP A token is not an integer';close(u);return
   endif
   if(.not.((s%z(rows)==-1.and.s%a(rows)==-1).or. &
       (s%z(rows)>0.and.s%a(rows)>=s%z(rows))))then
    status='BBP source composition token is invalid';close(u);return
   endif
   if(len_trim(tokens(6))/=4.or.index(tokens(6),'.')/=2.or. &
      verify(trim(tokens(6)),'0123456789.')/=0)then
    status='BBP Gamma token is not numeric';close(u);return
   endif
   read(tokens(6),*,iostat=ios)s%gamma(rows)
   if(ios/=0)then;status='BBP Gamma token is invalid';close(u);return;endif
   if(.not.ieee_is_finite(s%gamma(rows)))then
    status='BBP Gamma is nonfinite';close(u);return
   endif
   if(s%gamma(rows)<=0)then;status='BBP Gamma must be positive';close(u);return;endif
   s%rho(rows)=raw(1);s%eps(rows)=raw(1)*rscale
   s%p(rows)=raw(2)*pscale;s%nb(rows)=raw(3)*nscale
   s%er(rows)=half(1)*rscale;s%pr(rows)=half(2)*pscale;s%nr(rows)=half(3)*nscale
  enddo
  close(u)
  if(rows/=raw_count)then
   status='BBP primary table must contain exactly 35 numeric rows';return
  endif
  if(.not.all(ieee_is_finite(s%eps)).or..not.all(ieee_is_finite(s%p)).or. &
     .not.all(ieee_is_finite(s%nb)))then
   status='nonfinite converted BBP primary table';return
  endif
  if(any(s%eps<=s%er).or.any(s%p<=s%pr).or.any(s%nb<=s%nr))then
   status='BBP source rounding boxes must be strictly positive';return
  endif
  if(any(s%eps(2:)<=s%eps(:raw_count-1)).or.any(s%p(2:)<=s%p(:raw_count-1)).or. &
     any(s%nb(2:)<=s%nb(:raw_count-1)))then
   status='BBP source rho, pressure and density must be strictly increasing';return
  endif
  if(any(s%nb(2:)-s%nr(2:)<=s%nb(:raw_count-1)+s%nr(:raw_count-1)))then
   status='BBP density rounding intervals must be disjoint';return
  endif
  if(count(s%nb<=.08_wp)/=inner_count)then
   status='BBP declared n <= 0.08 fm^-3 domain must select exactly 26 rows';return
  endif
  status='ok'
 end subroutine

 subroutine tokenize(line,tokens,nt,status)
  character(*),intent(in)::line
  character(64),intent(out)::tokens(6)
  integer,intent(out)::nt
  character(*),intent(out)::status
  integer::i,j,last
  nt=0;tokens='';i=1;last=len_trim(line)
  status='BBP primary row requires exactly six whitespace-separated tokens'
  do while(i<=last)
   if(line(i:i)==' '.or.line(i:i)==achar(9))then;i=i+1;cycle;endif
   j=i
   do while(j<=last)
    if(line(j:j)==' '.or.line(j:j)==achar(9))exit
    j=j+1
   enddo
   nt=nt+1
   if(nt>6.or.j-i>64)return
   tokens(nt)=line(i:j-1);i=j
  enddo
  if(nt/=6)return
  status='ok'
 end subroutine

 subroutine decimal_token(token,value,half,status)
  character(*),intent(in)::token
  real(wp),intent(out)::value,half
  character(*),intent(out)::status
  integer::e,d,expo,ios,ndec
  status='BBP physical columns require finite decimal scientific tokens'
  value=0;half=0
  if(verify(trim(token),'+-0123456789.eEdD')/=0)return
  e=scan(token,'eEdD');d=index(token,'.')
  if(d/=2.or.e/=6.or.e>=len_trim(token))return
  if(verify(token(1:1),'123456789')/=0.or.verify(token(3:5),'0123456789')/=0)return
  if(verify(token(e+1:len_trim(token)),'+-0123456789')/=0)return
  read(token(e+1:),*,iostat=ios)expo
  if(ios/=0)return
  ndec=e-d-1
  if(expo < -200.or.expo>200)return
  read(token,*,iostat=ios)value
  if(ios/=0)return
  if(.not.ieee_is_finite(value))return
  half=.5_wp*10._wp**(expo-ndec)
  if(value<=half)return
  status='ok'
 end subroutine

 pure real(wp) function exprel(x)result(y)
  real(wp),intent(in)::x
  if(abs(x)<1e-4_wp)then
   y=1+x*(.5_wp+x*(1._wp/6+x*(1._wp/24+x*(1._wp/120+x/720))))
  else
   y=(exp(x)-1)/x
  endif
 end function

 subroutine build_integrals(s,power,integral,status)
  type(source_t),intent(in)::s
  real(wp),intent(out)::power(inner_count-1),integral(inner_count)
  character(*),intent(out)::status
  real(wp)::h
  integer::i
  integral=0
  do i=1,inner_count-1
   h=log(s%nb(i+1)/s%nb(i));power(i)=log(s%p(i+1)/s%p(i))/h
   integral(i+1)=integral(i)+(s%p(i)/s%nb(i))*h*exprel((power(i)-1)*h)
  enddo
  status='BBP pressure integral contains a nonfinite or nonpositive slope'
  if(.not.all(ieee_is_finite(power)).or..not.all(ieee_is_finite(integral)))return
  if(any(power<=0).or.any(integral(2:)<=integral(:inner_count-1)))return
  status='ok'
 end subroutine

 pure subroutine curve_state(s,power,integral,n,energy,pressure,mu)
  type(source_t),intent(in)::s
  real(wp),intent(in)::power(inner_count-1),integral(inner_count),n
  real(wp),intent(out)::energy,pressure,mu
  real(wp)::h,e
  integer::i
  i=1
  do while(i<inner_count-1)
   if(n<s%nb(i+1))exit
   i=i+1
  enddo
  h=log(n/s%nb(i))
  pressure=s%p(i)*exp(power(i)*h)
  e=s%eps(1)/s%nb(1)+integral(i)+(s%p(i)/s%nb(i))*h*exprel((power(i)-1)*h)
  energy=n*e
  if(n==s%nb(1))energy=s%eps(1)
  if(n==s%nb(inner_count))then
   pressure=s%p(inner_count);energy=n*(s%eps(1)/s%nb(1)+integral(inner_count))
  endif
  mu=(energy+pressure)/n
 end subroutine

 subroutine write_diagnostics(output,s,power,integral,derived,requested,generated,status)
  character(*),intent(in)::output
  type(source_t),intent(in)::s
  real(wp),intent(in)::power(inner_count-1),integral(inner_count)
  logical,intent(in)::derived
  integer,intent(in)::requested,generated
  character(*),intent(out)::status
  real(wp)::demin(inner_count-1),demax(inner_count-1),imin(inner_count-1),imax(inner_count-1)
  real(wp)::mu(inner_count),energy,pressure,mug,delta,maxround,maxrel,csmin,csmax,cs
  real(wp)::anchorlow,anchorhigh
  integer::i,j,u,ios,violations,first
  logical::bad(inner_count-1)
  mu=(s%eps(:inner_count)+s%p(:inner_count))/s%nb(:inner_count)
  do i=1,inner_count-1
   j=i+1
   demin(i)=(s%eps(j)-s%er(j))/(s%nb(j)+s%nr(j))-(s%eps(i)+s%er(i))/(s%nb(i)-s%nr(i))
   demax(i)=(s%eps(j)+s%er(j))/(s%nb(j)-s%nr(j))-(s%eps(i)-s%er(i))/(s%nb(i)+s%nr(i))
   imin(i)=(s%p(i)-s%pr(i))*(1/(s%nb(i)+s%nr(i))-1/(s%nb(j)-s%nr(j)))
   imax(i)=(s%p(j)+s%pr(j))*(1/(s%nb(i)-s%nr(i))-1/(s%nb(j)+s%nr(j)))
  enddo
  bad=demin>imax.or.demax<imin;violations=count(bad);first=0
  do i=1,inner_count-1
   if(bad(i))then;first=i;exit;endif
  enddo
  anchorlow=maxval((s%eps(:inner_count)-s%er(:inner_count))/s%nb(:inner_count)-integral)
  anchorhigh=minval((s%eps(:inner_count)+s%er(:inner_count))/s%nb(:inner_count)-integral)
  maxround=0;maxrel=0;csmin=huge(1._wp);csmax=0
  do i=1,inner_count
   call curve_state(s,power,integral,s%nb(i),energy,pressure,mug)
   maxround=max(maxround,abs(energy-s%eps(i))/s%er(i))
   maxrel=max(maxrel,abs(energy-s%eps(i))/s%eps(i))
  enddo
  ! On a polytropic segment r=P/(epsilon+P) obeys r'=r*(gamma-1-gamma*r).
  ! It is monotone or constant, so endpoint values bound c_s^2/c^2=gamma*r.
  do i=1,inner_count-1
   do j=i,i+1
    call curve_state(s,power,integral,s%nb(j),energy,pressure,mug)
    cs=power(i)*pressure/(energy+pressure);csmin=min(csmin,cs);csmax=max(csmax,cs)
   enddo
  enddo
  status='cannot create BBP primary source-values evidence'
  open(newunit=u,file=trim(output)//'/primary_source_values.csv',status='new',action='write',iostat=ios)
  if(ios/=0)return
  write(u,'(a)',err=900)'row,rho_g_cm3,n_b_fm3,epsilon_mev_fm3,pressure_mev_fm3,'// &
   'Z,A,Gamma,half_epsilon_rounding,half_pressure_rounding,half_n_rounding,in_inner_domain'
  do i=1,raw_count
   j=0;if(i<=inner_count)j=1
   write(u,'(i0,4(",",es25.16e3),2(",",i0),4(",",es25.16e3),",",i0)',err=900) &
    i,s%rho(i),s%nb(i),s%eps(i),s%p(i),s%z(i),s%a(i),s%gamma(i),s%er(i),s%pr(i),s%nr(i),j
  enddo
  close(u,iostat=ios);if(ios/=0)return
  status='cannot create BBP primary anchor-comparison evidence'
  open(newunit=u,file=trim(output)//'/primary_source_anchor_comparison.csv',status='new',action='write',iostat=ios)
  if(ios/=0)return
  write(u,'(a)',err=900)'row,n_b_fm3,raw_epsilon,pressure_integral_epsilon,pressure,'// &
   'raw_mu,pressure_integral_mu,epsilon_delta,epsilon_half_rounding,delta_over_half_rounding'
  do i=1,inner_count
   call curve_state(s,power,integral,s%nb(i),energy,pressure,mug)
   delta=energy-s%eps(i)
   write(u,'(i0,9(",",es25.16e3))',err=900) &
    i,s%nb(i),s%eps(i),energy,pressure,mu(i),mug,delta,s%er(i),delta/s%er(i)
  enddo
  close(u,iostat=ios);if(ios/=0)return
  status='cannot create BBP primary source diagnostics'
  open(newunit=u,file=trim(output)//'/primary_source_diagnostics.json',status='new',action='write',iostat=ios)
  if(ios/=0)return
  write(u,'(a)',err=900)'{"schema":"ns_cowling.bbp_primary.v1",'
  write(u,'(a)',err=900)'"source_doi":"10.1086/151216","source_table":"Table5_continued_p309",'
  write(u,'(a)',err=900)'"underlying_model_doi":"10.1016/0375-9474(71)90281-8",'
  write(u,'(a)',err=900)'"full_original_cldm_implemented":false,'
  if(derived)then
   write(u,'(a)',err=900)'"model":"BBP1971_pressure_integral",'
   write(u,'(a)',err=900)'"representation":"derived_pressure_primary_not_source_rounding_reconstruction",'
   write(u,'(a)',err=900)'"epsilon_method":"single_first_row_anchor_plus_analytic_integral_of_piecewise_polytropic_P(n)",'
   write(u,'(a)',err=900)'"source_epsilon_preserved":false,"source_n_P_knots_preserved":true,'
   write(u,'(a)',err=900)'"output_grid":"requested_log_nodes_union_all_original_inner_knots",'
   write(u,'(a)',err=900)'"firstlaw_continuous_by_construction":true,'
  else
   write(u,'(a)',err=900)'"model":"BBP1971_table5",'
   write(u,'(a)',err=900)'"representation":"primary_author_printed_table_raw_projection",'
   write(u,'(a)',err=900)'"epsilon_method":"supplied_rho_times_exact_SI_c_squared",'
   write(u,'(a)',err=900)'"source_epsilon_preserved":true,"source_n_P_knots_preserved":true,'
   write(u,'(a)',err=900)'"output_grid":"26_original_knots_no_resampling",'
   write(u,'(a)',err=900)'"firstlaw_continuous_by_construction":false,'
  endif
  write(u,'(a)',err=900)'"raw_source_rows":35,"selected_source_rows":26,"excluded_above_declared_domain":9,'
  write(u,'(a)',err=900)'"density_selection_max_fm3":0.08,"extrapolation":false,'
  write(u,'(a,i0,a,i0,a)',err=900)'"requested_n_points":',requested,',"generated_rows":',generated,','
  write(u,'(a)',err=900)'"rounding_assumption":"independent_closed_half_last_printed_decimal_digit_boxes",'
  write(u,'(a)',err=900)'"rounding_is_diagnostic_not_author_claim":true,'
  write(u,'(a,i0,a)',err=900)'"raw_nonincreasing_mu_intervals":',count(mu(2:)<=mu(:inner_count-1)),','
  write(u,'(a,i0,a)',err=900)'"rounding_incompatible_firstlaw_intervals":',violations,','
  write(u,'(a,i0,a)',err=900)'"first_incompatible_interval_one_based":',first,','
  if(violations>0)then
   write(u,'(a)',err=900)'"stable_firstlaw_within_printed_rounding":"REFUTED_FOR_SELECTED_SOURCE",'
  else
   write(u,'(a)',err=900)'"stable_firstlaw_within_printed_rounding":"NOT_REFUTED_BY_NECESSARY_BOUND_ONLY",'
  endif
  write(u,'(a,es25.16e3,a,es25.16e3,a)',err=900) &
   '"first_interval_delta_e_bounds_mev":[',demin(1),',',demax(1),'],'
  write(u,'(a,es25.16e3,a,es25.16e3,a)',err=900) &
   '"first_interval_pressure_integral_bounds_mev":[',imin(1),',',imax(1),'],'
  write(u,'(a,es25.16e3,a,es25.16e3,a)',err=900) &
   '"fixed_n_polytropic_rho_rounding_anchor_interval_mev":[',anchorlow,',',anchorhigh,'],'
  write(u,'(a,es25.16e3,a)',err=900)'"derived_max_abs_epsilon_residual_over_half_rounding":',maxround,','
  write(u,'(a,es25.16e3,a)',err=900)'"derived_max_relative_epsilon_residual":',maxrel,','
  write(u,'(a,es25.16e3,a,es25.16e3,a)',err=900) &
   '"derived_continuous_sound_speed_squared_over_c_squared_bounds":[',csmin,',',csmax,'],'
  write(u,'(a)',err=900)'"pressure_integral_anchor":"first_printed_rho_and_n_row_only",'
  write(u,'(a)',err=900)'"source_Gamma_role":"retained_as_source_diagnostic_not_interpolation_constraint",'
  write(u,'(a)',err=900)'"mass_header_mev":939.0,"mass_header_is_encoding_only":true,'
  write(u,'(a)',err=900)'"composition_available":false,"global_microphysical_reproduction_claim":false,'
  write(u,'(a)',err=900)'"physical_match_acceptance_asserted":false}'
  close(u,iostat=ios);if(ios/=0)return
  status='ok';return
900 continue
  close(u,iostat=ios)
  status='BBP diagnostic write failed'
 end subroutine
end module
