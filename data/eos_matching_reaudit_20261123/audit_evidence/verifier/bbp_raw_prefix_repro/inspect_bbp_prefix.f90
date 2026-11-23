program inspect_bbp_prefix
  use compose_eos_types
  use compose_eos_io, only: read_internal_csv
  use compose_eos_verify, only: verify_report_t, check_consistency
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(eos_data_t) :: eos, prefix, core
  type(verify_report_t) :: report
  character(len=2048) :: input_path, core_list_path, core_path, core_id
  character(len=512) :: status
  real(wp), allocatable :: p(:), eps(:), mu(:), mu_euler(:)
  real(wp) :: h1,h2,s12,s23,d1,d1_stable,dmu,dp,lower,upper
  real(wp) :: causal_lower,causal_upper,roundoff,mean_n
  real(wp) :: reencoded_mu(3),gauge_d1,gauge_min,gauge_max,pr1,pr2
  real(wp) :: target_q3(3),target_q4(3),target_q6(3),target_f(3),target_eps(3),target_muq(3)
  real(wp) :: gauge_dp,gauge_df,gauge_deps
  integer :: i,j,eligible,c11_fail,monotone_prefix,mu_decreases,first_decrease,integral_bad
  integer :: list_unit,ios,core_count,gauge_fail
  logical :: monotone
  call get_command_argument(1,input_path)
  call read_internal_csv(trim(input_path)//'/eos_data.csv',trim(input_path)//'/eos_data.nml',eos,status)
  if(status/='ok')then
    print '(a)',trim(status)
    stop 1
  endif
  if(eos%n<3)stop 2
  p=eos%nb*eos%q(1,:)
  eps=eos%nb*(1._wp+eos%q(7,:))*eos%mn
  mu=(1._wp+eos%q(3,:))*eos%mn
  mu_euler=(eps+p)/eos%nb
  print '(a)','scope=immutable source; exact first three retained knots; no source repair'
  print '(a,i0)','rows=',eos%n
  print '(a,i0)','Il=',eos%il
  print '(a,es26.17e3)','max_abs_muQ3_minus_muEuler=',maxval(abs(mu-mu_euler))
  print '(a)','first_three_rows: index,n,P,epsilon,mu_Q3,mu_Euler,Iphase'
  do i=1,3
    print '(i0,5(a,es26.17e3),a,i0)',i,',',eos%nb(i),',',p(i),',',eps(i),',',mu(i),',',mu_euler(i),',',eos%iphase(i)
  enddo
  h1=eos%nb(2)-eos%nb(1)
  h2=eos%nb(3)-eos%nb(2)
  s12=(mu(2)-mu(1))/h1
  s23=(mu(3)-mu(2))/h2
  d1=-(2._wp*h1+h2)/(h1*(h1+h2))*mu(1) &
       +(h1+h2)/(h1*h2)*mu(2)-h1/(h2*(h1+h2))*mu(3)
  d1_stable=((2._wp*h1+h2)*s12-h1*s23)/(h1+h2)
  print '(a,es26.17e3)','first_log_step=',log(eos%nb(2))-log(eos%nb(1))
  print '(a,es26.17e3)','second_log_step=',log(eos%nb(3))-log(eos%nb(2))
  print '(a,es26.17e3)','first_mu_secant=',s12
  print '(a,es26.17e3)','second_mu_secant=',s23
  print '(a,es26.17e3)','first_three_point_endpoint_derivative_native_form=',d1
  print '(a,es26.17e3)','first_three_point_endpoint_derivative_slope_form=',d1_stable
  print '(a,es26.17e3)','C11_threshold=',-1.e-10_wp
  print '(a,es26.17e3)','native_mask_max_log_step=',0.25_wp
  pr1=abs(p(2)-p(1))/max(abs(p(1)),abs(p(2)))
  pr2=abs(p(3)-p(2))/max(abs(p(2)),abs(p(3)))
  print '(a,es26.17e3)','first_pressure_relative_jump=',pr1
  print '(a,es26.17e3)','second_pressure_relative_jump=',pr2
  print '(a,es26.17e3)','thermodynamic_jump_pressure_threshold=',1.e-10_wp
  if(any(eos%iphase(1:3)/=eos%iphase(1)))stop 4
  if(max(log(eos%nb(2))-log(eos%nb(1)),log(eos%nb(3))-log(eos%nb(2)))>0.25_wp)stop 5
  if(min(pr1,pr2)<=1.e-10_wp)stop 6
  print '(a)','first_endpoint_stencil_cannot_be_boundary_masked=CONFIRMED'
  call get_command_argument(2,core_list_path)
  core_count=0;gauge_fail=0
  gauge_min=huge(1._wp);gauge_max=-huge(1._wp)
  if(len_trim(core_list_path)>0)then
    open(newunit=list_unit,file=trim(core_list_path),status='old',action='read',iostat=ios)
    if(ios/=0)stop 7
    print '(a)','affected_core_view_reencoding: source_view_id,mn,Il,Yq,dPdn,dFdn,dEPSdn,dMUdn'
    do
      read(list_unit,'(a)',iostat=ios)core_id
      if(ios<0)exit
      if(ios/=0)stop 8
      read(list_unit,'(a)',iostat=ios)core_path
      if(ios/=0)stop 9
      call read_internal_csv(trim(core_path)//'/eos_data.csv',trim(core_path)//'/eos_data.nml',core,status)
      if(status/='ok')stop 10
      do i=1,3
        target_q3(i)=((1._wp+eos%q(6,i))*eos%mn+eos%q(1,i))/core%mn-1._wp
        target_q6(i)=((1._wp+eos%q(6,i))*eos%mn)/core%mn-1._wp
        target_q4(i)=eos%q(4,i)*eos%mn/core%mn
      enddo
      if(core%il==0)target_q3=target_q3-core%charge_fraction_grid*target_q4
      reencoded_mu=(1._wp+target_q3)*core%mn
      target_muq=target_q4*core%mn
      if(core%il==0)reencoded_mu=reencoded_mu+core%charge_fraction_grid*target_muq
      target_f=(1._wp+target_q6)*core%mn
      target_eps=eos%nb(1:3)*target_f
      gauge_d1=first_derivative(reencoded_mu)
      gauge_dp=first_derivative(p(1:3))
      gauge_df=first_derivative(target_f)
      gauge_deps=first_derivative(target_eps)
      if(.not.all(ieee_is_finite([gauge_d1,gauge_dp,gauge_df,gauge_deps])))stop 12
      if(gauge_deps<=0._wp)stop 13
      core_count=core_count+1
      if(gauge_d1<-1.e-10_wp)gauge_fail=gauge_fail+1
      gauge_min=min(gauge_min,gauge_d1);gauge_max=max(gauge_max,gauge_d1)
      print '(a,a,es26.17e3,a,i0,5(a,es26.17e3))',trim(core_id),',',core%mn,',',core%il, &
          ',',core%charge_fraction_grid,',',gauge_dp,',',gauge_df,',',gauge_deps,',',gauge_d1
    enddo
    close(list_unit)
    print '(a,i0)','affected_core_views_checked=',core_count
    print '(a,i0)','C11_failed_affected_core_views=',gauge_fail
    print '(a,es26.17e3)','core_reencoding_min_endpoint_derivative=',gauge_min
    print '(a,es26.17e3)','core_reencoding_max_endpoint_derivative=',gauge_max
    if(core_count<1.or.gauge_fail/=core_count)stop 11
  endif
  print '(a)','intervals: i,delta_mu,delta_P,n1delta_mu,n2delta_mu,causal_lower,causal_upper,mean_n'
  integral_bad=0;mu_decreases=0;first_decrease=0
  do i=1,eos%n-1
    dmu=mu_euler(i+1)-mu_euler(i)
    dp=p(i+1)-p(i)
    if(dmu<=0._wp)then
      mu_decreases=mu_decreases+1
      if(first_decrease==0)first_decrease=i
      cycle
    endif
    lower=eos%nb(i)*dmu
    upper=eos%nb(i+1)*dmu
    causal_lower=lower*(mu_euler(i+1)+mu_euler(i))/(2._wp*mu_euler(i))
    causal_upper=upper*(mu_euler(i+1)+mu_euler(i))/(2._wp*mu_euler(i+1))
    roundoff=64._wp*epsilon(1._wp)*max(1._wp,abs(dp),abs(lower),abs(upper))
    if(dp<lower-roundoff.or.dp>upper+roundoff)integral_bad=integral_bad+1
    mean_n=dp/dmu
    if(i<=3) print '(i0,7(a,es26.17e3))',i,',',dmu,',',dp,',',lower,',',upper, &
        ',',causal_lower,',',causal_upper,',',mean_n
  enddo
  print '(a,i0)','mu_nonincreasing_edges=',mu_decreases
  print '(a,i0)','first_mu_nonincreasing_edge_left_1based=',first_decrease
  print '(a,i0)','positive_dmu_intervals_failing_necessary_integral_bound=',integral_bad
  eligible=0;c11_fail=0;monotone_prefix=0
  print '(a)','prefix: retained_rows,C11,min_dmu_dnb,derivative_points,masked,invalid'
  do j=3,eos%n
    eligible=eligible+1
    call allocate_eos(prefix,j)
    prefix%mn=eos%mn;prefix%mp=eos%mp;prefix%il=eos%il
    prefix%temperature_mev=eos%temperature_mev
    prefix%charge_fraction_grid=eos%charge_fraction_grid
    prefix%beta_equilibrated=eos%beta_equilibrated
    prefix%composition_available=.false.
    prefix%nb=eos%nb(:j);prefix%q=eos%q(:,:j);prefix%iphase=eos%iphase(:j)
    call check_consistency(prefix,report)
    if(.not.report%c11)c11_fail=c11_fail+1
    monotone=.true.
    do i=1,j-1
      roundoff=64._wp*epsilon(1._wp)*max(1._wp,abs(mu_euler(i)),abs(mu_euler(i+1)))
      if(mu_euler(i+1)<mu_euler(i)-roundoff)monotone=.false.
    enddo
    if(monotone)monotone_prefix=monotone_prefix+1
    print '(i0,a,l1,a,es26.17e3,3(a,i0))',j,',',report%c11,',',report%min_dmubdnb, &
        ',',report%derivative_points,',',report%derivative_masked,',',report%invalid_derivative_points
  enddo
  print '(a,i0)','eligible_retained_prefix_lengths=',eligible
  print '(a,i0)','native_C11_failed_prefix_lengths=',c11_fail
  print '(a,i0)','monotone_mu_prefix_lengths=',monotone_prefix
  if(d1>=-1.e-10_wp.or.c11_fail/=eligible)stop 3
  print '(a)','FIXED_PREFIX_C11_CERTIFICATE=CONFIRMED_FOR_THIS_SOURCE_AND_FROZEN_NATIVE_CONTRACT'
contains
  function first_derivative(v) result(result_d)
    real(wp),intent(in) :: v(3)
    real(wp) :: result_d
    result_d=-(2._wp*h1+h2)/(h1*(h1+h2))*v(1) &
        +(h1+h2)/(h1*h2)*v(2)-h1/(h2*(h1+h2))*v(3)
  end function first_derivative
end program inspect_bbp_prefix
