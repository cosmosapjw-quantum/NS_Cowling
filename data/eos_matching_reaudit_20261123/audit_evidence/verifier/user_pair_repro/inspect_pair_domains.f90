program inspect_pair_domains
  use compose_eos_types
  use compose_eos_io, only: read_internal_csv
  use, intrinsic::ieee_arithmetic, only:ieee_is_finite
  implicit none
  type(eos_data_t)::cr,co
  real(wp),allocatable::pc(:),fc(:),ec(:),mc(:),pk(:),fk(:),ek(:),mk(:)
  logical,allocatable::lc(:),rc(:),lk(:),rk(:)
  character(len=2048)::low,high
  character(len=512)::status
  real(wp)::n1,n2,p1,p2,f1,f2,mu1,mu2,average,lower,upper,roundoff,mu_mid
  integer::j,k,total,retained,ordered,causal,ordered_all,causal_all
  logical::order_ok,causal_ok
  call get_command_argument(1,low);call get_command_argument(2,high)
  call read_internal_csv(trim(low)//'/eos_data.csv',trim(low)//'/eos_data.nml',cr,status)
  if(status/='ok')then;print '(a)',trim(status);stop 1;endif
  call read_internal_csv(trim(high)//'/eos_data.csv',trim(high)//'/eos_data.nml',co,status)
  if(status/='ok')then;print '(a)',trim(status);stop 1;endif
  call inspect('MYN13',cr,pc,fc,ec,mc,lc,rc)
  call inspect('OMEG1',co,pk,fk,ek,mk,lk,rk)
  total=0;retained=0;ordered=0;causal=0;ordered_all=0;causal_all=0
  do j=cr%n,3,-1
    n1=cr%nb(j);p1=pc(j);f1=fc(j);mu1=mc(j)
    do k=1,co%n-2
      total=total+1
      n2=co%nb(k);p2=pk(k);f2=fk(k);mu2=mk(k)
      order_ok=n2>n1.and.p2>p1.and.mu2>mu1.and.f2>f1
      causal_ok=.false.
      if(order_ok)then
        ordered_all=ordered_all+1
        average=(p2-p1)/(mu2-mu1);mu_mid=.5_wp*(mu1+mu2)
        lower=n1*mu_mid/mu1;upper=n2*mu_mid/mu2
        roundoff=64._wp*epsilon(1._wp)*max(1._wp,abs(average),abs(lower),abs(upper))
        causal_ok=average>=lower-roundoff.and.average<=upper+roundoff
        if(causal_ok)causal_all=causal_all+1
      endif
      if(.not.lc(j).or..not.rk(k))cycle
      retained=retained+1
      if(.not.order_ok)cycle
      ordered=ordered+1
      if(causal_ok)causal=causal+1
    enddo
  enddo
  print '(a,i0)','total_endpoint_pairs=',total
  print '(a,i0)','crust_valid_prefix_endpoints_j_ge_3=',count(lc(3:))
  print '(a,i0)','core_valid_suffix_endpoints_k_le_n_minus_2=',count(rk(:co%n-2))
  print '(a,i0)','pairs_after_retained_source_gate=',retained
  print '(a,i0)','pairs_after_source_and_ordering_gate=',ordered
  print '(a,i0)','pairs_after_source_ordering_and_causal_integral_gate=',causal
  print '(a,i0)','ordering_pass_pairs_ignoring_source_gate_DIAGNOSTIC_ONLY=',ordered_all
  print '(a,i0)','causal_integral_pass_pairs_ignoring_source_gate_DIAGNOSTIC_ONLY=',causal_all
contains
  subroutine inspect(label,e,p,f,en,mu,prefix,suffix)
    character(len=*),intent(in)::label
    type(eos_data_t),intent(in)::e
    real(wp),allocatable,intent(out)::p(:),f(:),en(:),mu(:)
    logical,allocatable,intent(out)::prefix(:),suffix(:)
    logical,allocatable::good(:),edge(:)
    real(wp),allocatable::dmu(:),dp(:),de(:),cs(:)
    real(wp)::tol,roundoff
    integer::i,n,first,mu_strict_bad,mu_roundoff_bad,energy_bad,pressure_bad,causal_bad
    n=e%n;p=e%nb*e%q(1,:);f=(e%q(6,:)+1._wp)*e%mn
    en=e%nb*f;mu=f+e%q(1,:)
    dp=p(2:)-p(:n-1);de=en(2:)-en(:n-1);dmu=mu(2:)-mu(:n-1)
    allocate(prefix(n),suffix(n),good(n),edge(n-1),cs(n-1))
    tol=1.e-12_wp*maxval(abs(p))
    do i=1,n
      good(i)=all(ieee_is_finite([p(i),f(i),en(i),mu(i)])).and.min(p(i),f(i),en(i),mu(i))>0._wp
    enddo
    mu_strict_bad=0;mu_roundoff_bad=0;energy_bad=0;pressure_bad=0;causal_bad=0;cs=0._wp
    do i=1,n-1
      roundoff=64._wp*epsilon(1._wp)*max(1._wp,abs(mu(i)),abs(mu(i+1)))
      edge(i)=p(i+1)>=p(i)-tol.and.en(i+1)>en(i).and. &
        mu(i+1)>=mu(i)-roundoff.and.dp(i)<=(1._wp+1.e-6_wp)*de(i)
      if(dmu(i)<=0._wp)mu_strict_bad=mu_strict_bad+1
      if(dmu(i)<-roundoff)mu_roundoff_bad=mu_roundoff_bad+1
      if(de(i)<=0._wp)energy_bad=energy_bad+1
      if(dp(i)<-tol)pressure_bad=pressure_bad+1
      if(de(i)>0._wp)cs(i)=dp(i)/de(i)
      if(dp(i)>(1._wp+1.e-6_wp)*de(i))causal_bad=causal_bad+1
    enddo
    prefix(1)=good(1);suffix(n)=good(n)
    do i=2,n;prefix(i)=prefix(i-1).and.good(i).and.edge(i-1);enddo
    do i=n-1,1,-1;suffix(i)=suffix(i+1).and.good(i).and.edge(i);enddo
    print '(a)','source='//label
    print '(a,i0)','rows=',n
    print '(a,es25.16e3)','density_min_fm3=',e%nb(1)
    print '(a,es25.16e3)','density_max_fm3=',e%nb(n)
    print '(a,i0)','invalid_positive_finite_points=',count(.not.good)
    print '(a,i0)','strict_mu_nonincrease_edges=',mu_strict_bad
    print '(a,i0)','mu_decrease_beyond_native_roundoff_edges=',mu_roundoff_bad
    print '(a,i0)','energy_nonincrease_edges=',energy_bad
    print '(a,i0)','pressure_decrease_beyond_native_tolerance_edges=',pressure_bad
    print '(a,i0)','causal_secant_fail_edges=',causal_bad
    print '(a,es25.16e3)','min_delta_mu_MeV=',minval(dmu)
    print '(a,es25.16e3)','max_causal_secant=',maxval(cs)
    first=0
    do i=1,n-1
      if(dmu(i)<=0._wp)then;first=i;exit;endif
    enddo
    if(first>0)then
      i=first
      print '(a,i0,a,i0)','first_mu_nonincrease_rows_1based=',i,',',i+1
      print '(a,2(es25.16e3,1x))','first_mu_nonincrease_density_fm3=',e%nb(i),e%nb(i+1)
      print '(a,2(es25.16e3,1x))','first_mu_nonincrease_mu_MeV=',mu(i),mu(i+1)
      print '(a,es25.16e3)','first_mu_nonincrease_delta_MeV=',dmu(i)
    endif
    first=0
    do i=1,n-1
      if(dp(i)>(1._wp+1.e-6_wp)*de(i))then;first=i;exit;endif
    enddo
    if(first>0)then
      i=first
      print '(a,i0,a,i0)','first_causal_secant_failure_rows_1based=',i,',',i+1
      print '(a,2(es25.16e3,1x))','first_causal_secant_failure_density_fm3=',e%nb(i),e%nb(i+1)
      print '(a,es25.16e3)','first_causal_secant_failure_value=',cs(i)
    endif
    first=0
    do i=1,n
      if(.not.prefix(i))then;first=i;exit;endif
    enddo
    print '(a,i0)','first_invalid_prefix_endpoint_1based=',first
    first=0
    do i=n,1,-1
      if(.not.suffix(i))then;first=i;exit;endif
    enddo
    print '(a,i0)','highest_invalid_suffix_start_1based=',first
  end subroutine inspect
end program inspect_pair_domains
