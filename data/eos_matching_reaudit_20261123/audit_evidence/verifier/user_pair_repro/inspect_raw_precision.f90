program inspect_raw_precision
  use,intrinsic::iso_fortran_env,only:real64
  implicit none
  integer,parameter::wp=real64,nselected=182
  character(len=2048)::path,line
  character(len=128)::tokens(4)
  real(wp)::v(4),n(nselected),en(nselected),p(nselected)
  real(wp)::dn(nselected),de(nselected),dp(nselected),mu(nselected),lo(nselected),hi(nselected)
  real(wp)::scale,running_lower
  integer::u,ios,i,k,nread,ndown,nrobust,nmonotonebad,first_robust,first_monotone,worst
  call get_command_argument(1,path)
  open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
  if(ios/=0)stop 1
  nread=0
  do
    read(u,'(a)',iostat=ios)line
    if(ios/=0)exit
    read(line,*,iostat=ios)v
    if(ios/=0)cycle
    read(line,*,iostat=ios)tokens
    if(ios/=0)stop 2
    if(v(2)>0.06757_wp)exit
    nread=nread+1
    if(nread>nselected)stop 3
    n(nread)=v(2);en(nread)=v(3);p(nread)=v(4)
    dn(nread)=last_place(tokens(2));de(nread)=last_place(tokens(3));dp(nread)=last_place(tokens(4))
  enddo
  close(u)
  if(nread/=nselected)stop 4
  mu=(en+p)/n;ndown=count(mu(2:)<mu(:nselected-1));worst=minloc(mu(2:)-mu(:nselected-1),1)
  print '(a,i0)','selected_raw_rows=',nread
  print '(a,i0)','raw_mu_decreasing_edges=',ndown
  print '(a,i0,a,i0)','worst_raw_mu_drop_rows_1based=',worst,',',worst+1
  print '(a,2(es25.16e3,1x))','worst_raw_mu_drop_density_fm3=',n(worst),n(worst+1)
  print '(a,es25.16e3)','worst_raw_mu_drop_MeV=',mu(worst+1)-mu(worst)
  do k=1,2
    scale=.5_wp*real(k,wp)
    lo=(en-scale*de+p-scale*dp)/(n+scale*dn)
    hi=(en+scale*de+p+scale*dp)/(n-scale*dn)
    nrobust=0;first_robust=0;nmonotonebad=0;first_monotone=0
    running_lower=lo(1)
    do i=1,nselected-1
      if(hi(i+1)<lo(i))then
        nrobust=nrobust+1
        if(first_robust==0)first_robust=i
      endif
      if(hi(i+1)<running_lower)then
        nmonotonebad=nmonotonebad+1
        if(first_monotone==0)first_monotone=i+1
      endif
      running_lower=max(running_lower,lo(i+1))
    enddo
    print '(a,f4.1)','assumed_independent_last_digit_uncertainty_factor=',scale
    print '(a,i0)','decreasing_edges_not_explainable_within_printed_uncertainty=',nrobust
    print '(a,i0)','rows_violating_any_nondecreasing_mu_realization_within_uncertainty=',nmonotonebad
    print '(a,i0)','first_monotone_infeasibility_row_1based=',first_monotone
    print '(a,2(es25.16e3,1x))','first_row_mu_bounds_MeV=',lo(1),hi(1)
    print '(a,2(es25.16e3,1x))','second_row_mu_bounds_MeV=',lo(2),hi(2)
    if(first_robust>0)then
      i=first_robust
      print '(a,i0,a,i0)','first_robust_drop_rows_1based=',i,',',i+1
      print '(a,2(es25.16e3,1x))','first_robust_drop_density_fm3=',n(i),n(i+1)
      print '(a,2(es25.16e3,1x))','first_robust_drop_mu_MeV=',mu(i),mu(i+1)
      print '(a,es25.16e3)','first_robust_drop_bound_gap_MeV=',hi(i+1)-lo(i)
    endif
    print '(a,es25.16e3)','worst_drop_bound_gap_MeV=',hi(worst+1)-lo(worst)
  enddo
contains
  real(wp) function last_place(token)result(unit)
    character(len=*),intent(in)::token
    character(len=:),allocatable::s
    integer::ep,exponent,ending,dot,digits,j,stat
    s=trim(adjustl(token));ep=scan(s,'eEdD');exponent=0;ending=len(s)
    if(ep>0)then
      read(s(ep+1:),*,iostat=stat)exponent
      if(stat/=0)stop 5
      ending=ep-1
    endif
    dot=index(s(:ending),'.');digits=0
    if(dot>0)then
      do j=dot+1,ending
        if(s(j:j)>='0'.and.s(j:j)<='9')digits=digits+1
      enddo
    endif
    unit=10._wp**(exponent-digits)
  end function last_place
end program inspect_raw_precision
