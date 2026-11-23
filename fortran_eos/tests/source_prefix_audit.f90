program source_prefix_audit
  !! Finite-table necessary conditions only. No interpolation/continuous proof.
  !! Usage: source_prefix_audit INPUT.tsv OUTPUT.json
  !! Header: mn <MeV>; following rows: nb Q1 Q7 Iphase (whitespace separated).
  use, intrinsic :: iso_fortran_env, only: real64,real128,iostat_end,iostat_eor
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter::q=real128
  real(q),parameter::budget=64._q*real(epsilon(1._real64),q)
  real(q),allocatable::data(:,:),grown(:,:)
  integer,allocatable::phases(:),grown_phase(:)
  character(len=4096)::line
  character(len=1024)::input_path,output_path
  character(len=128)::tokens(4)
  character(len=160)::reason
  character(len=64)::status
  integer::unit,out,ios,n,capacity,line_number,nt,i,j,phase,pressure_bad,law_bad,first_bad,chars_read
  integer::nominal_pressure_bad,nominal_law_bad,coexistence_count
  real(q)::mass,mass_error,values(3)
  real(q)::nn(2),qq1(2),qq7(2),pp(2),ee(2),dn(2),dq1(2),dq7(2),pmin(2),pmax(2)
  real(q)::mumin(2),mumax(2),delta,delta_min,delta_max,integral_min,integral_max,ilo_nom,ihi_nom
  real(q)::a,b,qa,qb,dmin,dmax,quad_budget,quad_pressure,work_nominal
  logical::ok,header,pressure_failure,law_failure,coexistence,nominal_pf,nominal_lf
  if(command_argument_count()/=2)then
    print '(a)','usage: source_prefix_audit INPUT.tsv OUTPUT.json'
    stop 2
  endif
  call get_command_argument(1,input_path);call get_command_argument(2,output_path)
  line_number=0;n=0;header=.false.;capacity=256
  allocate(data(3,capacity),phases(capacity))
  open(newunit=unit,file=trim(input_path),status='old',action='read',iostat=ios)
  if(ios/=0)call invalid('input file open failed')
  do
    read(unit,'(a)',advance='no',size=chars_read,iostat=ios)line
    if(ios==iostat_end)exit
    line_number=line_number+1
    if(ios==0)call invalid('input line exceeds supported length')
    if(ios/=iostat_eor)call invalid('input line read failed')
    line=adjustl(line)
    if(len_trim(line)==0)cycle
    if(line(1:1)=='#')cycle
    call tokenize(line,tokens,nt,ok)
    if(.not.ok)call invalid('expected bounded whitespace-separated fields')
    if(.not.header)then
      if(nt/=2.or.trim(tokens(1))/='mn')call invalid('first data line must be: mn <MeV>')
      call parse_real(tokens(2),mass,ok)
      if(.not.ok)call invalid('mass token is malformed, nonfinite, or outside binary64 range')
      if(mass<=0)call invalid('reference mass must be positive')
      mass_error=budget*abs(mass);header=.true.;cycle
    endif
    if(nt/=4)call invalid('each row requires nb Q1 Q7 Iphase')
    do j=1,3
      call parse_real(tokens(j),values(j),ok)
      if(.not.ok)call invalid('numeric row token malformed, nonfinite, or outside binary64 range')
    enddo
    if(verify(trim(tokens(4)),'0123456789+-')/=0)call invalid('phase id must be an integer')
    read(tokens(4),*,iostat=ios)phase
    if(ios/=0)call invalid('phase id outside integer range')
    if(values(1)<=0.or.values(2)<=0.or.values(3)<=-1) &
      call invalid('nb, pressure, and total energy must be positive')
    if(n>0)then
      if(values(1)<=data(1,n))call invalid('nb must be strictly increasing')
      if(values(1)*(1-budget)<=data(1,n)*(1+budget)) &
        call invalid('neighbor density uncertainty intervals overlap; order unresolved')
    endif
    if(n==capacity)then
      if(capacity>=1048576)call invalid('more than 1048576 input rows are unsupported')
      allocate(grown(3,2*capacity),grown_phase(2*capacity))
      grown(:,:n)=data(:,:n);grown_phase(:n)=phases(:n)
      call move_alloc(grown,data);call move_alloc(grown_phase,phases);capacity=2*capacity
    endif
    n=n+1;data(:,n)=values;phases(n)=phase
  enddo
  close(unit)
  if(.not.header)call invalid('missing reference-mass header')
  if(n<2)call invalid('at least two finite rows are required')
  pressure_bad=0;law_bad=0;first_bad=0;nominal_pressure_bad=0;nominal_law_bad=0;coexistence_count=0
  do i=1,n-1
    call interval(i)
    if(pressure_failure)pressure_bad=pressure_bad+1
    if(law_failure)law_bad=law_bad+1
    if(nominal_pf)nominal_pressure_bad=nominal_pressure_bad+1
    if(nominal_lf)nominal_law_bad=nominal_law_bad+1
    if(coexistence)coexistence_count=coexistence_count+1
    if(first_bad==0.and.(pressure_failure.or.law_failure))first_bad=i
  enddo
  status='FINITE_TABLE_NECESSARY_BOUNDS_SATISFIED'
  if(first_bad>0)status='REFUTED_BY_FINITE_TABLE'
  open(newunit=out,file=trim(output_path),status='replace',action='write',iostat=ios)
  if(ios/=0)then
    print '(a)','output certificate open failed';stop 2
  endif
  write(out,'(a)')'{'
  write(out,'(a)')'"schema":"finite_cold_prefix_audit_v1",'
  write(out,'(a)')'"status":"'//trim(status)//'",'
  write(out,'(a)')'"input_valid":true,'
  write(out,'(a)')'"claim":"finite-table necessary condition only; no continuous first-law or global-stability proof",'
  write(out,'(a)')'"evaluation":"Fortran real128; P=nb*Q1; epsilon=nb*mn*(1+Q7)",'
  write(out,'(a)')'"uncertainty_policy":"64*epsilon(binary64): relative nb,mn,Q1; absolute max(1,abs(Q7))*factor for Q7 normalization",'
  write(out,'(a)')'"bound_policy":"correlated mn and nb are retained algebraically; extrema encompass every value in source boxes",'
  write(out,'(a)')'"tolerances":"no matching tolerance changed; only the explicitly reported serialization/roundoff model is used",'
  write(out,'(a,es28.18e4,a)')'"binary64_uncertainty_factor":',budget,','
  write(out,'(a,es28.18e4,a)')'"mn_MeV":',mass,','
  write(out,'(a,es28.18e4,a)')'"mn_halfwidth_MeV":',mass_error,','
  write(out,'(a,i0,a)')'"rows":',n,','
  write(out,'(a,i0,a)')'"pressure_decrease_count":',pressure_bad,','
  write(out,'(a,i0,a)')'"firstlaw_violation_count":',law_bad,','
  write(out,'(a,i0,a)')'"nominal_pressure_decrease_count":',nominal_pressure_bad,','
  write(out,'(a,i0,a)')'"nominal_firstlaw_violation_count":',nominal_law_bad,','
  write(out,'(a,i0,a)')'"coexistence_compatible_interval_count":',coexistence_count,','
  write(out,'(a,i0,a)')'"first_bad_interval_one_based":',first_bad,','
  write(out,'(a)')'"intervals":['
  do i=1,n-1
    call interval(i)
    write(out,'(a,i0,a)',advance='no')'{"i":',i,','
    call pair('nb',nn);call pair('Q1',qq1);call pair('Q7',qq7)
    write(out,'(a,i0,a,i0,a)',advance='no')'"phase":[',phases(i),',',phases(i+1),'],'
    call pair('P',pp);call pair('epsilon',ee)
    call pair('nb_halfwidth',dn);call pair('Q1_halfwidth',dq1);call pair('Q7_halfwidth',dq7)
    call pair('P1_bounds',[pmin(1),pmax(1)]);call pair('P2_bounds',[pmin(2),pmax(2)])
    call number('delta_e_over_n',delta)
    call pair('delta_e_over_n_bounds',[delta_min,delta_max])
    call pair('integral_bounds',[integral_min,integral_max])
    call pair('nominal_integral_bounds',[ilo_nom,ihi_nom])
    call number('quad_roundoff_budget_MeV',quad_budget)
    write(out,'(a)',advance='no')'"pressure_refuted":'//boolean(pressure_failure)//','
    write(out,'(a)',advance='no')'"firstlaw_refuted":'//boolean(law_failure)//','
    write(out,'(a)',advance='no')'"coexistence_compatible":'//boolean(coexistence)//'}'
    if(i<n-1)then
      write(out,'(a)')','
    else
      write(out,'(a)')''
    endif
  enddo
  write(out,'(a)',iostat=ios)']}';close(out,iostat=j)
  if(ios/=0.or.j/=0)then
    print '(a)','certificate write failed';stop 2
  endif
  print '(a)',trim(status)
contains
  subroutine tokenize(source,words,count,valid)
    character(len=*),intent(in)::source
    character(len=*),intent(out)::words(:)
    integer,intent(out)::count
    logical,intent(out)::valid
    integer::a,b,l
    words='';count=0;valid=.false.;l=len_trim(source);a=1
    do while(a<=l)
      if(source(a:a)==' '.or.source(a:a)==achar(9))then
        a=a+1;cycle
      endif
      b=a
      do while(b<=l)
        if(source(b:b)==' '.or.source(b:b)==achar(9))exit
        b=b+1
      enddo
      if(count==size(words).or.b-a>len(words(1)))return
      count=count+1;words(count)=source(a:b-1);a=b
    enddo
    valid=.true.
  end subroutine
  subroutine parse_real(token,value,valid)
    character(len=*),intent(in)::token
    real(q),intent(out)::value
    logical,intent(out)::valid
    integer::ep,power,stat
    value=0;valid=.false.
    if(len_trim(token)==0)return
    if(verify(trim(token),'0123456789.+-eEdD')/=0)return
    ep=scan(trim(token),'eEdD')
    if(ep>0)then
      if(len_trim(token)-ep>5.or.len_trim(token)==ep)return
      power=0;read(token(ep+1:),*,iostat=stat)power
      if(stat/=0.or.abs(power)>4000)return
    endif
    read(token,*,iostat=stat)value
    if(stat/=0)return
    if(.not.ieee_is_finite(value))return
    if(abs(value)>real(huge(1._real64),q))return
    if(value/=0)then
      if(abs(value)<scale(1._q,-1074)*0.49_q)return
    endif
    valid=.true.
  end subroutine
  subroutine interval(index)
    integer,intent(in)::index
    nn=data(1,index:index+1);qq1=data(2,index:index+1);qq7=data(3,index:index+1)
    dn=budget*abs(nn);dq1=budget*abs(qq1);dq7=budget*max(1._q,abs(qq7))
    pp=nn*qq1;ee=nn*mass*(1+qq7)
    pmin=(nn-dn)*(qq1-dq1);pmax=(nn+dn)*(qq1+dq1)
    quad_pressure=128*epsilon(1._q)*maxval(abs(pp))
    pmin=pmin-quad_pressure;pmax=pmax+quad_pressure
    pressure_failure=pmax(2)<pmin(1);nominal_pf=pp(2)<pp(1)
    delta=mass*(qq7(2)-qq7(1))
    dmin=qq7(2)-dq7(2)-qq7(1)-dq7(1)
    dmax=qq7(2)+dq7(2)-qq7(1)+dq7(1)
    delta_min=min((mass-mass_error)*dmin,(mass+mass_error)*dmin)
    delta_max=max((mass-mass_error)*dmax,(mass+mass_error)*dmax)
    ! P1*(1/n1-1/n2)=Q1_1*(1-n1/n2), retaining n/P correlation.
    integral_min=(qq1(1)-dq1(1))*(1-(nn(1)+dn(1))/(nn(2)-dn(2)))
    integral_max=(qq1(2)+dq1(2))*((nn(2)+dn(2))/(nn(1)-dn(1))-1)
    quad_budget=128*epsilon(1._q)*max(abs(delta_min),abs(delta_max), &
      abs(integral_min),abs(integral_max),maxval(abs(mass*(1+qq7))))
    delta_min=delta_min-quad_budget;delta_max=delta_max+quad_budget
    integral_min=integral_min-quad_budget;integral_max=integral_max+quad_budget
    law_failure=delta_max<integral_min.or.delta_min>integral_max
    ilo_nom=qq1(1)*(1-nn(1)/nn(2));ihi_nom=qq1(2)*(nn(2)/nn(1)-1)
    nominal_lf=delta<ilo_nom.or.delta>ihi_nom
    mumin=(mass-mass_error)*(1+qq7-dq7)+qq1-dq1
    mumax=(mass+mass_error)*(1+qq7+dq7)+qq1+dq1
    ! For physical rows, total specific energy is positive. Very small
    ! uncertain values can straddle zero; cover all four product corners.
    do j=1,2
      a=(mass-mass_error)*(1+qq7(j)-dq7(j))
      b=(mass+mass_error)*(1+qq7(j)-dq7(j))
      qa=(mass-mass_error)*(1+qq7(j)+dq7(j))
      qb=(mass+mass_error)*(1+qq7(j)+dq7(j))
      mumin(j)=min(a,b,qa,qb)+qq1(j)-dq1(j)-quad_budget
      mumax(j)=max(a,b,qa,qb)+qq1(j)+dq1(j)+quad_budget
    enddo
    coexistence=phases(index)/=phases(index+1).and. &
      max(pmin(1),pmin(2))<=min(pmax(1),pmax(2)).and. &
      max(mumin(1),mumin(2))<=min(mumax(1),mumax(2))
  end subroutine
  function boolean(value) result(word)
    logical,intent(in)::value
    character(len=:),allocatable::word
    if(value)then
      word='true'
    else
      word='false'
    endif
  end function
  subroutine pair(key,value)
    character(len=*),intent(in)::key
    real(q),intent(in)::value(2)
    write(out,'(a,es28.18e4,a,es28.18e4,a)',advance='no')'"'//trim(key)//'":[',value(1),',',value(2),'],'
  end subroutine
  subroutine number(key,value)
    character(len=*),intent(in)::key
    real(q),intent(in)::value
    write(out,'(a,es28.18e4,a)',advance='no')'"'//trim(key)//'":',value,','
  end subroutine
  subroutine invalid(message)
    character(len=*),intent(in)::message
    integer::u,stat
    open(newunit=u,file=trim(output_path),status='replace',action='write',iostat=stat)
    if(stat/=0)then
      print '(a)','invalid input and certificate open failed';stop 2
    endif
    write(u,'(a)')'{"schema":"finite_cold_prefix_audit_v1","status":"INVALID_INPUT",'
    write(u,'(a,i0,a,i0,a)')'"input_valid":false,"rows":',n,',"input_line_one_based":',line_number,','
    write(u,'(a)')'"message":"'//trim(message)//'"}'
    close(u)
    print '(a)',trim(message);stop 2
  end subroutine
end program source_prefix_audit
