! Independent Fortran C9 bulk scan regression with bounds checking.
! No Python numerical code or shared-thread DVODE execution.
program cowling_c9_bulk_scan_test
  use iso_c_binding, only: c_int, c_double
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use compose_eos_cowling_dvode, only: cowling_native_setup, &
      cowling_native_shot, cowling_native_scan
  implicit none
  integer,parameter :: n=1025,nf=6
  real(c_double),parameter :: pi=3.1415926535897932384626433832795_c_double
  real(c_double),parameter :: ckm=299792.458_c_double
  real(c_double) :: s(n),values(5,n),slopes(5,n)
  real(c_double) :: f(nf),results(nf),si(2),so(2),u(2,2),v(2,2),r,omega
  integer(c_int) :: init_status,flags(nf-1),nfev_total,scan_status
  integer(c_int) :: nfev,shot_status
  integer :: i,method,method_flag,expected_sum,expected_crossings(nf-1)
  do i=1,n
    s(i)=real(i-1,c_double)/real(n-1,c_double)
    values(:,i)=[0.1_c_double+12*s(i),-.3_c_double+.08_c_double*s(i), &
       .05_c_double+.12_c_double*s(i),.011_c_double+.004_c_double*s(i), &
       .12_c_double+.035_c_double*s(i)]
    slopes(:,i)=[12._c_double,.08_c_double,.12_c_double, &
                   .004_c_double,.035_c_double]
  enddo
  call cowling_native_setup(n,s,values,slopes,init_status)
  if(init_status/=0)stop 10
  f=[1800._c_double,2000._c_double,2300._c_double, &
     2600._c_double,3100._c_double,4100._c_double]
  si=[0._c_double,.5_c_double]
  so=[1._c_double,.5_c_double]
  do method=1,2
    method_flag=10
    if(method==2)method_flag=22
    call cowling_native_scan(nf,f,2_c_int,1.e-9_c_double,1.e-12_c_double, &
      .005_c_double,int(method_flag,c_int),.5_c_double,1._c_double, &
      results,flags,nfev_total,scan_status)
    if(scan_status/=0)then
      print *, 'SCAN_FAILED',scan_status,method
      stop 11
    endif
    if(.not.all(ieee_is_finite(results)))stop 12
    expected_sum=0
    expected_crossings=0
    do i=1,nf
      omega=2._c_double*pi*f(i)/ckm
      call cowling_native_shot(omega,2_c_int,1.e-9_c_double,1.e-12_c_double, &
        .005_c_double,int(method_flag,c_int),2_c_int,si,2_c_int,so, &
        u,v,r,nfev,shot_status)
      if(shot_status/=0)stop 13
      expected_sum=expected_sum+nfev
      if(abs(results(i)-r)>3.e-13_c_double)then
        print *, 'NUMERIC_MISMATCH',method,i,results(i),r
        stop 14
      endif
      if(i>1)then
        if(results(i-1)*results(i)<0._c_double)expected_crossings(i-1)=1
      endif
    enddo
    if(expected_sum/=nfev_total)stop 15
    if(any(flags/=expected_crossings))stop 16
  enddo
  f(3)=-1._c_double
  call cowling_native_scan(nf,f,2_c_int,1.e-9_c_double,1.e-12_c_double, &
    .005_c_double,10_c_int,.5_c_double,1._c_double, &
    results,flags,nfev_total,scan_status)
  if(scan_status==0.or.nfev_total/=0)stop 17
  if(any(results/=0._c_double).or.any(flags/=0))stop 18
  print *, 'FORTRAN_C9_BULK_SCAN_PASS_ADAMS_BDF_REJECT_BAD_INPUT'
end program cowling_c9_bulk_scan_test
