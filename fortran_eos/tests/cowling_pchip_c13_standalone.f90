! Self-contained Fortran only PCHIP verification on exact C12 star samples.
! The expected derivatives and field values were exported once from the
! independently source-locked SciPy reference, not generated at runtime.
program cowling_c13_pchip_standalone
  use iso_c_binding, only:c_int,c_double
  use, intrinsic :: ieee_arithmetic, only:ieee_is_finite
  use compose_eos_cowling_pchip, only: &
       cowling_native_pchip_slopes,cowling_native_pchip_values
  implicit none
  character(len=1024)::arg
  integer::u,ios,n,m,i
  real(c_double),allocatable::s(:),y(:,:),expected_d(:,:),d(:,:),q(:),expected_v(:,:),v(:,:)
  real(c_double)::max_scaled,max_abs,scaled,delta
  integer(c_int)::status
  call get_command_argument(1,arg)
  if(len_trim(arg)==0)stop 10
  open(newunit=u,file=trim(arg),action='read',status='old',iostat=ios)
  if(ios/=0)stop 11
  read(u,*,iostat=ios)n,m
  if(ios/=0.or.n<3.or.m<2)stop 12
  allocate(s(n),y(5,n),expected_d(5,n),d(5,n), &
      q(m),expected_v(5,m),v(5,m))
  do i=1,n
    read(u,*,iostat=ios)s(i),y(:,i),expected_d(:,i)
    if(ios/=0)stop 13
  enddo
  do i=1,m
    read(u,*,iostat=ios)q(i),expected_v(:,i)
    if(ios/=0)stop 14
  enddo
  close(u)
  call cowling_native_pchip_slopes(int(n,c_int),s,y,d,status)
  if(status/=0)stop 15
  max_scaled=maxval(abs(d-expected_d)/max(abs(expected_d),1.e-8_c_double))
  max_abs=maxval(abs(d-expected_d))
  if(any(abs(d-expected_d)>5.e-9_c_double + &
     2.e-10_c_double*abs(expected_d)))then
    print *, 'C13_BAD_PCHIP_D',max_scaled,max_abs
    stop 16
  endif
  call cowling_native_pchip_values(int(n,c_int),s,y,d,int(m,c_int),q,v,status)
  if(status/=0)stop 17
  if(any(abs(v-expected_v)>1.e-10_c_double + &
    2.e-10_c_double*abs(expected_v)))then
    print *, 'C13_BAD_PCHIP_VALUES'
    stop 18
  endif
  q(1)=s(1)-.01_c_double
  call cowling_native_pchip_values(int(n,c_int),s,y,d,int(m,c_int),q,v,status)
  if(status==0.or.any(v/=0._c_double))stop 19
  s(2)=s(1)
  call cowling_native_pchip_slopes(int(n,c_int),s,y,d,status)
  if(status==0.or.any(d/=0._c_double))stop 20
  print '("C13_FORTRAN_PCHIP_ORACLE_PASS N ",I5," QUERIES ",I5, &
      " MAX_SCALED_SLOPE_DELTA ",ES13.5)',n,m,max_scaled
end program cowling_c13_pchip_standalone
