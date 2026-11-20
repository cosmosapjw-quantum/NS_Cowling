! Independently exercise original EOS c_s² law with no Python runtime.
program cowling_source_exact_validation
  use iso_c_binding, only:c_double,c_int
  use, intrinsic::ieee_arithmetic,only:ieee_is_finite
  use compose_eos_cowling_source_exact,only: &
        cowling_native_source_eos_setup,cowling_native_source_cs2_eval
  implicit none
  character(len=1024)::path
  integer::u,ios,no,nu,jump,i,nprobe,sector_expected,sector_fail
  integer::side
  real(c_double),allocatable::outer(:,:),upper(:,:)
  real(c_double)::surface(4),bridge(7),h,expected,got,scale,maxerr
  integer(c_int)::status,region
  call get_command_argument(1,path)
  if(len_trim(path)==0)stop 10
  open(newunit=u,file=trim(path),status="old",action="read",iostat=ios)
  if(ios/=0)stop 11
  read(u,*,iostat=ios)no,nu,jump
  if(ios/=0.or.no<2.or.nu<4)stop 12
  allocate(outer(4,no),upper(3,nu))
  read(u,*,iostat=ios)surface
  if(ios/=0)stop 13
  read(u,*,iostat=ios)bridge
  if(ios/=0)stop 14
  do i=1,no
    read(u,*,iostat=ios)outer(:,i)
    if(ios/=0)stop 15
  enddo
  do i=1,nu
    read(u,*,iostat=ios)upper(:,i)
    if(ios/=0)stop 16
  enddo
  call cowling_native_source_eos_setup(int(no,c_int),outer, &
       int(nu,c_int),upper,surface,bridge,int(jump,c_int),status)
  if(status/=0)then
    print *,'EOS_SOURCE_SETUP_REJECTED',status
    stop 17
  endif
  read(u,*,iostat=ios)nprobe
  if(ios/=0.or.nprobe<100)stop 18
  maxerr=0._c_double
  sector_fail=0
  do i=1,nprobe
    read(u,*,iostat=ios)h,side,expected,sector_expected
    if(ios/=0)stop 19
    call cowling_native_source_cs2_eval(h,int(side,c_int),got,region,status)
    if(status/=0.or..not.ieee_is_finite(got))then
      print *,'SOURCE_EVAL_FAILED',i,h,side,status
      stop 20
    endif
    if(region/=sector_expected)then
      print *,'SOURCE_REGION_MISMATCH',i,h,region,sector_expected
      stop 21
    endif
    scale=max(abs(expected),1.e-5_c_double)
    maxerr=max(maxerr,abs(got-expected)/scale)
    if(abs(got-expected)>max(5.e-12_c_double,1.e-8_c_double*abs(expected)))then
      print *,'SOURCE_SPEED_MISMATCH',i,h,got,expected
      stop 22
    endif
  enddo
  close(u)
  call cowling_native_source_cs2_eval(-1._c_double,2_c_int,got,region,status)
  if(status==0.or.got/=0._c_double.or.region/=0)stop 23
  print '("FORTRAN_SOURCE_LAW_PASS PROBES ",I5," MAX_REL_ERROR ",ES13.5)', &
         nprobe,maxerr
end program cowling_source_exact_validation
