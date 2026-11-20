! Fortran-only characterization and process-isolated parallel work unit.
! Synthetic smooth geometry; results are numerical regression, not physical modes.
program cowling_c8_native_probe
  use iso_c_binding, only: c_double,c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use compose_eos_cowling_dvode, only: cowling_native_setup,cowling_native_shot
  implicit none
  integer, parameter :: n=8193, np=2
  real(c_double), parameter :: pi=3.1415926535897932384626433832795_c_double
  real(c_double), parameter :: ckm=299792.458_c_double
  real(c_double) :: s(n),v(5,n),dv(5,n),si(np),so(np),inner(2,np),outer(2,np)
  real(c_double) :: omega,residual,checksum,t0,t1
  integer(c_int) :: setup_status,status,nfev,nfev_total
  integer :: i,j,job,shots,ios,failures
  character(len=60) :: arg

  job=1;shots=500
  call get_command_argument(1,arg)
  if(len_trim(arg)>0)then
    read(arg,*,iostat=ios)job
    if(ios/=0)stop 3
  endif
  call get_command_argument(2,arg)
  if(len_trim(arg)>0)then
    read(arg,*,iostat=ios)shots
    if(ios/=0)stop 3
  endif
  if(job<1.or.job>1000.or.shots<1.or.shots>100000)stop 3

  do i=1,n
    s(i)=real(i-1,c_double)/real(n-1,c_double)
    ! Smooth synthetic coefficients. Native interpolation uses the exact
    ! same Hermite values and slopes before and after each optimization.
    v(:,i)=[0.1_c_double+12._c_double*s(i), &
            -0.3_c_double+0.08_c_double*s(i), &
             0.05_c_double+0.12_c_double*s(i), &
             0.011_c_double+0.004_c_double*s(i), &
             0.12_c_double+0.035_c_double*s(i)]
    dv(:,i)=[12._c_double,0.08_c_double,0.12_c_double, &
              0.004_c_double,0.035_c_double]
  enddo
  call cowling_native_setup(n,s,v,dv,setup_status)
  if(setup_status/=0)then
    print *, 'SETUP_FAILURE',setup_status
    stop 5
  endif

  si=[0._c_double,0.5_c_double]
  so=[1._c_double,0.5_c_double]
  nfev_total=0;checksum=0._c_double;failures=0
  call cpu_time(t0)
  do j=1,shots
    omega=2._c_double*pi*(1700._c_double + &
          10._c_double*real(mod(j+13*job,105),c_double))/ckm
    call cowling_native_shot(omega,2_c_int,1.e-9_c_double, &
         1.e-12_c_double,0.005_c_double,10_c_int, &
         np,si,np,so,inner,outer,residual,nfev,status)
    if(status/=0.or..not.ieee_is_finite(residual))then
      failures=failures+1
      exit
    endif
    nfev_total=nfev_total+nfev
    checksum=checksum+residual*real(mod(j,17)+1,c_double)
  enddo
  call cpu_time(t1)
  if(failures/=0)then
    print *, 'SHOOT_FAILURE',job,j,status,failures
    stop 7
  endif
  write(*,'("C8_JOB ",I4," SHOTS ",I7," NFEV ",I12," CHECKSUM ",ES25.16)') &
       job,shots,nfev_total,checksum
  write(*,'("C8_CPU_SECONDS ",ES14.6)') t1-t0
end program cowling_c8_native_probe
