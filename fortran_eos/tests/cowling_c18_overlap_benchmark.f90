! Isolated native-Fortran kernel microbenchmark, no DVODE and no Python.
! A synthetic two-phase profile probes interpolation/normalization throughput.
program cowling_c18_overlap_benchmark
  use iso_c_binding, only:c_double,c_int
  use omp_lib, only:omp_get_wtime,omp_get_max_threads
  use compose_eos_cowling_mode_overlap, only:cowling_native_mode_shape_overlap
  implicit none
  integer,parameter::n=360,split=180
  integer(c_int)::phase(n),status
  real(c_double)::r(n),w(n),v(n),total,inner,outer,check
  real(c_double)::t0,t1
  character(len=64)::arg
  integer::reps,j,it,ios,bad
  reps=20000
  call get_command_argument(1,arg)
  if(len_trim(arg)>0)then
    read(arg,*,iostat=ios)reps
    if(ios/=0.or.reps<1.or.reps>1000000)stop 2
  endif
  do j=1,n
    if(j<=split)then
      phase(j)=0_c_int
      r(j)=real(j-1,c_double)/real(split-1,c_double)
    else
      phase(j)=1_c_int
      r(j)=1._c_double+real(j-split-1,c_double)/real(n-split-1,c_double)
    endif
    w(j)=sin(real(j,c_double)/57._c_double)+.15_c_double
    v(j)=cos(real(j,c_double)/71._c_double)-.25_c_double
  enddo
  bad=0;check=0._c_double
  t0=omp_get_wtime()
  !$omp parallel do default(shared) private(total,inner,outer,status) &
  !$omp reduction(+:check,bad) schedule(static)
  do it=1,reps
    call cowling_native_mode_shape_overlap(int(n,c_int),phase,r,w, &
      int(n,c_int),phase,r,v,total,inner,outer,status)
    if(status/=0)then
      bad=bad+1
    else
      check=check+total+inner+outer
    endif
  enddo
  !$omp end parallel do
  t1=omp_get_wtime()
  if(bad/=0)stop 3
  write(*,'("N=",I0," REP=",I0," THREADS=",I0, &
     " WALL_S=",F14.7," CHECKSUM=",ES24.15)') &
     n,reps,omp_get_max_threads(),t1-t0,check
end program cowling_c18_overlap_benchmark
