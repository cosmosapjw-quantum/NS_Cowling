! C17 pure Fortran auditable onset and branch selector unit driver.
! No Python process, source contract values selected from original M2/M3.
program c17_branch_and_seidov_original_source
  use iso_c_binding, only:c_double,c_int
  use compose_eos_cowling_extended_branch, only: &
       cowling_native_maxwell_onset_threshold, &
       cowling_native_select_mode_branch_overlap
  implicit none
  real(c_double)::delta,delta_c,margin
  real(c_double)::f(4),high(4),low(4)
  integer(c_int)::flag,st,index,nw(4),nv(4),rootstatus(4)
  call cowling_native_maxwell_onset_threshold( &
    8.572084498151984_c_double,16.33826505200197_c_double, &
    .014307398278768153_c_double,delta,delta_c,margin,flag,st)
  if(st/=0.or.flag/=1)stop 10
  if(abs(delta-.9044754359602759_c_double)>1.e-13_c_double)stop 11
  if(abs(delta_c-.5016662875921976_c_double)>1.e-13_c_double)stop 12
  print '("C17_SOURCE_M2_DELTA ",F14.10," DELTA_CRIT ",F14.10, &
       " SMALL_CORE_ABOVE 1")',delta,delta_c
  call cowling_native_maxwell_onset_threshold( &
    76.94982713409907_c_double,89.34565739599732_c_double, &
    .5657052249953693_c_double,delta,delta_c,margin,flag,st)
  if(st/=0.or.flag/=0)stop 13
  if(abs(delta-.15991414732822803_c_double)>1.e-13_c_double)stop 14
  if(abs(delta_c-.507297959618915_c_double)>1.e-13_c_double)stop 15
  print '("C17_SOURCE_M3_DELTA ",F14.10," DELTA_CRIT ",F14.10, &
       " SMALL_CORE_ABOVE 0")',delta,delta_c
  f=[1200._c_double,1850._c_double,3730._c_double,4200._c_double]
  nw=[0,0,1,1]
  nv=[0,0,1,1]
  rootstatus=[1,0,0,0]
  high=[0.99_c_double,0.99_c_double,0.989_c_double,0.7_c_double]
  low=[0.99_c_double,0.99_c_double,0.999_c_double,0.999_c_double]
  call cowling_native_select_mode_branch_overlap(4_c_int,f,nw,nv,rootstatus,high,low, &
       1_c_int,3730._c_double,.12_c_double,.95_c_double,index,st)
  if(st/=0.or.index/=2)stop 16
  rootstatus(4)=0
  high(4)=.99_c_double
  f(4)=3750._c_double
  call cowling_native_select_mode_branch_overlap(4_c_int,f,nw,nv,rootstatus,high,low, &
       1_c_int,3730._c_double,.12_c_double,.95_c_double,index,st)
  if(st/=3.or.index/=-1)stop 17
  high(3)=.40_c_double
  rootstatus(4)=1
  call cowling_native_select_mode_branch_overlap(4_c_int,f,nw,nv,rootstatus,high,low, &
       1_c_int,3730._c_double,.12_c_double,.95_c_double,index,st)
  if(st/=5.or.index/=-1)stop 18
  print '(A)','C17_NATIVE_ORIGINAL_SOURCE_ONSET_AND_OVERLAP_SELECTOR_PASS'
end program c17_branch_and_seidov_original_source
