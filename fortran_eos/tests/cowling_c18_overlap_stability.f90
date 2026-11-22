! C18: native mode-shape overlap must be safe under finite amplitude scaling.
! The Maxwell radius may repeat only across the physical phase boundary.
program cowling_c18_overlap_stability
  use iso_c_binding, only:c_double,c_int
  use compose_eos_cowling_mode_overlap, only:cowling_native_mode_shape_overlap
  implicit none
  integer(c_int) :: phase(8),status
  real(c_double) :: r(8),w(8),scaled(8),r_bad(8),total,inner,outer
  real(c_double) :: base_total,base_inner,base_outer
  phase=[0_c_int,0_c_int,0_c_int,0_c_int,1_c_int,1_c_int,1_c_int,1_c_int]
  r=[0._c_double,.2_c_double,.4_c_double,.6_c_double, &
      .6_c_double,.8_c_double,1._c_double,1.2_c_double]
  w=[0._c_double,.7_c_double,1._c_double,.4_c_double, &
      .4_c_double,-.1_c_double,-.4_c_double,-.2_c_double]
  call cowling_native_mode_shape_overlap(8_c_int,phase,r,w, &
       8_c_int,phase,r,-w,total,inner,outer,status)
  if(status/=0.or.abs(total-1._c_double)>1.e-13_c_double)stop 10
  if(abs(inner-1._c_double)>1.e-13_c_double.or. &
     abs(outer-1._c_double)>1.e-13_c_double)stop 11
  base_total=total; base_inner=inner; base_outer=outer

  ! Finite input amplitudes must not overflow the quadratic norms.
  scaled=w*1.e200_c_double
  call cowling_native_mode_shape_overlap(8_c_int,phase,r,scaled, &
       8_c_int,phase,r,-scaled,total,inner,outer,status)
  if(status/=0)then
    print *, 'C18_RED_HUGE_FINITE_STATUS',status
    stop 12
  endif
  if(abs(total-base_total)>1.e-12_c_double.or. &
      abs(inner-base_inner)>1.e-12_c_double.or. &
      abs(outer-base_outer)>1.e-12_c_double)stop 13

  ! Finite tiny amplitudes must not disappear in squared norms.
  scaled=w*1.e-200_c_double
  call cowling_native_mode_shape_overlap(8_c_int,phase,r,scaled, &
       8_c_int,phase,r,-scaled,total,inner,outer,status)
  if(status/=0)then
    print *, 'C18_RED_TINY_FINITE_STATUS',status
    stop 14
  endif
  if(abs(total-base_total)>1.e-12_c_double.or. &
      abs(inner-base_inner)>1.e-12_c_double.or. &
      abs(outer-base_outer)>1.e-12_c_double)stop 15

  ! Unequal global amplitudes across different eigenfunction normalizations
  ! must cancel from all three correlations without overflow or underflow.
  call cowling_native_mode_shape_overlap(8_c_int,phase,r,w*1.e200_c_double, &
       8_c_int,phase,r,-w*1.e-200_c_double,total,inner,outer,status)
  if(status/=0.or.abs(total-base_total)>1.e-12_c_double.or. &
      abs(inner-base_inner)>1.e-12_c_double.or. &
      abs(outer-base_outer)>1.e-12_c_double)stop 17

  ! A repeated radius inside a phase is not a physical mesh.
  ! The duplicate at the Maxwell boundary (r(4)==r(5)) remains allowed.
  r_bad=r
  r_bad(3)=r_bad(2)
  call cowling_native_mode_shape_overlap(8_c_int,phase,r_bad,w, &
       8_c_int,phase,r,w,total,inner,outer,status)
  if(status==0)then
    print *, 'C18_RED_DUPLICATE_INTRAPHASE_ACCEPTED'
    stop 16
  endif
  print *, 'C18_OVERLAP_STABILITY_AND_PHASE_MESH_PASS'
end program cowling_c18_overlap_stability
