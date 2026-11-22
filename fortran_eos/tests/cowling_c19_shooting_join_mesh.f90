! C19 physical regression: C17 two-sided shooting joins duplicate radius.
! The duplicate is inside the high Maxwell phase, not a Maxwell jump.
! W sides differ by roundoff / finite matching error and must be compared
! without ever interpolating across the high/low physical phase interface.
program cowling_c19_shooting_join_mesh
  use iso_c_binding, only:c_double,c_int
  use compose_eos_cowling_mode_overlap, only:cowling_native_mode_shape_overlap
  implicit none
  integer(c_int)::phase(10),status
  real(c_double)::r(10),w(10),z(10),bad(10),bad_r(10)
  real(c_double)::total,high,low
  phase=[0_c_int,0_c_int,0_c_int,0_c_int,0_c_int, &
         1_c_int,1_c_int,1_c_int,1_c_int,1_c_int]
  r=[0._c_double,.2_c_double,.4_c_double,.4_c_double,.7_c_double, &
     .7_c_double,.8_c_double,.9_c_double,1._c_double,1.1_c_double]
  w=[0._c_double,.7_c_double,1._c_double,1._c_double+1.e-8_c_double, &
     .4_c_double,.4_c_double,.2_c_double,-.1_c_double,-.4_c_double,-.2_c_double]
  z=-w
  call cowling_native_mode_shape_overlap(10_c_int,phase,r,w, &
     10_c_int,phase,r,z,total,high,low,status)
  if(status/=0)then
    print *, 'C19_FAIL_C17_SHOOTING_JOIN_REJECTED',status
    stop 10
  endif
  if(abs(total-1._c_double)>5.e-13_c_double.or. &
     abs(high-1._c_double)>5.e-13_c_double.or. &
     abs(low-1._c_double)>5.e-13_c_double)stop 11

  ! The shared shooting match must remain well-defined under independently
  ! normalized tiny/huge eigenfunction amplitudes, too.
  call cowling_native_mode_shape_overlap(10_c_int,phase,r,w*1.e200_c_double, &
     10_c_int,phase,r,-w*1.e-200_c_double,total,high,low,status)
  if(status/=0.or.abs(total-1._c_double)>5.e-13_c_double.or. &
     abs(high-1._c_double)>5.e-13_c_double.or. &
     abs(low-1._c_double)>5.e-13_c_double)stop 15

  ! A large jump at a repeated radius cannot represent a continuous
  ! numerical shooting join: must fail closed.
  bad=w;bad(4)=bad(3)+1.e-2_c_double
  call cowling_native_mode_shape_overlap(10_c_int,phase,r,bad, &
     10_c_int,phase,r,w,total,high,low,status)
  if(status==0)stop 12

  ! Two different internal repeated radii are not one shooting junction.
  bad_r=r;bad_r(5)=bad_r(4)
  call cowling_native_mode_shape_overlap(10_c_int,phase,bad_r,w, &
     10_c_int,phase,r,w,total,high,low,status)
  if(status==0)stop 13

  ! Invalid downward step even inside one phase remains prohibited.
  bad_r=r;bad_r(4)=r(3)-1.e-3_c_double
  call cowling_native_mode_shape_overlap(10_c_int,phase,bad_r,w, &
     10_c_int,phase,r,w,total,high,low,status)
  if(status==0)stop 14
  ! Degenerate all-zero W shape remains a typed status=2 rather than a
  ! numerically spurious perfect overlap.
  z=0._c_double
  call cowling_native_mode_shape_overlap(10_c_int,phase,r,w, &
     10_c_int,phase,r,z,total,high,low,status)
  if(status/=2_c_int)stop 16
  print *, 'C19_REAL_SHOOTING_JOIN_MESH_AND_MAXWELL_SCOPE_PASS'
end program cowling_c19_shooting_join_mesh
