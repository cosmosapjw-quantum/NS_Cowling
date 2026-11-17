program test_matching_extrema_controls
  use compose_eos_types
  use compose_eos_matching
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  real(wp)::c(0:5),ymin,coefficients(0:4,1),a
  logical::ok
  character(len=512)::msg
  ! Exact stationary points: flat quartic, two minima, narrow negative
  ! pocket, and endpoint minimum. No uniform grid is used as an oracle.
  a=.37_wp;c=[a**4,-4._wp*a**3,6._wp*a*a,-4._wp*a,1._wp,0._wp]
  ymin=polynomial_minimum(c)
  if(abs(ymin)>2.e-13_wp)stop 3
  c=[.0156_wp,-.32_wp,1.32_wp,-2._wp,1._wp,0._wp]
  if(abs(polynomial_minimum(c)+.01_wp)>2.e-12_wp)stop 4
  a=.123456789_wp
  c=[100._wp*a*a+10._wp*a**4-1.e-6_wp,-200._wp*a-40._wp*a**3, &
     100._wp+60._wp*a*a,-40._wp*a,10._wp,0._wp]
  if(abs(polynomial_minimum(c)+1.e-6_wp)>5.e-12_wp)stop 5
  c=[1._wp,-2._wp,0._wp,0._wp,0._wp,0._wp]
  if(abs(polynomial_minimum(c)+1._wp)>1.e-14_wp)stop 6
  coefficients=0._wp;coefficients(0,1)=ieee_value(0._wp,ieee_quiet_nan)
  call certify_profile(coefficients,[0._wp,1._wp],.true.,0._wp,1._wp,939._wp,ok,msg)
  if(ok)stop 7
  print '(a)','matching polynomial extrema controls: PASS'
end program
