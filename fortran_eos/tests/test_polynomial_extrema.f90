program test_polynomial_extrema
  use compose_eos_types, only: wp
  use compose_eos_reconstruct, only: quartic_polynomial_minimum
  implicit none
  real(wp)::ymin,a,c0,c1,c2,c3,c4,grid_min,value
  integer::k
  a=0.37_wp
  c0=a**4;c1=-4._wp*a**3;c2=6._wp*a**2;c3=-4._wp*a;c4=1._wp
  call quartic_polynomial_minimum(c0,c1,c2,c3,c4,ymin)
  if(abs(ymin)>2.e-13_wp)error stop 'flat quartic minimum was not isolated'
  call quartic_polynomial_minimum(0.0156_wp,-0.32_wp,1.32_wp,-2._wp,1._wp,ymin)
  if(abs(ymin+0.01_wp)>2.e-12_wp)error stop 'two-minimum quartic failed'
  a=0.123456789_wp
  c0=100._wp*a*a+10._wp*a**4-1.e-6_wp
  c1=-200._wp*a-40._wp*a**3
  c2=100._wp+60._wp*a*a
  c3=-40._wp*a
  c4=10._wp
  grid_min=huge(1._wp)
  do k=0,256
    value=real(k,wp)/256._wp
    grid_min=min(grid_min,((((c4*value+c3)*value+c2)*value+c1)*value+c0))
  enddo
  if(grid_min<=0._wp)error stop 'adversarial narrow pocket unexpectedly sampled by old mesh'
  call quartic_polynomial_minimum(c0,c1,c2,c3,c4,ymin)
  if(abs(ymin+1.e-6_wp)>5.e-12_wp)error stop 'narrow negative quartic pocket was missed'
  print '(a)', 'quartic polynomial extrema: PASS'
end program test_polynomial_extrema
