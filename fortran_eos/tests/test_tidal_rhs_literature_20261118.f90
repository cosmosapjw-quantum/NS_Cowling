! Independent GR Hinderer Q(r) identity at a fixed interior stellar state.
! Old native code erroneously used dP/dE instead of dE/dP.
program test_tidal_rhs_literature_20261118
  use compose_eos_types, only: wp
  use compose_eos_tov, only: tov_rhs
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real(wp),parameter :: k=1.3243117038949245e-6_wp
  real(wp),parameter :: pi=acos(-1._wp)
  real(wp) :: p(3),eps(3),r,pm,m,y,fac,elam,nu,eg,dedp
  real(wp) :: dpdr,dmdr,dydr,expected,Q,wrongQ,wrong_rhs,ratio
  logical :: ok

  p=[1._wp,2._wp,4._wp]
  eps=[100._wp,104._wp,110._wp]
  r=12._wp
  pm=2._wp*k
  m=1.7_wp
  y=1.8_wp
  call tov_rhs(r,pm,m,y,p,eps,3,p(1)*k,dpdr,dmdr,dydr,ok)
  if(.not.ok.or..not.ieee_is_finite(dydr))stop 3

  fac=1._wp-2._wp*m/r
  elam=1._wp/fac
  nu=2._wp*(m+4._wp*pi*r**3*pm)/(r*(r-2._wp*m))
  dedp=(104._wp/2._wp)*log(110._wp/104._wp)/log(4._wp/2._wp)
  eg=104._wp*k
  Q=4._wp*pi*elam*(5._wp*eg+9._wp*pm+(eg+pm)*dedp) &
     -6._wp*elam/(r*r)-nu**2
  expected=-(y*y+y*elam*(1._wp+4._wp*pi*r**2*(pm-eg))+r*r*Q)/r
  wrongQ=4._wp*pi*elam*(5._wp*eg+9._wp*pm+(eg+pm)/dedp) &
     -6._wp*elam/(r*r)-nu**2
  wrong_rhs=-(y*y+y*elam*(1._wp+4._wp*pi*r**2*(pm-eg))+r*r*wrongQ)/r
  ratio=abs((dydr-expected)/expected)
  write(*,'(a,3es24.14)')'Q_PHYSICS_VALIDATION ',dydr,expected,ratio
  if(ratio>1.e-12_wp)then
    write(*,'(a)')'FAIL tidal dE/dP source-GR identity'
    stop 1
  endif
  if(abs(dydr-wrong_rhs)<1.e-6_wp*abs(dydr))then
    write(*,'(a)')'FAIL old inverse derivative bug cannot be detected by test'
    stop 2
  endif
  write(*,'(a)')'PASS tidal RHS independent Hinderer Q identity'
end program test_tidal_rhs_literature_20261118
