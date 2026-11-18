! Adversarial source-independent Hinderer Love normalization gate.
! The reference formula here is independently written from the paper and
! evaluated using the NATIVE SOLVER's beta and y_R, so failures localize
! to k2 algebra before testing the separate tidal ODE.
program test_tov_love_literature_20261118
  use compose_eos_types, only: wp
  use compose_eos_tov, only: tov_love_result_t,solve_tov_love_smooth
  implicit none
  integer,parameter::n=96
  real(wp):: pressure(n),epsilon(n),t,b,y,coef,denom,k2_expected,ratio
  integer::i
  type(tov_love_result_t)::result
  do i=1,n
    t=real(i-1,wp)/real(n-1,wp)
    epsilon(i)=100.1_wp*exp(log(1600._wp/100.1_wp)*t)
    pressure(i)=(epsilon(i)-100._wp)/3._wp
  enddo
  call solve_tov_love_smooth(n,pressure,epsilon,pressure(n/2),result)
  if(.not.result%converged)then
    write(*,'(a,a)')'FAIL native TOV baseline did not converge: ',trim(result%status)
    stop 2
  endif
  b=result%compactness
  y=result%y_surface
  coef=2._wp*b*(y-1._wp)-y+2._wp
  denom=2._wp*b*(4._wp*(y+1._wp)*b**4+(6._wp*y-4._wp)*b**3+ &
      (26._wp-22._wp*y)*b**2+3._wp*(5._wp*y-8._wp)*b-3._wp*y+6._wp) + &
      3._wp*(1._wp-2._wp*b)**2*coef*log(1._wp-2._wp*b)
  k2_expected=(8._wp/5._wp)*b**5*(1._wp-2._wp*b)**2*coef/denom
  ratio=result%k2/k2_expected
  write(*,'(a,3es25.16)')'NATIVE_LOVE_K2_NORMALIZATION ',result%k2,k2_expected,ratio
  if(abs(ratio-1._wp)>1.e-8_wp)then
    write(*,'(a)')'FAIL Hinderer analytic k2 normalization'
    stop 1
  endif
  write(*,'(a)')'PASS analytic Love k2 normalization'
end program
