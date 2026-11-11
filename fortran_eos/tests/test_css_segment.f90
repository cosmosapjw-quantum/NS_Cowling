program test_css_segment
  use compose_eos_types, only: wp
  use compose_eos_css
  implicit none
  type(css_segment_t) :: segment, zero_jump, invalid
  real(wp) :: pressure(3), epsilon(3), nb(3), mu_b(3)
  real(wp) :: pressure_fd(3), epsilon_fd(3), nb_fd(3), mu_fd(3)
  real(wp) :: de_dn, cs2_fd
  logical :: ok
  character(len=256) :: message

  call initialize_css_segment(1._wp, 940._wp, 1._wp, 941._wp, 0.25_wp, &
                              segment, ok, message)
  call require(ok, 'valid CSS initialization: '//trim(message))
  call require(segment%interface_class == css_interface_first_order, &
               'positive Delta epsilon must be first-order')
  call close_to(segment%mu_transition, 941._wp, 1.e-13_wp, 'mu transition')
  call close_to(segment%epsilon_high_transition, 1881._wp, 1.e-13_wp, 'epsilon plus')
  call close_to(segment%nb_high_transition, 2._wp, 1.e-13_wp, 'n plus')

  pressure = [1._wp, 2._wp, 4._wp]
  call evaluate_css_segment(segment, pressure, epsilon, nb, mu_b, ok, message)
  call require(ok, 'valid CSS branch evaluation: '//trim(message))
  call close_to(epsilon(1), 1881._wp, 1.e-13_wp, 'transition epsilon')
  call close_to(nb(1), 2._wp, 1.e-13_wp, 'transition density')
  call close_to(mu_b(1), 941._wp, 1.e-13_wp, 'chemical-potential continuity')
  cs2_fd = (pressure(3) - pressure(1)) / (epsilon(3) - epsilon(1))
  call close_to(cs2_fd, 0.25_wp, 1.e-12_wp, 'constant sound speed')

  pressure_fd = [1.99_wp, 2._wp, 2.01_wp]
  call evaluate_css_segment(segment, pressure_fd, epsilon_fd, nb_fd, mu_fd, ok, message)
  call require(ok, 'finite-difference CSS evaluation: '//trim(message))
  de_dn = (epsilon_fd(3) - epsilon_fd(1)) / (nb_fd(3) - nb_fd(1))
  call close_to(de_dn, mu_fd(2), 2.e-4_wp, 'first-law derivative d epsilon/d n = mu')

  call initialize_css_segment(1._wp, 940._wp, 1._wp, 0._wp, 0.25_wp, &
                              zero_jump, ok, message)
  call require(ok, 'zero-jump CSS initialization: '//trim(message))
  call require(zero_jump%interface_class == css_interface_continuous_unclassified, &
               'zero jump must remain thermodynamically unclassified')
  call close_to(zero_jump%nb_high_transition, 1._wp, 1.e-13_wp, 'zero-jump density continuity')

  pressure = [0.9_wp, 1._wp, 2._wp]
  call evaluate_css_segment(segment, pressure, epsilon, nb, mu_b, ok, message)
  call require(.not. ok, 'high-phase evaluator must reject P<P_t')

  call initialize_css_segment(1._wp, 940._wp, 1._wp, 1._wp, 1.01_wp, &
                              invalid, ok, message)
  call require(.not. ok, 'acausal CSS speed must be rejected')
  call initialize_css_segment(1._wp, 940._wp, 1._wp, -1._wp, 0.25_wp, &
                              invalid, ok, message)
  call require(.not. ok, 'negative CSS density jump must be rejected')

  print '(a)', 'CSS segment analytic checks: PASS'

contains

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(label)
      error stop 1
    endif
  end subroutine require

  subroutine close_to(value, expected, tolerance, label)
    real(wp), intent(in) :: value, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(value - expected) > tolerance * max(1._wp, abs(expected))) then
      write(*,'(a,es24.16,a,es24.16)') 'FAIL: '//trim(label)//' actual=', value, ' expected=', expected
      error stop 1
    endif
  end subroutine close_to

end program test_css_segment
