module compose_eos_css
  use compose_eos_types, only: wp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer, parameter, public :: css_interface_continuous_unclassified = 0
  integer, parameter, public :: css_interface_first_order = 1

  type, public :: css_segment_t
    real(wp) :: p_transition = 0._wp
    real(wp) :: epsilon_low = 0._wp
    real(wp) :: nb_low = 0._wp
    real(wp) :: delta_epsilon = 0._wp
    real(wp) :: cs2 = 0._wp
    real(wp) :: mu_transition = 0._wp
    real(wp) :: epsilon_high_transition = 0._wp
    real(wp) :: nb_high_transition = 0._wp
    integer :: interface_class = css_interface_continuous_unclassified
    logical :: initialized = .false.
  end type css_segment_t

  public :: initialize_css_segment, evaluate_css_segment

contains

  subroutine initialize_css_segment(p_transition, epsilon_low, nb_low, &
                                    delta_epsilon, cs2, segment, ok, message)
    real(wp), intent(in) :: p_transition, epsilon_low, nb_low
    real(wp), intent(in) :: delta_epsilon, cs2
    type(css_segment_t), intent(out) :: segment
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(wp) :: mu_t, epsilon_high, nb_high, jump_residual, scale

    ok = .false.
    message = 'invalid CSS segment input'
    segment%initialized = .false.

    if (.not. all(ieee_is_finite([p_transition, epsilon_low, nb_low, &
                                  delta_epsilon, cs2]))) then
      message = 'CSS parameters must be finite'
      return
    endif
    if (p_transition < 0._wp .or. epsilon_low <= 0._wp .or. nb_low <= 0._wp) then
      message = 'CSS transition state requires P_t>=0, epsilon_low>0 and n_low>0'
      return
    endif
    if (delta_epsilon < 0._wp) then
      message = 'CSS energy-density jump must be nonnegative'
      return
    endif
    if (cs2 <= 0._wp .or. cs2 > 1._wp) then
      message = 'CSS requires 0<cs2<=1 in units of c^2'
      return
    endif

    mu_t = (epsilon_low + p_transition) / nb_low
    epsilon_high = epsilon_low + delta_epsilon
    nb_high = (epsilon_high + p_transition) / mu_t
    if (.not. all(ieee_is_finite([mu_t, epsilon_high, nb_high]))) then
      message = 'CSS transition state overflowed'
      return
    endif
    if (mu_t <= 0._wp .or. epsilon_high <= 0._wp .or. nb_high <= 0._wp) then
      message = 'CSS transition state is not positive'
      return
    endif

    jump_residual = (epsilon_high - epsilon_low) - mu_t * (nb_high - nb_low)
    scale = max(1._wp, abs(delta_epsilon), abs(mu_t * (nb_high - nb_low)))
    if (abs(jump_residual) > 64._wp * epsilon(1._wp) * scale) then
      message = 'CSS first-order jump identity failed'
      return
    endif

    segment%p_transition = p_transition
    segment%epsilon_low = epsilon_low
    segment%nb_low = nb_low
    segment%delta_epsilon = delta_epsilon
    segment%cs2 = cs2
    segment%mu_transition = mu_t
    segment%epsilon_high_transition = epsilon_high
    segment%nb_high_transition = nb_high
    if (delta_epsilon > 0._wp) then
      segment%interface_class = css_interface_first_order
    else
      ! Continuity alone cannot distinguish a second-order transition from a crossover.
      segment%interface_class = css_interface_continuous_unclassified
    endif
    segment%initialized = .true.
    ok = .true.
    message = 'ok'
  end subroutine initialize_css_segment

  subroutine evaluate_css_segment(segment, pressure, epsilon, nb, mu_b, ok, message)
    type(css_segment_t), intent(in) :: segment
    real(wp), intent(in) :: pressure(:)
    real(wp), intent(out) :: epsilon(:), nb(:), mu_b(:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(wp) :: factor(size(pressure))
    integer :: n

    ok = .false.
    message = 'CSS segment is not initialized'
    n = size(pressure)
    if (.not. segment%initialized) return
    if (size(epsilon) /= n .or. size(nb) /= n .or. size(mu_b) /= n) then
      message = 'CSS pressure and output array shapes differ'
      return
    endif
    if (n < 1) then
      message = 'CSS pressure grid is empty'
      return
    endif
    if (.not. all(ieee_is_finite(pressure))) then
      message = 'CSS pressure grid contains nonfinite values'
      return
    endif
    if (any(pressure < segment%p_transition)) then
      message = 'CSS high-phase evaluation requested below P_t'
      return
    endif
    if (n > 1) then
      if (any(pressure(2:) <= pressure(:n-1))) then
        message = 'CSS high-phase pressure grid must be strictly increasing'
        return
      endif
    endif

    factor = 1._wp + (1._wp + segment%cs2) * &
      (pressure - segment%p_transition) / &
      (segment%cs2 * (segment%epsilon_high_transition + segment%p_transition))
    if (any(factor <= 0._wp) .or. .not. all(ieee_is_finite(factor))) then
      message = 'CSS thermodynamic density factor is nonpositive or nonfinite'
      return
    endif

    epsilon = segment%epsilon_high_transition + &
              (pressure - segment%p_transition) / segment%cs2
    nb = segment%nb_high_transition * factor ** (1._wp / (1._wp + segment%cs2))
    mu_b = (epsilon + pressure) / nb
    if (.not. all(ieee_is_finite(epsilon)) .or. &
        .not. all(ieee_is_finite(nb)) .or. &
        .not. all(ieee_is_finite(mu_b))) then
      message = 'CSS evaluation produced nonfinite thermodynamics'
      return
    endif
    if (any(epsilon <= 0._wp) .or. any(nb <= 0._wp) .or. any(mu_b <= 0._wp)) then
      message = 'CSS evaluation produced nonpositive thermodynamics'
      return
    endif

    ok = .true.
    message = 'ok'
  end subroutine evaluate_css_segment

end module compose_eos_css
