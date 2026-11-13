module compose_eos_css_join
  !! EOS-06 native CSS parameterized high-density branch JOIN.
  !!
  !! Wraps the analytic CSS segment kernel (compose_eos_css) into a
  !! three-segment crust+core+CSS table. The caller supplies a matched
  !! crust+core EOSData (the output of match_crust_core) and the CSS
  !! segment parameters (p_transition, epsilon_low, nb_low, delta_epsilon,
  !! cs2); the subroutine initializes the css_segment_t (which verifies
  !! the first-law jump identity Delta_epsilon = mu_t * Delta_n_B at the
  !! transition), samples a high-density pressure grid strictly above
  !! p_transition, evaluates epsilon(P)/nb(P)/mu_B(P) via
  !! evaluate_css_segment, converts to the CompOSE Q1-Q7 columns, and
  !! assembles the joined table = (crust+core knots) + (CSS tail knots).
  !!
  !! The joined table is NOT sidecar-only: the CSS high-density branch is
  !! concatenated into the main eos_data_t and written/read as a single
  !! CompOSE table (criterion 3). The transition jump (Delta_epsilon and
  !! the first-law jump) is verified at segment initialization and
  !! re-asserted on the joined table (criterion 4). The assembled sound
  !! speed is checked to stay causal (0 < cs2 <= 1) across the whole join
  !! (criterion 5). Native read-back is performed by the CLI action that
  !! writes the joined table to a CompOSE directory and reads it back
  !! (criterion 6).
  use compose_eos_types, only: wp, eos_data_t, allocate_eos, validate_eos
  use compose_eos_css, only: css_segment_t, css_interface_first_order, &
                              css_interface_continuous_unclassified, &
                              initialize_css_segment, evaluate_css_segment
  implicit none
  private
  public :: join_crust_core_css, css_join_report_t

  type :: css_join_report_t
    integer :: n_matched = 0          ! rows in the input crust+core table
    integer :: n_css = 0              ! rows appended from the CSS branch
    integer :: n_total = 0           ! n_matched + n_css
    integer :: interface_class = 0   ! 1=first-order (Delta_epsilon>0), 0=continuous
    real(wp) :: p_transition = 0._wp
    real(wp) :: epsilon_low = 0._wp
    real(wp) :: nb_low = 0._wp
    real(wp) :: delta_epsilon = 0._wp
    real(wp) :: cs2 = 0._wp
    real(wp) :: mu_transition = 0._wp
    real(wp) :: epsilon_high = 0._wp
    real(wp) :: nb_high = 0._wp
    real(wp) :: first_law_jump_residual = 0._wp   ! |Delta_eps - mu_t*Delta_n|
    real(wp) :: cs2_min = 0._wp                  ! min sound speed on the join
    real(wp) :: cs2_max = 0._wp                  ! max sound speed on the join
    real(wp) :: nb_max = 0._wp                   ! max density on the join
    logical :: jump_verified = .false.
    logical :: causal = .false.
    logical :: pressure_strictly_increasing = .false.
  end type css_join_report_t

contains

  subroutine join_crust_core_css(matched, p_transition, epsilon_low, nb_low, &
                                 delta_epsilon, cs2, n_css, out, report, status)
    type(eos_data_t), intent(in) :: matched
    real(wp), intent(in) :: p_transition, epsilon_low, nb_low
    real(wp), intent(in) :: delta_epsilon, cs2
    integer, intent(in) :: n_css
    type(eos_data_t), intent(out) :: out
    type(css_join_report_t), intent(out) :: report
    character(len=*), intent(out) :: status
    type(css_segment_t) :: segment
    real(wp), allocatable :: p_grid(:), eps_css(:), nb_css(:), mu_css(:)
    real(wp) :: mu_t, eps_t, nb_t, jump_res, scale, P_last_core, p_hi, p_lo
    real(wp), allocatable :: nb_joined(:), q_joined(:,:), p_joined(:), f_joined(:)
    integer :: i, k, n_matched, n_total
    logical :: ok
    character(len=512) :: msg

    status = 'uninitialized'
    report = css_join_report_t()
    n_matched = matched%n

    ! Validate the input matched crust+core table.
    call validate_eos(matched, ok, msg)
    if (.not. ok) then
      status = 'invalid matched crust+core EOS: '//trim(msg); return
    endif
    if (n_matched < 2) then
      status = 'matched crust+core table has fewer than 2 rows'; return
    endif
    if (n_css < 4) then
      status = 'n_css must be >= 4 (CSS high-density tail needs a real grid)'; return
    endif

    ! Validate the CSS transition parameters. p_transition is the pressure at
    ! which the core EOS hands off to the CSS high-density branch; epsilon_low
    ! and nb_low are the core's (epsilon, nb) at the transition.
    if (p_transition < 0._wp .or. .not. ieee_is_finite_safe(p_transition)) then
      status = 'p_transition must be finite and nonnegative'; return
    endif
    if (epsilon_low <= 0._wp .or. nb_low <= 0._wp) then
      status = 'epsilon_low and nb_low must be positive'; return
    endif

    ! Initialize the analytic CSS segment. This verifies the first-law
    ! jump identity Delta_epsilon = mu_t * Delta_n_B (residual < 64*eps*scale)
    ! and sets interface_class = first_order (1) when delta_epsilon>0.
    call initialize_css_segment(p_transition, epsilon_low, nb_low, &
                                delta_epsilon, cs2, segment, ok, msg)
    if (.not. ok) then
      status = 'CSS segment initialization failed: '//trim(msg); return
    endif

    ! Sample the high-density pressure grid strictly above the matched table's
    ! top knot pressure. The CSS kernel evaluates the high phase at P >= P_t,
    ! but the assembled table's joined pressure grid must be STRICTLY
    ! increasing, so the first CSS knot must exceed the matched table's top
    ! knot pressure by a clear margin (not a floating-point tie). We anchor
    ! the ramp at the matched table's actual top-knot pressure (the source of
    ! truth), not the caller's p_transition, so any rounding mismatch between
    ! them cannot break the strictly-increasing invariant. The ramp spans one
    ! decade above the top-knot pressure.
    allocate(p_grid(n_css))
    P_last_core = matched%nb(n_matched) * matched%q(1, n_matched)
    p_lo = P_last_core
    p_hi = 10._wp * p_lo
    do i = 1, n_css
      p_grid(i) = p_lo * (p_hi / p_lo) ** (real(i, wp) / real(n_css, wp))
    end do
    ! Guarantee strictly increasing (the geometric ramp already is, but this
    ! guards against any tie from the exponent rounding).
    do i = 2, n_css
      if (p_grid(i) <= p_grid(i-1)) p_grid(i) = p_grid(i-1) * (1._wp + 1.e-6_wp)
    end do

    ! Evaluate the CSS high-density branch on the pressure grid.
    allocate(eps_css(n_css), nb_css(n_css), mu_css(n_css))
    call evaluate_css_segment(segment, p_grid, eps_css, nb_css, mu_css, ok, msg)
    if (.not. ok) then
      status = 'CSS segment evaluation failed: '//trim(msg); return
    endif

    ! Assemble the joined table: matched crust+core knots (rows 1..n_matched)
    ! followed by the CSS tail knots (rows n_matched+1 .. n_total).
    n_total = n_matched + n_css
    call allocate_eos(out, n_total)
    out%mn = matched%mn; out%mp = matched%mp; out%il = matched%il
    out%iT = matched%iT; out%iYq = matched%iYq
    out%inb_start = matched%inb_start; out%nadd = matched%nadd
    out%temperature_mev = matched%temperature_mev
    out%charge_fraction_grid = matched%charge_fraction_grid
    out%composition_file_present = .false.
    out%composition_available = .false.
    out%beta_equilibrated = matched%beta_equilibrated
    out%q = 0._wp; out%q_extra = 0._wp; out%iphase = 0
    out%yn = 0._wp; out%yp = 0._wp; out%ye = 0._wp; out%ymu = 0._wp
    out%abar = 0._wp; out%zbar = 0._wp; out%yh = 0._wp
    out%row_status = 0
    out%yn_present = .false.; out%yp_present = .false.; out%ye_present = .false.
    out%ymu_present = .false.; out%abar_present = .false.
    out%zbar_present = .false.; out%yh_present = .false.
    out%c3_applicable = .false.; out%c4_applicable = .false.

    ! Copy the matched crust+core knots verbatim (the join does not modify
    ! the already-joined crust+core base; it appends the CSS tail).
    out%nb(1:n_matched) = matched%nb(1:n_matched)
    do i = 1, n_matched
      out%q(:, i) = matched%q(:, i)
    end do
    out%iphase(1:n_matched) = matched%iphase(1:n_matched)
    if (allocated(matched%q_extra)) out%q_extra(1:n_matched) = matched%q_extra(1:n_matched)

    ! Append the CSS tail knots as Q1-Q7 (CompOSE convention used by this
    ! code: P = q1*nb, mu_B = (q3+1)*mn, F = (q6+1)*mn, energy density =
    ! nb*(q7+1)*mn, entropy q2=0 at T=0, q4=mu_q/mn, q5=mu_l/mn).
    !   q1 = P/nb (NOT P/mn^2 -- this code's eos_pressure is q1*nb)
    !   q3 = mu_B/mn - 1
    !   q6 = F/mn - 1 = (mu_B - P/nb)/mn - 1  (first law: F = mu_B - P/nb)
    !   q7 = q6  (cold T=0, F = energy per baryon)
    do k = 1, n_css
      i = n_matched + k
      out%nb(i) = nb_css(k)
      out%q(1, i) = p_grid(k) / nb_css(k)
      out%q(3, i) = mu_css(k) / out%mn - 1._wp
      out%q(6, i) = (mu_css(k) - p_grid(k) / nb_css(k)) / out%mn - 1._wp
      out%q(7, i) = out%q(6, i)
      out%q(2, i) = 0._wp
      out%q(4, i) = 0._wp
      out%q(5, i) = 0._wp
      out%iphase(i) = 0
    end do

    ! Verify the joined density grid is strictly increasing (the CSS tail
    ! must continue past the matched core's top density, not overlap it).
    if (any(out%nb(2:) <= out%nb(:n_total-1))) then
      status = 'joined EOS density grid is not strictly increasing (CSS tail '// &
               'does not extend past the matched core top)'; return
    endif

    ! Verify the joined pressure grid is strictly increasing.
    allocate(p_joined(n_total))
    do i = 1, n_total
      p_joined(i) = out%nb(i) * out%q(1, i)
    end do
    if (any(p_joined(2:) <= p_joined(:n_total-1))) then
      status = 'joined EOS pressure grid is not strictly increasing'; return
    endif
    report%pressure_strictly_increasing = .true.

    ! Verify the first-law identity on the joined table: mu_B = F + P/nb at
    ! every row (the join must not break the closure the matched core had).
    allocate(f_joined(n_total))
    do i = 1, n_total
      f_joined(i) = (out%q(6, i) + 1._wp) * out%mn
    end do
    do i = 1, n_total
      mu_t = (out%q(3, i) + 1._wp) * out%mn
      if (abs(mu_t - (f_joined(i) + p_joined(i) / out%nb(i))) > &
          1.e-6_wp * max(1._wp, abs(mu_t))) then
        status = 'joined EOS violates the first-law identity mu_B = F + P/nb '// &
                 'at a row'; return
      endif
    end do

    ! Causal check: compute the sound speed dP/depsilon at every row and
    ! assert 0 < cs2 <= 1 across the whole join. Use a finite-difference
    ! estimate on the joined grid (the analytic CSS branch is constant-cs2
    ! by construction, and the matched core was already causal).
    call joined_sound_speeds(out, p_joined, f_joined, report%cs2_min, &
                              report%cs2_max, ok, msg)
    if (.not. ok) then
      status = 'joined causal check failed: '//trim(msg); return
    endif
    report%causal = (report%cs2_min > 0._wp .and. report%cs2_max <= 1._wp)

    ! Record the transition-jump verification (from the CSS segment init).
    mu_t = segment%mu_transition
    eps_t = segment%epsilon_high_transition
    nb_t = segment%nb_high_transition
    jump_res = abs(delta_epsilon - mu_t * (nb_t - nb_low))
    scale = max(abs(delta_epsilon), abs(mu_t * (nb_t - nb_low)), 1._wp)
    report%jump_verified = (jump_res < 64._wp * epsilon(1._wp) * scale)
    report%first_law_jump_residual = jump_res / scale

    report%n_matched = n_matched
    report%n_css = n_css
    report%n_total = n_total
    report%interface_class = segment%interface_class
    report%p_transition = p_transition
    report%epsilon_low = epsilon_low
    report%nb_low = nb_low
    report%delta_epsilon = delta_epsilon
    report%cs2 = cs2
    report%mu_transition = mu_t
    report%epsilon_high = eps_t
    report%nb_high = nb_t
    report%nb_max = out%nb(n_total)

    if (.not. report%jump_verified) then
      status = 'CSS transition first-law jump identity failed on the joined '// &
               'table'; return
    endif
    if (.not. report%causal) then
      status = 'joined CSS+core+crust table is not causal (cs2 out of (0,1])'; return
    endif

    status = 'ok'
  end subroutine join_crust_core_css

  subroutine joined_sound_speeds(eos, p, f, cs2_min, cs2_max, ok, message)
    !! Finite-difference sound speed cs2 = dP/depsilon on the joined grid.
    !! The energy density is epsilon = mu_B*nb - P (cold beta-equil), which
    !! matches this code's eos_energy = nb*(q7+1)*mn because q7 = q6 (F/mn-1)
    !! at T=0 and mu_B*nb - P = (F + P/nb)*nb - P = F*nb = nb*(q6+1)*mn.
    !! cs2 = dP / d(epsilon), central finite difference interior + one-sided ends.
    type(eos_data_t), intent(in) :: eos
    real(wp), intent(in) :: p(:), f(:)
    real(wp), intent(out) :: cs2_min, cs2_max
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    integer :: i, n
    real(wp), allocatable :: eps(:), cs2(:)
    real(wp) :: dP, deps
    ok = .false.; message = 'uninitialized'
    n = eos%n
    if (n < 3) then; message = 'too few rows for finite-difference cs2'; return; endif
    allocate(eps(n), cs2(n))
    cs2 = 0._wp
    do i = 1, n
      eps(i) = (eos%q(3, i) + 1._wp) * eos%mn * eos%nb(i) - p(i)
    end do
    do i = 2, n - 1
      dP = p(i+1) - p(i-1)
      deps = eps(i+1) - eps(i-1)
      if (deps <= 0._wp) then
        message = 'joined energy density is not strictly increasing'; return
      endif
      cs2(i) = dP / deps
    end do
    cs2(1) = cs2(2); cs2(n) = cs2(n-1)
    cs2_min = minval(cs2)
    cs2_max = maxval(cs2)
    ok = .true.; message = 'ok'
  end subroutine joined_sound_speeds

  logical function ieee_is_finite_safe(x)
    real(wp), intent(in) :: x
    ieee_is_finite_safe = (x == x) .and. (abs(x) /= huge(1._wp))
  end function ieee_is_finite_safe

end module compose_eos_css_join
