! E390-06 (M6): native TOV + l=2 tidal Love solver with sharp Maxwell
! phase-transition support, ported from scripts/m6/m6_tov_love_solver.py.
!
! Physics (Takatsy-Kovacs 2020 PRD 102 028501, all equations literature-checked
! in the m6 Python skeleton):
!   1. TOV equations in G=c=1 geometric units (Eq. 7-8):
!        m'(r) = 4 pi r^2 epsilon(r)
!        p'(r) = -(eps + p)(m + 4 pi r^3 p)/(r^2 - 2 m r)
!   2. l=2 tidal Love y(r) ODE (Eq. 2):
!        r y'(r) + y^2 + r^2 Q(r) + y e^lambda [1 + 4 pi r^2 (p - eps)] = 0
!      with Q(r) (Eq. 3) and metric e^lambda = (1 - 2m/r)^-1 (Eq. 4),
!      nu'(r) = 2(m + 4 pi r^3 p)/(r^2 - 2mr) (Eq. 5).
!   3. Phase-transition tidal y-jump (Eq. 11, pressure-corrected Postnikov):
!        Delta y = -4 pi r_d^3 Delta eps / (m_d + 4 pi r_d^3 p_d)
!   4. Love number k2 (Thorne 1971 / Eq. 1) and Lambda = (2/3) k2 / beta^5.
!
! Units: geometric (km, km^-2) inside the integrator; the EOS barotrope is
! supplied as P(eps) tables in MeV/fm^3 and converted. The module supports a
! smooth single-phase barotrope and a two-branch sharp Maxwell barotrope.
module compose_eos_tov
  use compose_eos_types, only: wp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: tov_love_result_t
  public :: solve_tov_love_smooth
  public :: solve_tov_love_maxwell

  ! Geometric unit conversion: MeV/fm^3 -> km^-2 in G=c=1.
  ! _MEVFM3_TO_KM2 = G * MeV / (1e-45 fm^3->m^3 * c^4) * 1e6 (m^2->km^2)
  ! G=6.67430e-11, MeV=1.602176634e-13 J, c=299792458, fm^3=1e-45 m^3
  real(wp), parameter :: mevfm3_to_km2 = 1.3243117038949245e-6_wp
  real(wp), parameter :: msun_geom_km = 1.47662504_wp  ! G*Msun/c^2 in km
  real(wp), parameter :: pi = acos(-1._wp)
  real(wp), parameter :: log_min = -700._wp   ! safe lower bound for log()

  type :: tov_love_result_t
    real(wp) :: r_km = 0._wp          ! stellar radius
    real(wp) :: m_msun = 0._wp        ! gravitational mass
    real(wp) :: compactness = 0._wp   ! M/R (dimensionless, geometric)
    real(wp) :: k2 = 0._wp            ! l=2 Love number
    real(wp) :: lambda_tidal = 0._wp  ! dimensionless tidal deformability
    real(wp) :: y_surface = 0._wp     ! y(R)
    real(wp) :: central_pressure = 0._wp
    real(wp) :: central_eps = 0._wp
    logical  :: converged = .false.
    integer  :: n_steps = 0
    character(len=64) :: status = ''
    ! sharp-Maxwell only: transition diagnostics
    logical  :: has_transition = .false.
    real(wp) :: r_transition_km = 0._wp
    real(wp) :: m_transition_msun = 0._wp
    real(wp) :: p_transition_mevfm3 = 0._wp
    real(wp) :: delta_eps_mevfm3 = 0._wp
    real(wp) :: delta_y = 0._wp
  end type tov_love_result_t

contains

  !------------------------------------------------------------------
  ! Monotone log-log linear interpolation of eps(P) and d(eps)/dP for one
  ! phase branch. P must be strictly increasing. Returns eps at P_mev and
  ! the local derivative; out-of-range queries return the endpoint.
  subroutine branch_interp(P_arr, eps_arr, n, P_mev, eps_out, dedp_out, ok)
    integer, intent(in) :: n
    real(wp), intent(in) :: P_arr(n), eps_arr(n)
    real(wp), intent(in) :: P_mev
    real(wp), intent(out) :: eps_out, dedp_out
    logical, intent(out) :: ok
    integer :: i, j
    real(wp) :: lp, le_l, le_r, lr_l, lr_r, t, le, dle_dlr

    ok = .false.; eps_out = 0._wp; dedp_out = 0._wp
    if (n < 2) return
    if (P_mev <= P_arr(1)) then
      eps_out = eps_arr(1); dedp_out = (eps_arr(2) - eps_arr(1))/(P_arr(2) - P_arr(1))
      ok = .true.; return
    endif
    if (P_mev >= P_arr(n)) then
      eps_out = eps_arr(n); dedp_out = (eps_arr(n) - eps_arr(n-1))/(P_arr(n) - P_arr(n-1))
      ok = .true.; return
    endif
    ! find bracket
    j = 1
    do i = 1, n-1
      if (P_arr(i) <= P_mev .and. P_mev < P_arr(i+1)) then
        j = i; exit
      endif
    end do
    lp = log(P_mev)
    le_l = log(eps_arr(j)); le_r = log(eps_arr(j+1))
    lr_l = log(P_arr(j));   lr_r = log(P_arr(j+1))
    t = (lp - lr_l)/(lr_r - lr_l)
    le = le_l + t*(le_r - le_l)
    eps_out = exp(le)
    dle_dlr = (le_r - le_l)/(lr_r - lr_l)
    dedp_out = exp(le)/P_mev * dle_dlr
    ok = .true.
  end subroutine branch_interp

  !------------------------------------------------------------------
  ! TOV + y RHS at radius r given (p_geom, m, y) and a branch interpolator.
  ! Returns (dpdr, dmdr, dydr) in geometric units. Sets ok=.false. on the
  ! trapped-surface or surface-cutoff condition.
  subroutine tov_rhs(r, p_geom, m, y, P_arr, eps_arr, n, &
                     surface_cutoff, dpdr, dmdr, dydr, ok)
    real(wp), intent(in) :: r, p_geom, m, y
    integer, intent(in) :: n
    real(wp), intent(in) :: P_arr(n), eps_arr(n)
    real(wp), intent(in) :: surface_cutoff
    real(wp), intent(out) :: dpdr, dmdr, dydr
    logical, intent(out) :: ok
    real(wp) :: p_mev, eps_v, dedp_v, eg, pg, fac, elam, nu_prime, cs2, Q

    ok = .false.; dpdr = 0._wp; dmdr = 0._wp; dydr = 0._wp
    if (r <= 0._wp) return
    if (p_geom <= 0._wp) then
      ok = .true.; return   ! surface reached
    endif
    if (surface_cutoff > 0._wp .and. p_geom <= surface_cutoff) then
      ok = .true.; return   ! freeze at surface cutoff
    endif
    p_mev = p_geom/mevfm3_to_km2
    call branch_interp(P_arr, eps_arr, n, p_mev, eps_v, dedp_v, ok)
    if (.not. ok .or. dedp_v <= 0._wp) return
    eg = eps_v*mevfm3_to_km2
    pg = p_geom
    fac = 1._wp - 2._wp*m/r
    if (fac <= 0._wp) return   ! trapped surface
    elam = 1._wp/fac
    nu_prime = 2._wp*(m + 4._wp*pi*r**3*pg)/(r*(r - 2._wp*m))
    cs2 = 1._wp/dedp_v
    Q = (4._wp*pi*elam*(5._wp*eg + 9._wp*pg + (eg + pg)*cs2) &
         - 6._wp*elam/(r*r) - nu_prime**2)
    dpdr = -(eg + pg)*(m + 4._wp*pi*r**3*pg)/(r*(r - 2._wp*m))
    dmdr = 4._wp*pi*r*r*eg
    dydr = -(y*y + y*elam*(1._wp + 4._wp*pi*r**2*(pg - eg)) + r*r*Q)/r
    ok = ieee_is_finite(dpdr + dmdr + dydr)
  end subroutine tov_rhs

  !------------------------------------------------------------------
  ! Central Taylor expansion for the initial condition at small r0.
  ! state = (p_geom, m, y). (Takatsy-Kovacs 2020 near-center expansion.)
  subroutine central_expansion(pc_mev, ec_mev, cs2c, r0, p0, m0, y0)
    real(wp), intent(in) :: pc_mev, ec_mev, cs2c
    real(wp), intent(out) :: r0, p0, m0, y0
    real(wp) :: pcg, ecg
    r0 = 1.0e-4_wp   ! km
    pcg = pc_mev*mevfm3_to_km2
    ecg = ec_mev*mevfm3_to_km2
    p0 = pcg - (2._wp*pi/3._wp)*(ecg + pcg)*(ecg + 3._wp*pcg)*r0*r0
    m0 = (4._wp*pi/3._wp)*ecg*r0**3
    y0 = 2._wp - (4._wp*pi/7._wp)*(5._wp*ecg + 9._wp*pcg + (ecg + pcg)/cs2c)*r0*r0
  end subroutine central_expansion

  !------------------------------------------------------------------
  ! Love number k2 and Lambda from compactness beta and surface y.
  ! (Thorne 1971 / Takatsy-Kovacs 2020 Eq. 1.)
  subroutine love_number_k2(beta, yR, k2, lambda_tidal, ok)
    real(wp), intent(in) :: beta, yR
    real(wp), intent(out) :: k2, lambda_tidal
    logical, intent(out) :: ok
    real(wp) :: b, t1, t2, t3, den
    ok = .false.; k2 = 0._wp; lambda_tidal = 0._wp
    if (.not. (beta > 0._wp .and. beta < 0.5_wp)) return
    b = beta
    t1 = (1._wp - 2._wp*b)**2 * b**5 * (2._wp*b*(yR - 1._wp) - yR + 2._wp)
    t2 = 2._wp*b*(4._wp*(yR + 1._wp)*b**4 + (6._wp*yR - 4._wp)*b**3 &
                 + (26._wp - 22._wp*yR)*b**2 + 3._wp*(5._wp*yR - 8._wp)*b - 3._wp*yR + 6._wp)
    t3 = 3._wp*(1._wp - 2._wp*b)**2 * (2._wp*b*(yR - 1._wp) - yR + 2._wp)*log(1._wp - 2._wp*b)
    den = t2 + t3
    if (abs(den) < 1.0e-30_wp) return
    k2 = (t1/2._wp)/den
    lambda_tidal = (2._wp/3._wp)*k2/b**5
    ok = ieee_is_finite(k2 + lambda_tidal)
  end subroutine love_number_k2

  !------------------------------------------------------------------
  ! RK45 (Dormand-Prince) adaptive integrator for the TOV+y system on one
  ! branch. Integrates from r0 outward until p<=surface_cutoff or the
  ! trapped-surface condition. Returns the final (r, p, m, y) and n_steps.
  subroutine integrate_branch(r0, p0, m0, y0, P_arr, eps_arr, n, &
                              dr_max, p_surface, r_out, m_out, y_out, n_steps, ok, status)
    real(wp), intent(in) :: r0, p0, m0, y0
    integer, intent(in) :: n
    real(wp), intent(in) :: P_arr(n), eps_arr(n)
    real(wp), intent(in) :: dr_max, p_surface
    real(wp), intent(out) :: r_out, m_out, y_out
    integer, intent(out) :: n_steps
    logical, intent(out) :: ok
    character(len=*), intent(out) :: status
    ! Dormand-Prince RK45 coefficients
    real(wp), parameter :: a21=1._wp/5._wp
    real(wp), parameter :: a31=3._wp/40._wp, a32=9._wp/40._wp
    real(wp), parameter :: a41=44._wp/45._wp, a42=-56._wp/15._wp, a43=32._wp/9._wp
    real(wp), parameter :: a51=19372._wp/6561._wp, a52=-25360._wp/2187._wp, &
                          a53=64448._wp/6561._wp, a54=-212._wp/729._wp
    real(wp), parameter :: a61=9017._wp/3168._wp, a62=-355._wp/33._wp, &
                          a63=46732._wp/5247._wp, a54_=49._wp/176._wp, a65=-5103._wp/18656._wp
    real(wp), parameter :: b1=35._wp/384._wp, b3=500._wp/1113._wp, b4=125._wp/192._wp, &
                          b5=-2187._wp/6784._wp, b6=11._wp/84._wp
    real(wp), parameter :: e1=71._wp/57600._wp, e3=-71._wp/16695._wp, e4=71._wp/1920._wp, &
                          e5=-17253._wp/339200._wp, e6=22._wp/525._wp
    real(wp) :: r, p, m, y, h, dp, dm, dy, r2, p2, m2, y2
    real(wp) :: k1p, k1m, k1y, k2p, k2m, k2y, k3p, k3m, k3y
    real(wp) :: k4p, k4m, k4y, k5p, k5m, k5y, k6p, k6m, k6y
    real(wp) :: err, tol, h_new
    logical :: step_ok
    integer :: i
    integer, parameter :: max_steps = 200000

    ok = .false.; status = ''
    r = r0; p = p0; m = m0; y = y0; n_steps = 0
    tol = 1.0e-6_wp
    h = min(dr_max, 5.0e-3_wp)
    do i = 1, max_steps
      ! k1
      call tov_rhs(r, p, m, y, P_arr, eps_arr, n, p_surface, k1p, k1m, k1y, step_ok)
      if (.not. step_ok) exit
      ! k2
      call tov_rhs(r + 0.5_wp*h, p + 0.5_wp*h*k1p, m + 0.5_wp*h*k1m, &
                   y + 0.5_wp*h*k1y, P_arr, eps_arr, n, p_surface, k2p, k2m, k2y, step_ok)
      if (.not. step_ok) exit
      ! k3
      call tov_rhs(r + 0.5_wp*h, p + h*(a31*k1p + a32*k2p), m + h*(a31*k1m + a32*k2m), &
                   y + h*(a31*k1y + a32*k2y), P_arr, eps_arr, n, p_surface, k3p, k3m, k3y, step_ok)
      if (.not. step_ok) exit
      ! k4
      call tov_rhs(r + h, p + h*(a41*k1p + a42*k2p + a43*k3p), &
                   m + h*(a41*k1m + a42*k2m + a43*k3m), &
                   y + h*(a41*k1y + a42*k2y + a43*k3y), P_arr, eps_arr, n, p_surface, &
                   k4p, k4m, k4y, step_ok)
      if (.not. step_ok) exit
      ! k5
      call tov_rhs(r + h, p + h*(a51*k1p + a52*k2p + a53*k3p + a54*k4p), &
                   m + h*(a51*k1m + a52*k2m + a53*k3m + a54*k4m), &
                   y + h*(a51*k1y + a52*k2y + a53*k3y + a54*k4y), &
                   P_arr, eps_arr, n, p_surface, k5p, k5m, k5y, step_ok)
      if (.not. step_ok) exit
      ! k6
      call tov_rhs(r + h, p + h*(a61*k1p + a62*k2p + a63*k3p + a54_*k4p + a65*k5p), &
                   m + h*(a61*k1m + a62*k2m + a63*k3m + a54_*k4m + a65*k5m), &
                   y + h*(a61*k1y + a62*k2y + a63*k3y + a54_*k4y + a65*k5y), &
                   P_arr, eps_arr, n, p_surface, k6p, k6m, k6y, step_ok)
      if (.not. step_ok) exit
      ! 5th-order step
      r2 = r + h
      p2 = p + h*(b1*k1p + b3*k3p + b4*k4p + b5*k5p + b6*k6p)
      m2 = m + h*(b1*k1m + b3*k3m + b4*k4m + b5*k5m + b6*k6m)
      y2 = y + h*(b1*k1y + b3*k3y + b4*k4y + b5*k5y + b6*k6y)
      ! embedded error estimate (4th vs 5th order)
      err = abs(h*(e1*k1p + e3*k3p + e4*k4p + e5*k5p + e6*k6p))
      err = err + abs(h*(e1*k1m + e3*k3m + e4*k4m + e5*k5m + e6*k6m))
      err = err + abs(h*(e1*k1y + e3*k3y + e4*k4y + e5*k5y + e6*k6y))
      if (err < 1.0e-30_wp) err = 1.0e-30_wp
      ! accept the step
      r = r2; p = p2; m = m2; y = y2
      n_steps = i
      ! adapt: increase or decrease h within [dr_max/1000, dr_max]
      h_new = h*0.9_wp*(tol/err)**0.2_wp
      h = min(dr_max, max(h_new, dr_max/1000._wp))
      ! surface reached?
      if (p <= p_surface .or. p <= 0._wp) then
        ok = .true.; exit
      endif
      ! trapped surface?
      if (1._wp - 2._wp*m/r <= 0._wp) then
        status = 'trapped surface'; ok = .false.; exit
      endif
    end do
    r_out = r; m_out = m; y_out = y
    if (i >= max_steps .and. .not. ok) status = 'max steps reached'
    if (ok .and. len_trim(status) == 0) status = 'ok'
  end subroutine integrate_branch

  !------------------------------------------------------------------
  ! Solve TOV + l=2 Love on a smooth (single-phase) barotrope.
  ! P_arr, eps_arr are MeV/fm^3 (strictly increasing P). central_pressure
  ! in MeV/fm^3. Returns R(km), M(Msun), C, k2, Lambda.
  subroutine solve_tov_love_smooth(n, P_arr, eps_arr, central_pressure, result)
    integer, intent(in) :: n
    real(wp), intent(in) :: P_arr(n), eps_arr(n)
    real(wp), intent(in) :: central_pressure
    type(tov_love_result_t), intent(out) :: result
    real(wp) :: pc, ec, dedp_c, cs2c, r0, p0, m0, y0
    real(wp) :: p_surface, r_out, m_out, y_out, beta
    integer :: n_steps
    logical :: ok
    character(len=64) :: status

    result = tov_love_result_t()
    result%central_pressure = central_pressure
    pc = central_pressure
    if (pc <= P_arr(1) .or. pc > P_arr(n)) then
      result%status = 'central pressure outside EOS range'; return
    endif
    call branch_interp(P_arr, eps_arr, n, pc, ec, dedp_c, ok)
    if (.not. ok .or. dedp_c <= 0._wp) then
      result%status = 'EOS interp failed at center'; return
    endif
    cs2c = 1._wp/dedp_c
    if (.not. (cs2c > 0._wp .and. cs2c <= 1._wp + 1.0e-5_wp)) then
      result%status = 'central cs2 unphysical'; return
    endif
    result%central_eps = ec
    call central_expansion(pc, ec, cs2c, r0, p0, m0, y0)
    p_surface = P_arr(1)*mevfm3_to_km2   ! lowest tabulated pressure (geometric)
    call integrate_branch(r0, p0, m0, y0, P_arr, eps_arr, n, 0.12_wp, p_surface, &
                          r_out, m_out, y_out, n_steps, ok, status)
    result%n_steps = n_steps
    result%status = trim(status)
    if (.not. ok) return
    result%r_km = r_out
    result%m_msun = m_out/msun_geom_km
    if (r_out > 0._wp) then
      beta = m_out/r_out
      result%compactness = beta
      call love_number_k2(beta, y_out, result%k2, result%lambda_tidal, ok)
      if (.not. ok) then
        result%status = 'Love number singular'; return
      endif
    endif
    result%y_surface = y_out
    result%converged = .true.
    result%status = 'ok'
  end subroutine solve_tov_love_smooth

  !------------------------------------------------------------------
  ! Solve TOV + l=2 Love on a two-branch sharp Maxwell barotrope.
  ! high_P/high_eps: the high-density branch (P > P_transition).
  ! low_P/low_eps:   the low-density branch (P <= P_transition).
  ! P_transition is the transition pressure; delta_eps is the energy-density
  ! jump at constant P_transition. Integrate high branch from center to
  ! P_transition, apply the Takatsy-Kovacs y-jump, integrate low branch to P=0.
  subroutine solve_tov_love_maxwell(nh, high_P, high_eps, nl, low_P, low_eps, &
                                     P_transition, delta_eps, central_pressure, result)
    integer, intent(in) :: nh, nl
    real(wp), intent(in) :: high_P(nh), high_eps(nh)
    real(wp), intent(in) :: low_P(nl), low_eps(nl)
    real(wp), intent(in) :: P_transition, delta_eps, central_pressure
    type(tov_love_result_t), intent(out) :: result
    real(wp) :: pc, ec, dedp_c, cs2c, r0, p0, m0, y0
    real(wp) :: p_surface, r_t, m_t, y_t, r_out, m_out, y_out, beta
    real(wp) :: p_trans_geom, y_after_jump, m_t_msun
    integer :: n_steps1, n_steps2
    logical :: ok
    character(len=64) :: status

    result = tov_love_result_t()
    result%central_pressure = central_pressure
    result%has_transition = .true.
    result%p_transition_mevfm3 = P_transition
    result%delta_eps_mevfm3 = delta_eps
    pc = central_pressure
    if (pc <= P_transition) then
      result%status = 'central pressure below transition'; return
    endif
    if (P_transition < low_P(1) .or. P_transition > high_P(nh)) then
      result%status = 'transition pressure outside branches'; return
    endif
    call branch_interp(high_P, high_eps, nh, pc, ec, dedp_c, ok)
    if (.not. ok .or. dedp_c <= 0._wp) then
      result%status = 'high-branch interp failed at center'; return
    endif
    cs2c = 1._wp/dedp_c
    if (.not. (cs2c > 0._wp .and. cs2c <= 1._wp + 1.0e-5_wp)) then
      result%status = 'central cs2 unphysical'; return
    endif
    result%central_eps = ec
    call central_expansion(pc, ec, cs2c, r0, p0, m0, y0)
    p_trans_geom = P_transition*mevfm3_to_km2
    ! integrate high branch from center to P_transition
    call integrate_branch(r0, p0, m0, y0, high_P, high_eps, nh, 0.12_wp, p_trans_geom, &
                          r_t, m_t, y_t, n_steps1, ok, status)
    result%n_steps = n_steps1
    if (.not. ok) then
      result%status = 'high-branch integration: '//trim(status); return
    endif
    result%r_transition_km = r_t
    result%m_transition_msun = m_t/msun_geom_km
    m_t_msun = m_t/msun_geom_km
    ! Takatsy-Kovacs y-jump (Eq. 11, pressure-corrected):
    !   Delta y = -4 pi r_d^3 delta_eps / (m_d + 4 pi r_d^3 p_d)
    ! delta_eps in MeV/fm^3 -> geometric km^-2 for the jump.
    result%delta_y = -4._wp*pi*r_t**3*(delta_eps*mevfm3_to_km2) / &
                      (m_t + 4._wp*pi*r_t**3*p_trans_geom)
    y_after_jump = y_t + result%delta_y
    ! integrate low branch from r_t with the jumped y to P=0
    ! the low branch starts at P_transition; eps jumps down by delta_eps at
    ! the same pressure, so the low-branch eps at P_transition is
    ! high_eps(P_transition) - delta_eps.
    p_surface = low_P(1)*mevfm3_to_km2
    call integrate_branch(r_t, p_trans_geom, m_t, y_after_jump, &
                          low_P, low_eps, nl, 0.12_wp, p_surface, &
                          r_out, m_out, y_out, n_steps2, ok, status)
    result%n_steps = result%n_steps + n_steps2
    if (.not. ok) then
      result%status = 'low-branch integration: '//trim(status); return
    endif
    result%r_km = r_out
    result%m_msun = m_out/msun_geom_km
    if (r_out > 0._wp) then
      beta = m_out/r_out
      result%compactness = beta
      call love_number_k2(beta, y_out, result%k2, result%lambda_tidal, ok)
      if (.not. ok) then
        result%status = 'Love number singular'; return
      endif
    endif
    result%y_surface = y_out
    result%converged = .true.
    result%status = 'ok'
  end subroutine solve_tov_love_maxwell

end module compose_eos_tov
