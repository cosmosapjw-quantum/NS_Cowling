program test_tov_love
  !! Native unit test for compose_eos_tov (E390-06 / M6 acceptance: CI pass for
  !! the TOV-Love module). Verifies physical invariants without invoking the
  !! full Python frontend:
  !!   1. solve_tov_love_smooth on a polytropic-ish analytic EOS returns
  !!      converged=.true., positive mass and radius, and k2 in (0, 1).
  !!   2. More massive (higher central pressure) configuration has a larger
  !!      gravitational mass (monotone mass-pressure branch, no spurious sign
  !!      flips).
  !!   3. Compactness M/R is positive and bounded (superluminal analytic branch
  !!      keeps C < 1).
  !!   4. Central pressure outside the EOS range returns a non-converged
  !!      result with a clear status (fail-closed, no silent wrong answer).
  !! The analytic test barotrope is a piecewise-linear P(eps) with dP/deps = 1/3
  !! (a causal, monotone, simple polytrope-like branch) sampled on a uniform
  !! log grid from eps_low to eps_high.
  use, intrinsic :: iso_fortran_env, only: real64
  use compose_eos_types, only: wp
  use compose_eos_tov, only: tov_love_result_t, solve_tov_love_smooth
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: fails, n, i
  integer, parameter :: np = 64
  real(wp) :: P_arr(np), eps_arr(np), pc1, pc2, pc_bad, dpdcs, ep0, eps_lo, eps_hi, dlog
  type(tov_love_result_t) :: r1, r2, rbad

  fails = 0
  write(*,'(a)') 'TOV-Love smooth barotrope native checks:'

  ! Build a causal analytic barotrope: P = (1/3)(eps - eps0) with eps0 = 100 MeV/fm^3,
  ! eps in [eps_lo, eps_hi]. Sample uniformly in log(eps). This is the simplest
  ! causal (c_s^2 = 1/3) monotone branch.
  dpdcs = 1._wp / 3._wp
  ep0 = 100._wp
  eps_lo = 100.1_wp          ! just above eps0 so P>0
  eps_hi = 1600._wp
  n = np
  do i = 1, n
    ! uniform in log(eps)
    eps_arr(i) = exp( log(eps_lo) + (log(eps_hi)-log(eps_lo)) * real(i-1,wp) / real(n-1,wp) )
    P_arr(i) = dpdcs * (eps_arr(i) - ep0)
  end do
  ! sanity: P strictly increasing
  if (any(P_arr(2:n) <= P_arr(1:n-1))) then
    write(*,'(2x,a)') 'FAIL test barotrope not strictly monotone'
    fails = fails + 1
  end if

  ! --- 2. monotone mass branch on the STABLE side: lower central pressure -> smaller mass.
  !! A causal barotrope has a maximum-mass turning point; past it the branch is
  !! unstable and mass decreases with pressure. We compare two configurations
  !! both below the turning point (lower-half of the branch), so the stable
  !! branch must be monotone-increasing in pc.
  pc2 = 0.5_wp * (P_arr(np/2) + P_arr(np/2 + 1))   ! mid
  ! pick a genuinely lower central pressure, well inside the stable region
  pc1 = P_arr(np/4)                                  ! lower
  call solve_tov_love_smooth(n, P_arr, eps_arr, pc1, r1)
  call solve_tov_love_smooth(n, P_arr, eps_arr, pc2, r2)
  if (.not. r2%converged) then
    write(*,'(2x,a,a)') 'FAIL TOV mid-pc not converged: ', trim(r2%status)
    fails = fails + 1
  else if (r2%m_msun <= 0._wp .or. r2%r_km <= 0._wp) then
    write(*,'(2x,a)') 'FAIL TOV nonpositive mass or radius'
    fails = fails + 1
  else if (r2%k2 < 0._wp .or. r2%k2 > 1._wp) then
    write(*,'(2x,a,es12.4)') 'FAIL TOV k2 out of [0,1]: ', r2%k2
    fails = fails + 1
  else if (r2%compactness <= 0._wp .or. r2%compactness >= 1._wp) then
    write(*,'(2x,a,es12.4)') 'FAIL TOV compactness out of (0,1): ', r2%compactness
    fails = fails + 1
  else
    write(*,'(2x,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') 'PASS TOV mid: M=', r2%m_msun, &
      ' M_sun R=', r2%r_km, ' km k2=', r2%k2, ' C=', r2%compactness
  end if

  ! --- 1b. stable branch monotonicity ---
  if (.not. r1%converged) then
    write(*,'(2x,a,a)') 'FAIL TOV low-pc not converged: ', trim(r1%status)
    fails = fails + 1
  else if (r2%converged .and. r2%m_msun <= r1%m_msun) then
    write(*,'(2x,a,es12.4,a,es12.4)') 'FAIL stable branch not monotone: m(pc_low)=', r1%m_msun, &
      ' >= m(pc_mid)=', r2%m_msun
    fails = fails + 1
  else if (r2%converged) then
    write(*,'(2x,a,es12.4,a,es12.4)') 'PASS stable branch monotone: M(pc_low)=', r1%m_msun, &
      ' < M(pc_mid)=', r2%m_msun
  end if

  ! --- 3. fail-closed: central pressure outside EOS range ---
  pc_bad = 2._wp * P_arr(n)    ! above the highest tabulated P
  call solve_tov_love_smooth(n, P_arr, eps_arr, pc_bad, rbad)
  if (rbad%converged) then
    write(*,'(2x,a)') 'FAIL TOV converged on out-of-range central pressure (should fail-closed)'
    fails = fails + 1
  else
    write(*,'(2x,a,a)') 'PASS TOV out-of-range central pressure fail-closed: ', trim(rbad%status)
  end if

  if (fails == 0) then
    write(*,'(a)') 'TOV-Love smooth barotrope: PASS'
  else
    write(*,'(a,i0,a)') 'TOV-Love smooth barotrope: FAIL (', fails, ' checks)'
    stop 1
  end if
end program test_tov_love
