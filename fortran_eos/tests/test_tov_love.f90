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
  use compose_eos_types, only: wp,eos_data_t,allocate_eos
  use compose_eos_tov, only: tov_love_result_t, solve_tov_love_smooth, solve_tov_love_maxwell
  use compose_eos_verify, only: verify_report_t,check_consistency
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

  ! --- 4. sharp-Maxwell TOV-Love branch-domain and jump-contract tests ---
  block
    real(wp) :: high_P(4),high_eps(4),low_P(4),low_eps(4),high_bad(4),low_bad(4)
    type(tov_love_result_t) :: rmaxwell

    ! Two causal log-log branches with a positive one-sided energy-density
    ! jump at Pt=1 MeV/fm^3. Both branches include Pt in their interpolation
    ! domains; epsilon_high(Pt)-epsilon_low(Pt)=50 MeV/fm^3.
    high_P=[1.0_wp,1.2_wp,1.5_wp,2.0_wp]
    high_eps=[153.0_wp,153.6_wp,154.5_wp,156.0_wp]
    low_P=[0.01_wp,0.1_wp,0.5_wp,1.0_wp]
    low_eps=[100.03_wp,100.3_wp,101.5_wp,103.0_wp]
    call solve_tov_love_maxwell(4,high_P,high_eps,4,low_P,low_eps,1.0_wp,50.0_wp,1.8_wp,rmaxwell)
    if(.not.rmaxwell%converged)then
      write(*,'(2x,a,a)')'FAIL valid Maxwell TOV-Love fixture: ',trim(rmaxwell%status)
      fails=fails+1
    else if(rmaxwell%delta_y>=0._wp.or.rmaxwell%r_transition_km<=0._wp)then
      write(*,'(2x,a)')'FAIL Maxwell TOV-Love transition diagnostics'
      fails=fails+1
    else
      write(*,'(2x,a,es12.4,a,es12.4)')'PASS Maxwell TOV-Love jump: R_t=', &
        rmaxwell%r_transition_km,' km Delta-y=',rmaxwell%delta_y
    end if

    ! Old code accepted these out-of-domain Pt values because it checked the
    ! wrong branch endpoints and branch_interp silently clamped them.
    high_bad=[1.1_wp,1.3_wp,1.6_wp,2.0_wp]
    call solve_tov_love_maxwell(4,high_bad,high_eps,4,low_P,low_eps,1.0_wp,50.0_wp,1.8_wp,rbad)
    if(rbad%converged)then
      write(*,'(2x,a)')'FAIL Maxwell Pt below high branch minimum was not rejected'
      fails=fails+1
    else
      write(*,'(2x,a)')'PASS Maxwell Pt outside high branch rejected'
    end if

    low_bad=[0.01_wp,0.1_wp,0.5_wp,0.8_wp]
    call solve_tov_love_maxwell(4,high_P,high_eps,4,low_bad,low_eps,1.0_wp,50.0_wp,1.8_wp,rbad)
    if(rbad%converged)then
      write(*,'(2x,a)')'FAIL Maxwell Pt above low branch maximum was not rejected'
      fails=fails+1
    else
      write(*,'(2x,a)')'PASS Maxwell Pt outside low branch rejected'
    end if

    call solve_tov_love_maxwell(4,high_P,high_eps,4,low_P,low_eps,low_P(1),50.0_wp,1.8_wp,rbad)
    if(rbad%converged)then
      write(*,'(2x,a)')'FAIL Maxwell Pt at low-branch surface cutoff was not rejected'
      fails=fails+1
    else
      write(*,'(2x,a)')'PASS Maxwell Pt at low-branch surface cutoff rejected'
    end if

    call solve_tov_love_maxwell(4,high_P,high_eps,4,low_P,low_eps,1.0_wp,49.0_wp,1.8_wp,rbad)
    if(rbad%converged)then
      write(*,'(2x,a)')'FAIL inconsistent supplied delta-epsilon was not rejected'
      fails=fails+1
    else
      write(*,'(2x,a)')'PASS inconsistent delta-epsilon rejected'
    end if

    call solve_tov_love_maxwell(4,high_P,high_eps,4,low_P,low_eps,1.0_wp,50.0_wp,2.1_wp,rbad)
    if(rbad%converged)then
      write(*,'(2x,a)')'FAIL Maxwell central pressure above high branch was not rejected'
      fails=fails+1
    else
      write(*,'(2x,a)')'PASS Maxwell central pressure outside high branch rejected'
    end if
  end block

  call check_thermodynamic_jump_consistency(fails)

  if (fails == 0) then
    write(*,'(a)') 'TOV-Love smooth + Maxwell branches: PASS'
  else
    write(*,'(a,i0,a)') 'TOV-Love smooth + Maxwell branches: FAIL (', fails, ' checks)'
    stop 1
  end if
contains

  subroutine check_thermodynamic_jump_consistency(failures)
    integer,intent(inout)::failures
    type(eos_data_t)::eos
    type(verify_report_t)::report
    real(wp)::nvals(11),k,m,p,eps,mu,delta_n,delta_eps
    integer::i
    ! The 0.30 -> 0.33 transition density ratio is below exp(0.25), so
    ! the derivative mask must come from the thermodynamic jump detector,
    ! not the generic maximum-log-step guard.
    nvals=[0.10_wp,0.12_wp,0.15_wp,0.18_wp,0.21_wp,0.25_wp, &
           0.30_wp,0.33_wp,0.38_wp,0.44_wp,0.50_wp]
    call allocate_eos(eos,size(nvals))
    eos%mn=939._wp;eos%mp=939._wp;eos%il=1
    eos%temperature_mev=0._wp;eos%beta_equilibrated=.true.
    eos%composition_available=.false.;eos%composition_file_present=.false.
    do i=1,size(nvals)
      if(i<=7)then
        k=1._wp/(0.30_wp**2)
        m=950._wp-2._wp*k*0.30_wp
      else
        k=1._wp/(0.33_wp**2)
        m=950._wp-2._wp*k*0.33_wp
      endif
      p=k*nvals(i)**2
      eps=m*nvals(i)+p
      mu=m+2._wp*k*nvals(i)
      eos%nb(i)=nvals(i)
      eos%q(1,i)=p/nvals(i)
      eos%q(2,i)=0._wp
      eos%q(3,i)=mu/eos%mn-1._wp
      eos%q(4,i)=0._wp;eos%q(5,i)=0._wp
      eos%q(6,i)=eps/(nvals(i)*eos%mn)-1._wp
      eos%q(7,i)=eos%q(6,i)
    enddo

    call check_consistency(eos,report)
    delta_n=eos%nb(8)-eos%nb(7)
    delta_eps=eos%nb(8)*(eos%q(7,8)+1._wp)*eos%mn- &
              eos%nb(7)*(eos%q(7,7)+1._wp)*eos%mn
    if(.not.report%overall.or.report%thermodynamic_jump_boundaries/=1.or. &
       report%thermodynamic_jump_failures/=0.or.report%derivative_masked<2.or. &
       abs(delta_eps-950._wp*delta_n)>1.e-9_wp*max(abs(delta_eps),1._wp))then
      write(*,'(a,l1,a,i0,a,i0,a,i0)')'FAIL valid thermo-only Maxwell consistency overall=', &
        report%overall,' jump_boundaries=',report%thermodynamic_jump_boundaries, &
        ' jump_failures=',report%thermodynamic_jump_failures,' derivative_masked=',report%derivative_masked
      failures=failures+1
    else
      write(*,'(a)')'PASS thermo-only Maxwell jump detected and derivative stencil masked'
    endif

    ! CompOSE Il=0 stores Q3=mu_Euler-Yq*mu_q. The jump detector must
    ! recover mu_Euler=(epsilon+P)/n_B, not compare the raw Q3 field.
    eos%il=0
    eos%charge_fraction_grid=0.1_wp
    eos%q(4,:)=0.2_wp
    eos%q(3,:)=eos%q(3,:)-eos%charge_fraction_grid*eos%q(4,:)
    call check_consistency(eos,report)
    if(.not.report%overall.or.report%thermodynamic_jump_boundaries/=1.or. &
       report%thermodynamic_jump_failures/=0)then
      write(*,'(a,l1,a,i0,a,i0)')'FAIL Il=0 thermo-only Maxwell jump overall=', &
        report%overall,' jump_boundaries=',report%thermodynamic_jump_boundaries, &
        ' jump_failures=',report%thermodynamic_jump_failures
      failures=failures+1
    else
      write(*,'(a)')'PASS Il=0 charge-shifted Q3 uses Euler mu for jump closure'
    endif

    ! Corrupt only the high-side energy row. The transition is still at common
    ! pressure, but the Euler and latent-energy identities must reject it.
    eos%q(7,8)=eos%q(7,8)+1.e-3_wp
    call check_consistency(eos,report)
    if(report%thermodynamic_jump_failures/=1.or.report%overall)then
      write(*,'(a,i0,a,l1)')'FAIL corrupted Maxwell jump gate failures=', &
        report%thermodynamic_jump_failures,' overall=',report%overall
      failures=failures+1
    else
      write(*,'(a)')'PASS corrupted Maxwell latent-energy jump rejected'
    endif
  end subroutine check_thermodynamic_jump_consistency

end program test_tov_love
