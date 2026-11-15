program test_ddrmf_homogeneous
  !! Native unit test for compose_eos_rmf (E390-03 / M3 acceptance: CI pass for
  !! the DD-RMF module). Verifies physical invariants of the homogeneous DD-RMF
  !! solver without invoking the full Python frontend:
  !!   1. load_ddrmf_model succeeds for DD2/DDME2/TW and fails for an unknown id.
  !!   2. solve_ddrmf_homogeneous at nuclear saturation (rho_b=rho_sat, yp=0.5,
  !!      symmetric nuclear matter) returns status='ok', finite energy and
  !!      pressure, and rho_n+rho_p == rho_b (number conservation). Saturation
  !!      pressure is near zero by calibration (DD2 has a small residual).
  !!   3. solve_ddrmf_beta_equilibrium at rho_sat returns status='ok', a proton
  !!      fraction in (0, 0.5), and charge balance rho_p == rho_e + rho_mu.
  !!   4. Couplings and meson masses loaded are positive (calibration sanity).
  use, intrinsic :: iso_fortran_env, only: real64
  use compose_eos_types, only: wp
  use compose_eos_rmf, only: ddrmf_model_t, ddrmf_state_t, ddrmf_beta_result_t, &
                             load_ddrmf_model, solve_ddrmf_homogeneous, &
                             solve_ddrmf_beta_equilibrium
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: fails
  type(ddrmf_model_t) :: model
  type(ddrmf_state_t) :: st
  type(ddrmf_beta_result_t) :: beta
  character(len=120) :: status
  real(wp) :: rho_b, charge_resid
  integer :: i

  fails = 0
  write(*,'(a)') 'DD-RMF homogeneous + beta-equilibrium native checks:'

  ! --- 1. load known models + 4. coupling/meson sanity ---
  do i = 1, 3
    select case (i)
    case (1); call load_ddrmf_model('GPPVA_DD2', model, status)
    case (2); call load_ddrmf_model('GPPVA_DDME2', model, status)
    case (3); call load_ddrmf_model('GPPVA_TW', model, status)
    end select
    if (trim(status) /= 'ok') then
      write(*,'(2x,a,i0,a)') 'FAIL load model ', i, ' -> '//trim(status)
      fails = fails + 1
    else if (min(model%ms, model%mw, model%mr, model%gs0, model%gw0, model%gr0, &
                 model%rho_sat) <= 0._wp) then
      write(*,'(2x,a,i0)') 'FAIL model ', i, ' has nonpositive meson/coupling/sat'
      fails = fails + 1
    else
      write(*,'(2x,a,a)') 'PASS load ', trim(model%id)
    end if
  end do

  ! unknown id must fail (status starts with 'unknown')
  call load_ddrmf_model('NOT_A_MODEL', model, status)
  if (index(status, 'unknown') == 1) then
    write(*,'(2x,a)') 'PASS load rejects unknown model id'
  else
    write(*,'(2x,a)') 'FAIL load did not reject unknown model id'
    fails = fails + 1
  end if

  ! --- 2. homogeneous SNM at saturation ---
  call load_ddrmf_model('GPPVA_DD2', model, status)
  rho_b = model%rho_sat
  call solve_ddrmf_homogeneous(model, rho_b, 0.5_wp, st, status)
  if (trim(status) /= 'ok') then
    write(*,'(2x,a,a)') 'FAIL solve_ddrmf_homogeneous SNM -> ', trim(status)
    fails = fails + 1
  else if (.not. ieee_is_finite(st%eps_b) .or. .not. ieee_is_finite(st%p_b)) then
    write(*,'(2x,a)') 'FAIL SNM energy/pressure not finite'
    fails = fails + 1
  else if (abs(st%rho_n + st%rho_p - rho_b) > 1.e-9_wp) then
    write(*,'(2x,a,es12.4)') 'FAIL SNM number conservation residual: ', st%rho_n+st%rho_p-rho_b
    fails = fails + 1
  else
    write(*,'(2x,a,es12.4,a,es12.4)') 'PASS SNM at sat: |P|=', abs(st%p_b), ' rho_n+rho_p=', st%rho_n+st%rho_p
  end if

  ! --- 3. beta-equilibrium at saturation ---
  call solve_ddrmf_beta_equilibrium(model, rho_b, beta, status)
  if (trim(status) /= 'ok') then
    write(*,'(2x,a,a)') 'FAIL solve_ddrmf_beta_equilibrium -> ', trim(status)
    fails = fails + 1
  else if (beta%yp <= 0._wp .or. beta%yp > 0.5_wp) then
    write(*,'(2x,a,es12.4)') 'FAIL beta yp out of (0,0.5]: ', beta%yp
    fails = fails + 1
  else
    charge_resid = beta%rho_e + beta%rho_mu - beta%baryons%rho_p
    if (abs(charge_resid) > 1.e-6_wp * max(1._wp, beta%baryons%rho_p)) then
      write(*,'(2x,a,es12.4)') 'FAIL charge balance residual: ', charge_resid
      fails = fails + 1
    else
      write(*,'(2x,a,es12.4,a,es12.4)') 'PASS beta: yp=', beta%yp, ' charge_resid=', charge_resid
    end if
  end if

  if (fails == 0) then
    write(*,'(a)') 'DD-RMF homogeneous + beta-equilibrium: PASS'
  else
    write(*,'(a,i0,a)') 'DD-RMF homogeneous + beta-equilibrium: FAIL (', fails, ' checks)'
    stop 1
  end if
end program test_ddrmf_homogeneous
