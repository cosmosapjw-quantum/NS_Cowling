program test_cld_cell
  !! Native unit test for compose_eos_cld (E390-03 / M3 acceptance: CI pass for
  !! the CLD module). Verifies physical invariants of the CLD cell
  !! reconstruction at representative inner-crust densities without invoking the
  !! full Python frontend:
  !!   1. cld_reconstruct_point returns a valid result (f, r_d, rho_i all in
  !!      variational bounds) for DD2 at three inner-crust densities.
  !!   2. The cluster volume fraction f is in (0, 1).
  !!   3. Cluster proton density and gas proton density are non-negative.
  !!   4. Cell-averaged proton fraction is in (0, 0.5] (beta-equilibrated crust).
  !!   5. Independent full-cell proton derivative vanishes, and reported
  !!      energy excludes constituent baryon rest energies.
  !! The lattice/surface terms are explicitly an exploratory proxy, not RPW
  !! or the source-calibrated Scurto2024 functions.
  use, intrinsic :: iso_fortran_env, only: real64
  use compose_eos_types, only: wp
  use compose_eos_cld, only: cld_cell_result_t, cld_reconstruct_point
  use compose_eos_rmf, only: ddrmf_model_t, ddrmf_state_t, load_ddrmf_model, solve_ddrmf_homogeneous
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: fails, i
  type(cld_cell_result_t) :: res, plus_cell, minus_cell
  real(wp) :: grad(3), delta, gas_diff
  type(ddrmf_state_t) :: gas
  real(wp), parameter :: sigma_0 = 1.0_wp   ! MeV fm^-2 (proxy parameter)
  real(wp), parameter :: b = 0.59_wp       ! surface tension stiffness
  real(wp) :: rho_b(3), d_eps, step, total, rest
  type(ddrmf_model_t) :: model
  character(len=80) :: model_status
  character(len=120) :: note

  fails = 0
  write(*,'(a)') 'CLD cell reconstruction native checks:'

  ! three representative inner-crust baryon densities (fm^-3)
  rho_b = [0.01_wp, 0.03_wp, 0.055_wp]
  call load_ddrmf_model('GPPVA_DD2', model, model_status)

  do i = 1, 3
    call cld_reconstruct_point('GPPVA_DD2', rho_b(i), sigma_0, b, res)
    note = ' '
    if(.not.res%converged) then
      write(*,'(2x,a,i0)') 'FAIL nominal CLD optimizer not converged at point ',i
      fails=fails+1
    endif
    ! Derivative of the complete declared cell energy, independently varying
    ! cluster protons at fixed f, R, rho_I, rho_II (charge neutrality fixes n_e).
    step = min(1.e-7_wp, res%rho_ip*1.e-3_wp)
    d_eps = (cell_total(res, res%rho_ip+step)-cell_total(res, res%rho_ip-step))/(2._wp*step*res%f)
    if (abs(d_eps) > 2.e-4_wp) then
      write(*,'(2x,a,i0,a,es14.6)') 'FAIL full proton stationarity point ', i, ': MeV=', d_eps
      fails = fails + 1
    endif
    if (abs(abs(d_eps)-res%beta_resid) > 2.e-4_wp) then
      write(*,'(2x,a,i0)') 'FAIL reported beta residual is not full cell derivative, point ', i
      fails = fails + 1
    endif
    write(*,'(2x,a,es14.6,a,es14.6)') 'Full cell proton derivative MeV=',d_eps,' binding E/A=',res%e_binding_per_baryon
    plus_cell=res; minus_cell=res; delta=1.e-5_wp*res%f
    plus_cell%f=res%f+delta; minus_cell%f=res%f-delta
    plus_cell%rho_ii=(rho_b(i)-plus_cell%f*res%rho_i)/(1._wp-plus_cell%f)
    minus_cell%rho_ii=(rho_b(i)-minus_cell%f*res%rho_i)/(1._wp-minus_cell%f)
    grad(1)=res%f*(cell_total(plus_cell,res%rho_ip)-cell_total(minus_cell,res%rho_ip))/(2._wp*delta*rho_b(i))
    plus_cell=res; minus_cell=res; delta=1.e-5_wp
    plus_cell%r_d=res%r_d*exp(delta); minus_cell%r_d=res%r_d*exp(-delta)
    grad(2)=(cell_total(plus_cell,res%rho_ip)-cell_total(minus_cell,res%rho_ip))/(2._wp*delta*rho_b(i))
    plus_cell=res; minus_cell=res; delta=1.e-5_wp*res%rho_i
    plus_cell%rho_i=res%rho_i+delta; minus_cell%rho_i=res%rho_i-delta
    plus_cell%rho_ii=(rho_b(i)-res%f*plus_cell%rho_i)/(1._wp-res%f)
    minus_cell%rho_ii=(rho_b(i)-res%f*minus_cell%rho_i)/(1._wp-res%f)
    grad(3)=res%rho_i*(cell_total(plus_cell,res%rho_ip)-cell_total(minus_cell,res%rho_ip))/(2._wp*delta*rho_b(i))
    write(*,'(2x,a,i0,a,3es18.9,a,i0)') 'GRAD ',i,' f/logR/rhoI scaled MeV/baryon ',grad,' iterations ',res%n_iter
    if(maxval(abs(grad))>1.e-4_wp) then
      write(*,'(2x,a,i0)') 'FAIL full variational stationarity at point ',i
      fails=fails+1
    endif
    total = cell_total(res, res%rho_ip)
    rest = model%mn*(rho_b(i)-res%rho_e)+model%mp*res%rho_e
    if (abs(res%e_cell+rest-total) > 1.e-8_wp) then
      write(*,'(2x,a,i0,a,es14.6)') 'FAIL binding energy includes baryon rest energy, point ', i, ': E/A=', res%e_binding_per_baryon
      fails = fails + 1
    endif
    if (.not. ieee_is_finite(res%f) .or. .not. ieee_is_finite(res%r_d) .or. &
        .not. ieee_is_finite(res%rho_i) .or. .not. ieee_is_finite(res%e_cell)) then
      write(*,'(2x,a,i0,a,es9.2)') 'FAIL CLD point ', i, ' not finite at rho_b=', rho_b(i)
      fails = fails + 1
      cycle
    end if
    ! bounds: cluster volume fraction in (0,1), radius in [2,40] fm,
    ! cluster baryon density positive and below 0.5 fm^-3.
    if (res%f <= 0._wp .or. res%f >= 1._wp) then
      write(*,'(2x,a,i0,a,es9.2)') 'FAIL CLD f out of (0,1) at point ', i, ': f=', res%f
      fails = fails + 1
    else if (res%r_d < 2._wp .or. res%r_d > 40._wp) then
      write(*,'(2x,a,i0,a,es9.2)') 'FAIL CLD r_d out of [2,40] at point ', i, ': r_d=', res%r_d
      fails = fails + 1
    else if (res%rho_i <= 0._wp .or. res%rho_i > 0.5_wp) then
      write(*,'(2x,a,i0,a,es9.2)') 'FAIL CLD rho_i out of (0,0.5] at point ', i, ': rho_i=', res%rho_i
      fails = fails + 1
    else if (res%rho_ip < 0._wp .or. res%rho_iip < 0._wp) then
      write(*,'(2x,a,i0)') 'FAIL CLD negative cluster/gas proton density at point ', i
      fails = fails + 1
    else if (res%yp_cell <= 0._wp .or. res%yp_cell > 0.5_wp) then
      write(*,'(2x,a,i0,a,es9.2)') 'FAIL CLD yp_cell out of (0,0.5] at point ', i, ': yp=', res%yp_cell
      fails = fails + 1
    else
      write(*,'(2x,a,i0,a,es9.2,a,es9.2,a,es9.2,a,es9.2)') 'PASS CLD point ', i, &
        ' rho_b=', rho_b(i), ' f=', res%f, ' r_d=', res%r_d, ' yp=', res%yp_cell
    end if
  end do

  ! unknown model id must return a zeroed/invalid result (f stays at default 0
  ! or finite but the reconstruction should not crash). We only require it not
  ! to crash and to leave f in [0,1].
  call cld_reconstruct_point('NOT_A_MODEL', 0.03_wp, sigma_0, b, res)
  if (ieee_is_finite(res%f) .and. res%f >= 0._wp .and. res%f <= 1._wp) then
    write(*,'(2x,a)') 'PASS CLD unknown-model does not crash'
  else
    write(*,'(2x,a)') 'FAIL CLD unknown-model crashed or out of bounds'
    fails = fails + 1
  end if

  call cld_reconstruct_point('GPPVA_DD2', -0.01_wp, sigma_0, b, res)
  if(res%converged.or.index(res%status,'invalid')/=1) fails=fails+1
  call cld_reconstruct_point('GPPVA_DD2', 0.01_wp, -1._wp, b, res)
  if(res%converged.or.index(res%status,'invalid')/=1) fails=fails+1

  if (fails == 0) then
    write(*,'(a)') 'CLD cell reconstruction: PASS'
  else
    write(*,'(a,i0,a)') 'CLD cell reconstruction: FAIL (', fails, ' checks)'
    stop 1
  end if
contains
  function cell_total(r, nip) result(e)
    type(cld_cell_result_t), intent(in) :: r
    real(wp), intent(in) :: nip
    real(wp) :: e, y, ne, k, ef, el, phi, surf
    real(wp), parameter :: hc=197.3269804_wp, pi=acos(-1._wp), me=0.51099895_wp
    type(ddrmf_state_t) :: si, sg
    character(len=80) :: status
    call solve_ddrmf_homogeneous(model, r%rho_i, nip/r%rho_i, si, status)
    if(status/='ok') stop 2
    call solve_ddrmf_homogeneous(model, r%rho_ii, 0._wp, sg, status)
    if(status/='ok') stop 2
    y=nip/r%rho_i; ne=r%f*nip
    k=hc*(3._wp*pi*pi*ne)**(1._wp/3._wp); ef=sqrt(k*k+me*me)
    el=(k*(2._wp*k*k+me*me)*ef-me**4*log((k+ef)/me))/(8._wp*pi*pi*hc**3)
    phi=0.5_wp*sqrt(max(0._wp,2._wp-3._wp*r%f))
    surf=3._wp*r%f*sigma_0*(y*(1._wp-y))**b/r%r_d
    e=r%f*si%eps_b+(1._wp-r%f)*sg%eps_b+el+surf+ &
       2._wp*pi*(hc/137.036_wp)*r%f*phi*r%r_d**2*nip**2
  end function cell_total
end program test_cld_cell
