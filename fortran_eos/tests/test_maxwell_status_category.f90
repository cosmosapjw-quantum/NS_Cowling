program test_maxwell_status_category
  !! Native unit test for maxwell_status_category (§3/§4 structured status
  !! taxonomy). Verifies the Fortran-side mapping from a v1 kernel status
  !! string to a §3 terminal-outcome identifier matches the expected category
  !! for a representative sample of every category, and that the mapping is
  !! data-driven from the status string only (never the model name).
  use, intrinsic :: iso_fortran_env, only: real64
  use compose_eos_types, only: wp
  use compose_eos_phase_transition, only: maxwell_status_category, MAXWELL_STATUS_UNKNOWN
  implicit none
  integer :: fails, i
  character(len=160) :: statuses(12)
  character(len=40)  :: expected(12)
  character(len=40)  :: got

  fails = 0
  ! representative v1 kernel status strings, one per §3 category
  statuses(1) = 'phase tables must use the same neutron and proton rest-mass conventions'
  expected(1) = 'MASS_CONVENTION_FAIL'
  statuses(2) = 'Maxwell construction requires beta-equilibrium lepton potential Q5 approximately zero'
  expected(2) = 'BETA_EQUILIBRIUM_FAIL'
  statuses(3) = 'provider Q3 chemical potential disagrees with the Euler value (relative tolerance 5e-4)'
  expected(3) = 'FIRST_LAW_IDENTITY_FAIL'
  statuses(4) = 'P(mu_B) interpolation violates the causal bound c_s^2/c^2<=1'
  expected(4) = 'CAUSALITY_FAIL'
  statuses(5) = 'the two phases have no overlapping baryon-chemical-potential interval'
  expected(5) = 'NO_COMMON_MU_SUPPORT'
  statuses(6) = 'no pressure equality was found in the common mu_B interval'
  expected(6) = 'NO_COEXISTENCE_CROSSING'
  statuses(7) = 'multiple pressure-equality roots found; v1 accepts exactly one isolated transition'
  expected(7) = 'MULTIPLE_COEXISTENCE_CROSSINGS'
  statuses(8) = 'pressure equality does not connect a distinct denser high-density phase'
  expected(8) = 'NONPOSITIVE_DENSITY_JUMP'
  statuses(9) = 'P(mu_B) interpolation violates positive density or mechanical stability'
  expected(9) = 'INTERPOLANT_SHAPE_FAIL'
  statuses(10) = 'v1 Maxwell construction is thermo-only; composition-bearing inputs are rejected, not discarded'
  expected(10) = 'COMPOSITION_UNSUPPORTED'
  statuses(11) = 'some unrecognized status string'
  expected(11) = MAXWELL_STATUS_UNKNOWN
  statuses(12) = 'phase tables must use the same CompOSE lepton-content and charge-fraction conventions'
  expected(12) = 'SLICE_CONVENTION_MISMATCH'

  write(*,'(a)') 'maxwell_status_category native checks:'
  do i = 1, 12
    got = maxwell_status_category(statuses(i))
    if (trim(got) /= trim(expected(i))) then
      write(*,'(2x,a,i0,2a)') 'FAIL ', i, ' expected ', trim(expected(i))
      write(*,'(6x,a)') 'got '//trim(got)//' for: '//trim(statuses(i))
      fails = fails + 1
    else
      write(*,'(2x,a,i0,2a)') 'PASS ', i, ' -> ', trim(got)
    end if
  end do

  ! phase-prefix stripping: 'low-density phase: <interp>' must map to
  ! INTERPOLANT_SHAPE_FAIL, not be lost
  got = maxwell_status_category('low-density phase: P(mu_B) interpolation violates positive density or mechanical stability')
  if (trim(got) /= 'INTERPOLANT_SHAPE_FAIL') then
    write(*,'(2x,2a)') 'FAIL prefix-strip: got ', trim(got)
    fails = fails + 1
  else
    write(*,'(2x,a)') 'PASS phase-prefix stripped'
  end if

  ! model-name independence: embedding a model name does not change category
  got = maxwell_status_category('SPG_M1 P(mu_B) interpolation violates positive density or mechanical stability')
  if (trim(got) /= 'INTERPOLANT_SHAPE_FAIL') then
    write(*,'(2x,2a)') 'FAIL model-name independence: got ', trim(got)
    fails = fails + 1
  else
    write(*,'(2x,a)') 'PASS model-name independence'
  end if

  if (fails == 0) then
    write(*,'(a)') 'ALL maxwell_status_category native checks PASSED'
    stop 0
  else
    write(*,'(a,i0)') 'maxwell_status_category native checks FAILED: ', fails
    stop 1
  end if
end program test_maxwell_status_category
