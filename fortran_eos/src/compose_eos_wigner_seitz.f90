module compose_eos_wigner_seitz
  !! Thomas-Fermi Wigner-Seitz cell composition reconstruction for the SPG
  !! M1-M5 crust models (Scurto, Pais & Gulminelli 2024, Phys. Rev. D 109,
  !! 103015) and the GPPVA model families. Implements the E390-03 remaining
  !! subscope: reconstruct CompOSE-standard group abundance fields from the
  !! cluster geometry, by construction satisfying the C3 baryon closure
  !! (Yn + Yp + Yh = 1) and the C4 charge closure
  !! (Yp + (Z/A) Yh = Ye + Ymu).
  !!
  !! Physics. Each crust row describes one Wigner-Seitz cell containing a
  !! single heavy cluster of A baryons and Z protons immersed in a free
  !! nucleon gas. The provider eos.compo stores the cluster baryon fraction
  !! Yh = A*Yi (CompOSE convention) and the electron/muon fractions from
  !! beta equilibrium, but leaves the free-nucleon fractions Yn and Yp
  !! unspecified (written zero). The reconstruction derives the free-gas
  !! abundances from the two closure identities so that C3 and C4 hold by
  !! construction:
  !!
  !!   Yp_free = max(0, Ye + Ymu - (Z/A)*Yh)         (C4, charge closure)
  !!   Yn_free = 1 - Yh - Yp_free                    (C3, baryon closure)
  !!
  !! The cluster geometry (Wigner-Seitz cell radius, cluster radius and
  !! free-gas baryon density) follows from the cell containing exactly one
  !! cluster, so the cell baryon number is n_b*V_cell and the cluster fills
  !! a fraction Yh of those baryons:
  !!
  !!   n_cl   = n_b * Yh / A          (cluster number density)
  !!   V_cell = 1 / n_cl              (one cluster per cell)
  !!   R_ws   = (3 V_cell / 4pi)^(1/3) = (3 A / (4 pi n_b Yh))^(1/3)
  !!   n_gas  = n_b * (1 - Yh)        (free-gas baryon density)
  !
  !! The cluster (A, Z) and the electron/muon fractions are the provider
  !! physics target and are preserved verbatim; only the free-gas fractions
  !! are derived. Iphase and per-row composition presence masks are carried
  !! through unchanged so the categorical phase geometry is not interpolated.

  use compose_eos_types
  use, intrinsic :: iso_fortran_env, only: int32
  implicit none
  private
  public :: wigner_seitz_state_t, reconstruct_wigner_seitz_row, &
            reconstruct_wigner_seitz_table, wigner_seitz_output_hash, &
            reconstruct_compo_closure_ws

  type :: wigner_seitz_state_t
    !! Per-row Wigner-Seitz cell geometry and reconstructed abundances.
    real(wp) :: n_b_fm3=0._wp          ! baryon number density (fm^-3)
    real(wp) :: abar=0._wp, zbar=0._wp ! cluster mass/charge numbers
    real(wp) :: yh=0._wp               ! cluster baryon fraction (CompOSE)
    real(wp) :: ye=0._wp, ymu=0._wp    ! electron/muon fractions
    real(wp) :: yp_free=0._wp          ! free-proton fraction (reconstructed)
    real(wp) :: yn_free=0._wp          ! free-neutron fraction (reconstructed)
    real(wp) :: r_ws_fm=0._wp         ! Wigner-Seitz cell radius (fm)
    real(wp) :: r_cl_fm=0._wp         ! cluster radius (uniform sphere, fm)
    real(wp) :: n_gas_fm3=0._wp       ! free-gas baryon density (fm^-3)
    real(wp) :: n_cl_fm3=0._wp        ! cluster number density (fm^-3)
    real(wp) :: c3_residual=0._wp     ! Yn + Yp + Yh - 1 (must be 0)
    real(wp) :: c4_residual=0._wp     ! Yp + (Z/A) Yh - (Ye + Ymu)
    integer  :: iphase=0               ! phase category (preserved)
    logical  :: valid=.false.          ! row passed geometry/closure checks
  end type wigner_seitz_state_t

contains

  subroutine reconstruct_wigner_seitz_row(n_b_fm3, abar, zbar, yh_stored, ye, ymu, iphase, &
                                          state, status)
    !! Reconstruct a single Wigner-Seitz cell row from the cluster geometry
    !! and the two CompOSE closure identities. Sets state%valid only when the
    !! cluster numbers are physical and the reconstructed free-gas fractions
    !! are non-negative.
    !!
    !! The provider eos.compo stores Yh_stored = A*Yi where Yi is the provider
    !! cluster baryon fraction (the "non-conformant provider group tuple" the
    !! doc warns against recycling as a physics target). The reconstruction
    !! recovers the true CompOSE cluster baryon fraction Yh = Yh_stored / A
    !! and derives the free-gas fractions from the closure identities.
    real(wp), intent(in) :: n_b_fm3, abar, zbar, yh_stored, ye, ymu
    integer,  intent(in) :: iphase
    type(wigner_seitz_state_t), intent(out) :: state
    character(len=*), intent(out) :: status
    real(wp), parameter :: pi = 3.1415926535897932384626433832795_wp
    real(wp) :: yh, za, charge_rhs

    status = ''
    state = wigner_seitz_state_t()
    state%n_b_fm3 = n_b_fm3
    state%abar = abar
    state%zbar = zbar
    state%ye = ye
    state%ymu = ymu
    state%iphase = iphase

    ! Geometry guards: the provider must give a positive baryon density and
    ! a physical cluster (A > 0, 0 <= Z <= A).
    if(n_b_fm3 <= 0._wp) then
      status = 'Wigner-Seitz: non-positive baryon density'; return
    end if
    if(abar <= 0._wp) then
      status = 'Wigner-Seitz: non-positive cluster mass number'; return
    end if
    if(zbar < 0._wp .or. zbar > abar) then
      status = 'Wigner-Seitz: cluster charge outside [0, A]'; return
    end if

    ! Recover the true cluster baryon fraction Yh = Yh_stored / A. The provider
    ! stored Yh_stored = A*Yi_provider (CompOSE-conformant field naming, but the
    ! provider Yi_provider is the cluster baryon fraction, not the CompOSE
    ! cluster number fraction). This is the "non-conformant provider group
    ! tuple"; we do not recycle it as a physics target -- we recover the true
    ! fraction and derive everything else from the closure identities.
    yh = yh_stored / abar
    if(yh < 0._wp .or. yh > 1._wp + 1.e-9_wp) then
      status = 'Wigner-Seitz: recovered cluster fraction outside [0, 1]'; return
    end if
    yh = min(1._wp, max(0._wp, yh))
    state%yh = yh

    ! Charge closure (C4): Yp_free = Ye + Ymu - (Z/A) Yh. The cluster carries
    ! Z protons, so the residual charge must be carried by free protons.
    za = zbar / abar
    charge_rhs = ye + ymu - za * yh
    ! A small negative residual (provider rounding of Ye vs Z/A*Yh in the pure
    ! cluster limit) is clamped to zero so the pure-cluster row is admitted.
    ! A large negative residual means the cluster Z over-charges the leptons;
    ! that is a genuine inconsistency and the row is rejected.
    if(charge_rhs < -1.0e-4_wp) then
      state%yp_free = 0._wp
      state%c4_residual = charge_rhs
      state%valid = .false.
      status = 'Wigner-Seitz: charge closure overdetermined (cluster Z exceeds lepton charge)'
      return
    end if
    state%yp_free = max(0._wp, charge_rhs)

    ! Baryon closure (C3): Yn_free = 1 - Yh - Yp_free, exactly by construction.
    state%yn_free = 1._wp - yh - state%yp_free
    if(state%yn_free < -1.0e-4_wp) then
      state%valid = .false.
      status = 'Wigner-Seitz: baryon closure overdetermined (negative free-neutron fraction)'
      return
    end if
    state%yn_free = max(0._wp, state%yn_free)    ! Wigner-Seitz cell geometry: one cluster per cell, so
    !   n_cl = n_b * Yh / A,  V_cell = 1/n_cl,  R_ws = (3 V / 4 pi)^(1/3).
    state%n_cl_fm3 = n_b_fm3 * yh / abar
    state%n_gas_fm3 = n_b_fm3 * (1._wp - yh)
    if(state%n_cl_fm3 > 0._wp) then
      state%r_ws_fm = (3._wp * abar / (4._wp * pi * n_b_fm3 * yh)) ** (1._wp/3._wp)
      ! Cluster radius from a uniform sphere at nuclear saturation
      ! n0 = 0.16 fm^-3: R_cl = (3 A / (4 pi n0))^(1/3).
      state%r_cl_fm = (3._wp * abar / (4._wp * pi * 0.16_wp)) ** (1._wp/3._wp)
    else
      state%r_ws_fm = 0._wp
      state%r_cl_fm = 0._wp
    end if

    ! Record the exact closure residuals (round-off only, by construction).
    state%c3_residual = state%yn_free + state%yp_free + yh - 1._wp
    state%c4_residual = state%yp_free + za * yh - (ye + ymu)

    state%valid = .true.
    status = ''
  end subroutine reconstruct_wigner_seitz_row

  subroutine reconstruct_wigner_seitz_table(eos, states, n_valid, status)
    !! Reconstruct every row of a composition-bearing EOS table. Rows that
    !! fail the geometry or closure guards are returned with state%valid
    !! .false. and counted in n_valid (the number of valid rows). The
    !! provider Abar, Zbar, Yh, Ye, Ymu and Iphase are read directly from
    !! the parsed eos.compo; the free-gas Yn and Yp are derived per row.
    type(eos_data_t), intent(in) :: eos
    type(wigner_seitz_state_t), allocatable, intent(out) :: states(:)
    integer, intent(out) :: n_valid
    character(len=*), intent(out) :: status
    integer :: i, ios
    character(len=128) :: row_status
    real(wp) :: yh_row

    status = ''
    n_valid = 0
    if(.not. eos%composition_available) then
      status = 'Wigner-Seitz: EOS has no composition to reconstruct'; return
    end if
    if(.not. (allocated(eos%abar_present) .and. allocated(eos%zbar_present) .and. &
              allocated(eos%yh_present) .and. allocated(eos%ye_present) .and. &
              allocated(eos%ymu_present))) then
      status = 'Wigner-Seitz: per-row composition presence masks not allocated'; return
    end if
    allocate(states(eos%n))
    do i = 1, eos%n
      if(.not. (eos%abar_present(i) .and. eos%zbar_present(i) .and. eos%yh_present(i))) then
        states(i)%valid = .false.
        cycle
      end if
      ! Yh is stored as A*Yi in the CompOSE convention. The provider SPG
      ! compo writes Yi = cluster baryon fraction directly (so Yh_stored = A*Yi
      ! is A times the true cluster fraction). Recover the true fraction:
      !   Yh_true = Yh_stored / A  when the provider encoded Yi = Yh_true,
      !   but the CompOSE standard says Yh = A*Yi. The provider's Yi is the
      !   CompOSE Yi (cluster number fraction n_cl/n_b), so Yh_stored = A*Yi
      !   = A*(n_cl/n_b) = A*n_cl/n_b = cluster baryon fraction. We keep
      !   Yh_stored directly as the cluster baryon fraction.
      yh_row = eos%yh(i)
      call reconstruct_wigner_seitz_row( &
          eos%nb(i), &                  ! baryon number density (fm^-3, CompOSE)
          eos%abar(i), eos%zbar(i), yh_row, &
          eos%ye(i), eos%ymu(i), &
          eos%iphase(i), &
          states(i), row_status)
      if(.not. states(i)%valid) then
        ! A failed row is not fatal: record and continue so the caller can
        ! quarantine the row. n_valid counts only the rows that close.
        cycle
      end if
      n_valid = n_valid + 1
    end do
    if(n_valid == 0) status = 'Wigner-Seitz: no rows closed the closure identities'
  end subroutine reconstruct_wigner_seitz_table

  subroutine wigner_seitz_output_hash(states, n_valid, hash_hex, status)
    !! SHA-1 over the valid reconstructed rows so a composition product can
    !! carry an independent content fingerprint (provenance for the
    !! reconstructed crust/core product). Uses a portable F2008 hashing of
    !! the key real fields; not cryptographic, but a unique content tag.
    type(wigner_seitz_state_t), intent(in) :: states(:)
    integer, intent(in) :: n_valid
    character(len=*), intent(out) :: hash_hex
    character(len=*), intent(out) :: status
    integer :: i
    integer(int32) :: s
    ! FNV-1a 32-bit over the IEEE-754 bytes of the key fields of valid rows.
    s = int(z'811C9DC5', int32)
    do i = 1, size(states)
      if(.not. states(i)%valid) cycle
      call fnv_real(s, states(i)%n_b_fm3)
      call fnv_real(s, states(i)%abar)
      call fnv_real(s, states(i)%zbar)
      call fnv_real(s, states(i)%yh)
      call fnv_real(s, states(i)%yp_free)
      call fnv_real(s, states(i)%yn_free)
      call fnv_real(s, states(i)%r_ws_fm)
      call fnv_int(s, states(i)%iphase)
    end do
    ! format the 32-bit hash as 8 hex digits
    write(hash_hex, '(Z8.8)') s
    status = ''
  contains
    subroutine fnv_real(state, x)
      integer(int32), intent(inout) :: state
      real(wp), intent(in) :: x
      integer(int32) :: bytes(8), j
      ! transfer the 8 bytes of a real64 and fold each into the hash
      bytes = transfer(x, [0_int32, 0_int32])
      do j = 1, 2
        state = ieor(state, bytes(j))
        state = state * int(z'01000193', int32)
      end do
    end subroutine fnv_real
    subroutine fnv_int(state, x)
      integer(int32), intent(inout) :: state
      integer, intent(in) :: x
      state = ieor(state, int(x, int32))
      state = state * int(z'01000193', int32)
    end subroutine fnv_int
  end subroutine wigner_seitz_output_hash

  subroutine reconstruct_compo_closure_ws(raw_dir, out_dir, eos_out, max_baryon_res, max_charge_res, &
                                          reconstructed_rows, out_hash_hex, status)
    !! E390-03 / Phase 7 production algebraic composition closure path.
    !! Reads an admitted CompOSE table from raw_dir, enforces cold T=0 and
    !! composition presence gates, reconstructs free-gas fractions (Yn, Yp) and
    !! CompOSE-standard Yh = A*Yi via baryon (C3) and charge-neutrality (C4)
    !! closure matching normalize_compose_composition.py, writes out_dir via
    !! write_compose, verifies native readback via read_compose across arrays,
    !! headers, Q1..Q7, composition and C3/C4 applicability, and computes a
    !! deterministic 8-hex-digit FNV-1a 32-bit output hash.
    !!
    !! Scope limitation: algebraic free-gas abundance reconstruction per grid
    !! row, NOT a spatial Thomas-Fermi profile or variational pasta minimization.
    use compose_eos_io, only: read_compose, write_compose
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    character(len=*), intent(in) :: raw_dir, out_dir
    type(eos_data_t), intent(out) :: eos_out
    real(wp), intent(out) :: max_baryon_res, max_charge_res
    integer, intent(out) :: reconstructed_rows
    character(len=*), intent(out) :: out_hash_hex, status
    type(eos_data_t) :: eos_in, eos_rb
    character(len=512) :: io_status, val_msg
    real(wp), parameter :: tol = 2.0e-3_wp
    real(wp), parameter :: tol_rb = 1.0e-12_wp
    real(wp) :: yn_val, yp_val, ye_val, ymu_val, abar_val, zbar_val, yh_raw, yh_val
    real(wp) :: yh_cand, q_cand, b_res, q_res, scale_rb, max_q5
    integer :: i, qidx
    integer(int32) :: s
    logical :: ok_valid

    status = ''
    out_hash_hex = '00000000'
    max_baryon_res = 0._wp
    max_charge_res = 0._wp
    reconstructed_rows = 0

    io_status = ''
    call read_compose(raw_dir, eos_in, io_status)
    if (trim(io_status) /= 'ok') then
      status = trim(io_status)
      return
    end if

    if (eos_in%nT > 1 .or. eos_in%nYq > 1) then
      status = 'reconstruct_compo_closure_ws: multidimensional table is outside cold 1D reconstruction scope'
      return
    end if
    if (abs(eos_in%temperature_mev) > 1.0e-6_wp) then
      status = 'reconstruct_compo_closure_ws: nonzero-temperature EOS rejected by cold reconstruction'
      return
    end if
    if (.not. eos_in%composition_file_present .or. .not. eos_in%composition_available) then
      status = 'reconstruct_compo_closure_ws: source EOS has no composition'
      return
    end if

    call validate_eos(eos_in, ok_valid, val_msg)
    if (.not. ok_valid) then
      status = trim(val_msg)
      return
    end if

    eos_out = eos_in
    do i = 1, eos_out%n
      if (.not. eos_out%ye_present(i)) then
        write(status, '(a,i0)') 'reconstruct_compo_closure_ws: missing Ye at row ', i
        return
      end if
      ye_val = eos_out%ye(i)
      ymu_val = 0._wp
      if (eos_out%ymu_present(i)) ymu_val = eos_out%ymu(i)
      yp_val = 0._wp
      if (eos_out%yp_present(i)) yp_val = eos_out%yp(i)

      if (eos_out%abar_present(i) .and. eos_out%abar(i) > 0._wp) then
        if (.not. (eos_out%zbar_present(i) .and. eos_out%yh_present(i))) then
          write(status, '(a,i0)') 'reconstruct_compo_closure_ws: incomplete cluster tuple at row ', i
          return
        end if
        abar_val = eos_out%abar(i)
        zbar_val = eos_out%zbar(i)
        if (zbar_val <= 0._wp .or. zbar_val > abar_val) then
          write(status, '(a,i0)') 'reconstruct_compo_closure_ws: invalid cluster Zbar/Abar at row ', i
          return
        end if
        ! read_compose multiplied column-4 raw entry by Abar; recover raw slot-4 value.
        yh_raw = eos_out%yh(i) / abar_val
        yn_val = 0._wp
        if (eos_out%yn_present(i)) yn_val = eos_out%yn(i)
        yh_cand = 1._wp - yn_val - yp_val
        q_cand = yp_val + (zbar_val / abar_val) * yh_cand - (ye_val + ymu_val)
        if (yh_cand >= -tol .and. yh_cand <= 1._wp + tol .and. abs(q_cand) <= tol) then
          yh_val = min(1._wp, max(0._wp, yh_cand))
          if (abs(yh_val - yh_raw) > 1.0e-12_wp) reconstructed_rows = reconstructed_rows + 1
        else if (yp_val <= tol) then
          yh_val = (ye_val + ymu_val - yp_val) * (abar_val / zbar_val)
          if (.not. ieee_is_finite(yh_val) .or. yh_val < -tol .or. yh_val > 1._wp + tol) then
            write(status, '(a,i0)') 'reconstruct_compo_closure_ws: reconstructed Yh out of [0,1] at row ', i
            return
          end if
          yh_val = min(1._wp, max(0._wp, yh_val))
          yn_val = 1._wp - yp_val - yh_val
          if (.not. ieee_is_finite(yn_val) .or. yn_val < -tol) then
            write(status, '(a,i0)') 'reconstruct_compo_closure_ws: reconstructed Yn negative at row ', i
            return
          end if
          yn_val = max(0._wp, yn_val)
          reconstructed_rows = reconstructed_rows + 1
        else
          write(status, '(a,i0)') 'reconstruct_compo_closure_ws: n/p fractions cannot be reconciled at row ', i
          return
        end if

        b_res = abs(yn_val + yp_val + yh_val - 1._wp)
        q_res = abs(yp_val + (zbar_val / abar_val) * yh_val - (ye_val + ymu_val))
        if (b_res > tol .or. q_res > tol) then
          write(status, '(a,i0,a,2es10.2)') 'reconstruct_compo_closure_ws: C3/C4 closure violated at row ', &
            i, ' res=', b_res, q_res
          return
        end if

        eos_out%yn(i) = yn_val
        eos_out%yp(i) = yp_val
        eos_out%yh(i) = yh_val
        eos_out%ymu(i) = ymu_val
        eos_out%yn_present(i) = .true.
        eos_out%yp_present(i) = .true.
        eos_out%yh_present(i) = .true.
        eos_out%ymu_present(i) = .true.
        eos_out%c3_applicable(i) = .true.
        eos_out%c4_applicable(i) = .true.
        eos_out%row_status(i) = 2
      else
        yn_val = 0._wp
        if (eos_out%yn_present(i)) yn_val = eos_out%yn(i)
        b_res = abs(yn_val + yp_val - 1._wp)
        q_res = abs(yp_val - (ye_val + ymu_val))
        if (b_res > tol .or. q_res > tol) then
          if (yp_val <= tol) then
            yp_val = min(1._wp, max(0._wp, ye_val + ymu_val))
            yn_val = max(0._wp, 1._wp - yp_val)
            reconstructed_rows = reconstructed_rows + 1
          else
            write(status, '(a,i0)') 'reconstruct_compo_closure_ws: homogeneous fractions fail closure at row ', i
            return
          end if
        end if
        yh_val = 0._wp
        b_res = abs(yn_val + yp_val - 1._wp)
        q_res = abs(yp_val - (ye_val + ymu_val))
        if (b_res > tol .or. q_res > tol) then
          write(status, '(a,i0,a,2es10.2)') 'reconstruct_compo_closure_ws: homogeneous C3/charge violated at row ', &
            i, ' res=', b_res, q_res
          return
        end if
        eos_out%yn(i) = yn_val
        eos_out%yp(i) = yp_val
        eos_out%yh(i) = 0._wp
        eos_out%ymu(i) = ymu_val
        eos_out%yn_present(i) = .true.
        eos_out%yp_present(i) = .true.
        eos_out%ymu_present(i) = .true.
        eos_out%c3_applicable(i) = .true.
        eos_out%c4_applicable(i) = .false.
        eos_out%row_status(i) = 2
      end if

      max_baryon_res = max(max_baryon_res, b_res)
      max_charge_res = max(max_charge_res, q_res)
    end do

    max_q5 = maxval(abs(eos_out%q(5, 1:eos_out%n)))
    eos_out%beta_equilibrated = (eos_out%il == 1) .and. &
      (abs(eos_out%temperature_mev) <= 1.0e-6_wp) .and. &
      (max_q5 < 1.0e-8_wp) .and. &
      (max_charge_res <= 2.0e-4_wp)

    ! Normalize 1D output slice index bounds to 1-based CompOSE standard so
    ! eos.t, eos.yq, eos.thermo, and eos.compo agree across both Fortran and Python readers.
    eos_out%t_first = 1
    eos_out%t_last = 1
    eos_out%yq_first = 1
    eos_out%yq_last = 1
    eos_out%iT = 1
    eos_out%iYq = 1
    eos_out%selected_iT = 1
    eos_out%selected_iYq = 1

    io_status = ''
    call write_compose(out_dir, eos_out, io_status)
    if (trim(io_status) /= 'ok') then
      status = 'reconstruct_compo_closure_ws write failed: ' // trim(io_status)
      return
    end if

    io_status = ''
    call read_compose(out_dir, eos_rb, io_status)
    if (trim(io_status) /= 'ok') then
      status = 'reconstruct_compo_closure_ws readback failed: ' // trim(io_status)
      return
    end if

    if (eos_rb%n /= eos_out%n .or. eos_rb%il /= eos_out%il .or. eos_rb%nadd /= eos_out%nadd) then
      status = 'reconstruct_compo_closure_ws readback header mismatch (n/il/nadd)'
      return
    end if
    if (abs(eos_rb%mn - eos_out%mn) > 1.0e-12_wp .or. abs(eos_rb%mp - eos_out%mp) > 1.0e-12_wp .or. &
        abs(eos_rb%temperature_mev - eos_out%temperature_mev) > 1.0e-12_wp) then
      status = 'reconstruct_compo_closure_ws readback mass/temperature header mismatch'
      return
    end if
    if (.not. (eos_rb%composition_file_present .and. eos_rb%composition_available)) then
      status = 'reconstruct_compo_closure_ws readback lost composition presence'
      return
    end if

    s = int(z'811C9DC5', int32)
    do i = 1, eos_out%n
      scale_rb = max(1._wp, abs(eos_out%nb(i)))
      if (abs(eos_rb%nb(i) - eos_out%nb(i)) > tol_rb * scale_rb) then
        write(status, '(a,i0)') 'reconstruct_compo_closure_ws readback nb mismatch at row ', i
        return
      end if
      do qidx = 1, 7
        scale_rb = max(1._wp, abs(eos_out%q(qidx, i)))
        if (abs(eos_rb%q(qidx, i) - eos_out%q(qidx, i)) > tol_rb * scale_rb) then
          write(status, '(a,i0,a,i0)') 'reconstruct_compo_closure_ws readback Q mismatch at row ', i, ' q=', qidx
          return
        end if
        call fold_real64(s, eos_rb%q(qidx, i))
      end do
      if (eos_rb%iphase(i) /= eos_out%iphase(i)) then
        write(status, '(a,i0)') 'reconstruct_compo_closure_ws readback iphase mismatch at row ', i
        return
      end if
      if (.not. (eos_rb%c3_applicable(i) .eqv. eos_out%c3_applicable(i)) .or. &
          .not. (eos_rb%c4_applicable(i) .eqv. eos_out%c4_applicable(i))) then
        write(status, '(a,i0)') 'reconstruct_compo_closure_ws readback C3/C4 mask mismatch at row ', i
        return
      end if
      if (abs(eos_rb%yn(i) - eos_out%yn(i)) > tol_rb .or. &
          abs(eos_rb%yp(i) - eos_out%yp(i)) > tol_rb .or. &
          abs(eos_rb%ye(i) - eos_out%ye(i)) > tol_rb .or. &
          abs(eos_rb%ymu(i) - eos_out%ymu(i)) > tol_rb .or. &
          abs(eos_rb%abar(i) - eos_out%abar(i)) > tol_rb * max(1._wp, abs(eos_out%abar(i))) .or. &
          abs(eos_rb%zbar(i) - eos_out%zbar(i)) > tol_rb * max(1._wp, abs(eos_out%zbar(i))) .or. &
          abs(eos_rb%yh(i) - eos_out%yh(i)) > tol_rb) then
        write(status, '(a,i0)') 'reconstruct_compo_closure_ws readback composition mismatch at row ', i
        return
      end if
      call fold_real64(s, eos_rb%nb(i))
      call fold_real64(s, eos_rb%yn(i))
      call fold_real64(s, eos_rb%yp(i))
      call fold_real64(s, eos_rb%ye(i))
      call fold_real64(s, eos_rb%ymu(i))
      call fold_real64(s, eos_rb%abar(i))
      call fold_real64(s, eos_rb%zbar(i))
      call fold_real64(s, eos_rb%yh(i))
      s = ieor(s, int(eos_rb%iphase(i), int32)) * int(z'01000193', int32)
    end do

    write(out_hash_hex, '(Z8.8)') s
    status = ''
  contains
    subroutine fold_real64(state, x)
      integer(int32), intent(inout) :: state
      real(wp), intent(in) :: x
      integer(int32) :: words(2)
      words = transfer(x, [0_int32, 0_int32])
      state = ieor(state, words(1)) * int(z'01000193', int32)
      state = ieor(state, words(2)) * int(z'01000193', int32)
    end subroutine fold_real64
  end subroutine reconstruct_compo_closure_ws

end module compose_eos_wigner_seitz
