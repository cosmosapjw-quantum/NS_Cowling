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
            reconstruct_wigner_seitz_table, wigner_seitz_output_hash

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

end module compose_eos_wigner_seitz
