! E390-03: native Compressible Liquid Drop (CLD) cell reconstruction for the
! inner crust, ported from python_frontend/compose_eos/cld_cell.py.
!
! Exploratory proxy CLD functional, NOT the fitted surface/lattice functions
! of Scurto 2024 Eqs. 24-25. This module operates on the validated
! DD-RMF homogeneous-matter solver (module compose_eos_rmf):
!   E_cell = f*E(rho_I,rho_Ip) + (1-f)*E(rho_II,rho_IIp) + E_Coul + E_surf + E_L
! minimised over (f, R_d, rho_I); the cluster proton density rho_Ip is solved
! from stationarity of the full cell energy including finite-size terms; the gas is the
! standard low-density pure-neutron approximation (rho_IIp = 0).
!
! All units: densities fm^-3, energies MeV, energy densities MeV/fm^3,
! lengths fm.  This module reuses solve_ddrmf_homogeneous from compose_eos_rmf
! and therefore inherits the existing DD-RMF beta-equilibrium, pressure-
! identity and scalar-field residuals that the homogeneous solver enforces.
module compose_eos_cld
  use compose_eos_types, only: wp
  use compose_eos_rmf, only: ddrmf_model_t, ddrmf_state_t, load_ddrmf_model, &
                             solve_ddrmf_homogeneous
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: cld_cell_result_t
  public :: cld_reconstruct_model
  public :: cld_reconstruct_point

  ! Physical constants (kept in sync with compose_eos_rmf).
  real(wp), parameter :: hbarc = 197.3269804_wp
  real(wp), parameter :: pi    = acos(-1._wp)
  real(wp), parameter :: alpha_fs = 1._wp/137.036_wp
  real(wp), parameter :: e2    = alpha_fs*hbarc        ! e^2 in MeV*fm
  real(wp), parameter :: m_e   = 0.51099895_wp

  ! Variational bounds (mirror python_frontend/compose_eos/cld_cell.py).
  real(wp), parameter :: f_lo  = 1.e-4_wp
  real(wp), parameter :: f_hi  = 0.95_wp
  real(wp), parameter :: rd_lo = 2._wp
  real(wp), parameter :: rd_hi = 40._wp
  real(wp), parameter :: rhoI_lo = 0.01_wp
  real(wp), parameter :: rhoI_hi = 0.5_wp
  integer, parameter :: nm_max_iter = 500

  type :: cld_cell_result_t
    real(wp) :: rho_b = 0._wp
    real(wp) :: f = 0._wp              ! cluster volume fraction
    real(wp) :: r_d = 0._wp            ! cluster radius (fm)
    real(wp) :: rho_i = 0._wp          ! cluster baryon density
    real(wp) :: rho_ip = 0._wp        ! cluster proton density
    real(wp) :: rho_ii = 0._wp        ! gas baryon density
    real(wp) :: rho_iip = 0._wp       ! gas proton density
    real(wp) :: rho_e = 0._wp         ! electron density
    real(wp) :: yp_cell = 0._wp       ! cell-averaged proton fraction
    real(wp) :: e_cell = 0._wp        ! baryon-rest-subtracted energy density
    real(wp) :: eps_total = 0._wp     ! physical energy density, including rest masses
    real(wp) :: e_binding_per_baryon = 0._wp
    real(wp) :: e_coul = 0._wp
    real(wp) :: e_surf = 0._wp
    real(wp) :: e_l = 0._wp           ! electron energy density
    real(wp) :: mu_n = 0._wp
    real(wp) :: mu_p = 0._wp
    real(wp) :: mu_e = 0._wp
    real(wp) :: beta_resid = 0._wp
    logical  :: converged = .false.
    integer  :: n_iter = 0
    character(len=64) :: status = ''
  end type cld_cell_result_t

contains

  !------------------------------------------------------------------
  ! Proxy lattice function; NOT the Ravenhall-Pethick-Wilson sphere function.
  pure function coulomb_lattice_phi(u) result(phi)
    real(wp), intent(in) :: u
    real(wp) :: phi
    if (u <= 0._wp .or. u >= 2._wp/3._wp) then
      phi = 0._wp
    else
      phi = 0.5_wp*sqrt(2._wp - 3._wp*u)
    endif
  end function coulomb_lattice_phi

  ! Proxy surface tension (not a source-calibrated fit): sigma = sigma_0 * (yp*(1-yp))^b
  pure function surface_tension(yp_i, sigma_0, b) result(sigma)
    real(wp), intent(in) :: yp_i, sigma_0, b
    real(wp) :: sigma
    if (yp_i <= 0._wp .or. yp_i >= 1._wp) then
      sigma = 0._wp
    else
      sigma = sigma_0*(yp_i*(1._wp - yp_i))**b
    endif
  end function surface_tension

  ! Relativistic degenerate electron energy density (MeV/fm^3) at T=0.
  pure function electron_energy_density(rho_e) result(e_l)
    real(wp), intent(in) :: rho_e
    real(wp) :: e_l, kfe_mev, e, integral
    if (rho_e <= 0._wp) then
      e_l = 0._wp; return
    endif
    kfe_mev = hbarc*(3._wp*pi*pi*rho_e)**(1._wp/3._wp)
    e = sqrt(kfe_mev*kfe_mev + m_e*m_e)
    integral = 0.125_wp*(kfe_mev*(2._wp*kfe_mev*kfe_mev + m_e*m_e)*e - &
                         m_e**4*log((kfe_mev + e)/m_e))
    e_l = (2._wp/(2._wp*pi*pi))*integral/hbarc**3
  end function electron_energy_density

  ! Electron chemical potential (Fermi energy, MeV) at T=0.
  pure function electron_chemical_potential(rho_e) result(mu_e)
    real(wp), intent(in) :: rho_e
    real(wp) :: mu_e, kfe_mev
    if (rho_e <= 0._wp) then
      mu_e = m_e; return
    endif
    kfe_mev = hbarc*(3._wp*pi*pi*rho_e)**(1._wp/3._wp)
    mu_e = sqrt(kfe_mev*kfe_mev + m_e*m_e)
  end function electron_chemical_potential

  !------------------------------------------------------------------
  ! Conditional proton stationarity for the declared proxy functional.
  ! The derivative of epsilon_cell/f includes electron, Coulomb and surface
  ! terms. A finite scan brackets local minima; it is not a proof of global
  ! uniqueness for arbitrary sigma_0,b. Endpoint solutions are not admitted
  ! as beta-equilibrated cells by the residual gate below.
  subroutine proton_energy_residual(model,rho_i,f,r_d,sigma_0,b,rho_ip,energy,residual,ok)
    type(ddrmf_model_t), intent(in) :: model
    real(wp), intent(in) :: rho_i,f,r_d,sigma_0,b,rho_ip
    real(wp), intent(out) :: energy,residual
    logical, intent(out) :: ok
    type(ddrmf_state_t) :: st
    character(len=80) :: status
    real(wp) :: y, ysafe, phi, sigma, dsigma, ne
    y=rho_ip/rho_i; ysafe=max(1.e-12_wp,min(1._wp-1.e-12_wp,y))
    call solve_ddrmf_homogeneous(model,rho_i,y,st,status)
    ok=status=='ok'
    energy=huge(1._wp); residual=huge(1._wp)
    if(.not.ok) return
    phi=coulomb_lattice_phi(f); sigma=surface_tension(y,sigma_0,b)
    dsigma=sigma_0*b*(ysafe*(1._wp-ysafe))**(b-1._wp)*(1._wp-2._wp*ysafe)
    ne=f*rho_ip
    energy=st%eps_b+electron_energy_density(ne)/f+ &
      2._wp*pi*e2*phi*r_d*r_d*rho_ip*rho_ip+3._wp*sigma/r_d
    residual=st%mu_p-st%mu_n+electron_chemical_potential(ne)+ &
      4._wp*pi*e2*phi*r_d*r_d*rho_ip+3._wp*dsigma/(r_d*rho_i)
    ok=ieee_is_finite(energy).and.ieee_is_finite(residual)
  end subroutine proton_energy_residual

  subroutine solve_beta_equilibrium_cluster(model,rho_i,f,r_d,sigma_0,b,rho_ip,rho_e,mu_e,beta_resid)
    type(ddrmf_model_t), intent(in) :: model
    real(wp), intent(in) :: rho_i,f,r_d,sigma_0,b
    real(wp), intent(out) :: rho_ip,rho_e,mu_e,beta_resid
    real(wp) :: left,right,dl,dr,lo,hi,mid,dm,energy,best_energy,best_residual
    logical :: ok, ok_left
    real(wp) :: nodes(151), temp
    integer :: j,it,k,l
    do j=1,129
      nodes(j)=real(j-1,wp)/128._wp
    enddo
    do j=1,11
      nodes(129+2*j-1)=10._wp**(-j)
      nodes(129+2*j)=1._wp-10._wp**(-j)
    enddo
    do k=2,size(nodes)
      temp=nodes(k); l=k-1
      do while(l>=1)
        if(nodes(l)<=temp) exit
        nodes(l+1)=nodes(l); l=l-1
      enddo
      nodes(l+1)=temp
    enddo
    call proton_energy_residual(model,rho_i,f,r_d,sigma_0,b,0._wp,best_energy,dl,ok_left)
    rho_ip=0._wp; best_residual=huge(1._wp)
    left=0._wp
    do j=2,size(nodes)
      right=rho_i*nodes(j)
      call proton_energy_residual(model,rho_i,f,r_d,sigma_0,b,right,energy,dr,ok)
      if(ok.and.ok_left.and.dl<0._wp.and.dr>=0._wp) then
        lo=left; hi=right
        do it=1,64
          mid=0.5_wp*(lo+hi)
          call proton_energy_residual(model,rho_i,f,r_d,sigma_0,b,mid,energy,dm,ok)
          if(.not.ok) exit
          if(abs(dm)<=1.e-9_wp.or.hi-lo<=1.e-14_wp*rho_i) exit
          if(dm>0._wp) then
            hi=mid
          else
            lo=mid
          endif
        enddo
        if(ok.and.abs(dm)<1.e-6_wp.and.energy<best_energy) then
          rho_ip=mid; best_energy=energy; best_residual=abs(dm)
        endif
      endif
      left=right
      call proton_energy_residual(model,rho_i,f,r_d,sigma_0,b,left,energy,dl,ok_left)
    enddo
    ! Include the other composition endpoint in the conditional energy comparison.
    if(ok_left.and.energy<best_energy) then
      rho_ip=rho_i; best_residual=huge(1._wp)
    endif
    rho_e=f*rho_ip; mu_e=electron_chemical_potential(rho_e)
    beta_resid=best_residual
  end subroutine solve_beta_equilibrium_cluster

  !------------------------------------------------------------------
  ! CLD total energy density (including baryon rest masses) at fixed rho_b and variational point
  ! (f, R_d, rho_i).  rho_ii is set by baryon conservation; rho_ip by beta-
  ! equilibrium at the cluster; rho_iip = 0 (pure neutron gas).
  function cld_cell_energy(model, rho_b, f, r_d, rho_i, sigma_0, b, &
                            rho_ip, rho_ii, rho_e, e_coul, e_surf, e_l) result(e_cell)
    type(ddrmf_model_t), intent(in) :: model
    real(wp), intent(in)  :: rho_b, f, r_d, rho_i, sigma_0, b
    real(wp), intent(out) :: rho_ip, rho_ii, rho_e, e_coul, e_surf, e_l
    real(wp) :: e_cell
    real(wp) :: rho_iip, phi, yp_i, sigma, rho_e_tmp, mu_e, beta_resid
    type(ddrmf_state_t) :: st_i, st_ii
    character(len=80) :: st_status
    real(wp) :: e_i, e_ii

    rho_iip = 0._wp
    if (abs(f - 1._wp) < 1.e-12_wp) then
      rho_ii = 0._wp
    else
      rho_ii = (rho_b - f*rho_i)/(1._wp - f)
    endif
    rho_ii = max(0._wp, rho_ii)

    ! Full proton variation, including finite-size terms.
    call solve_beta_equilibrium_cluster(model,rho_i,f,r_d,sigma_0,b,rho_ip,rho_e_tmp,mu_e,beta_resid)
    if(beta_resid>1.e-6_wp) then
      e_cell=huge(1._wp); return
    endif
    rho_e = max(0._wp, rho_e_tmp)

    ! Native eps_b is TOTAL baryonic energy, including constituent rest masses.
    call solve_ddrmf_homogeneous(model, rho_i, rho_ip/rho_i, st_i, st_status)
    if (st_status /= 'ok') then
      e_cell = huge(1._wp); return
    endif
    e_i = st_i%eps_b
    if (rho_ii > 0._wp) then
      call solve_ddrmf_homogeneous(model, rho_ii, 0._wp, st_ii, st_status)
      if (st_status /= 'ok') then
        e_cell = huge(1._wp); return
      endif
      e_ii = st_ii%eps_b
    else
      e_ii = 0._wp
    endif

    ! Proxy Coulomb term: E_Coul = 2*pi*e^2 * f * Phi * R_d^2 * (rho_ip - rho_iip)^2
    phi = coulomb_lattice_phi(f)
    e_coul = 2._wp*pi*e2*f*phi*r_d*r_d*(rho_ip - rho_iip)**2
    ! Proxy surface term: E_surf = 3*f*sigma/R_d
    yp_i = rho_ip/rho_i
    sigma = surface_tension(yp_i, sigma_0, b)
    if (r_d > 0._wp) then
      e_surf = 3._wp*f*sigma/r_d
    else
      e_surf = 0._wp
    endif
    ! electron background
    e_l = electron_energy_density(rho_e)

    e_cell = f*e_i + (1._wp - f)*e_ii + e_coul + e_surf + e_l
    if (.not. ieee_is_finite(e_cell)) e_cell = huge(1._wp)
  end function cld_cell_energy

  !------------------------------------------------------------------
  ! Objective for the Nelder-Mead: x = (f, log R_d, rho_i).  Mirrors the
  ! python _cld_objective clamp + feasibility logic.
  subroutine cld_objective(model, rho_b, sigma_0, b, x, e_cell, info_ok)
    type(ddrmf_model_t), intent(in) :: model
    real(wp), intent(in)  :: rho_b, sigma_0, b, x(3)
    real(wp), intent(out) :: e_cell
    logical, intent(out) :: info_ok
    real(wp) :: f, log_rd, r_d, rho_i, rho_ii, rho_ip, rho_e, e_coul, e_surf, e_l

    info_ok = .false.
    e_cell = huge(1._wp)
    f = max(f_lo, min(f_hi, x(1)))
    r_d = exp(max(log(rd_lo), min(log(rd_hi), x(2))))
    rho_i = max(rhoI_lo, min(rhoI_hi, x(3)))
    if (f < 1._wp) then
      rho_ii = (rho_b - f*rho_i)/(1._wp - f)
    else
      rho_ii = 0._wp
    endif
    if (rho_ii < 0._wp .or. rho_ii > rho_i) return
    e_cell = cld_cell_energy(model, rho_b, f, r_d, rho_i, sigma_0, b, &
                             rho_ip, rho_ii, rho_e, e_coul, e_surf, e_l)
    if (e_cell >= huge(1._wp)) return
    info_ok = .true.
  end subroutine cld_objective

  !------------------------------------------------------------------
  ! Lightweight Nelder-Mead in 3 variables, no external dependency.
  subroutine nelder_mead(model, rho_b, sigma_0, b, x0, x_best, e_best, n_iter)
    type(ddrmf_model_t), intent(in) :: model
    real(wp), intent(in)  :: rho_b, sigma_0, b, x0(3)
    real(wp), intent(out) :: x_best(3), e_best
    integer, intent(out) :: n_iter
    integer, parameter :: max_iter = nm_max_iter
    real(wp), parameter :: tol = 1.e-10_wp, x_tol = 1.e-6_wp
    real(wp), parameter :: alpha = 1._wp, gamma = 2._wp, rho_c = 0.5_wp, sigma_c = 0.5_wp
    real(wp) :: simplex(4, 3), fvals(4), centroid(3), xr(3), xe(3), xc(3)
    real(wp) :: fr, fe, fc, steps(3)
    logical :: ok(4), okr, oke, okc
    integer :: it, i, j, order(4)
    real(wp) :: e_tmp

    steps = [0.02_wp, log(1.5_wp), 0.02_wp*model%rho_sat]
    simplex(1, :) = x0
    do i = 1, 3
      simplex(i+1, :) = x0
      simplex(i+1, i) = x0(i) + steps(i)
    end do
    do i = 1, 4
      call cld_objective(model, rho_b, sigma_0, b, simplex(i, :), fvals(i), ok(i))
      if (.not. ok(i)) fvals(i) = huge(1._wp)
    end do

    do it = 1, max_iter
      ! sort vertices by ascending f
      order = sort_indices(fvals)
      simplex = simplex(order, :)
      fvals = fvals(order)
      ok = ok(order)
      if (abs(fvals(4) - fvals(1)) < tol .and. fvals(1) < huge(1._wp) &
          .and.maxval(abs(simplex(2:4,:)-spread(simplex(1,:),1,3)))<x_tol) exit
      if (fvals(1) >= huge(1._wp)) exit   ! no feasible point at all
      centroid = 0._wp
      do j = 1, 3
        centroid = centroid + simplex(j, :)
      end do
      centroid = centroid/3._wp
      ! reflection
      xr = centroid + alpha*(centroid - simplex(4, :))
      call cld_objective(model, rho_b, sigma_0, b, xr, fr, okr)
      if (.not. okr) fr = huge(1._wp)
      if (fvals(1) <= fr .and. fr < fvals(3)) then
        simplex(4, :) = xr; fvals(4) = fr; ok(4) = okr
        cycle
      endif
      if (fr < fvals(1)) then
        xe = centroid + gamma*(xr - centroid)
        call cld_objective(model, rho_b, sigma_0, b, xe, fe, oke)
        if (.not. oke) fe = huge(1._wp)
        if (fe < fr) then
          simplex(4, :) = xe; fvals(4) = fe; ok(4) = oke
        else
          simplex(4, :) = xr; fvals(4) = fr; ok(4) = okr
        endif
        cycle
      endif
      ! contraction
      xc = centroid + rho_c*(simplex(4, :) - centroid)
      call cld_objective(model, rho_b, sigma_0, b, xc, fc, okc)
      if (.not. okc) fc = huge(1._wp)
      if (fc < fvals(4)) then
        simplex(4, :) = xc; fvals(4) = fc; ok(4) = okc
        cycle
      endif
      ! shrink toward the best vertex
      do j = 2, 4
        simplex(j, :) = simplex(1, :) + sigma_c*(simplex(j, :) - simplex(1, :))
        call cld_objective(model, rho_b, sigma_0, b, simplex(j, :), fvals(j), ok(j))
        if (.not. ok(j)) fvals(j) = huge(1._wp)
      end do
    end do

    order = sort_indices(fvals)
    x_best = simplex(order(1), :)
    e_best = fvals(order(1))
    n_iter = it
  end subroutine nelder_mead

  pure function sort_indices(v) result(idx)
    real(wp), intent(in) :: v(4)
    integer :: idx(4), i, j, tmp
    do i = 1, 4; idx(i) = i; end do
    do i = 1, 3
      do j = i+1, 4
        if (v(idx(j)) < v(idx(i))) then
          tmp = idx(i); idx(i) = idx(j); idx(j) = tmp
        endif
      end do
    end do
  end function sort_indices

  !------------------------------------------------------------------
  ! Minimise the CLD cell energy at a single total baryon density rho_b.
  ! Two-stage: coarse grid scan, then Nelder-Mead refinement.
  subroutine cld_reconstruct_point(model_id, rho_b, sigma_0, b, result)
    character(len=*), intent(in) :: model_id
    real(wp), intent(in)  :: rho_b, sigma_0, b
    type(cld_cell_result_t), intent(out) :: result
    type(ddrmf_model_t) :: model
    type(ddrmf_state_t) :: st_i
    character(len=80) :: status
    real(wp) :: rho_sat
    real(wp) :: f_grid(7), rd_grid(5), rhoi_grid(5)
    real(wp) :: best_e, best_x(3), x0(3), e_grid, e_nm, x_nm(3)
    real(wp) :: f, r_d, rho_i, rho_ip, rho_ii, rho_e, e_coul, e_surf, e_l
    real(wp), parameter :: rho_min = 0.001_wp, rho_max = 0.08_wp
    integer :: n_iter_grid, n_iter_nm, if_, ir, iri
    logical :: ok_grid

    result = cld_cell_result_t()
    result%rho_b = rho_b
    if(.not.all(ieee_is_finite([rho_b,sigma_0,b])).or.rho_b<=0._wp.or.sigma_0<0._wp.or.b<=0._wp) then
      result%status='invalid density or surface parameters'; return
    endif
    call load_ddrmf_model(trim(model_id), model, status)
    if (status /= 'ok') then
      result%status = 'model load failed'; return
    endif
    rho_sat = model%rho_sat
    f_grid    = [0.005_wp, 0.01_wp, 0.03_wp, 0.08_wp, 0.15_wp, 0.25_wp, 0.4_wp]
    rd_grid   = [4._wp, 6._wp, 9._wp, 13._wp, 20._wp]
    rhoi_grid = [0.5_wp*rho_sat, 0.7_wp*rho_sat, 0.9_wp*rho_sat, &
                 1.1_wp*rho_sat, 1.3_wp*rho_sat]

    best_e = huge(1._wp); best_x = [0._wp, 0._wp, 0._wp]
    do if_ = 1, size(f_grid)
      do ir = 1, size(rd_grid)
        do iri = 1, size(rhoi_grid)
          x0 = [f_grid(if_), log(rd_grid(ir)), rhoi_grid(iri)]
          call cld_objective(model, rho_b, sigma_0, b, x0, e_grid, ok_grid)
          if (ok_grid .and. e_grid < best_e) then
            best_e = e_grid; best_x = x0
          endif
        end do
      end do
    end do
    if (best_e >= huge(1._wp)) then
      result%status = 'no feasible point'; result%n_iter = 0; return
    endif

    call nelder_mead(model, rho_b, sigma_0, b, best_x, x_nm, e_nm, n_iter_nm)
    if (e_nm < best_e) then
      best_x = x_nm; best_e = e_nm
    endif
    n_iter_grid = 0  ! grid iteration count not tracked separately (mirror python)
    result%n_iter = n_iter_nm

    ! extract final values with the same clamps
    f = max(f_lo, min(f_hi, best_x(1)))
    r_d = exp(max(log(rd_lo), min(log(rd_hi), best_x(2))))
    rho_i = max(rhoI_lo, min(rhoI_hi, best_x(3)))
    result%eps_total = cld_cell_energy(model, rho_b, f, r_d, rho_i, sigma_0, b, &
                                    rho_ip, rho_ii, rho_e, e_coul, e_surf, e_l)
    result%f = f; result%r_d = r_d; result%rho_i = rho_i; result%rho_ip = rho_ip
    result%rho_ii = rho_ii; result%rho_iip = 0._wp; result%rho_e = rho_e
    result%e_cell=result%eps_total-model%mn*(rho_b-rho_e)-model%mp*rho_e
    result%e_coul = e_coul; result%e_surf = e_surf; result%e_l = e_l
    if (rho_b > 0._wp) then
      result%e_binding_per_baryon = result%e_cell/rho_b
      result%yp_cell = (f*rho_ip + (1._wp - f)*0._wp)/rho_b
    endif
    ! cell mu_n/mu_p from the cluster phase
    call solve_ddrmf_homogeneous(model, rho_i, rho_ip/rho_i, st_i, status)
    if (status == 'ok') then
      result%mu_n = st_i%mu_n; result%mu_p = st_i%mu_p
    endif
    call solve_beta_equilibrium_cluster(model,rho_i,f,r_d,sigma_0,b,rho_ip,rho_e,result%mu_e,result%beta_resid)
    result%converged = ieee_is_finite(result%eps_total).and.result%eps_total<huge(1._wp) &
      .and.result%beta_resid<=1.e-6_wp.and.n_iter_nm<=nm_max_iter
    result%status = 'ok'
    if(.not.result%converged) result%status='proxy optimizer or stationarity not converged'
  end subroutine cld_reconstruct_point

  !------------------------------------------------------------------
  ! Reconstruct one model over a log-spaced inner-crust density grid and emit
  ! a JSON array (one object per density point) to stdout.  This is the
  ! exploratory entry point used by the cld-reconstruct CLI subcommand.
  subroutine cld_reconstruct_model(model_id, sigma_0, b, rho_min, rho_max, n_points)
    character(len=*), intent(in) :: model_id
    real(wp), intent(in) :: sigma_0, b, rho_min, rho_max
    integer, intent(in) :: n_points
    integer :: i, n_eff
    real(wp) :: rho_b, log_min, log_max
    type(cld_cell_result_t) :: r

    n_eff = max(1, n_points)
    log_min = log(rho_min); log_max = log(rho_max)
    write(*, '(a)', advance='no') '['
    do i = 1, n_eff
      if (n_eff == 1) then
        rho_b = rho_min
      else
        rho_b = exp(log_min + (log_max - log_min)*real(i-1, wp)/real(n_eff-1, wp))
      endif
      call cld_reconstruct_point(model_id, rho_b, sigma_0, b, r)
      if (i > 1) write(*, '(a)', advance='no') ','
      call emit_point_json(model_id, r)
    end do
    write(*, '(a)') ']'
  end subroutine cld_reconstruct_model

  subroutine emit_point_json(model_id, r)
    character(len=*), intent(in) :: model_id
    type(cld_cell_result_t), intent(in) :: r
    write(*, '(a)', advance='no') &
      '{"model":"'//trim(model_id)// &
      '","rho_b_fm3":'//trim(fmt(r%rho_b))// &
      ',"f":'//trim(fmt(r%f))// &
      ',"R_d_fm":'//trim(fmt(r%r_d))// &
      ',"rho_I_fm3":'//trim(fmt(r%rho_i))// &
      ',"rho_Ip_fm3":'//trim(fmt(r%rho_ip))// &
      ',"rho_II_fm3":'//trim(fmt(r%rho_ii))// &
      ',"rho_IIp_fm3":'//trim(fmt(r%rho_iip))// &
      ',"rho_e_fm3":'//trim(fmt(r%rho_e))// &
      ',"yp_cell":'//trim(fmt(r%yp_cell))// &
      ',"E_cell_MeV_fm3":'//trim(fmt(r%e_cell))// &
      ',"epsilon_total_MeV_fm3":'//trim(fmt(r%eps_total))// &
      ',"physical_scope":"exploratory_proxy_CLD","source_calibrated":false'// &
      ',"E_binding_per_baryon_MeV":'//trim(fmt(r%e_binding_per_baryon))// &
      ',"E_coul_MeV_fm3":'//trim(fmt(r%e_coul))// &
      ',"E_surf_MeV_fm3":'//trim(fmt(r%e_surf))// &
      ',"E_L_MeV_fm3":'//trim(fmt(r%e_l))// &
      ',"mu_n_MeV":'//trim(fmt(r%mu_n))// &
      ',"mu_p_MeV":'//trim(fmt(r%mu_p))// &
      ',"mu_e_MeV":'//trim(fmt(r%mu_e))// &
      ',"beta_residual_MeV":'//trim(fmt(r%beta_resid))// &
      ',"converged":'//merge('true ','false',r%converged)// &
      ',"n_iter":'//trim(int_fmt(r%n_iter))// &
      ',"status":"'//trim(r%status)//'"}'
  end subroutine emit_point_json

  ! Compact real-to-string (no leading/trailing whitespace, scientific).
  pure function fmt(x) result(s)
    real(wp), intent(in) :: x
    character(len=32) :: s, tmp
    if (.not. ieee_is_finite(x)) then
      if (x > 0._wp) then
        s = '1.0e99'
      else
        s = 'null'
      endif
      return
    endif
    write(tmp, '(es24.16e3)') x
    s = adjustl(tmp)
  end function fmt

  pure function int_fmt(i) result(s)
    integer, intent(in) :: i
    character(len=16) :: s, tmp
    write(tmp, '(i0)') i
    s = adjustl(tmp)
  end function int_fmt

end module compose_eos_cld
