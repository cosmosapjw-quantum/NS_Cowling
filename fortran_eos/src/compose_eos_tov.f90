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
  public :: solve_tov_love_maxwell, solve_tov_love_layered, solve_tov_layered_target
  ! Read-only low-level scientific validation seam for tidal RHS.
  public :: tov_rhs

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
    integer :: transition_count = 0  ! count crossed by the layered integrator
    real(wp) :: surface_pressure = 0._wp, surface_epsilon = 0._wp, surface_jump_scale = 0._wp
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
    ! Find the strictly monotone EOS pressure bracket in O(log n).
    ! This preserves the previous right-side convention when P==P_arr(j).
    block
      integer :: lo, hi, mid
      lo = 1; hi = n
      do while(hi-lo>1)
        mid = lo + (hi-lo)/2
        if (P_mev < P_arr(mid)) then
          hi = mid
        else
          lo = mid
        endif
      enddo
      j = lo
    end block
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
  ! Validate one continuous P(epsilon) branch before any interpolation.
  ! branch_interp uses log-log interpolation and clamps out-of-range queries;
  ! callers that model a phase interface must prove the complete query range
  ! is covered before entering the integrator.
  subroutine validate_maxwell_branch(P_arr,eps_arr,n,label,ok,status)
    integer,intent(in)::n
    real(wp),intent(in)::P_arr(n),eps_arr(n)
    character(len=*),intent(in)::label
    logical,intent(out)::ok
    character(len=*),intent(out)::status
    real(wp)::alpha,cs2_left,cs2_right,log_dp,log_de
    integer::i
    ok=.false.;status=trim(label)//' branch requires at least two rows'
    if(n<2)return
    status=trim(label)//' branch contains nonfinite P or epsilon'
    if(.not.all(ieee_is_finite(P_arr)).or..not.all(ieee_is_finite(eps_arr)))return
    status=trim(label)//' branch requires positive pressure and energy density'
    if(any(P_arr<=0._wp).or.any(eps_arr<=0._wp))return
    do i=1,n-1
      status=trim(label)//' branch pressure must be strictly increasing'
      if(P_arr(i+1)<=P_arr(i))return
      status=trim(label)//' branch energy density must be strictly increasing'
      if(eps_arr(i+1)<=eps_arr(i))return
      log_dp=log(P_arr(i+1)/P_arr(i))
      log_de=log(eps_arr(i+1)/eps_arr(i))
      if(log_dp<=0._wp.or.log_de<=0._wp)return
      alpha=log_de/log_dp
      if(.not.ieee_is_finite(alpha).or.alpha<=0._wp)return
      ! For the log-log linear segment, c_s^2=P/(epsilon*alpha).
      ! P/epsilon is monotone within each segment, so endpoint checks bound
      ! the full segment without an arbitrary sampling grid.
      cs2_left=P_arr(i)/(eps_arr(i)*alpha)
      cs2_right=P_arr(i+1)/(eps_arr(i+1)*alpha)
      status=trim(label)//' branch is mechanically unstable or acausal'
      if(.not.all(ieee_is_finite([cs2_left,cs2_right])))return
      if(min(cs2_left,cs2_right)<=0._wp.or.max(cs2_left,cs2_right)>1._wp+1.e-5_wp)return
    enddo
    ok=.true.;status='ok'
  end subroutine validate_maxwell_branch

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
    ! Hinderer tidal Q contains d(epsilon)/dP = 1/c_s^2, not c_s^2.
    Q = (4._wp*pi*elam*(5._wp*eg + 9._wp*pg + (eg + pg)*dedp_v) &
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
    ! Hinderer 2008 Eq. (23): (8/5) C^5 (1-2C)^2 (...) / denominator.
    k2 = (8._wp/5._wp)*t1/den
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
                          e5=-17253._wp/339200._wp, e6=22._wp/525._wp, e7=-1._wp/40._wp
    real(wp) :: r, p, m, y, h, r2, p2, m2, y2
    real(wp) :: k1p, k1m, k1y, k2p, k2m, k2y, k3p, k3m, k3y
    real(wp) :: k4p, k4m, k4y, k5p, k5m, k5y, k6p, k6m, k6y, k7p, k7m, k7y
    real(wp) :: err, tol, h_new, errp, errm, erry, fraction
    character(len=32) :: native_integrator
    integer :: env_status
    logical :: step_ok
    integer :: i
    integer, parameter :: max_steps = 200000

    ! Corrected DVODE Adams is the production default; RK45 remains
    ! independently available with genuine trial-step rejection.
    native_integrator='dvode_adams'
    call get_environment_variable('NS_COWLING_TOV_INTEGRATOR',native_integrator,status=env_status)
    if(env_status==1)native_integrator='dvode_adams'
    if(trim(native_integrator)=='dvode_adams'.or.trim(native_integrator)=='dvode_bdf')then
      call integrate_branch_dvode(r0,p0,m0,y0,P_arr,eps_arr,n,dr_max,p_surface, &
                                  trim(native_integrator),r_out,m_out,y_out,n_steps,ok,status)
      return
    else if(trim(native_integrator)/='rk45')then
      ok=.false.;status='unsupported NS_COWLING_TOV_INTEGRATOR';return
    endif
    ok = .false.; status = ''
    r = r0; p = p0; m = m0; y = y0; n_steps = 0
    tol = 1._wp  ! normalized embedded error threshold; trial steps may be rejected
    h = min(dr_max, 5.0e-3_wp)
    do i = 1, max_steps
      ! k1
      call tov_rhs(r, p, m, y, P_arr, eps_arr, n, p_surface, k1p, k1m, k1y, step_ok)
      if (.not. step_ok) exit
      ! k2
      ! DOPRI c2=1/5, NOT 1/2: stage-time errors invalidate the tableau.
      call tov_rhs(r + a21*h, p + a21*h*k1p, m + a21*h*k1m, &
                   y + a21*h*k1y, P_arr, eps_arr, n, p_surface, k2p, k2m, k2y, step_ok)
      if (.not. step_ok) exit
      ! k3
      call tov_rhs(r + (a31+a32)*h, p + h*(a31*k1p + a32*k2p), m + h*(a31*k1m + a32*k2m), &
                   y + h*(a31*k1y + a32*k2y), P_arr, eps_arr, n, p_surface, k3p, k3m, k3y, step_ok)
      if (.not. step_ok) exit
      ! k4
      call tov_rhs(r + (a41+a42+a43)*h, p + h*(a41*k1p + a42*k2p + a43*k3p), &
                   m + h*(a41*k1m + a42*k2m + a43*k3m), &
                   y + h*(a41*k1y + a42*k2y + a43*k3y), P_arr, eps_arr, n, p_surface, &
                   k4p, k4m, k4y, step_ok)
      if (.not. step_ok) exit
      ! k5: DOPRI c5=8/9; previously evaluated at t+h in error.
      call tov_rhs(r + (a51+a52+a53+a54)*h, p + h*(a51*k1p + a52*k2p + a53*k3p + a54*k4p), &
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
      ! The true Dormand-Prince 5(4) embedded estimate includes k7.
      ! Local error is normalized by separately scaled pressure, mass and y.
      call tov_rhs(r2,p2,m2,y2,P_arr,eps_arr,n,p_surface,k7p,k7m,k7y,step_ok)
      if (.not.step_ok) then
        h=max(h/2._wp,1.e-10_wp)
        if (h<=1.e-10_wp) then
          status='RK45 stage invalid at minimum step';exit
        endif
        cycle
      endif
      errp=abs(h*(e1*k1p+e3*k3p+e4*k4p+e5*k5p+e6*k6p+e7*k7p))/ &
           (1.e-15_wp+2.e-11_wp*max(abs(p),abs(p2)))
      errm=abs(h*(e1*k1m+e3*k3m+e4*k4m+e5*k5m+e6*k6m+e7*k7m))/ &
           (1.e-13_wp+2.e-11_wp*max(abs(m),abs(m2)))
      erry=abs(h*(e1*k1y+e3*k3y+e4*k4y+e5*k5y+e6*k6y+e7*k7y))/ &
           (1.e-12_wp+2.e-11_wp*max(abs(y),abs(y2)))
      err=max(errp,errm,erry)
      if (.not.ieee_is_finite(err)) then
        status='nonfinite native RK45 local error';exit
      endif
      h_new=h*min(4._wp,max(.2_wp,.9_wp/max(err,1.e-12_wp)**.2_wp))
      if (err>tol) then
        ! Reject the trial in ALL state components. Old code accepted here.
        h=max(1.e-10_wp,min(dr_max,h_new))
        if(h<=1.e-10_wp) then
          status='RK45 tolerance unattainable at minimum step';exit
        endif
        cycle
      endif
      if(p2<=p_surface)then
        ! Interpolate the pressure-defined endpoint within the accepted step.
        ! Do not apply an O(h) endpoint overshoot to M, R or tidal y.
        fraction=max(0._wp,min(1._wp,(p-p_surface)/(p-p2)))
        m=m+fraction*(m2-m)
        y=y+fraction*(y2-y)
        r=r+fraction*h
        p=p_surface
        n_steps=n_steps+1
        ok=.true.;exit
      endif
      r=r2;p=p2;m=m2;y=y2
      n_steps=n_steps+1
      h=max(1.e-10_wp,min(dr_max,h_new))
      if (1._wp-2._wp*m/r<=0._wp) then
        status='trapped surface';ok=.false.;exit
      endif
    end do
    r_out = r; m_out = m; y_out = y
    if (i >= max_steps .and. .not. ok) status = 'max steps reached'
    if (ok .and. len_trim(status) == 0) status = 'ok'
  end subroutine integrate_branch

  ! The DVODE package has internal COMMON workspace. Each solver process uses
  ! one instance at a time. Never call concurrently from OpenMP threads.
  ! RK45 and DVODE use the same physical tov_rhs and source EOS arrays.
  subroutine dvode_rhs(neq,t,y,ydot,rpar,ipar)
    integer,intent(in)::neq
    real(wp),intent(in)::t,y(neq),rpar(*)
    real(wp),intent(out)::ydot(neq)
    integer,intent(inout)::ipar(*)
    integer::n
    logical::valid
    n=ipar(1)
    ydot=0._wp
    if(neq/=3.or.n<2)then
      ipar(2)=1
      return
    endif
    call tov_rhs(t,y(1),y(2),y(3),rpar(2:1+n),rpar(2+n:1+2*n),n,rpar(1), &
                 ydot(1),ydot(2),ydot(3),valid)
    if(.not.valid)then
      ipar(2)=1
      ydot=0._wp
    endif
  end subroutine dvode_rhs

  subroutine dvode_jac(neq,t,y,ml,mu,pd,nrowpd,rpar,ipar)
    integer,intent(in)::neq,ml,mu,nrowpd
    real(wp),intent(in)::t,y(neq),rpar(*)
    integer,intent(in)::ipar(*)
    real(wp),intent(out)::pd(nrowpd,neq)
    ! Only MF=10 (Adams) or MF=22 (BDF internal Jacobian) is supported.
    ! DVODE never calls this user callback under those method flags.
    pd=0._wp
  end subroutine dvode_jac

  subroutine integrate_branch_dvode(r0,p0,m0,y0,P_arr,eps_arr,n, &
                                    dr_max,p_surface,algorithm, &
                                    r_out,m_out,y_out,n_steps,ok,status)
    integer,intent(in)::n
    real(wp),intent(in)::r0,p0,m0,y0,P_arr(n),eps_arr(n),dr_max,p_surface
    character(len=*),intent(in)::algorithm
    real(wp),intent(out)::r_out,m_out,y_out
    integer,intent(out)::n_steps
    logical,intent(out)::ok
    character(len=*),intent(out)::status

    integer,parameter::neq=3,lrw=200,liw=70,max_steps=200000
    real(wp)::y(neq),old_y(neq),t,t_old,tout,atol(neq),rtol
    real(wp)::rwork(lrw),rpar(1+2*n),dky(neq),mid,left,right
    integer::iwork(liw),ipar(2),istate,iopt,itol,itask,mf,iflag,i
    integer::step_counter,env_status,ios
    character(len=48)::rtol_text
    external::dvode,dvindy

    ok=.false.;status='dvode not started'
    r_out=r0;m_out=m0;y_out=y0;n_steps=0
    if(n<2.or.dr_max<=0._wp.or.p_surface<=0._wp.or.p0<=p_surface)return
    select case(trim(algorithm))
    case('dvode_adams')
      mf=10
    case('dvode_bdf')
      mf=22
    case default
      status='unknown native DVODE method';return
    end select
    y=[p0,m0,y0];t=r0
    rpar(1)=p_surface
    rpar(2:1+n)=P_arr
    rpar(2+n:1+2*n)=eps_arr
    ipar=[n,0];iwork=0;rwork=0._wp
    rtol=2.e-9_wp
    rtol_text=''
    call get_environment_variable('NS_COWLING_TOV_RTOL',rtol_text,status=env_status)
    if(env_status==0)then
      read(rtol_text,*,iostat=ios)rtol
      if(ios/=0)then
        status='DVODE invalid relative tolerance text';return
      endif
      if(.not.ieee_is_finite(rtol).or.rtol<1.e-12_wp.or.rtol>1.e-5_wp)then
        status='DVODE relative tolerance out of [1e-12,1e-5]';return
      endif
    else if(env_status/=1)then
      status='DVODE environment tolerance unavailable';return
    endif
    ! Pressure, mass and tidal-y component tolerances are all scaled with
    ! the same requested relative target, not silently frozen at baseline.
    atol=[5.e-5_wp,5.e-3_wp,5.e-2_wp]*rtol
    istate=1;itol=2;itask=2;iopt=1
    rwork(6)=dr_max
    rwork(5)=min(0.005_wp,dr_max)
    iwork(6)=100000
    do i=1,max_steps
      t_old=t;old_y=y
      tout=t+dr_max
      call dvode(dvode_rhs,neq,y,t,tout,itol,rtol,atol,itask, &
          istate,iopt,rwork,lrw,iwork,liw,dvode_jac,mf,rpar,ipar)
      if(ipar(2)/=0)then
        status='DVODE physical RHS left admissible state';exit
      endif
      if(istate/=2)then
        write(status,'(a,i0)')'DVODE ISTATE failure ',istate
        exit
      endif
      if(.not.all(ieee_is_finite(y)).or..not.ieee_is_finite(t))then
        status='DVODE nonfinite state';exit
      endif
      if(t<=t_old.or.t>50._wp)then
        status='DVODE nonmonotone/out-of-range stellar radius';exit
      endif
      if(y(2)<=0._wp.or.2._wp*y(2)>=t)then
        status='DVODE trapped star or nonpositive mass';exit
      endif
      n_steps=i
      if(y(1)<=p_surface)then
        if(old_y(1)<=p_surface)then
          status='DVODE crossing was not bracketed';exit
        endif
        ! DVINDY evaluates the dense output from DVODE's *accepted* last
        ! internal step. Pressure-defined surface root, no O(h) overshoot.
        left=t_old;right=t
        do step_counter=1,60
          mid=0.5_wp*(left+right)
          call dvindy(mid,0,rwork(21),neq,dky,iflag)
          if(iflag/=0.or..not.all(ieee_is_finite(dky)))then
            status='DVODE dense surface interpolation failed'
            return
          endif
          if(dky(1)>p_surface)then
            left=mid
          else
            right=mid
          endif
          if(right-left<1.e-12_wp*max(1._wp,abs(mid)))exit
        enddo
        mid=0.5_wp*(left+right)
        call dvindy(mid,0,rwork(21),neq,dky,iflag)
        if(iflag/=0)then
          status='DVODE final surface interpolation failed';return
        endif
        if(abs(dky(1)-p_surface)>1.e-5_wp*max(p_surface,1.e-14_wp))then
          status='DVODE surface residual outside tolerance';return
        endif
        r_out=mid;m_out=dky(2);y_out=dky(3)
        ok=.true.;status='ok'
        return
      endif
      if(istate/=2)exit
    enddo
    r_out=t;m_out=y(2);y_out=y(3)
    if(i>=max_steps.and..not.ok)status='DVODE maximum step count'
  end subroutine integrate_branch_dvode


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
    real(wp) :: eps_high_t,eps_low_t,dedp_high_t,dedp_low_t,jump_from_tables,jump_scale
    integer :: n_steps1, n_steps2
    logical :: ok,high_ok,low_ok
    character(len=128) :: branch_status
    character(len=64) :: status

    result = tov_love_result_t()
    result%central_pressure = central_pressure
    result%has_transition = .true.
    result%p_transition_mevfm3 = P_transition
    result%delta_eps_mevfm3 = delta_eps
    pc = central_pressure

    ! Validate every input branch before branch_interp is allowed to clamp an
    ! out-of-range query to an endpoint. The high branch must span Pt..Pc and
    ! the low branch must span P_surface..Pt.
    if(nh<2.or.nl<2)then
      result%status='Maxwell TOV requires at least two rows in each phase branch';return
    endif
    if(.not.all(ieee_is_finite([P_transition,delta_eps,central_pressure])))then
      result%status='Maxwell TOV transition inputs must be finite';return
    endif
    if(P_transition<=0._wp.or.delta_eps<=0._wp.or.central_pressure<=P_transition)then
      result%status='Maxwell TOV requires Pt>0, delta_eps>0 and Pc>Pt';return
    endif
    call validate_maxwell_branch(high_P,high_eps,nh,'high',high_ok,branch_status)
    if(.not.high_ok)then
      result%status=trim(branch_status);return
    endif
    call validate_maxwell_branch(low_P,low_eps,nl,'low',low_ok,branch_status)
    if(.not.low_ok)then
      result%status=trim(branch_status);return
    endif
    if(P_transition<high_P(1).or.P_transition>high_P(nh).or. &
       central_pressure<high_P(1).or.central_pressure>high_P(nh))then
      result%status='Maxwell Pt and Pc must lie inside the high-density branch range';return
    endif
    if(P_transition<=low_P(1).or.P_transition>low_P(nl))then
      result%status='Maxwell Pt must lie inside the low-density branch range';return
    endif
    call branch_interp(high_P,high_eps,nh,P_transition,eps_high_t,dedp_high_t,ok)
    if(.not.ok.or.dedp_high_t<=0._wp)then
      result%status='high-density branch interpolation failed at Pt';return
    endif
    call branch_interp(low_P,low_eps,nl,P_transition,eps_low_t,dedp_low_t,ok)
    if(.not.ok.or.dedp_low_t<=0._wp)then
      result%status='low-density branch interpolation failed at Pt';return
    endif
    jump_from_tables=eps_high_t-eps_low_t
    jump_scale=max(abs(jump_from_tables),abs(delta_eps),tiny(1._wp))
    if(jump_from_tables<=0._wp)then
      result%status='Maxwell energy-density ordering must have a positive jump';return
    endif
    if(abs(delta_eps-jump_from_tables)>1.e-8_wp*jump_scale)then
      result%status='supplied delta_eps does not match the one-sided EOS jump at Pt';return
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

  ! Equilibrium l=2 tide for a cold barotrope with any number of explicitly
  ! labelled density jumps. Sampled coordinates stop at the finite source
  ! pressure; k2 and Lambda there are finite-surface proxies.
  subroutine solve_tov_love_layered(p, e, phase, pc, max_dlogp, result)
    real(wp), intent(in) :: p(:), e(:), pc, max_dlogp
    integer, intent(in) :: phase(:)
    type(tov_love_result_t), intent(out) :: result
    integer :: i, center, last_continuous, n
    real(wp) :: ratio, width, slope, xc, z(3), t, goal, h, ptop, etop, ec, dedp
    real(wp) :: p0, e0, u, dy, k1(3), k2(3), k3(3), k4(3), ps, es, pm, em, pe, ee
    real(wp), parameter :: jump_tol=32*epsilon(1._wp)
    logical :: ok

    result=tov_love_result_t()
    result%central_pressure=pc
    n=size(p)
    result%status='invalid layered table, step, or central pressure'
    if(n<2.or.size(e)/=n.or.size(phase)/=n)return
    if(.not.all(ieee_is_finite(p)).or..not.all(ieee_is_finite(e)))return
    if(.not.ieee_is_finite(pc).or..not.ieee_is_finite(max_dlogp))return
    if(max_dlogp<1.e-5_wp.or.max_dlogp>.1_wp.or.any(p<=0).or.any(e<=0))return
    result%status='central pressure outside tabulated pressure range'
    if(pc<=p(1).or.pc>p(n))return
    center=0
    last_continuous=0
    do i=1,n-1
      result%status='unordered pressure or invalid labelled density jump'
      ratio=p(i+1)/p(i)
      if(phase(i)/=phase(i+1).and.abs(ratio-1._wp)<=jump_tol)then
        if(e(i+1)<=e(i))return
        ! The central state at an interface is its limit of zero inner-core radius.
        if(pc>=min(p(i),p(i+1)).and.pc<=max(p(i),p(i+1)))then
          if(last_continuous==0)return
          center=last_continuous
        endif
      else
        if(p(i+1)<=p(i))return
        last_continuous=i
        if(center==0.and.pc<=p(i+1))center=i
      endif
    enddo
    result%status='no continuous EOS segment at center'
    if(center==0)return
    ! Use log-pressure distance x=ln(Pc/P), integration variables
    ! u=r^2, v=m/r, y=rH'/H. Near x=0, m=(4pi/3)e_c r^3.
    slope=log(e(center+1)/e(center))/log(p(center+1)/p(center))
    if(slope<=0.or..not.ieee_is_finite(slope))return
    ec=e(center)*(pc/p(center))**slope
    result%central_eps=ec
    p0=pc*mevfm3_to_km2
    e0=ec*mevfm3_to_km2
    dedp=slope*ec/pc
    xc=min(1.e-8_wp,0.25_wp*log(pc/p(center)))
    if(xc<=0)return
    u=xc*p0/(2*pi*(e0+p0)*(e0/3+p0))
    z=[u,4*pi*e0*u/3,2._wp-4*pi*(e0/3+11*p0+(e0+p0)*dedp)*u/7]
    t=xc
    do i=center,1,-1
      ratio=p(i+1)/p(i)
      if(phase(i)/=phase(i+1).and.abs(ratio-1._wp)<=jump_tol)then
        ! Pressure-corrected Takatsy--Kovacs jump, traversed outwards.
        dy=-4*pi*z(1)*((e(i+1)-e(i))*mevfm3_to_km2)/(z(2)+4*pi*z(1)*p(i)*mevfm3_to_km2)
        z(3)=z(3)+dy
        result%transition_count=result%transition_count+1
        cycle
      endif
      ptop=min(pc,p(i+1))
      slope=log(e(i+1)/e(i))/log(ratio)
      etop=e(i)*(ptop/p(i))**slope
      result%status='unstable or acausal active log-linear EOS segment'
      if(slope<=0.or..not.ieee_is_finite(slope))return
      ! In a log-linear cell c_s^2/c^2=P/(slope*epsilon) is monotone.
      if(max(p(i)/e(i),ptop/etop)>slope*(1._wp+1.e-5_wp))return
      goal=log(pc/p(i))
      width=log(pc/ptop)
      do while(t<goal)
        h=min(max_dlogp,goal-t,.2_wp*t)
        ps=p0*exp(-t)
        es=etop*mevfm3_to_km2*exp(-slope*(t-width))
        pm=ps*exp(-h/2); em=es*exp(-slope*h/2)
        pe=ps*exp(-h); ee=es*exp(-slope*h)
        k1=layer_rhs(z,ps,es,slope)
        k2=layer_rhs(z+h*k1/2,pm,em,slope)
        k3=layer_rhs(z+h*k2/2,pm,em,slope)
        k4=layer_rhs(z+h*k3,pe,ee,slope)
        z=z+h*(k1+2*k2+2*k3+k4)/6
        result%status='nonfinite or trapped layered stellar solution'
        if(.not.all(ieee_is_finite(z)).or.z(1)<=0.or.z(2)<=0.or.2*z(2)>=1)return
        result%n_steps=result%n_steps+1
        if(h>=goal-t)then
          t=goal
        else
          t=t+h
        endif
      enddo
    enddo
    result%r_km=sqrt(z(1));result%compactness=z(2)
    result%m_msun=z(2)*result%r_km/msun_geom_km
    result%y_surface=z(3)
    result%has_transition=result%transition_count>0
    result%surface_pressure=p(1);result%surface_epsilon=e(1)
    ! Magnitude of the local y correction for a drop by epsilon_cut.
    ! Report this sensitivity separately; a finite-P endpoint has no vacuum continuation.
    result%surface_jump_scale=4*pi*z(1)*e(1)*mevfm3_to_km2/(z(2)+4*pi*z(1)*p(1)*mevfm3_to_km2)
    call love_number_k2(z(2),z(3),result%k2,result%lambda_tidal,ok)
    result%status='invalid finite-surface tidal response'
    if(.not.ok.or.result%k2<=0)return
    result%converged=.true.;result%status='ok'
  end subroutine solve_tov_love_layered

  ! Find the first sampled ascending mass crossing in a caller-specified
  ! pressure range. No continuation past a failed or missing EOS interval.
  subroutine solve_tov_layered_target(p,e,phase,mass,low,high,max_dlogp,result)
    real(wp),intent(in)::p(:),e(:),mass,low,high,max_dlogp
    integer,intent(in)::phase(:)
    type(tov_love_result_t),intent(out)::result
    integer,parameter::scan_intervals=32, refinements=45
    real(wp)::previous,current,previous_mass,mid
    integer::i,j
    result=tov_love_result_t()
    result%status='invalid target mass or pressure search interval'
    if(size(p)<2.or..not.all(ieee_is_finite([mass,low,high])))return
    if(mass<=0.or.low<=p(1).or.high<=low.or.low>=p(size(p)))return
    previous=low
    call solve_tov_love_layered(p,e,phase,low,max_dlogp,result)
    if(.not.result%converged)return
    if(abs(result%m_msun-mass)<=1.e-8_wp)return
    previous_mass=result%m_msun
    do i=1,scan_intervals
      current=low*(min(high,p(size(p)))/low)**(real(i,wp)/scan_intervals)
      if(i==scan_intervals)current=min(high,p(size(p)))
      call solve_tov_love_layered(p,e,phase,current,max_dlogp,result)
      if(.not.result%converged)return
      if(previous_mass<mass.and.result%m_msun>=mass)then
        do j=1,refinements
          mid=sqrt(previous*current)
          call solve_tov_love_layered(p,e,phase,mid,max_dlogp,result)
          if(.not.result%converged)return
          if(abs(result%m_msun-mass)<=1.e-8_wp)return
          if(result%m_msun<mass)then
            previous=mid
          else
            current=mid
          endif
        enddo
        result%converged=.false.;result%status='target mass residual above tolerance'
        return
      endif
      previous_mass=result%m_msun
      previous=current
    enddo
    result%converged=.false.;result%status='no sampled ascending mass crossing in supplied interval'
  end subroutine solve_tov_layered_target

  pure function layer_rhs(z,p,e,slope) result(d)
    real(wp), intent(in) :: z(3),p,e,slope
    real(wp) :: d(3),u,v,y,den,g,k,tidal_f,tidal_q
    u=z(1);v=z(2);y=z(3)
    den=1-2*v
    g=v+4*pi*u*p
    k=p*den/((e+p)*g)
    tidal_f=(1+4*pi*u*(p-e))/den
    tidal_q=(4*pi*u*(5*e+9*p+(e+p)*(slope*e/p))-6)/den-4*(g/den)**2
    d=[2*u*k,(4*pi*u*e-v)*k,-(y*y+y*tidal_f+tidal_q)*k]
  end function layer_rhs

end module compose_eos_tov
