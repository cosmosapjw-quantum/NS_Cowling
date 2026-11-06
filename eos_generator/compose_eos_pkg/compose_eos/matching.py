"""Core--crust matching for cold, one-dimensional CompOSE barotropes.

The module constructs smooth single-valued bridges in either the canonical
F_b(n_b) or grand-canonical P(mu_b) representation. The historical method token
``maxwell`` is retained for API compatibility, but currently denotes a smooth
near-pressure-plateau bridge: it does NOT construct an exact Maxwell coexistence
point, latent-heat density jump, or a Seidov-applicable first-order transition.
The other methods are interpolation regularity choices, not phase-transition
orders. ``find_convex_transition_window`` is a candidate-screening heuristic;
all completed EOS tables must pass independent consistency and causality gates.

Composition is propagated only when the core input passes the CompOSE baryon and
charge closure checks. By default, missing/unusable core composition remains
absent. An explicit ``literature_beta_eq`` option supplies a model-dependent
parabolic symmetry-energy approximation for registered APR/BBB core models only.

The TOV helper in this module supports strictly monotone smooth P(epsilon) tables.
It rejects out-of-range central pressures and missed pressure events. It does not
implement the tidal jump condition for a genuine first-order density discontinuity.
Target-mass observables are evaluated at an actual bracketed root of M(P_c), never
manufactured by linear interpolation of neighbouring stellar observables.
"""
from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.optimize import brentq

from .eos_data import EOSData
from .consistency import check_consistency

HBARC = 197.3269804
ME = 0.51099895
MMU = 105.6583755
R_SCHW_SUN = 1.4766250385079258
CONV_MEV_FM3 = 8.96526390e-7
FOUR_PI = 4.0 * np.pi

METHODS = ("maxwell", "css", "cubic_F", "cubic_P", "quintic_F", "quintic_P")

CRUST_NAMES = (
    "SPG_M1", "SPG_M2", "SPG_M3", "SPG_M4", "SPG_M5",
    "GPPVA_DD2", "GPPVA_DDME2", "GPPVA_FSU2", "GPPVA_FSU2H",
    "GPPVA_FSU2R", "GPPVA_NL3wrL55", "GPPVA_TM1e", "GPPVA_TW",
)

CORE_NAMES = ("APR_APR", "BBB_BHF2", "GM_GM1", "RG_SLY4", "RG_KDE0v1")


def get_core_start_index(co_name: str, d_co: EOSData) -> int:
    """Return the first index of the pure homogeneous core in ``d_co``."""
    if d_co.composition_available and d_co.Aav.shape == d_co.nb.shape:
        hom = np.where(d_co.is_homogeneous & (d_co.nb >= 0.03))[0]
        if hom.size:
            return int(hom[0])
    idx = np.flatnonzero(d_co.nb >= 0.045)
    if not idx.size:
        raise ValueError("core EOS has no point at or above the homogeneous-core density threshold 0.045 fm^-3")
    return int(idx[0])


def get_gauge_shifts(
    cr_name: str, co_name: str, d_cr: EOSData, d_co: EOSData
) -> Tuple[float, float]:
    """Return ``(dm_crust, dm_core)`` (in MeV) to align both tables to the
    common nucleon rest-mass zero point ``d_co.mn``.
    """
    # A gauge shift may align different rest-mass conventions; it must not be
    # used as an empirical fit parameter. Historical +4.96005 MeV (GM1 core)
    # and +2.45 MeV (FSU2 crust) offsets are deliberately no longer applied:
    # their provenance is not established by the source data sheets.
    dm_cr = float(d_co.mn - d_cr.mn)
    dm_co = 0.0
    return dm_cr, dm_co


def find_convex_transition_window(
    d_cr: EOSData, d_co: EOSData, cr_name: str = "", co_name: str = ""
) -> Tuple[int, int]:
    """Select a candidate transition window using a two-plane secant heuristic.

    The secant inequalities are a screening heuristic only. They are not a
    sufficient theorem for convexity, causality, or stability of the completed
    bridge; every returned candidate must pass the downstream numerical gates.
    """
    icore = get_core_start_index(co_name, d_co)
    dm_cr, dm_co = get_gauge_shifts(cr_name, co_name, d_cr, d_co)
    mub_cr = d_cr.mub + dm_cr
    F_cr = d_cr.F + dm_cr
    mub_co = d_co.mub + dm_co
    F_co = d_co.F + dm_co

    for i1 in range(len(d_cr.nb) - 1, max(10, len(d_cr.nb) - 150), -1):
        for i2 in range(icore, len(d_co.nb)):
            n1, p1, m1, f1 = float(d_cr.nb[i1]), float(d_cr.p[i1]), float(mub_cr[i1]), float(F_cr[i1])
            n2, p2, m2, f2 = float(d_co.nb[i2]), float(d_co.p[i2]), float(mub_co[i2]), float(F_co[i2])
            if n2 > n1 * 1.10 and p2 > p1 * 1.10 and m2 > m1 + 0.2 and f2 > f1 + 0.2:
                n_bar = (p2 - p1) / (m2 - m1)
                q1_bar = (f2 - f1) / np.log(n2 / n1)
                th_p = (n_bar - n1) / (n2 - n1)
                th_f = (q1_bar - (p1 / n1)) / ((p2 / n2) - (p1 / n1))
                if 0.12 < th_p < 0.88 and 0.12 < th_f < 0.88:
                    return i1, i2
    raise RuntimeError(f"No valid convex transition window found for ({cr_name}, {co_name})")


def build_monotone_hermite_bridge(
    y0: float,
    y1: float,
    dy0: float,
    dy1: float,
    target_avg: float,
    match_slopes: bool,
    u_grid: np.ndarray,
    *,
    canonical_dx: float | None = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, float, str]:
    """Construct a smooth, spike-free Hermite bridge :math:`y(u)` on :math:`u \\in [0, 1]`
    satisfying :math:`y(0)=y_0`, :math:`y(1)=y_1`, and
    :math:`\\int_0^1 y(u)\\,du = \\text{target\\_avg}` to machine precision.

    - When ``match_slopes=False`` (3rd-order :math:`C^1` Hermite potential):
      Uses the exact single-segment cubic Hermite potential polynomial
      (quadratic :math:`y(u)`) whenever :math:`\\theta \\in (1/3, 2/3)` lies inside
      the Hyman--Fritsch--Carlson monotone convexity cone, and a 1-knot monotone
      cubic Hermite spline otherwise.
    - When ``match_slopes=True`` (5th-order :math:`C^2` Hermite potential):
      Preserves the **exact unclipped** boundary derivatives :math:`y'(0) = dy_0`
      and :math:`y'(1) = dy_1` (so :math:`\\Delta c_{s,\\mathrm{L}}^2 = 0` and
      :math:`\\Delta c_{s,\\mathrm{R}}^2 = 0` to machine precision :math:`\\sim 10^{-16}`),
      using the exact single-segment 5th-order Hermite potential polynomial
      (quartic :math:`y(u)`) when strictly convex and spike-free, or an adaptive
      boundary-layer 2-knot :math:`C^2` Hermite spline otherwise.
    """
    dy_tot = y1 - y0
    theta = (target_avg - y0) / dy_tot
    u_test = np.linspace(0.0, 1.0, 401)

    if not match_slopes:
        # Check if exact single-segment cubic Hermite potential Phi_3(u)
        # (quadratic y_2(u) = Phi_3'(u)) is strictly convex on [0, 1]:
        y_quad_test = (
            y0 * (1.0 - u_test) * (1.0 - 3.0 * u_test)
            + y1 * u_test * (3.0 * u_test - 2.0)
            + 6.0 * target_avg * u_test * (1.0 - u_test)
        )
        dy_quad_test = (
            y0 * (-4.0 + 6.0 * u_test)
            + y1 * (6.0 * u_test - 2.0)
            + 6.0 * target_avg * (1.0 - 2.0 * u_test)
        )
        if canonical_dx is not None:
            convex_ok = bool(np.min(y_quad_test + dy_quad_test / canonical_dx) >= 0.15 * y0 and (0.28 <= theta <= 0.72))
        else:
            convex_ok = bool(np.min(dy_quad_test) >= 0.15 * dy_tot and (1.0 / 3.0 < theta < 2.0 / 3.0))

        if convex_ok:
            y_u = (
                y0 * (1.0 - u_grid) * (1.0 - 3.0 * u_grid)
                + y1 * u_grid * (3.0 * u_grid - 2.0)
                + 6.0 * target_avg * u_grid * (1.0 - u_grid)
            )
            dy_u = (
                y0 * (-4.0 + 6.0 * u_grid)
                + y1 * (6.0 * u_grid - 2.0)
                + 6.0 * target_avg * (1.0 - 2.0 * u_grid)
            )
            int_u = (
                y0 * u_grid * (1.0 - u_grid) ** 2
                + y1 * (u_grid ** 2) * (u_grid - 1.0)
                + target_avg * (u_grid ** 2) * (3.0 - 2.0 * u_grid)
            )
            d0_exact = float(-4.0 * y0 - 2.0 * y1 + 6.0 * target_avg)
            d1_exact = float(2.0 * y0 + 4.0 * y1 - 6.0 * target_avg)
            return y_u, dy_u, int_u, d0_exact, d1_exact, "exact_single_cubic"

        um = float(np.clip(1.0 - theta, 0.18, 0.82))

        def _make_c1_sp(ym: float) -> CubicHermiteSpline:
            sec1 = (ym - y0) / um
            sec2 = (y1 - ym) / (1.0 - um)
            s0 = 0.85 * sec1
            s1 = 0.85 * sec2
            sm = (2.0 * sec1 * sec2 / max(1e-12, sec1 + sec2)) if (sec1 > 0 and sec2 > 0) else 0.5 * (sec1 + sec2)
            return CubicHermiteSpline([0.0, um, 1.0], [y0, ym, y1], [s0, sm, s1])

        ym_sol = brentq(
            lambda ym: float(_make_c1_sp(ym).integrate(0.0, 1.0) - target_avg),
            y0 + 1e-5 * dy_tot, y1 - 1e-5 * dy_tot, xtol=1e-14
        )
        sp = _make_c1_sp(ym_sol)
        anti = sp.antiderivative()
        return (
            sp(u_grid),
            sp.derivative()(u_grid),
            anti(u_grid),
            float(sp.derivative()(0.0)),
            float(sp.derivative()(1.0)),
            "monotone_1knot_cubic",
        )

    # 5th-Order / C^2 Hermite bridge with EXACT unclipped dy0 and dy1:
    # 1. Test exact single-segment 5th-order Hermite potential Phi_5(u) (quartic y_4(u)):
    c1 = y0
    c2 = 0.5 * dy0
    A0 = target_avg - c1 - c2
    A1 = y1 - c1 - 2.0 * c2
    A2 = dy1 - 2.0 * c2
    c3 = 10.0 * A0 - 4.0 * A1 + 0.5 * A2
    c4 = -15.0 * A0 + 7.0 * A1 - A2
    c5 = 6.0 * A0 - 3.0 * A1 + 0.5 * A2

    y_q_test = c1 + 2.0 * c2 * u_test + 3.0 * c3 * (u_test ** 2) + 4.0 * c4 * (u_test ** 3) + 5.0 * c5 * (u_test ** 4)
    dy_q_test = 2.0 * c2 + 6.0 * c3 * u_test + 12.0 * c4 * (u_test ** 2) + 20.0 * c5 * (u_test ** 3)
    if canonical_dx is not None:
        dpdn_test = y_q_test + dy_q_test / canonical_dx
        dpdn_0 = y0 + dy0 / canonical_dx
        dpdn_1 = y1 + dy1 / canonical_dx
        quintic_ok = bool(
            np.min(dpdn_test) >= 0.35 * min(dpdn_0, dpdn_1)
            and np.max(dpdn_test) <= 1.8 * max(dpdn_0, dpdn_1)
            and np.min(dy_q_test) >= -0.1 * dy_tot
        )
    else:
        quintic_ok = bool(
            np.min(dy_q_test) >= 0.75 * min(dy0, dy1, dy_tot)
            and np.max(dy_q_test) <= 2.2 * max(dy0, dy1, dy_tot)
        )

    if quintic_ok:
        y_u = c1 + 2.0 * c2 * u_grid + 3.0 * c3 * (u_grid ** 2) + 4.0 * c4 * (u_grid ** 3) + 5.0 * c5 * (u_grid ** 4)
        dy_u = 2.0 * c2 + 6.0 * c3 * u_grid + 12.0 * c4 * (u_grid ** 2) + 20.0 * c5 * (u_grid ** 3)
        int_u = c1 * u_grid + c2 * (u_grid ** 2) + c3 * (u_grid ** 3) + c4 * (u_grid ** 4) + c5 * (u_grid ** 5)
        return y_u, dy_u, int_u, float(dy0), float(dy1), "exact_single_quintic"

    # 2. Diagonal-knot boundary-layer C^2 Hermite spline preserving exact dy0 and dy1:
    ua = max(1e-4, min(0.02, 0.025 * dy_tot / max(abs(dy0), 1e-12)))
    ub_w = max(1e-4, min(0.02, 0.025 * dy_tot / max(abs(dy1), 1e-12)))
    ub = 1.0 - ub_w
    ya = y0 + np.clip(0.85 * max(dy0, 0.0) * ua, 0.002 * dy_tot, 0.025 * dy_tot)
    yb = y1 - np.clip(0.85 * max(dy1, 0.0) * ub_w, 0.002 * dy_tot, 0.025 * dy_tot)

    def _make_c2_sp(s: float) -> CubicHermiteSpline:
        um = ua + (0.96 - 0.92 * s) * (ub - ua)
        ym = ya + s * (yb - ya)
        knots_u = [0.0, ua, um, ub, 1.0]
        knots_y = [y0, ya, ym, yb, y1]
        secs = np.diff(knots_y) / np.diff(knots_u)
        slopes = [dy0]
        for k in range(3):
            sl, sr = secs[k], secs[k + 1]
            sh = (2.0 * sl * sr / (sl + sr)) if (sl > 0 and sr > 0) else 0.5 * (sl + sr)
            slopes.append(sh)
        slopes.append(dy1)
        return CubicHermiteSpline(knots_u, knots_y, slopes)

    s_sol = brentq(
        lambda s: float(_make_c2_sp(s).integrate(0.0, 1.0) - target_avg),
        1e-4, 1.0 - 1e-4, xtol=1e-14
    )
    sp = _make_c2_sp(s_sol)
    anti = sp.antiderivative()
    return (
        sp(u_grid),
        sp.derivative()(u_grid),
        anti(u_grid),
        float(sp.derivative()(0.0)),
        float(sp.derivative()(1.0)),
        "adaptive_2knot_c2_hermite",
    )


def build_maxwell_plateau_bridge(
    x1: float,
    x2: float,
    P1: float,
    P2: float,
    F1: float,
    F2: float,
    mu1: float,
    mu2: float,
    cs2_1: float,
    cs2_2: float,
    u_grid: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, float]:
    """Construct a smooth near-constant-pressure bridge, not a Maxwell transition.

    The interior pressure is deliberately allowed to vary by approximately
    3 percent of the endpoint pressure difference. There is no density jump,
    latent-heat discontinuity, or exact equal-pressure coexistence construction.
    The integrated first-law condition is enforced for the smooth barotrope.
    The public ``maxwell`` method token is retained only for compatibility.
    """
    Dx = x2 - x1
    x_br = x1 + u_grid * Dx
    nb_br = np.exp(x_br)
    ua, ub = 0.22, 0.78
    dP_plat = 0.015 * (P2 - P1)
    dP_du_1 = float(np.exp(x1) * mu1 * cs2_1 * Dx)
    dP_du_2 = float(np.exp(x2) * mu2 * cs2_2 * Dx)

    xg, wg = np.polynomial.legendre.leggauss(64)
    ug = 0.5 * (xg + 1.0)
    wg = 0.5 * wg
    nb_g = np.exp(x1 + ug * Dx)

    def _make_P_sp(Pt: float) -> CubicHermiteSpline:
        Pa = Pt - dP_plat
        Pb = Pt + dP_plat
        sec1 = (Pa - P1) / ua
        sec_p = (Pb - Pa) / (ub - ua)
        sec2 = (P2 - Pb) / (1.0 - ub)
        s0 = min(dP_du_1, 2.2 * max(1e-12, sec1))
        s1 = min(dP_du_2, 2.2 * max(1e-12, sec2))
        sa = min(1.5 * sec_p, 1.5 * max(1e-12, sec1))
        sb = min(1.5 * sec_p, 1.5 * max(1e-12, sec2))
        return CubicHermiteSpline([0.0, ua, ub, 1.0], [P1, Pa, Pb, P2], [s0, sa, sb, s1])

    def _int_err(Pt: float) -> float:
        sp = _make_P_sp(Pt)
        return float(Dx * np.sum(wg * (sp(ug) / nb_g)) - (F2 - F1))

    Pt_sol = brentq(_int_err, P1 + 1.05 * dP_plat, P2 - 1.05 * dP_plat, xtol=1e-14)
    sp_P = _make_P_sp(Pt_sol)
    P_br = sp_P(u_grid)
    dP_du = sp_P.derivative()(u_grid)
    Q1_br = P_br / nb_br

    F_br = np.empty_like(u_grid)
    xg16, wg16 = np.polynomial.legendre.leggauss(16)
    for k, uk in enumerate(u_grid):
        uk_g = 0.5 * uk * (xg16 + 1.0)
        wk_g = 0.5 * uk * wg16
        F_br[k] = F1 + Dx * float(np.sum(wk_g * (sp_P(uk_g) / np.exp(x1 + uk_g * Dx))))
    mu_br = F_br + Q1_br
    cs2_br = (dP_du / Dx) / (nb_br * mu_br)
    cs2_left = float((sp_P.derivative()(0.0) / Dx) / (np.exp(x1) * mu1))
    cs2_right = float((sp_P.derivative()(1.0) / Dx) / (np.exp(x2) * mu2))
    return nb_br, P_br, F_br, mu_br, Q1_br, cs2_br, float(Pt_sol), cs2_left, cs2_right


def _literature_core_symmetry_energy(co_name: str, nb_arr: np.ndarray) -> np.ndarray:
    """Return the nuclear symmetry energy S(nb) [MeV] from primary literature
    for homogeneous cores that lack an input eos.compo (APR_APR, BBB_BHF2):
      - APR_APR : Akmal, Pandharipande & Ravenhall (1998), Phys. Rev. C 58, 1804,
                  A18 + delta v + UIX* parameterization (S(0.16) = 32.6 MeV).
      - BBB_BHF2: Baldo, Bombaci & Burgio (1997), Astron. Astrophys. 328, 274,
                  Paris + TBF BHF parameterization (S(0.16) = 30.0 MeV).
    """
    base = co_name.replace("_fixed", "")
    u = nb_arr / 0.16
    if base == "APR_APR":
        return 13.0 * (u ** (2.0 / 3.0)) + 19.6 * (u ** 0.62)
    if base == "BBB_BHF2":
        return 13.0 * (u ** (2.0 / 3.0)) + 17.0 * (u ** 0.72)
    raise ValueError(
        f"No verified symmetry-energy fit is registered for {base!r}; "
        "use strict_thermo_only or supply model-specific composition."
    )


_BETA_EQ_CACHE: Dict[Tuple[str, int, float], Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}

def _solve_beta_eq_core_fractions(co_name: str, nb_co: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Approximate npe-mu beta equilibrium using a parabolic symmetry-energy fit.

    This is a model-dependent approximation, not an exact reconstruction of
    the microscopic core composition. It is available only for APR/BBB fits.
    """
    cache_key = (co_name.replace("_fixed", ""), len(nb_co), float(nb_co[0]))
    if cache_key in _BETA_EQ_CACHE:
        return _BETA_EQ_CACHE[cache_key]
    S_nb = _literature_core_symmetry_energy(co_name, nb_co)
    Yp_co = np.zeros_like(nb_co)
    Ye_co = np.zeros_like(nb_co)
    Ymu_co = np.zeros_like(nb_co)
    muq_co = np.zeros_like(nb_co)

    for k, nb_v in enumerate(nb_co):
        S_v = float(S_nb[k])

        def _lepton_split(yp_v: float) -> Tuple[float, float, float]:
            nq = nb_v * yp_v
            if nq <= 0.0:
                return 0.0, 0.0, ME
            kfe0 = HBARC * (3.0 * np.pi ** 2 * nq) ** (1.0 / 3.0)
            mue0 = float(np.sqrt(ME ** 2 + kfe0 ** 2))
            if mue0 <= MMU:
                return yp_v, 0.0, mue0
            ymu_v = float(brentq(
                lambda ym: np.sqrt(ME ** 2 + (HBARC * (3.0 * np.pi ** 2 * nb_v * max(0.0, yp_v - ym)) ** (1.0 / 3.0)) ** 2)
                - np.sqrt(MMU ** 2 + (HBARC * (3.0 * np.pi ** 2 * nb_v * ym) ** (1.0 / 3.0)) ** 2),
                0.0, yp_v, xtol=1e-14
            ))
            ye_v = yp_v - ymu_v
            kfe = HBARC * (3.0 * np.pi ** 2 * nb_v * ye_v) ** (1.0 / 3.0)
            return ye_v, ymu_v, float(np.sqrt(ME ** 2 + kfe ** 2))

        yp_sol = float(brentq(
            lambda yp: 4.0 * (1.0 - 2.0 * yp) * S_v - _lepton_split(yp)[2],
            1e-4, 0.49, xtol=1e-12
        ))
        ye_sol, ymu_sol, mue_sol = _lepton_split(yp_sol)
        Yp_co[k] = yp_sol
        Ye_co[k] = ye_sol
        Ymu_co[k] = ymu_sol
        muq_co[k] = -mue_sol

    Yn_co = 1.0 - Yp_co
    res_tup = (Yn_co, Yp_co, Ye_co, Ymu_co, muq_co)
    _BETA_EQ_CACHE[cache_key] = res_tup
    return res_tup


def _endpoint_sound_speed_squared(data: EOSData, index: int, side: str, points: int = 6) -> float:
    """Estimate c_s^2=dP/deps from one-sided PCHIP derivatives in log(n_b).

    No clipping is applied. Invalid/unstable/acausal endpoint data are rejected
    rather than silently altered to make a bridge fit.
    """
    if side == "left":
        lo = max(0, index - points + 1)
        ids = np.arange(lo, index + 1)
    elif side == "right":
        ids = np.arange(index, min(data.n, index + points))
    else:
        raise ValueError("side must be 'left' or 'right'")
    if ids.size < 3:
        raise ValueError(f"not enough points for one-sided {side} endpoint derivative")
    x = np.log(np.asarray(data.nb[ids], dtype=float))
    p = np.asarray(data.p[ids], dtype=float)
    eps = np.asarray(data.eps[ids], dtype=float)
    if np.any(p <= 0) or np.any(eps <= 0) or np.any(np.diff(x) <= 0):
        raise ValueError(f"invalid endpoint data for {side} sound-speed estimate")
    lp = PchipInterpolator(x, np.log(p), extrapolate=False)
    le = PchipInterpolator(x, np.log(eps), extrapolate=False)
    x0 = float(x[-1] if side == "left" else x[0])
    j = -1 if side == "left" else 0
    dlogp = float(lp.derivative()(x0))
    dloge = float(le.derivative()(x0))
    if not np.isfinite(dlogp + dloge) or dloge <= 0:
        raise ValueError(f"non-positive/non-finite dlog(eps)/dlog(n) at {side} endpoint")
    cs2 = float((p[j] / eps[j]) * (dlogp / dloge))
    if not np.isfinite(cs2) or cs2 <= 0.0 or cs2 >= 1.0:
        raise ValueError(f"{side} endpoint c_s^2={cs2:.8g} is outside the stable causal interval (0,1)")
    return cs2


def match_crust_core(
    d_cr: EOSData,
    d_co: EOSData,
    cr_name: str = "",
    co_name: str = "",
    method: str = "quintic_P",
    n_bridge: int = 40,
    *,
    core_compo_mode: str = "strict_thermo_only",
) -> Tuple[EOSData, Dict[str, object]]:
    """Join ``d_cr`` (crust) and ``d_co`` (core) using ``method``.

    Parameters
    ----------
    core_compo_mode : {"strict_thermo_only", "literature_beta_eq"}
        The default is ``strict_thermo_only``: absent or unusable core composition
        remains absent in the matched table. ``literature_beta_eq`` is an explicit
        opt-in parabolic-symmetry-energy approximation for APR_APR/BBB_BHF2 only.
    """
    if method not in METHODS:
        raise ValueError(f"Unknown matching method {method!r}; choose from {METHODS}")
    if core_compo_mode not in ("literature_beta_eq", "strict_thermo_only"):
        raise ValueError("core_compo_mode must be 'strict_thermo_only' or 'literature_beta_eq'")

    i1, i2 = find_convex_transition_window(d_cr, d_co, cr_name, co_name)
    mn = d_co.mn
    mp = d_co.mp
    dm_cr, dm_co = get_gauge_shifts(cr_name, co_name, d_cr, d_co)

    nb_cr = d_cr.nb[:i1 + 1]
    F_cr = d_cr.F[:i1 + 1] + dm_cr
    Q1_cr = d_cr.Q1[:i1 + 1]
    mub_cr = F_cr + Q1_cr
    P_cr = nb_cr * Q1_cr

    nb_co = d_co.nb[i2:]
    F_co = d_co.F[i2:] + dm_co
    Q1_co = d_co.Q1[i2:]
    mub_co = F_co + Q1_co
    P_co = nb_co * Q1_co

    n1, P1, F1, mu1, Q1_1 = float(nb_cr[-1]), float(P_cr[-1]), float(F_cr[-1]), float(mub_cr[-1]), float(Q1_cr[-1])
    n2, P2, F2, mu2, Q1_2 = float(nb_co[0]), float(P_co[0]), float(F_co[0]), float(mub_co[0]), float(Q1_co[0])

    cs2_1 = _endpoint_sound_speed_squared(d_cr, i1, "left")
    cs2_2 = _endpoint_sound_speed_squared(d_co, i2, "right")

    x1, x2 = np.log(n1), np.log(n2)
    Dx = x2 - x1
    Dmu = mu2 - mu1
    n_bar = (P2 - P1) / Dmu
    q1_bar = (F2 - F1) / Dx

    u_grid = np.linspace(0.0, 1.0, n_bridge + 2)[1:-1]
    extra_meta: Dict[str, object] = {}

    if method == "maxwell":
        nb_br, P_br, F_br, mu_br, Q1_br, cs2_br, Pt_sol, cs2_left_br, cs2_right_br = build_maxwell_plateau_bridge(
            x1, x2, P1, P2, F1, F2, mu1, mu2, cs2_1, cs2_2, u_grid
        )
        extra_meta["P_plateau_reference"] = Pt_sol
        extra_meta["pressure_span_bridge"] = float(np.max(P_br) - np.min(P_br))
        extra_meta["pressure_span_fraction_of_endpoint_delta"] = float((np.max(P_br) - np.min(P_br)) / max(abs(P2 - P1), 1e-30))
        extra_meta["spline_mode"] = "smooth_near_pressure_plateau"

    elif method in ("cubic_P", "quintic_P"):
        mu_br = mu1 + u_grid * Dmu
        dn_du_1 = (n1 / (mu1 * cs2_1)) * Dmu
        dn_du_2 = (n2 / (mu2 * cs2_2)) * Dmu
        nb_br, dnb_du, int_nb, d0_p, d1_p, sp_mode = build_monotone_hermite_bridge(
            n1, n2, dn_du_1, dn_du_2, n_bar,
            match_slopes=(method == "quintic_P"),
            u_grid=u_grid,
            canonical_dx=None,
        )
        P_br = P1 + Dmu * int_nb
        Q1_br = P_br / nb_br
        F_br = mu_br - Q1_br
        cs2_br = (nb_br * Dmu) / (mu_br * dnb_du)
        cs2_left_br = float((n1 * Dmu) / (mu1 * d0_p))
        cs2_right_br = float((n2 * Dmu) / (mu2 * d1_p))
        extra_meta["spline_mode"] = sp_mode

    elif method == "css":
        def _css_knot(s: float) -> Tuple[float, float]:
            mu_m_v = mu1 * ((mu2 / mu1) ** (1.0 - s))
            nm_v = n1 + s * (n2 - n1)
            return mu_m_v, nm_v

        def _css_int(s: float) -> float:
            mu_m_v, nm_v = _css_knot(s)
            csA2 = np.log(mu_m_v / mu1) / np.log(nm_v / n1)
            csB2 = np.log(mu2 / mu_m_v) / np.log(n2 / nm_v)
            gA, gB = 1.0 / csA2, 1.0 / csB2
            IA = (n1 * mu1 / (gA + 1.0)) * ((mu_m_v / mu1) ** (gA + 1.0) - 1.0)
            IB = (nm_v * mu_m_v / (gB + 1.0)) * ((mu2 / mu_m_v) ** (gB + 1.0) - 1.0)
            return float((IA + IB) - (P2 - P1))

        s_sol = brentq(_css_int, 1e-6, 1.0 - 1e-6, xtol=1e-14)
        mu_m, nm_sol = _css_knot(s_sol)
        csA2 = float(np.log(mu_m / mu1) / np.log(nm_sol / n1))
        csB2 = float(np.log(mu2 / mu_m) / np.log(n2 / nm_sol))
        gA, gB = 1.0 / csA2, 1.0 / csB2
        P_m = P1 + (n1 * mu1 / (gA + 1.0)) * ((mu_m / mu1) ** (gA + 1.0) - 1.0)
        mu_br = mu1 * ((mu2 / mu1) ** u_grid)
        nb_br = np.empty_like(mu_br)
        P_br = np.empty_like(mu_br)
        cs2_br = np.empty_like(mu_br)
        for k, m_v in enumerate(mu_br):
            if m_v <= mu_m:
                nb_br[k] = n1 * ((m_v / mu1) ** gA)
                P_br[k] = P1 + (n1 * mu1 / (gA + 1.0)) * ((m_v / mu1) ** (gA + 1.0) - 1.0)
                cs2_br[k] = csA2
            else:
                nb_br[k] = nm_sol * ((m_v / mu_m) ** gB)
                P_br[k] = P_m + (nm_sol * mu_m / (gB + 1.0)) * ((m_v / mu_m) ** (gB + 1.0) - 1.0)
                cs2_br[k] = csB2
        Q1_br = P_br / nb_br
        F_br = mu_br - Q1_br
        cs2_left_br = csA2
        cs2_right_br = csB2
        extra_meta["csA2"] = csA2
        extra_meta["csB2"] = csB2
        extra_meta["spline_mode"] = "2segment_css"

    elif method in ("cubic_F", "quintic_F"):
        x_br = x1 + u_grid * Dx
        nb_br = np.exp(x_br)
        dQ1_du_1 = (mu1 * cs2_1 - Q1_1) * Dx
        dQ1_du_2 = (mu2 * cs2_2 - Q1_2) * Dx
        Q1_br, dQ1_du, int_Q1, d0_f, d1_f, sp_mode = build_monotone_hermite_bridge(
            Q1_1, Q1_2, dQ1_du_1, dQ1_du_2, q1_bar,
            match_slopes=(method == "quintic_F"),
            u_grid=u_grid,
            canonical_dx=Dx,
        )
        F_br = F1 + Dx * int_Q1
        P_br = nb_br * Q1_br
        mu_br = F_br + Q1_br
        cs2_br = (Q1_br + dQ1_du / Dx) / mu_br
        cs2_left_br = float((Q1_1 + d0_f / Dx) / mu1)
        cs2_right_br = float((Q1_2 + d1_f / Dx) / mu2)
        extra_meta["spline_mode"] = sp_mode

    nb_full = np.concatenate([nb_cr, nb_br, nb_co])
    Q1_full = np.concatenate([Q1_cr, Q1_br, Q1_co])
    F_full = np.concatenate([F_cr, F_br, F_co])
    mub_full = F_full + Q1_full
    N_full = len(nb_full)

    # C^2 quintic smootherstep sigma(u) = 6 u^5 - 15 u^4 + 10 u^3 (matches manual Section 2.2.2)
    s_u = (u_grid ** 3) * (10.0 - 15.0 * u_grid + 6.0 * (u_grid ** 2))

    core_composition_usable = False
    core_composition_rejection = "missing"
    if d_co.composition_available:
        core_comp_report = check_consistency(d_co, verbose=False)
        core_composition_usable = bool(
            core_comp_report["C3"]["status"] == "CHECKED"
            and core_comp_report["C4"]["status"] == "CHECKED"
            and core_comp_report["C3"]["pass_"]
            and core_comp_report["C4"]["pass_"]
        )
        if core_composition_usable:
            core_composition_rejection = "validated_baryon_and_charge_closure"
        else:
            core_composition_rejection = (
                "rejected_input_composition:"
                f"C3={core_comp_report['C3']['max_abs']:.3e},"
                f"C4={core_comp_report['C4']['max_abs']:.3e}"
            )
    else:
        core_composition_rejection = "input_eos_compo_absent"

    if core_composition_usable:
        compo_out_available = bool(d_cr.composition_available)
    elif core_compo_mode == "literature_beta_eq":
        if not d_cr.composition_available:
            raise ValueError("literature_beta_eq reconstruction also requires validated crust composition")
        # The approximation is registered only for APR and BBB. Unsupported
        # cores raise in _literature_core_symmetry_energy instead of guessing.
        compo_out_available = True
    else:
        compo_out_available = False

    if compo_out_available:
        Yh_cr, Yp_cr, Yn_cr = d_cr.Yh[:i1 + 1], d_cr.Yp[:i1 + 1], d_cr.Yn[:i1 + 1]
        Ye_cr, Ymu_cr = d_cr.Ye[:i1 + 1], d_cr.Ymu[:i1 + 1]
        A_cr, Z_cr, Iph_cr = d_cr.Aav[:i1 + 1], d_cr.Zav[:i1 + 1], d_cr.Iphase[:i1 + 1]

        Yh_co = np.zeros(len(nb_co))
        A_co = np.zeros(len(nb_co))
        Z_co = np.zeros(len(nb_co))
        if core_composition_usable:
            Iph_co = d_co.Iphase[i2:].copy()
            Yp_co = d_co.Yp[i2:].copy()
            Ye_co = d_co.Ye[i2:].copy()
            Ymu_co = d_co.Ymu[i2:].copy()
            Yn_co = 1.0 - Yp_co
            muq_co = d_co.muq[i2:].copy()
            compo_prov = "tabulated_crust_and_core"
        else:
            Iph_co = np.full(len(nb_co), int(Iph_cr[-1]) + 1 if int(Iph_cr[-1]) != 6 else 6, dtype=int)
            Yn_co, Yp_co, Ye_co, Ymu_co, muq_co = _solve_beta_eq_core_fractions(co_name, nb_co)
            compo_prov = "tabulated_crust + approximate_parabolic_symmetry_energy_beta_eq_core"

        Yh_br = Yh_cr[-1] * ((1.0 - s_u) ** 2)
        Z_br = Z_cr[-1] * (1.0 - 0.2 * s_u)
        za1 = Z_cr[-1] / A_cr[-1] if A_cr[-1] > 0 else 0.30
        A_br = Z_br / za1
        Yq1 = Yp_cr[-1] + za1 * Yh_cr[-1]
        Yq2 = Yp_co[0]
        Yq_br = (1.0 - s_u) * Yq1 + s_u * Yq2
        Yp_br = np.clip(Yq_br - za1 * Yh_br, 0.0, 1.0)
        Yh_br = (Yq_br - Yp_br) / za1
        Yn_br = 1.0 - Yp_br - Yh_br
        Ye_br = np.zeros(n_bridge)
        Ymu_br = np.zeros(n_bridge)
        for k in range(n_bridge):
            nq = nb_br[k] * Yq_br[k]
            kfe0 = HBARC * (3.0 * np.pi ** 2 * nq) ** (1.0 / 3.0)
            if np.sqrt(ME ** 2 + kfe0 ** 2) <= MMU:
                Ye_br[k] = Yq_br[k]
            else:
                ymu_v = float(brentq(
                    lambda ym: np.sqrt(ME ** 2 + (HBARC * (3.0 * np.pi ** 2 * nb_br[k] * max(0.0, Yq_br[k] - ym)) ** (1.0 / 3.0)) ** 2)
                    - np.sqrt(MMU ** 2 + (HBARC * (3.0 * np.pi ** 2 * nb_br[k] * ym) ** (1.0 / 3.0)) ** 2),
                    0.0, Yq_br[k], xtol=1e-14
                ))
                Ymu_br[k] = ymu_v
                Ye_br[k] = Yq_br[k] - ymu_v
        # Inherit the crust boundary phase ID (Finding F-06: never hardcode 7)
        Iph_br = np.full(n_bridge, int(Iph_cr[-1]), dtype=int)

        Yn_full = np.concatenate([Yn_cr, Yn_br, Yn_co])
        Yp_full = np.concatenate([Yp_cr, Yp_br, Yp_co])
        Ye_full = np.concatenate([Ye_cr, Ye_br, Ye_co])
        Ymu_full = np.concatenate([Ymu_cr, Ymu_br, Ymu_co])
        A_full = np.concatenate([A_cr, A_br, A_co])
        Z_full = np.concatenate([Z_cr, Z_br, Z_co])
        Yh_full = np.concatenate([Yh_cr, Yh_br, Yh_co])
        Iph_full = np.concatenate([Iph_cr, Iph_br, Iph_co]).astype(int)

        # Preserve exact provider-tabulated mu_q in untouched crust and core
        muq_cr = d_cr.muq[:i1 + 1]
        kFe_br = HBARC * (3.0 * np.pi ** 2 * np.maximum(0.0, nb_br * Ye_br)) ** (1.0 / 3.0)
        muq_br_fermi = -np.sqrt(ME ** 2 + kFe_br ** 2)
        # Smoothly blend endpoint lattice offsets across the bridge
        off1 = float(muq_cr[-1]) - (-float(np.sqrt(ME ** 2 + (HBARC * (3.0 * np.pi ** 2 * max(0.0, n1 * Ye_cr[-1])) ** (1.0 / 3.0)) ** 2)))
        off2 = float(muq_co[0]) - (-float(np.sqrt(ME ** 2 + (HBARC * (3.0 * np.pi ** 2 * max(0.0, n2 * Ye_co[0])) ** (1.0 / 3.0)) ** 2)))
        muq_br = muq_br_fermi + (1.0 - s_u) * off1 + s_u * off2
        muq_full = np.concatenate([muq_cr, muq_br, muq_co])
    else:
        Yn_full = Yp_full = Ye_full = Ymu_full = A_full = Z_full = Yh_full = np.zeros(N_full)
        Iph_full = np.zeros(N_full, dtype=int)
        muq_cr = d_cr.muq[:i1 + 1]
        muq_co = d_co.muq[i2:]
        muq_br = (1.0 - s_u) * float(muq_cr[-1]) + s_u * float(muq_co[0])
        muq_full = np.concatenate([muq_cr, muq_br, muq_co])
        compo_prov = "none_strict_thermo_only"

    Q2_full = np.zeros(N_full)
    Q3_full = mub_full / mn - 1.0
    Q4_full = muq_full / mn
    Q5_full = np.zeros(N_full)
    Q6_full = F_full / mn - 1.0
    Q7_full = Q6_full.copy()
    q_extra = np.zeros(N_full)

    matched = EOSData(
        nb=nb_full, mn=mn, mp=mp,
        Q1=Q1_full, Q2=Q2_full, Q3=Q3_full, Q4=Q4_full, Q5=Q5_full, Q6=Q6_full, Q7=Q7_full,
        q_extra=q_extra, Iphase=Iph_full,
        Yn=Yn_full, Yp=Yp_full, Ye=Ye_full, Ymu=Ymu_full,
        Aav=A_full, Zav=Z_full, Yh=Yh_full,
        Il=int(d_co.Il), iT=int(d_co.iT), iYq=int(d_co.iYq), inb_start=1,
        temperature_mev=0.0, charge_fraction_grid=0.0, q_extra_count=0,
        composition_available=compo_out_available,
        beta_equilibrated=True,
    )

    eps1, eps2 = n1 * F1, n2 * F2
    # The bridges are smooth single-valued barotropes with no density jump.
    # Seidov's criterion is therefore not applicable to these constructions.
    seidov_ratio = None
    regularity_class = {
        "maxwell": "smooth_near_pressure_plateau_no_phase_jump",
        "css": "piecewise_constant_sound_speed_continuous_barotrope",
        "cubic_F": "cubic_free_energy_bridge",
        "cubic_P": "cubic_pressure_chemical_potential_bridge",
        "quintic_F": "quintic_C2_free_energy_bridge",
        "quintic_P": "quintic_C2_pressure_chemical_potential_bridge",
    }[method]

    info = dict(
        crust=cr_name, core=co_name, method=method,
        i1=i1, i2=i2, n_bridge=n_bridge,
        n1=n1, n2=n2, dnb=n2 - n1,
        P1=P1, P2=P2, dP=P2 - P1,
        mu1=mu1, mu2=mu2, dmu=Dmu,
        eps1=eps1, eps2=eps2, deps=eps2 - eps1,
        seidov_ratio=seidov_ratio,
        seidov_applicable=False,
        is_exact_maxwell=False,
        phase_transition_order=None,
        phase_transition_class="not_assigned_no_density_jump",
        density_jump_fm3=0.0,
        construction_kind=regularity_class,
        regularity_class=regularity_class,
        cs2_1=cs2_1, cs2_2=cs2_2,
        cs2_left_br=cs2_left_br, cs2_right_br=cs2_right_br,
        dcs2_left=abs(cs2_left_br - cs2_1),
        dcs2_right=abs(cs2_right_br - cs2_2),
        cs2_min_bridge=float(np.min(cs2_br)),
        cs2_max_bridge=float(np.max(cs2_br)),
        pt_order=None,
        pt_class="not_assigned_no_density_jump",
        dm_crust=dm_cr, dm_core=dm_co,
        composition_provenance=compo_prov,
        core_composition_input_status=core_composition_rejection,
        core_composition_input_usable=bool(core_composition_usable),
        **extra_meta,
    )
    return matched, info


def solve_star_3stage(
    eos_data: EOSData,
    P1_MeV: float,
    P2_MeV: float,
    Pc_MeV: float,
    *,
    method_ivp: str = "DOP853",
    rtol: float = 1e-7,
) -> Tuple[float, float, float, float, float]:
    """Integrate TOV plus l=2 tidal perturbations for a smooth barotropic EOS.

    This implementation requires a strictly monotone tabulated P(epsilon) and
    does not apply density-jump tidal matching for first-order phase transitions.
    Central pressures are never clipped or extrapolated. Each integration stage
    must terminate at its requested pressure event or the calculation fails.
    """
    p_grid = np.asarray(eos_data.p, dtype=float)
    eps_grid = np.asarray(eos_data.eps, dtype=float)
    if p_grid.ndim != 1 or eps_grid.shape != p_grid.shape or p_grid.size < 4:
        raise ValueError("TOV EOS must contain matching one-dimensional pressure/energy arrays")
    if (np.any(~np.isfinite(p_grid)) or np.any(~np.isfinite(eps_grid))
            or np.any(p_grid <= 0) or np.any(eps_grid <= 0)
            or np.any(np.diff(p_grid) <= 0) or np.any(np.diff(eps_grid) <= 0)):
        raise ValueError("TOV solver requires finite, positive, strictly increasing P and epsilon")
    if not np.all(np.isfinite([P1_MeV, P2_MeV, Pc_MeV, rtol])) or rtol <= 0:
        raise ValueError("TOV pressures and rtol must be finite; rtol must be positive")
    P_surf_MeV = max(float(p_grid[0]), 1e-14)
    if not (P_surf_MeV < P1_MeV < P2_MeV < Pc_MeV <= float(p_grid[-1])):
        raise ValueError(
            "TOV pressure ordering/range invalid: require p_surface < P1 < P2 < Pc <= p_table_max; "
            f"got {P_surf_MeV:g}, {P1_MeV:g}, {P2_MeV:g}, {Pc_MeV:g}, {p_grid[-1]:g}"
        )

    lnP = np.log(p_grid)
    lneps = np.log(eps_grid)
    lneps_sp = PchipInterpolator(lnP, lneps, extrapolate=False)
    dlneps_sp = lneps_sp.derivative()
    lnP_min, lnP_max = float(lnP[0]), float(lnP[-1])
    lnPc = float(np.log(Pc_MeV))
    eps_c = float(np.exp(lneps_sp(lnPc)))
    dedp_c = float(dlneps_sp(lnPc) * (eps_c / Pc_MeV))
    if not np.isfinite(dedp_c) or dedp_c <= 0:
        raise ValueError("central d epsilon / dP is non-positive or non-finite")

    r0 = 1e-3
    P0 = Pc_MeV - (2.0 * np.pi / 3.0) * R_SCHW_SUN * CONV_MEV_FM3 * (eps_c + Pc_MeV) * (eps_c + 3.0 * Pc_MeV) * (r0 ** 2)
    m0 = (FOUR_PI / 3.0) * CONV_MEV_FM3 * eps_c * (r0 ** 3)
    y0 = 2.0 - (FOUR_PI / 7.0) * R_SCHW_SUN * CONV_MEV_FM3 * (5.0 * eps_c + 9.0 * Pc_MeV + (eps_c + Pc_MeV) * dedp_c) * (r0 ** 2)
    if not np.isfinite(P0 + m0 + y0) or P0 <= P2_MeV:
        raise ValueError("central expansion starts below the requested core/bridge boundary")

    def make_rhs(p_stop: float):
        def rhs(r: float, state: list[float]) -> list[float]:
            P_v, m_v, y_v = state
            P_eval = max(float(P_v), p_stop)
            lnp = float(np.clip(np.log(P_eval), lnP_min, lnP_max))
            eps_v = float(np.exp(lneps_sp(lnp)))
            dedp_v = float(dlneps_sp(lnp) * (eps_v / P_eval))
            if not np.isfinite(eps_v + dedp_v) or dedp_v <= 0:
                raise FloatingPointError("non-positive/non-finite EOS derivative during TOV integration")
            m_km = m_v * R_SCHW_SUN
            P_km = P_eval * CONV_MEV_FM3 * R_SCHW_SUN
            eps_km = eps_v * CONV_MEV_FM3 * R_SCHW_SUN
            elam = 1.0 / max(1e-8, 1.0 - 2.0 * m_km / r)
            dnu_dr = 2.0 * elam * (m_km + FOUR_PI * (r ** 3) * P_km) / (r ** 2)
            dP_dr = -0.5 * (eps_v + P_eval) * dnu_dr
            dm_dr = FOUR_PI * (r ** 2) * (eps_v * CONV_MEV_FM3)
            F_r = elam * (1.0 - FOUR_PI * (r ** 2) * (eps_km - P_km))
            Q_r = FOUR_PI * elam * (5.0 * eps_km + 9.0 * P_km + (eps_km + P_km) * dedp_v) - 6.0 * elam / (r ** 2) - (dnu_dr ** 2)
            dy_dr = -(y_v ** 2 + y_v * F_r + (r ** 2) * Q_r) / r
            return [dP_dr, dm_dr, dy_dr]
        return rhs

    def run_stage(r_start: float, state_start: list[float], p_target: float, max_step: float):
        if state_start[0] <= p_target:
            raise RuntimeError(f"TOV stage starts at P={state_start[0]:.8g} <= target {p_target:.8g}")
        ev = lambda r, state: state[0] - p_target
        ev.terminal = True
        ev.direction = -1
        sol = solve_ivp(make_rhs(p_target), [r_start, 35.0], state_start,
                        method=method_ivp, events=ev, rtol=rtol,
                        atol=1e-10, max_step=max_step)
        if not sol.success:
            raise RuntimeError(f"TOV integration failed before P={p_target:.8g}: {sol.message}")
        if not sol.t_events or len(sol.t_events[0]) == 0:
            raise RuntimeError(f"TOV integration did not reach pressure event P={p_target:.8g} before 35 km")
        state = sol.y_events[0][0]
        return float(sol.t_events[0][0]), [float(state[0]), float(state[1]), float(state[2])]

    r1, s1 = run_stage(r0, [P0, m0, y0], P2_MeV, max_step=0.15)
    r2, s2 = run_stage(r1, s1, P1_MeV, max_step=0.002)
    R_km, s3 = run_stage(r2, s2, P_surf_MeV, max_step=0.05)
    M_sun, y_R = s3[1], s3[2]
    C = M_sun * R_SCHW_SUN / R_km
    if not (0.01 < C < 0.48):
        raise RuntimeError(f"unphysical compactness C={C:.8g} after TOV integration")
    t1 = (8.0 / 5.0) * (C ** 5) * ((1.0 - 2.0 * C) ** 2) * (2.0 + 2.0 * C * (y_R - 1.0) - y_R)
    t2 = 2.0 * C * (6.0 - 3.0 * y_R + 3.0 * C * (5.0 * y_R - 8.0))
    t3 = 4.0 * (C ** 3) * (13.0 - 11.0 * y_R + C * (3.0 * y_R - 2.0) + 2.0 * (C ** 2) * (1.0 + y_R))
    t4 = 3.0 * ((1.0 - 2.0 * C) ** 2) * (2.0 - y_R + 2.0 * C * (y_R - 1.0)) * np.log(1.0 - 2.0 * C)
    denom = t2 + t3 + t4
    if not np.isfinite(denom) or abs(denom) < 1e-30:
        raise RuntimeError("tidal Love-number denominator is singular")
    k2 = t1 / denom
    Lam = (2.0 / 3.0) * k2 / (C ** 5)
    return float(M_sun), float(R_km), float(k2), float(Lam), float(R_km - r1)


def solve_star_at_target_mass(
    eos_data: EOSData,
    P1_MeV: float,
    P2_MeV: float,
    Pc_guess_MeV: float,
    M_target_Msun: float = 1.4,
    *,
    dM_dlnPc_hint: float | None = None,
    mass_tolerance: float = 1e-6,
) -> Tuple[float, float, float, float, float, float]:
    """Find a stable-branch central pressure with an actual TOV root solve.

    Returns the observables evaluated at the brentq root, never interpolated from
    neighbouring stars. If the requested mass cannot be bracketed on the
    monotonically increasing branch inside the tabulated EOS range, raises.
    """
    pmin = max(float(eos_data.p[0]), float(P2_MeV) * (1.0 + 1e-8))
    pmax = float(eos_data.p[-1])
    if not np.all(np.isfinite([P1_MeV, P2_MeV, Pc_guess_MeV, M_target_Msun, mass_tolerance])):
        raise ValueError("target-mass inputs must be finite")
    if M_target_Msun <= 0 or mass_tolerance <= 0 or not (P1_MeV < P2_MeV < pmax):
        raise ValueError("invalid target mass, tolerance, or matching pressure bounds")
    if not (pmin < Pc_guess_MeV <= pmax):
        raise ValueError(f"Pc guess {Pc_guess_MeV:g} is outside stable tabulated bounds ({pmin:g},{pmax:g}]")

    ln_guess = float(np.log(Pc_guess_MeV))
    M0, R0, k20, Lam0, dR0 = solve_star_3stage(eos_data, P1_MeV, P2_MeV, Pc_guess_MeV)
    f0 = M0 - M_target_Msun
    if abs(f0) <= mass_tolerance:
        return M0, R0, k20, Lam0, dR0, Pc_guess_MeV

    direction = 1.0 if f0 < 0.0 else -1.0
    step = 0.025
    if dM_dlnPc_hint is not None and np.isfinite(dM_dlnPc_hint) and abs(dM_dlnPc_hint) > 1e-6:
        step = min(0.20, max(0.005, abs(f0 / dM_dlnPc_hint)))
    previous_ln, previous_mass = ln_guess, M0
    bracket = None
    ln_lo, ln_hi = float(np.log(pmin)), float(np.log(pmax))
    for _ in range(32):
        trial_ln = previous_ln + direction * step
        trial_ln = min(ln_hi, max(ln_lo, trial_ln))
        if abs(trial_ln - previous_ln) < 1e-12:
            break
        trial_pc = float(np.exp(trial_ln))
        trial = solve_star_3stage(eos_data, P1_MeV, P2_MeV, trial_pc)
        trial_mass = float(trial[0])
        f_trial = trial_mass - M_target_Msun
        f_previous = previous_mass - M_target_Msun
        if f_trial == 0.0 or f_trial * f_previous < 0.0:
            bracket = (min(previous_ln, trial_ln), max(previous_ln, trial_ln))
            break
        # Refuse to cross a mass maximum or reverse direction on the stable branch.
        if direction > 0.0 and trial_mass <= previous_mass:
            break
        if direction < 0.0 and trial_mass >= previous_mass:
            break
        previous_ln, previous_mass = trial_ln, trial_mass
        step *= 1.7
    if bracket is None:
        raise ValueError(
            f"target mass {M_target_Msun:.6f} Msun is not bracketed on the stable branch "
            f"within Pc in [{pmin:.8g},{pmax:.8g}] MeV; latest M={previous_mass:.6f} Msun"
        )

    root_ln = brentq(
        lambda ln_pc: solve_star_3stage(eos_data, P1_MeV, P2_MeV, float(np.exp(ln_pc)))[0] - M_target_Msun,
        bracket[0], bracket[1], xtol=5e-11, rtol=5e-11, maxiter=100
    )
    Pc_root = float(np.exp(root_ln))
    M, R, k2, Lam, dR = solve_star_3stage(eos_data, P1_MeV, P2_MeV, Pc_root)
    if abs(M - M_target_Msun) > mass_tolerance:
        raise RuntimeError(
            f"bracketed target-mass solve did not converge: |M-M_target|={abs(M-M_target_Msun):.3e} Msun"
        )
    return float(M), float(R), float(k2), float(Lam), float(dR), Pc_root


__all__ = [
    "METHODS",
    "CRUST_NAMES",
    "CORE_NAMES",
    "get_core_start_index",
    "get_gauge_shifts",
    "find_convex_transition_window",
    "build_monotone_hermite_bridge",
    "build_maxwell_plateau_bridge",
    "match_crust_core",
    "solve_star_3stage",
    "solve_star_at_target_mass",
]
