r"""Reconstruction of a cold EOS table by filling a density gap or master rebuild.

Upgraded in the adversarial audit to enforce BOTH the Euler relation (C2:
:math:`p = n_b \mu_b - f`) AND the differential First Law / Maxwell identities
(C7: :math:`dp/dn_b = n_b d\mu_b/dn_b`, C8: :math:`p = n_b^2 dF/dn_b`, C9:
:math:`d\varepsilon/dn_b = \mu_b`) to exact analytical precision across filled
density gaps and in low-density outer crusts.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator, CubicHermiteSpline
from scipy.integrate import quad, cumulative_trapezoid

from .eos_data import EOSData

_HBARC = 197.3269804   # MeV fm
_ME = 0.51099895       # electron mass, MeV
_MMU = 105.6583745     # muon mass, MeV
_ALPHA = 1.0 / 137.035999084


def find_density_gap(nb: np.ndarray, min_ratio: float = 5.0) -> int:
    """Return the index of the largest logarithmic step in ``nb``."""
    dln = np.diff(np.log(nb))
    igap = int(dln.argmax())
    if nb[igap + 1] / nb[igap] < min_ratio:
        raise ValueError(f"no density gap found (max ratio "
                         f"{nb[igap+1]/nb[igap]:.2f} < {min_ratio})")
    return igap


def fill_density_gap(data: EOSData, n_fill: int = 60,
                     fit_points: int = 12, *, gap_phase: int | None = None) -> EOSData:
    r"""Return a new :class:`EOSData` with the density gap filled using a true
    First-Law + Euler-exact master function :math:`F(\ln n_b)`.

    Unlike independent PCHIP interpolation of :math:`Q_1` and :math:`\mu_b`
    (which violates the First Law :math:`dF/d\ln n_b = Q_1` and the work
    integral :math:`\Delta F = \int Q_1 d\ln n_b` across the gap), this method
    constructs :math:`F(x)` (:math:`x = \ln n_b`) so that:

    1. :math:`F(x_0) = F_0,\; F(x_1) = F_1` and :math:`F'(x_0) = Q_1(x_0),\;
       F'(x_1) = Q_1(x_1)`;
    2. :math:`Q_1(x) = F'(x)` integrates to :math:`F(x_1) - F(x_0)` to machine
       precision (:math:`\sim 10^{-16}`);
    3. :math:`\mu_b(x) = F(x) + Q_1(x)` and :math:`p(x) = n_b Q_1(x)`, so C2,
       C7, C8, and C9 hold simultaneously.

    For crust gap fills where nuclear charge :math:`Z` is tabulated on both
    sides, :math:`Z` and :math:`Y_e` are interpolated with monotone PCHIP and
    cell mass number is derived as :math:`A = Z / (Y_e + Y_\mu)` (preventing
    unphysical :math:`Z` overshoot across outer-to-pasta transitions).
    """
    nb = data.nb
    N = data.n
    igap = find_density_gap(nb)
    if n_fill < 1:
        raise ValueError("n_fill must be >= 1")
    left_phase = int(data.Iphase[igap]) if data.composition_available else 0
    right_phase = int(data.Iphase[igap + 1]) if data.composition_available else 0
    if data.composition_available and left_phase != right_phase and gap_phase is None:
        raise ValueError(
            f"density gap crosses phase IDs {left_phase}->{right_phase}; "
            "pass gap_phase explicitly instead of silently assigning a phase"
        )
    if gap_phase is None:
        gap_phase = left_phase

    nb_fill = np.exp(np.linspace(np.log(nb[igap]), np.log(nb[igap + 1]),
                                 n_fill + 2)[1:-1])
    nb_new = np.concatenate([nb[:igap + 1], nb_fill, nb[igap + 1:]])
    Nnew = N + n_fill
    lo, hi = igap + 1, igap + 1 + n_fill

    fit_idx = np.concatenate([
        np.arange(max(0, igap - fit_points), igap + 1),
        np.arange(igap + 1, min(N, igap + 1 + fit_points)),
    ])
    xfig = np.log(nb[fit_idx])
    xfill = np.log(nb_fill)

    x0, x1 = float(np.log(nb[igap])), float(np.log(nb[igap + 1]))
    Dx = x1 - x0
    F0 = float(data.F[igap])
    F1 = float(data.F[igap + 1])
    Q1_0 = float(data.Q1[igap])
    Q1_1 = float(data.Q1[igap + 1])

    # Build PCHIP for Q1(x) and add the C^1 quartic bump phi(u) = 30 u^2 (1-u)^2
    # so that Integral_{x0}^{x1} Q1(x) dx == F1 - F0 identically.
    # If Q1_fill would become non-positive (as in giant 191x gaps like GM1 where
    # Q1_0 ~ 1.88 and Q1_1 ~ 1.93 while <Q1> = 0.795), use CubicHermiteSpline
    # directly on ([x0, x1], [F0, F1], [Q1_0, Q1_1]).
    q1_fit = np.asarray(data.Q1[fit_idx], dtype=float)
    if np.any(~np.isfinite(q1_fit)) or np.any(q1_fit <= 0):
        raise ValueError("log-Q1 gap interpolation requires finite positive Q1 anchor values")
    logQ_pc = PchipInterpolator(xfig, np.log(q1_fit))
    def Q1_pc(x):
        return np.exp(logQ_pc(x))
    def I_pc_to(x):
        return float(quad(lambda xx: float(Q1_pc(xx)), x0, float(x),
                          epsabs=1e-11, epsrel=1e-10, limit=200)[0])
    I_pc = I_pc_to(x1)
    Delta_Q1 = (F1 - F0 - I_pc) / Dx
    u = (xfill - x0) / Dx
    phi = 30.0 * u ** 2 * (1.0 - u) ** 2
    int_phi = u ** 3 * (10.0 - 15.0 * u + 6.0 * u ** 2)

    Q1_cand = Q1_pc(xfill) + Delta_Q1 * phi
    F_cand = np.array([F0 + I_pc_to(xf) for xf in xfill]) + (F1 - F0 - I_pc) * int_phi
    P_cand = nb_fill * Q1_cand

    def _candidate_is_physical(x_dense, F_dense, Q_dense):
        # At T=0, epsilon=n F and mu_b=F+Q1.  With x=ln(n),
        # dP/dn=Q1+dQ1/dx and d epsilon/dn=F+Q1.
        dQdx = np.gradient(Q_dense, x_dense, edge_order=2)
        mu_dense = F_dense + Q_dense
        dpdn = Q_dense + dQdx
        cs2_dense = np.divide(dpdn, mu_dense,
                              out=np.full_like(dpdn, np.nan), where=mu_dense > 0)
        return bool(
            np.all(np.isfinite(F_dense)) and np.all(np.isfinite(Q_dense))
            and np.all(Q_dense > 0) and np.all(mu_dense > 0)
            and np.all(dpdn > 0) and np.all(cs2_dense >= -1e-8)
            and np.all(cs2_dense < 1.0 + 1e-6)
        )

    x_dense = np.linspace(x0, x1, max(512, 8 * (n_fill + 2)))
    u_dense = (x_dense - x0) / Dx
    int_phi_dense = u_dense ** 3 * (10 - 15 * u_dense + 6 * u_dense ** 2)
    base_dense = Q1_pc(x_dense)
    F_dense = F0 + cumulative_trapezoid(base_dense, x_dense, initial=0.0) + (F1 - F0 - I_pc) * int_phi_dense
    Q_dense = base_dense + Delta_Q1 * (30 * u_dense ** 2 * (1 - u_dense) ** 2)
    if _candidate_is_physical(x_dense, F_dense, Q_dense):
        Q1_fill = Q1_cand
        F_fill = F_cand
    else:
        chs = CubicHermiteSpline([x0, x1], [F0, F1], [Q1_0, Q1_1])
        F_dense_fallback = chs(x_dense)
        Q_dense_fallback = chs.derivative()(x_dense)
        if not _candidate_is_physical(x_dense, F_dense_fallback, Q_dense_fallback):
            raise ValueError(
                "no positive, mechanically stable and causal First-Law gap interpolant "
                "satisfies the endpoint constraints; refusing to emit an unchecked EOS"
            )
        F_fill = chs(xfill)
        Q1_fill = chs.derivative()(xfill)

    mub_fill = F_fill + Q1_fill

    def _splice(orig):
        out = np.empty(Nnew)
        out[:igap + 1] = orig[:igap + 1]
        out[hi:] = orig[igap + 1:]
        return out

    Q1n = _splice(data.Q1); Q1n[lo:hi] = Q1_fill
    Q3n = _splice(data.Q3); Q3n[lo:hi] = mub_fill / data.mn - 1.0
    Q6n = _splice(data.Q6); Q6n[lo:hi] = F_fill / data.mn - 1.0
    Q4n = _splice(data.Q4)
    Q4n[lo:hi] = PchipInterpolator(xfig, data.Q4[fit_idx])(xfill)
    Q2n = np.zeros(Nnew)
    Q5n = np.zeros(Nnew)
    Q7n = Q6n.copy()
    # Set q_extra in the gap to match the post-gap phase convention
    qn = _splice(data.q_extra)
    qn[lo:hi] = 1.0 if (data.q_extra.size and np.any(data.q_extra != 0)) else 0.0

    # For GPPVA tables where the original post-gap point has an ~8.8e-5 offset
    # between F_orig and mub_orig - Q1_orig, enforce F = mub - Q1 on the
    # original points only if needed for C2/C8 continuity, or keep original
    # points untouched as designed.

    def _pchip(arr):
        pc = PchipInterpolator(xfig, arr[fit_idx])
        out = _splice(arr)
        out[lo:hi] = pc(xfill)
        return out

    Ye_n = _pchip(data.Ye)
    Zav_n = _pchip(data.Zav)
    Yh_n = _pchip(data.Yh)
    Ymu_n = _pchip(data.Ymu)
    Yp_n = _pchip(data.Yp)
    Iphase_n = np.concatenate([data.Iphase[:igap + 1],
                               np.full(n_fill, int(gap_phase)),
                               data.Iphase[igap + 1:]]).astype(int)

    # Derive A from interpolated Z and (Ye + Ymu) where nuclei exist so Z stays
    # strictly monotone between its outer-crust and pasta endpoints (no Z overshoot)
    Aav_n = _splice(data.Aav)
    lept_gap = Ye_n[lo:hi] + Ymu_n[lo:hi]
    if np.all(lept_gap > 1e-8) and np.all(Zav_n[lo:hi] > 0):
        Aav_n[lo:hi] = Zav_n[lo:hi] / lept_gap
        # Keep Yh_n smoothly PCHIP-interpolated from 1.0 (outer crust) to pasta Yh
    else:
        Aav_n = _pchip(data.Aav)
        Zav_n[lo:hi] = lept_gap * Aav_n[lo:hi]

    Yn_n = np.zeros(Nnew)
    Yn_n[:igap + 1] = data.Yn[:igap + 1]
    Yn_n[hi:] = data.Yn[igap + 1:]

    return EOSData(nb=nb_new, mn=data.mn, mp=data.mp,
                   Q1=Q1n, Q2=Q2n, Q3=Q3n, Q4=Q4n, Q5=Q5n, Q6=Q6n, Q7=Q7n,
                   q_extra=qn, Iphase=Iphase_n,
                   Yn=Yn_n, Yp=Yp_n, Ye=Ye_n, Ymu=Ymu_n,
                   Aav=Aav_n, Zav=Zav_n, Yh=Yh_n,
                   Il=data.Il, iT=data.iT, iYq=data.iYq, inb_start=data.inb_start,
                   temperature_mev=data.temperature_mev,
                   charge_fraction_grid=data.charge_fraction_grid,
                   q_extra_count=data.q_extra_count,
                   composition_available=data.composition_available,
                   beta_equilibrated=data.beta_equilibrated)


def reconstruct_master(data: EOSData, smooth: bool = True,
                       smooth_sigma: float = 0.05,
                       fix_dissolution: bool = True,
                       dissolution_tol: float = 0.015) -> EOSData:
    r"""Rebuild a table from the free energy per baryon as master function.

    Upgraded to:
      1. Use ``dissolution_tol = 0.015`` (or detect an anomalous drop in :math:`A`
         or :math:`Y_e` at the last crust point) so that the anomalous remnant
         nucleus at ``idx = 142`` of ``RG_SLY4`` is properly dissolved into
         homogeneous matter instead of corrupting :math:`Y_e[142]` and :math:`Q_4[142]`.
      2. Clip ``Yh_f`` to :math:`[0, 1]` before computing ``Yn_f`` so baryon
         conservation holds to machine precision even when raw outer-crust rows
         have :math:`Y_h = 1.0000007`.
      3. Enforce exact relativistic electron + Wigner-Seitz lattice pressure in
         the low-density outer crust (:math:`Y_n = 0,\; A > 0`) so that spline
         boundary run-off never distorts low-density crust pressures or makes
         :math:`Q_1(n_b)` non-monotonic below neutron drip.
      4. Enforce exact muon beta-equilibrium (:math:`\mu_\mu = \mu_e`) where
         :math:`\mu_e \ge m_\mu`.
      5. Recompute ``q_extra`` as :math:`c_s^2 = dP/d\varepsilon` when the
         table stores sound speed squared in its extra thermo column.
    """
    from scipy.interpolate import UnivariateSpline

    nb = data.nb
    N = data.n
    mn = data.mn
    lognb = np.log(nb)

    Yp_f = np.clip(data.Yp.copy().astype(float), 0.0, 1.0)
    Ymu_f = np.clip(data.Ymu.copy().astype(float), 0.0, 1.0)
    Yh_f = np.clip(data.Yh.copy().astype(float), 0.0, 1.0)
    Aav_f = data.Aav.copy().astype(float)
    Zav_f = data.Zav.copy().astype(float)
    Iphase_f = data.Iphase.copy()

    crust = data.has_nucleus
    core = data.is_homogeneous

    # Dissolve an anomalous remnant nucleus at the crust-core boundary
    if fix_dissolution and crust.any() and core.any():
        last_crust = int(np.where(crust)[0][-1])
        ZA = np.where(Aav_f > 0, Zav_f / Aav_f, 0.0)
        resid = abs(Yp_f[last_crust] + ZA[last_crust] * Yh_f[last_crust]
                    - (data.Ye[last_crust] + Ymu_f[last_crust]))
        # Also detect if Ye = ZA*Yh at last_crust drops anomalously (>50% drop
        # relative to both neighbours, as in already-processed RG_SLY4_fixed)
        lo = last_crust - 1
        hi = last_crust + 1
        anom_drop = False
        if lo >= 0 and hi < N and core[hi]:
            ye_cand = Yp_f[last_crust] + ZA[last_crust] * Yh_f[last_crust] - Ymu_f[last_crust]
            ye_prev = Yp_f[lo] + ZA[lo] * Yh_f[lo] - Ymu_f[lo]
            ye_next = Yp_f[hi] - Ymu_f[hi]
            if ye_cand < 0.5 * min(ye_prev, ye_next):
                anom_drop = True
        if (resid > dissolution_tol or anom_drop) and last_crust + 2 < N:
            if core[hi] and lo >= 0:
                Yq_lo = (Zav_f[lo] / Aav_f[lo]) * Yh_f[lo]
                Yq_hi = Yp_f[hi]
                t = (lognb[last_crust] - lognb[lo]) / (lognb[hi] - lognb[lo])
                Yq = Yq_lo * (1.0 - t) + Yq_hi * t
                Iphase_f[last_crust] = Iphase_f[hi]
                Yh_f[last_crust] = 0.0
                Aav_f[last_crust] = 0.0
                Zav_f[last_crust] = 0.0
                Yp_f[last_crust] = Yq
                Ymu_f[last_crust] = 0.0

    # Baryon conservation: Yn = 1 - Yp - Yh
    Yn_f = np.clip(1.0 - Yp_f - Yh_f, 0.0, 1.0)
    Yh_f = np.clip(1.0 - Yp_f - Yn_f, 0.0, 1.0)

    # Charge neutrality: Ye + Ymu = Yp + (Z/A) Yh
    ZA = np.where(Aav_f > 0, Zav_f / Aav_f, 0.0)
    Yq_tot = np.clip(Yp_f + ZA * Yh_f, 0.0, 1.0)

    # Solve coupled beta-equilibrium mu_e(n_e) = mu_mu(n_mu) with n_e + n_mu = nb * Yq_tot
    # where mu_e >= m_mu, so that Y_mu and Y_e are self-consistent to machine precision!
    Ye_f = np.zeros(N)
    Ymu_f_new = np.zeros(N)
    from scipy.optimize import brentq
    for i in range(N):
        nq = nb[i] * Yq_tot[i]
        if nq <= 0:
            continue
        kFe0 = (3.0 * np.pi ** 2 * nq) ** (1.0 / 3.0)
        mue0 = np.sqrt(_ME ** 2 + (_HBARC * kFe0) ** 2)
        if mue0 <= _MMU:
            Ye_f[i] = Yq_tot[i]
            Ymu_f_new[i] = 0.0
        else:
            def _diff_mu(ymu_val):
                ye_val = max(Yq_tot[i] - ymu_val, 0.0)
                kfe = (3.0 * np.pi ** 2 * (nb[i] * ye_val)) ** (1.0 / 3.0)
                kfmu = (3.0 * np.pi ** 2 * (nb[i] * ymu_val)) ** (1.0 / 3.0)
                return np.sqrt(_ME ** 2 + (_HBARC * kfe) ** 2) - np.sqrt(_MMU ** 2 + (_HBARC * kfmu) ** 2)
            ymu_eq = brentq(_diff_mu, 0.0, Yq_tot[i], xtol=1e-14)
            Ymu_f_new[i] = ymu_eq
            Ye_f[i] = Yq_tot[i] - ymu_eq
    Ymu_f = Ymu_f_new

    # Charge chemical potential mu_q = -mu_e (including Wigner-Seitz lattice shift in crust if desired,
    # or pure Fermi gas with exact muon consistency; we keep exact Fermi gas for C1-C11 compatibility)
    n_e = Ye_f * nb
    kFe = (3.0 * np.pi ** 2 * np.maximum(n_e, 0.0)) ** (1.0 / 3.0)
    mue = np.sqrt(_ME ** 2 + (_HBARC * kFe) ** 2)
    muq_f = -mue

    # Master function F(nb); derive P, mu_b, eps
    F = (data.Q6 + 1.0) * mn
    if smooth:
        s = N * smooth_sigma ** 2
        spl = UnivariateSpline(lognb, F, k=4, s=s)
        F_sm = spl(lognb)
        Q1_sm = spl.derivative()(lognb)
    else:
        F_sm = F.copy()
        Q1_sm = np.gradient(F, lognb)

    # Fix low-density outer-crust spline run-off (where Yn == 0 and A > 0):
    outer = np.where((Aav_f > 0) & (Yn_f == 0))[0]
    if len(outer) > 5:
        i_drip = int(outer[-1])
        i_join = min(N - 1, i_drip + max(8, int(0.025 * N)))
        x_e = _HBARC * kFe / _ME
        pref = _ME ** 4 / (24.0 * np.pi ** 2 * _HBARC ** 3)
        P_el = pref * (x_e * (2.0 * x_e ** 2 - 3.0) * np.sqrt(1.0 + x_e ** 2) + 3.0 * np.arcsinh(x_e))
        eps_lat = -0.9 * (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0) * _ALPHA * np.maximum(Zav_f, 0.0) ** (2.0 / 3.0) * _HBARC * kFe * n_e
        P_lat = (1.0 / 3.0) * eps_lat
        Q1_phys = (P_el + P_lat) / nb

        idx_bridge = np.concatenate([np.arange(0, i_drip + 1), np.arange(i_join, N)])
        Q1_bridge_vals = np.concatenate([Q1_phys[:i_drip + 1], Q1_sm[i_join:]])
        q1_pc = PchipInterpolator(lognb[idx_bridge], Q1_bridge_vals)
        for i in range(i_join):
            Q1_sm[i] = float(q1_pc(lognb[i]))
            F_sm[i] = float(F_sm[i_join] - q1_pc.integrate(lognb[i], lognb[i_join]))

    P_sm = nb * Q1_sm
    mub_sm = F_sm + Q1_sm
    eps_sm = nb * F_sm

    Q1n = Q1_sm
    Q2n = np.zeros(N)
    Q3n = mub_sm / mn - 1.0
    Q4n = muq_f / mn
    Q5n = np.zeros(N)
    Q6n = F_sm / mn - 1.0
    Q7n = Q6n.copy()

    # If the table carries an extra column (like cs^2 in RG_KDE0v1), recompute cs^2 = dP/deps
    deps_dnb = np.gradient(eps_sm, nb)
    dP_dnb = np.gradient(P_sm, nb)
    cs2_sm = np.where(np.abs(deps_dnb) > 0, dP_dnb / deps_dnb, 0.0)
    if data.q_extra.size == N and np.any(data.q_extra != 0):
        qn = cs2_sm.copy()
    else:
        qn = np.zeros(N)

    return EOSData(nb=nb, mn=mn, mp=data.mp,
                   Q1=Q1n, Q2=Q2n, Q3=Q3n, Q4=Q4n, Q5=Q5n, Q6=Q6n, Q7=Q7n,
                   q_extra=qn, Iphase=Iphase_f,
                   Yn=Yn_f, Yp=Yp_f, Ye=Ye_f, Ymu=Ymu_f,
                   Aav=Aav_f, Zav=Zav_f, Yh=Yh_f,
                   Il=data.Il, iT=data.iT, iYq=data.iYq, inb_start=data.inb_start,
                   temperature_mev=data.temperature_mev,
                   charge_fraction_grid=data.charge_fraction_grid,
                   q_extra_count=data.q_extra_count,
                   composition_available=data.composition_available,
                   beta_equilibrated=data.beta_equilibrated)


__all__ = ["fill_density_gap", "find_density_gap", "reconstruct_master"]
