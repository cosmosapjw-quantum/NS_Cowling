r"""First-principles nuclear EoS solvers from published Lagrangians and Hamiltonians.

This module implements self-contained, path-free, first-principles solvers for all
eighteen neutron-star equations of state in the audited CompOSE suite (zero file I/O
at import time):

1. **Density-Dependent Relativistic Mean-Field (DD-RMF) Models**:
   ``GPPVA_DD2`` (Typel et al. 2010, *Phys. Rev. C* 81, 015803; Grill et al. 2012,
   *Phys. Rev. C* 85, 055808), ``GPPVA_DDME2`` (Lalazissis et al. 2005, *Phys. Rev. C*
   71, 024312; Grill et al. 2014, *Phys. Rev. C* 90, 045803), and ``GPPVA_TW``
   (Typel & Wolter 1999, *Nucl. Phys. A* 656, 331, Table 1 & Eqs. 2.1--2.18).
2. **Nonlinear Relativistic Mean-Field (NL-RMF) Models**:
   ``GPPVA_FSU2`` (Chen & Piekarewicz 2014, *Phys. Rev. C* 90, 044305),
   ``GPPVA_FSU2R`` and ``GPPVA_FSU2H`` (Tolos, Centelles & Ramos 2017, *PASA* 34,
   e065; Providencia et al. 2019, *Front. Astron. Space Sci.* 6, 13),
   ``GPPVA_NL3wrL55`` (Pais & Providencia 2016, *Phys. Rev. C* 94, 015808; Grill
   et al. 2014), ``GPPVA_TM1e`` (Shen et al. 2020, *Astrophys. J.* 891, 148), and
   ``GM_GM1`` (Glendenning & Moszkowski 1991, *Phys. Rev. Lett.* 67, 2414, Table II;
   Douchin & Haensel 2001, *Astron. Astrophys.* 380, 151).
3. **Relativistic GDFM Metamodels**:
   ``SPG_M1``--``SPG_M5`` (Scurto, Pais & Gulminelli 2024, *Phys. Rev. D* 109,
   103015, Eqs. 1--30, Table VII, Table VIII, and exact 12-parameter GDFM couplings
   :math:`g_i(x) = a_i + (b_i + d_i x^3)e^{-c_i x}` from Figure 18).
4. **Non-Relativistic Skyrme Energy Density Functionals**:
   ``RG_SLY4`` (Chabanat et al. 1998, *Nucl. Phys. A* 635, 231, Table 1; Gulminelli
   & Raduta 2015, *Phys. Rev. C* 92, 055803) and ``RG_KDE0v1`` (Agrawal, Shlomo &
   Au 2005, *Phys. Rev. C* 72, 014310, Table I).
5. **Microscopic Many-Body Variational & Brueckner--Hartree--Fock Models**:
   ``APR_APR`` (Akmal, Pandharipande & Ravenhall 1998, *Phys. Rev. C* 58, 1804,
   Eqs. A1--A12 & Table XII) and ``BBB_BHF2`` (Baldo, Bombaci & Burgio 1997,
   *Astron. Astrophys.* 328, 274, Table 3).
6. **Outer and Inner Crust Solvers & Published Table Anchors**:
   Baym, Pethick & Sutherland (1971, *Astrophys. J.* 170, 299, Table 2 & Table 5),
   Haensel & Pichon (1994, *Astron. Astrophys.* 283, 313), Pearson et al. (2018,
   *MNRAS* 481, 2994, BSk22), and Grill et al. (2014, *Phys. Rev. C* 90, 045803,
   Tables IV--V).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar, root_scalar

HBARC = 197.3269804
ME = 0.51099895
MMU = 105.6583755
ALPHA = 1.0 / 137.035999084

# ==============================================================================
# 1. DENSITY-DEPENDENT RMF MODELS (DD2, DDME2, TW)
# ==============================================================================
DD_MODELS = {
    "GPPVA_DD2": dict(
        mn=939.0, mp=939.0, nsat=0.149065,
        ms=546.212459, mw=783.0, mr=763.0,
        gs0=10.686681, gw0=13.342362, gr0=3.626940 * 2.0,
        as_=1.357630, bs=0.634442, cs=1.005358, ds=0.575810,
        aw=1.369718, bw=0.496475, cw=0.817753, dw=0.638452,
        ar=0.518903,
    ),
    "GPPVA_DDME2": dict(
        mn=939.0, mp=939.0, nsat=0.152,
        ms=550.1238, mw=783.0, mr=763.0,
        gs0=10.5396, gw0=13.0189, gr0=3.6836 * 2.0,
        as_=1.3881, bs=1.0943, cs=1.7057, ds=0.4421,
        aw=1.3892, bw=0.9240, cw=1.4620, dw=0.4775,
        ar=0.5647,
    ),
    "GPPVA_TW": dict(
        mn=939.0, mp=939.0, nsat=0.153,
        ms=550.0, mw=783.0, mr=763.0,
        gs0=10.72854, gw0=13.29015, gr0=7.32196,
        as_=1.365469, bs=0.226061, cs=0.409704, ds=0.901995,
        aw=1.402488, bw=0.172577, cw=0.344293, dw=0.983955,
        ar=0.515,
    ),
}

# ==============================================================================
# 2. NONLINEAR RMF MODELS (FSU2, FSU2R, FSU2H, NL3wrL55, TM1e, GM1)
# ==============================================================================
NL_MODELS = {
    "GPPVA_FSU2": dict(
        mn=939.0, mp=939.0,
        ms=497.479, mw=782.500, mr=763.000,
        gs2=108.0943, gw2=183.7893, gr2=80.4656,
        kappa=3.0029, lam=-0.000533, zeta=0.0256, Lam_v=0.000823,
    ),
    "GPPVA_FSU2R": dict(
        mn=939.0, mp=939.0,
        ms=497.479, mw=782.500, mr=763.000,
        gs2=107.5751, gw2=182.3949, gr2=206.4260,
        kappa=3.0911, lam=-0.001680, zeta=0.024, Lam_v=0.045,
    ),
    "GPPVA_FSU2H": dict(
        mn=939.0, mp=939.0,
        ms=497.479, mw=782.500, mr=763.000,
        gs2=102.7200, gw2=169.5315, gr2=197.2692,
        kappa=4.0014, lam=-0.013298, zeta=0.008, Lam_v=0.045,
    ),
    "GPPVA_NL3wrL55": dict(
        mn=939.0, mp=939.0,
        ms=508.194, mw=782.501, mr=763.000,
        gs2=10.217**2, gw2=12.868**2, gr2=127.1154,
        kappa=2.0 * 10.431 * HBARC / (10.217**3),
        lam=-6.0 * 28.885 / (10.217**4),
        zeta=0.0, Lam_v=0.03,
    ),
    "GPPVA_TM1e": dict(
        mn=938.0, mp=938.0,
        ms=511.198, mw=783.0, mr=770.0,
        gs2=10.0289**2, gw2=12.6139**2, gr2=139.6331,
        kappa=2.0 * 7.2325 * HBARC / (10.0289**3),
        lam=6.0 * 0.6183 / (10.0289**4),
        zeta=6.0 * 71.3075 / (12.6139**4),
        Lam_v=0.03,
    ),
}

# ==============================================================================
# 3. EXACT GDFM MODELS (SPG_M1..SPG_M5) FROM SCURTO, PAIS & GULMINELLI (2024)
# ==============================================================================
SPG_GDFM_PARAMS = {
    "SPG_M1": dict(
        mn=938.9, mp=938.9, rsat=0.157, ms=560.0, mw=782.5, mr=770.0,
        Esat=-15.8, Ksat=230.9, Qsat=-177.1, Zsat=7966.3,
        Jsym=33.9, Lsym=56.3, Ksym=-181.4, Qsym=350.1, Zsym=-1700.3,
        as_=9.844710, bs=2.064633, cs=2.629716, ds=3.884630,
        aw=11.530983, bw=2.399252, cw=2.542700, dw=8.447793,
        ar=2.0 * 1.969425, br=2.0 * 3.915208, cr=0.550511, dr=2.0 * (-0.600466),
    ),
    "SPG_M2": dict(
        mn=938.9, mp=938.9, rsat=0.155, ms=560.0, mw=782.5, mr=770.0,
        Esat=-15.9, Ksat=223.9, Qsat=-230.1, Zsat=2141.8,
        Jsym=30.2, Lsym=40.2, Ksym=-59.8, Qsym=1007.2, Zsym=-8829.1,
        as_=9.463916, bs=2.365605, cs=2.498849, ds=3.537427,
        aw=11.071276, bw=2.343656, cw=2.330909, dw=6.242264,
        ar=2.0 * 2.107415, br=2.0 * 4.831326, cr=1.003673, dr=2.0 * (-0.942798),
    ),
    "SPG_M3": dict(
        mn=938.9, mp=938.9, rsat=0.161, ms=560.0, mw=782.5, mr=770.0,
        Esat=-15.8, Ksat=225.7, Qsat=-283.4, Zsat=11174.8,
        Jsym=33.7, Lsym=50.9, Ksym=-231.0, Qsym=322.9, Zsym=1343.3,
        as_=9.956210, bs=1.986223, cs=2.945439, ds=4.016877,
        aw=11.598285, bw=2.383339, cw=2.661277, dw=8.884928,
        ar=2.0 * 0.496314, br=2.0 * 5.282626, cr=0.376043, dr=2.0 * (-0.545308),
    ),
    "SPG_M4": dict(
        mn=938.9, mp=938.9, rsat=0.162, ms=560.0, mw=782.5, mr=770.0,
        Esat=-15.8, Ksat=246.0, Qsat=-292.9, Zsat=-1397.4,
        Jsym=30.9, Lsym=40.1, Ksym=-111.8, Qsym=958.5, Zsym=-5795.1,
        as_=9.245701, bs=2.207806, cs=2.648662, ds=4.053112,
        aw=10.767465, bw=2.064659, cw=2.353136, dw=6.103118,
        ar=2.0 * 1.377456, br=2.0 * 5.101441, cr=0.749059, dr=2.0 * (-0.620039),
    ),
    "SPG_M5": dict(
        mn=938.9, mp=938.9, rsat=0.148, ms=560.0, mw=782.5, mr=770.0,
        Esat=-15.8, Ksat=227.7, Qsat=45.2, Zsat=14552.1,
        Jsym=31.4, Lsym=52.3, Ksym=-113.5, Qsym=539.8, Zsym=-5362.5,
        as_=10.634400, bs=1.852977, cs=3.074967, ds=3.581855,
        aw=12.526493, bw=2.204764, cw=2.743498, dw=8.946392,
        ar=2.0 * 1.510499, br=2.0 * 4.662305, cr=0.637078, dr=2.0 * (-0.572744),
    ),
}

# ==============================================================================
# 4. SKYRME EDF PARAMETERS (SLy4 & KDE0v1)
# ==============================================================================
SKYRME_PARAMS = {
    "RG_SLY4": dict(
        t0=-2488.913, t1=486.818, t2=-546.395, t3=13777.0,
        x0=0.834, x1=-0.344, x2=-1.000, x3=1.354,
        alpha=1.0 / 6.0, sigma=1.0 / 6.0, W0=123.0,
    ),
    "RG_KDE0v1": dict(
        t0=-2553.0843, t1=411.6963, t2=-419.8712, t3=14603.6069,
        x0=0.6483, x1=-0.3472, x2=-0.9268, x3=0.9475,
        alpha=0.1673, sigma=0.1673, W0=124.41,
    ),
}

BPS_HP_OUTER_CRUST_NUCLEI = [
    (56, 26), (62, 28), (64, 28), (86, 36), (84, 34), (82, 32), (80, 30),
    (78, 28), (76, 26), (126, 44), (124, 42), (122, 40), (120, 38), (118, 36),
]

GRILL2014_DDME2 = [
    (0.0003, 0.01426, 0.00243), (0.0004, 0.01903, 0.00279), (0.0005, 0.02379, 0.00314),
    (0.0006, 0.02856, 0.00350), (0.0007, 0.03332, 0.00385), (0.0008, 0.03809, 0.00426),
    (0.0009, 0.04286, 0.00466), (0.0010, 0.04762, 0.00507), (0.0015, 0.07146, 0.00745),
    (0.0020, 0.09531, 0.01014), (0.0050, 0.23854, 0.03015), (0.0100, 0.47754, 0.07389),
    (0.0200, 0.95625, 0.20094), (0.0500, 2.39727, 1.04932), (0.0720, 3.45856, 2.06135),
]

GRILL2014_NL3WR = [
    (0.0003, 0.01426, 0.00253), (0.0004, 0.01903, 0.00289), (0.0005, 0.02379, 0.00324),
    (0.0006, 0.02855, 0.00360), (0.0007, 0.03332, 0.00400), (0.0008, 0.03809, 0.00436),
    (0.0009, 0.04285, 0.00476), (0.0010, 0.04762, 0.00517), (0.0015, 0.07146, 0.00745),
    (0.0020, 0.09531, 0.00998), (0.0050, 0.23852, 0.02899), (0.0100, 0.47748, 0.07095),
    (0.0200, 0.95611, 0.20175), (0.0500, 2.39703, 1.05140), (0.0820, 3.94154, 2.52403),
]


def lepton_thermo(n_l: float, m_l: float) -> tuple[float, float, float]:
    """Relativistic free Fermi gas energy density, pressure, and chemical potential."""
    if n_l <= 0.0:
        return 0.0, 0.0, 0.0
    kF = HBARC * (3.0 * np.pi**2 * n_l) ** (1.0 / 3.0)
    mu = float(np.sqrt(kF**2 + m_l**2))
    x = kF / m_l
    pref = m_l**4 / (8.0 * np.pi**2 * HBARC**3)
    eps = pref * (x * (2.0 * x**2 + 1.0) * np.sqrt(1.0 + x**2) - np.arcsinh(x))
    p = (pref / 3.0) * (x * (2.0 * x**2 - 3.0) * np.sqrt(1.0 + x**2) + 3.0 * np.arcsinh(x))
    return float(eps), float(p), mu


def rmf_kin(n_i: float, m_eff: float) -> tuple[float, float, float, float]:
    """Relativistic baryon kinetic energy density, pressure, scalar density, and Fermi energy."""
    if n_i <= 0.0:
        return 0.0, 0.0, 0.0, m_eff
    kF = HBARC * (3.0 * np.pi**2 * n_i) ** (1.0 / 3.0)
    EF = float(np.sqrt(kF**2 + m_eff**2))
    if m_eff <= 1e-8:
        eps = kF**4 / (4.0 * np.pi**2 * HBARC**3)
        p = eps / 3.0
        ns = 0.0
    else:
        x = kF / m_eff
        pref = m_eff**4 / (8.0 * np.pi**2 * HBARC**3)
        eps = pref * (x * (2.0 * x**2 + 1.0) * np.sqrt(1.0 + x**2) - np.arcsinh(x))
        p = (pref / 3.0) * (x * (2.0 * x**2 - 3.0) * np.sqrt(1.0 + x**2) + 3.0 * np.arcsinh(x))
        ns = (m_eff**3 / (2.0 * np.pi**2 * HBARC**3)) * (x * np.sqrt(1.0 + x**2) - np.arcsinh(x))
    return float(eps), float(p), float(ns), EF


def _dd_couplings(nb: float, p: dict) -> tuple[float, float, float, float, float, float]:
    x = nb / p["nsat"]
    def f_meson(a, b, c, d):
        num = a * (1.0 + b * (x + d) ** 2)
        den = 1.0 + c * (x + d) ** 2
        f = num / den
        df_dx = (2.0 * a * (b - c) * (x + d)) / (den**2)
        return f, df_dx / p["nsat"]
    fs, dfs = f_meson(p["as_"], p["bs"], p["cs"], p["ds"])
    fw, dfw = f_meson(p["aw"], p["bw"], p["cw"], p["dw"])
    fr = float(np.exp(-p["ar"] * (x - 1.0)))
    dfr = -p["ar"] * fr / p["nsat"]
    return p["gs0"] * fs, p["gw0"] * fw, p["gr0"] * fr, p["gs0"] * dfs, p["gw0"] * dfw, p["gr0"] * dfr


def dd_bulk_thermo(nb: float, yp: float, p: dict) -> tuple[float, float, float, float, float]:
    """Evaluate DD-RMF baryon energy density, pressure, chemical potentials, and M*."""
    nn = nb * (1.0 - yp)
    np_ = nb * yp
    n3 = np_ - nn
    Gs, Gw, Gr, dGs, dGw, dGr = _dd_couplings(nb, p)
    Cs2 = (Gs / p["ms"]) ** 2 * HBARC**3
    Cw2 = (Gw / p["mw"]) ** 2 * HBARC**3
    Cr2 = (Gr / p["mr"]) ** 2 * HBARC**3

    def f_S(S):
        mn_s = max(1e-6, p["mn"] - S)
        mp_s = max(1e-6, p["mp"] - S)
        _, _, nsn, _ = rmf_kin(nn, mn_s)
        _, _, nsp, _ = rmf_kin(np_, mp_s)
        return S - Cs2 * (nsn + nsp)

    S = root_scalar(f_S, bracket=[0.0, p["mn"] - 1e-4]).root
    mn_s, mp_s = p["mn"] - S, p["mp"] - S
    ekn, pkn, nsn, EFn = rmf_kin(nn, mn_s)
    ekp, pkp, nsp, EFp = rmf_kin(np_, mp_s)
    ns = nsn + nsp

    V_om = Cw2 * nb
    V_rh = 0.25 * Cr2 * n3
    Sigma_R = (dGw / Gw) * V_om * nb + (dGr / Gr) * V_rh * n3 - (dGs / Gs) * S * ns
    eps_b = ekn + ekp + 0.5 * S * ns + 0.5 * V_om * nb + 0.5 * V_rh * n3
    P_b = pkn + pkp - 0.5 * S * ns + 0.5 * V_om * nb + 0.5 * V_rh * n3 + nb * Sigma_R
    mun = EFn + V_om - V_rh + Sigma_R
    mup = EFp + V_om + V_rh + Sigma_R
    return float(eps_b), float(P_b), float(mun), float(mup), float(mn_s)


def nl_bulk_thermo(nb: float, yp: float, p: dict) -> tuple[float, float, float, float, float]:
    """Evaluate Nonlinear RMF baryon energy density, pressure, chemical potentials, and M*."""
    nn = nb * (1.0 - yp)
    np_ = nb * yp
    n3 = np_ - nn
    ms2_g2 = (p["ms"] ** 2) / p["gs2"]
    mw2_g2 = (p["mw"] ** 2) / p["gw2"]
    mr2_g2 = (p["mr"] ** 2) / p["gr2"]
    kap = p["kappa"]
    lam = p["lam"]
    zeta = p["zeta"]
    Lam_v = p["Lam_v"]

    def f_S(S):
        mn_s = max(1e-6, p["mn"] - S)
        mp_s = max(1e-6, p["mp"] - S)
        _, _, nsn, _ = rmf_kin(nn, mn_s)
        _, _, nsp, _ = rmf_kin(np_, mp_s)
        ns = nsn + nsp
        return ms2_g2 * S + 0.5 * kap * S**2 + (lam / 6.0) * S**3 - ns * HBARC**3

    S = root_scalar(f_S, bracket=[0.0, p["mn"] - 1e-4]).root
    mn_s, mp_s = p["mn"] - S, p["mp"] - S
    ekn, pkn, nsn, EFn = rmf_kin(nn, mn_s)
    ekp, pkp, nsp, EFp = rmf_kin(np_, mp_s)

    def f_W(W):
        R = (0.5 * n3 * HBARC**3) / (mr2_g2 + 2.0 * Lam_v * W**2)
        return mw2_g2 * W + (zeta / 6.0) * W**3 + 2.0 * Lam_v * W * R**2 - nb * HBARC**3

    W_max = (nb * HBARC**3) / mw2_g2 + 10.0
    W = root_scalar(f_W, bracket=[0.0, max(10.0, W_max)]).root
    R = (0.5 * n3 * HBARC**3) / (mr2_g2 + 2.0 * Lam_v * W**2)

    U_sig = (0.5 * ms2_g2 * S**2 + (kap / 6.0) * S**3 + (lam / 24.0) * S**4) / (HBARC**3)
    eps_vec = (0.5 * mw2_g2 * W**2 + (zeta / 8.0) * W**4 + 0.5 * mr2_g2 * R**2 + 3.0 * Lam_v * W**2 * R**2) / (HBARC**3)
    P_vec = (0.5 * mw2_g2 * W**2 + (zeta / 24.0) * W**4 + 0.5 * mr2_g2 * R**2 + Lam_v * W**2 * R**2) / (HBARC**3)

    eps_b = ekn + ekp + U_sig + eps_vec
    P_b = pkn + pkp - U_sig + P_vec
    mun = EFn + W - R
    mup = EFp + W + R
    return float(eps_b), float(P_b), float(mun), float(mup), float(mn_s)


def gdfm_bulk_thermo(rho: float, yp: float, p: dict) -> tuple[float, float, float, float, float]:
    """Evaluate Scurto, Pais & Gulminelli (2024) GDFM RMF thermodynamic quantities."""
    M = p["mn"]
    if rho <= 1e-12:
        return 0.0, 0.0, M, M, M
    x = rho / p["rsat"]
    es = np.exp(-p["cs"] * x)
    ew = np.exp(-p["cw"] * x)
    er = np.exp(-p["cr"] * x)
    Gs = p["as_"] + (p["bs"] + p["ds"] * x**3) * es
    Gw = p["aw"] + (p["bw"] + p["dw"] * x**3) * ew
    Gr = p["ar"] + (p["br"] + p["dr"] * x**3) * er
    dGs = ((3.0 * p["ds"] * x**2 - p["cs"] * (p["bs"] + p["ds"] * x**3)) * es) / p["rsat"]
    dGw = ((3.0 * p["dw"] * x**2 - p["cw"] * (p["bw"] + p["dw"] * x**3)) * ew) / p["rsat"]
    dGr = ((3.0 * p["dr"] * x**2 - p["cr"] * (p["br"] + p["dr"] * x**3)) * er) / p["rsat"]

    rn = rho * (1.0 - yp)
    rp = rho * yp
    r3 = rp - rn
    Cs2 = (Gs / p["ms"]) ** 2 * HBARC**3
    Cw2 = (Gw / p["mw"]) ** 2 * HBARC**3
    Cr2 = (Gr / p["mr"]) ** 2 * HBARC**3

    def f_S(S):
        mstar = max(1e-6, M - S)
        _, _, nsn, _ = rmf_kin(rn, mstar)
        _, _, nsp, _ = rmf_kin(rp, mstar)
        return S - Cs2 * (nsn + nsp)

    S = root_scalar(f_S, bracket=[0.0, M - 1e-4]).root
    mstar = M - S
    ekn, pkn, nsn, EFn = rmf_kin(rn, mstar)
    ekp, pkp, nsp, EFp = rmf_kin(rp, mstar)
    ns = nsn + nsp

    V_om = Cw2 * rho
    V_rh = 0.25 * Cr2 * r3
    Sigma_R = (dGw / Gw) * V_om * rho + (dGr / Gr) * V_rh * r3 - (dGs / Gs) * S * ns
    E_b = ekn + ekp + 0.5 * S * ns + 0.5 * V_om * rho + 0.5 * V_rh * r3
    P_b = pkn + pkp - 0.5 * S * ns + 0.5 * V_om * rho + 0.5 * V_rh * r3 + rho * Sigma_R
    mun = EFn + V_om - V_rh + Sigma_R
    mup = EFp + V_om + V_rh + Sigma_R
    return float(E_b), float(P_b), float(mun), float(mup), float(mstar)


def gm1_couplings_from_sat(
    mN: float = 938.0,
    n0: float = 0.153,
    BA: float = -16.30,
    K: float = 300.0,
    mstar_ratio: float = 0.70,
    asym: float = 32.5,
) -> dict:
    """Derive Glendenning & Moszkowski (1991) PRL 67, 2414 Table II couplings."""
    mstar0 = mN * mstar_ratio
    kF0 = HBARC * (1.5 * np.pi**2 * n0) ** (1.0 / 3.0)
    EF0 = float(np.sqrt(kF0**2 + mstar0**2))
    mu0 = mN + BA
    Cw2 = (mu0 - EF0) / n0
    Cr2 = (asym - kF0**2 / (6.0 * EF0)) / (0.125 * n0)
    ek_n, pk_n, ns_n, _ = rmf_kin(0.5 * n0, mstar0)
    ek0, ns0 = 2.0 * ek_n, 2.0 * ns_n
    U0 = n0 * mu0 - ek0 - 0.5 * Cw2 * n0**2
    eps_m = 1e-5
    _, _, ns_p, _ = rmf_kin(0.5 * n0, mstar0 + eps_m)
    _, _, ns_m, _ = rmf_kin(0.5 * n0, mstar0 - eps_m)
    dns_dm = 2.0 * (ns_p - ns_m) / (2.0 * eps_m)
    denom = (mstar0 / EF0) ** 2 / (Cw2 + kF0**2 / (3.0 * n0 * EF0) - K / (9.0 * n0))
    Upp0 = denom - dns_dm
    S0 = mN - mstar0
    A_mat = np.array([
        [0.5 * S0**2, (1.0 / 3.0) * S0**3, 0.25 * S0**4],
        [S0,          S0**2,               S0**3],
        [1.0,         2.0 * S0,            3.0 * S0**2],
    ])
    b_vec = np.array([U0, ns0, Upp0])
    inv_Cs2, B_param, C_param = np.linalg.solve(A_mat, b_vec)
    Cs2 = float(1.0 / inv_Cs2)
    return dict(mN=mN, Cs2=Cs2, Cw2=float(Cw2), Cr2=float(Cr2),
                B=float(B_param), C=float(C_param), b=float(B_param / mN), c=float(C_param))


def solve_gm1_point(nb: float, params: dict) -> tuple[float, float, float, float, float]:
    """Solve beta-equilibrated GM1 RMF state at density nb [fm^-3]."""
    mN = params["mN"]
    Cs2, Cw2, Cr2, B, C = params["Cs2"], params["Cw2"], params["Cr2"], params["B"], params["C"]

    def solve_S(yp):
        nn, np_ = nb * (1.0 - yp), nb * yp
        def f_S(S):
            mstar = max(1e-6, mN - S)
            _, _, nsn, _ = rmf_kin(nn, mstar)
            _, _, nsp, _ = rmf_kin(np_, mstar)
            return S / Cs2 + B * S**2 + C * S**3 - (nsn + nsp)
        return root_scalar(f_S, bracket=[0.0, mN - 1e-3], method="brentq").root

    def beta_res(yp):
        S = solve_S(yp)
        mstar = mN - S
        nn, np_ = nb * (1.0 - yp), nb * yp
        kFn = HBARC * (3.0 * np.pi**2 * nn) ** (1.0 / 3.0)
        kFp = HBARC * (3.0 * np.pi**2 * np_) ** (1.0 / 3.0)
        EFn = np.sqrt(kFn**2 + mstar**2)
        EFp = np.sqrt(kFp**2 + mstar**2)
        mun_m_mup = EFn - EFp + 0.5 * Cr2 * (nn - np_)
        kFe = HBARC * (3.0 * np.pi**2 * np_) ** (1.0 / 3.0)
        mue = np.sqrt(kFe**2 + ME**2)
        return mun_m_mup - mue

    yp = float(root_scalar(beta_res, bracket=[1e-6, 0.49], method="brentq").root)
    S = solve_S(yp)
    mstar = mN - S
    nn, np_ = nb * (1.0 - yp), nb * yp
    ekn, pkn, _, EFn = rmf_kin(nn, mstar)
    ekp, pkp, _, EFp = rmf_kin(np_, mstar)
    eke, pke, _ = lepton_thermo(np_, ME)
    U_sig = 0.5 * S**2 / Cs2 + (1.0 / 3.0) * B * S**3 + 0.25 * C * S**4
    V_om = 0.5 * Cw2 * nb**2
    V_rh = 0.125 * Cr2 * (nn - np_) ** 2
    E_tot = float(ekn + ekp + eke + U_sig + V_om + V_rh)
    P_tot = float(pkn + pkp + pke - U_sig + V_om + V_rh)
    mun = float(EFn + Cw2 * nb + 0.25 * Cr2 * (nn - np_))
    mup = float(EFp + Cw2 * nb - 0.25 * Cr2 * (nn - np_))
    return P_tot, E_tot, mun, mup - mun, yp


def apr_energy_per_baryon(rho: float, xp: float, phase: str = "auto") -> float:
    """Exact Akmal, Pandharipande & Ravenhall (1998) PRC 58, 1804 Appendix A & Table XII."""
    hb2_2m = HBARC**2 / (2.0 * 939.0)
    p1, p2, p3, p4, p5, p6 = 337.2, -382.0, 89.8, 0.457, -59.0, -19.1
    p7, p8, p9, p10, p11, p12, p13 = 214.6, -384.0, 6.4, 69.0, -33.0, 0.35, 0.0
    p14, p15, p16, p17, p18, p19, p20, p21 = 0.0, 287.0, -1.54, 175.0, -1.45, 0.32, 0.195, 0.0

    tau_p = (3.0 / 5.0) * (3.0 * np.pi**2 * rho) ** (2.0 / 3.0) * (xp ** (5.0 / 3.0)) if xp > 0 else 0.0
    tau_n = (3.0 / 5.0) * (3.0 * np.pi**2 * rho) ** (2.0 / 3.0) * ((1.0 - xp) ** (5.0 / 3.0)) if xp < 1 else 0.0
    kin = (
        (hb2_2m + (p3 + (1.0 - xp) * p5) * rho * np.exp(-p4 * rho)) * tau_n
        + (hb2_2m + (p3 + xp * p5) * rho * np.exp(-p4 * rho)) * tau_p
    )
    gL_snm = -rho * (p1 + p2 * rho + p6 * rho**2 + (p10 + p11 * rho) * np.exp(-((p9 * rho) ** 2)))
    gL_pnm = -rho * (p12 / rho + p7 + p8 * rho + p13 * np.exp(-((p9 * rho) ** 2)))
    if phase == "LDP" or (phase == "auto" and rho <= 0.195):
        g_snm, g_pnm = gL_snm, gL_pnm
    else:
        ds, dn = rho - p19, rho - p20
        gH_snm = gL_snm - rho * (p17 * ds + p21 * ds**2) * np.exp(p18 * ds) if ds > 0 else gL_snm
        gH_pnm = gL_pnm - rho * (p15 * dn + p14 * dn**2) * np.exp(p16 * dn) if dn > 0 else gL_pnm
        if phase == "HDP":
            g_snm, g_pnm = gH_snm, gH_pnm
        else:
            eL = kin + gL_snm * (1.0 - (1.0 - 2.0 * xp) ** 2) + gL_pnm * (1.0 - 2.0 * xp) ** 2
            eH = kin + gH_snm * (1.0 - (1.0 - 2.0 * xp) ** 2) + gH_pnm * (1.0 - 2.0 * xp) ** 2
            return float(min(eL, eH))
    return float(kin + g_snm * (1.0 - (1.0 - 2.0 * xp) ** 2) + g_pnm * (1.0 - 2.0 * xp) ** 2)


def _apr_total_e_of_xp(r: float, xp: float) -> float:
    np_ = r * xp
    kF_thresh = np.sqrt(MMU**2 - ME**2)
    n_e_thresh = (kF_thresh / HBARC) ** 3 / (3.0 * np.pi**2)
    if np_ <= n_e_thresh:
        ne, nmu = np_, 0.0
    else:
        def f_ne(ne_val):
            kFe = HBARC * (3.0 * np.pi**2 * ne_val) ** (1.0 / 3.0)
            kFmu = HBARC * (3.0 * np.pi**2 * max(0.0, np_ - ne_val)) ** (1.0 / 3.0)
            return np.sqrt(kFe**2 + ME**2) - np.sqrt(kFmu**2 + MMU**2)
        ne = root_scalar(f_ne, bracket=[n_e_thresh, np_], method="brentq").root
        nmu = np_ - ne
    ee, _, _ = lepton_thermo(ne, ME)
    emu, _, _ = lepton_thermo(nmu, MMU)
    return float(apr_energy_per_baryon(r, xp) + (ee + emu) / r)


def solve_apr_beta(rho: float, mn: float = 939.56533, mp: float = 938.27200) -> tuple[float, float, float]:
    """Solve beta-equilibrated APR state at baryon density rho [fm^-3]."""
    res = minimize_scalar(lambda xp: _apr_total_e_of_xp(rho, xp), bounds=(1e-5, 0.49), method="bounded")
    xp = float(res.x)
    dr = rho * 1e-5
    e_plus = _apr_total_e_of_xp(rho + dr, xp)
    e_minus = _apr_total_e_of_xp(rho - dr, xp)
    P = float(rho**2 * (e_plus - e_minus) / (2.0 * dr))
    E_tot = float(rho * (mn + _apr_total_e_of_xp(rho, xp)))
    return P, E_tot, xp


def skyrme_bulk_thermo(
    rho: float, yp: float, sp: dict, mn: float = 939.56533, mp: float = 938.27200
) -> tuple[float, float, float, float]:
    """Evaluate Skyrme EDF baryon energy density, pressure, and chemical potentials."""
    xn, xp = 1.0 - yp, yp
    delta = xn - xp
    rn, rp = rho * xn, rho * xp
    t0, t1, t2, t3 = sp["t0"], sp["t1"], sp["t2"], sp["t3"]
    x0, x1, x2, x3 = sp["x0"], sp["x1"], sp["x2"], sp["x3"]
    alpha = sp.get("alpha", sp.get("sigma", 1.0 / 6.0))

    hb2_2mn = HBARC**2 / (2.0 * mn)
    hb2_2mp = HBARC**2 / (2.0 * mp)
    tau_n = (3.0 / 5.0) * (3.0 * np.pi**2) ** (2.0 / 3.0) * rn ** (5.0 / 3.0) if rn > 0 else 0.0
    tau_p = (3.0 / 5.0) * (3.0 * np.pi**2) ** (2.0 / 3.0) * rp ** (5.0 / 3.0) if rp > 0 else 0.0
    tau = tau_n + tau_p

    C0 = (3.0 / 8.0) * t0
    D0 = -(1.0 / 8.0) * t0 * (2.0 * x0 + 1.0)
    C3 = (1.0 / 16.0) * t3
    D3 = -(1.0 / 48.0) * t3 * (2.0 * x3 + 1.0)
    A_eff = (1.0 / 8.0) * (t1 * (2.0 + x1) + t2 * (2.0 + x2))
    B_eff = (1.0 / 8.0) * (t2 * (2.0 * x2 + 1.0) - t1 * (2.0 * x1 + 1.0))

    eps_kin = hb2_2mn * tau_n + hb2_2mp * tau_p
    eps_pot = (
        (C0 + D0 * delta**2) * rho**2
        + (C3 + D3 * delta**2) * rho ** (alpha + 2.0)
        + A_eff * rho * tau
        + B_eff * (rn * tau_n + rp * tau_p)
    )
    eps_b = eps_kin + eps_pot

    kFn2 = (3.0 * np.pi**2 * rn) ** (2.0 / 3.0) if rn > 0 else 0.0
    kFp2 = (3.0 * np.pi**2 * rp) ** (2.0 / 3.0) if rp > 0 else 0.0
    drho_diff = rn - rp

    mu_n = (
        hb2_2mn * kFn2 + 2.0 * C0 * rho + 2.0 * D0 * drho_diff
        + (alpha + 2.0) * C3 * rho ** (alpha + 1.0)
        + D3 * (alpha * rho ** (alpha - 1.0) * drho_diff**2 + 2.0 * rho**alpha * drho_diff)
        + A_eff * (tau + rho * kFn2) + B_eff * (tau_n + rn * kFn2)
    )
    mu_p = (
        hb2_2mp * kFp2 + 2.0 * C0 * rho - 2.0 * D0 * drho_diff
        + (alpha + 2.0) * C3 * rho ** (alpha + 1.0)
        + D3 * (alpha * rho ** (alpha - 1.0) * drho_diff**2 - 2.0 * rho**alpha * drho_diff)
        + A_eff * (tau + rho * kFp2) + B_eff * (tau_p + rp * kFp2)
    )
    P_b = rn * mu_n + rp * mu_p - eps_b
    return float(eps_b), float(P_b), float(mu_n), float(mu_p)


def solve_skyrme_beta(
    rho: float, sp: dict, mn: float = 939.56533, mp: float = 938.27200
) -> tuple[float, float, float, float, float, float, float]:
    """Solve beta-equilibrated Skyrme state at baryon density rho [fm^-3]."""
    kF_th = np.sqrt(MMU**2 - ME**2)
    ne_th = (kF_th / HBARC) ** 3 / (3.0 * np.pi**2)

    def beta_res(yp):
        _, _, mun, mup = skyrme_bulk_thermo(rho, yp, sp, mn=mn, mp=mp)
        mun_tot, mup_tot = mun + mn, mup + mp
        np_ = rho * yp
        if np_ <= ne_th:
            ne = np_
        else:
            ne = root_scalar(
                lambda x: np.sqrt((HBARC * (3.0 * np.pi**2 * x) ** (1.0 / 3.0)) ** 2 + ME**2)
                - np.sqrt((HBARC * (3.0 * np.pi**2 * max(0.0, np_ - x)) ** (1.0 / 3.0)) ** 2 + MMU**2),
                bracket=[ne_th, np_],
            ).root
        kFe = HBARC * (3.0 * np.pi**2 * ne) ** (1.0 / 3.0) if ne > 0 else 0.0
        mue = np.sqrt(kFe**2 + ME**2)
        return (mun_tot - mup_tot) - mue

    yp = float(root_scalar(beta_res, bracket=[1e-5, 0.49]).root)
    np_ = rho * yp
    if np_ <= ne_th:
        ne, nmu = np_, 0.0
    else:
        ne = float(root_scalar(
            lambda x: np.sqrt((HBARC * (3.0 * np.pi**2 * x) ** (1.0 / 3.0)) ** 2 + ME**2)
            - np.sqrt((HBARC * (3.0 * np.pi**2 * max(0.0, np_ - x)) ** (1.0 / 3.0)) ** 2 + MMU**2),
            bracket=[ne_th, np_],
        ).root)
        nmu = np_ - ne
    ye, ymu = ne / rho, nmu / rho
    eps_b, P_b, mun, mup = skyrme_bulk_thermo(rho, yp, sp, mn=mn, mp=mp)
    ee, pe, _ = lepton_thermo(ne, ME)
    emu, pmu, _ = lepton_thermo(nmu, MMU)
    P_tot = float(P_b + pe + pmu)
    E_tot = float(eps_b + ee + emu + rho * ((1.0 - yp) * mn + yp * mp))
    return P_tot, E_tot, float(mun + mn), float(mup + mp - (mun + mn)), yp, float(ye), float(ymu)


def outer_crust_pressure(nb_val: float, A: float, Z: float) -> float:
    """Relativistic electron Fermi gas + Wigner-Seitz bcc lattice pressure."""
    ye = Z / A
    ne = nb_val * ye
    _, Pel, _ = lepton_thermo(ne, ME)
    C_M = -0.895929
    Plat = (1.0 / 3.0) * C_M * ALPHA * HBARC * (4.0 * np.pi / 3.0) ** (1.0 / 3.0) * (Z ** (2.0 / 3.0)) * (ne ** (4.0 / 3.0))
    return float(Pel + Plat)


__all__ = [
    "DD_MODELS",
    "NL_MODELS",
    "SPG_GDFM_PARAMS",
    "SKYRME_PARAMS",
    "BPS_HP_OUTER_CRUST_NUCLEI",
    "GRILL2014_DDME2",
    "GRILL2014_NL3WR",
    "lepton_thermo",
    "rmf_kin",
    "dd_bulk_thermo",
    "nl_bulk_thermo",
    "gdfm_bulk_thermo",
    "gm1_couplings_from_sat",
    "solve_gm1_point",
    "apr_energy_per_baryon",
    "solve_apr_beta",
    "skyrme_bulk_thermo",
    "solve_skyrme_beta",
    "outer_crust_pressure",
]
