"""Exploratory spherical-droplet CLD proxy for inner-crust experiments.

The bulk energies come from rmf_solver.  The explicitly chosen interface
functions below are historical proxies, not Scurto et al. (2024) eqs. 24-25.
In particular, sigma_0 and b have not been fitted to measured nuclear masses.
Source-faithful crust predictions require those functions and fitted inputs.

The total cell energy is minimized over f, R_d, rho_I and rho_Ip, with pure
neutron gas and uniform electrons.  The conditional rho_Ip solve includes
surface/Coulomb derivatives and the full nucleon rest masses.  Units: fm^-3,
MeV, MeV fm^-3 and fm.  Output E_cell excludes nucleon rest masses only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from .rmf_solver import (
    HBARC, M_E, M_MU, M_N, M_P, RMFModel, solve_homogeneous, energy_per_baryon,
    HomogeneousResult,
)

ALPHA = 1.0 / 137.036   # fine-structure constant
E2 = ALPHA * HBARC      # e^2 in MeV*fm = 1.439965 MeV*fm


def coulomb_lattice_phi(u: float) -> float:
    """Historical proxy Phi(u)=sqrt(2-3u)/2 on 0<u<2/3.

    This is not the spherical WS result (2-3*u**(1/3)+u)/5.  Kept explicitly
    for reproducibility of the exploratory functional; no source calibration.
    """
    if u <= 0 or u >= 2.0 / 3.0:
        return 0.0
    return 0.5 * math.sqrt(2.0 - 3.0 * u)


def surface_tension(yp_I: float, sigma_0: float, b: float) -> float:
    """Historical proxy sigma=sigma_0*[yp*(1-yp)]**b.

    This is not Scurto eq. 25 or a fitted Lattimer-Swesty surface tension.
    b is an exponent in this proxy, not the source's denominator parameter.
    """
    if yp_I <= 0 or yp_I >= 1:
        return 0.0
    # Empirical form: sigma peaks at y_p=0.5, vanishes at y_p=0,1
    return sigma_0 * (yp_I * (1.0 - yp_I)) ** b


def electron_energy_density(rho_e: float, mass: float = M_E) -> float:
    """Relativistic degenerate electron energy density (MeV/fm^3) at T=0.
    E_L = (1/(pi^2 hbar^3)) ∫_0^kFe k^2 sqrt(k^2 + m_e^2 c^4) dk  (in natural units)
    In our convention: rho_e in fm^-3, returns MeV/fm^3.
    """
    if rho_e <= 0:
        return 0.0
    kFe = (3.0 * math.pi ** 2 * rho_e) ** (1.0 / 3.0)  # fm^-1
    kFe_MeV = kFe * HBARC  # MeV
    E = math.sqrt(kFe_MeV ** 2 + mass ** 2)
    # integral in MeV^4 (degeneracy 2 for spin-1/2)
    I = 0.125 * (kFe_MeV * (2 * kFe_MeV ** 2 + mass ** 2) * E
                 - mass ** 4 * math.log((kFe_MeV + E) / mass))
    e_MeV4 = (2.0 / (2.0 * math.pi ** 2)) * I
    return e_MeV4 / HBARC ** 3


def electron_chemical_potential(rho_e: float, mass: float = M_E) -> float:
    """Electron chemical potential (Fermi energy, MeV) at T=0."""
    if rho_e <= 0:
        return mass
    kFe = (3.0 * math.pi ** 2 * rho_e) ** (1.0 / 3.0) * HBARC  # MeV
    return math.sqrt(kFe ** 2 + mass ** 2)


@dataclass
class CLDCellResult:
    """Result of the CLD cell minimization at fixed total baryon density rho_B."""
    rho_B: float            # total baryon density fm^-3
    f: float                # cluster volume fraction
    R_d: float              # cluster radius fm
    rho_I: float            # cluster baryon density fm^-3
    rho_Ip: float           # cluster proton density fm^-3
    rho_II: float           # gas baryon density fm^-3
    rho_IIp: float          # gas proton density fm^-3
    rho_e: float            # electron density fm^-3
    y_p_cell: float         # cell-averaged proton fraction
    E_cell: float           # cell energy density MeV/fm^3 (incl. electrons)
    E_binding_per_baryon: float  # (E_cell - M*rho_B)/rho_B MeV
    E_coul: float           # Coulomb energy density MeV/fm^3
    E_surf: float           # surface energy density MeV/fm^3
    E_L: float              # electron energy density MeV/fm^3
    P_cell: float           # cell pressure MeV/fm^3
    mu_n: float; mu_p: float; mu_e: float
    converged: bool
    n_iter: int
    message: str = ""
    # Neutron chemical potential of the uniform gas phase.  At the variational
    # optimum the cell mu_B equals this value; the cluster *bulk* mu_n differs
    # by the surface and Coulomb derivative terms.
    mu_n_gas: float = float("nan")


def cld_cell_energy(model: RMFModel, rho_B: float, f: float, R_d: float,
                    rho_I: float, rho_Ip: float, rho_II: float, rho_IIp: float,
                    sigma_0: float, b: float) -> tuple[float, dict]:
    """Compute the CLD cell energy density at fixed total density rho_B and
    the four variational variables (f, R_d, rho_I, rho_Ip); the gas density
    rho_II, rho_IIp are constrained by the baryon-number conservation
    (rho_B = f*rho_I + (1-f)*rho_II) and the cell charge neutrality
    (f*rho_Ip + (1-f)*rho_IIp = rho_e).

    Returns (E_cell_MeV_fm3, info_dict)."""
    # baryon conservation gives rho_II from rho_B, f, rho_I:
    if abs(f) < 1e-12:
        rho_II_eff = rho_B
    elif abs(f - 1) < 1e-12:
        rho_II_eff = 0.0
    else:
        rho_II_eff = (rho_B - f * rho_I) / (1.0 - f)
    # clamp to physical
    rho_II_eff = max(0.0, rho_II_eff)
    # gas proton density from charge neutrality: f*rho_Ip + (1-f)*rho_IIp = rho_e
    # and beta-equilibrium: rho_e is determined by mu_e = mu_n - mu_p at the gas.
    # For the energy evaluation we still need rho_IIp; use the cell-averaged
    # proton fraction from the cluster and assume beta-equilibrated gas:
    # we treat rho_IIp as a free variable minimized alongside, but for a first
    # implementation set rho_IIp from the gas neutron-rich limit (rho_IIp ~ 0 at
    # low density) and let the outer minimizer vary it. Here we just use the
    # provided rho_IIp.
    rho_e = f * rho_Ip + (1.0 - f) * rho_IIp
    rho_e = max(0.0, rho_e)

    # homogeneous-matter energies (binding) for cluster and gas
    rI = solve_homogeneous(model, rho_I, rho_Ip)
    rII = solve_homogeneous(model, rho_II_eff, rho_IIp)
    E_I = rI.binding_energy       # MeV/fm^3 (full minus rest mass)
    E_II = rII.binding_energy
    # electron background
    E_L = electron_energy_density(rho_e)
    # Coulomb geometry structure with the declared proxy Phi: E_Coul = 2*pi*e^2 * f * Phi * R_d^2 * (rho_Ip - rho_IIp)^2
    Phi = coulomb_lattice_phi(f)
    rho_p_diff = rho_Ip - rho_IIp
    E_coul = 2.0 * math.pi * E2 * f * Phi * R_d ** 2 * rho_p_diff ** 2
    # Surface geometry structure with the declared proxy sigma: E_surf = (3*f*sigma)/R_d
    yp_I = rho_Ip / rho_I if rho_I > 0 else 0.0
    sigma = surface_tension(yp_I, sigma_0, b)
    E_surf = (3.0 * f * sigma) / R_d if R_d > 0 else 0.0
    # cell energy density (binding + rest mass recovered as M*rho_B for nucleons)
    # E_cell = f*E_I_full + (1-f)*E_II_full + E_coul + E_surf + E_L
    # where E_I_full = E_I + M_P*rho_Ip + M_N*rho_In (full including rest mass)
    # and the total rest mass of the cell = M_P*(f*rho_Ip + (1-f)*rho_IIp) + M_N*(rho_B - that)
    # so E_cell_binding = f*E_I + (1-f)*E_II + E_coul + E_surf + E_L  (already binding)
    E_cell = f * E_I + (1.0 - f) * E_II + E_coul + E_surf + E_L

    info = {
        "rho_II": rho_II_eff, "rho_e": rho_e, "E_I": E_I, "E_II": E_II,
        "E_coul": E_coul, "E_surf": E_surf, "E_L": E_L, "sigma": sigma,
        "Phi": Phi, "rI_converged": rI.converged, "rII_converged": rII.converged,
        "Mstar_I_p": rI.Mstar_p, "Mstar_II_p": rII.Mstar_p,
    }
    return E_cell, info


def solve_beta_equilibrium_cluster(model: RMFModel, rho_I: float,
                                   rho_IIp_gas: float, f: float, *, R_d: float,
                                   sigma_0: float, b: float) -> tuple[float, float, float, float]:
    """Solve d eps_cell/d rho_Ip=0 at fixed f,R_d,rho_I and pure neutron gas.

    Dividing by f gives mu_p(I)-mu_n(I)+mu_e + 4*pi*e^2*Phi*R_d^2*rho_Ip
    + 3*sigma'(yp)/(R_d*rho_I).  Full homogeneous chemical potentials include
    the neutron/proton rest-mass difference.  The last two terms are essential
    even when bulk beta equilibrium is very accurately solved.
    """
    rho_IIp_gas = 0.0  # pure neutron gas (low-density inner-crust approximation)
    def residual(rho_Ip):
        rho_Ip = max(1e-12*rho_I, min(rho_Ip, (1.0-1e-12)*rho_I))
        rI = solve_homogeneous(model, rho_I, rho_Ip)
        rho_e = f * rho_Ip + (1.0 - f) * rho_IIp_gas
        rho_e = max(0.0, rho_e)
        mu_e = electron_chemical_potential(rho_e) if rho_e > 0 else 0.511
        yp = rho_Ip / rho_I
        dsigma = b * surface_tension(yp, sigma_0, b) * (1.0 - 2.0*yp) / (yp*(1.0-yp))
        interface = (4.0 * math.pi * E2 * coulomb_lattice_phi(f) * R_d**2 * rho_Ip
                     + 3.0*dsigma/(R_d*rho_I))
        return (rI.mu_n - rI.mu_p) - mu_e - interface
    if b < 1.0:
        # Fractional powers have singular endpoint derivatives: the endpoint
        # sign bracket alone can select a maximum near yp=1.  Isolate observed
        # sign changes on a uniform grid plus logarithmic endpoint refinement,
        # then select by the FULL energy.  This is a numerical search, not a
        # proof of global uniqueness for arbitrary RMF/surface parameters.
        from scipy.optimize import brentq
        ys = sorted(set([i/64.0 for i in range(65)] +
                        [10.0**(-k) for k in range(2, 13)] +
                        [1.0-10.0**(-k) for k in range(2, 13)]))
        vals = [residual(y*rho_I) for y in ys]
        candidates = [0.0, rho_I]
        for y0, y1, v0, v1 in zip(ys[:-1], ys[1:], vals[:-1], vals[1:]):
            if v0 == 0.0:
                candidates.append(y0*rho_I)
            if v0*v1 < 0.0:
                root = brentq(lambda y: residual(y*rho_I), y0, y1,
                              xtol=5e-324, rtol=1e-14)
                candidates.append(root*rho_I)
        def energy(p):
            r = solve_homogeneous(model, rho_I, p)
            y = p/rho_I
            return (f*(r.binding_energy + (model.mp-model.mn)*p)
                    + electron_energy_density(f*p)
                    + 2*math.pi*E2*f*coulomb_lattice_phi(f)*R_d**2*p**2
                    + 3*f*surface_tension(y, sigma_0, b)/R_d)
        p = min(candidates, key=energy)
        ne = f*p
        return p, ne, electron_chemical_potential(ne), abs(residual(p))
    r_lo = residual(0.0)
    r_hi = residual(rho_I)
    if r_lo * r_hi > 0:
        # no sign change: return endpoint with smaller |residual|
        if abs(r_lo) < abs(r_hi):
            return 0.0, f * 0.0 + (1 - f) * rho_IIp_gas, 0.511, abs(r_lo)
        rho_Ip = rho_I; rho_e = f * rho_Ip + (1 - f) * rho_IIp_gas
        return rho_Ip, max(0.0, rho_e), electron_chemical_potential(max(0.0, rho_e)), abs(r_hi)
    lo, hi = 0.0, rho_I
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        r_mid = residual(mid)
        if r_mid == 0 or (hi - lo) < 1e-12:
            break
        if r_lo * r_mid < 0:
            hi = mid; r_hi = r_mid
        else:
            lo = mid; r_lo = r_mid
    rho_Ip = 0.5 * (lo + hi)
    rho_e = max(0.0, f * rho_Ip + (1.0 - f) * rho_IIp_gas)
    mu_e = electron_chemical_potential(rho_e) if rho_e > 0 else 0.511
    return rho_Ip, rho_e, mu_e, abs(residual(rho_Ip))


def _cld_objective(model: RMFModel, rho_B: float, sigma_0: float, b: float,
                   x: tuple[float, float, float]) -> tuple[float, dict]:
    """Objective for the Nelder-Mead: x = (f, log R_d, rho_I).
    The cluster proton density rho_Ip solves the full cell beta-equilibrium
    condition including interface terms, with the gas treated as pure neutron
    (rho_IIp = 0) - the standard low-density inner-crust approximation.
    Returns (E_cell, info)."""
    f, log_Rd, rho_I = x
    f = max(1e-4, min(0.95, f))
    R_d = math.exp(max(math.log(2.0), min(math.log(40.0), log_Rd)))
    rho_I = max(0.01, min(0.5, rho_I))
    if f < 1:
        rho_II = (rho_B - f * rho_I) / (1.0 - f)
    else:
        rho_II = 0.0
    if rho_II < 0 or rho_II > rho_I:
        return float("inf"), {}
    # beta-equilibrium at the cluster: solve for rho_Ip
    rho_Ip, rho_e, mu_e, beta_resid = solve_beta_equilibrium_cluster(model, rho_I, 0.0, f,
        R_d=R_d, sigma_0=sigma_0, b=b)
    rho_IIp = 0.0  # pure neutron gas
    try:
        E, info = cld_cell_energy(model, rho_B, f, R_d, rho_I, rho_Ip, rho_II, rho_IIp,
                                  sigma_0, b)
    except (ValueError, ZeroDivisionError):
        return float("inf"), {}
    if not math.isfinite(E):
        return float("inf"), {}
    info["rho_Ip"] = rho_Ip; info["rho_IIp"] = rho_IIp
    info["rho_e"] = rho_e; info["mu_e"] = mu_e; info["beta_resid"] = beta_resid
    if not info["rI_converged"] or not info["rII_converged"] or beta_resid > 1e-6:
        return float("inf"), {}
    # At fixed n_B the mn*n_B term is constant.  The proton rest-mass
    # correction is not: varying charge changes the species mixture.
    return E + (model.mp - model.mn) * rho_e, info


def _nelder_mead(model, rho_B, sigma_0, b, x0, max_iter=500, tol=1e-10):
    """Lightweight Nelder-Mead (no scipy) in 3 variables (f, log R_d, rho_I).
    Returns (best_x, best_E, best_info, n_iter, converged)."""
    n = 3
    simplex = [list(x0)]
    steps = [0.02, math.log(1.5), 0.02 * model.rho_sat]
    for i in range(n):
        xi = list(x0); xi[i] += steps[i]; simplex.append(xi)
    fvals = [None] * len(simplex)
    infos = [None] * len(simplex)
    for i in range(len(simplex)):
        fvals[i], infos[i] = _cld_objective(model, rho_B, sigma_0, b, tuple(simplex[i]))
    alpha, gamma, rho_coef, sigma_coef = 1.0, 2.0, 0.5, 0.5
    converged = False
    for it in range(max_iter):
        order = sorted(range(len(simplex)), key=lambda k: fvals[k])
        simplex = [simplex[k] for k in order]
        fvals = [fvals[k] for k in order]
        infos = [infos[k] for k in order]
        spread = max(abs(simplex[j][i] - simplex[0][i]) for j in range(1, n+1) for i in range(n))
        if abs(fvals[-1] - fvals[0]) < tol and fvals[0] != float("inf") and spread < 1e-6:
            converged = True
            break
        centroid = [sum(simplex[j][i] for j in range(n)) / n for i in range(n)]
        xr = [centroid[i] + alpha * (centroid[i] - simplex[-1][i]) for i in range(n)]
        fr, ir = _cld_objective(model, rho_B, sigma_0, b, tuple(xr))
        if fvals[0] <= fr < fvals[-2]:
            simplex[-1], fvals[-1], infos[-1] = xr, fr, ir
            continue
        if fr < fvals[0]:
            xe = [centroid[i] + gamma * (xr[i] - centroid[i]) for i in range(n)]
            fe, ie = _cld_objective(model, rho_B, sigma_0, b, tuple(xe))
            if fe < fr:
                simplex[-1], fvals[-1], infos[-1] = xe, fe, ie
            else:
                simplex[-1], fvals[-1], infos[-1] = xr, fr, ir
            continue
        xc = [centroid[i] + rho_coef * (simplex[-1][i] - centroid[i]) for i in range(n)]
        fc, ic = _cld_objective(model, rho_B, sigma_0, b, tuple(xc))
        if fc < fvals[-1]:
            simplex[-1], fvals[-1], infos[-1] = xc, fc, ic
            continue
        for j in range(1, len(simplex)):
            simplex[j] = [simplex[0][i] + sigma_coef * (simplex[j][i] - simplex[0][i]) for i in range(n)]
            fvals[j], infos[j] = _cld_objective(model, rho_B, sigma_0, b, tuple(simplex[j]))
    bi = min(range(len(simplex)), key=lambda k: fvals[k])
    return simplex[bi], fvals[bi], infos[bi], it + 1, converged


def minimize_cld_cell(model: RMFModel, rho_B: float, sigma_0: float, b: float,
                     verbose: bool = False) -> Optional[CLDCellResult]:
    """Minimize the CLD cell energy at fixed total baryon density rho_B.

    Two-stage optimization (3 free variables: f, R_d, rho_I; the cluster proton
    density rho_Ip solves full cell beta equilibrium with interface terms; the
    gas is treated as pure neutron, rho_IIp = 0):
      1. Coarse grid scan over (f, R_d, rho_I) for a good starting point.
      2. Nelder-Mead simplex refinement.

    The gas density rho_II is set by baryon conservation. beta-equilibrium is
    enforced for the full cell energy, including interface terms; the pure-
    neutron gas approximation is the standard low-density inner-crust limit.
    """
    if not (math.isfinite(rho_B) and rho_B > 0.0 and math.isfinite(sigma_0)
            and sigma_0 >= 0.0 and math.isfinite(b) and b > 0.0):
        raise ValueError("CLD requires finite positive density/exponent and nonnegative surface coefficient")
    rho_sat = model.rho_sat
    f_grid = [0.005, 0.01, 0.03, 0.08, 0.15, 0.25, 0.4]
    Rd_grid = [4.0, 6.0, 9.0, 13.0, 20.0]
    rhoI_grid = [0.5 * rho_sat, 0.7 * rho_sat, 0.9 * rho_sat, 1.1 * rho_sat, 1.3 * rho_sat]

    best_E = float("inf"); best_x = None; best_info = None
    for f in f_grid:
        for R_d in Rd_grid:
            for rho_I in rhoI_grid:
                x = (f, math.log(R_d), rho_I)
                E, info = _cld_objective(model, rho_B, sigma_0, b, x)
                if E < best_E:
                    best_E = E; best_x = x; best_info = info
    if best_x is None:
        return CLDCellResult(rho_B, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, False, 0, "no feasible point")
    x_nm, E_nm, info_nm, n_iter_nm, optimizer_converged = _nelder_mead(model, rho_B, sigma_0, b, best_x)
    if E_nm < best_E:
        best_x, best_E, best_info = x_nm, E_nm, info_nm
    f, log_Rd, rho_I = best_x
    f = max(1e-4, min(0.95, f))
    R_d = math.exp(max(math.log(2.0), min(math.log(40.0), log_Rd)))
    rho_I = max(0.01, min(0.5, rho_I))
    # final beta-equilibrium solve
    rho_Ip, rho_e, mu_e, beta_resid = solve_beta_equilibrium_cluster(model, rho_I, 0.0, f,
        R_d=R_d, sigma_0=sigma_0, b=b)
    rho_IIp = 0.0
    if f < 1:
        rho_II = (rho_B - f * rho_I) / (1.0 - f)
    else:
        rho_II = 0.0
    rho_II = max(0.0, rho_II)
    E_cell, info = cld_cell_energy(model, rho_B, f, R_d, rho_I, rho_Ip, rho_II, rho_IIp,
                                   sigma_0, b)
    info["rho_Ip"] = rho_Ip; info["rho_IIp"] = rho_IIp
    info["rho_e"] = rho_e; info["mu_e"] = mu_e; info["beta_resid"] = beta_resid
    # cell mu_n/mu_p from the cluster phase (equilibrating)
    rI = solve_homogeneous(model, rho_I, rho_Ip)
    mu_n_cell = rI.mu_n; mu_p_cell = rI.mu_p
    # Pressure is a property of the minimized energy curve, not of one cell:
    # reconstruct_cld.thermodynamic_table derives it from eps(n_B).
    P_cell = float("nan")
    yp_cell = (f * rho_Ip + (1 - f) * rho_IIp) / rho_B if rho_B > 0 else 0
    return CLDCellResult(
        rho_B=rho_B, f=f, R_d=R_d, rho_I=rho_I, rho_Ip=rho_Ip,
        rho_II=rho_II, rho_IIp=rho_IIp, rho_e=rho_e, y_p_cell=yp_cell,
        E_cell=E_cell, E_binding_per_baryon=E_cell / rho_B if rho_B > 0 else 0,
        E_coul=info["E_coul"], E_surf=info["E_surf"], E_L=info["E_L"],
        P_cell=P_cell, mu_n=mu_n_cell, mu_p=mu_p_cell, mu_e=mu_e,
        mu_n_gas=solve_homogeneous(model, rho_II, rho_IIp).mu_n if rho_II > 0 else float("nan"),
        converged=(optimizer_converged and info["rI_converged"] and info["rII_converged"]
                   and beta_resid <= 1e-6), n_iter=n_iter_nm,
        message=f"proxy CLD; Nelder-Mead converged={optimizer_converged}; full cell beta resid={beta_resid:.3e}",
    )
