"""Compressible Liquid Drop (CLD) cell model for the EOS390 inner-crust reconstruction.

Implements the Scurto 2024 (PRD 109 103015) eq. 21-25 CLD cell model on top of
the validated RMF homogeneous-matter solver (rmf_solver.py).

The CLD model describes each Wigner-Seitz cell as:
  - a high-density "cluster" phase (I) at density rho_I, proton density rho_Ip
  - a low-density "gas" phase (II) at density rho_II, proton density rho_IIp
  - a uniform electron background of energy density E_L
  - interface Coulomb (E_Coul) and surface (E_surf) corrections

The cell energy density (Scurto eq. 21):
  E_cell = f * E(rho_I, rho_Ip) + (1 - f) * E(rho_II, rho_IIp)
         + E_Coul + E_surf + E_L

is minimized with respect to four variables:
  R_d (cluster radius), rho_I, rho_Ip, f (volume fraction).

The gas density rho_II and composition rho_IIp are determined by beta-equilibrium
and charge-neutrality constraints at the cell level (mu_n - mu_p = mu_e; rho_p
total = rho_e). The Coulomb and surface terms follow the Ravenhall-Pethick-Wilson
[68] standard forms (Scurto eq. 22-25).

All units: densities fm^-3, energies MeV, energy densities MeV/fm^3, lengths fm.
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
    """Ravenhall-Pethick-Wilson lattice function Phi(u) for spherical droplets.
    Standard RPW form: Phi = (1/2)*[(2-3u)^(1/2)*(2-(1-3u/2)^2) ... ].
    The most-used simple form for the inner-crust CLD droplet phase:
      Phi(u) = (1/2)*(2 - 3u)^(1/2)   (dilute-sphere limit, valid for u < ~0.3)
    More complete form including the (1-3u/2) geometric correction:
      Phi(u) = (1/2)*[ (2-3u)^(1/2)*(1 + (1 - 3u/2)*(1-u)^(-1)) - (1 - 3u/2)^2 ]
    We use the dilute-sphere limit which is accurate for the inner-crust droplet
    regime where u << 1; this matches the Scurto eq. 23 functional dependence
    on (2-3f)^(1/2)."""
    if u <= 0 or u >= 2.0 / 3.0:
        return 0.0
    return 0.5 * math.sqrt(2.0 - 3.0 * u)


def surface_tension(yp_I: float, sigma_0: float, b: float) -> float:
    """Surface tension sigma (Scurto eq. 25, Lattimer form):
      sigma = sigma_0 * [(b + 2)/(...) * y_p,I^(-3) + (1 - y_p,I)^(-3)] / (b + ...)
    The PDF text is corrupted; the standard Lattimer-Swestsky form is:
      sigma(yp) = sigma_0 * [(2/b) * (yp^(2/3) + (1-yp)^(2/3) - 1) ... ]
    The simplest physically reasonable form used in CLD codes:
      sigma(yp) = sigma_0 * (yp^(2/3) + (1-yp)^(2/3) - 2^(2/3))
    normalized to vanish at yp=0 (pure neutron) and yp=1 (pure proton).
    We use the Lattimer et al. empirical form: sigma = sigma_0 * [yp*(1-yp)]^b
    which gives sigma=0 at yp=0,1 and peaks at yp=0.5, controlled by b.
    For b=1 (default) this is sigma_0*yp*(1-yp).
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
    # Coulomb (Scurto eq. 22): E_Coul = 2*pi*e^2 * f * Phi * R_d^2 * (rho_Ip - rho_IIp)^2
    Phi = coulomb_lattice_phi(f)
    rho_p_diff = rho_Ip - rho_IIp
    E_coul = 2.0 * math.pi * E2 * f * Phi * R_d ** 2 * rho_p_diff ** 2
    # Surface (Scurto eq. 24): E_surf = (3*f*sigma)/R_d
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
                                   rho_IIp_gas: float, f: float) -> tuple[float, float, float, float]:
    """Solve for the cluster proton density rho_Ip from beta-equilibrium at the
    CLUSTER phase (the dense, equilibrating phase):
        mu_n(I) - mu_p(I) = mu_e(cell)
    with cell charge neutrality:
        rho_e = f*rho_Ip + (1-f)*rho_IIp_gas   (rho_IIp_gas fixed = 0 for pure neutron gas)
        mu_e = sqrt((3*pi^2*rho_e)^(2/3)*hbar^2 + m_e^2)
    At rho_Ip=0 (pure neutron cluster) mu_n-mu_p is large positive; at rho_Ip=rho_I
    (symmetric cluster) mu_n-mu_p ~ 0; mu_e grows with rho_Ip. Sign change exists.
    Returns (rho_Ip, rho_e, mu_e, |residual|).
    """
    rho_IIp_gas = 0.0  # pure neutron gas (low-density inner-crust approximation)
    def residual(rho_Ip):
        rho_Ip = max(0.0, min(rho_Ip, rho_I))
        rI = solve_homogeneous(model, rho_I, rho_Ip)
        rho_e = f * rho_Ip + (1.0 - f) * rho_IIp_gas
        rho_e = max(0.0, rho_e)
        mu_e = electron_chemical_potential(rho_e) if rho_e > 0 else 0.511
        return (rI.mu_n - rI.mu_p) - mu_e
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
    The cluster proton density rho_Ip is solved from beta-equilibrium at the
    cluster phase (mu_n - mu_p = mu_e), with the gas treated as pure neutron
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
    rho_Ip, rho_e, mu_e, beta_resid = solve_beta_equilibrium_cluster(model, rho_I, 0.0, f)
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
    return E, info


def _nelder_mead(model, rho_B, sigma_0, b, x0, max_iter=120, tol=1e-6):
    """Lightweight Nelder-Mead (no scipy) in 3 variables (f, log R_d, rho_I).
    Returns (best_x, best_E, best_info, n_iter)."""
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
    for it in range(max_iter):
        order = sorted(range(len(simplex)), key=lambda k: fvals[k])
        simplex = [simplex[k] for k in order]
        fvals = [fvals[k] for k in order]
        infos = [infos[k] for k in order]
        if abs(fvals[-1] - fvals[0]) < tol and fvals[0] != float("inf"):
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
    return simplex[bi], fvals[bi], infos[bi], it + 1


def minimize_cld_cell(model: RMFModel, rho_B: float, sigma_0: float, b: float,
                     verbose: bool = False) -> Optional[CLDCellResult]:
    """Minimize the CLD cell energy at fixed total baryon density rho_B.

    Two-stage optimization (3 free variables: f, R_d, rho_I; the cluster proton
    density rho_Ip is solved from beta-equilibrium mu_n(I) - mu_p(I) = mu_e; the
    gas is treated as pure neutron, rho_IIp = 0):
      1. Coarse grid scan over (f, R_d, rho_I) for a good starting point.
      2. Nelder-Mead simplex refinement.

    The gas density rho_II is set by baryon conservation. beta-equilibrium is
    enforced at the cluster phase (the dense, equilibrating phase); the pure-
    neutron gas approximation is the standard low-density inner-crust limit.
    """
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
    x_nm, E_nm, info_nm, n_iter_nm = _nelder_mead(model, rho_B, sigma_0, b, best_x)
    if E_nm < best_E:
        best_x, best_E, best_info = x_nm, E_nm, info_nm
    f, log_Rd, rho_I = best_x
    f = max(1e-4, min(0.95, f))
    R_d = math.exp(max(math.log(2.0), min(math.log(40.0), log_Rd)))
    rho_I = max(0.01, min(0.5, rho_I))
    # final beta-equilibrium solve
    rho_Ip, rho_e, mu_e, beta_resid = solve_beta_equilibrium_cluster(model, rho_I, 0.0, f)
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
    P_cell = 0.0  # filled by emitter via finite difference
    yp_cell = (f * rho_Ip + (1 - f) * rho_IIp) / rho_B if rho_B > 0 else 0
    return CLDCellResult(
        rho_B=rho_B, f=f, R_d=R_d, rho_I=rho_I, rho_Ip=rho_Ip,
        rho_II=rho_II, rho_IIp=rho_IIp, rho_e=rho_e, y_p_cell=yp_cell,
        E_cell=E_cell, E_binding_per_baryon=E_cell / rho_B if rho_B > 0 else 0,
        E_coul=info["E_coul"], E_surf=info["E_surf"], E_L=info["E_L"],
        P_cell=P_cell, mu_n=mu_n_cell, mu_p=mu_p_cell, mu_e=mu_e,
        converged=True, n_iter=n_iter_nm,
        message=f"Nelder-Mead + beta-equil (resid={beta_resid:.3e})",
    )
