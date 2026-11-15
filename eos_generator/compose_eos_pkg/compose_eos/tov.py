"""Causality-aware barotropic interpolation and TOV integration.

All EOS arrays are in CompOSE natural units (MeV, fm^-3). TOV integration is
performed in SI units. The interpolant is a monotone PCHIP representation of
log(P) as a function of log(epsilon), so its derivative defines a single
barotropic sound speed consistently with the pressure used in the TOV RHS.

This module deliberately does not extrapolate below/above the tabulated EOS.
The stellar surface is the lowest tabulated pressure; the returned radius is
therefore a finite-table radius unless an outer-crust/vacuum extension is part
of the input table.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.integrate import cumulative_simpson

MEV_FM3_TO_PA = 1.602176634e32
G = 6.67430e-11
C = 2.99792458e8
MSUN = 1.98847e30


@dataclass(frozen=True)
class TOVConfig:
    rtol: float = 2.0e-8
    max_radius_m: float = 3.0e6
    max_step_m: float = 2.0e4
    initial_radius_m: float = 1.0
    causal_tolerance: float = 1.0e-8

    def __post_init__(self):
        if not (0.0 < self.rtol < 1.0):
            raise ValueError("rtol must lie in (0, 1)")
        if self.max_radius_m <= 0 or self.max_step_m <= 0 or self.initial_radius_m <= 0:
            raise ValueError("radii and max_step must be positive")
        if self.causal_tolerance < 0:
            raise ValueError("causal_tolerance must be non-negative")


@dataclass(frozen=True)
class CausalDomain:
    status: str
    nb_max: float
    eps_max: float
    pressure_max: float
    cs2_at_limit: float
    reason: str
    crossing_index: int | None = None
    cs2_certified_max: float | None = None

    def to_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class TOVStarResult:
    status: str
    central_nb: float
    mass_solar: float | None
    radius_km: float | None
    surface_event: bool
    surface_pressure_mev_fm3: float
    nfev: int
    message: str

    def to_dict(self):
        return asdict(self)


class BarotropicEOS:
    """Validated monotone PCHIP interpolation for a cold barotropic EOS."""

    def __init__(self, nb, pressure, energy_density):
        self.nb = np.asarray(nb, dtype=float)
        self.pressure = np.asarray(pressure, dtype=float)
        self.energy_density = np.asarray(energy_density, dtype=float)
        arrays = (self.nb, self.pressure, self.energy_density)
        if any(a.ndim != 1 for a in arrays) or len({a.size for a in arrays}) != 1:
            raise ValueError("nb, pressure and energy_density must be same-length 1D arrays")
        if self.nb.size < 4:
            raise ValueError("at least four EOS points are required for robust PCHIP interpolation")
        if any(np.any(~np.isfinite(a)) for a in arrays):
            raise ValueError("EOS arrays must contain only finite values")
        if any(np.any(a <= 0.0) for a in arrays):
            raise ValueError("logarithmic EOS interpolation requires positive values")
        if any(np.any(np.diff(a) <= 0.0) for a in arrays):
            raise ValueError("nb, pressure and energy_density must be strictly increasing")

        self._ln_nb = np.log(self.nb)
        self._ln_p = np.log(self.pressure)
        self._ln_eps = np.log(self.energy_density)
        self._lnp_of_lneps = PchipInterpolator(self._ln_eps, self._ln_p, extrapolate=False)
        self._lneps_of_lnp = PchipInterpolator(self._ln_p, self._ln_eps, extrapolate=False)
        self._lnp_of_lnnb = PchipInterpolator(self._ln_nb, self._ln_p, extrapolate=False)
        self._lneps_of_lnnb = PchipInterpolator(self._ln_nb, self._ln_eps, extrapolate=False)
        self._lnnb_of_lneps = PchipInterpolator(self._ln_eps, self._ln_nb, extrapolate=False)

    @staticmethod
    def _check_range(values, lo, hi, label):
        x = np.asarray(values, dtype=float)
        if np.any(~np.isfinite(x)) or np.any(x < lo) or np.any(x > hi):
            raise ValueError(f"{label} outside tabulated range [{lo:.12g}, {hi:.12g}]")

    def pressure_of_energy(self, eps):
        self._check_range(eps, self.energy_density[0], self.energy_density[-1], "energy density")
        return np.exp(self._lnp_of_lneps(np.log(eps)))

    def energy_of_pressure(self, pressure):
        self._check_range(pressure, self.pressure[0], self.pressure[-1], "pressure")
        return np.exp(self._lneps_of_lnp(np.log(pressure)))

    def pressure_of_nb(self, nb):
        self._check_range(nb, self.nb[0], self.nb[-1], "baryon density")
        return np.exp(self._lnp_of_lnnb(np.log(nb)))

    def energy_of_nb(self, nb):
        self._check_range(nb, self.nb[0], self.nb[-1], "baryon density")
        return np.exp(self._lneps_of_lnnb(np.log(nb)))

    def nb_of_energy(self, eps):
        self._check_range(eps, self.energy_density[0], self.energy_density[-1], "energy density")
        return np.exp(self._lnnb_of_lneps(np.log(eps)))

    def cs2(self, eps):
        """Return dP/d epsilon from the same log-log PCHIP barotrope."""
        self._check_range(eps, self.energy_density[0], self.energy_density[-1], "energy density")
        logeps = np.log(eps)
        p = np.exp(self._lnp_of_lneps(logeps))
        dlnp_dlneps = self._lnp_of_lneps.derivative()(logeps)
        return (p / np.asarray(eps, dtype=float)) * dlnp_dlneps

    def cs2_of_pressure(self, pressure):
        """dP/d epsilon of the interpolant the TOV RHS integrates, ln eps(ln P)."""
        self._check_range(pressure, self.pressure[0], self.pressure[-1], "pressure")
        x = np.log(np.asarray(pressure, dtype=float))
        z = self._lneps_of_lnp(x)
        return np.exp(x - z) / self._lneps_of_lnp.derivative()(x)

    def cs2_on_grid(self):
        return np.asarray(self.cs2(self.energy_density), dtype=float)


def _rhs_interval_polys(eos: BarotropicEOS, j: int):
    """Interval j of the RHS interpolant z(x)=ln eps(x), x=ln P, s in [0, 1].

    Returns (x0, h, z_s, dz_s) as numpy Polynomials in s, with dz/dx = dz_s/h.
    """
    from numpy.polynomial import Polynomial
    spline = eos._lneps_of_lnp
    x0 = float(spline.x[j]); h = float(spline.x[j + 1] - spline.x[j])
    c3, c2, c1, c0 = (float(v) for v in spline.c[:, j])
    z_s = Polynomial([c0, c1 * h, c2 * h * h, c3 * h ** 3])
    return x0, h, z_s, z_s.deriv()


def _poly_range(poly, a, b):
    """Exact (min, max) of a polynomial of degree <= 3 on [a, b]."""
    cand = [a, b]
    d = poly.deriv().trim(tol=0.0)
    if d.degree() >= 1:
        for r in d.roots():
            if abs(r.imag) <= 1e-12 * max(1.0, abs(r.real)) and a < r.real < b:
                cand.append(float(r.real))
    vals = [float(poly(t)) for t in cand]
    return min(vals), max(vals)


def _certify_rhs_causality(eos: BarotropicEOS, tolerance: float, min_width: float = 1e-12):
    """Ordered branch-and-bound certificate of 0 < c_s^2 <= 1 + tolerance.

    On the RHS interpolant c_s^2 = exp(x - z(x)) / z'(x).  On any sub-interval
    w = x - z is a cubic and z' a quadratic, so their exact ranges give rigorous
    bounds  exp(min w)/max z' <= c_s^2 <= exp(max w)/min z'.  Sub-intervals are
    processed left to right and bisected until certified or a violating point is
    found.  A leaf of width ``min_width`` (in the unit interval coordinate) that
    is neither certified nor violating at its endpoints is recorded as
    unresolved and the scan continues; if no genuine violation follows, the
    domain is cut fail-closed at the first unresolved leaf.

    Returns (violation_x, certified_max, x_safe, undecided) with x = ln P.
    """
    bound = 1.0 + tolerance
    certified_max = -np.inf
    unresolved = None
    nint = len(eos._lneps_of_lnp.x) - 1
    for j in range(nint):
        x0, h, z_s, dz_s = _rhs_interval_polys(eos, j)
        w = Polynomial_from_shift(x0, h) - z_s

        def cs2_at(s):
            dz = float(dz_s(s)) / h
            return np.exp(float(w(s))) / dz if dz > 0.0 else np.inf

        stack = [(0.0, 1.0)]
        while stack:
            a, b = stack.pop()
            wmin, wmax = _poly_range(w, a, b)
            dmin, dmax = _poly_range(dz_s, a, b)
            dmin /= h; dmax /= h
            upper = np.exp(wmax) / dmin if dmin > 0.0 else np.inf
            if upper <= bound:
                certified_max = max(certified_max, upper)
                continue
            if cs2_at(a) > bound:
                return x0 + h * a, certified_max, x0 + h * a, False
            if b - a <= min_width:
                if cs2_at(b) > bound:
                    return x0 + h * b, certified_max, x0 + h * a, False
                certified_max = max(certified_max, cs2_at(a), cs2_at(b))
                if unresolved is None:
                    unresolved = x0 + h * a
                continue
            mid = 0.5 * (a + b)
            stack.append((mid, b))   # right half processed after left half
            stack.append((a, mid))
    if unresolved is not None:
        return unresolved, certified_max, unresolved, True
    return None, certified_max, float(eos._lneps_of_lnp.x[-1]), False


def Polynomial_from_shift(x0, h):
    from numpy.polynomial import Polynomial
    return Polynomial([x0, h])


def find_causal_domain(eos: BarotropicEOS, tolerance: float = 1e-8,
                       samples_per_interval: int = 16) -> CausalDomain:
    r"""Return the connected low-density domain with 0 < c_s^2 <= 1 + tolerance.

    Certification is performed on the interpolant the TOV right-hand side
    actually integrates, ln eps(ln P), with exact interval bounds (no sampling):
    see :func:`_certify_rhs_causality`.  The first violation is bracketed by a
    certified point and root-solved for c_s^2 = 1 (or 1 + tolerance when the
    certified side already lies in (1, 1 + tolerance]).  ``samples_per_interval``
    is accepted for configuration compatibility and is not used for any gate.
    """
    if tolerance < 0:
        raise ValueError("tolerance must be non-negative")
    if samples_per_interval < 2:
        raise ValueError("samples_per_interval must be at least 2")
    violation, cert_max, x_safe, undecided = _certify_rhs_causality(eos, tolerance)
    if violation is None:
        return CausalDomain(
            status="FULL_TABLE_CAUSAL", nb_max=float(eos.nb[-1]),
            eps_max=float(eos.energy_density[-1]), pressure_max=float(eos.pressure[-1]),
            cs2_at_limit=float(eos.cs2_of_pressure(eos.pressure[-1])),
            reason=("Certified 0 < c_s^2 <= 1 + tol on the TOV RHS interpolant "
                    "ln eps(ln P) by exact interval bounds; no sampling."),
            cs2_certified_max=float(cert_max))
    x_first = float(eos._lneps_of_lnp.x[0])
    if violation <= x_first:
        raise ValueError("EOS is mechanically unstable or acausal at the first tabulated point")
    cs2 = lambda x: float(eos.cs2_of_pressure(np.exp(x)))
    target = 1.0
    if undecided:
        x_cross = x_safe
    else:
        # Cut at the c_s^2 = 1 crossing connected to the violation: step back
        # from the last certified point until c_s^2 <= 1, then root-solve.
        x_hi, x_lo, step = violation, x_safe, 1e-12 * max(1.0, abs(x_safe))
        while cs2(x_lo) > 1.0 and x_lo > x_first:
            x_hi = x_lo
            x_lo = max(x_first, x_lo - step)
            step *= 2.0
        if cs2(x_lo) > 1.0:
            target = 1.0 + tolerance
            x_lo = x_first
        x_cross = float(brentq(lambda x: cs2(x) - target, x_lo, x_hi,
                               xtol=1e-15 * max(1.0, abs(x_hi)), rtol=1e-15))
    p_cross = float(np.exp(x_cross))
    eps_cross = float(eos.energy_of_pressure(p_cross))
    nb_cross = float(eos.nb_of_energy(min(max(eps_cross, eos.energy_density[0]), eos.energy_density[-1])))
    source_interval = max(0, min(len(eos.pressure) - 2,
                                 int(np.searchsorted(eos.pressure, p_cross, side="right") - 1)))
    reason = ("First c_s^2 crossing of the TOV RHS interpolant ln eps(ln P), located by "
              "exact interval bounds; higher-density branch excluded.")
    if undecided:
        reason += " Bound undecided at minimum width; cut fail-closed at the last certified point."
    elif target != 1.0:
        reason += " No c_s^2 <= 1 point precedes the violation; root solved at c_s^2 = 1 + tol."
    return CausalDomain(
        status="ACASUAL_CUTOFF", nb_max=nb_cross, eps_max=eps_cross,
        pressure_max=p_cross, cs2_at_limit=cs2(x_cross), reason=reason,
        crossing_index=source_interval, cs2_certified_max=float(cert_max))


def solve_tov_star(eos: BarotropicEOS, central_nb: float,
                   config: TOVConfig | None = None,
                   causal_domain: CausalDomain | None = None) -> TOVStarResult:
    """Integrate one TOV star; stop at the finite-table surface pressure."""
    cfg = config or TOVConfig()
    if causal_domain is not None and central_nb > causal_domain.nb_max * (1.0 + 2e-12):
        raise ValueError(f"central density {central_nb:g} lies outside causal domain ending at {causal_domain.nb_max:g}")
    pc_mev = float(eos.pressure_of_nb(central_nb))
    ec_mev = float(eos.energy_of_nb(central_nb))
    pc = pc_mev * MEV_FM3_TO_PA
    surface_p_mev = float(eos.pressure[0])
    surface_p = max(1.0, surface_p_mev * MEV_FM3_TO_PA)
    if pc <= surface_p:
        raise ValueError("central pressure must exceed the tabulated surface pressure")

    ec = ec_mev * MEV_FM3_TO_PA / C**2
    # Near-origin expansion. Choose r0 so the initial pressure decrement is
    # small compared with Pc, rather than relying on a fixed central seed.
    coeff = (2.0 * np.pi * G / 3.0) * (ec + pc / C**2) * (ec + 3.0 * pc / C**2)
    if not np.isfinite(coeff) or coeff <= 0.0:
        raise ValueError("invalid central TOV expansion coefficient")
    r0 = min(cfg.initial_radius_m, max(1e-4, np.sqrt(1e-8 * pc / coeff)))
    p0 = pc - coeff * r0**2
    if p0 <= surface_p:
        raise ValueError("central expansion step is too large relative to central pressure")
    m0 = (4.0 * np.pi / 3.0) * ec * r0**3

    def rhs(r, state):
        mass, pressure_si = state
        # Clamp intermediate Runge-Kutta stages to the table surface. This
        # keeps the RHS continuous through the terminal event and avoids the
        # discontinuous zero-RHS branch that caused step-size failures.
        p_eval = max(float(pressure_si), surface_p)
        eps_mev = float(eos.energy_of_pressure(p_eval / MEV_FM3_TO_PA))
        rho = eps_mev * MEV_FM3_TO_PA / C**2
        compactness_factor = 1.0 - 2.0 * G * mass / (r * C**2)
        if compactness_factor <= 0.0:
            raise FloatingPointError("TOV integration reached a trapped-surface condition")
        dm_dr = 4.0 * np.pi * r**2 * rho
        dp_dr = (-G * (rho + p_eval / C**2)
                 * (mass + 4.0 * np.pi * r**3 * p_eval / C**2)
                 / (r**2 * compactness_factor))
        return [dm_dr, dp_dr]

    def surface_event(r, state):
        return state[1] - surface_p
    surface_event.terminal = True
    surface_event.direction = -1

    try:
        sol = solve_ivp(
            rhs, (r0, cfg.max_radius_m), [m0, p0], method="RK45",
            rtol=cfg.rtol, atol=np.array([1e-10 * MSUN, 1e-12 * pc]),
            max_step=cfg.max_step_m, events=surface_event)
    except Exception as exc:
        return TOVStarResult("INTEGRATION_ERROR", float(central_nb), None, None,
                              False, surface_p_mev, 0, f"{type(exc).__name__}: {exc}")
    event_hit = bool(sol.t_events and len(sol.t_events[0]) == 1)
    if not sol.success or not event_hit:
        status = "NO_SURFACE_EVENT" if sol.success else "SOLVER_FAILURE"
        return TOVStarResult(status, float(central_nb), None, None, False,
                              surface_p_mev, int(sol.nfev), str(sol.message))
    r_surface = float(sol.t_events[0][0])
    m_surface = float(sol.y_events[0][0][0])
    if not (np.isfinite(r_surface) and np.isfinite(m_surface) and r_surface > 0 and m_surface > 0):
        return TOVStarResult("NONFINITE_OUTPUT", float(central_nb), None, None,
                              True, surface_p_mev, int(sol.nfev), "non-finite surface state")
    return TOVStarResult("SUCCESS", float(central_nb), m_surface / MSUN,
                         r_surface / 1000.0, True, surface_p_mev,
                         int(sol.nfev), str(sol.message))



def first_law_nb_consistency(eos: BarotropicEOS, oversampling: int = 12) -> dict:
    r"""Compare tabulated n_b(epsilon) with the zero-T first-law integral.

    For a cold one-component barotrope, d epsilon = (epsilon + P) d n_b/n_b,
    hence d ln n_b/d ln epsilon = epsilon/(epsilon+P).  The integral is
    anchored at the lowest-density table point; its shape residual is therefore
    independent of the arbitrary normalization of baryon number density.
    """
    if oversampling < 2:
        raise ValueError("oversampling must be at least 2")
    n_grid = max(256, int(eos.nb.size * oversampling))
    logeps = np.linspace(float(np.log(eos.energy_density[0])),
                         float(np.log(eos.energy_density[-1])), n_grid)
    eps = np.exp(logeps)
    # exp(log(endpoint)) can exceed the exact endpoint by one ulp. Restore the
    # tabulated endpoints before range-checked interpolation.
    eps[0] = eos.energy_density[0]
    eps[-1] = eos.energy_density[-1]
    pressure = np.asarray(eos.pressure_of_energy(eps), dtype=float)
    integrand = eps / (eps + pressure)
    integral = cumulative_simpson(integrand, x=logeps, initial=0.0)
    lognb_pred = np.log(eos.nb[0]) + integral
    lognb_actual = np.log(np.asarray(eos.nb_of_energy(eps), dtype=float))
    residual = lognb_actual - lognb_pred
    # Exclude only the anchored point; report log residual to preserve scale
    # symmetry and a direct multiplicative interpretation.
    r = residual[1:]
    return {
        "status": "CHECKED",
        "oversampling": int(oversampling),
        "integration_points": int(n_grid),
        "max_abs_log_ratio": float(np.max(np.abs(r))),
        "rms_log_ratio": float(np.sqrt(np.mean(r * r))),
        "max_fractional_mismatch": float(np.expm1(np.max(np.abs(r)))),
        "formula": "d ln(n_b)/d ln(epsilon) = epsilon/(epsilon + P)",
    }


def stable_branch_summary(results: Iterable[TOVStarResult]) -> dict:
    """Summarize the sampled stable branch up to the first mass maximum."""
    results = list(results)
    good = [r for r in results if r.status == "SUCCESS"]
    if len(good) < 3:
        return {"status": "INSUFFICIENT_SUCCESSFUL_STARS", "n_success": len(good)}
    nc = np.asarray([r.central_nb for r in good])
    mass = np.asarray([r.mass_solar for r in good])
    radius = np.asarray([r.radius_km for r in good])
    order = np.argsort(nc)
    nc, mass, radius = nc[order], mass[order], radius[order]
    imax = int(np.argmax(mass))
    branch_mass = mass[:imax + 1]
    branch_radius = radius[:imax + 1]
    branch_nc = nc[:imax + 1]
    r14 = None
    crosses = np.flatnonzero((branch_mass[:-1] - 1.4) * (branch_mass[1:] - 1.4) <= 0.0)
    if crosses.size:
        j = int(crosses[0])
        if branch_mass[j + 1] != branch_mass[j]:
            frac = (1.4 - branch_mass[j]) / (branch_mass[j + 1] - branch_mass[j])
            r14 = float(branch_radius[j] + frac * (branch_radius[j + 1] - branch_radius[j]))
    upper_limited = bool(imax == len(mass) - 1)
    return {
        "status": "SAMPLED_SEQUENCE",
        "n_success": len(good),
        "n_fail": sum(r.status != "SUCCESS" for r in results),
        "Mmax_sampled_solar": float(mass[imax]),
        "Mmax_peak_classification": ("UPPER_DOMAIN_LIMITED" if upper_limited else "INTERIOR_PEAK_RESOLVED"),
        "Mmax_is_global_maximum_claimed": False,
        "R_at_Mmax_km": float(radius[imax]),
        "nc_at_Mmax_fm3": float(nc[imax]),
        "R14_km_linear_on_stable_branch": r14,
        "R14_in_sampled_stable_branch": r14 is not None,
    }


def solve_tov_sequence(eos: BarotropicEOS, n_points: int = 24,
                       central_nb_min: float = 0.08,
                       config: TOVConfig | None = None,
                       causal_domain: CausalDomain | None = None) -> dict:
    """Run a log-spaced central-density sequence up to the causal cutoff."""
    if n_points < 4:
        raise ValueError("n_points must be at least 4")
    domain = causal_domain or find_causal_domain(eos, (config or TOVConfig()).causal_tolerance)
    lo = max(float(eos.nb[0]) * (1.0 + 1e-10), float(central_nb_min))
    hi = float(domain.nb_max) * (1.0 - 1e-9)
    if not (lo < hi):
        raise ValueError(f"central-density interval is empty: [{lo:g}, {hi:g}]")
    centers = np.geomspace(lo, hi, n_points)
    stars = [solve_tov_star(eos, float(nc), config=config, causal_domain=domain)
             for nc in centers]
    summary = stable_branch_summary(stars)
    return {
        "causal_domain": domain.to_dict(),
        "config": asdict(config or TOVConfig()),
        "central_nb_grid_fm3": centers.tolist(),
        "stars": [s.to_dict() for s in stars],
        "summary": summary,
        "surface_definition": "lowest tabulated pressure; radius is finite-table radius unless the table includes an outer-crust/vacuum extension",
    }


__all__ = [
    "BarotropicEOS", "TOVConfig", "CausalDomain", "TOVStarResult",
    "find_causal_domain", "first_law_nb_consistency", "solve_tov_star", "solve_tov_sequence",
    "stable_branch_summary",
]
