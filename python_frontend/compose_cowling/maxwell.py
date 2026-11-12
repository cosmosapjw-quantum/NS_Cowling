"""Phase-aware TOV and Cowling integration for a sharp Maxwell transition.

The existing native Maxwell builder supplies a two-branch, phase-labelled EOS.
This module is an independent SciPy reference/integration path: it never sends a
phase-discontinuous table through the smooth single-valued Fortran table reader.

Units at the public EOS boundary are CompOSE-native (MeV fm^-3, fm^-3, MeV);
stellar structure and perturbations are integrated in SI units.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from compose_eos import EOSData
from .direct_barotrope import DirectBarotrope, EOSState

MEV_FM3_TO_J_M3 = 1.602176634e32
MEV_TO_J = 1.602176634e-13
FM3_TO_M3 = 1.0e45
C_SI = 299792458.0
G_SI = 6.67430e-11
M_SUN_NOMINAL_KG = 1.3271244e20 / G_SI

Conversion = Literal["slow", "rapid"]


@dataclass(frozen=True)
class MaxwellInterface:
    pressure_pa: float
    chemical_potential_j: float
    energy_inner_j_m3: float
    energy_outer_j_m3: float
    density_inner_m3: float
    density_outer_m3: float
    radius_m: float

    @property
    def energy_jump_j_m3(self) -> float:
        """Positive high-density-core minus low-density-envelope jump."""
        return self.energy_inner_j_m3 - self.energy_outer_j_m3

    @property
    def density_jump_m3(self) -> float:
        """Positive high-density-core minus low-density-envelope jump."""
        return self.density_inner_m3 - self.density_outer_m3


class MaxwellEOS:
    """Two single-phase providers separated by one exact Maxwell jump."""

    def __init__(self, transition_table: EOSData):
        if not isinstance(transition_table, EOSData):
            raise TypeError("transition_table must be an EOSData instance")
        tags = np.asarray(transition_table.Iphase, dtype=int)
        if tags.ndim != 1 or tags.size != transition_table.n:
            raise ValueError("invalid phase labels")
        changes = np.flatnonzero(np.diff(tags) != 0)
        if changes.size != 1 or set(np.unique(tags)) != {1, 2}:
            raise ValueError("exactly one ordered low(1)-to-high(2) Maxwell interface is required")
        cut = int(changes[0])
        if tags[cut] != 1 or tags[cut + 1] != 2:
            raise ValueError("phase labels must be ordered low-density then high-density")

        def subtable(mask: np.ndarray) -> DirectBarotrope:
            return DirectBarotrope(
                np.asarray(transition_table.p[mask], dtype=float),
                np.asarray(transition_table.eps[mask], dtype=float),
                np.asarray(transition_table.nb[mask], dtype=float),
                np.asarray(transition_table.mub[mask], dtype=float),
            )

        self.low = subtable(tags == 1)
        self.high = subtable(tags == 2)
        p_lo = float(self.low.pressure[-1])
        p_hi = float(self.high.pressure[0])
        scale = max(abs(p_lo), abs(p_hi), 1.0)
        if abs(p_lo - p_hi) > 5e-10 * scale:
            raise ValueError("phase providers do not share one transition pressure")
        self.transition_pressure_mev_fm3 = 0.5 * (p_lo + p_hi)
        self.transition_pressure_pa = self.transition_pressure_mev_fm3 * MEV_FM3_TO_J_M3
        low_t = self.low.evaluate(p_lo)
        high_t = self.high.evaluate(p_hi)
        mu_scale = max(abs(low_t.chemical_potential), abs(high_t.chemical_potential), 1.0)
        if abs(low_t.chemical_potential - high_t.chemical_potential) > 2e-8 * mu_scale:
            raise ValueError("phase providers violate chemical-potential continuity at the interface")
        if high_t.number_density <= low_t.number_density or high_t.energy <= low_t.energy:
            raise ValueError("Maxwell interface must have a positive density and energy jump")
        closure = (high_t.energy - low_t.energy) - low_t.chemical_potential * (
            high_t.number_density - low_t.number_density
        )
        if abs(closure) > 2e-8 * max(abs(high_t.energy - low_t.energy), 1.0):
            raise ValueError("Maxwell latent-energy identity failed after phase-provider construction")
        self.mu_transition_mev = 0.5 * (low_t.chemical_potential + high_t.chemical_potential)
        delta_eps = high_t.energy - low_t.energy
        delta_n = high_t.number_density - low_t.number_density
        seidov_ratio = delta_eps / low_t.energy
        seidov_critical = 0.5 + 1.5 * self.transition_pressure_mev_fm3 / low_t.energy
        self.interface_thermo = {
            "pressure_pa": self.transition_pressure_pa,
            "pressure_mev_fm3": self.transition_pressure_mev_fm3,
            "chemical_potential_mev": self.mu_transition_mev,
            "energy_low_mev_fm3": low_t.energy,
            "energy_high_mev_fm3": high_t.energy,
            "density_low_fm3": low_t.number_density,
            "density_high_fm3": high_t.number_density,
            "latent_energy_mev_fm3": delta_eps,
            "density_jump_fm3": delta_n,
            "maxwell_identity_relative_residual": abs(closure) / max(abs(delta_eps), 1.0),
            # Seidov is an onset-stability diagnostic for rapid phase conversion,
            # not a substitute for constructing the full stable TOV sequence.
            "seidov_ratio": seidov_ratio,
            "seidov_critical_ratio": seidov_critical,
            "seidov_onset_stable_rapid_diagnostic": seidov_ratio < seidov_critical,
            "provider_first_law_residual_low": self.low.raw_first_law_max_relative_residual,
            "provider_first_law_residual_high": self.high.raw_first_law_max_relative_residual,
        }

    def evaluate(self, pressure_pa: float, side: Literal["low", "high"]) -> dict[str, float]:
        p = float(pressure_pa)
        if not np.isfinite(p) or p <= 0:
            raise ValueError("pressure must be finite and positive")
        if side == "low":
            provider = self.low
        elif side == "high":
            provider = self.high
        else:
            raise ValueError("side must be 'low' or 'high'")
        s: EOSState = provider.evaluate(p / MEV_FM3_TO_J_M3)
        return {
            "pressure_pa": p,
            "energy_j_m3": s.energy * MEV_FM3_TO_J_M3,
            "number_density_m3": s.number_density * FM3_TO_M3,
            "chemical_potential_j": s.chemical_potential * MEV_TO_J,
            "chi": s.sound_speed_sq,
        }

    @property
    def pressure_cutoff_pa(self) -> float:
        return float(self.low.pressure[0] * MEV_FM3_TO_J_M3)

    @property
    def high_pressure_max_pa(self) -> float:
        return float(self.high.pressure[-1] * MEV_FM3_TO_J_M3)


@dataclass(frozen=True)
class TOVSegment:
    side: Literal["high", "low"]
    solution: object  # dense solve_ivp solution in q=ln(P_ref/P), state [r,m,nu]
    r_start_m: float
    r_end_m: float
    pressure_reference_pa: float
    q_of_r: object


@dataclass(frozen=True)
class MaxwellStar:
    eos: MaxwellEOS
    inner: TOVSegment
    outer: TOVSegment
    r0_m: float
    interface: MaxwellInterface
    radius_m: float
    mass_kg: float
    mass_solar: float
    p_c_pa: float
    p_cut_pa: float
    nu_shift: float
    rtol: float
    accepted_steps: int

    def background(self, r_m: float, side: Literal["inner", "outer"] | None = None) -> dict[str, float]:
        r = float(r_m)
        if side is None:
            if abs(r - self.interface.radius_m) <= 1e-10 * max(self.interface.radius_m, 1.0):
                raise ValueError("side='inner' or 'outer' is required exactly at the Maxwell interface")
            side = "inner" if r < self.interface.radius_m else "outer"
        if side == "inner":
            seg, eos_side = self.inner, "high"
        elif side == "outer":
            seg, eos_side = self.outer, "low"
        else:
            raise ValueError("side must be 'inner' or 'outer'")
        tol = 1e-10 * max(abs(seg.r_start_m), abs(seg.r_end_m), 1.0)
        if r < seg.r_start_m - tol or r > seg.r_end_m + tol:
            raise ValueError("radius lies outside requested background segment")
        r = min(max(r, seg.r_start_m), seg.r_end_m)
        q = float(seg.q_of_r(r))
        q_lo, q_hi = sorted((float(seg.solution.t[0]), float(seg.solution.t[-1])))
        q = float(np.clip(q, q_lo, q_hi))
        r_sol, m_kg, nu_raw = np.asarray(seg.solution.sol(q), dtype=float)
        pressure_pa = seg.pressure_reference_pa * np.exp(-q)
        provider = self.eos.high if eos_side == "high" else self.eos.low
        p_lo = float(provider.pressure[0] * MEV_FM3_TO_J_M3)
        p_hi = float(provider.pressure[-1] * MEV_FM3_TO_J_M3)
        scale = max(abs(p_lo), abs(p_hi), 1.0)
        if pressure_pa < p_lo and p_lo - pressure_pa <= 2e-12 * scale:
            pressure_pa = p_lo
        elif pressure_pa > p_hi and pressure_pa - p_hi <= 2e-12 * scale:
            pressure_pa = p_hi
        therm = self.eos.evaluate(float(pressure_pa), eos_side)
        compact = 2.0 * G_SI * m_kg / (r * C_SI**2)
        if compact <= 0.0 or compact >= 1.0:
            raise ArithmeticError("invalid Schwarzschild compactness in background query")
        lam = -0.5 * np.log1p(-compact)
        nu = float(nu_raw) + self.nu_shift
        nu_prime = (G_SI / C_SI**2) * (
            m_kg / r**2 + 4.0 * np.pi * r * pressure_pa / C_SI**2
        ) / (1.0 - compact)
        return {
            "r_m": r, "m_kg": float(m_kg), **therm,
            "nu": nu, "lambda": float(lam), "dnu_dr": float(nu_prime),
        }


def solve_maxwell_tov(
    eos: MaxwellEOS,
    central_pressure_pa: float,
    *,
    pressure_cutoff_pa: float | None = None,
    r0_m: float = 10.0,
    r_limit_m: float = 1.0e6,
    rtol: float = 2.0e-10,
    atol_mass_kg: float = 1.0e10,
    atol_nu: float = 1.0e-13,
    max_step_log_pressure: float = 0.02,
) -> MaxwellStar:
    """Integrate a hybrid TOV model in two phases with an exact Pt restart."""
    if not isinstance(eos, MaxwellEOS):
        raise TypeError("eos must be a MaxwellEOS built from native Maxwell output")
    pc = float(central_pressure_pa)
    if not np.isfinite(pc) or pc <= eos.transition_pressure_pa or pc >= eos.high_pressure_max_pa:
        raise ValueError("central pressure must be above Pt and inside the high-phase EOS domain")
    pcut = eos.pressure_cutoff_pa if pressure_cutoff_pa is None else float(pressure_cutoff_pa)
    if not np.isfinite(pcut) or pcut <= 0 or pcut >= eos.transition_pressure_pa:
        raise ValueError("surface pressure cutoff must satisfy 0 < Pcut < Pt")
    if not (np.isfinite(r0_m) and r0_m > 0 and np.isfinite(r_limit_m) and r_limit_m > r0_m):
        raise ValueError("invalid radial integration limits")
    if not (0 < rtol < 1 and np.isfinite(max_step_log_pressure) and max_step_log_pressure > 0):
        raise ValueError("invalid ODE tolerances")

    center = eos.evaluate(pc, "high")
    eps_c = center["energy_j_m3"]
    alpha = (2.0 * np.pi * G_SI / 3.0) * ((eps_c + pc) / C_SI**2) * ((eps_c + 3.0 * pc) / C_SI**2)
    if alpha <= 0 or not np.isfinite(alpha):
        raise ArithmeticError("invalid TOV center-series coefficient")
    r0 = float(r0_m)
    # The leading center expansion must be a small correction to Pc.
    r0 = min(r0, np.sqrt(5.0e-3 * pc / alpha))
    drop = alpha * r0**2
    p0 = pc - drop
    if not (p0 > eos.transition_pressure_pa and p0 < pc):
        raise ValueError("center start is too close to Pt for the requested central pressure")
    m0 = (4.0 * np.pi / 3.0) * (eps_c / C_SI**2) * r0**3
    # Integrate in q=ln(P_ref/P), not r.  The transition and surface are exact
    # coordinate endpoints, so explicit RK stages never evaluate either phase
    # beyond its own pressure domain.
    y0 = np.array([r0, m0, 0.0], dtype=float)
    atol = np.array([1.0e-7, atol_mass_kg, atol_nu], dtype=float)

    def q_rhs(side: Literal["high", "low"], p_ref: float):
        def fun(q: float, y: np.ndarray) -> np.ndarray:
            r, m, _nu = y
            p = p_ref * np.exp(-q)
            if r <= 0 or p <= 0:
                raise ValueError("TOV state left the positive-r, positive-P domain")
            s = eos.evaluate(float(p), side)
            compact = 2.0 * G_SI * m / (r * C_SI**2)
            if compact >= 1.0:
                raise ArithmeticError("TOV integration encountered a trapped-surface denominator")
            phi = (G_SI / C_SI**2) * (
                m / r**2 + 4.0 * np.pi * r * p / C_SI**2
            ) / (1.0 - compact)
            dmdr = 4.0 * np.pi * r**2 * s["energy_j_m3"] / C_SI**2
            dpdr = -(s["energy_j_m3"] + p) * phi
            if dpdr >= 0 or not np.isfinite(dpdr):
                raise ArithmeticError("TOV pressure gradient is not negative")
            drdq = -p / dpdr
            return np.array([drdq, dmdr * drdq, phi * drdq])
        return fun

    def radius_limit_event(_q: float, y: np.ndarray) -> float:
        return r_limit_m - float(y[0])
    radius_limit_event.terminal = True
    radius_limit_event.direction = -1

    q_transition = float(np.log(p0 / eos.transition_pressure_pa))
    inner_sol = solve_ivp(
        q_rhs("high", p0), (0.0, q_transition), y0, method="DOP853",
        rtol=rtol, atol=atol, max_step=max_step_log_pressure,
        events=radius_limit_event, dense_output=True,
    )
    if not inner_sol.success:
        raise RuntimeError(f"inner TOV segment failed before Pt: {inner_sol.message}")
    if inner_sol.t[-1] < q_transition * (1.0 - 1e-12):
        raise RuntimeError("inner TOV segment reached r_limit_m before the phase interface")
    rt, mt, nut = np.asarray(inner_sol.y[:, -1], dtype=float)
    if not (np.isfinite(rt) and rt > r0 and np.isfinite(mt) and mt > 0):
        raise ArithmeticError("invalid TOV state at Maxwell interface")

    q_surface = float(np.log(eos.transition_pressure_pa / pcut))
    interface_y = np.array([rt, mt, nut], dtype=float)
    outer_sol = solve_ivp(
        q_rhs("low", eos.transition_pressure_pa), (0.0, q_surface), interface_y,
        method="DOP853", rtol=rtol, atol=atol, max_step=max_step_log_pressure,
        events=radius_limit_event, dense_output=True,
    )
    if not outer_sol.success:
        raise RuntimeError(f"outer TOV segment failed before surface cutoff: {outer_sol.message}")
    if outer_sol.t[-1] < q_surface * (1.0 - 1e-12):
        raise RuntimeError("outer TOV segment reached r_limit_m before surface cutoff")
    radius, ms, nu_raw = np.asarray(outer_sol.y[:, -1], dtype=float)
    compact = 2.0 * G_SI * ms / (radius * C_SI**2)
    if compact <= 0 or compact >= 1:
        raise ArithmeticError("invalid surface compactness")
    nu_shift = 0.5 * np.log1p(-compact) - nu_raw

    def make_segment(side: Literal["high", "low"], sol, q0: float, q1: float,
                     p_ref: float) -> TOVSegment:
        q_grid = np.linspace(q0, q1, max(1200, 8 * len(sol.t)))
        r_grid = np.asarray(sol.sol(q_grid)[0], dtype=float)
        if np.any(np.diff(r_grid) <= 0):
            raise ArithmeticError("TOV radius is not strictly increasing in q")
        q_map = PchipInterpolator(r_grid, q_grid, extrapolate=False)
        return TOVSegment(side, sol, float(r_grid[0]), float(r_grid[-1]), p_ref, q_map)

    inner_seg = make_segment("high", inner_sol, 0.0, q_transition, p0)
    outer_seg = make_segment("low", outer_sol, 0.0, q_surface, eos.transition_pressure_pa)

    # One-sided thermodynamics at the exact same geometric interface.
    inner_bg = _background_from_state(eos, rt, eos.transition_pressure_pa,
                                      float(nut) + nu_shift, mt, "high")
    outer_bg = _background_from_state(eos, rt, eos.transition_pressure_pa,
                                      float(nut) + nu_shift, mt, "low")
    interface = MaxwellInterface(
        pressure_pa=eos.transition_pressure_pa,
        chemical_potential_j=eos.mu_transition_mev * MEV_TO_J,
        energy_inner_j_m3=inner_bg["energy_j_m3"],
        energy_outer_j_m3=outer_bg["energy_j_m3"],
        density_inner_m3=inner_bg["number_density_m3"],
        density_outer_m3=outer_bg["number_density_m3"],
        radius_m=rt,
    )
    if abs(inner_bg["m_kg"] - outer_bg["m_kg"]) > 1e-12 * max(abs(ms), 1.0):
        raise ArithmeticError("mass is discontinuous at Maxwell interface")
    if abs(inner_bg["nu"] - outer_bg["nu"]) > 1e-12:
        raise ArithmeticError("metric potential is discontinuous at Maxwell interface")
    return MaxwellStar(
        eos=eos,
        inner=inner_seg,
        outer=outer_seg,
        r0_m=r0,
        interface=interface,
        radius_m=radius,
        mass_kg=float(ms),
        mass_solar=float(ms / M_SUN_NOMINAL_KG),
        p_c_pa=pc,
        p_cut_pa=pcut,
        nu_shift=float(nu_shift),
        rtol=rtol,
        accepted_steps=int(len(inner_sol.t) + len(outer_sol.t) - 2),
    )


def _background_from_state(eos: MaxwellEOS, r: float, p: float, nu: float,
                            m: float, side: str) -> dict[str, float]:
    m, p = float(m), float(p)
    therm = eos.evaluate(p, side)  # type: ignore[arg-type]
    compact = 2.0 * G_SI * m / (r * C_SI**2)
    lam = -0.5 * np.log1p(-compact)
    nu_prime = (G_SI / C_SI**2) * (
        m / r**2 + 4.0 * np.pi * r * p / C_SI**2
    ) / (1.0 - compact)
    return {"r_m": r, "m_kg": m, **therm, "nu": nu, "lambda": float(lam), "dnu_dr": float(nu_prime)}


def _cowling_rhs(star: MaxwellStar, r: float, y: np.ndarray, omega: float, side: str) -> np.ndarray:
    bg = star.background(r, "inner" if side == "high" else "outer")
    w, v = map(float, y)
    chi, nu, lam, dnu = bg["chi"], bg["nu"], bg["lambda"], bg["dnu_dr"]
    expo = lam - 2.0 * nu
    if not np.isfinite(expo) or abs(expo) > 700:
        raise ArithmeticError("Cowling metric exponential is outside floating-point range")
    a = (omega * r / C_SI)**2 * np.exp(expo)
    qterm = a * v + dnu * w
    return np.array([qterm / chi - 6.0 * np.exp(lam) * v,
                     2.0 * dnu * v - np.exp(lam) * w / r**2])


def _interface_transfer(
    star: MaxwellStar, omega: float, y_inner: np.ndarray, conversion: Conversion
) -> tuple[np.ndarray, dict[str, float]]:
    rt = star.interface.radius_m
    bi = star.background(rt, "inner")
    bo = star.background(rt, "outer")
    w_i, v_i = map(float, y_inner)
    a = (omega * rt / C_SI)**2 * np.exp(bi["lambda"] - 2.0 * bi["nu"])
    if a <= 0 or not np.isfinite(a):
        raise ArithmeticError("invalid interface frequency coefficient")
    q_i = a * v_i + bi["dnu_dr"] * w_i
    h_i = bi["energy_j_m3"] + bi["pressure_pa"]
    h_o = bo["energy_j_m3"] + bo["pressure_pa"]
    if h_i <= 0 or h_o <= 0:
        raise ArithmeticError("nonpositive interface enthalpy density")
    q_o = (h_i / h_o) * q_i  # [Delta P]=0, since exp(-lambda)/r^2 is continuous.
    if conversion == "slow":
        # [xi^r]=0 -> [W]=0 in the present W convention.
        w_o = w_i
    elif conversion == "rapid":
        # [xi^r - Delta P/P'_0]=0 -> [W - qterm/nu']=0.
        dnu = bi["dnu_dr"]
        if dnu <= 0 or abs(dnu - bo["dnu_dr"]) > 2e-12 * max(abs(dnu), 1.0):
            raise ArithmeticError("metric derivative is not continuous at the interface")
        w_o = w_i + (q_o - q_i) / dnu
    else:
        raise ValueError("conversion must be 'slow' or 'rapid'")
    v_o = (q_o - bo["dnu_dr"] * w_o) / a
    y_outer = np.array([w_o, v_o], dtype=float)

    delta_i = -h_i * np.exp(-bi["lambda"]) * q_i / rt**2
    delta_o = -h_o * np.exp(-bo["lambda"]) * q_o / rt**2
    pprime_i = -h_i * bi["dnu_dr"]
    pprime_o = -h_o * bo["dnu_dr"]
    xi_i = np.exp(-bi["lambda"]) * w_i / rt**2
    xi_o = np.exp(-bo["lambda"]) * w_o / rt**2
    isobar_i = xi_i - delta_i / pprime_i
    isobar_o = xi_o - delta_o / pprime_o
    diag = {
        "w_inner": w_i, "v_inner": v_i, "w_outer": w_o, "v_outer": v_o,
        "delta_p_inner_pa": delta_i, "delta_p_outer_pa": delta_o,
        "delta_p_relative_jump": abs(delta_o - delta_i) / max(abs(delta_i), abs(delta_o), 1.0),
        "radial_displacement_jump_m": xi_o - xi_i,
        "isobar_displacement_jump_m": isobar_o - isobar_i,
        "slow_condition_residual": abs(xi_o - xi_i) / max(abs(xi_o), abs(xi_i), 1.0),
        "rapid_condition_residual": abs(isobar_o - isobar_i) / max(abs(isobar_o), abs(isobar_i), 1.0),
    }
    if diag["delta_p_relative_jump"] > 2e-10:
        raise ArithmeticError("Lagrangian pressure perturbation is not continuous")
    required = "slow_condition_residual" if conversion == "slow" else "rapid_condition_residual"
    if diag[required] > 2e-10:
        raise ArithmeticError(f"{conversion} conversion interface condition failed")
    return y_outer, diag


@dataclass(frozen=True)
class MaxwellCowlingShot:
    frequency_hz: float
    conversion: Conversion
    residual: float
    radius_inner_m: np.ndarray
    w_inner_m3: np.ndarray
    v_inner_m2: np.ndarray
    radius_outer_m: np.ndarray
    w_outer_m3: np.ndarray
    v_outer_m2: np.ndarray
    interface_diagnostics: dict[str, float]
    surface_diagnostics: dict[str, float]
    accepted_steps: int


def shoot_maxwell_cowling(
    star: MaxwellStar,
    frequency_hz: float,
    *,
    conversion: Conversion = "slow",
    ell: int = 2,
    rtol: float = 2.0e-9,
    atol_w: float = 1.0e-10,
    atol_v: float = 1.0e-10,
    max_step_m: float = 100.0,
    output_points_per_segment: int = 600,
) -> MaxwellCowlingShot:
    """Shoot a real-frequency, polar Cowling mode across a sharp interface."""
    if ell != 2:
        raise NotImplementedError("the current native Cowling equations are quadrupolar ell=2 only")
    if conversion not in ("slow", "rapid"):
        raise ValueError("conversion must be 'slow' or 'rapid'")
    f = float(frequency_hz)
    if not np.isfinite(f) or f <= 0:
        raise ValueError("frequency must be finite and positive")
    if output_points_per_segment < 20 or max_step_m <= 0 or not (0 < rtol < 1):
        raise ValueError("invalid Cowling integration controls")
    omega = 2.0 * np.pi * f
    r0, rt, radius = star.r0_m, star.interface.radius_m, star.radius_m
    y0 = np.array([r0**3, -(r0**2)/2.0], dtype=float)
    inner = solve_ivp(
        lambda r, y: _cowling_rhs(star, r, y, omega, "high"),
        (r0, rt), y0, method="DOP853", rtol=rtol, atol=[atol_w, atol_v],
        max_step=max_step_m, dense_output=True,
    )
    if not inner.success:
        raise RuntimeError(f"inner Cowling segment failed: {inner.message}")
    y_t_inner = np.asarray(inner.sol(rt), dtype=float)
    y_t_outer, idiag = _interface_transfer(star, omega, y_t_inner, conversion)
    outer = solve_ivp(
        lambda r, y: _cowling_rhs(star, r, y, omega, "low"),
        (rt, radius), y_t_outer, method="DOP853", rtol=rtol, atol=[atol_w, atol_v],
        max_step=max_step_m, dense_output=True,
    )
    if not outer.success:
        raise RuntimeError(f"outer Cowling segment failed: {outer.message}")
    ri = np.linspace(r0, rt, output_points_per_segment)
    ro = np.linspace(rt, radius, output_points_per_segment)
    yi = np.asarray(inner.sol(ri), dtype=float)
    yo = np.asarray(outer.sol(ro), dtype=float)
    surface = star.background(radius, "outer")
    w_s, v_s = map(float, yo[:, -1])
    b1 = radius * surface["dnu_dr"]
    b2 = (omega * radius / C_SI)**2 * np.exp(surface["lambda"] - 2.0 * surface["nu"])
    z = np.array([w_s / radius**3, v_s / radius**2])
    b = np.array([b1, b2])
    norm_b, norm_z = np.linalg.norm(b), np.linalg.norm(z)
    if norm_b <= 0 or norm_z <= 0 or not np.all(np.isfinite([norm_b, norm_z])):
        raise ArithmeticError("invalid surface residual normalization")
    residual = float(np.dot(b / norm_b, z / norm_z))
    surface_diag = {
        "boundary_unscaled": float(b1 * z[0] + b2 * z[1]),
        "boundary_normalized_residual": residual,
        "surface_pressure_pa": surface["pressure_pa"],
    }
    return MaxwellCowlingShot(
        frequency_hz=f, conversion=conversion, residual=residual,
        radius_inner_m=ri, w_inner_m3=yi[0], v_inner_m2=yi[1],
        radius_outer_m=ro, w_outer_m3=yo[0], v_outer_m2=yo[1],
        interface_diagnostics=idiag, surface_diagnostics=surface_diag,
        accepted_steps=int(len(inner.t) + len(outer.t) - 2),
    )


def _sampled_segment_nodes(values: np.ndarray, relative_zero_tol: float = 128.0 * np.finfo(float).eps) -> int:
    """Count resolved sign transitions in one continuous phase segment only.

    Samples smaller than a segment-relative threshold are treated as zeros.
    The function never compares the last inner-phase sample to the first
    outer-phase sample, so a finite interface jump is not a radial node.
    """
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1 or arr.size == 0 or np.any(~np.isfinite(arr)):
        raise ValueError("node-count segment must be a nonempty finite 1-D array")
    tol = float(relative_zero_tol)
    if not np.isfinite(tol) or tol < 0.0 or tol >= 0.1:
        raise ValueError("relative_zero_tol must satisfy 0 <= tol < 0.1")
    scale = float(np.max(np.abs(arr)))
    if scale <= np.finfo(float).tiny:
        return 0
    resolved = arr[np.abs(arr) > tol * scale]
    if resolved.size < 2:
        return 0
    return int(np.count_nonzero(np.signbit(resolved[1:]) != np.signbit(resolved[:-1])))


def _sampled_phase_nodes(
    w_inner: np.ndarray,
    w_outer: np.ndarray,
    v_inner: np.ndarray,
    v_outer: np.ndarray,
    relative_zero_tol: float = 128.0 * np.finfo(float).eps,
) -> dict[str, int]:
    """Return sampled W/V node counts separately in each sharp-interface phase."""
    arrays = {
        "w_inner": np.asarray(w_inner, dtype=float),
        "w_outer": np.asarray(w_outer, dtype=float),
        "v_inner": np.asarray(v_inner, dtype=float),
        "v_outer": np.asarray(v_outer, dtype=float),
    }
    if any(x.ndim != 1 or x.size < 2 for x in arrays.values()):
        raise ValueError("each phase-segment profile must contain at least two samples")
    result = {
        "w_nodes_inner_sampled": _sampled_segment_nodes(arrays["w_inner"], relative_zero_tol),
        "w_nodes_outer_sampled": _sampled_segment_nodes(arrays["w_outer"], relative_zero_tol),
        "v_nodes_inner_sampled": _sampled_segment_nodes(arrays["v_inner"], relative_zero_tol),
        "v_nodes_outer_sampled": _sampled_segment_nodes(arrays["v_outer"], relative_zero_tol),
    }
    result["w_nodes_total_sampled"] = (
        result["w_nodes_inner_sampled"] + result["w_nodes_outer_sampled"]
    )
    result["v_nodes_total_sampled"] = (
        result["v_nodes_inner_sampled"] + result["v_nodes_outer_sampled"]
    )
    return result


def scan_maxwell_cowling_modes(
    star: MaxwellStar,
    frequency_grid_hz: np.ndarray,
    *,
    conversion: Conversion = "slow",
    rtol: float = 2e-9,
    max_step_m: float = 100.0,
) -> list[dict[str, object]]:
    """Find real-frequency roots from sign changes on a caller-supplied grid.

    The caller must inspect node counts and mode overlaps before assigning f/p1
    labels; roots are returned unlabelled to avoid order-based mode identity.
    """
    grid = np.asarray(frequency_grid_hz, dtype=float)
    if grid.ndim != 1 or grid.size < 3 or np.any(~np.isfinite(grid)) or np.any(grid <= 0) or np.any(np.diff(grid) <= 0):
        raise ValueError("frequency_grid_hz must be a strictly increasing positive 1-D grid")
    residuals = np.empty(grid.size)
    for i, f in enumerate(grid):
        residuals[i] = shoot_maxwell_cowling(star, float(f), conversion=conversion,
                                             rtol=rtol, max_step_m=max_step_m,
                                             output_points_per_segment=80).residual
    roots: list[dict[str, object]] = []
    for i in range(grid.size - 1):
        if residuals[i] == 0:
            root_hz = float(grid[i])
        elif residuals[i] * residuals[i + 1] < 0:
            root_hz = float(brentq(
                lambda f: shoot_maxwell_cowling(
                    star, float(f), conversion=conversion, rtol=rtol,
                    max_step_m=max_step_m, output_points_per_segment=80,
                ).residual,
                float(grid[i]), float(grid[i + 1]), xtol=1e-7, rtol=1e-10,
            ))
        else:
            continue
        shot = shoot_maxwell_cowling(star, root_hz, conversion=conversion,
                                     rtol=rtol, max_step_m=max_step_m)
        # Count W and V separately in each phase. The interface transfer is
        # one-sided, and the finite jump must never be counted as an interior node.
        nodes = _sampled_phase_nodes(
            shot.w_inner_m3, shot.w_outer_m3, shot.v_inner_m2, shot.v_outer_m2
        )
        roots.append({
            "frequency_hz": root_hz,
            "residual": shot.residual,
            **nodes,
            "conversion": conversion,
            "interface_diagnostics": shot.interface_diagnostics,
            "shot": shot,
        })
    return roots
