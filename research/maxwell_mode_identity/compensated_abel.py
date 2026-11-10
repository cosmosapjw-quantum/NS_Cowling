"""Bulk and interface compensation for a specified 2x2 Cowling system.

This is an algebraic Abel-Liouville diagnostic, not a canonical symplectic
current. It is valid only when the supplied ODE matrix has the documented
trace and no additional composition-gradient terms have been omitted.
"""
from __future__ import annotations
import math
from collections.abc import Iterable
import numpy as np

def _finite_real_array(value: object, name: str) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite real values")
    return arr

def closure_drift_density(nu_prime_per_m: object, chi_ad: object, chi_eq: object) -> np.ndarray | float:
    """Return d(log Omega)/dr = nu'*(1/chi_ad - 1/chi_eq).

    The formula assumes the exact matrix trace tr(M)=nu'*(2+1/chi_ad)
    and the equilibrium first law h'/h=-(1+1/chi_eq)*nu'.
    """
    g = _finite_real_array(nu_prime_per_m, "nu_prime_per_m")
    ca = _finite_real_array(chi_ad, "chi_ad")
    ce = _finite_real_array(chi_eq, "chi_eq")
    try:
        g, ca, ce = np.broadcast_arrays(g, ca, ce)
    except ValueError as exc:
        raise ValueError("nu_prime_per_m, chi_ad and chi_eq are not broadcast-compatible") from exc
    if np.any(ca <= 0.0) or np.any(ce <= 0.0):
        raise ValueError("chi_ad and chi_eq must be strictly positive")
    out = g * (1.0 / ca - 1.0 / ce)
    return float(out) if out.ndim == 0 else out

def integrated_closure_drift(r_m: object, nu_prime_per_m: object, chi_ad: object, chi_eq: object) -> float:
    """Trapezoid-integrate the closure drift over one smooth phase segment."""
    r = _finite_real_array(r_m, "r_m")
    g = _finite_real_array(nu_prime_per_m, "nu_prime_per_m")
    ca = _finite_real_array(chi_ad, "chi_ad")
    ce = _finite_real_array(chi_eq, "chi_eq")
    if r.ndim != 1 or r.size < 2 or np.any(np.diff(r) <= 0.0):
        raise ValueError("r_m must be one-dimensional, strictly increasing and length >= 2")
    try:
        g, ca, ce = np.broadcast_arrays(g, ca, ce)
    except ValueError as exc:
        raise ValueError("drift arrays are not broadcast-compatible") from exc
    if g.ndim != 1 or g.size != r.size:
        raise ValueError("drift arrays must be one-dimensional and match r_m length")
    density = np.asarray(closure_drift_density(g, ca, ce), dtype=float)
    trapz = getattr(np, "trapezoid", np.trapz)
    result = float(trapz(density, r))
    if not math.isfinite(result):
        raise ArithmeticError("integrated closure drift is not finite")
    return result

def interface_multiplier(h_inner_j_m3: float, h_outer_j_m3: float, nu_inner: float,
                          nu_outer: float, determinant_transfer: complex | float) -> complex | float:
    """Return J=Omega_outer/Omega_inner for y_outer=T*y_inner."""
    hi, ho, ni, no = map(float, (h_inner_j_m3, h_outer_j_m3, nu_inner, nu_outer))
    det_t = complex(determinant_transfer)
    if not all(map(math.isfinite, (hi, ho, ni, no))):
        raise ValueError("h and lapse values must be finite")
    if hi <= 0.0 or ho <= 0.0:
        raise ValueError("one-sided enthalpy densities must be positive")
    if not (math.isfinite(det_t.real) and math.isfinite(det_t.imag)) or abs(det_t) == 0.0:
        raise ValueError("determinant_transfer must be finite and nonzero")
    value = (ho / hi) * math.exp(-(no - ni)) * det_t
    if not (math.isfinite(value.real) and math.isfinite(value.imag)):
        raise ArithmeticError("interface multiplier overflowed or became non-finite")
    if not np.iscomplexobj(determinant_transfer) and value.imag == 0.0:
        return float(value.real)
    return value

def compensate_concomitant(omega: complex | float, bulk_drift_integral: float,
                            crossed_interface_multipliers: Iterable[complex | float] = ()) -> complex | float:
    """Return K=Omega*exp(-integral(delta))*product(J_k**-1).

    The bulk integral runs from the fixed reference point to the evaluation
    point, and interface multipliers are ordered from reference point outward.
    Each crossed interface must appear exactly once.
    """
    om = complex(omega)
    integral = float(bulk_drift_integral)
    if not (math.isfinite(om.real) and math.isfinite(om.imag) and math.isfinite(integral)):
        raise ValueError("omega and bulk_drift_integral must be finite")
    product = 1.0 + 0.0j
    for j, factor in enumerate(crossed_interface_multipliers):
        fj = complex(factor)
        if not (math.isfinite(fj.real) and math.isfinite(fj.imag)) or abs(fj) == 0.0:
            raise ValueError(f"interface multiplier {j} must be finite and nonzero")
        product *= fj
    value = om * math.exp(-integral) / product
    if not (math.isfinite(value.real) and math.isfinite(value.imag)):
        raise ArithmeticError("compensated concomitant is not finite")
    if not np.iscomplexobj(omega) and value.imag == 0.0:
        return float(value.real)
    return value
