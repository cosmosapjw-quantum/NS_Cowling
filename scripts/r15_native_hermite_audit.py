"""Read-only audit of solver-native cubic-Hermite P(mu) EOS intervals.

The algebra and interval acceptance reproduce ns_mu_bridge.f90:
P(t)=a0+a1 t+a2 t^2+a3 t^3,
n=dP/dmu, epsilon=mu*n-P, chi=n/(mu*dn/dmu).
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np
from numpy.polynomial import Polynomial


def _finite_positive(name: str, value: float) -> float:
    x = float(value)
    if not math.isfinite(x) or x <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return x


def _chi_value(t: float, *, mu0: float, h: float, a1: float, a2: float, a3: float) -> float:
    mu = mu0 + h * t
    q = a1 + 2.0 * a2 * t + 3.0 * a3 * t * t
    r = 2.0 * a2 + 6.0 * a3 * t
    return h * q / (mu * r)


def audit_native_interval(
    mu0_j: float, p0_pa: float, n0_m3: float,
    mu1_j: float, p1_pa: float, n1_m3: float,
) -> dict[str, Any]:
    """Audit one native cubic-Hermite P(mu) cell using the production formulas.

    Whole-interval checks match the Fortran constructor:
      min(P'') = min(P''(0), P''(1)) because P'' is linear;
      min(mu P'' - P') = min(endpoint values) because its derivative is
      mu P''' and P''' is constant.
    """
    mu0 = _finite_positive("mu0_j", mu0_j)
    p0 = _finite_positive("p0_pa", p0_pa)
    n0 = _finite_positive("n0_m3", n0_m3)
    mu1 = _finite_positive("mu1_j", mu1_j)
    p1 = _finite_positive("p1_pa", p1_pa)
    n1 = _finite_positive("n1_m3", n1_m3)
    if mu1 <= mu0 or p1 <= p0:
        raise ValueError("native interval requires strictly increasing mu and pressure")

    h = mu1 - mu0
    a0 = p0
    a1 = h * n0
    a2 = 3.0 * (p1 - p0) - h * (2.0 * n0 + n1)
    a3 = 2.0 * (p0 - p1) + h * (n0 + n1)

    d2_low = 2.0 * a2 / (h * h)
    d2_high = (2.0 * a2 + 6.0 * a3) / (h * h)
    causal_margin_low = mu0 * d2_low - n0
    causal_margin_high = mu1 * d2_high - n1

    min_d2p = min(d2_low, d2_high)
    min_causal_margin = min(causal_margin_low, causal_margin_high)
    stable = min_d2p > 0.0
    causal = stable and min_causal_margin >= 0.0

    epsilon_low = mu0 * n0 - p0
    epsilon_high = mu1 * n1 - p1
    if min(epsilon_low, epsilon_high) <= 0.0:
        status = "REJECT_NONPOSITIVE_ENERGY_DENSITY"
    elif not stable:
        status = "REJECT_NONPOSITIVE_COMPRESSIBILITY"
    elif not causal:
        status = "REJECT_ACAUSAL_NATIVE_CELL"
    else:
        status = "PASS_NATIVE_CELL_ADMISSIBLE"

    chi_min = None
    chi_max = None
    chi_min_t = None
    chi_max_t = None
    stationary_t: list[float] = []
    chi_candidates: list[tuple[float, float]] = []

    if stable:
        # chi(t)=h*q(t)/(mu(t)*r(t)); its derivative numerator is
        # mu*r^2 - h*q*r - mu*q*r', a polynomial of degree <= 3.
        t_poly = Polynomial([0.0, 1.0])
        mu_poly = Polynomial([mu0, h])
        q_poly = Polynomial([a1, 2.0 * a2, 3.0 * a3])
        r_poly = Polynomial([2.0 * a2, 6.0 * a3])
        rprime_poly = Polynomial([6.0 * a3])
        derivative_numerator = (
            mu_poly * r_poly * r_poly
            - h * q_poly * r_poly
            - mu_poly * q_poly * rprime_poly
        )
        for root in derivative_numerator.roots():
            real = float(np.real(root))
            imag = abs(float(np.imag(root)))
            if imag <= 1.0e-10 and 1.0e-12 < real < 1.0 - 1.0e-12:
                if not any(abs(real - previous) <= 1.0e-10 for previous in stationary_t):
                    stationary_t.append(real)
        stationary_t.sort()
        for t in [0.0, *stationary_t, 1.0]:
            chi_candidates.append((t, _chi_value(t, mu0=mu0, h=h, a1=a1, a2=a2, a3=a3)))
        chi_min_t, chi_min = min(chi_candidates, key=lambda item: item[1])
        chi_max_t, chi_max = max(chi_candidates, key=lambda item: item[1])

    return {
        "status": status,
        "stable_interval": stable,
        "causal_interval": causal,
        "interval_extrema_checked": True,
        "mu_span_j": h,
        "a0_pa": a0,
        "a1_pa": a1,
        "a2_pa": a2,
        "a3_pa": a3,
        "d2p_low_pa_per_j2": d2_low,
        "d2p_high_pa_per_j2": d2_high,
        "min_d2p_pa_per_j2": min_d2p,
        "causal_margin_low_m3": causal_margin_low,
        "causal_margin_high_m3": causal_margin_high,
        "min_causal_margin_m3": min_causal_margin,
        "chi_low": _chi_value(0.0, mu0=mu0, h=h, a1=a1, a2=a2, a3=a3) if stable else None,
        "chi_high": _chi_value(1.0, mu0=mu0, h=h, a1=a1, a2=a2, a3=a3) if stable else None,
        "chi_min": chi_min,
        "chi_min_t": chi_min_t,
        "chi_max": chi_max,
        "chi_max_t": chi_max_t,
        "chi_stationary_t": stationary_t,
        "chi_stationary_count": len(stationary_t),
        "epsilon_low_pa": epsilon_low,
        "epsilon_high_pa": epsilon_high,
        "secant_dlnp_dlnn": math.log(p1 / p0) / math.log(n1 / n0),
        "density_ratio": n1 / n0,
        "pressure_ratio": p1 / p0,
        "energy_density_ratio": epsilon_high / epsilon_low,
        "mu_change_fraction": mu1 / mu0 - 1.0,
    }
