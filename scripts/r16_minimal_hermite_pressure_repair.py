"""Minimal per-cell pressure-increment repair for frozen cubic-Hermite EOS admission."""
from __future__ import annotations

import math
from typing import Any


def minimal_admissible_delta_p(
    *, mu0: float, n0: float, mu1: float, n1: float,
    delta_p_source: float, interior_margin_fraction: float = 1.0e-10,
) -> dict[str, Any]:
    """Project a source pressure increment to the nearest Fortran-admissible interval.

    Holding mu0,mu1,n0,n1 fixed, the exact native-Hermite conditions
    P''(0)>0, P''(1)>0, G(0)>=0, G(1)>=0 reduce to delta_p in [L,U],
    where the causal margins are the stricter bounds:
      L = h(2*n0+n1)/3 + h^2*n0/(6*mu0)
      U = h(n0+2*n1)/3 - h^2*n1/(6*mu1).
    A tiny interior offset is used only when clipping, so strict P'' and
    nonnegative causal margins survive floating-point evaluation.
    """
    mu0, n0, mu1, n1, delta_p_source = map(
        float, (mu0, n0, mu1, n1, delta_p_source)
    )
    values = [mu0, n0, mu1, n1, delta_p_source, interior_margin_fraction]
    if any(not math.isfinite(x) for x in values) or min(values[:5]) <= 0.0:
        raise ValueError("cell input values must be finite and positive")
    if mu1 <= mu0:
        raise ValueError("mu endpoints must be strictly ordered")
    if not (0.0 < interior_margin_fraction < 0.01):
        raise ValueError("interior margin fraction must be in (0,0.01)")

    h = mu1 - mu0
    lower_stability = h * (2.0 * n0 + n1) / 3.0
    upper_stability = h * (n0 + 2.0 * n1) / 3.0
    lower_causality = lower_stability + h * h * n0 / (6.0 * mu0)
    upper_causality = upper_stability - h * h * n1 / (6.0 * mu1)
    lower = lower_causality
    upper = upper_causality
    width = upper - lower
    if not math.isfinite(width) or width <= 0.0:
        raise ValueError(
            f"empty cubic-Hermite admissibility interval: lower={lower}, upper={upper}"
        )

    margin = interior_margin_fraction * width
    if delta_p_source < lower:
        repaired = lower + margin
    elif delta_p_source > upper:
        repaired = upper - margin
    else:
        repaired = delta_p_source
    return {
        "delta_p_source": delta_p_source,
        "delta_p_lower": lower,
        "delta_p_upper": upper,
        "delta_p_width": width,
        "delta_p_repaired": repaired,
        "delta_p_adjustment": repaired - delta_p_source,
        "was_adjusted": repaired != delta_p_source,
        "lower_stability_bound": lower_stability,
        "upper_stability_bound": upper_stability,
        "lower_causality_bound": lower_causality,
        "upper_causality_bound": upper_causality,
        "interior_margin_fraction": interior_margin_fraction,
        "interior_margin": margin,
    }


def audit_hermite_cell(
    *, mu0: float, p0: float, n0: float,
    mu1: float, p1: float, n1: float,
) -> dict[str, Any]:
    """Exact source-matched endpoint conditions for a cubic-Hermite P(mu) cell."""
    mu0, p0, n0, mu1, p1, n1 = map(
        float, (mu0, p0, n0, mu1, p1, n1)
    )
    if any(not math.isfinite(x) or x <= 0.0 for x in (mu0,p0,n0,mu1,p1,n1)):
        raise ValueError("cell endpoints must be finite and positive")
    if mu1 <= mu0 or p1 <= p0:
        raise ValueError("cell endpoints must be strictly ordered")
    h = mu1 - mu0
    delta = p1 - p0
    a2 = 3.0 * delta - h * (2.0 * n0 + n1)
    a3 = 2.0 * (p0 - p1) + h * (n0 + n1)
    d20 = 2.0 * a2 / (h*h)
    d21 = (2.0*a2 + 6.0*a3)/(h*h)
    g0 = mu0*d20 - n0
    g1 = mu1*d21 - n1
    stable = min(d20,d21) > 0.0
    causal = stable and min(g0,g1) >= 0.0
    return {
        "status":"PASS_NATIVE_CELL_ADMISSIBLE" if causal else (
            "REJECT_NONPOSITIVE_COMPRESSIBILITY" if not stable else "REJECT_ACAUSAL_NATIVE_CELL"
        ),
        "stable_interval":stable,
        "causal_interval":causal,
        "d2p_low":d20,
        "d2p_high":d21,
        "causal_margin_low":g0,
        "causal_margin_high":g1,
        "min_d2p":min(d20,d21),
        "min_causal_margin":min(g0,g1),
    }
