"""Reproduce the frozen q-isotonic projection and audit all 234 EOS tables."""
from __future__ import annotations

import math
from typing import Any, Sequence


def _finite_positive(values: Sequence[float], name: str) -> list[float]:
    out = [float(x) for x in values]
    if any(not math.isfinite(x) or x <= 0.0 for x in out):
        raise ValueError(f"{name} values must be finite and positive")
    return out


def pav_increasing(values: Sequence[float], weights: Sequence[float]) -> tuple[list[float], list[dict[str, Any]]]:
    """Weighted pool-adjacent-violators algorithm for a nondecreasing sequence."""
    y = [float(x) for x in values]
    w = _finite_positive(weights, "PAV weight")
    if len(y) != len(w) or len(y) < 2:
        raise ValueError("PAV values and weights must have equal length >= 2")
    if any(not math.isfinite(x) for x in y):
        raise ValueError("PAV values must be finite")

    blocks: list[dict[str, float | int]] = []
    for i, (value, weight) in enumerate(zip(y, w)):
        blocks.append({"lo": i, "hi": i + 1, "weight": weight, "weighted_sum": weight * value})
        while len(blocks) >= 2:
            a, b = blocks[-2], blocks[-1]
            mean_a = float(a["weighted_sum"]) / float(a["weight"])
            mean_b = float(b["weighted_sum"]) / float(b["weight"])
            if mean_a <= mean_b:
                break
            blocks[-2:] = [{
                "lo": int(a["lo"]),
                "hi": int(b["hi"]),
                "weight": float(a["weight"]) + float(b["weight"]),
                "weighted_sum": float(a["weighted_sum"]) + float(b["weighted_sum"]),
            }]

    projected = [0.0] * len(y)
    normalized_blocks: list[dict[str, Any]] = []
    for block in blocks:
        lo, hi = int(block["lo"]), int(block["hi"])
        mean = float(block["weighted_sum"]) / float(block["weight"])
        for i in range(lo, hi):
            projected[i] = mean
        normalized_blocks.append({
            "lo_zero_based": lo,
            "hi_exclusive_zero_based": hi,
            "width": hi - lo,
            "weight": float(block["weight"]),
            "mean": mean,
        })
    return projected, normalized_blocks


def hermite_cell_diagnostics(
    mu0: float, p0: float, n0: float,
    mu1: float, p1: float, n1: float,
) -> dict[str, Any]:
    """Whole-cell endpoint curvature and causal-margin diagnostics for cubic P(mu)."""
    if not (mu1 > mu0 and p1 > p0 and mu0 > 0.0 and mu1 > 0.0):
        raise ValueError("Hermite cell mu and pressure must be positive and strictly increasing")
    h = mu1 - mu0
    dp_over_h = (p1 - p0) / h
    d2_low = (6.0 * dp_over_h - 4.0 * n0 - 2.0 * n1) / h
    d2_high = (2.0 * n0 + 4.0 * n1 - 6.0 * dp_over_h) / h
    g_low = mu0 * d2_low - n0
    g_high = mu1 * d2_high - n1
    stable = d2_low > 0.0 and d2_high > 0.0
    causal = stable and g_low >= 0.0 and g_high >= 0.0
    return {
        "d2p_low": d2_low,
        "d2p_high": d2_high,
        "min_d2p": min(d2_low, d2_high),
        "causal_margin_low": g_low,
        "causal_margin_high": g_high,
        "min_causal_margin": min(g_low, g_high),
        "stable": stable,
        "causal": causal,
    }


def reconstruct_projection(
    mu: Sequence[float],
    n_source: Sequence[float],
    p_source: Sequence[float],
    eps_source: Sequence[float],
) -> dict[str, Any]:
    """Reproduce project_isotonic.py: q-PAV, n=mu*q, exact integral P, eps=mu*n-P."""
    mu_v = _finite_positive(mu, "mu")
    n_v = _finite_positive(n_source, "n_source")
    p_v = _finite_positive(p_source, "p_source")
    eps_v = _finite_positive(eps_source, "eps_source")
    if not (len(mu_v) == len(n_v) == len(p_v) == len(eps_v)) or len(mu_v) < 3:
        raise ValueError("source EOS arrays must have equal length >= 3")

    h = [b - a for a, b in zip(mu_v, mu_v[1:])]
    if any(x <= 0.0 for x in h) or any(b <= a for a, b in zip(p_v, p_v[1:])):
        raise ValueError("mu and source pressure must be strictly increasing")

    q_source = [n / m for n, m in zip(n_v, mu_v)]
    weights = [0.0] * len(mu_v)
    weights[0] = 0.5 * h[0]
    weights[-1] = 0.5 * h[-1]
    for i in range(1, len(mu_v) - 1):
        weights[i] = 0.5 * (h[i - 1] + h[i])

    q_projected, blocks = pav_increasing(q_source, weights)
    if any(x <= 0.0 for x in q_projected) or any(
        b < a - 1.0e-14 * max(q_projected) for a, b in zip(q_projected, q_projected[1:])
    ):
        raise ValueError("isotonic projection failed positive monotonic q")

    n_projected = [m * q for m, q in zip(mu_v, q_projected)]
    p_projected = [p_v[0]]
    for i, step in enumerate(h):
        dq = q_projected[i + 1] - q_projected[i]
        increment = step * (
            (mu_v[i] + 0.5 * step) * q_projected[i]
            + (0.5 * mu_v[i] + step / 3.0) * dq
        )
        p_projected.append(p_projected[-1] + increment)
    eps_projected = [
        m * n - p for m, n, p in zip(mu_v, n_projected, p_projected)
    ]
    if any(x <= 0.0 for x in p_projected + eps_projected) or any(
        b <= a for a, b in zip(p_projected, p_projected[1:])
    ) or any(b <= a for a, b in zip(eps_projected, eps_projected[1:])):
        raise ValueError("projected thermodynamics are not positive and monotone")

    source_cell_diagnostics = [
        hermite_cell_diagnostics(mu_v[i], p_v[i], n_v[i], mu_v[i + 1], p_v[i + 1], n_v[i + 1])
        for i in range(len(mu_v) - 1)
    ]
    projected_cell_diagnostics = [
        hermite_cell_diagnostics(
            mu_v[i], p_projected[i], n_projected[i],
            mu_v[i + 1], p_projected[i + 1], n_projected[i + 1],
        )
        for i in range(len(mu_v) - 1)
    ]
    source_nonstable_cells = sum(not x["stable"] for x in source_cell_diagnostics)
    source_noncausal_cells = sum(not x["causal"] for x in source_cell_diagnostics)
    projected_nonstable_cells = sum(not x["stable"] for x in projected_cell_diagnostics)
    projected_noncausal_cells = sum(not x["causal"] for x in projected_cell_diagnostics)
    source_boundary_cell = source_cell_diagnostics[216] if len(source_cell_diagnostics) > 216 else None
    projected_boundary_cell = projected_cell_diagnostics[216] if len(projected_cell_diagnostics) > 216 else None

    source_pressure_secants = [
        (p_v[i + 1] - p_v[i]) / h[i] for i in range(len(h))
    ]
    projected_pressure_secants = [
        (p_projected[i + 1] - p_projected[i]) / h[i] for i in range(len(h))
    ]
    source_secant_change_pct = [None] + [
        100.0 * (source_pressure_secants[i] / source_pressure_secants[i - 1] - 1.0)
        for i in range(1, len(source_pressure_secants))
    ]
    projected_secant_change_pct = [None] + [
        100.0 * (projected_pressure_secants[i] / projected_pressure_secants[i - 1] - 1.0)
        for i in range(1, len(projected_pressure_secants))
    ]
    source_secant_decrease_indices = [
        i + 1 for i in range(1, len(source_pressure_secants))
        if source_pressure_secants[i] < source_pressure_secants[i - 1] * (1.0 - 1.0e-12)
    ]
    projected_secant_decrease_indices = [
        i + 1 for i in range(1, len(projected_pressure_secants))
        if projected_pressure_secants[i] < projected_pressure_secants[i - 1] * (1.0 - 1.0e-12)
    ]

    source_q_violation_count = sum(
        b < a for a, b in zip(q_source, q_source[1:])
    )
    pooled_blocks = [b for b in blocks if b["width"] > 1]
    pooled_knots = sum(int(b["width"]) for b in pooled_blocks)
    pooled_reduction = sum(int(b["width"]) - 1 for b in pooled_blocks)

    def max_relative(actual: Sequence[float], expected: Sequence[float]) -> float:
        return max(
            abs(a - e) / max(abs(e), 1.0e-300)
            for a, e in zip(actual, expected)
        )

    source_first_law = max(
        abs(e + p - m * n) / max(abs(e + p), 1.0e-300)
        for m, n, p, e in zip(mu_v, n_v, p_v, eps_v)
    )
    projected_first_law = max(
        abs(e + p - m * n) / max(abs(e + p), 1.0e-300)
        for m, n, p, e in zip(mu_v, n_projected, p_projected, eps_projected)
    )

    row217, row218 = 216, 217
    has_r16_boundary = len(mu_v) > row218
    if has_r16_boundary:
        block217 = next(b for b in blocks if b["lo_zero_based"] <= row217 < b["hi_exclusive_zero_based"])
        block218 = next(b for b in blocks if b["lo_zero_based"] <= row218 < b["hi_exclusive_zero_based"])
    else:
        block217 = block218 = None

    def boundary_value(values: Sequence[float], index: int):
        return values[index] if has_r16_boundary else None

    return {
        "row_count": len(mu_v),
        "source_cell_count": len(source_cell_diagnostics),
        "source_nonstable_cell_count": source_nonstable_cells,
        "source_noncausal_cell_count": source_noncausal_cells,
        "projected_nonstable_cell_count": projected_nonstable_cells,
        "projected_noncausal_cell_count": projected_noncausal_cells,
        "source_min_d2p_all_cells": min(x["min_d2p"] for x in source_cell_diagnostics),
        "source_min_causal_margin_all_cells": min(x["min_causal_margin"] for x in source_cell_diagnostics),
        "projected_min_d2p_all_cells": min(x["min_d2p"] for x in projected_cell_diagnostics),
        "projected_min_causal_margin_all_cells": min(x["min_causal_margin"] for x in projected_cell_diagnostics),
        "source_row217_218_cell": source_boundary_cell,
        "projected_row217_218_cell": projected_boundary_cell,
        "source_cell_diagnostics": source_cell_diagnostics,
        "projected_cell_diagnostics": projected_cell_diagnostics,
        "source_q_monotone": source_q_violation_count == 0,
        "source_q_violation_count": source_q_violation_count,
        "source_pressure_secants_dpdmu": source_pressure_secants,
        "projected_pressure_secants_dpdmu": projected_pressure_secants,
        "source_secant_change_pct": source_secant_change_pct,
        "projected_secant_change_pct": projected_secant_change_pct,
        "source_secant_decrease_count": len(source_secant_decrease_indices),
        "source_secant_decrease_lower_cell_indices": source_secant_decrease_indices,
        "projected_secant_decrease_count": len(projected_secant_decrease_indices),
        "projected_secant_decrease_lower_cell_indices": projected_secant_decrease_indices,
        "source_secant_worst_relative_drop_pct": min(
            (x for x in source_secant_change_pct if x is not None), default=0.0
        ),
        "pav_block_count": len(blocks),
        "pav_pooled_block_count": len(pooled_blocks),
        "pav_pooled_knots": pooled_knots,
        "pav_knot_reduction": pooled_reduction,
        "pav_row217_block_width": block217["width"] if has_r16_boundary else None,
        "pav_row218_block_width": block218["width"] if has_r16_boundary else None,
        "max_abs_delta_n_over_n": max_relative(n_projected, n_v),
        "max_abs_delta_p_over_p": max_relative(p_projected, p_v),
        "max_abs_delta_p_over_pmax": max(abs(a - b) for a, b in zip(p_projected, p_v)) / max(p_v),
        "max_abs_delta_eps_over_eps": max_relative(eps_projected, eps_v),
        "max_abs_delta_eps_over_mu_n": max(
            abs(a - b) / max(abs(m * n), 1.0e-300)
            for a, b, m, n in zip(eps_projected, eps_v, mu_v, n_v)
        ),
        "source_first_law_max_rel": source_first_law,
        "projected_first_law_max_rel": projected_first_law,
        "row217_n_source": boundary_value(n_v, row217),
        "row218_n_source": boundary_value(n_v, row218),
        "row217_mu_source": boundary_value(mu_v, row217),
        "row218_mu_source": boundary_value(mu_v, row218),
        "row217_p_source": boundary_value(p_v, row217),
        "row218_p_source": boundary_value(p_v, row218),
        "row217_p_projected": boundary_value(p_projected, row217),
        "row218_p_projected": boundary_value(p_projected, row218),
        "row217_eps_source": boundary_value(eps_v, row217),
        "row218_eps_source": boundary_value(eps_v, row218),
        "row217_eps_projected": boundary_value(eps_projected, row217),
        "row218_eps_projected": boundary_value(eps_projected, row218),
        "row218_pressure_endpoint_delta_pct": (
            100.0 * (p_projected[row218] / p_v[row218] - 1.0)
            if has_r16_boundary else None
        ),
        "row218_energy_endpoint_delta_pct": (
            100.0 * (eps_projected[row218] / eps_v[row218] - 1.0)
            if has_r16_boundary else None
        ),
        "row217_218_source_delta_p": p_v[row218] - p_v[row217] if has_r16_boundary else None,
        "row217_218_projected_delta_p": p_projected[row218] - p_projected[row217] if has_r16_boundary else None,
        "row217_218_delta_p_increment_delta_pct": (
            100.0 * (
                (p_projected[row218] - p_projected[row217])
                / (p_v[row218] - p_v[row217]) - 1.0
            ) if has_r16_boundary else None
        ),
        "projected_n": n_projected,
        "projected_p": p_projected,
        "projected_eps": eps_projected,
        "projected_q": q_projected,
        "q_source": q_source,
        "q_projected": q_projected,
        "n_projected": n_projected,
        "p_projected": p_projected,
        "eps_projected": eps_projected,
        "pav_blocks": blocks,
    }
