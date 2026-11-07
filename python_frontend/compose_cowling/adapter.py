"""Thermodynamically consistent q-linear adapter for the native P(mu_B) Hermite solver.

The upstream EOS generator remains authoritative for source and joined EOS knots.
This module constructs a separate, explicitly versioned interpolation closure used
by the Cowling solver. It never overwrites the input table.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any
import numpy as np

ENERGY_FACTOR = 1.602176634e32
DENSITY_FACTOR = 1.0e45

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def pav_increasing(values: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, list[dict[str, float]]]:
    """Weighted pool-adjacent-violators projection onto nondecreasing sequences."""
    blocks: list[list[float]] = []
    for i, (value, weight) in enumerate(zip(values, weights)):
        if not (math.isfinite(float(value)) and math.isfinite(float(weight)) and weight > 0):
            raise ValueError("PAV requires finite values and positive finite weights")
        blocks.append([float(i), float(i + 1), float(weight), float(weight * value)])
        while len(blocks) >= 2:
            left, right = blocks[-2:]
            if left[3] / left[2] <= right[3] / right[2]:
                break
            blocks[-2:] = [[left[0], right[1], left[2] + right[2], left[3] + right[3]]]
    projected = np.empty_like(values, dtype=float)
    report: list[dict[str, float]] = []
    for lo, hi, weight, weighted_sum in blocks:
        projected[int(lo):int(hi)] = weighted_sum / weight
        report.append({"start": lo, "stop": hi, "weight": weight, "mean": weighted_sum / weight})
    return projected, report

def build_adapter(
    mu: np.ndarray,
    n_source: np.ndarray,
    p_source: np.ndarray,
    eps_source: np.ndarray,
    *,
    enforce_monotone_q: bool = True,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Construct n=mu*q, P=P0+integral(n dmu), eps=mu*n-P.

    q is linear in each chemical-potential interval. This makes n quadratic and P
    exactly cubic, matching the EOS solver's endpoint-Hermite P(mu_B) representation.
    """
    mu = np.asarray(mu, dtype=float)
    n0 = np.asarray(n_source, dtype=float)
    p0 = np.asarray(p_source, dtype=float)
    e0 = np.asarray(eps_source, dtype=float)
    arrays = (mu, n0, p0, e0)
    if any(a.ndim != 1 for a in arrays) or len(mu) < 3 or any(len(a) != len(mu) for a in arrays):
        raise ValueError("EOS columns must be equally sized 1-D arrays with at least three rows")
    if any(not np.all(np.isfinite(a)) for a in arrays):
        raise ValueError("nonfinite EOS input")
    if np.any(mu <= 0) or np.any(n0 <= 0) or np.any(p0 <= 0) or np.any(e0 <= 0):
        raise ValueError("mu, n, P and epsilon must be strictly positive")
    if np.any(np.diff(mu) <= 0) or np.any(np.diff(n0) <= 0) or np.any(np.diff(p0) <= 0) or np.any(np.diff(e0) <= 0):
        raise ValueError("mu, n, P and epsilon must be strictly increasing")
    identity = e0 + p0 - mu * n0
    identity_rel = np.abs(identity) / np.maximum(np.abs(e0 + p0), 1e-300)
    if float(np.max(identity_rel)) > 2e-12:
        raise ValueError(f"source knots fail the cold first-law identity: {np.max(identity_rel):.3e}")

    h = np.diff(mu)
    q0 = n0 / mu
    weights = np.empty_like(mu)
    weights[0] = 0.5 * h[0]
    weights[-1] = 0.5 * h[-1]
    weights[1:-1] = 0.5 * (h[:-1] + h[1:])
    tol = 64.0 * np.finfo(float).eps * max(float(np.max(np.abs(q0))), 1e-300)
    monotone_source = bool(np.all(np.diff(q0) >= -tol))
    if not monotone_source and not enforce_monotone_q:
        raise ValueError("q=n/mu is not monotone; explicit PAV repair is disabled")
    q, blocks = pav_increasing(q0, weights) if enforce_monotone_q else (q0.copy(), [])
    if np.any(q <= 0) or np.any(np.diff(q) < -tol):
        raise ValueError("positive monotone q could not be constructed")

    dq_dmu = np.diff(q) / h
    if np.any(dq_dmu < -tol):
        raise ValueError("negative q slope violates the causal closure")
    n = mu * q
    p = np.empty_like(mu)
    p[0] = p0[0]
    for i, width in enumerate(h):
        delta_q = q[i + 1] - q[i]
        # Exact integral of mu*q(mu), with q linear over this native interval.
        p[i + 1] = p[i] + width * (
            (mu[i] + 0.5 * width) * q[i]
            + (0.5 * mu[i] + width / 3.0) * delta_q
        )
    eps = mu * n - p
    if np.any(p <= 0) or np.any(eps <= 0) or np.any(np.diff(p) <= 0) or np.any(np.diff(eps) <= 0):
        raise ValueError("constructed EOS is not positive and strictly monotone")

    # Exact thermodynamic and causality checks for the q-linear closure.
    cs2_left = q[:-1] / (q[:-1] + mu[:-1] * dq_dmu)
    cs2_right = q[1:] / (q[1:] + mu[1:] * dq_dmu)
    if not np.all(np.isfinite(np.r_[cs2_left, cs2_right])) or np.any(cs2_left <= 0) or np.any(cs2_right <= 0):
        raise ValueError("nonpositive or nonfinite sound speed")
    if np.any(cs2_left > 1.0 + 128 * np.finfo(float).eps) or np.any(cs2_right > 1.0 + 128 * np.finfo(float).eps):
        raise ValueError("causal bound violated")

    # Verify that the exact q-linear pressure equals the cubic Hermite P(mu) at
    # interior points. The derivative of this polynomial is n=mu*q.
    t = np.array([0.0, 0.125, 0.25, 0.5, 0.75, 0.875, 1.0])[None, :]
    h2 = h[:, None]
    h00 = 2*t**3 - 3*t**2 + 1
    h10 = t**3 - 2*t**2 + t
    h01 = -2*t**3 + 3*t**2
    h11 = t**3 - t**2
    ph = h00*p[:-1, None] + h10*h2*n[:-1, None] + h01*p[1:, None] + h11*h2*n[1:, None]
    # Evaluate exact integral at the same local coordinates.
    mu_t = mu[:-1, None] + h2*t
    q_t = q[:-1, None] + (q[1:] - q[:-1])[:, None]*t
    p_exact = p[:-1, None] + h2 * (
        (mu[:-1, None] + 0.5*h2*t) * q[:-1, None] * t
        + (0.5*mu[:-1, None] + h2*t/3.0) * (q[1:] - q[:-1])[:, None] * t**2
    )
    # The expression above is the exact integral over the partial interval.
    hermite_rel = float(np.max(np.abs(ph - p_exact) / np.maximum(np.abs(p_exact), 1e-300)))
    if hermite_rel > 2e-12:
        raise ValueError(f"q-linear pressure is not represented by the native cubic Hermite polynomial: {hermite_rel:.3e}")

    eps_closure = eps + p - mu*n
    metrics: dict[str, Any] = {
        "n_rows": int(len(mu)),
        "source_first_law_max_relative": float(np.max(identity_rel)),
        "adapter_first_law_max_relative": float(np.max(np.abs(eps_closure) / np.maximum(np.abs(eps+p), 1e-300))),
        "source_q_monotone": monotone_source,
        "pav_pooled_blocks": int(sum(1 for b in blocks if b["stop"] - b["start"] > 1)),
        "pav_pooled_knots": int(sum(int(b["stop"] - b["start"]) for b in blocks if b["stop"] - b["start"] > 1)),
        "max_abs_delta_n_over_n_source": float(np.max(np.abs(n - n0) / n0)),
        "rms_delta_n_over_n_source": float(np.sqrt(np.mean(((n - n0) / n0)**2))),
        "max_abs_delta_P_over_Pmax": float(np.max(np.abs(p - p0)) / np.max(p0)),
        "median_abs_delta_P_over_Pmax": float(np.median(np.abs(p - p0)) / np.max(p0)),
        "max_abs_delta_epsilon_over_mu_n": float(np.max(np.abs(eps - e0) / (mu*n0))),
        "min_cs2": float(min(np.min(cs2_left), np.min(cs2_right))),
        "max_cs2": float(max(np.max(cs2_left), np.max(cs2_right))),
        "min_q_slope": float(np.min(dq_dmu)),
        "max_q_slope": float(np.max(dq_dmu)),
        "max_relative_hermite_pressure_residual": hermite_rel,
    }
    return {
        "mu_b_mev": mu,
        "n_b_fm3": n,
        "pressure_mev_fm3": p,
        "energy_mev_fm3": eps,
        "q_n_over_mu": q,
        "cs2_cell_left": cs2_left,
        "cs2_cell_right": cs2_right,
    }, metrics

def build_adapter_csv(source: Path, destination: Path, *, enforce_monotone_q: bool = True) -> dict[str, Any]:
    source, destination = Path(source), Path(destination)
    raw = np.genfromtxt(source, delimiter=",", names=True, dtype=None, encoding="utf-8")
    mu = np.asarray(raw["mu_b_mev"], dtype=float)
    n0 = np.asarray(raw["n_b_fm3"], dtype=float)
    p0 = np.asarray(raw["pressure_mev_fm3"], dtype=float)
    e0 = np.asarray(raw["energy_mev_fm3"], dtype=float)
    adapter, metrics = build_adapter(mu, n0, p0, e0, enforce_monotone_q=enforce_monotone_q)
    columns = list(raw.dtype.names)
    values = {name: np.asarray(raw[name]).copy() for name in columns}
    values["n_b_fm3"] = adapter["n_b_fm3"]
    values["n_b_m3"] = adapter["n_b_fm3"] * DENSITY_FACTOR
    values["pressure_mev_fm3"] = adapter["pressure_mev_fm3"]
    values["pressure_pa"] = adapter["pressure_mev_fm3"] * ENERGY_FACTOR
    values["energy_mev_fm3"] = adapter["energy_mev_fm3"]
    values["energy_j_m3"] = adapter["energy_mev_fm3"] * ENERGY_FACTOR
    values["q1_p_over_nb_mev"] = adapter["pressure_mev_fm3"] / adapter["n_b_fm3"]
    # The physical nucleon mass may differ by source model; infer it from the
    # source identity mu_B = m_n(1+Q3) row by row, then retain source mass convention.
    mn_values = mu / (1.0 + np.asarray(raw["q3_mu_b_over_mn_minus_1"], dtype=float))
    mn = float(np.median(mn_values))
    if not np.allclose(mn_values, mn, rtol=2e-11, atol=0):
        raise ValueError("source Q3 does not encode one consistent neutron mass")
    values["q3_mu_b_over_mn_minus_1"] = mu / mn - 1.0
    values["q7_e_over_nb_mn_minus_1"] = adapter["energy_mev_fm3"] / (adapter["n_b_fm3"] * mn) - 1.0
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        for i in range(len(mu)):
            writer.writerow([format(float(values[k][i]), ".17g") if np.issubdtype(np.asarray(values[k]).dtype, np.number) else values[k][i] for k in columns])
    receipt = {
        "schema": "compose-cowling-q-linear-adapter-v1",
        "adapter_kind": "solver_interpolation_closure_not_source_eos_replacement",
        "source_path": str(source),
        "source_sha256": sha256(source),
        "adapter_path": str(destination),
        "adapter_sha256": sha256(destination),
        "neutron_mass_mev": mn,
        "construction": "weighted_PAV_q_if_needed; piecewise_linear_q(mu_B); n=mu_B*q; P=P0+integral(n dmu_B); epsilon=mu_B*n-P",
        "metrics": metrics,
    }
    sidecar = destination.with_suffix(".adapter.json")
    sidecar.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt
