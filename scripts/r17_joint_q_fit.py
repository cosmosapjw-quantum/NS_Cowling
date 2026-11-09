"""Constrained joint pressure-density reconstruction sensitivity pilot.

Lambda values are sensitivity parameters, not inferred statistical precisions.
"""
from __future__ import annotations

import math
from typing import Any, Sequence

import numpy as np


def _isotonic_projection(values: np.ndarray, lower_bound: float) -> np.ndarray:
    """Euclidean projection onto nondecreasing sequences with a positive floor."""
    blocks: list[list[float | int]] = []
    for i, value in enumerate(np.asarray(values, dtype=float)):
        blocks.append([i, i + 1, 1.0, float(value)])
        while len(blocks) >= 2:
            a, b = blocks[-2:]
            if float(a[3]) / float(a[2]) <= float(b[3]) / float(b[2]):
                break
            blocks[-2:] = [[
                int(a[0]), int(b[1]),
                float(a[2]) + float(b[2]),
                float(a[3]) + float(b[3]),
            ]]
    result = np.empty(len(values), dtype=float)
    for lo, hi, weight, weighted_sum in blocks:
        result[int(lo):int(hi)] = float(weighted_sum) / float(weight)
    return np.maximum(result, lower_bound)


def pressure_operator(mu: Sequence[float]) -> np.ndarray:
    """Return A such that P(q)=P0+Aq for cellwise-linear q and n=mu*q."""
    mu = np.asarray(mu, dtype=float)
    if mu.ndim != 1 or len(mu) < 3 or not np.all(np.isfinite(mu)) or np.any(mu <= 0.0):
        raise ValueError("mu must be a finite positive one-dimensional array with >=3 knots")
    h = np.diff(mu)
    if np.any(h <= 0.0):
        raise ValueError("mu must be strictly increasing")
    operator = np.zeros((len(mu), len(mu)), dtype=float)
    for i, step in enumerate(h):
        c0 = step * (0.5 * mu[i] + step / 6.0)
        c1 = step * (0.5 * mu[i] + step / 3.0)
        operator[i + 1:, i] += c0
        operator[i + 1:, i + 1] += c1
    return operator


def fit_monotone_q(
    mu: Sequence[float],
    n_source: Sequence[float],
    p_source: Sequence[float],
    *,
    lambda_density: float,
    weight_coordinate: str = "mu",
    max_iterations: int = 50000,
    tolerance: float = 1.0e-12,
) -> dict[str, Any]:
    """Solve a convex q-fit under positivity and nondecreasing-q constraints.

    Minimize
      lambda * sum(w_i * ((mu_i*q_i)/n_i - 1)^2)
      + (1-lambda) * sum(w_i * ((P_i(q)-P_i_source)/P_i_source)^2)
    subject to q_i >= q_floor and q[i+1] >= q[i].

    Trapezoidal chemical-potential weights are normalized to sum to one.
    Lambda is a sensitivity parameter, not a statistical confidence weight.
    """
    mu_v = np.asarray(mu, dtype=float)
    n_v = np.asarray(n_source, dtype=float)
    p_v = np.asarray(p_source, dtype=float)
    if not (mu_v.ndim == n_v.ndim == p_v.ndim == 1 and len(mu_v) == len(n_v) == len(p_v) and len(mu_v) >= 3):
        raise ValueError("mu, n and pressure must be equal-length one-dimensional arrays with >=3 knots")
    if not all(np.all(np.isfinite(x)) for x in (mu_v, n_v, p_v)):
        raise ValueError("source arrays must be finite")
    if not all(np.all(x > 0.0) for x in (mu_v, n_v, p_v)):
        raise ValueError("source arrays must be positive")
    if any(np.any(np.diff(x) <= 0.0) for x in (mu_v, n_v, p_v)):
        raise ValueError("mu, density and source pressure must be strictly increasing")

    lam = float(lambda_density)
    if not math.isfinite(lam) or not (0.0 < lam <= 1.0):
        raise ValueError("lambda_density must lie in (0, 1]")
    if max_iterations < 1 or not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("iteration limit and tolerance must be positive")

    h = np.diff(mu_v)
    q_source = n_v / mu_v
    q_scale = float(np.median(q_source))
    if not math.isfinite(q_scale) or q_scale <= 0.0:
        raise ValueError("source q scale is not positive")

    if weight_coordinate == "mu":
        quadrature_coordinate = mu_v
    elif weight_coordinate == "log_n":
        quadrature_coordinate = np.log(n_v)
    else:
        raise ValueError("weight_coordinate must be 'mu' or 'log_n'")
    quadrature_steps = np.diff(quadrature_coordinate)
    if np.any(quadrature_steps <= 0.0):
        raise ValueError("objective quadrature coordinate must be strictly increasing")
    weights = np.empty(len(mu_v), dtype=float)
    weights[0], weights[-1] = 0.5 * quadrature_steps[0], 0.5 * quadrature_steps[-1]
    weights[1:-1] = 0.5 * (quadrature_steps[:-1] + quadrature_steps[1:])
    weights /= np.sum(weights)

    operator = pressure_operator(mu_v)
    density_operator = mu_v * q_scale / n_v
    pressure_scale = np.maximum(p_v, p_v[0])
    pressure_operator_scaled = operator * q_scale / pressure_scale[:, None]
    pressure_target = (p_v - p_v[0]) / pressure_scale

    hessian = lam * np.diag(weights * density_operator**2)
    hessian += (1.0 - lam) * (
        pressure_operator_scaled.T @ (weights[:, None] * pressure_operator_scaled)
    )
    linear = lam * weights * density_operator
    linear += (1.0 - lam) * (
        pressure_operator_scaled.T @ (weights * pressure_target)
    )
    # grad J = 2(Hx-b), so L is twice the largest eigenvalue of H.
    lipschitz = 2.0 * float(np.linalg.eigvalsh(hessian)[-1])
    if not math.isfinite(lipschitz) or lipschitz <= 0.0:
        raise RuntimeError("joint-fit Hessian has invalid Lipschitz bound")

    q_floor = max(float(np.min(q_source)) * 1.0e-10, np.finfo(float).tiny)
    x_floor = q_floor / q_scale

    def project(x: np.ndarray) -> np.ndarray:
        return _isotonic_projection(x, x_floor)

    def components(x: np.ndarray) -> tuple[float, float]:
        density_residual = density_operator * x - 1.0
        pressure_residual = pressure_operator_scaled @ x - pressure_target
        return (
            float(np.dot(weights, density_residual**2)),
            float(np.dot(weights, pressure_residual**2)),
        )

    def objective(x: np.ndarray) -> float:
        jn, jp = components(x)
        return lam * jn + (1.0 - lam) * jp

    def gradient(x: np.ndarray) -> np.ndarray:
        return 2.0 * (hessian @ x - linear)

    x_source = q_source / q_scale
    x = project(x_source)
    y = x.copy()
    acceleration = 1.0
    kkt_residual = math.inf
    for iteration in range(1, max_iterations + 1):
        x_next = project(y - gradient(y) / lipschitz)
        acceleration_next = 0.5 * (1.0 + math.sqrt(1.0 + 4.0 * acceleration**2))
        y_next = x_next + ((acceleration - 1.0) / acceleration_next) * (x_next - x)
        kkt_residual = float(np.linalg.norm(
            x_next - project(x_next - gradient(x_next) / lipschitz), ord=np.inf
        ))
        if kkt_residual <= tolerance * (1.0 + float(np.linalg.norm(x_next, ord=np.inf))):
            x = x_next
            break
        x, y, acceleration = x_next, y_next, acceleration_next
    else:
        raise RuntimeError(
            f"joint q fit did not converge in {max_iterations} iterations; KKT={kkt_residual:.3e}"
        )

    q_fit = q_scale * x
    n_fit = mu_v * q_fit
    p_fit = np.empty_like(mu_v)
    p_fit[0] = p_v[0]
    for i, step in enumerate(h):
        dq = q_fit[i + 1] - q_fit[i]
        p_fit[i + 1] = p_fit[i] + step * (
            (mu_v[i] + 0.5 * step) * q_fit[i]
            + (0.5 * mu_v[i] + step / 3.0) * dq
        )
    eps_fit = mu_v * n_fit - p_fit
    jn, jp = components(x)
    source_jn, source_jp = components(project(x_source))

    density_relative = n_fit / n_v - 1.0
    pressure_relative = p_fit / p_v - 1.0
    q_monotone = bool(np.all(np.diff(q_fit) >= -1.0e-13 * q_scale))
    positive = bool(
        np.all(q_fit > 0.0) and np.all(n_fit > 0.0)
        and np.all(p_fit > 0.0) and np.all(eps_fit > 0.0)
    )
    pressure_increasing = bool(np.all(np.diff(p_fit) > 0.0))
    energy_increasing = bool(np.all(np.diff(eps_fit) > 0.0))

    return {
        "lambda_density": lam,
        "weight_coordinate": weight_coordinate,
        "objective_before_source_q": lam * source_jn + (1.0 - lam) * source_jp,
        "objective_after_fit": lam * jn + (1.0 - lam) * jp,
        "density_objective_before": source_jn,
        "pressure_objective_before": source_jp,
        "density_objective_after": jn,
        "pressure_objective_after": jp,
        "objective_reduction_fraction": (
            (lam * source_jn + (1.0 - lam) * source_jp
             - lam * jn - (1.0 - lam) * jp)
            / max(lam * source_jn + (1.0 - lam) * source_jp, 1.0e-300)
        ),
        "iterations": iteration,
        "projected_gradient_kkt_residual": kkt_residual,
        "q_monotone": q_monotone,
        "positive_thermodynamics": positive,
        "pressure_strictly_increasing": pressure_increasing,
        "energy_density_strictly_increasing": energy_increasing,
        "density_rms_relative_change": float(np.sqrt(np.mean(density_relative**2))),
        "density_max_abs_relative_change": float(np.max(np.abs(density_relative))),
        "pressure_rms_relative_change": float(np.sqrt(np.mean(pressure_relative**2))),
        "pressure_max_abs_relative_change": float(np.max(np.abs(pressure_relative))),
        "energy_rms_relative_change": float(
            np.sqrt(np.mean((eps_fit / (mu_v * n_v) - 1.0)**2))
        ),
        "source_knot_count": int(len(mu_v)),
        "n_source": n_v.tolist(),
        "n_fitted": n_fit.tolist(),
        "p_source": p_v.tolist(),
        "p_fitted": p_fit.tolist(),
        "epsilon_fitted": eps_fit.tolist(),
        "mu_b": mu_v.tolist(),
        "q_source": q_source.tolist(),
        "q_fitted": q_fit.tolist(),
    }
