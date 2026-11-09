"""Analytic pressure-preserving, causal q=n/mu cell interpolation family."""
from __future__ import annotations

import math


def cell_pressure_increment_bounds(
    mu0: float, mu1: float, q0: float, q1: float
) -> tuple[float, float]:
    """Bounds on integral(mu*q dmu) for any nondecreasing q joining q0 to q1."""
    mu0, mu1, q0, q1 = map(float, (mu0, mu1, q0, q1))
    if not all(math.isfinite(x) for x in (mu0, mu1, q0, q1)):
        raise ValueError("cell endpoints must be finite")
    if mu0 <= 0.0 or mu1 <= mu0 or q0 <= 0.0 or q1 <= 0.0:
        raise ValueError("mu and q endpoints must be positive with increasing mu")
    if q1 < q0:
        raise ValueError("q endpoints are decreasing; no nondecreasing q interpolant exists")
    area = 0.5 * (mu1 * mu1 - mu0 * mu0)
    return q0 * area, q1 * area


def _shape_integrals(k: float) -> tuple[float, float]:
    """A0=integral S_k(t)dt and A1=integral t*S_k(t)dt on [0,1]."""
    k = float(k)
    if not math.isfinite(k):
        raise ValueError("shape parameter must be finite")
    if abs(k) < 1.0e-2:
        k2 = k * k
        k3 = k2 * k
        k4 = k2 * k2
        k5 = k4 * k
        k6 = k3 * k3
        k7 = k6 * k
        k8 = k4 * k4
        k9 = k8 * k
        # Taylor expansions avoid cancellation in the closed forms near k=0.
        a0 = (
            0.5 + k / 12.0 - k3 / 720.0 + k5 / 30240.0
            - k7 / 1209600.0 + k9 / 47900160.0
        )
        a1 = (
            1.0 / 3.0 + k / 24.0 - k2 / 720.0
            - k3 / 1440.0 + k4 / 30240.0 + k5 / 60480.0
            - k6 / 1209600.0 - k7 / 2419200.0
            + k8 / 47900160.0 + k9 / 95800320.0
        )
        return a0, a1

    if k > 0.0:
        denominator = -math.expm1(-k)
        exp_minus = math.exp(-k) if k < 745.0 else 0.0
        a0 = 1.0 / denominator - 1.0 / k
        a1 = (
            0.5 - (1.0 - (1.0 + k) * exp_minus) / (k * k)
        ) / denominator
        return a0, a1

    a = -k
    if a > 500.0:
        # The omitted exp(-a) terms are below double-precision relevance.
        a0 = 1.0 / a
        a1 = 1.0 / a - 1.0 / (a * a)
        return a0, a1

    denominator = math.expm1(a)
    exp_a = math.exp(a)
    a0 = 1.0 / a - 1.0 / denominator
    a1 = (
        (exp_a * (a - 1.0) + 1.0) / (a * a) - 0.5
    ) / denominator
    return a0, a1


def qexp_pressure_increment(
    mu0: float, mu1: float, q0: float, q1: float, k: float
) -> float:
    """Exact integral for q(t)=q0+(q1-q0)*(1-exp(-k*t))/(1-exp(-k))."""
    mu0, mu1, q0, q1 = map(float, (mu0, mu1, q0, q1))
    if not (mu0 > 0.0 and mu1 > mu0 and q0 > 0.0 and q1 >= q0):
        raise ValueError("invalid increasing-mu, nondecreasing-q cell")
    h = mu1 - mu0
    dq = q1 - q0
    a0, a1 = _shape_integrals(k)
    return h * (q0 * (mu0 + 0.5 * h) + dq * (mu0 * a0 + h * a1))


def qexp_shape_derivatives(k: float) -> tuple[float, float]:
    """Return dS/dt at t=0 and t=1; both are positive for finite k."""
    k = float(k)
    if not math.isfinite(k):
        raise ValueError("shape parameter must be finite")
    if abs(k) < 1.0e-4:
        k2 = k * k
        k4 = k2 * k2
        return (
            1.0 + 0.5 * k + k2 / 12.0 - k4 / 720.0,
            1.0 - 0.5 * k + k2 / 12.0 - k4 / 720.0,
        )
    if k > 0.0:
        d0 = k / (-math.expm1(-k))
        d1 = k / math.expm1(k) if k < 700.0 else 0.0
        return d0, d1
    a = -k
    d0 = a / math.expm1(a) if a < 700.0 else 0.0
    d1 = a / (-math.expm1(-a))
    return d0, d1


def qexp_shape_parameter(
    mu0: float,
    mu1: float,
    q0: float,
    q1: float,
    source_pressure_increment: float,
    *,
    regularization_fraction: float = 0.0,
) -> dict:
    """Find k matching a source pressure increment, or its nearest regularized bound.

    q is nondecreasing iff q1>=q0. For any continuous monotone q between
    endpoint values, the pressure increment lies strictly between
    q0*integral(mu dmu) and q1*integral(mu dmu), unless q0==q1.
    A positive regularization_fraction keeps the target mean away from those
    limiting step-function bounds and therefore gives finite, smooth endpoint
    slopes for the exponential shape family.
    """
    mu0, mu1, q0, q1, source_dp = map(
        float, (mu0, mu1, q0, q1, source_pressure_increment)
    )
    if not all(math.isfinite(x) for x in (mu0, mu1, q0, q1, source_dp)):
        raise ValueError("cell values must be finite")
    if mu0 <= 0.0 or mu1 <= mu0 or q0 <= 0.0 or q1 <= 0.0 or source_dp <= 0.0:
        raise ValueError("mu, q and source pressure increment must be positive")
    if not (0.0 <= regularization_fraction < 0.5):
        raise ValueError("regularization_fraction must lie in [0, 0.5)")

    lower, upper = cell_pressure_increment_bounds(mu0, mu1, q0, q1)
    area = 0.5 * (mu1 * mu1 - mu0 * mu0)
    dq = q1 - q0
    tol = 1.0e-12 * max(upper - lower, abs(source_dp), 1.0e-300)

    if dq == 0.0:
        exact = abs(source_dp - lower) <= tol
        return {
            "status": "PASS_CONSTANT_Q" if exact else "INFEASIBLE_CONSTANT_Q",
            "feasible": exact,
            "realizable_by_continuous_q": exact,
            "k": 0.0 if exact else None,
            "lower_bound": lower,
            "upper_bound": upper,
            "source_pressure_increment": source_dp,
            "target_pressure_increment": lower,
            "source_mean_q": source_dp / area,
            "target_mean_q": q0,
            "regularization_applied": False,
            "qprime_lower": 0.0 if exact else None,
            "qprime_upper": 0.0 if exact else None,
            "chi_lower": 1.0 if exact else None,
            "chi_upper": 1.0 if exact else None,
        }

    source_mean_q = source_dp / area
    # Strict feasibility for a continuous q curve with q0<q1.
    feasible = (source_dp > lower + tol) and (source_dp < upper - tol)
    if source_dp <= lower + tol:
        status = "BELOW_MONOTONE_Q_ENVELOPE"
    elif source_dp >= upper - tol:
        status = "ABOVE_MONOTONE_Q_ENVELOPE"
    else:
        status = "SOURCE_INCREMENT_FEASIBLE"

    # At either envelope boundary the admissible monotone q curve is only
    # approached by an endpoint step. No finite k in this smooth family reaches
    # that limit. Do not silently pretend that a finite shape parameter exists.
    if not feasible and regularization_fraction == 0.0:
        target_dp = lower if source_dp <= lower + tol else upper
        return {
            "status": status,
            "feasible": False,
            "realizable_by_continuous_q": False,
            "k": None,
            "limiting_shape": "endpoint_step_limit",
            "lower_bound": lower,
            "upper_bound": upper,
            "source_pressure_increment": source_dp,
            "target_pressure_increment": target_dp,
            "source_mean_q": source_mean_q,
            "target_mean_q": target_dp / area,
            "regularization_applied": False,
            "qprime_lower": None,
            "qprime_upper": None,
            "chi_lower": None,
            "chi_upper": None,
        }

    margin = regularization_fraction * dq
    target_mean_q = min(max(source_mean_q, q0 + margin), q1 - margin)
    target_dp = area * target_mean_q
    regularized = abs(target_dp - source_dp) > tol

    linear_dp = qexp_pressure_increment(mu0, mu1, q0, q1, 0.0)
    if abs(target_dp - linear_dp) <= tol:
        k = 0.0
    else:
        lo, hi = (-1.0, 0.0) if target_dp < linear_dp else (0.0, 1.0)
        # q-exponential approaches its step-function envelope only as |k|→∞.
        # The analytic integral uses asymptotic forms beyond |k|=500, so allow
        # a wider bracket for targets close to the regularized envelope.
        shape_limit = 1.0e6
        while qexp_pressure_increment(mu0, mu1, q0, q1, lo) > target_dp and lo > -shape_limit:
            lo = max(-shape_limit, 2.0 * lo)
        while qexp_pressure_increment(mu0, mu1, q0, q1, hi) < target_dp and hi < shape_limit:
            hi = min(shape_limit, 2.0 * hi)
        if (
            qexp_pressure_increment(mu0, mu1, q0, q1, lo) > target_dp + tol
            or qexp_pressure_increment(mu0, mu1, q0, q1, hi) < target_dp - tol
        ):
            raise RuntimeError("could not bracket finite q-exponential shape parameter")
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            value = qexp_pressure_increment(mu0, mu1, q0, q1, mid)
            if value < target_dp:
                lo = mid
            else:
                hi = mid
        k = 0.5 * (lo + hi)

    slope_t0, slope_t1 = qexp_shape_derivatives(k)
    h = mu1 - mu0
    qprime0 = dq * slope_t0 / h
    qprime1 = dq * slope_t1 / h
    chi0 = q0 / (q0 + mu0 * qprime0)
    chi1 = q1 / (q1 + mu1 * qprime1)
    return {
        "status": status,
        "feasible": feasible,
        "realizable_by_continuous_q": True,
        "k": k,
        "lower_bound": lower,
        "upper_bound": upper,
        "source_pressure_increment": source_dp,
        "target_pressure_increment": target_dp,
        "source_mean_q": source_mean_q,
        "target_mean_q": target_mean_q,
        "regularization_applied": regularized,
        "linear_q_pressure_increment": linear_dp,
        "qprime_lower": qprime0,
        "qprime_upper": qprime1,
        "chi_lower": chi0,
        "chi_upper": chi1,
        "whole_cell_stable": q0 > 0.0 and q1 > 0.0 and qprime0 >= 0.0 and qprime1 >= 0.0,
        "whole_cell_causal": q0 > 0.0 and q1 > 0.0 and qprime0 >= 0.0 and qprime1 >= 0.0,
    }
