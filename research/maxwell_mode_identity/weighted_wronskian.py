"""Weighted Abel-Liouville concomitant for the relativistic Cowling W,V system.

This module encodes an algebraic determinant identity. It does not claim that
this reduced quantity is the canonical covariant symplectic current.
"""
from __future__ import annotations
import math
import numpy as np


def _state(value: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(value)
    if arr.shape[-1:] != (2,):
        raise ValueError(f"{name} must have final dimension 2 in [W,V] order")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def weighted_wronskian(state_1, state_2, enthalpy_j_m3, nu):
    """Return h exp(-nu) det([state_1,state_2]) without complex conjugation."""
    y1, y2 = _state(state_1, "state_1"), _state(state_2, "state_2")
    h, metric_nu = np.asarray(enthalpy_j_m3), np.asarray(nu)
    if np.any(~np.isfinite(h)) or np.any(h <= 0.0):
        raise ValueError("enthalpy_j_m3 must be finite and strictly positive")
    if np.any(~np.isfinite(metric_nu)):
        raise ValueError("nu must be finite")
    determinant = y1[..., 0] * y2[..., 1] - y1[..., 1] * y2[..., 0]
    result = h * np.exp(-metric_nu) * determinant
    return result.item() if np.ndim(result) == 0 else result


def interface_invariant_residual(
    state_1_inner, state_2_inner, state_1_outer, state_2_outer, *,
    enthalpy_inner_j_m3, enthalpy_outer_j_m3, nu_inner, nu_outer,
    relative_floor=1e-300, determinant_transfer=None
):
    """Compare one-sided concomitants and return the general jump diagnostics.

    If the inner concomitant is zero or below relative_floor, the observed
    ratio is mathematically undefined. Return None plus an explicit status
    instead of a NaN, because NaN is not valid strict JSON.
    """
    omega_i = complex(weighted_wronskian(
        state_1_inner, state_2_inner, enthalpy_inner_j_m3, nu_inner))
    omega_o = complex(weighted_wronskian(
        state_1_outer, state_2_outer, enthalpy_outer_j_m3, nu_outer))
    floor = float(relative_floor)
    if not math.isfinite(floor) or floor <= 0.0:
        raise ValueError("relative_floor must be finite and positive")
    scale = max(abs(omega_i), abs(omega_o), floor)
    h_ratio = enthalpy_inner_j_m3 / enthalpy_outer_j_m3
    nu_jump = nu_outer - nu_inner
    if abs(omega_i) > floor:
        observed = omega_o / omega_i
        observed_status = "DEFINED"
    else:
        observed = None
        observed_status = "UNDEFINED_ZERO_DENOMINATOR"
    predicted = None
    if determinant_transfer is not None:
        predicted = (1.0 / h_ratio) * math.exp(-nu_jump) * determinant_transfer
    return {
        "omega_inner": omega_i,
        "omega_outer": omega_o,
        "absolute_jump": abs(omega_o - omega_i),
        "relative_jump": abs(omega_o - omega_i) / scale,
        "enthalpy_ratio_inner_over_outer": h_ratio,
        "metric_nu_jump_outer_minus_inner": nu_jump,
        "determinant_transfer": determinant_transfer,
        "predicted_omega_ratio": predicted,
        "observed_omega_ratio": observed,
        "observed_omega_ratio_status": observed_status,
    }


def surface_boundary_basis_wv(
    radius_m, nu_prime_per_m, omega_rad_s, nu, lambda_, *,
    c_m_s=299_792_458.0
):
    """Return one [W,V] vector spanning the current ell=2 finite-cutoff surface line.

    The current boundary form is b1*(W/R^3)+b2*(V/R^2)=0, with
    b1=R*nu_prime and b2=(omega*R/c)^2 exp(lambda-2nu).
    """
    vals = (radius_m, nu_prime_per_m, omega_rad_s, nu, lambda_, c_m_s)
    if not all(math.isfinite(float(x)) for x in vals):
        raise ValueError("surface boundary inputs must be finite")
    r, nu_p, omega, metric_nu, metric_lambda, c = map(float, vals)
    if r <= 0.0 or omega <= 0.0 or c <= 0.0 or nu_p <= 0.0:
        raise ValueError("radius, frequency, c and nu_prime must be positive")
    b1 = r * nu_p
    b2 = (omega * r / c) ** 2 * math.exp(metric_lambda - 2.0 * metric_nu)
    if not math.isfinite(b2):
        raise ArithmeticError("surface frequency coefficient overflowed")
    return np.array([-b2 / r**2, b1 / r**3], dtype=float)


def two_sided_matching_concomitant(center_state_wv, surface_state_wv,
                                    enthalpy_j_m3, nu):
    """Raw weighted two-sided matching determinant at one matching radius.

    This is the Abel-Liouville concomitant. Do not normalize it by Euclidean
    vector norms when using complex frequency or claiming analytic dependence.
    """
    return weighted_wronskian(center_state_wv, surface_state_wv,
                              enthalpy_j_m3, nu)


def normalized_matching_determinant(center_state_wv, surface_line_state_wv,
                                     radius_m):
    """Real-frequency determinant angle for z=(W/r^3,V/r^2); diagnostic only."""
    yc = _state(center_state_wv, "center_state_wv")
    ys = _state(surface_line_state_wv, "surface_line_state_wv")
    r = float(radius_m)
    if not math.isfinite(r) or r <= 0.0:
        raise ValueError("radius_m must be finite and positive")
    if yc.ndim != 1 or ys.ndim != 1:
        raise ValueError("matching determinant expects two single [W,V] vectors")
    zc = np.array([yc[0]/r**3, yc[1]/r**2], dtype=np.result_type(yc, ys, float))
    zs = np.array([ys[0]/r**3, ys[1]/r**2], dtype=np.result_type(yc, ys, float))
    denom = float(np.linalg.norm(zc) * np.linalg.norm(zs))
    if denom <= np.finfo(float).tiny:
        raise ArithmeticError("matching state is numerically zero")
    value = (zc[0]*zs[1] - zc[1]*zs[0]) / denom
    return value.item() if hasattr(value, "item") else value
