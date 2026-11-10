"""Symbolically derived sharp-Maxwell Cowling interface transfer maps.

Convention: Q=A V+nu_prime W and DeltaP=-(h exp(-lambda)/r_t**2) Q.
W has units m^3, V m^2, A is dimensionless, nu_prime has units m^-1.
Maps send the high-density inner state [W_i,V_i] to the low-density outer
state [W_o,V_o]. This helper does not use EOS composition data.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np


@dataclass(frozen=True)
class InterfaceMap:
    conversion: str
    alpha: float
    A: float
    nu_prime: float
    matrix: np.ndarray
    determinant: float


def _validate(alpha: float, A: float, nu_prime: float) -> tuple[float, float, float]:
    vals = tuple(float(x) for x in (alpha, A, nu_prime))
    if not all(math.isfinite(x) for x in vals):
        raise ValueError("alpha, A and nu_prime must be finite")
    alpha, A, nu_prime = vals
    if alpha <= 0.0:
        raise ValueError("alpha=h_inner/h_outer must be positive")
    if A <= 0.0:
        raise ValueError("A must be positive for a positive-frequency interface shot")
    if nu_prime <= 0.0:
        raise ValueError("nu_prime must be positive for a regular TOV interface")
    return alpha, A, nu_prime


def interface_map(conversion: str, alpha: float, A: float, nu_prime: float) -> InterfaceMap:
    """Return exact slow/rapid transfer matrices.

    alpha=h_inner/h_outer. At cold Maxwell coexistence h=epsilon+P=mu_B*n_B,
    so alpha=n_inner/n_outer. Slow conversion imposes W_o=W_i; rapid
    conversion imposes continuity of W-Q/nu_prime; both impose [DeltaP]=0.
    """
    alpha, A, nu_prime = _validate(alpha, A, nu_prime)
    if conversion == "slow":
        matrix = np.array([[1.0, 0.0],
                           [(alpha - 1.0) * nu_prime / A, alpha]], dtype=float)
    elif conversion == "rapid":
        matrix = np.array([[alpha, (alpha - 1.0) * A / nu_prime],
                           [0.0, 1.0]], dtype=float)
    else:
        raise ValueError("conversion must be 'slow' or 'rapid'")
    matrix.setflags(write=False)
    return InterfaceMap(conversion, alpha, A, nu_prime, matrix, float(np.linalg.det(matrix)))


def apply_interface_map(conversion: str, alpha: float, A: float, nu_prime: float,
                        state_inner: np.ndarray) -> np.ndarray:
    state = np.asarray(state_inner, dtype=float)
    if state.shape != (2,) or not np.all(np.isfinite(state)):
        raise ValueError("state_inner must be a finite length-2 [W,V] vector")
    return interface_map(conversion, alpha, A, nu_prime).matrix @ state


def interface_residuals(conversion: str, alpha: float, A: float, nu_prime: float,
                        state_inner: np.ndarray, state_outer: np.ndarray) -> dict[str, float]:
    alpha, A, nu_prime = _validate(alpha, A, nu_prime)
    yi, yo = np.asarray(state_inner, dtype=float), np.asarray(state_outer, dtype=float)
    if yi.shape != (2,) or yo.shape != (2,) or not np.all(np.isfinite(np.r_[yi, yo])):
        raise ValueError("inner and outer states must be finite length-2 [W,V] vectors")
    wi, vi = yi
    wo, vo = yo
    qi, qo = A * vi + nu_prime * wi, A * vo + nu_prime * wo
    pressure_scale = max(abs(alpha * qi), abs(qo), np.finfo(float).tiny)
    result = {
        "q_scaling_residual": abs(qo - alpha * qi) / pressure_scale,
        "pressure_continuity_residual": abs(qo - alpha * qi) / pressure_scale,
    }
    if conversion == "slow":
        scale = max(abs(wi), abs(wo), np.finfo(float).tiny)
        result["slow_radial_displacement_residual"] = abs(wo - wi) / scale
    elif conversion == "rapid":
        xi_i, xi_o = wi - qi / nu_prime, wo - qo / nu_prime
        scale = max(abs(xi_i), abs(xi_o), np.finfo(float).tiny)
        result["rapid_isobar_displacement_residual"] = abs(xi_o - xi_i) / scale
    else:
        raise ValueError("conversion must be 'slow' or 'rapid'")
    return result
