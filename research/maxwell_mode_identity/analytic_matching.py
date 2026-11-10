"""Holomorphic matching-determinant and sharp-interface audit helpers.

The functions in this module are algebraic/numerical utilities.  They do not
assert that a particular production Cowling RHS, endpoint frame, or interface
law has the hypotheses required by an analytic characteristic determinant.
Those hypotheses must be checked against the actual source before use.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np


def _matrix(value: Any, name: str) -> np.ndarray:
    arr = np.asarray(value, dtype=complex)
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1] or arr.shape[0] == 0:
        raise ValueError(f"{name} must be a nonempty square matrix")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def _full_rank(matrix: np.ndarray, name: str) -> None:
    # Numerical rank is deliberately not used as a hidden tolerance gate.
    # Exact singularity of the represented floating matrix is checked through
    # slogdet; conditioning belongs in a separate diagnostic.
    sign, logabs = np.linalg.slogdet(matrix)
    if sign == 0 or not np.isfinite(logabs):
        raise ValueError(f"{name} must be nonsingular")


def matching_determinant(center_frame: Any, surface_frame: Any) -> complex:
    """Return det([Phi_c, Phi_s]) using a bilinear determinant (no conjugation).

    The frame matrices must have n rows, p and n-p columns respectively.
    This computes the raw finite-dimensional matching determinant; holomorphy
    is a property of the parameter-dependent inputs, not this pointwise call.
    """
    pc = np.asarray(center_frame, dtype=complex)
    ps = np.asarray(surface_frame, dtype=complex)
    if pc.ndim != 2 or ps.ndim != 2 or pc.shape[0] != ps.shape[0]:
        raise ValueError("endpoint frames must be rank-2 arrays with the same row count")
    n = pc.shape[0]
    if n == 0 or pc.shape[1] + ps.shape[1] != n:
        raise ValueError("center/surface frame column counts must sum to the state dimension")
    if pc.shape[1] == 0 or ps.shape[1] == 0:
        raise ValueError("both endpoint subspaces must have positive dimension")
    if not np.all(np.isfinite(pc)) or not np.all(np.isfinite(ps)):
        raise ValueError("endpoint frames must contain only finite values")
    return complex(np.linalg.det(np.column_stack((pc, ps))))


def basis_change_determinant(
    center_frame: Any,
    surface_frame: Any,
    center_basis_change: Any,
    surface_basis_change: Any,
) -> tuple[complex, complex]:
    """Apply invertible column-basis changes and return (new D, holomorphic unit).

    Pointwise identity: D_tilde = D det(A_c) det(A_s).  To preserve a
    holomorphic characteristic determinant on a domain, A_c and A_s must be
    holomorphic and invertible throughout that domain, not merely nonzero on
    its contour boundary.
    """
    pc = np.asarray(center_frame, dtype=complex)
    ps = np.asarray(surface_frame, dtype=complex)
    ac = _matrix(center_basis_change, "center_basis_change")
    ass = _matrix(surface_basis_change, "surface_basis_change")
    if pc.ndim != 2 or ps.ndim != 2:
        raise ValueError("endpoint frames must be rank-2 arrays")
    if ac.shape != (pc.shape[1], pc.shape[1]):
        raise ValueError("center basis change has incompatible shape")
    if ass.shape != (ps.shape[1], ps.shape[1]):
        raise ValueError("surface basis change has incompatible shape")
    _full_rank(ac, "center_basis_change")
    _full_rank(ass, "surface_basis_change")
    unit = complex(np.linalg.det(ac) * np.linalg.det(ass))
    new_d = matching_determinant(pc @ ac, ps @ ass)
    return new_d, unit


def gauge_transform_system(matrix_y: Any, gauge_g: Any, gauge_g_prime: Any) -> np.ndarray:
    """Transform y=G(r,lambda) z: M_z=G^-1 M_y G-G^-1 G'.

    G must be invertible.  For a holomorphic characteristic determinant, G
    must additionally be holomorphic in lambda and piecewise C1 in r.
    """
    m = _matrix(matrix_y, "matrix_y")
    g = _matrix(gauge_g, "gauge_g")
    gp = _matrix(gauge_g_prime, "gauge_g_prime")
    if g.shape != m.shape or gp.shape != m.shape:
        raise ValueError("matrix_y, gauge_g and gauge_g_prime must have equal square shapes")
    _full_rank(g, "gauge_g")
    # Solve rather than form G^{-1} explicitly.
    return np.linalg.solve(g, m @ g - gp)


def gauge_transform_interface(transfer_y: Any, gauge_inner: Any, gauge_outer: Any) -> np.ndarray:
    """Transform y_o=T_y y_i into z_o=T_z z_i, T_z=G_o^-1 T_y G_i."""
    t = _matrix(transfer_y, "transfer_y")
    gi = _matrix(gauge_inner, "gauge_inner")
    go = _matrix(gauge_outer, "gauge_outer")
    if gi.shape != t.shape or go.shape != t.shape:
        raise ValueError("transfer and one-sided gauges must have equal square shapes")
    _full_rank(gi, "gauge_inner")
    _full_rank(go, "gauge_outer")
    return np.linalg.solve(go, t @ gi)


def gauge_determinant_factor(gauge_at_match: Any) -> complex:
    """Return det(G) in D_y=det(G) D_z at a common matching point."""
    g = _matrix(gauge_at_match, "gauge_at_match")
    _full_rank(g, "gauge_at_match")
    return complex(np.linalg.det(g))


def smooth_match_radius_factor(trace_m: complex, delta_r: float) -> complex:
    """Return exp(integral tr(M)dr) for constant trace over a smooth interval."""
    tr = complex(trace_m)
    dr = float(delta_r)
    if not (math.isfinite(tr.real) and math.isfinite(tr.imag) and math.isfinite(dr)):
        raise ValueError("trace_m and delta_r must be finite")
    return complex(np.exp(tr * dr))


@dataclass(frozen=True)
class InterfaceMatchingForm:
    """Research representation of a 2-state inner/outer interface relation."""
    kind: str
    theta: float
    kappa: float
    A: complex
    nu_prime: float
    C_inner: np.ndarray | None
    C_outer: np.ndarray | None
    T: np.ndarray | None
    determinant: complex | None
    implicit_to_transfer_factor: complex | None
    scope: str = "ALGEBRAIC_BOUNDARY_CONSTRAINT_ONLY_NOT_FINITE_RATE_KINETICS"


def interface_matching_form(theta: float, kappa: float, A: complex, nu_prime: float) -> InterfaceMatchingForm:
    """Build the LA10 constraint blocks or a regular transfer map.

    Convention: inner is the high-density side, outer the low-density side,
    [X]_i^o=X_o-X_i, Q=A V+nu_prime W, Q_o=kappa Q_i.

    For 0<=theta<1, the implicit kinematic constraint is
      (1-theta)*nu_prime*[W] - theta*A*[V] = 0.
    For A!=0 and kappa!=1 it is a graph z_o=T z_i.  At A=0 its block form
    remains finite while the graph transfer has a pole in T_21.  For theta=1
    we use the regular analytic endpoint extension [V]=0, equivalent to the
    rapid constraint for A!=0.  For kappa=1 there is no physical jump and the
    exact identity map is returned before any rank-degenerate block form.
    """
    th = float(theta)
    kap = float(kappa)
    nu_p = float(nu_prime)
    aval = complex(A)
    if not all(map(math.isfinite, (th, kap, nu_p, aval.real, aval.imag))):
        raise ValueError("theta, kappa, A and nu_prime must be finite")
    if not 0.0 <= th <= 1.0:
        raise ValueError("theta must lie in [0,1]")
    if kap <= 0.0 or nu_p == 0.0:
        raise ValueError("kappa must be positive and nu_prime nonzero")

    delta = kap - 1.0
    eye = np.eye(2, dtype=complex)
    if kap == 1.0:
        return InterfaceMatchingForm("transfer", th, kap, aval, nu_p, None, None, eye, 1.0+0.0j, 1.0+0.0j)

    # The rapid endpoint admits the regularized [V]=0 form, including A=0.
    if th == 1.0:
        ci = np.array([[-kap*nu_p, -kap*aval], [0.0, -1.0]], dtype=complex)
        co = np.array([[nu_p, aval], [0.0, 1.0]], dtype=complex)
        t = np.array([[kap, delta*aval/nu_p], [0.0, 1.0]], dtype=complex)
        return InterfaceMatchingForm("transfer", th, kap, aval, nu_p, ci, co, t,
                                     complex(np.linalg.det(t)), complex(-nu_p))

    # Intermediate/slow algebraic family: retain implicit constraints at A=0.
    ci = np.array([
        [-kap*nu_p, -kap*aval],
        [-(1.0-th)*nu_p, th*aval],
    ], dtype=complex)
    co = np.array([
        [nu_p, aval],
        [(1.0-th)*nu_p, -th*aval],
    ], dtype=complex)
    if aval == 0.0:
        return InterfaceMatchingForm("implicit", th, kap, aval, nu_p, ci, co, None, None, None)

    # With A!=0, det(Co)=-A*nu_prime, hence the transfer graph exists.
    t = np.linalg.solve(co, -ci)
    return InterfaceMatchingForm("transfer", th, kap, aval, nu_p, ci, co, t,
                                 complex(np.linalg.det(t)), complex(aval*nu_p))


def implicit_interface_determinant(form: InterfaceMatchingForm, state_inner: Any, state_outer: Any) -> complex:
    """Return det[C_i y_i, C_o y_o] for a one-dimensional 2D matching problem.

    A zero is equivalent to the existence of nonzero endpoint amplitudes
    satisfying the two interface constraints, provided the endpoint vectors
    are nonzero and the resulting 2x2 amplitude matrix is interpreted with
    its actual rank. This is not a general n-state implementation.
    """
    if form.C_inner is None or form.C_outer is None:
        raise ValueError("this form has no implicit constraint blocks")
    yi = np.asarray(state_inner, dtype=complex)
    yo = np.asarray(state_outer, dtype=complex)
    if yi.shape != (2,) or yo.shape != (2,):
        raise ValueError("state_inner and state_outer must be length-2 vectors")
    if not np.all(np.isfinite(yi)) or not np.all(np.isfinite(yo)):
        raise ValueError("interface states must be finite")
    cols = np.column_stack((form.C_inner @ yi, form.C_outer @ yo))
    return complex(np.linalg.det(cols))


def transfer_matching_determinant(transfer: Any, state_inner: Any, state_outer: Any) -> complex:
    """Return det(T y_inner, y_outer), the graph-map matching scalar."""
    t = _matrix(transfer, "transfer")
    yi = np.asarray(state_inner, dtype=complex)
    yo = np.asarray(state_outer, dtype=complex)
    if t.shape != (2, 2) or yi.shape != (2,) or yo.shape != (2,):
        raise ValueError("transfer must be 2x2 and states must be length-2")
    if not np.all(np.isfinite(yi)) or not np.all(np.isfinite(yo)):
        raise ValueError("states must be finite")
    return complex(np.linalg.det(np.column_stack((t @ yi, yo))))


def dimensionless_interpolation_transfer(form: InterfaceMatchingForm, length_scale_m: float) -> np.ndarray:
    """Return LA10 T in the dimensionless state (W/L^3,V/L^2)."""
    if form.kind != "transfer" or form.T is None:
        raise ValueError("a graph transfer is unavailable for this interface form")
    L = float(length_scale_m)
    if not math.isfinite(L) or L <= 0.0:
        raise ValueError("length_scale_m must be finite and positive")
    scale = np.diag([L**3, L**2]).astype(complex)
    return np.linalg.solve(scale, form.T @ scale)


def transfer_condition_number_diagnostic(form: InterfaceMatchingForm, length_scale_m: float) -> float:
    """Return a 2-norm condition number for real-frequency diagnostics only.

    This uses singular values and conjugation, so it is deliberately NOT a
    holomorphic characteristic function and must never replace the raw
    matching determinant on a complex-frequency contour.
    """
    t_hat = dimensionless_interpolation_transfer(form, length_scale_m)
    if abs(form.A.imag) > 0.0 or form.A.real <= 0.0:
        raise ValueError("condition number diagnostic is restricted to positive real A")
    value = float(np.linalg.cond(t_hat, 2))
    if not math.isfinite(value):
        raise ArithmeticError("transfer condition number is not finite")
    return value
