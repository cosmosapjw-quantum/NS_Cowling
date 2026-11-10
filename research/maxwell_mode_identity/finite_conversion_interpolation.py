"""Algebraic slow/rapid kinematic interpolation; NOT a finite-rate law."""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np

@dataclass(frozen=True)
class AlgebraicInterpolationMap:
    theta: float
    enthalpy_ratio_inner_over_outer: float
    A: float
    nu_prime_per_m: float
    matrix: np.ndarray
    determinant: float
    model_scope: str = "ALGEBRAIC_BOUNDARY_INTERPOLATION_ONLY_NOT_FINITE_RATE_KINETICS"

def _validated_parameters(theta: float, kappa: float, A: float, nu_prime_per_m: float):
    values = tuple(float(x) for x in (theta, kappa, A, nu_prime_per_m))
    if not all(math.isfinite(x) for x in values):
        raise ValueError("theta, kappa, A and nu_prime_per_m must be finite")
    theta, kappa, A, nu_prime_per_m = values
    if not 0.0 <= theta <= 1.0:
        raise ValueError("theta must lie in [0,1]; extrapolation is outside this contract")
    if kappa <= 0.0 or A <= 0.0 or nu_prime_per_m <= 0.0:
        raise ValueError("kappa, A and nu_prime_per_m must be strictly positive")
    return theta, kappa, A, nu_prime_per_m

def algebraic_interpolation_map(theta: float, kappa: float, A: float, nu_prime_per_m: float):
    """Map inner [W,V] to outer [W,V]; theta is NOT a rate or timescale."""
    theta, kappa, A, nu_prime_per_m = _validated_parameters(theta,kappa,A,nu_prime_per_m)
    delta=kappa-1.0
    matrix=np.array([[1.0+theta*delta,theta*delta*A/nu_prime_per_m],
                     [(1.0-theta)*delta*nu_prime_per_m/A,1.0+(1.0-theta)*delta]],dtype=float)
    matrix.setflags(write=False)
    return AlgebraicInterpolationMap(theta,kappa,A,nu_prime_per_m,matrix,float(np.linalg.det(matrix)))

def interpolation_constraint_residual(theta,kappa,A,nu_prime_per_m,state_inner,state_outer):
    theta,kappa,A,nu_prime_per_m=_validated_parameters(theta,kappa,A,nu_prime_per_m)
    yi=np.asarray(state_inner,dtype=float); yo=np.asarray(state_outer,dtype=float)
    if yi.shape!=(2,) or yo.shape!=(2,) or not np.all(np.isfinite(np.r_[yi,yo])):
        raise ValueError("states must be finite length-2 [W,V] vectors")
    wi,vi=yi; wo,vo=yo
    qi=A*vi+nu_prime_per_m*wi; qo=A*vo+nu_prime_per_m*wo
    tiny=np.finfo(float).tiny
    pressure_scale=max(abs(qo),abs(kappa*qi),tiny)
    slow_jump=wo-wi
    rapid_jump=(wo-qo/nu_prime_per_m)-(wi-qi/nu_prime_per_m)
    mixed=(1.0-theta)*slow_jump+theta*rapid_jump
    kin_scale=max(abs((1.0-theta)*slow_jump),abs(theta*rapid_jump),abs(wi),abs(wo),tiny)
    return {"pressure_residual":abs(qo-kappa*qi)/pressure_scale,
            "interpolated_kinematic_residual":abs(mixed)/kin_scale}
