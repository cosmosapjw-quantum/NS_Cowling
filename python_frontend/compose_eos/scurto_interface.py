"""Spherical interface functions from Scurto et al. 2024, Eqs. 24–25.

These functions alone are not a calibrated CLD model.  In particular, b is
not the exponent of the historical proxy, and model-specific mass-fit inputs
and the rest of the functional are required before making a source claim.
The numerical support bounds below are conservative float64 implementation
limits, not physical constraints on a calibrated parameter set.
See https://arxiv.org/html/2402.15548v1.
"""
import math


def spherical_phi(f: float) -> float:
    """Spherical Wigner-Seitz geometry factor for 0 <= f <= 1."""
    if not math.isfinite(f) or not 0.0 <= f <= 1.0:
        raise ValueError("volume fraction must be finite and in [0, 1]")
    x = f ** (1.0 / 3.0)
    one_minus_x = (1.0-f)/(1.0+x+x*x)
    return one_minus_x**2 * (x + 2.0) / 5.0


def _terms(y: float, sigma_0: float, b: float):
    if not (math.isfinite(y) and 0.0 < y < 1.0 and
            math.isfinite(sigma_0) and sigma_0 >= 0.0 and
            math.isfinite(b) and b > -16.0):
        raise ValueError("require 0 < y < 1, sigma_0 >= 0, finite b > -16")
    if not (1e-8 <= y <= 1.0-1e-8 and sigma_0 <= 1e6 and b <= 1e100):
        raise ValueError("parameters outside supported float64 range")
    q = y*(1.0-y)
    k = b+16.0
    scale = max(k, 1.0)
    # d = y^3+(1-y)^3+b*q^3, written without cancellation
    # near b=-16 and y=1/2, then scaled for large b.
    q3 = q**3
    d = ((1.0-2.0*y)**2 * (1.0+q+4.0*q*q))/scale + (k/scale)*q3
    if d == 0.0:
        raise ValueError("surface parameters outside floating-point range")
    return q, q3, k/scale, scale, d


def surface_sigma(y: float, sigma_0: float, b: float) -> float:
    """Surface tension in MeV/fm² when sigma_0 is given in MeV/fm²."""
    _, q3, k_scaled, _, d = _terms(y, sigma_0, b)
    return sigma_0 * (k_scaled*q3/d)


def surface_sigma_prime(y: float, sigma_0: float, b: float) -> float:
    """Derivative of surface tension with respect to dense-phase proton fraction."""
    q, _, k_scaled, scale, d = _terms(y, sigma_0, b)
    numerator = 3.0*q*q*(1.0-2.0*y)*(1.0-2.0*q)
    return sigma_0 * (k_scaled * (numerator/scale) / d / d)
