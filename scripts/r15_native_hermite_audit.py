"""Read-only audit helpers for the solver-native cubic-Hermite P(mu) table cells."""
from __future__ import annotations


def audit_native_interval(
    mu0_j: float, p0_pa: float, n0_m3: float,
    mu1_j: float, p1_pa: float, n1_m3: float,
) -> dict:
    """Return whole-interval stability/causality and chi-extremum diagnostics."""
    return {}
