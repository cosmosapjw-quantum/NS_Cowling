"""Fortran-backed Python frontend for TOV and relativistic Cowling calculations."""
from .adapter import build_adapter, build_adapter_csv
from .direct_barotrope import DirectBarotrope, audit_csv
from .solver import solve_tov_star, solve_tov_sequence, solve_cowling_spectrum

__all__ = [
    "build_adapter",
    "build_adapter_csv",
    "DirectBarotrope",
    "audit_csv",
    "solve_tov_star",
    "solve_tov_sequence",
    "solve_cowling_spectrum",
]

