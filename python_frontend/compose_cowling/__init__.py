"""TOV/Cowling APIs: native Fortran smooth-EOS path plus sharp-Maxwell SciPy reference."""
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


# Explicitly separate sharp-interface Maxwell support from the smooth EOS path.
from .maxwell import MaxwellEOS, MaxwellStar, MaxwellCowlingShot, solve_maxwell_tov, shoot_maxwell_cowling, scan_maxwell_cowling_modes

__all__ += [
    "MaxwellEOS",
    "MaxwellStar",
    "MaxwellCowlingShot",
    "solve_maxwell_tov",
    "shoot_maxwell_cowling",
    "scan_maxwell_cowling_modes",
]
