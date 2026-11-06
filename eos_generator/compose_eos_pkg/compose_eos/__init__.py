"""compose_eos v3.3.0 — CompOSE EoS Parser, Thermodynamic Auditor, Gap Reconstructor,
First-Principles Literature Solver, Causal TOV Runner, and Universal Crust–Core Matcher.
"""
from .eos_data import EOSData
from .io import read_compose, write_compose
from .consistency import check_consistency
from .reconstruct import fill_density_gap, find_density_gap, reconstruct_master
from .plotting import (
    make_figures,
    plot_pressure,
    plot_sound_speed,
    plot_composition,
    plot_gibbs_duhem,
    plot_mass_radius,
)
from .tov import (
    BarotropicEOS,
    TOVConfig,
    CausalDomain,
    TOVStarResult,
    find_causal_domain,
    first_law_nb_consistency,
    solve_tov_star,
    solve_tov_sequence,
    stable_branch_summary,
)
from .cross_verify import cross_verify_eos
from .matching import (
    METHODS,
    CRUST_NAMES,
    CORE_NAMES,
    get_core_start_index,
    get_gauge_shifts,
    find_convex_transition_window,
    build_monotone_hermite_bridge,
    build_maxwell_plateau_bridge,
    match_crust_core,
    solve_star_3stage,
    solve_star_at_target_mass,
)

__version__ = "3.3.0"

__all__ = ["build_maxwell_css_eos",
    "MaxwellEOS", "GibbsResult", "ChargeResolvedPhase", "LeptonSector",
    "build_maxwell_eos", "write_maxwell_eos", "read_maxwell_eos",
    "tidal_y_jump", "solve_tov_with_phase_transition", "gibbs_mixed_phase",
    "write_gibbs_eos", "read_gibbs_eos", "causality_profile", "causal_completion",
    "EOSData",
    "read_compose",
    "write_compose",
    "check_consistency",
    "fill_density_gap",
    "find_density_gap",
    "reconstruct_master",
    "make_figures",
    "plot_pressure",
    "plot_sound_speed",
    "plot_composition",
    "plot_gibbs_duhem",
    "plot_mass_radius",
    "BarotropicEOS",
    "TOVConfig",
    "CausalDomain",
    "TOVStarResult",
    "find_causal_domain",
    "first_law_nb_consistency",
    "solve_tov_star",
    "solve_tov_sequence",
    "stable_branch_summary",
    "cross_verify_eos",
    "METHODS",
    "CRUST_NAMES",
    "CORE_NAMES",
    "get_core_start_index",
    "get_gauge_shifts",
    "find_convex_transition_window",
    "build_monotone_hermite_bridge",
    "build_maxwell_plateau_bridge",
    "match_crust_core",
    "solve_star_3stage",
    "solve_star_at_target_mass",
]

from .phase_transitions import (MaxwellEOS, GibbsResult, ChargeResolvedPhase, LeptonSector, build_maxwell_eos, write_maxwell_eos, read_maxwell_eos, tidal_y_jump, solve_tov_with_phase_transition, gibbs_mixed_phase, write_gibbs_eos, read_gibbs_eos)
from .causality import causality_profile, causal_completion

from .phase_transitions import build_maxwell_css_eos
