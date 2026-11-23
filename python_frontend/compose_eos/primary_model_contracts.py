"""Literature identities and I/O contracts; no physical calculations."""
from copy import deepcopy

_BASE = "data/eos_primary_crust_generation_20261123"
_PRIMARY = {
    "BPS1971_sequence": {
        "region": "outer_crust_only",
        "representation": "BPS electron gas and bcc lattice; fixed-pressure Gibbs minimization over the 12 published equilibrium nuclei",
        "source_equivalence": "PRIMARY_EQUATIONS_RESTRICTED_PUBLISHED_SEQUENCE",
        "composition": "original Table 2 nuclei; modern fundamental constants; explicit one-sided phase transitions",
        "source_path": _BASE + "/research_bps/bps1971_table2.csv",
        "doi": ["10.1086/151216"],
        "minimum_points": 16,
        "default_points": 257,
        "grid_policy": "target_pressure_grid_plus_all_phase_endpoints_and_interiors",
        "diagnostics": ["primary_source_diagnostics.json", "bps_states.csv", "bps_transitions.csv"],
    },
    "BBP1971_table5": {
        "region": "inner_crust_only",
        "representation": "Exact published BBP thermodynamic knots in BPS Table 5, through the declared n<=0.08 fm^-3 domain",
        "source_equivalence": "PRIMARY_PUBLISHED_TABLE;THERMODYNAMIC_INCONSISTENCY_PRESERVED",
        "composition": "thermodynamic-only export; no inferred free/cluster inventory",
        "source_path": _BASE + "/research_bbp/BPS1971_TABLE5_BBP_CGS.dat",
        "doi": ["10.1086/151216", "10.1016/0375-9474(71)90281-8"],
        "default_points": 26,
        "grid_policy": "exact_26_primary_knots_no_resampling",
        "diagnostics": ["primary_source_diagnostics.json", "primary_source_values.csv", "primary_source_anchor_comparison.csv"],
    },
    "BBP1971_pressure_integral": {
        "region": "inner_crust_only",
        "representation": "Derived piecewise-polytropic primary pressure table with one original energy anchor and analytic cold-first-law integration",
        "source_equivalence": "DERIVED_PRESSURE_PRIMARY_REPRESENTATION_NOT_FULL_TABLE_OR_CLDM_REPRODUCTION",
        "composition": "thermodynamic-only export; no inferred free/cluster inventory",
        "source_path": _BASE + "/research_bbp/BPS1971_TABLE5_BBP_CGS.dat",
        "doi": ["10.1086/151216", "10.1016/0375-9474(71)90281-8"],
        "default_points": 257,
        "grid_policy": "log_density_grid_plus_all_26_primary_knots",
        "diagnostics": ["primary_source_diagnostics.json", "primary_source_values.csv", "primary_source_anchor_comparison.csv"],
    },
    "NV1974_author_fit": {
        "region": "inner_crust_only",
        "representation": "Negele 1974 Table III analytic total-energy representation, differentiated analytically in Fortran",
        "source_equivalence": "AUTHOR_ANALYTIC_REPRESENTATION_NOT_1973_HARTREE_FOCK_RERUN",
        "composition": "thermodynamic-only analytic representation; original Table II retained as validation data",
        "source_path": _BASE + "/research_nv/nv1974_primary_source.dat",
        "doi": ["10.1017/S0074180900099885", "10.1016/0375-9474(73)90349-7"],
        "default_points": 257,
        "grid_policy": "log_density_grid_in_original_table_domain",
        "diagnostics": ["NV_PRIMARY_GENERATION.json", "nv_primary_states.csv", "nv_primary_source_comparison.csv"],
    },
}
PRIMARY_IDS = tuple(_PRIMARY)

def primary_crust_catalog():
    return {name: dict(family="primary_crust", **deepcopy(contract))
            for name, contract in _PRIMARY.items()}
