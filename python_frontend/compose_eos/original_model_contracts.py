"""Primary-paper identities and native I/O contracts, with no EOS arithmetic."""
from copy import deepcopy

_BASE = "data/eos_original_papers_20261123"
_BBP_SHA = "f5d15dbc97c8dae24557b55844c2bfc550984aaa50e514e8d9412e11a05ca495"
_NV_SHA = "5611ae66700d983e7f66881839802c2204e16f91cafbf19c670fbd178042bb5b"
_MYERS_SHA = "23d5031429ac8d635fca19db993972d1527b3be03a05f6647042faea01078561"
_BBP = dict(
    region="inner_crust_only", source_path=_BASE + "/research_bbp/BBP1971_ORIGINAL_PARAMETERS.dat",
    doi=["10.1016/0375-9474(71)90281-8"], source_document_sha256=_BBP_SHA,
    default_points=257, requires_composition=True,
    composition="original spherical CLDM neutron gas and nuclear clusters",
    stationary_search_scope="finite multistart stationary candidates in declared coordinate box; no global spherical-branch or pasta minimum proof",
    coexistence_scope="selected gas-free/gas-bearing branch pair near neutron drip; no full convex-hull proof",
    diagnostics=["bbp_original_diagnostics.json", "bbp_original_states.csv"],
)
_MYERS = dict(
    region="outer_crust_only", source_path=_BASE + "/research_myers/myers1966_source.csv",
    doi=["10.1016/S0029-5582(66)80001-9", "10.1086/151216"], source_document_sha256=_MYERS_SHA,
    default_points=257, minimum_points=17, requires_composition=True,
    composition="finite integer mass-candidate search with original deformation formula, BPS electrons and bcc lattice",
    row_policy="phase_complete_pressure_grid", native_grid_receipt="MYERS_BPS_GENERATION.json",
    explicit_coexistence_chain=True,
    diagnostics=["myers_mass_catalog.csv", "myers_bps_states.csv", "myers_bps_transitions.csv", "MYERS_BPS_GENERATION.json"],
)
_ORIGINAL = {
    "BBP1971_original_published": dict(
        **_BBP, representation="Original BBP spherical CLDM and published gas-surface derivative omission",
        source_equivalence="ORIGINAL_FUNCTIONAL_WITH_PUBLISHED_NONVARIATIONAL_CLOSURE",
        derivative_closure="published_omits_surface_neutron_gas_energy_derivative",
        row_policy="exact_requested", grid_policy="log_density_grid_in_declared_original_model_domain",
    ),
    "BBP1971_original_variational": dict(
        **_BBP, representation="Original BBP energy functional with all composition derivatives retained",
        source_equivalence="ORIGINAL_FUNCTIONAL_FULL_VARIATIONAL_POINTWISE_BRANCH",
        derivative_closure="complete_variational",
        row_policy="exact_requested", grid_policy="log_density_grid_raw_pointwise_minima",
    ),
    "BBP1971_variational_maxwell": dict(
        **{**_BBP, "diagnostics": _BBP["diagnostics"] + ["bbp_original_transitions.csv"]},
        representation="Original BBP variational functional with native common-tangent phase coexistence",
        source_equivalence="ORIGINAL_FUNCTIONAL_VARIATIONAL_COMMON_TANGENT_REPRESENTATION",
        derivative_closure="complete_variational_with_maxwell_convexification",
        row_policy="coexistence_density_grid", native_grid_receipt="bbp_original_diagnostics.json",
        explicit_coexistence_chain=True,
        grid_policy="log_density_grid_outside_coexistence_plus_exact_one_sided_endpoints",
    ),
    "NV1973_author_fit": dict(
        region="inner_crust_only",
        representation="Negele-Vautherin1973 Table4 published analytic total-energy fit and exact analytic derivatives",
        source_equivalence="ORIGINAL_1973_AUTHOR_FIT_IDENTICAL_FUNCTION_TO_1974_REPRINT",
        source_path=_BASE + "/research_nv/nv1973_original_source.dat",
        doi=["10.1016/0375-9474(73)90349-7"], source_document_sha256=_NV_SHA,
        semantic_equivalent_model="NV1974_author_fit", full_hartree_fock_rerun=False,
        distinct_new_physical_eos=False,
        composition="thermodynamic-only author fit; original Table3 and bulk functional retained as separate diagnostics",
        default_points=257, requires_composition=False,
        row_policy="exact_requested", grid_policy="log_density_grid_in_original_table_domain",
        diagnostics=["nv1973_states.csv", "nv1973_table3_comparison.csv", "nv1973_uniform_reconstruction.csv", "NV_ORIGINAL_GENERATION.json"],
    ),
    "BPS_Myers1966_magic14": dict(
        **_MYERS, representation="Myers-Swiatecki1966 mass/deformation formula with magic14 and finite BPS Gibbs minimization",
        source_equivalence="PUBLISHED_MASS_FORMULA_FINITE_CANDIDATE_GENERATION",
        magic_choice=14,
    ),
    "BPS_Myers1966_magic20": dict(
        **_MYERS, representation="Myers-Swiatecki1966 mass/deformation formula with magic20 and finite BPS Gibbs minimization",
        source_equivalence="PUBLISHED_MASS_FORMULA_FINITE_CANDIDATE_GENERATION",
        magic_choice=20,
    ),
}
ORIGINAL_IDS = tuple(_ORIGINAL)

def original_crust_catalog():
    return {name: dict(family="primary_crust", **deepcopy(contract))
            for name, contract in _ORIGINAL.items()}
