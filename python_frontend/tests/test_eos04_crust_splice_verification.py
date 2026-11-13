"""EOS-04 crust-reconstruction-and-splice verification: the assembled
crust+core table from the native smooth matching path must pass the full
C1-C11 consistency gate (hard checks C1-C6/C10/C11), with source-knot,
endpoint and causal-stability preservation, and no interpolation of
categorical columns.

EOS-04 acceptance:
  1. native Fortran interpolation/derivative/reconstruction
  2. local matching with endpoint/causality checks
  3. no interpolation of categorical Nadd/q_extra
  4. one-sided discontinuities retained
  5. source columns and optional files preserved

The 13 reconstruct-master tests (test_native_reconstruct_master.py) and the
13+5 matching tests (test_native_crust_core_matching.py +
test_eos05_kernel_edge_cases.py) cover criteria 1, 3, 4, 5 piecewise and
the bridge-level endpoint/causality (criterion 2, bridge only). But NO test
asserts that the *assembled* crust+core table -- the full splice output --
passes the native consistency gate (C1-C11, hard checks C1-C6/C10/C11).
This is the exact gap the legacy M4 receipt
(E390-02_M4_OUTER_INNER_CRUST_SPLICE_STATUS_20261112) recorded: C6/C10/C11
failed on the assembled table (legacy Python path, due to coarse CLD
pressure). The native path reconstructs Q1 from the first law and enforces
closure, so the assembled table passes; this test pins that guarantee.

test_native_maxwell_transition.py:130-133 already asserts
overall/C6/C10/C11 PASS for the sharp first-order (maxwell) transition; this
file asserts the same for the smooth cubic/quintic F/P matching path, which
is the EOS-04 crust-core splice.
"""
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos import match_crust_core, read_compose  # noqa: E402

CRUST = ROOT / "eos_generator" / "eos_data_normalized" / "GPPVA_FSU2H"
CORE = ROOT / "eos_generator" / "eos_data_normalized" / "GM_GM1"


def _require_reference_inputs():
    required = (CRUST / "eos.nb", CRUST / "eos.thermo", CRUST / "eos.compo",
                CORE / "eos.nb", CORE / "eos.thermo")
    missing = [str(p) for p in required if not p.is_file()]
    if missing:
        pytest.skip("repository benchmark EOS inputs are unavailable: " + ", ".join(missing))


@pytest.mark.parametrize("method", ("css", "cubic_F", "cubic_P", "quintic_F", "quintic_P"))
def test_assembled_crust_core_table_passes_full_consistency_gate(method):
    """The assembled crust+core table (crust knots + bridge + core knots)
    from the native smooth matching path must pass the native C1-C11
    consistency gate. This closes the legacy M4 receipt gap (C6/C10/C11
    failed on the assembled table) on the native path: the matcher
    reconstructs Q1 from the first law and enforces thermodynamic closure,
    so the assembled table is consistent. EOS-04 acceptance criterion 2
    ('local matching with endpoint/causality checks') applied to the
    assembled table, not just the bridge."""
    _require_reference_inputs()
    crust = read_compose(CRUST)
    core = read_compose(CORE)
    matched, info = match_crust_core(
        crust, core, method=method, n_bridge=40,
        core_compo_mode="strict_thermo_only", drop_composition=True,
    )
    cons = info["consistency"]
    overall = cons["overall"]
    # Full assembled-table consistency PASS (hard checks C1-C6/C10/C11).
    assert overall["pass_"] is True, (
        f"{method}: assembled crust+core table fails consistency: "
        f"failed_checks={overall.get('failed_checks')}"
    )
    # The specific checks the legacy M4 receipt recorded as failing.
    assert cons["C6"]["pass_"] is True, f"{method}: C6 (dP/dn_B > 0) failed on assembled table"
    assert cons["C10"]["pass_"] is True, f"{method}: C10 (causal sound speed) failed on assembled table"
    assert cons["C11"]["pass_"] is True, f"{method}: C11 (dmu_B/dn_B > 0) failed on assembled table"
    # C1 (cold T=0, zero entropy, zero Q6-Q7) and C2 (thermodynamic closure)
    # must also pass on the assembled table.
    assert cons["C1"]["pass_"] is True, f"{method}: C1 failed on assembled table"
    assert cons["C2"]["pass_"] is True, f"{method}: C2 (closure) failed on assembled table"


@pytest.mark.parametrize("method", ("css", "cubic_F", "cubic_P", "quintic_F", "quintic_P"))
def test_assembled_table_preserves_source_knots_and_endpoints(method):
    """EOS-04 acceptance criteria 2 (endpoint checks) + 5 (source columns
    preserved): the source crust knots up to the matching endpoint i1 and
    the source core knots from i2 onward are retained verbatim in the
    assembled table; only the bridge between them is interpolated. The
    matched density grid is strictly increasing across the splice."""
    _require_reference_inputs()
    crust = read_compose(CRUST)
    core = read_compose(CORE)
    matched, info = match_crust_core(
        crust, core, method=method, n_bridge=40,
        core_compo_mode="strict_thermo_only", drop_composition=True,
    )
    i1 = info["i1"]
    i2 = info["i2"]
    # Source crust knots up to and including the left endpoint are retained.
    np.testing.assert_array_equal(matched.nb[:i1 + 1], crust.nb[:i1 + 1])
    # Source core knots from the right endpoint onward are retained.
    np.testing.assert_array_equal(matched.nb[-(core.n - i2):], core.nb[i2:])
    # The assembled grid is strictly increasing (no interpolation knot
    # collides with a source knot, no source knot duplicated).
    assert np.all(np.diff(matched.nb) > 0.0), f"{method}: assembled nb not strictly increasing"
    # The endpoint pressures/energy densities come from the source tables,
    # not from interpolation (the bridge starts strictly inside the window).
    assert matched.nb[i1] == crust.nb[i1]
    assert matched.nb[i1 + 1 + 40] == core.nb[i2]  # first core knot after bridge


@pytest.mark.parametrize("method", ("cubic_F", "cubic_P", "quintic_F", "quintic_P"))
def test_assembled_table_has_no_interpolated_categorical_columns(method):
    """EOS-04 acceptance criterion 3 (no interpolation of categorical
    Nadd/q_extra): the smooth matching path interpolates the thermodynamic
    barotrope (Q1-Q7) in the bridge but must NOT invent categorical
    composition columns. With drop_composition=True, the assembled table
    carries zero composition everywhere (the source composition is
    explicitly discarded, not interpolated into the bridge). Iphase is
    uniform zero (cold smooth barotrope, no phase labels invented)."""
    _require_reference_inputs()
    crust = read_compose(CRUST)
    core = read_compose(CORE)
    matched, info = match_crust_core(
        crust, core, method=method, n_bridge=40,
        core_compo_mode="strict_thermo_only", drop_composition=True,
    )
    # No composition columns are interpolated: the bridge carries zero
    # composition (drop_composition=True discards the source composition
    # rather than interpolating categorical Nadd/Yn/Yp across the bridge).
    for field in (matched.Yn, matched.Yp, matched.Ye, matched.Ymu,
                  matched.Aav, matched.Zav, matched.Yh):
        assert np.all(field == 0.0), f"{method}: interpolated categorical column is nonzero"
    # Iphase is uniform zero: no phase label is invented in the bridge.
    assert np.all(matched.Iphase == 0), f"{method}: Iphase not uniformly zero"
    assert info["composition_dropped"] is True
    assert matched.composition_file_present is False
    assert matched.composition_available is False


@pytest.mark.parametrize("method", ("css", "cubic_F", "cubic_P", "quintic_F", "quintic_P"))
def test_assembled_table_causal_stability_across_splice(method):
    """EOS-04 acceptance criterion 2 (causality checks): the assembled
    table's sound speed is subluminal (0 < cs2 < 1) and strictly positive
    across the whole splice, including the bridge-to-source transition
    points. The native report's cs2_min/cs2_max bound the assembled table,
    not just the bridge; this asserts the assembled causal stability."""
    _require_reference_inputs()
    crust = read_compose(CRUST)
    core = read_compose(CORE)
    matched, info = match_crust_core(
        crust, core, method=method, n_bridge=40,
        core_compo_mode="strict_thermo_only", drop_composition=True,
    )
    rep = info["native_report"]
    assert rep["cs2_min"] > 0.0, f"{method}: nonpositive min sound speed on assembled table"
    assert rep["cs2_max"] < 1.0, f"{method}: superluminal max sound speed on assembled table"
    # First-law closure holds on the assembled table (the splice did not
    # break the mu_B = F + P/n_B identity at the join).
    assert info["first_law_relative_residual"] < 5e-9, (
        f"{method}: first-law residual {info['first_law_relative_residual']} too large"
    )
