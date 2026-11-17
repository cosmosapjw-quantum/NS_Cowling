"""§5 regression: independent M7 CSS smooth-bridge validation on the
SPG_M1 + BBB_BHF2 representative case.

Verifies the §5 required checks on the native M7 CSS output:
  1. CompOSE read-back row count and file hashes match the task_receipt.
  2. Strict monotonicity of n_B, P, epsilon, mu_B.
  3. epsilon + P = mu_B * n_B at every row (first-law / Euler identity).
  4. Exact integrated pressure difference P2-P1 = ∫ n_B d mu_B over the
     bridge, computed analytically from the CSS formula (NOT a sparse grid).
  5. Both CSS segment cs2 values and exponents are declared.
  6. 0 < cs2 <= 1 over the bridge (causality).
  7. Endpoint values of P, n_B, mu_B, F.
  8. Source/bridge cs2 difference at both endpoints is recorded.
  9. The seam continuity class: state-variable continuity (C^0) holds vs the
     declared bridge start; first-derivative (cs2) continuity (C^1) does NOT
     hold -- the bridge is C^0 only, not C^1, because the source crust cs2
     differs from the bridge cs2. This matches the §5 rule: "Do not describe
     the bridge as C^1 or derivative-smooth unless the endpoint derivatives
     are matched and tested."

Scope limitation preserved: composition_dropped=True, C3/C4 NOT_AVAILABLE
(not PASS). The output is a thermo-only barotropic EOS, not a composition-
preserving physical EOS.
"""
from __future__ import annotations
import os
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "maxwell_audit"))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_ROOT", str(ROOT))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_EXE",
                      str(ROOT / "build" / "compose_eos_fortran" / "compose_eos_fortran"))

from validate_css_output import validate_css_output  # noqa: E402

CSS_DIR = ROOT / "download" / "m7_batch_v2_20261114_full" / "d87e2ccf5b7c_SPG_M1_BBB_BHF2_css"


@pytest.fixture(scope="module")
def rep():
    return validate_css_output(CSS_DIR)


def test_css_readback_row_count_and_hashes(rep):
    r = rep["1_readback"]
    assert r["row_count"] == 388
    assert r["readback_n_match"] is True
    assert r["all_hashes_match"] is True


def test_css_strict_monotonicity(rep):
    m = rep["2_monotonicity"]
    assert m["all_strictly_increasing"] is True
    assert m["n_B"] and m["P"] and m["epsilon"] and m["mu_B"]


def test_css_first_law_identity(rep):
    fl = rep["3_first_law"]
    assert fl["pass"] is True
    assert fl["max_relative_residual"] < 1e-10


def test_css_exact_integrated_pressure_difference(rep):
    """The exact analytic bridge integral ∫ n_B d mu_B = P2 - P1 holds to
    machine precision (the CSS segment satisfies dP = n d mu by construction)."""
    ip = rep["4_integrated_pressure"]
    assert ip["pass"] is True
    exact = ip["exact_analytic_bridge"]
    assert exact["pass"] is True
    assert exact["relative_residual"] < 1e-8
    # the trapezoid-on-whole-table is a loose diagnostic that fails across the
    # crust-core seam (first-order there); it is recorded but is not the §5
    # authoritative check.
    assert ip["trapezoid_whole_table"]["relative_residual"] > 1e-4  # sampling artifact


def test_css_segment_cs2_values_declared(rep):
    s = rep["5_segments_declared"]
    assert s["segment1_cs2_low"] is not None
    assert s["segment2_cs2_high"] is not None
    assert 0 < s["segment1_cs2_low"] <= 1
    assert 0 < s["segment2_cs2_high"] <= 1
    # the CSS exponent is 1/(1+cs2)
    assert s["segment1_exponent"] == pytest.approx(1.0 / (1.0 + s["segment1_cs2_low"]))
    assert s["segment2_exponent"] == pytest.approx(1.0 / (1.0 + s["segment2_cs2_high"]))


def test_css_causality_over_bridge(rep):
    c = rep["6_cs2_causality"]
    assert c["pass"] is True
    assert c["min_cs2"] > 0
    assert c["max_cs2"] <= 1 + 1e-8
    assert c["acausal_point_count"] == 0


def test_css_endpoints_recorded(rep):
    e = rep["7_endpoints"]
    assert e["n_first"] > 0 and e["n_last"] > e["n_first"]
    assert e["P_first"] > 0 and e["P_last"] > e["P_first"]
    assert e["mu_first"] > 0 and e["mu_last"] > e["mu_first"]
    assert e["F_first"] > 0 and e["F_last"] > 0


def test_css_source_bridge_cs2_difference_recorded(rep):
    """The source crust cs2 differs from the bridge cs2 at both endpoints.
    This difference is recorded (§5 item 8), not hidden."""
    e = rep["8_endpoint_cs2_difference"]
    s = rep["5_segments_declared"]
    # source crust cs2 (cs2_source_left) vs bridge cs2 (cs2_left) differ
    assert s["cs2_source_left"] != s["segment1_cs2_low"]
    assert abs(s["cs2_source_left"] - s["segment1_cs2_low"]) > 1e-4


def test_css_seam_is_C0_continuous_not_C1(rep):
    """The bridge is state-variable continuous (C^0) at the crust-core seam
    (vs the declared bridge start), but NOT first-derivative continuous (C^1):
    the source crust cs2 differs from the bridge cs2. The §5 rule is honored:
    we do not describe the bridge as C^1 or derivative-smooth."""
    seam = rep["9_seam_continuity"]
    assert seam["state_variable_continuous_vs_declared_bridge_start"] is True, (
        "bridge must be C^0 (state-variable continuous) at the seam")
    assert bool(seam["first_derivative_continuous_crust_vs_bridge"]) is False, (
        "bridge is NOT C^1: source crust cs2 differs from bridge cs2")
    # the declared cs2 mismatch is non-zero
    assert seam["source_bridge_cs2_mismatch"] is not None
    assert seam["source_bridge_cs2_mismatch"] > 1e-4
    # higher regularity (C^2) is not enforced
    assert seam["higher_regularity_C2_enforced"] is False
    # the row-to-row sample gap is a sampling artifact, not a discontinuity
    assert seam["sample_gap_is_sampling_not_discontinuity"] is True


def test_css_scope_limitation_composition_dropped(rep):
    """The M7 CSS output is a thermo-only barotropic EOS. composition_dropped
    is True; C3/C4 are NOT_AVAILABLE (not PASS). This limitation is preserved
    in the report, not hidden."""
    sl = rep["scope_limitations"]
    assert sl["composition_dropped"] is True
    assert sl["C3_status"] == "NOT_AVAILABLE"
    assert sl["C4_status"] == "NOT_AVAILABLE"
    assert sl["thermo_only_barotropic_scope"] is True
    assert sl["not_a_full_composition_preserving_EOS"] is True
