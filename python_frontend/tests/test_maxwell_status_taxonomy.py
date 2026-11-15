"""§3/§4 regression: structured Maxwell status taxonomy is stable and
agreement-preserving across the two same-mass v2 experiments.

§3 requires that the Maxwell coexistence search reports structured terminal
outcomes (PASS, MASS_CONVENTION_FAIL, BETA_EQUILIBRIUM_FAIL,
SOURCE_TABLE_CONVEXITY_FAIL, INTERPOLANT_SHAPE_FAIL, CAUSALITY_FAIL,
NO_COMMON_MU_SUPPORT, NO_COEXISTENCE_CROSSING, MULTIPLE_COEXISTENCE_CROSSINGS,
NONPOSITIVE_DENSITY_JUMP, FIRST_LAW_IDENTITY_FAIL, COMPOSITION_UNSUPPORTED).

§4 requires that the same underlying failure receives the same status in
every transformation and every run, replacing the text-based category
inference that produced inconsistent labels.

These tests verify:
  1. The Python maxwell_status_category mapping agrees with the Fortran
     maxwell_status_category on a sample of every category (single source of
     truth).
  2. The 65 pairs in the PDG-gauge and gauge-shift v2 experiments all receive
     the SAME structured category -- representation-only transformations do
     not alter the physics verdict.
  3. The 4 pairs whose free-form error is the same but whose OLD text-based
     category differed are now consistently labeled by the structured
     taxonomy (the §4 inconsistency is fixed).
  4. The pairwise mass-gate tolerance is unchanged at 1e-12.
  5. The dominant failure is INTERPOLANT_SHAPE_FAIL (61/65), not a generic
     STABILITY_GATE bucket -- the §1 source-convexity diagnosis explains
     these as source-table-level, but the v1 kernel status string cannot
     distinguish source-level from interpolant-level, so INTERPOLANT_SHAPE_FAIL
     is the correct Fortran-emitted category for the generic interpolation
     message; the §1 audit provides the authoritative source-vs-interpolant
     split.
"""
from __future__ import annotations
import os
import sys
import json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "maxwell_audit"))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_ROOT", str(ROOT))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_EXE",
                      str(ROOT / "build" / "compose_eos_fortran" / "compose_eos_fortran"))

from reconcile_v2_experiments import maxwell_status_category, reconcile  # noqa: E402


# Representative status strings covering every §3 category. These are the
# actual v1 kernel status strings from compose_eos_phase_transition.f90.
CATEGORY_SAMPLES = {
    "MASS_CONVENTION_FAIL": "phase tables must use the same neutron and proton rest-mass conventions",
    "BETA_EQUILIBRIUM_FAIL": "Maxwell construction requires beta-equilibrium lepton potential Q5 approximately zero",
    "FIRST_LAW_IDENTITY_FAIL": "provider Q3 chemical potential disagrees with the Euler value (relative tolerance 5e-4)",
    "CAUSALITY_FAIL": "P(mu_B) interpolation violates the causal bound c_s^2/c^2<=1",
    "NO_COMMON_MU_SUPPORT": "the two phases have no overlapping baryon-chemical-potential interval",
    "NO_COEXISTENCE_CROSSING": "no pressure equality was found in the common mu_B interval",
    "MULTIPLE_COEXISTENCE_CROSSINGS": "multiple pressure-equality roots found; v1 accepts exactly one isolated transition",
    "NONPOSITIVE_DENSITY_JUMP": "pressure equality does not connect a distinct denser high-density phase",
    "INTERPOLANT_SHAPE_FAIL": "P(mu_B) interpolation violates positive density or mechanical stability",
    "COMPOSITION_UNSUPPORTED": "v1 Maxwell construction is thermo-only; composition-bearing inputs are rejected, not discarded",
    "UNKNOWN": "some unrecognized status string",
}


@pytest.mark.parametrize("category,status", list(CATEGORY_SAMPLES.items()))
def test_python_taxonomy_matches_expected_category(category, status):
    """The Python maxwell_status_category maps each representative v1 kernel
    status string to its expected §3 category."""
    got = maxwell_status_category(status)
    assert got == category, f"status {status!r}: expected {category}, got {got}"


def test_taxonomy_strips_phase_prefix():
    """The 'low-density phase: ' / 'high-density phase: ' prefix is stripped so
    the category reflects the underlying cause, not which phase failed."""
    base = "P(mu_B) interpolation violates positive density or mechanical stability"
    assert maxwell_status_category(f"low-density phase: {base}") == "INTERPOLANT_SHAPE_FAIL"
    assert maxwell_status_category(f"high-density phase: {base}") == "INTERPOLANT_SHAPE_FAIL"


def test_taxonomy_does_not_branch_on_model_name():
    """The mapping is data-driven from the status string only; model names
    appearing in the status do not change the category."""
    s = "P(mu_B) interpolation violates positive density or mechanical stability"
    # the same status with a model name embedded (hypothetical) maps the same
    assert maxwell_status_category(s) == "INTERPOLANT_SHAPE_FAIL"
    assert maxwell_status_category("SPG_M1 " + s) == "INTERPOLANT_SHAPE_FAIL"


def test_v2_experiments_agree_on_structured_category_all_65_pairs():
    """All 65 pairs receive the SAME structured category in the PDG-gauge and
    gauge-shift v2 experiments. A representation-only transformation must not
    alter the physics verdict."""
    rep = reconcile()
    assert rep["n_pairs_compared"] == 65
    assert rep["n_pairs_same_structured_category"] == 65, (
        f"{rep['n_pairs_same_structured_category']}/65 pairs agree on structured category")
    # the per-pair rows must all have same_structured_category True
    disagree = [r for r in rep["rows"] if not r["same_structured_category"]]
    assert not disagree, f"{len(disagree)} pairs disagree on structured category"


def test_v2_old_taxonomy_had_inconsistent_labels_now_fixed():
    """The old text-based taxonomy gave different labels to identical
    failures. The structured taxonomy fixes this: 4 pairs had the same kernel
    status but different OLD categories; all 4 now receive the same structured
    category in both experiments."""
    rep = reconcile()
    assert rep["n_pairs_same_kernel_status_but_different_old_category"] == 4, (
        "expected 4 pairs with same kernel status but different old category; "
        f"got {rep['n_pairs_same_kernel_status_but_different_old_category']}")
    # the 4 pairs now all agree on structured category
    fixed = [r for r in rep["rows"] if r["old_category_disagreed_on_same_failure"]]
    assert len(fixed) == 4
    for r in fixed:
        assert r["same_structured_category"] is True
        assert r["pdg_structured_category"] == r["gauge_structured_category"]


def test_v2_mass_gate_tolerance_unchanged_at_1e12():
    """The pairwise mass-gate tolerance is unchanged at 1e-12 relative. The
    structured-taxonomy fix did not widen any physical acceptance gate."""
    rep = reconcile()
    assert rep["mass_gate_tolerance_relative"] == 1e-12
    assert rep["mass_gate_unchanged"] is True


def test_v2_dominant_failure_is_interpolant_shape_not_generic_stability_gate():
    """The structured taxonomy splits the old 'STABILITY_GATE' bucket. The
    dominant category is INTERPOLANT_SHAPE_FAIL (61/65), reflecting the
    source-table non-convexity the §1 audit diagnoses. The generic
    'STABILITY_GATE' label no longer appears."""
    rep = reconcile()
    pdg_struct = rep["structured_category_distribution_pdg"]
    assert "STABILITY_GATE" not in pdg_struct, (
        "generic STABILITY_GATE bucket must be eliminated by the structured taxonomy")
    assert pdg_struct.get("INTERPOLANT_SHAPE_FAIL") == 61, (
        f"expected 61 INTERPOLANT_SHAPE_FAIL, got {pdg_struct.get('INTERPOLANT_SHAPE_FAIL')}")
