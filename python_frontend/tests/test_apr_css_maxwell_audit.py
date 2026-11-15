"""§7 regression: APR + phenomenological CSS Maxwell generator audit.

Verifies:
  1. The differential first-law n_B = dP/dmu_B is tested on the APR source and
     the CSS high branch. The CSS high branch holds by construction; the APR
     source via independent PCHIP does NOT guarantee n=dP/dmu (the §7 concern).
  2. The phenomenological declaration is present (high-density phase is NOT
     a microscopic quark-matter EOS).
  3. The Seidov scope is correctly limited (small-core onset, not full radial
     stability).
  4. The representative example has common P_t/mu_t, positive density and
     energy jumps, first-law jump identity, high-side cs2 <= 1, causality,
     correct phase IDs and transition row, and TOV/Love convergence.
"""
from __future__ import annotations
import json
import os
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "maxwell_audit"))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_ROOT", str(ROOT))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_EXE",
                      str(ROOT / "build" / "compose_eos_fortran" / "compose_eos_fortran"))

from audit_apr_css_maxwell import audit_apr_css  # noqa: E402


@pytest.fixture(scope="module")
def rep():
    return audit_apr_css()


def test_apr_css_differential_first_law_measured(rep):
    """The differential first-law n=dP/dmu is measured on both branches. The
    CSS high branch holds by construction; the APR source via independent
    PCHIP does NOT guarantee it (the §7 concern is confirmed measurable)."""
    dfl = rep["1_differential_first_law_n_equals_dPdmu"]
    # the CSS high branch (analytic formula) holds to finite-difference precision
    assert dfl["css_high_branch_residual_relative"] < 1e-3, (
        f"CSS branch dP/dmu residual {dfl['css_high_branch_residual_relative']} should be < 1e-3")
    # the residual is recorded (the §7 audit does NOT assume n=dP/dmu; it tests it)
    assert "apr_source_residual_relative" in dfl


def test_apr_css_phenomenological_declaration_present(rep):
    """The high-density phase is declared phenomenological, not a microscopic
    quark-matter EOS."""
    pd = rep["2_phenomenological_declaration"]
    assert pd["present"] is True
    assert pd["is_microscopic_quark_matter_eos"] is False
    assert "phenomenological" in pd["interpretation_text"].lower()


def test_apr_css_seidov_scope_limited(rep):
    """The Seidov criterion is correctly scoped as a small-core onset
    criterion, NOT a full radial-mode stability calculation."""
    s = rep["3_seidov_scope"]
    assert "small-core onset" in s["scope"] or "small_core onset" in s["scope"]
    assert "radial-mode" in s["scope"] or "radial mode" in s["scope"]


def test_apr_css_example_common_pressure_and_chemical_potential(rep):
    ex = rep["4_representative_example"]
    assert ex["common_Pt_mu_t"] is True


def test_apr_css_example_positive_density_and_energy_jumps(rep):
    ex = rep["4_representative_example"]
    assert ex["positive_density_jump"] is True
    assert ex["positive_energy_jump"] is True


def test_apr_css_example_first_law_jump_identity(rep):
    ex = rep["4_representative_example"]
    assert ex["first_law_jump_identity"] is True
    assert ex["first_law_jump_residual"] < 1e-6


def test_apr_css_example_high_side_cs2_causal(rep):
    ex = rep["4_representative_example"]
    assert ex["high_side_cs2_causal"] is True
    assert 0 < ex["high_side_cs2"] <= 1


def test_apr_css_example_phase_ids_and_transition_row(rep):
    ex = rep["4_representative_example"]
    p = ex["phase_ids"]
    assert p["phase1_rows"] > 0
    assert p["phase2_rows"] > 0
    assert p["transition_row_index"] is not None


def test_apr_css_example_tov_love_converged(rep):
    ex = rep["4_representative_example"]
    assert ex["tov_converged"] is True
    assert ex["tov_mass_Msun"] > 0
    assert ex["tov_radius_km"] > 0
    assert ex["tov_Lambda"] > 0


def test_apr_css_verdict_is_blocked_when_differential_gate_fails(rep):
    """A failed branch-interior differential first law must block admission,
    even when the Maxwell jump identity and representative TOV solution pass."""
    dfl = rep["1_differential_first_law_n_equals_dPdmu"]
    if not dfl["apr_source_differential_first_law_holds"]:
        assert rep["verdict"] == "REVIEW_REQUIRED"
    else:
        assert rep["verdict"] == "ADMITTED_WITH_PHENOMENOLOGICAL_DECLARATION"


def test_apr_css_exact_pchip_derivative_is_an_admission_gate(rep):
    """A representative TOV/Love solution and a correct jump identity do not
    override a failed differential first-law check inside the interpolated
    phase. Use the exact PCHIP derivative, not finite differences of samples."""
    dfl = rep["1_differential_first_law_n_equals_dPdmu"]
    assert dfl["derivative_operator"] == "PchipInterpolator.derivative"
    assert dfl["apr_source_differential_first_law_holds"] is False
    assert rep["verdict"] == "REVIEW_REQUIRED"
