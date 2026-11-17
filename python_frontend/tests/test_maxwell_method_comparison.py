"""§2 regression: no interpolation method rescues the SPG_M1 Iphase-7 crust
segment within the source-error budget -> BLOCKED_SOURCE_RECONSTRUCTION_REQUIRED.

Confirms the §2 finding: the source-table non-convexity (§1) is not repairable
by changing the interpolant. Cubic Hermite, PCHIP, constrained Hermite, and the
thermodynamic-potential reconstruction from epsilon(n_B) all fail to produce a
stable, causal, thermodynamically consistent EOS on SPG_M1 Iphase-7 within the
source-error budget. The segment is BLOCKED_SOURCE_RECONSTRUCTION_REQUIRED.

The contrast case (BBB_BHF2 core, all 83 intervals source-convex) is admitted
by the existing cubic Hermite v1 kernel -- confirming the failure is specific
to the non-convex crust source, not a blanket rejection.
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

from method_comparison import compare_methods_on_segment  # noqa: E402


def test_spg_m1_iphase7_is_blocked_source_reconstruction_required():
    """All four methods fail on SPG_M1 Iphase-7 within the source budget."""
    rep = compare_methods_on_segment("SPG_M1",
                                     ROOT / "eos_generator" / "eos_data" / "SPG_M1")
    assert rep["overall_verdict"] == "BLOCKED_SOURCE_RECONSTRUCTION_REQUIRED", (
        f"SPG_M1 Iphase-7 should be BLOCKED, got {rep['overall_verdict']}")
    # each method individually rejected
    for name in ("M1_cubic_hermite", "M2_pchip_derived_n",
                 "M3_constrained_hermite", "M4_thermo_from_epsilon"):
        assert rep["methods"][name]["verdict"] == "REJECTED", (
            f"{name} should be REJECTED on SPG_M1 Iphase-7, "
            f"got {rep['methods'][name]['verdict']}")


def test_spg_m1_m4_rejected_for_P_residual_not_first_law():
    """M4's first-law closes by construction (P=n*mu-eps at nodes), but the
    reconstructed P deviates 15% from source P -- far beyond the 1e-9 budget.
    The verdict must be REJECTED for the P residual, not for first-law."""
    rep = compare_methods_on_segment("SPG_M1",
                                     ROOT / "eos_generator" / "eos_data" / "SPG_M1")
    m4 = rep["methods"]["M4_thermo_from_epsilon"]
    assert m4["verdict"] == "REJECTED"
    assert m4["first_law_max_residual_at_nodes"] < 1e-6, (
        "M4 first-law should close by construction")
    assert m4["P_residual_within_budget"] is False, (
        "M4 P residual should exceed the 1e-9 budget")
    assert m4["residuals_vs_source"]["P_max_rel"] > 1e-9, (
        f"M4 P_max_rel {m4['residuals_vs_source']['P_max_rel']} should exceed 1e-9")


def test_spg_m1_m1_cubic_hermite_violates_dndmu_on_most_intervals():
    """The cubic Hermite (v1 kernel interpolant) has negative dn/dmu on 123
    of 147 intervals -- confirming the §1 symptom is reproduced by the actual
    v1 interpolant, not just by the source-secant check."""
    rep = compare_methods_on_segment("SPG_M1",
                                     ROOT / "eos_generator" / "eos_data" / "SPG_M1")
    m1 = rep["methods"]["M1_cubic_hermite"]
    assert m1["n_intervals_violate_dndmu"] == 123, (
        f"cubic Hermite should violate dn/dmu on 123 intervals, "
        f"got {m1['n_intervals_violate_dndmu']}")


def test_bbb_bhf2_core_admitted_by_cubic_hermite():
    """The BBB_BHF2 core (all 83 intervals source-convex) is admitted by the
    cubic Hermite v1 kernel. This confirms the BLOCKED verdict on SPG_M1 is
    specific to the non-convex crust, not a blanket rejection of the method."""
    rep = compare_methods_on_segment("BBB_BHF2",
                                     ROOT / "eos_generator" / "eos_data" / "BBB_BHF2")
    assert rep["overall_verdict"] == "ADMITTED_BY_SOME_METHOD"
    assert rep["methods"]["M1_cubic_hermite"]["verdict"] == "ADMITTED"
    assert rep["n_rejected_source_convexity"] == 0, (
        "BBB_BHF2 core should have zero source-convexity failures")


def test_bbb_bhf2_core_m4_rejected_for_P_residual():
    """Even on the clean core, M4 (thermo from epsilon PCHIP) is rejected:
    the PCHIP derivative of epsilon does not reproduce the source mu_B, so
    P=n*mu_B-eps deviates 144% from source P. This proves the §2 rule: do not
    assume a PCHIP derivative reproduces the source chemical potential."""
    rep = compare_methods_on_segment("BBB_BHF2",
                                     ROOT / "eos_generator" / "eos_data" / "BBB_BHF2")
    m4 = rep["methods"]["M4_thermo_from_epsilon"]
    assert m4["verdict"] == "REJECTED"
    assert m4["P_residual_within_budget"] is False
