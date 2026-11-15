"""§1 regression test: SPG_M1 Iphase-7 first interval is a SOURCE_TABLE_CONVEXITY_FAIL.

This is the canonical reproduction of the Maxwell stability-gate failure. It
asserts that the first interval of the last/highest-density contiguous Iphase
segment (Iphase 7, 148 rows for SPG_M1) is diagnosed as a source-table-level
convexity inconsistency, NOT a generic Hermite-interpolant overshoot. The
diagnosis is derived from the interval data (the source secant lies above
both endpoint densities), not from the model name.

The fixture values are the ones the audit brief specified (§1):

  n0   = 0.002511886432  fm^-3
  n1   = 0.002570395783  fm^-3
  mu0  = 941.1233955964  MeV
  mu1  = 941.1475485761  MeV
  P0   = 0.0030998574392 MeV fm^-3
  P1   = 0.0031724602093 MeV fm^-3
  dP/dmu_B ~ 0.003005955  fm^-3  (above both endpoint densities)

This secant lies above max(n0, n1), so the source table itself is not convex
in (mu_B, P): no interpolant that reproduces the source endpoints can have
dn_B/dmu_B >= 0 throughout the interval. The cubic Hermite on this interval
reproduces a negative dn_B/dmu_B (min approx -0.113 fm^-3/MeV), which is the
symptom the v1 kernel currently reports as a generic "P(mu_B) interpolation
violates positive density or mechanical stability" -- conflating the
source-level root cause with an interpolant-level symptom.
"""
from __future__ import annotations
import os
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "maxwell_audit"))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_ROOT", str(ROOT))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_EXE",
                      str(ROOT / "build" / "compose_eos_fortran" / "compose_eos_fortran"))

from source_consistency_audit import audit_segment, select_last_iphase_segment  # noqa: E402
from compose_eos import read_compose  # noqa: E402

SPG_M1_DIR = ROOT / "eos_generator" / "eos_data" / "SPG_M1"


def _first_interval():
    d = read_compose(str(SPG_M1_DIR))
    il = np.asarray(d.Iphase, int)
    seg = select_last_iphase_segment(il)
    diags = audit_segment("SPG_M1", SPG_M1_DIR, seg)
    assert len(diags) >= 1, "Iphase-7 segment produced no intervals"
    return seg, diags[0]


def test_spg_m1_iphase7_segment_is_148_rows():
    seg, _ = _first_interval()
    label, lo, hi = seg
    assert label == 7, f"expected last Iphase=7, got {label}"
    assert hi - lo == 148, f"expected 148 rows in Iphase-7, got {hi - lo}"


def test_spg_m1_iphase7_first_interval_reproduces_fixture_values():
    """The first Iphase-7 interval reproduces the §1 fixture values to 10 sig figs."""
    _, iv = _first_interval()
    # Fixture values from the audit brief, §1.
    assert iv.n0 == pytest.approx(0.002511886432, rel=1e-9)
    assert iv.n1 == pytest.approx(0.002570395783, rel=1e-9)
    assert iv.mu0 == pytest.approx(941.1233955964, rel=1e-9)
    assert iv.mu1 == pytest.approx(941.1475485761, rel=1e-9)
    assert iv.P0 == pytest.approx(0.0030998574392, rel=1e-8)
    assert iv.P1 == pytest.approx(0.0031724602093, rel=1e-8)
    assert iv.secant == pytest.approx(0.003005955, rel=1e-6)


def test_spg_m1_iphase7_first_interval_secant_lies_above_endpoint_densities():
    """The source secant dP/dmu_B lies ABOVE both endpoint densities, which is
    the defining data-derived condition for SOURCE_TABLE_CONVEXITY_FAIL."""
    _, iv = _first_interval()
    assert iv.secant > iv.n_max, (
        f"secant {iv.secant} should lie above max endpoint density {iv.n_max}")
    assert iv.n_max == max(iv.n0, iv.n1)


def test_spg_m1_iphase7_first_interval_is_source_table_convexity_fail():
    """The diagnosis is SOURCE_TABLE_CONVEXITY_FAIL, not a generic Hermite
    overshoot. This is derived from the interval data, not from the model
    name (the audit code never branches on the model string)."""
    _, iv = _first_interval()
    assert iv.category == "SOURCE_TABLE_CONVEXITY_FAIL", (
        f"expected SOURCE_TABLE_CONVEXITY_FAIL, got {iv.category}: {iv.reason}")


def test_spg_m1_iphase7_first_interval_hermite_has_negative_dndmu():
    """The cubic Hermite on this interval reproduces a negative dn_B/dmu_B
    (the symptom the v1 kernel currently reports generically). This confirms
    the interval would be rejected by the interpolant-level check too, but the
    ROOT CAUSE is source-level (the secant already violates convexity)."""
    _, iv = _first_interval()
    assert iv.hermite_min_dndmu < 0, (
        f"Hermite dn/dmu minimum should be negative, got {iv.hermite_min_dndmu}")


def test_spg_m1_iphase7_first_interval_thermo_identity_passes():
    """The Q3/Euler and first-law residuals pass their tolerances at this
    interval. This proves the failure is NOT a thermo-identity or
    mass-convention problem: the source is algebraically self-consistent,
    it is just not convex in (mu_B, P)."""
    _, iv = _first_interval()
    assert iv.euler_residual_0 < 5e-4, (
        f"Q3/Euler residual at lo = {iv.euler_residual_0}")
    assert iv.euler_residual_1 < 5e-4, (
        f"Q3/Euler residual at hi = {iv.euler_residual_1}")
    assert iv.first_law_residual_0 < 1e-9, (
        f"first-law residual at lo = {iv.first_law_residual_0}")
    assert iv.first_law_residual_1 < 1e-9, (
        f"first-law residual at hi = {iv.first_law_residual_1}")


def test_spg_m1_iphase7_has_many_source_convexity_failures_not_one():
    """The convexity failure is not isolated to the first interval: 99 of 147
    intervals in the Iphase-7 segment are source-table convexity failures.
    This rules out a single bad row and confirms a systematic source-table
    non-convexity in the SPG_M1 inner-crust branch."""
    d = read_compose(str(SPG_M1_DIR))
    il = np.asarray(d.Iphase, int)
    seg = select_last_iphase_segment(il)
    diags = audit_segment("SPG_M1", SPG_M1_DIR, seg)
    cats = [dd.category for dd in diags]
    n_source_fail = cats.count("SOURCE_TABLE_CONVEXITY_FAIL")
    n_hermite_fail = cats.count("HERMITE_INTERPOLANT_SHAPE_FAIL")
    n_pass = cats.count("PASS")
    assert n_source_fail == 99, f"expected 99 source-convexity fails, got {n_source_fail}"
    assert n_hermite_fail == 24, f"expected 24 hermite-shape fails, got {n_hermite_fail}"
    assert n_pass == 24, f"expected 24 passing intervals, got {n_pass}"


def test_diagnosis_is_data_derived_not_model_name_hardcoded():
    """The audit module does not branch on the model name to assign a category.
    Verify by checking the classification if/elif chain only: no 'SPG_M1' or
    model-name string appears in the category assignment. The category is
    computed from secant_in_range and the Hermite extrema only."""
    import inspect
    import source_consistency_audit as m
    src = inspect.getsource(m.audit_segment)
    # isolate the classification if/elif chain (from 'category = "PASS"' to
    # the 'diags.append(' call that closes the per-interval block).
    i0 = src.index('category = "PASS"')
    i1 = src.index('diags.append(', i0)
    classify_block = src[i0:i1]
    # the category field and the data-derived conditions are the only allowed
    # dependencies in this block.
    for forbidden in ("SPG_M1", "model", "iphase_label", "model_dir"):
        # `model` appears only in IntervalDiag(model=model,...) at the end,
        # not in the if/elif conditions; exclude that line from the check.
        condition_lines = [ln for ln in classify_block.splitlines()
                           if "category" in ln or "if " in ln or "elif " in ln]
        for ln in condition_lines:
            assert forbidden not in ln, (
                f"category condition must not depend on '{forbidden}': {ln.strip()}")
    # the only data-derived conditions used are secant_ok, hermite extrema,
    # euler/first-law residuals.
    assert "secant_ok" in classify_block
    assert "herm_min_n" in classify_block
    assert "herm_min_dndmu" in classify_block
