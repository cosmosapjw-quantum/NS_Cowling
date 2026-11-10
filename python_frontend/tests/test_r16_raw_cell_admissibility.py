from __future__ import annotations

import pytest

from scripts.r16_raw_cell_admissibility import audit_hermite_cell


def test_strictly_stable_and_causal_hermite_cell_passes():
    result = audit_hermite_cell(
        mu0=2.0, p0=0.5, n0=1.0,
        mu1=3.0, p1=2.0, n1=2.0,
    )
    assert result["stable_interval"] is True
    assert result["causal_interval"] is True
    assert result["status"] == "PASS_NATIVE_CELL_ADMISSIBLE"
    assert result["d2p_low"] == pytest.approx(1.0)
    assert result["d2p_high"] == pytest.approx(1.0)
    assert result["causal_margin_low"] == pytest.approx(1.0)
    assert result["causal_margin_high"] == pytest.approx(1.0)


def test_stable_but_acausal_cell_is_rejected():
    result = audit_hermite_cell(
        mu0=0.5, p0=0.1, n0=1.0,
        mu1=1.5, p1=1.6, n1=2.0,
    )
    assert result["stable_interval"] is True
    assert result["causal_interval"] is False
    assert result["status"] == "REJECT_ACAUSAL_NATIVE_CELL"
    assert result["causal_margin_low"] == pytest.approx(-0.5)
    assert result["causal_margin_high"] == pytest.approx(-0.5)


def test_negative_curvature_cell_is_rejected_before_causality():
    result = audit_hermite_cell(
        mu0=2.0, p0=0.5, n0=2.0,
        mu1=3.0, p1=1.5, n1=1.0,
    )
    assert result["stable_interval"] is False
    assert result["causal_interval"] is False
    assert result["status"] == "REJECT_NONPOSITIVE_COMPRESSIBILITY"
    assert min(result["d2p_low"], result["d2p_high"]) < 0.0


def test_cell_requires_ordered_positive_finite_endpoints():
    with pytest.raises(ValueError, match="strictly ordered"):
        audit_hermite_cell(
            mu0=2.0, p0=0.5, n0=1.0,
            mu1=2.0, p1=2.0, n1=2.0,
        )
    with pytest.raises(ValueError, match="finite and positive"):
        audit_hermite_cell(
            mu0=2.0, p0=0.5, n0=1.0,
            mu1=3.0, p1=float("nan"), n1=2.0,
        )
