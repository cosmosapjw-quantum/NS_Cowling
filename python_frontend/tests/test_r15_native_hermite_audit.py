from __future__ import annotations

import pytest

from scripts.r15_native_hermite_audit import audit_native_interval


def test_convex_causal_interval_returns_exact_endpoint_chi_extrema():
    # h=1, P(t)=0.5+t+0.5*t^2, mu=2+t, n=1+t, dn/dmu=1.
    result = audit_native_interval(
        mu0_j=2.0, p0_pa=0.5, n0_m3=1.0,
        mu1_j=3.0, p1_pa=2.0, n1_m3=2.0,
    )

    assert result["status"] == "PASS_NATIVE_CELL_ADMISSIBLE"
    assert result["stable_interval"] is True
    assert result["causal_interval"] is True
    assert result["d2p_low_pa_per_j2"] == pytest.approx(1.0, rel=0.0, abs=1e-14)
    assert result["d2p_high_pa_per_j2"] == pytest.approx(1.0, rel=0.0, abs=1e-14)
    assert result["chi_min"] == pytest.approx(0.5, rel=0.0, abs=1e-14)
    assert result["chi_max"] == pytest.approx(2.0 / 3.0, rel=0.0, abs=1e-14)
    assert result["epsilon_low_pa"] == pytest.approx(1.5, rel=0.0, abs=1e-14)
    assert result["epsilon_high_pa"] == pytest.approx(4.0, rel=0.0, abs=1e-14)


def test_stable_but_acausal_interval_is_rejected_by_endpoint_margin():
    # P''=1 > 0, but mu*P''-P'=-0.5 at both endpoints, so chi>1 throughout.
    result = audit_native_interval(
        mu0_j=0.5, p0_pa=0.1, n0_m3=1.0,
        mu1_j=1.5, p1_pa=1.6, n1_m3=2.0,
    )

    assert result["stable_interval"] is True
    assert result["causal_interval"] is False
    assert result["status"] == "REJECT_ACAUSAL_NATIVE_CELL"
    assert result["chi_min"] > 1.0


def test_zero_curvature_interval_is_rejected_as_not_strictly_stable():
    result = audit_native_interval(
        mu0_j=2.0, p0_pa=0.5, n0_m3=1.0,
        mu1_j=3.0, p1_pa=1.5, n1_m3=1.0,
    )

    assert result["stable_interval"] is False
    assert result["status"] == "REJECT_NONPOSITIVE_COMPRESSIBILITY"
