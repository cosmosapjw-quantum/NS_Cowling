from __future__ import annotations

import pytest

from scripts.r16_pressure_area_feasibility import audit_pressure_area_feasibility


def test_integral_inside_monotone_q_bounds_is_feasible():
    # Integral_mu0^mu1 mu dmu = 1.5; q0=1, q1=2, qbar=1.5.
    result = audit_pressure_area_feasibility(
        mu0=1.0, p0=1.0, n0=1.0,
        mu1=2.0, p1=3.25, n1=4.0,
    )
    assert result["q_low"] == pytest.approx(1.0)
    assert result["q_high"] == pytest.approx(2.0)
    assert result["q_pressure_average"] == pytest.approx(1.5)
    assert result["status"] == "FEASIBLE_MONOTONE_CAUSAL_INTERPOLANT"


def test_pressure_increment_below_lower_bound_is_infeasible():
    result = audit_pressure_area_feasibility(
        mu0=1.0, p0=1.0, n0=1.0,
        mu1=2.0, p1=1.75, n1=4.0,
    )
    assert result["status"] == "INFEASIBLE_PRESSURE_INCREMENT_BELOW_Q_BOUNDS"
    assert result["lower_margin"] < 0.0


def test_pressure_increment_above_upper_bound_is_infeasible():
    result = audit_pressure_area_feasibility(
        mu0=1.0, p0=1.0, n0=1.0,
        mu1=2.0, p1=4.25, n1=4.0,
    )
    assert result["status"] == "INFEASIBLE_PRESSURE_INCREMENT_ABOVE_Q_BOUNDS"
    assert result["upper_margin"] < 0.0


def test_nonmonotone_q_endpoints_have_no_causal_barotropic_interpolant():
    result = audit_pressure_area_feasibility(
        mu0=1.0, p0=1.0, n0=2.0,
        mu1=2.0, p1=2.5, n1=2.0,
    )
    assert result["status"] == "INFEASIBLE_NONMONOTONE_Q_ENDPOINTS"
    assert result["q_high"] < result["q_low"]
