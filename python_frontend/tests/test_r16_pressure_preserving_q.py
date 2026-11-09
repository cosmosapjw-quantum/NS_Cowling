from __future__ import annotations

import math

import pytest

from scripts.r16_pressure_preserving_q import (
    cell_pressure_increment_bounds,
    qexp_pressure_increment,
    qexp_shape_parameter,
    qexp_shape_derivatives,
)


def test_qexponential_linear_limit_matches_q_linear_integral():
    mu0, mu1 = 2.0, 3.0
    q0, q1 = 0.5, 1.0
    h = mu1 - mu0

    value = qexp_pressure_increment(mu0, mu1, q0, q1, k=0.0)
    expected = h * (
        (mu0 + h / 2.0) * q0
        + (mu0 / 2.0 + h / 3.0) * (q1 - q0)
    )
    assert value == pytest.approx(expected, rel=0.0, abs=1e-14)


def test_qexponential_shape_parameter_matches_a_feasible_source_increment():
    mu0, mu1 = 939.56, 941.96
    q0, q1 = 2.86e-7, 2.12e-6
    lower, upper = cell_pressure_increment_bounds(mu0, mu1, q0, q1)
    target = lower + 0.30 * (upper - lower)

    result = qexp_shape_parameter(mu0, mu1, q0, q1, target)
    assert result["feasible"] is True
    assert result["k"] < 0.0
    assert result["target_pressure_increment"] == pytest.approx(target, rel=0.0, abs=1e-15)
    assert qexp_pressure_increment(mu0, mu1, q0, q1, result["k"]) == pytest.approx(
        target, rel=2e-12, abs=1e-15
    )

    slope0, slope1 = qexp_shape_derivatives(result["k"])
    assert slope0 > 0.0 and slope1 > 0.0
    assert result["chi_lower"] > 0.0 and result["chi_lower"] <= 1.0
    assert result["chi_upper"] > 0.0 and result["chi_upper"] <= 1.0


def test_qexponential_solver_uses_positive_and_negative_shapes_on_each_side_of_linear():
    mu0, mu1 = 2.0, 3.0
    q0, q1 = 0.5, 1.0
    linear = qexp_pressure_increment(mu0, mu1, q0, q1, 0.0)
    lower, upper = cell_pressure_increment_bounds(mu0, mu1, q0, q1)

    low_target = lower + 0.35 * (upper-lower)
    high_target = lower + 0.65 * (upper-lower)
    low = qexp_shape_parameter(mu0, mu1, q0, q1, low_target)
    high = qexp_shape_parameter(mu0, mu1, q0, q1, high_target)

    assert low["k"] < 0.0
    assert high["k"] > 0.0
    assert low["target_pressure_increment"] < linear < high["target_pressure_increment"]


def test_pressure_increment_outside_monotone_q_envelope_is_marked_infeasible():
    mu0, mu1 = 2.0, 3.0
    q0, q1 = 0.5, 1.0
    lower, upper = cell_pressure_increment_bounds(mu0, mu1, q0, q1)

    below = qexp_shape_parameter(mu0, mu1, q0, q1, lower * 0.9)
    above = qexp_shape_parameter(mu0, mu1, q0, q1, upper * 1.1)

    assert below["feasible"] is False
    assert below["status"] == "BELOW_MONOTONE_Q_ENVELOPE"
    assert below["k"] is None
    assert below["realizable_by_continuous_q"] is False
    assert above["feasible"] is False
    assert above["status"] == "ABOVE_MONOTONE_Q_ENVELOPE"
    assert above["k"] is None
    assert above["realizable_by_continuous_q"] is False


def test_bounded_regularization_produces_a_finite_shape_for_outside_increment():
    mu0, mu1 = 2.0, 3.0
    q0, q1 = 0.5, 1.0
    lower, upper = cell_pressure_increment_bounds(mu0, mu1, q0, q1)

    result = qexp_shape_parameter(
        mu0, mu1, q0, q1, lower * 0.9,
        regularization_fraction=1.0e-3,
    )
    assert result["feasible"] is False
    assert result["target_pressure_increment"] > lower
    assert math.isfinite(result["k"])
    assert result["k"] < 0.0


def test_constant_q_cell_cannot_match_a_different_source_pressure_increment():
    result = qexp_shape_parameter(
        2.0, 3.0, 0.5, 0.5, 1.2,
    )
    assert result["feasible"] is False
    assert result["status"] == "INFEASIBLE_CONSTANT_Q"
    assert result["k"] is None
    assert result["realizable_by_continuous_q"] is False
