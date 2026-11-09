from __future__ import annotations

import numpy as np
import pytest

from scripts.r16_projection_audit import hermite_cell_diagnostics
from scripts.r17_joint_q_fit import fit_monotone_q


def _nonconvex_pressure_fixture():
    mu = np.asarray([1.0, 2.0, 3.0, 4.0, 5.0])
    n = mu.copy()  # q=n/mu is already positive and nondecreasing.
    p = np.asarray([0.2, 1.4, 2.0, 4.6, 7.0])  # adjacent secants decrease.
    return mu, n, p


def test_lambda_one_is_the_existing_q_preserving_exact_integral_baseline():
    mu, n, p = _nonconvex_pressure_fixture()
    result = fit_monotone_q(mu, n, p, lambda_density=1.0)

    assert result["q_monotone"] is True
    assert result["positive_thermodynamics"] is True
    assert result["pressure_strictly_increasing"] is True
    assert result["energy_density_strictly_increasing"] is True
    assert result["n_fitted"] == pytest.approx(n, rel=0.0, abs=1e-13)
    assert result["density_rms_relative_change"] == pytest.approx(0.0, abs=1e-13)
    assert result["objective_after_fit"] == pytest.approx(0.0, abs=1e-13)


def test_pressure_density_tradeoff_improves_pressure_fit_without_losing_cell_admissibility():
    mu, n, p = _nonconvex_pressure_fixture()
    baseline = fit_monotone_q(mu, n, p, lambda_density=1.0)
    candidate = fit_monotone_q(mu, n, p, lambda_density=0.5)

    assert candidate["objective_after_fit"] < candidate["objective_before_source_q"]
    assert candidate["pressure_rms_relative_change"] < baseline["pressure_rms_relative_change"]
    assert candidate["density_rms_relative_change"] > 0.0
    assert candidate["q_monotone"] is True
    assert candidate["positive_thermodynamics"] is True

    for i in range(len(mu) - 1):
        cell = hermite_cell_diagnostics(
            float(mu[i]), float(candidate["p_fitted"][i]), float(candidate["n_fitted"][i]),
            float(mu[i + 1]), float(candidate["p_fitted"][i + 1]), float(candidate["n_fitted"][i + 1]),
        )
        assert cell["stable"] is True
        # A q plateau has exactly zero analytic causal margin; allow endpoint
        # Hermite roundoff at 1e-12 relative to the endpoint density scale.
        causal_roundoff = 1.0e-12 * max(
            abs(float(candidate["n_fitted"][i])),
            abs(float(candidate["n_fitted"][i + 1])),
        )
        assert cell["min_causal_margin"] >= -causal_roundoff


def test_joint_fit_rejects_invalid_objective_weight_and_source_domain():
    mu, n, p = _nonconvex_pressure_fixture()
    with pytest.raises(ValueError, match="lambda_density"):
        fit_monotone_q(mu, n, p, lambda_density=0.0)
    with pytest.raises(ValueError, match="strictly increasing"):
        fit_monotone_q([1.0, 1.0, 2.0], [1.0, 2.0, 3.0], [0.2, 1.4, 2.0], lambda_density=0.5)


def test_log_density_weighted_fit_is_a_declared_coordinate_sensitivity():
    mu, n, p = _nonconvex_pressure_fixture()
    result = fit_monotone_q(
        mu, n, p, lambda_density=0.5, weight_coordinate="log_n"
    )

    assert result["weight_coordinate"] == "log_n"
    assert result["q_monotone"] is True
    assert result["positive_thermodynamics"] is True
    assert result["objective_after_fit"] <= result["objective_before_source_q"] + 1e-12


def test_joint_fit_reports_scale_relative_projected_gradient_residual():
    mu, n, p = _nonconvex_pressure_fixture()
    result = fit_monotone_q(mu, n, p, lambda_density=0.01, tolerance=1e-12)

    assert result["projected_gradient_kkt_residual_relative"] <= 1e-12
    assert result["projected_gradient_kkt_residual_relative"] == pytest.approx(
        result["projected_gradient_kkt_residual"] /
        (1.0 + max(abs(x) for x in result["q_fitted"]) / np.median(n / mu)),
        rel=1e-12,
        abs=1e-15,
    )
