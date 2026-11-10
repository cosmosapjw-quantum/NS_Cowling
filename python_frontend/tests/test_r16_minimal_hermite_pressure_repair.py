from __future__ import annotations

import pytest

from scripts.r16_minimal_hermite_pressure_repair import (
    minimal_admissible_delta_p,
    audit_hermite_cell,
)


def test_already_admissible_pressure_increment_is_preserved_exactly():
    result = minimal_admissible_delta_p(
        mu0=2.0, n0=1.0, mu1=3.0, n1=2.0,
        delta_p_source=1.5,
    )
    assert result["was_adjusted"] is False
    assert result["delta_p_repaired"] == pytest.approx(1.5, rel=0.0, abs=0.0)
    cell = audit_hermite_cell(
        mu0=2.0, p0=0.5, n0=1.0,
        mu1=3.0, p1=0.5 + result["delta_p_repaired"], n1=2.0,
    )
    assert cell["status"] == "PASS_NATIVE_CELL_ADMISSIBLE"


def test_too_small_pressure_increment_is_clipped_to_minimal_admissible_bound():
    result = minimal_admissible_delta_p(
        mu0=2.0, n0=1.0, mu1=3.0, n1=2.0,
        delta_p_source=0.7,
    )
    assert result["was_adjusted"] is True
    assert result["delta_p_repaired"] > result["delta_p_lower"]
    assert result["delta_p_repaired"] < result["delta_p_upper"]
    cell = audit_hermite_cell(
        mu0=2.0, p0=0.5, n0=1.0,
        mu1=3.0, p1=0.5 + result["delta_p_repaired"], n1=2.0,
    )
    assert cell["status"] == "PASS_NATIVE_CELL_ADMISSIBLE"


def test_too_large_pressure_increment_is_clipped_to_maximal_admissible_bound():
    result = minimal_admissible_delta_p(
        mu0=2.0, n0=1.0, mu1=3.0, n1=2.0,
        delta_p_source=2.0,
    )
    assert result["was_adjusted"] is True
    cell = audit_hermite_cell(
        mu0=2.0, p0=0.5, n0=1.0,
        mu1=3.0, p1=0.5 + result["delta_p_repaired"], n1=2.0,
    )
    assert cell["status"] == "PASS_NATIVE_CELL_ADMISSIBLE"


def test_no_nonempty_cubic_hermite_admissibility_interval_is_rejected():
    with pytest.raises(ValueError, match="empty"):
        minimal_admissible_delta_p(
            mu0=1.0, n0=1.0, mu1=2.0, n1=1.0000000001,
            delta_p_source=1.0,
        )


def test_fortran_cell_acceptance_is_unit_consistent_in_si():
    mev_to_j = 1.602176634e-13
    mev_fm3_to_pa = 1.602176634e32
    fm3_inverse_to_m3_inverse = 1.0e45
    result = audit_hermite_cell(
        mu0=2.0 * mev_to_j,
        p0=0.5 * mev_fm3_to_pa,
        n0=1.0 * fm3_inverse_to_m3_inverse,
        mu1=3.0 * mev_to_j,
        p1=2.0 * mev_fm3_to_pa,
        n1=2.0 * fm3_inverse_to_m3_inverse,
    )
    assert result["status"] == "PASS_NATIVE_CELL_ADMISSIBLE"
    assert result["stable_interval"] is True
    assert result["causal_interval"] is True
