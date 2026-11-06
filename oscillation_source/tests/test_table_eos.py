"""T14 thermodynamically coupled table representation tests.

These tests are new in T14.  They intentionally do not re-run pre-T14 tests.
"""
from pathlib import Path
import os
import pytest

from tools.runtime import CaseRunner

ROOT = Path(__file__).resolve().parents[1]
SYNTHETIC = ROOT / "fixtures/tables/hermite_consistent.csv"
COMPOSE_134 = ROOT / "data/normalized/compose_134_si.csv"


@pytest.fixture(scope="module")
def table_case(tmp_path_factory):
    runner = CaseRunner(
        ROOT,
        tmp_path_factory.mktemp("t14-table"),
        workspace=os.environ.get("NS_WORKSPACE", "solutions"),
    )
    artifact = runner.build_fpm()

    def call(case_id, parameters, timeout=20):
        return runner.run_artifact(artifact, case_id, parameters, timeout=timeout)

    return call


def test_consistent_hermite_table_is_admitted_and_recovers_same_representation(table_case):
    r = table_case("table_diagnostics", {"table_path": str(SYNTHETIC)})
    assert r["status"] == "ok"
    assert r["row_count"] == 3
    assert r["interval_count"] == 2
    assert r["mu_strict_increasing"] is True
    assert r["pressure_strict_increasing"] is True
    assert r["rejected_intervals"] == 0
    assert r["fully_admitted"] is True
    assert r["interpolation_id"] == "cubic_Hermite_P_of_mu_raw_n_v1"

    s = table_case(
        "table_state",
        {"table_path": str(SYNTHETIC), "query_kind": "mu", "query_value": 1.3e-10},
    )
    assert s["status"] == "ok"
    assert s["mu_j"] == pytest.approx(1.3e-10, rel=3e-14)
    assert s["pressure_pa"] == pytest.approx(9.0e32, rel=3e-13)
    assert s["n_b_m3"] == pytest.approx(6.0e43, rel=3e-13)
    assert s["energy_j_m3"] == pytest.approx(6.9e33, rel=4e-13)
    assert s["chi"] == pytest.approx(3.0 / 13.0, rel=4e-13)


def test_pressure_inverse_uses_the_same_hermite_polynomial(table_case):
    r = table_case(
        "table_state",
        {"table_path": str(SYNTHETIC), "query_kind": "pressure", "query_value": 9.0e32},
    )
    assert r["status"] == "ok"
    assert r["mu_j"] == pytest.approx(1.3e-10, rel=2e-12)
    assert r["pressure_pa"] == pytest.approx(9.0e32, rel=3e-12)
    assert r["n_b_m3"] == pytest.approx(6.0e43, rel=3e-12)


def test_table_queries_do_not_extrapolate(table_case):
    lo = table_case(
        "table_state",
        {"table_path": str(SYNTHETIC), "query_kind": "mu", "query_value": 1.1e-10},
    )
    assert lo["status"] == "out_of_domain"
    assert lo["result_valid"] is False
    assert "pressure_pa" not in lo

    hi = table_case(
        "table_state",
        {"table_path": str(SYNTHETIC), "query_kind": "pressure", "query_value": 4.0e33},
    )
    assert hi["status"] == "out_of_domain"
    assert hi["result_valid"] is False


def test_compose_134_native_knots_preserve_admission_obstruction(table_case):
    r = table_case("table_diagnostics", {"table_path": str(COMPOSE_134)})
    assert r["status"] == "ok"
    assert r["row_count"] == 197
    assert r["interval_count"] == 196
    assert r["mu_strict_increasing"] is True
    assert r["pressure_strict_increasing"] is True
    assert r["admitted_intervals"] == 124
    assert r["rejected_intervals"] == 72
    assert r["first_rejected_interval"] == 1
    assert r["fully_admitted"] is False
    # Independent raw-row checks, not identities of the fitted cubic itself.
    assert r["raw_energy_closure_max_rel"] == pytest.approx(1.1327080889264811e-4, rel=5e-11)
    assert r["raw_trapezoid_first_law_max_rel"] == pytest.approx(1.261585717787261e-1, rel=5e-11)


def test_unadmitted_real_table_cannot_be_promoted_to_eos_provider(table_case):
    r = table_case(
        "table_state",
        {"table_path": str(COMPOSE_134), "query_kind": "mu", "query_value": 2.0e-10},
    )
    assert r["status"] == "matching_inconsistent"
    assert r["result_valid"] is False
    assert "pressure_pa" not in r
    assert "full table" in r["detail"].lower()
