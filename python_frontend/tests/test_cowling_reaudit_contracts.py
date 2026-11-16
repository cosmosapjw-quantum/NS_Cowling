"""Adversarial contracts at the phase solver / serialized census boundary."""
from pathlib import Path
import copy
import importlib.util
import json
from types import SimpleNamespace

import numpy as np
import pytest

from compose_cowling.maxwell import MaxwellEOS, MEV_FM3_TO_J_M3
import compose_cowling.maxwell as maxwell
from test_maxwell_star import quadratic_phase

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("survivor_stars_audit", ROOT / "scripts/run_css_survivor_stars.py")
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)


def transition():
    low = quadratic_phase(rest_mass=939., stiffness=1., n_min=.3, n_max=1., count=8)
    high = quadratic_phase(rest_mass=940., stiffness=.25, n_min=2., n_max=2.8, count=8)
    table = copy.deepcopy(low)
    for key, value in vars(low).items():
        if isinstance(value, np.ndarray) and value.shape == low.nb.shape:
            setattr(table, key, np.concatenate((value, getattr(high, key))))
    table.Iphase = np.r_[np.ones(8, dtype=int), np.full(8, 2, dtype=int)]
    return table


def test_admitted_roundoff_pressure_mismatch_evaluates_both_interface_sides():
    table = transition()
    # Preserve the exact source Euler identity after the small P perturbation.
    table.Q1[8] *= 1 + 5e-11
    table.Q3[8] = (table.eps[8] + table.p[8]) / (table.mn * table.nb[8]) - 1
    eos = MaxwellEOS(table)
    low = eos.evaluate(eos.transition_pressure_pa, "low")
    high = eos.evaluate(eos.transition_pressure_pa, "high")
    assert low["pressure_pa"] == high["pressure_pa"]
    assert low["energy_j_m3"] < high["energy_j_m3"]


def test_pressure_jump_gate_is_relative_even_below_one_mev():
    table = transition()
    table.nb *= .01
    # Q1 and Q7 preserve per-baryon values, so P and eps scale with n_B.
    table.Q1[8] *= 1 + 2e-9
    table.Q3[8] = (table.eps[8] + table.p[8]) / (table.mn * table.nb[8]) - 1
    with pytest.raises(ValueError, match="transition pressure"):
        MaxwellEOS(table)


def test_fractional_phase_labels_are_not_truncated_to_integers():
    table = transition()
    table.Iphase = table.Iphase.astype(float) + .25
    with pytest.raises(ValueError, match="phase label"):
        MaxwellEOS(table)


def test_frequency_scan_includes_last_exact_endpoint_root(monkeypatch):
    def shot(star, frequency_hz, **kwargs):
        ones = np.ones(3)
        return SimpleNamespace(residual=frequency_hz - 3.,
            w_inner_m3=ones, w_outer_m3=ones, v_inner_m2=ones, v_outer_m2=ones,
            interface_diagnostics={})
    monkeypatch.setattr(maxwell, "shoot_maxwell_cowling", shot)
    roots = maxwell.scan_maxwell_cowling_modes(None, np.array([1., 2., 3.]))
    assert [r["frequency_hz"] for r in roots] == [3.]


def test_survivor_diagnostics_remain_machine_readable():
    result = driver._plain({"residual": np.float64(.01), "nodes": np.array([1, 2])})
    assert result == {"residual": .01, "nodes": [1, 2]}


def test_survivor_driver_reports_failure_exit_code(tmp_path):
    (tmp_path / "matching_manifest.json").write_text(json.dumps({
        "survivors": [{"pair_id": "missing", "match_report": {}}]
    }))
    target = tmp_path / "stars.json"
    assert driver.main(tmp_path, target) == 1
    record = json.loads(target.read_text())
    assert record["survivors"][0]["status"] == "ERROR"

def _survivor_files(tmp_path):
    import hashlib
    pair = tmp_path / "example" / "compose"
    pair.mkdir(parents=True)
    hashes = {}
    for name in ("eos.nb", "eos.t", "eos.yq", "eos.thermo"):
        (pair / name).write_text(name + "\n")
        hashes[name] = hashlib.sha256((pair / name).read_bytes()).hexdigest()
    survivor = {"pair_id": "example", "method": "css",
        "admission": {"status": "ADMITTED", "schema": "ns-cowling.match-report-admission.v1",
                      "method": "css", "failed_checks": [], "checks": {"example": {"pass": True}}},
        "output_sha256": hashes, "match_report": {"status": "ok"}, "source_sha256": {},
        "fortran_executable_sha256": "test-double", "readback_validation": {},
        "consistency_report": {}, "combined_n": 4, "readback_n": 4,
        "readback_status": "ok", "consistency_status": "PASS"}
    receipt = pair.parent / "pair_report.json"
    receipt.write_text(json.dumps(dict(survivor, survived=True, stage="complete")))
    survivor["pair_report_sha256"] = hashlib.sha256(receipt.read_bytes()).hexdigest()
    return survivor


def test_survivor_input_hashes_are_verified_before_stellar_solves(tmp_path):
    survivor = _survivor_files(tmp_path)
    assert driver._validated_survivor_dir(tmp_path, survivor) == (tmp_path / "example" / "compose").resolve()
    (tmp_path / "example" / "compose" / "eos.thermo").write_text("changed\n")
    with pytest.raises(ValueError, match="sha256"):
        driver._validated_survivor_dir(tmp_path, survivor)


@pytest.mark.parametrize("field,value", [("method", "maxwell"),
    ("admission", {"status": "REJECTED"}), ("output_sha256", {})])
def test_survivor_input_contract_is_enforced(tmp_path, field, value):
    survivor = _survivor_files(tmp_path)
    survivor[field] = value
    with pytest.raises(ValueError):
        driver._validated_survivor_dir(tmp_path, survivor)

@pytest.mark.parametrize("field,value", [
    ("css_transition_low_index", 7.5), ("css_transition_high_index", 8.5),
    ("n", 16.5), ("css_transition_low_index", True)])
def test_css_report_requires_exact_integer_indices(monkeypatch, field, value):
    import compose_cowling.css_adapter as adapter
    table = transition()
    report = {"construction_kind": adapter.CSS_KIND, "css_transition_low_index": 7,
              "css_transition_high_index": 8, "n": 16}
    report[field] = value
    monkeypatch.setattr(adapter, "read_compose", lambda path: table)
    with pytest.raises(ValueError, match="index"):
        adapter.phase_table_from_css_match("unused", report)

def test_survivor_receipt_is_not_overwritten(tmp_path):
    target = tmp_path / "stars.json"
    target.write_text("preserved\n")
    (tmp_path / "matching_manifest.json").write_text('{"survivors": []}')
    with pytest.raises(FileExistsError):
        driver.main(tmp_path, target)
    assert target.read_text() == "preserved\n"


def test_mass_bracket_does_not_cross_failed_or_decreasing_sequence_segments():
    a = {"p_c_MeV_fm3": 1., "mass_solar": 1.2}
    b = {"p_c_MeV_fm3": 2., "mass_solar": 1.6}
    assert driver._first_rising_bracket([a, b], 1.4) == (a, b)
    failure = {"p_c_MeV_fm3": 1.5, "error": "integration failed"}
    assert driver._first_rising_bracket([a, failure, b], 1.4) is None
    turn = {"p_c_MeV_fm3": 1.5, "mass_solar": 1.1}
    assert driver._first_rising_bracket([a, turn, b], 1.4) is None


def test_stellar_receipt_records_executed_source_identity(tmp_path):
    import hashlib
    (tmp_path / "matching_manifest.json").write_text('{"survivors": []}')
    target = tmp_path / "stars.json"
    assert driver.main(tmp_path, target) == 0
    receipt = json.loads(target.read_text())
    relative = "python_frontend/compose_cowling/css_barotrope.py"
    assert receipt["executed_source_sha256"][relative] == hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
    assert isinstance(receipt["git_tracked_tree_dirty"], bool)
    assert receipt["runtime"]["python"]


def test_mode_scan_rejects_a_sign_jump_without_a_residual_zero(monkeypatch):
    def shot(star, frequency_hz, **kwargs):
        ones = np.ones(3)
        return SimpleNamespace(residual=-1. if frequency_hz < 2.1 else 1.,
            w_inner_m3=ones, w_outer_m3=ones, v_inner_m2=ones, v_outer_m2=ones,
            interface_diagnostics={})
    monkeypatch.setattr(maxwell, "shoot_maxwell_cowling", shot)
    with pytest.raises(ArithmeticError, match="residual"):
        maxwell.scan_maxwell_cowling_modes(None, np.array([1., 2., 3.]))
