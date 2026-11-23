"""Fail-closed orchestration contracts for the fresh native EOS census."""
from pathlib import Path
from types import SimpleNamespace
import importlib.util
import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "eos_reaudit_runner", ROOT / "scripts/run_eos_all_pairs_20261123.py")
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


def source(identifier, role, available=True):
    return {"source_id": identifier, "role": role, "available": available,
            "status": "AVAILABLE" if available else "MISSING_INPUT",
            "matching_views": []}


def test_cartesian_census_keeps_blocked_sources_and_all_methods():
    entries = [source("c1", "core"), source("c2", "core"),
               source("k1", "crust"), source("k2", "crust", False),
               source("k3", "crust")]
    cells = runner.enumerate_cells(entries)
    assert len(cells) == 2 * 3 * 7
    assert len({c[0] for c in cells}) == len(cells)
    assert {(c[1]["source_id"], c[2]["source_id"], c[3]) for c in cells} == {
        (k, c, m) for k in ("k1", "k2", "k3") for c in ("c1", "c2")
        for m in runner.METHODS}


def test_duplicate_source_identity_fails_before_dispatch():
    with pytest.raises(ValueError, match="duplicate"):
        runner.enumerate_cells([source("x", "core"), source("x", "crust")])


def test_missing_source_is_explicit_blocker_without_physical_computation(monkeypatch):
    monkeypatch.setattr(runner, "try_candidate",
                        lambda *args: pytest.fail("must not compute unavailable source"))
    result = runner.run_cell(("cell_1", source("k", "crust", False),
                              source("c", "core"), "css"))
    assert result["status"] == "BLOCKED_SOURCE"
    assert result["native_candidate_attempts"] == 0
    assert result["failure_category"] == "data_provenance"


@pytest.mark.parametrize("returncode,stdout,accepted", [
    (0, '{"status":"ok","operation":"match-core"}', True),
    (3, '{"status":"ok","operation":"match-core"}', False),
    (0, '{"status":"error","operation":"match-core"}', False),
    (0, '{"status":"ok","operation":"write"}', False),
    (0, "unstructured success text", False),
])
def test_native_protocol_requires_agreeing_exit_status_and_operation(
        monkeypatch, returncode, stdout, accepted):
    monkeypatch.setattr(runner, "SETTINGS", {"exe": "/pinned/native", "timeout": 1})
    monkeypatch.setattr(runner.subprocess, "run",
                        lambda *a, **k: SimpleNamespace(
                            returncode=returncode, stdout=stdout, stderr=""))
    assert runner.invoke("match-core")["ok"] is accepted


def test_views_require_explicit_method_permission_and_minimum_rows():
    entry = {"matching_views": [
        {"id": "supported", "methods": ["css"], "rows": 3},
        {"id": "unapproved", "methods": [], "rows": 100},
        {"id": "short", "methods": ["css"], "rows": 2},
        {"id": "other", "methods": ["maxwell"], "rows": 10},
    ]}
    assert [v["id"] for v in runner.allowed_views(entry, "css")] == ["supported"]


def test_post_assembly_rejection_is_not_reported_as_exhaustive_no_bridge():
    call = {"payload": {"window_search_exhausted": False},
            "message": "Assembled consistency failed"}
    outcome = runner.failure("consistency", call)
    assert outcome["failure_category"] == "native_consistency_gate"
    assert "exhaustive" not in outcome["failure_category"]


def test_existing_interface_receipt_contract_is_admitted_not_generic_pass():
    receipt = {"schema": "ns-cowling.match-report-admission.v1", "method": "css",
               "status": "ADMITTED", "failed_checks": [],
               "checks": {"independent_interface": {"pass": True}}}
    assert runner.interface_admitted(receipt, "css")
    assert not runner.interface_admitted({**receipt, "status": "PASS"}, "css")
    assert not runner.interface_admitted(receipt, "maxwell")
    assert not runner.interface_admitted({**receipt, "checks": {}}, "css")
    assert not runner.interface_admitted(
        {**receipt, "checks": {"bad": {"pass": False}}}, "css")


def test_prepared_view_bounds_cannot_escape_full_selected_domain():
    entry = {"source_id": "x", "available": True, "full_selected_rows": 5,
             "matching_views": [{"start": 0, "stop": 6, "rows": 6, "methods": ["css"]}]}
    with pytest.raises(ValueError, match="escapes"):
        runner.validate_catalog([entry])
    entry["matching_views"][0].update(stop=5, rows=5)
    runner.validate_catalog([entry])


def test_private_numeric_outputs_cannot_enter_public_repository(tmp_path):
    with pytest.raises(ValueError, match="Private"):
        runner.validate_private_location(ROOT / "private", tmp_path / "public")
    with pytest.raises(ValueError, match="Private"):
        runner.validate_private_location(tmp_path / "public" / "private", tmp_path / "public")
    runner.validate_private_location(tmp_path / "private", tmp_path / "public")
