"""Counterexamples for the explicit current-generator intake."""
import copy
import json
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from current_generator_manifest import (build_manifest, validate_manifest,
                                        validate_source, source_dependency_mismatches)


def _allowance(tmp_path):
    path = tmp_path / "allowance"
    path.mkdir()
    (path / "ALLOWANCE_UNUSED.json").write_bytes(
        (ROOT / "data/current_generator_intake_20261117/ALLOWANCE_UNUSED.json").read_bytes())
    return path


def test_matrix_and_legacy_separation():
    manifest = build_manifest(ROOT)
    assert set(manifest) == {"schema_version", "base_commit", "sources", "cells"}
    assert len(manifest["cells"]) == 390
    assert sum(c["requested"] for c in manifest["cells"]) == 390
    assert sum(c["eligible"] for c in manifest["cells"]) == 1
    assert sum(c["attempted"] for c in manifest["cells"]) == 0
    assert len({c["id"] for c in manifest["cells"]}) == 390
    assert {c["method"] for c in manifest["cells"]} == {
        "cubic_F", "cubic_P", "quintic_F", "quintic_P", "css", "maxwell"}
    assert all(not c["admitted"] for c in manifest["cells"])
    crust = next(s for s in manifest["sources"] if s["id"] == "crust:SPG_M1")
    core = next(s for s in manifest["sources"] if s["id"] == "core:GM_GM1")
    assert crust["phase_segments"] == [{"iphase": 1, "row_start": 0, "row_end": 128}]
    assert core["phase_segments"] == [{"iphase": 6, "row_start": 0, "row_end": 128}]
    assert crust["generation_executable_identity"]["sha256"] == "7f0a4d3849746b5c7d3a0dc3a4adc5189819c264163a417369f69ddb401c4894"
    assert crust["generation_executable_identity"]["readback_sha256"] == "2bc5e45f6d8ccbcaf66eacc366c805918a1afb1f4d8d44f49e674f03114c1b50"
    assert core["generation_executable_identity"]["status"] == "NOT_INDEPENDENTLY_OBSERVED"
    assert core["generation_executable_identity"]["sha256"] is None
    assert crust["density_segments"][0]["row_start"] == 0
    assert crust["density_segments"][0]["row_end"] == 128
    assert crust["density_segments"][0]["n_b_min"] < crust["density_segments"][0]["n_b_max"]
    assert all(c["status"] == "HOLD" for c in manifest["cells"] if c["method"].startswith(("cubic", "quintic")))
    validate_manifest(manifest)
    from match_core_crust import CRUST_MODELS, CORE_MODELS, METHODS
    assert len(CRUST_MODELS) * len(CORE_MODELS) * len(METHODS) == 80


def test_missing_and_duplicate_cells_rejected():
    manifest = build_manifest(ROOT)
    manifest["cells"].pop()
    with pytest.raises(ValueError):
        validate_manifest(manifest)
    manifest = build_manifest(ROOT)
    manifest["cells"][1] = copy.deepcopy(manifest["cells"][0])
    with pytest.raises(ValueError):
        validate_manifest(manifest)
    manifest = build_manifest(ROOT)
    extra = next(c for c in manifest["cells"] if c["id"] == "SPG_M1__GM_GM1__maxwell")
    extra["eligible"] = True
    with pytest.raises(ValueError, match="only minimum"):
        validate_manifest(manifest)


def test_source_seals_and_tamper(tmp_path):
    manifest = build_manifest(ROOT)
    crust = next(s for s in manifest["sources"] if s["id"] == "crust:SPG_M1")
    core = next(s for s in manifest["sources"] if s["id"] == "core:GM_GM1")
    validate_source(ROOT, crust)
    validate_source(ROOT, core)
    mismatches = source_dependency_mismatches(ROOT, core)
    assert mismatches
    assert all(set(m) == {"path", "expected_sha256", "observed_sha256"} for m in mismatches)
    minimum = next(c for c in manifest["cells"] if c["id"] == "SPG_M1__GM_GM1__css")
    assert minimum["reason_code"] == "HISTORICAL_GENERATOR_DEPENDENCY_VERSION_DIFFERENCE"
    assert minimum["eligible"] is True
    bad = copy.deepcopy(crust)
    bad["request_sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_source(ROOT, bad)
    bad = copy.deepcopy(core)
    bad["artifact_sha256"]["eos.nb"] = "0" * 64
    with pytest.raises(ValueError):
        validate_source(ROOT, bad)
    bad = copy.deepcopy(core)
    bad["support_sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_source(ROOT, bad)
    bad = copy.deepcopy(crust)
    bad["request_path"] = "data/current_generator_intake_20261117/missing_request.nml"
    with pytest.raises(ValueError):
        validate_source(ROOT, bad)
    bad = copy.deepcopy(crust)
    bad["phase_segments"][0]["iphase"] = 9
    with pytest.raises(ValueError, match="segments differ"):
        validate_source(ROOT, bad)


def test_diagnostic_is_sealed_and_reuse_cannot_promote():
    manifest = build_manifest(ROOT)
    core = next(s for s in manifest["sources"] if s["id"] == "core:GM_GM1")
    assert "sound_speed.csv" in core["artifact_sha256"]
    bad = copy.deepcopy(core)
    del bad["artifact_sha256"]["sound_speed.csv"]
    with pytest.raises(ValueError):
        validate_source(ROOT, bad)
    cell = next(c for c in manifest["cells"] if c["id"] == "SPG_M1__GM_GM1__css")
    cell["reuse_mode"] = "historical"
    cell["admitted"] = True
    with pytest.raises(ValueError):
        validate_manifest(manifest)


def test_execution_one_call_and_preserved_failure(tmp_path, monkeypatch):
    import match_core_crust as matching
    manifest_path = tmp_path / "intake.json"
    manifest_path.write_text(json.dumps(build_manifest(ROOT)))
    exe = tmp_path / "native"
    exe.write_bytes(b"test executable identity")
    monkeypatch.setattr(matching, "executable", lambda: exe)
    calls = []

    def fake_pipeline(method, crust, core, output):
        calls.append((method, crust, core, output, matching.N_BRIDGE))
        return {"survived": False, "stage": "css", "matching_kernel_invoked": True,
                "error": "preserved raw failure"}

    monkeypatch.setattr(matching, "run_pipeline", fake_pipeline)
    out = tmp_path / "run"
    allowance = _allowance(tmp_path)
    monkeypatch.setattr(matching, "CURRENT_GENERATOR_ALLOWANCE_ROOT", allowance)
    receipt = matching.execute_current_generator_manifest(manifest_path, out, root=ROOT)
    assert len(calls) == 1
    assert calls[0][0] == "css" and calls[0][4] == 40
    assert calls[0][1] == ROOT / next(s["compose_path"] for s in build_manifest(ROOT)["sources"] if s["id"] == "crust:SPG_M1")
    assert receipt["attempted"] is True and receipt["admitted"] is False
    assert receipt["pipeline_invoked"] is True and receipt["matching_kernel_invoked"] is True
    assert receipt["result"]["error"] == "preserved raw failure"
    assert receipt["current_executable_sha256"] == matching.sha256_file(exe)
    assert json.loads((out / "MATCH_EXECUTION.json").read_text()) == receipt
    census = json.loads((out / "MATCH_RESULT_CENSUS.json").read_text())
    assert len(census["cells"]) == 390
    assert sum(c["attempted"] for c in census["cells"]) == 1
    assert sum(c["admitted"] for c in census["cells"]) == 0
    assert matching.N_BRIDGE == 8
    with pytest.raises(FileExistsError):
        matching.execute_current_generator_manifest(manifest_path, out, root=ROOT)
    assert len(calls) == 1
    second = matching.execute_current_generator_manifest(manifest_path, tmp_path / "other_run",
                                                         root=ROOT)
    assert second["status"] == "HOLD_ALLOWANCE_ALREADY_CLAIMED"
    assert second["attempted"] is False and second["pipeline_invoked"] is False
    assert len(calls) == 1


def test_executable_lookup_failure_does_not_attempt(tmp_path, monkeypatch):
    import match_core_crust as matching
    manifest_path = tmp_path / "intake.json"
    manifest_path.write_text(json.dumps(build_manifest(ROOT)))
    allowance = _allowance(tmp_path)
    monkeypatch.setattr(matching, "CURRENT_GENERATOR_ALLOWANCE_ROOT", allowance)
    monkeypatch.setattr(matching, "executable", lambda: (_ for _ in ()).throw(RuntimeError("no binary")))
    monkeypatch.setattr(matching, "run_pipeline", lambda *args: pytest.fail("pipeline called"))
    receipt = matching.execute_current_generator_manifest(manifest_path, tmp_path / "run", root=ROOT)
    assert receipt["execution_started"] is True and receipt["source_verified"] is True
    assert receipt["pipeline_invoked"] is False and receipt["attempted"] is False
    assert receipt["matching_kernel_invoked"] is False
    assert not (allowance / "runtime/CLAIM.json").exists()
    census = json.loads((tmp_path / "run/MATCH_RESULT_CENSUS.json").read_text())
    assert sum(c["attempted"] for c in census["cells"]) == 0


def test_pipeline_read_failure_and_native_failure_kernel_boundary(tmp_path, monkeypatch):
    import match_core_crust as matching
    exe = tmp_path / "native"
    exe.write_bytes(b"test executable identity")
    monkeypatch.setattr(matching, "executable", lambda: exe)
    for side in ("crust", "core"):
        directory = tmp_path / side
        directory.mkdir()
        for name in matching.REQUIRED_SOURCE_FILES:
            (directory / name).write_text("mock input\n")
    calls = []

    def read_failure(action, *args):
        calls.append(action)
        return {"status": "error", "message": "read refused"}

    monkeypatch.setattr(matching, "run_native", read_failure)
    early = matching.run_pipeline("css", tmp_path / "crust", tmp_path / "core",
                                  tmp_path / "early")
    assert early["stage"] == "read_crust"
    assert early["matching_kernel_invoked"] is False
    assert calls == ["read"]

    calls.clear()

    def native_failure(action, *args):
        calls.append(action)
        return {"status": "ok"} if action == "read" else {"status": "error", "message": "native refused"}

    monkeypatch.setattr(matching, "run_native", native_failure)
    failed = matching.run_pipeline("css", tmp_path / "crust", tmp_path / "core",
                                   tmp_path / "failed")
    assert failed["stage"] == "css"
    assert failed["matching_kernel_invoked"] is True
    assert calls == ["read", "read", "match-core"]


def test_execution_read_stage_does_not_count_attempt(tmp_path, monkeypatch):
    import match_core_crust as matching
    manifest_path = tmp_path / "intake.json"
    manifest_path.write_text(json.dumps(build_manifest(ROOT)))
    exe = tmp_path / "native"
    exe.write_bytes(b"test executable identity")
    monkeypatch.setattr(matching, "executable", lambda: exe)
    monkeypatch.setattr(matching, "CURRENT_GENERATOR_ALLOWANCE_ROOT", _allowance(tmp_path))
    calls = []

    def prekernel_result(*args):
        calls.append(args)
        return {"survived": False, "stage": "read_crust", "matching_kernel_invoked": False}

    monkeypatch.setattr(matching, "run_pipeline", prekernel_result)
    receipt = matching.execute_current_generator_manifest(manifest_path, tmp_path / "run", root=ROOT)
    assert len(calls) == 1
    assert receipt["pipeline_invoked"] is True and receipt["matching_kernel_invoked"] is False
    assert receipt["attempted"] is False
    census = json.loads((tmp_path / "run/MATCH_RESULT_CENSUS.json").read_text())
    assert sum(cell["attempted"] for cell in census["cells"]) == 0


def test_missing_malformed_or_changed_pinned_source_is_390_hold(tmp_path):
    copy_root = tmp_path / "repo"
    paths = [
        "data/research_matching_20261117/spg/native_129_repaired/SPG_M1",
        "data/research_matching_20261117/integration/real_pairs_v2/core/GM_GM1",
    ]
    for rel in paths:
        shutil.copytree(ROOT / rel, copy_root / rel)
    receipt = "data/research_matching_20261117/spg/native_129_repaired/RUN_RECEIPT.json"
    (copy_root / receipt).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / receipt, copy_root / receipt)
    nb = copy_root / paths[0] / "generated/eos.nb"
    original = nb.read_bytes()
    nb.unlink()
    missing = build_manifest(copy_root)
    assert len(missing["cells"]) == 390
    assert sum(c["eligible"] for c in missing["cells"]) == 0
    assert all(c["status"] == "HOLD" for c in missing["cells"])
    assert any("eos.nb" in str(ref) for ref in next(c for c in missing["cells"] if c["id"] == "SPG_M1__GM_GM1__css")["evidence_refs"])
    nb.write_bytes(original)
    dispatch = copy_root / paths[1] / "DISPATCH.json"
    original_dispatch = dispatch.read_bytes()
    dispatch.write_text("{ malformed")
    malformed = build_manifest(copy_root)
    assert len(malformed["cells"]) == 390 and sum(c["eligible"] for c in malformed["cells"]) == 0
    dispatch.write_bytes(original_dispatch)
    nb.write_bytes(original + b"\n")
    changed = build_manifest(copy_root)
    assert len(changed["cells"]) == 390 and sum(c["eligible"] for c in changed["cells"]) == 0
    assert next(c for c in changed["cells"] if c["id"] == "SPG_M1__GM_GM1__css")["reason_code"] == "BYTE_IDENTITY_MISMATCH"


def test_execution_rejects_tamper_and_eligibility(tmp_path, monkeypatch):
    import match_core_crust as matching
    monkeypatch.setattr(matching, "run_pipeline", lambda *args: pytest.fail("matcher called"))
    manifest = build_manifest(ROOT)
    cases = []
    bad = copy.deepcopy(manifest)
    bad["base_commit"] = "0" * 40
    cases.append(bad)
    bad = copy.deepcopy(manifest)
    next(c for c in bad["cells"] if c["eligible"])["eligible"] = False
    cases.append(bad)
    bad = copy.deepcopy(manifest)
    next(c for c in bad["cells"] if c["id"] == "SPG_M1__GM_GM1__maxwell")["eligible"] = True
    cases.append(bad)
    bad = copy.deepcopy(manifest)
    next(s for s in bad["sources"] if s["id"] == "core:GM_GM1")["support_sha256"] = "0" * 64
    cases.append(bad)
    for index, case in enumerate(cases):
        path = tmp_path / f"bad{index}.json"
        path.write_text(json.dumps(case))
        with pytest.raises(ValueError):
            matching.execute_current_generator_manifest(path, tmp_path / f"out{index}", root=ROOT)
        assert not (tmp_path / f"out{index}").exists()


def test_execution_exception_receipt(tmp_path, monkeypatch):
    import match_core_crust as matching
    manifest_path = tmp_path / "intake.json"
    manifest_path.write_text(json.dumps(build_manifest(ROOT)))
    exe = tmp_path / "native"
    exe.write_bytes(b"test executable identity")
    monkeypatch.setattr(matching, "executable", lambda: exe)

    def failed_match(*args):
        raise RuntimeError("first failure retained")

    monkeypatch.setattr(matching, "run_pipeline", failed_match)
    out = tmp_path / "run"
    monkeypatch.setattr(matching, "CURRENT_GENERATOR_ALLOWANCE_ROOT", _allowance(tmp_path))
    receipt = matching.execute_current_generator_manifest(manifest_path, out, root=ROOT)
    assert receipt["pipeline_invoked"] is True and receipt["attempted"] is False
    assert receipt["admitted"] is False
    assert receipt["status"] == "HOLD_EXECUTION_EXCEPTION"
    assert "first failure retained" in (out / "MATCH_EXECUTION.json").read_text()
    census = json.loads((out / "MATCH_RESULT_CENSUS.json").read_text())
    assert sum(c["attempted"] for c in census["cells"]) == 0
    assert sum(c["admitted"] for c in census["cells"]) == 0


def test_execution_admitted_result_still_publication_hold(tmp_path, monkeypatch):
    import match_core_crust as matching
    manifest_path = tmp_path / "intake.json"
    manifest_path.write_text(json.dumps(build_manifest(ROOT)))
    exe = tmp_path / "native"
    exe.write_bytes(b"test executable identity")
    monkeypatch.setattr(matching, "executable", lambda: exe)
    calls = []

    def admitted_match(*args):
        calls.append(args)
        return {"survived": True, "matching_kernel_invoked": True,
                "admission": {"status": "ADMITTED"}}

    monkeypatch.setattr(matching, "run_pipeline", admitted_match)
    out = tmp_path / "run"
    monkeypatch.setattr(matching, "CURRENT_GENERATOR_ALLOWANCE_ROOT", _allowance(tmp_path))
    receipt = matching.execute_current_generator_manifest(manifest_path, out, root=ROOT)
    census = json.loads((out / "MATCH_RESULT_CENSUS.json").read_text())
    selected = next(c for c in census["cells"] if c["id"] == "SPG_M1__GM_GM1__css")
    assert len(calls) == 1
    assert receipt["admitted"] is True and selected["admitted"] is True
    assert sum(c["attempted"] for c in census["cells"]) == 1
    assert sum(c["admitted"] for c in census["cells"]) == 1
    assert census["execution"]["publication_status"] == "HOLD_INDEPENDENT_REVIEW"
    assert json.loads(manifest_path.read_text()) == build_manifest(ROOT)
