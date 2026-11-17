"""Adversarial cold-source and receipt-binding regressions, 2026-11-17."""
from __future__ import annotations

import csv
import hashlib
import json

import pytest

from scripts import match_report_admission as adm
from scripts import run_css_survivor_stars as stars
from pipeline.tests.test_match_report_admission import (
    HEADER, META, MN, _css_case, _maxwell_case,
)


@pytest.mark.parametrize("method,factory", [("css", _css_case), ("maxwell", _maxwell_case)])
@pytest.mark.parametrize("column,value", [(2, 0.1), (5, 0.1), (6, 0.1)])
def test_cold_serialized_columns_cannot_disappear(method, factory, column, value):
    report, raw = factory()
    raw[0][column] = value
    result = getattr(adm, "admit_" + method)(report, adm.rows_from_table(raw, META), META)
    assert result["status"] == "REJECTED", result
    assert "serialized_cold_contract" in result["failed_checks"]


@pytest.mark.parametrize("method,factory", [("css", _css_case), ("maxwell", _maxwell_case)])
@pytest.mark.parametrize("field,value", [("status", "error"), ("operation", "read")])
def test_receipt_identity_must_match_native_operation(method, factory, field, value):
    report, raw = factory()
    report[field] = value
    result = getattr(adm, "admit_" + method)(report, adm.rows_from_table(raw, META), META)
    assert result["status"] == "REJECTED", result
    assert "native_report_identity" in result["failed_checks"]


def test_boolean_pressure_crossing_is_not_an_integer_count():
    report, raw = _maxwell_case()
    report["pressure_crossings"] = True
    result = adm.admit_maxwell(report, adm.rows_from_table(raw, META), META)
    assert result["status"] == "REJECTED", result


def _write_rows(tmp_path, raw, *, temperature=0.0):
    csv_path, nml = tmp_path / "eos_data.csv", tmp_path / "eos_data.nml"
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        writer.writerows(raw)
    nml.write_text(f"&metadata\nmn={MN}\nil=1\ncharge_fraction_grid=0\n"
                   f"temperature_mev={temperature}\nbeta_equilibrated=T\n/\n")
    return csv_path, nml


def test_warm_serialized_metadata_is_rejected(tmp_path):
    report, raw = _css_case()
    paths = _write_rows(tmp_path, raw, temperature=1.0)
    result = adm.admit_match_report("css", report, *paths)
    assert result["status"] == "REJECTED", result


def test_maxwell_phase_boundary_must_match_serialized_output(tmp_path):
    report, raw = _maxwell_case()
    for i, row in enumerate(raw):
        row[9] = 1 if i < report["low_phase_rows"] else 2
    paths = _write_rows(tmp_path, raw)
    assert adm.admit_match_report("maxwell", report, *paths)["status"] == "ADMITTED"
    raw[0][9] = 9  # native output has an explicit 1/2 phase contract, unlike provider IDs
    paths = _write_rows(tmp_path, raw)
    result = adm.admit_match_report("maxwell", report, *paths)
    assert result["status"] == "REJECTED", result


@pytest.mark.parametrize("tamper", ["erase", "receipt_missing", "both_fields_missing", "bridge_collision"])
def test_css_phase_retention_receipt_is_bound_to_serialized_rows(tmp_path, tamper):
    report, raw = _css_case()
    hi, core = report["css_transition_high_index"], report["css_core_bridge_index"]
    for i, row in enumerate(raw):
        row[9] = 17 if hi <= i < core else 23
    report.update(bridge_phase_id=17, output_phase_ids=[int(row[9]) for row in raw])
    assert adm.admit_match_report("css", report, *_write_rows(tmp_path, raw))["status"] == "ADMITTED"
    if tamper == "erase":
        raw[-1][9] = 0
    elif tamper == "receipt_missing":
        del report["output_phase_ids"]
    elif tamper == "both_fields_missing":
        del report["output_phase_ids"]
        del report["bridge_phase_id"]
        for row in raw:
            row[9] = 0
    else:
        report["bridge_phase_id"] = 23
    result = adm.admit_match_report("css", report, *_write_rows(tmp_path, raw))
    assert "serialized_phase_layout" in result["failed_checks"]


def _survivor(tmp_path):
    directory = tmp_path / "case" / "compose"
    directory.mkdir(parents=True)
    hashes = {}
    for name in ("eos.nb", "eos.t", "eos.yq", "eos.thermo"):
        (directory / name).write_text(name + "\n")
        hashes[name] = hashlib.sha256((directory / name).read_bytes()).hexdigest()
    return {"pair_id": "case", "method": "css", "admission": {"status": "ADMITTED"},
            "output_sha256": hashes}


def test_preserved_phase_erased_apr_baseline_cannot_enter_new_stellar_solver(monkeypatch):
    from pathlib import Path
    from compose_eos import EOSData
    import compose_cowling.css_adapter as adapter
    root = Path(__file__).resolve().parents[2] / "data/research_matching_20261117/phase_retention_repair/attempts"
    report = next(json.loads(line) for line in (root / "overbroad_guard/baseline.log").read_text().splitlines()
                  if line.startswith("{"))
    internal = root / "phase_erased_baseline/internal"
    result = adm.admit_match_report("css", report, internal / "eos_data.csv", internal / "eos_data.nml")
    assert "serialized_phase_layout" in result["failed_checks"]
    data = EOSData.from_internal_directory(internal)
    assert not any(data.Iphase)  # actual old phase-erased output, unchanged
    monkeypatch.setattr(adapter, "read_compose", lambda _: data)
    with pytest.raises(ValueError, match="phase layout"):
        adapter.css_match_maxwell_eos(internal, report)


def test_status_only_survivor_cannot_replace_pair_receipt(tmp_path):
    survivor = _survivor(tmp_path)
    with pytest.raises(ValueError, match="receipt|admission"):
        stars._validated_survivor_dir(tmp_path, survivor)


def test_survivor_cannot_change_match_report_after_pair_seal(tmp_path):
    survivor = _survivor(tmp_path)
    survivor.update(match_report={"status": "ok", "method": "css"})
    receipt = dict(survivor, survived=True, stage="complete")
    path = tmp_path / "case" / "pair_report.json"
    path.write_text(json.dumps(receipt))
    survivor["pair_report_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    survivor["match_report"] = {"status": "error", "method": "css"}
    with pytest.raises(ValueError, match="receipt|admission"):
        stars._validated_survivor_dir(tmp_path, survivor)


def test_direct_pipeline_missing_source_is_not_a_physics_attempt(tmp_path, monkeypatch):
    from scripts import match_core_crust as pipeline
    monkeypatch.setattr(pipeline, "executable", lambda: (_ for _ in ()).throw(RuntimeError("must not execute")))
    result = pipeline.run_pipeline("css", tmp_path / "absent_crust", tmp_path / "absent_core", tmp_path / "attempt")
    assert result["failure_category"] == "data_provenance", result
    assert result["outcome"] == "NOT_EVALUATED_SOURCE_INCOMPLETE", result
    assert not (tmp_path / "attempt").exists()


def test_actual_generated_tf_phase_support_rejected_as_input_contract(tmp_path):
    from pathlib import Path
    from scripts import match_core_crust as pipeline
    root = Path(__file__).resolve().parents[2]
    crust = root / "data/research_native_20261116/public_cohort/GPPVA_FSU2R/attempt_001/output/generated"
    core = root / "eos_generator/eos_data_normalized/GM_GM1"
    result = pipeline.run_pipeline("css", crust, core, tmp_path / "tf_pair")
    assert result["survived"] is False, result
    assert result["stage"] == "source_support", result
    assert result["failure_category"] == "input_contract", result
    assert "cross-phase" in result["error"], result
    support = result["source_scope"]["crust"]
    assert support["declaration"]["cross_phase_interpolation"] == "UNAPPROVED"
    assert len(support["support_sha256"]) == 64
