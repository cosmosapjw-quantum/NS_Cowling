"""M7 acceptance gates: method ontology, native-only PASS, and immutable outputs."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "m7"))
sys.path.insert(0, str(ROOT / "python_frontend"))

import m7_batch_orchestrator as m7  # noqa: E402
from m7_manifest_generator import METHODS, METHOD_PHYSICS  # noqa: E402


def test_method_ids_distinguish_legacy_smooth_bridge_from_sharp_maxwell():
    legacy = "first_law_constrained_hermite_pressure_bridge"
    assert legacy in METHODS
    assert "maxwell" not in METHODS
    assert METHOD_PHYSICS[legacy]["kind"] == "legacy_smooth_interpolation"
    assert METHOD_PHYSICS[legacy]["method_identity"] == "LEGACY_NPPHB_FLCHPB"
    assert "first-law-constrained Hermite pressure bridge" in METHOD_PHYSICS[legacy]["description"]
    assert METHOD_PHYSICS["maxwell"]["kind"] == "sharp_maxwell_phase_transition"
    assert METHOD_PHYSICS["maxwell"]["method_identity"] == "SHARP_MAXWELL_LOCAL_NEUTRALITY"
    assert METHOD_PHYSICS["maxwell"]["physical_classification"] == "FIRST_ORDER_SHARP"


def test_python_oracle_and_placeholder_methods_are_blocked():
    with pytest.raises(m7.BackendUnavailable, match="native Fortran"):
        m7.run_first_law_constrained_hermite_pressure_bridge_task("crust", "core", "out")
    with pytest.raises(m7.BackendUnavailable, match="native Fortran"):
        m7.run_maxwell_task("crust", "core", "out")
    with pytest.raises(m7.BackendUnavailable, match="native Fortran"):
        m7.run_css_task("crust", "core", "out")
    for method in ("cubic_F", "cubic_P", "quintic_F", "quintic_P"):
        with pytest.raises(m7.BackendUnavailable, match="native Fortran"):
            m7.run_smooth_task(method, "crust", "core", "out")

def test_run_task_does_not_promote_python_results_to_pass(tmp_path, monkeypatch):
    crust = tmp_path / "crust"
    core = tmp_path / "core"
    out = tmp_path / "new-output"
    crust.mkdir()
    core.mkdir()
    out_task = {
        "task_id": "unit-test",
        "crust_model": "C",
        "core_model": "H",
        "method": "first_law_constrained_hermite_pressure_bridge",
        "method_kind": "legacy_smooth_interpolation",
        "method_reference": "project-defined",
        "acceptance_node": "E390-02C",
        "output_folder": str(out),
        "status": "PENDING",
        "crust_eos_path": str(crust),
        "core_eos_path": str(core),
    }
    monkeypatch.setattr(
        m7, "run_first_law_constrained_hermite_pressure_bridge_task",
        lambda *_: {
            "output_files": [],
            "execution_backend": "python_oracle",
            "complete_combined_eos": True,
            "compose_readback": True,
        },
    )
    receipt = m7.run_task(out_task)
    assert receipt["status"] == "BLOCKED"
    assert receipt["blocked_reason"]
    assert receipt.get("execution_backend") != "native_fortran"


def test_m7_smooth_runner_rejects_sharp_maxwell_method_id(tmp_path):
    crust = tmp_path / "crust"
    core = tmp_path / "core"
    out = tmp_path / "sharp-output"
    crust.mkdir()
    core.mkdir()
    task = {
        "task_id": "sharp-maxwell-in-smooth-matrix",
        "crust_model": "C",
        "core_model": "H",
        "method": "maxwell",
        "method_kind": "sharp_maxwell_phase_transition",
        "method_reference": "sharp local-neutrality coexistence",
        "acceptance_node": "M5",
        "output_folder": str(out),
        "status": "PENDING",
        "crust_eos_path": str(crust),
        "core_eos_path": str(core),
    }
    receipt = m7.run_task(task)
    assert receipt["status"] == "BLOCKED"
    assert "separate M5 sharp-Maxwell workflow" in receipt["blocked_reason"]


def test_m7_smooth_runner_rejects_a_historical_v1_manifest(tmp_path):
    manifest = tmp_path / "v1.json"
    manifest.write_text(
        '{"schema":"ns-cowling.eos390-m7-batch-manifest","schema_version":"1.0.0","matrix":{"methods":["maxwell"]},"tasks":[]}',
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Historical v1 manifests"):
        m7.run_batch(str(manifest))


def test_run_task_refuses_to_delete_or_overwrite_an_existing_output_folder(tmp_path):
    existing = tmp_path / "existing"
    existing.mkdir()
    sentinel = existing / "keep.txt"
    sentinel.write_text("immutable", encoding="utf-8")
    task = {
        "task_id": "unit-test",
        "crust_model": "C",
        "core_model": "H",
        "method": "first_law_constrained_hermite_pressure_bridge",
        "method_kind": "legacy_smooth_interpolation",
        "method_reference": "project-defined",
        "acceptance_node": "E390-02C",
        "output_folder": str(existing),
        "status": "PENDING",
        "crust_eos_path": str(tmp_path / "missing-crust"),
        "core_eos_path": str(tmp_path / "missing-core"),
    }
    with pytest.raises(FileExistsError):
        m7.run_task(task)
    assert sentinel.read_text(encoding="utf-8") == "immutable"
