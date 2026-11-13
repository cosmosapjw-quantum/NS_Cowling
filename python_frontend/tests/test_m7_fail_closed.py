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


def test_historical_maxwell_token_is_not_a_sharp_first_order_transition():
    assert "maxwell" in METHODS
    assert METHOD_PHYSICS["maxwell"]["kind"] == "legacy_smooth_interpolation"
    assert "NPPHB" in METHOD_PHYSICS["maxwell"]["description"]
    assert METHOD_PHYSICS["maxwell"]["physical_classification"] == "UNRESOLVED"


def test_python_oracle_and_placeholder_methods_are_blocked():
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
        "method": "maxwell",
        "method_kind": "legacy_smooth_interpolation",
        "method_reference": "project-defined",
        "acceptance_node": "E390-02C",
        "output_folder": str(out),
        "status": "PENDING",
        "crust_eos_path": str(crust),
        "core_eos_path": str(core),
    }
    monkeypatch.setattr(
        m7, "run_maxwell_task",
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


def test_run_task_refuses_to_delete_or_overwrite_an_existing_output_folder(tmp_path):
    existing = tmp_path / "existing"
    existing.mkdir()
    sentinel = existing / "keep.txt"
    sentinel.write_text("immutable", encoding="utf-8")
    task = {
        "task_id": "unit-test",
        "crust_model": "C",
        "core_model": "H",
        "method": "maxwell",
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
