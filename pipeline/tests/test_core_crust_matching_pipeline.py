from __future__ import annotations

import csv
from pathlib import Path

import pytest

from scripts import match_core_crust as pipeline


def _write_internal_csv(
    path: Path, *, q1_second: float = 2.0, yn_present_second: int = 0
) -> None:
    header = [
        "nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "q_extra",
        "Iphase", "Yn", "Yp", "Ye", "Ymu", "Abar", "Zbar", "Yh",
        "row_status", "yn_present", "yp_present", "ye_present", "ymu_present",
        "abar_present", "zbar_present", "yh_present", "c3_applicable",
        "c4_applicable",
    ]
    first = [
        0.1, 1.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 7,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
    ]
    second = [
        0.2, q1_second, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 6,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,
        yn_present_second, 0, 0, 0, 0, 0, 0, 0, 0,
    ]
    assert len(header) == len(first) == len(second)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerow(first)
        writer.writerow(second)


def test_internal_readback_comparison_detects_changed_thermodynamic_column(tmp_path):
    expected = tmp_path / "expected.csv"
    actual = tmp_path / "actual.csv"
    _write_internal_csv(expected, q1_second=2.0)
    _write_internal_csv(actual, q1_second=2.0001)

    report = pipeline.verify_internal_readback(
        expected, actual, expected_n=2, actual_n=2
    )

    assert report["status"] == "FAIL"
    assert any("Q1" in error for error in report["errors"])
    assert report["max_abs_residuals"]["Q1"] > 0.0


def test_internal_readback_comparison_checks_row_count_before_contents(tmp_path):
    expected = tmp_path / "expected.csv"
    actual = tmp_path / "actual.csv"
    _write_internal_csv(expected)
    _write_internal_csv(actual)

    report = pipeline.verify_internal_readback(
        expected, actual, expected_n=2, actual_n=1
    )

    assert report["status"] == "FAIL"
    assert any("row count" in error.lower() for error in report["errors"])


def test_internal_readback_comparison_checks_presence_masks_exactly(tmp_path):
    expected = tmp_path / "expected.csv"
    actual = tmp_path / "actual.csv"
    _write_internal_csv(expected, yn_present_second=0)
    _write_internal_csv(actual, yn_present_second=1)

    report = pipeline.verify_internal_readback(
        expected, actual, expected_n=2, actual_n=2
    )

    assert report["status"] == "FAIL"
    assert any("yn_present" in error for error in report["errors"])


def test_internal_readback_comparison_rejects_nonfinite_payload(tmp_path):
    expected = tmp_path / "expected.csv"
    actual = tmp_path / "actual.csv"
    _write_internal_csv(expected)
    _write_internal_csv(actual, q1_second=float("nan"))

    report = pipeline.verify_internal_readback(
        expected, actual, expected_n=2, actual_n=2
    )

    assert report["status"] == "FAIL"
    assert any("non-finite" in error for error in report["errors"])


def _write_internal_metadata(
    path: Path, *, mn: float = 939.0, beta_equilibrated: bool = True
) -> None:
    path.write_text(
        "&metadata\n"
        f" mn={mn:.16E}\n"
        " mp=939.0000000000000E+00\n"
        " il=1\n iT=1\n iYq=1\n inb_start=1\n"
        " nadd=0\n nadd_max=0\n naxis_t=1\n naxis_yq=1\n"
        " temperature_mev=0.0000000000000E+00\n"
        " charge_fraction_grid=0.0000000000000E+00\n"
        " composition_file_present=F\n"
        " composition_available=F\n"
        f" beta_equilibrated={'T' if beta_equilibrated else 'F'}\n"
        "/\n",
        encoding="utf-8",
    )


def test_internal_metadata_readback_detects_rest_mass_change(tmp_path):
    expected = tmp_path / "expected.nml"
    actual = tmp_path / "actual.nml"
    _write_internal_metadata(expected, mn=939.0)
    _write_internal_metadata(actual, mn=939.1)

    report = pipeline.verify_internal_metadata_readback(expected, actual)

    assert report["status"] == "FAIL"
    assert any("metadata mn" in error for error in report["errors"])


def test_internal_metadata_readback_checks_discrete_flags_and_indices(tmp_path):
    expected = tmp_path / "expected.nml"
    actual = tmp_path / "actual.nml"
    _write_internal_metadata(expected, beta_equilibrated=True)
    _write_internal_metadata(actual, beta_equilibrated=False)

    report = pipeline.verify_internal_metadata_readback(expected, actual)

    assert report["status"] == "FAIL"
    assert any("beta_equilibrated" in error for error in report["errors"])


def test_census_main_refuses_nonempty_output_root_without_deleting_existing_files(tmp_path, monkeypatch):
    occupied = tmp_path / "previous_census"
    occupied.mkdir()
    sentinel = occupied / "previous_manifest.json"
    sentinel.write_text('{"preserve": true}', encoding="utf-8")
    monkeypatch.setattr(pipeline, "OUT_ROOT", occupied)

    with pytest.raises(FileExistsError, match="not empty"):
        pipeline.main()

    assert sentinel.read_text(encoding="utf-8") == '{"preserve": true}'


def test_pair_runner_preserves_previous_work_and_output(tmp_path, monkeypatch):
    out = tmp_path / "pair"
    (out / "work").mkdir(parents=True)
    (out / "compose").mkdir()
    trace = out / "work" / "previous_trace.txt"
    trace.write_text("preserve this evidence")
    output = out / "compose" / "eos.thermo"
    output.write_text("previous eos")
    with pytest.raises(FileExistsError, match="not empty"):
        pipeline.run_pipeline("css", tmp_path / "crust", tmp_path / "core", out)
    assert trace.read_text() == "preserve this evidence"
    assert output.read_text() == "previous eos"


def test_pair_runner_rejects_unknown_method_before_writing(tmp_path):
    out = tmp_path / "pair"
    with pytest.raises(ValueError, match="unsupported matching method"):
        pipeline.run_pipeline("typo", tmp_path / "crust", tmp_path / "core", out)
    assert not out.exists()


def test_pair_runner_classifies_missing_backend_as_environment(tmp_path, monkeypatch):
    from compose_eos.native import NativeExecutionError

    def unavailable():
        raise NativeExecutionError("compiler unavailable")

    monkeypatch.setattr(pipeline, "executable", unavailable)
    result = pipeline.run_pipeline("css", tmp_path / "crust", tmp_path / "core", tmp_path / "pair")
    assert result["survived"] is False
    assert result["stage"] == "executable"
    assert result["failure_category"] == "environment"


def test_pipeline_does_not_mark_css_survived_when_readback_count_differs(tmp_path, monkeypatch):
    exe = tmp_path / "compose_eos_fortran"
    exe.write_bytes(b"test executable identity")
    monkeypatch.setattr(pipeline, "executable", lambda: exe)

    crust_dir = tmp_path / "crust"
    core_dir = tmp_path / "core"
    out_dir = tmp_path / "out"
    crust_dir.mkdir()
    core_dir.mkdir()

    actions = iter([
        {"status": "ok", "n": 3},
        {"status": "ok", "n": 4},
        {"status": "ok", "n": 6, "construction_kind": "css_single_speed_with_crust_density_jump"},
        {"status": "ok"},
        {"status": "ok", "n": 5},
    ])
    monkeypatch.setattr(pipeline, "run_native", lambda action, *args: next(actions))

    result = pipeline.run_pipeline("css", crust_dir, core_dir, out_dir)

    assert result["survived"] is False
    assert result["stage"] == "readback"
    assert result["readback_validation"]["status"] == "FAIL"


def test_pipeline_rejects_native_consistency_failure_after_verified_readback(tmp_path, monkeypatch):
    exe = tmp_path / "compose_eos_fortran"
    exe.write_bytes(b"test executable identity")
    monkeypatch.setattr(pipeline, "executable", lambda: exe)
    monkeypatch.setattr(
        pipeline, "verify_internal_readback",
        lambda *args, **kwargs: {
            "status": "PASS", "errors": [],
            "max_abs_residuals": {}, "max_rel_residuals": {},
        },
    )
    monkeypatch.setattr(
        pipeline, "verify_internal_metadata_readback",
        lambda *args, **kwargs: {
            "status": "PASS", "errors": [],
            "max_abs_residuals": {}, "max_rel_residuals": {},
        },
    )

    crust_dir = tmp_path / "crust"
    core_dir = tmp_path / "core"
    out_dir = tmp_path / "out"
    crust_dir.mkdir()
    core_dir.mkdir()

    def fake_run_native(action, *args):
        if action == "match-core":
            return {"status": "ok", "n": 6,
                    "construction_kind": "css_single_speed_with_crust_density_jump"}
        if action == "read":
            return {"status": "ok", "n": 6}
        if action == "verify":
            Path(args[1]).write_text(
                '{"schema":"compose-eos-fortran-verification-v1",'
                '"C1":{"pass":true},"C2":{"pass":true},'
                '"C6":{"pass":true},"C10":{"pass":false},'
                '"C11":{"pass":true},"failed_checks":"C10"}',
                encoding="utf-8",
            )
            return {"status": "ok", "overall": "FAIL", "failed_checks": "C10"}
        return {"status": "ok"}

    monkeypatch.setattr(pipeline, "run_native", fake_run_native)
    result = pipeline.run_pipeline("css", crust_dir, core_dir, out_dir)

    assert result["survived"] is False
    assert result["stage"] == "consistency"
    assert result["consistency_status"] == "FAIL"
    assert "C10" in result["error"]


def test_css_runner_roundtrips_real_native_benchmark_pair(tmp_path):
    """Exercise actual Fortran match -> CompOSE write -> native readback."""
    root = Path(__file__).resolve().parents[2]
    crust = root / "eos_generator" / "eos_data_normalized" / "GPPVA_FSU2H"
    core = root / "eos_generator" / "eos_data_normalized" / "GM_GM1"
    required = (crust / "eos.nb", crust / "eos.thermo",
                core / "eos.nb", core / "eos.thermo")
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        pytest.skip("native CSS runner benchmark inputs unavailable: " + ", ".join(missing))
    try:
        exe = pipeline.executable()
    except (OSError, RuntimeError, FileNotFoundError) as exc:
        pytest.skip(f"native Fortran executable unavailable: {exc}")
    if not Path(exe).is_file():
        pytest.skip(f"native Fortran executable does not exist: {exe}")

    result = pipeline.run_pipeline("css", crust, core, tmp_path / "css_roundtrip")

    assert result["survived"] is True, result
    assert result["stage"] == "complete"
    assert result["readback_status"] == "ok"
    assert result["combined_n"] == result["readback_n"]
    assert result["readback_validation"]["status"] == "PASS"
    assert result["readback_validation"]["expected_rows"] == result["combined_n"]
    assert result["readback_validation"]["actual_rows"] == result["readback_n"]
    assert result["readback_validation"]["metadata_status"] == "PASS"
    assert result["consistency_status"] == "PASS"
    consistency = result["consistency_report"]
    assert consistency["C1"]["pass"] is True
    assert consistency["C2"]["pass"] is True
    assert consistency["C6"]["pass"] is True
    assert consistency["C10"]["pass"] is True
    assert consistency["C11"]["pass"] is True
    metadata = result["readback_validation"]["metadata"]
    assert metadata["expected_values"]["naxis_t"] == 1
    assert metadata["expected_values"]["naxis_yq"] == 1
    assert metadata["actual_values"]["naxis_t"] == 1
    assert metadata["actual_values"]["naxis_yq"] == 1
    assert result["output_sha256"]
    assert result["source_sha256"]["crust"]
    assert result["source_sha256"]["core"]
    admission = result["admission"]
    assert admission["status"] == "ADMITTED", admission["failed_checks"]
    for name in ("construction_kind", "css_sound_speed_range", "css_single_sound_speed",
                 "crust_transition_pressure", "crust_transition_euler_mu",
                 "crust_transition_positive_jump", "crust_transition_latent_identity",
                 "css_segment_single_speed", "css_core_continuity",
                 "first_law_relative_residual", "report_matches_rows"):
        assert admission["checks"][name]["pass"] is True, name


def _passing_readback(monkeypatch):
    ok = {"status": "PASS", "errors": [], "max_abs_residuals": {}, "max_rel_residuals": {}}
    monkeypatch.setattr(pipeline, "verify_internal_readback", lambda *a, **k: dict(ok, errors=[]))
    monkeypatch.setattr(pipeline, "verify_internal_metadata_readback", lambda *a, **k: dict(ok, errors=[]))


def _fake_native_all_ok(action, *args):
    if action == "verify":
        import json
        Path(args[1]).write_text(json.dumps(_valid_consistency_report()), encoding="utf-8")
        return {"status": "ok", "overall": "PASS"}
    if action in ("match-core", "maxwell-transition"):
        return {"status": "ok", "n": 6,
                "construction_kind": "css_single_speed_with_crust_density_jump"}
    return {"status": "ok", "n": 6}


def _valid_consistency_report():
    return {
        "schema": "compose-eos-fortran-verification-v1", "n": 6,
        "overall": {"status": "PASS", "pass": True}, "failed_checks": "",
        **{name: {"pass": True} for name in ("C1", "C2", "C6", "C10", "C11")},
        **{name: {"status": "NOT_AVAILABLE", "pass": None} for name in ("C3", "C4")},
        "C3_applicable_points": 0, "C4_applicable_points": 0,
        "C5": {"applicable": False, "pass": None},
    }


@pytest.mark.parametrize("corruption", ["failed_gate", "missing_gate", "wrong_count", "false_unavailable"])
def test_pipeline_checks_consistency_file_despite_pass_stdout(tmp_path, monkeypatch, corruption):
    import json
    exe = tmp_path / "exe"
    exe.write_bytes(b"x")
    monkeypatch.setattr(pipeline, "executable", lambda: exe)
    _passing_readback(monkeypatch)
    monkeypatch.setattr(pipeline, "admit_match_report", lambda *a, **k: {"status": "ADMITTED"})

    def native(action, *args):
        if action != "verify":
            return _fake_native_all_ok(action, *args)
        report = _valid_consistency_report()
        if corruption == "failed_gate":
            report["C10"]["pass"] = False
        elif corruption == "missing_gate":
            del report["C11"]
        elif corruption == "wrong_count":
            report["n"] = 5
        else:
            report["C3_applicable_points"] = 1
        Path(args[1]).write_text(json.dumps(report))
        return {"status": "ok", "overall": "PASS"}

    monkeypatch.setattr(pipeline, "run_native", native)
    result = pipeline.run_pipeline("css", tmp_path / "crust", tmp_path / "core", tmp_path / "out")
    assert result["survived"] is False
    assert result["stage"] == "consistency"
    assert result["consistency_status"] == "INVALID_REPORT"


def test_pipeline_requires_independent_admission_after_native_pass(tmp_path, monkeypatch):
    exe = tmp_path / "exe"
    exe.write_bytes(b"x")
    monkeypatch.setattr(pipeline, "executable", lambda: exe)
    _passing_readback(monkeypatch)
    monkeypatch.setattr(pipeline, "run_native", _fake_native_all_ok)
    for d in ("crust", "core"):
        (tmp_path / d).mkdir()

    result = pipeline.run_pipeline("css", tmp_path / "crust", tmp_path / "core", tmp_path / "out")

    # Native status ok + verifier PASS alone must not create a survivor.
    assert result["survived"] is False
    assert result["stage"] == "admission"
    assert result["admission"]["status"] == "REJECTED"
    assert result["failure_category"] == "admission"


def test_pipeline_admits_only_when_admission_layer_admits(tmp_path, monkeypatch):
    exe = tmp_path / "exe"
    exe.write_bytes(b"x")
    monkeypatch.setattr(pipeline, "executable", lambda: exe)
    _passing_readback(monkeypatch)
    monkeypatch.setattr(pipeline, "run_native", _fake_native_all_ok)
    monkeypatch.setattr(pipeline, "admit_match_report",
                        lambda *a, **k: {"status": "ADMITTED", "failed_checks": []})
    for d in ("crust", "core"):
        (tmp_path / d).mkdir()
    result = pipeline.run_pipeline("css", tmp_path / "crust", tmp_path / "core", tmp_path / "out")
    assert result["survived"] is True
    assert result["admission"]["status"] == "ADMITTED"


def test_native_failure_keeps_stage_and_category(tmp_path, monkeypatch):
    from compose_eos.native import NativeExecutionError

    exe = tmp_path / "exe"
    exe.write_bytes(b"x")
    monkeypatch.setattr(pipeline, "executable", lambda: exe)

    def fake(action, *args):
        if action == "maxwell-transition":
            raise NativeExecutionError(
                "phase tables must use the same neutron and proton rest-mass conventions",
                category="MASS_CONVENTION_FAIL", operation=action)
        return {"status": "ok", "n": 5}

    monkeypatch.setattr(pipeline, "run_native", fake)
    for d in ("crust", "core"):
        (tmp_path / d).mkdir()
    result = pipeline.run_pipeline("maxwell", tmp_path / "crust", tmp_path / "core", tmp_path / "out")
    assert result["survived"] is False
    assert result["stage"] == "maxwell"
    assert result["native_category"] == "MASS_CONVENTION_FAIL"
    assert result["failure_category"] == "implementation_convention"


def test_census_counts_missing_sources_as_provenance_and_preserves_failed_reports(
        tmp_path, monkeypatch):
    crust_root, core_root, out_root = tmp_path / "crust", tmp_path / "core", tmp_path / "runs" / "r1"
    for name in ("eos.nb", "eos.t", "eos.yq", "eos.thermo"):
        (crust_root / "COMPLETE").mkdir(parents=True, exist_ok=True)
        (crust_root / "COMPLETE" / name).write_text("x", encoding="utf-8")
        (core_root / "CORE").mkdir(parents=True, exist_ok=True)
        (core_root / "CORE" / name).write_text("x", encoding="utf-8")
    (crust_root / "PARTIAL").mkdir()
    (crust_root / "PARTIAL" / "eos.thermo").write_text("x", encoding="utf-8")
    monkeypatch.setattr(pipeline, "CRUST_MODELS", ["COMPLETE", "PARTIAL", "ABSENT"])
    monkeypatch.setattr(pipeline, "CORE_MODELS", ["CORE"])
    monkeypatch.setattr(pipeline, "CRUST_ROOT", crust_root)
    monkeypatch.setattr(pipeline, "CORE_ROOT", core_root)
    monkeypatch.setattr(pipeline, "OUT_ROOT", out_root)
    monkeypatch.setattr(pipeline, "executable", lambda: tmp_path / "exe")

    calls = []

    def fake_pipeline(method, crust_dir, core_dir, out_dir):
        calls.append((crust_dir.name, method))
        (out_dir / "work").mkdir(parents=True, exist_ok=True)
        (out_dir / "work" / "trace.txt").write_text("kept", encoding="utf-8")
        return {"survived": False, "method": method, "stage": method,
                "error": "no pressure equality was found in the common mu_B interval",
                "failure_category": "physical"}

    monkeypatch.setattr(pipeline, "run_pipeline", fake_pipeline)
    assert pipeline.main() == 0

    import json
    manifest = json.loads((out_root / "matching_manifest.json").read_text(encoding="utf-8"))
    assert manifest["configured_attempts"] == 6
    assert manifest["evaluated_attempts"] == 2
    assert manifest["not_evaluated_source_incomplete"] == 4
    assert sorted(calls) == [("COMPLETE", "css"), ("COMPLETE", "maxwell")]
    cats = manifest["failures_by_category"]
    assert cats["data_provenance"] == 4
    assert cats["physical"] == 2
    partial = [r for r in manifest["all_results"] if r["crust_model"] == "PARTIAL"]
    assert all(r["stage"] == "source" and r["outcome"] == "NOT_EVALUATED_SOURCE_INCOMPLETE"
               for r in partial)
    assert "eos.nb" in partial[0]["missing_source_files"]
    # Failed attempts are preserved with a per-pair JSON report and hash.
    for r in manifest["all_results"]:
        report = out_root / r["pair_id"] / "pair_report.json"
        assert report.is_file()
        assert r["pair_report_sha256"] == pipeline.sha256_file(report)
    assert (out_root / "COMPLETE__CORE__css" / "work" / "trace.txt").read_text() == "kept"
    assert manifest["by_method"]["css"]["evaluated"] == 1
    assert manifest["by_method"]["maxwell"]["survivors"] == 0


def test_census_flags_byte_identical_source_models(tmp_path, monkeypatch):
    crust_root, core_root, out_root = tmp_path / "crust", tmp_path / "core", tmp_path / "run"
    for model, payload in (("A", "same"), ("B", "same"), ("C", "different")):
        (crust_root / model).mkdir(parents=True)
        for name in ("eos.nb", "eos.t", "eos.yq", "eos.thermo"):
            (crust_root / model / name).write_text(payload + name, encoding="utf-8")
    (core_root / "K").mkdir(parents=True)
    for name in ("eos.nb", "eos.t", "eos.yq", "eos.thermo"):
        (core_root / "K" / name).write_text(name, encoding="utf-8")
    monkeypatch.setattr(pipeline, "CRUST_MODELS", ["A", "B", "C"])
    monkeypatch.setattr(pipeline, "CORE_MODELS", ["K"])
    monkeypatch.setattr(pipeline, "CRUST_ROOT", crust_root)
    monkeypatch.setattr(pipeline, "CORE_ROOT", core_root)
    monkeypatch.setattr(pipeline, "OUT_ROOT", out_root)
    monkeypatch.setattr(pipeline, "executable", lambda: tmp_path / "exe")

    def fake_pipeline(method, crust_dir, core_dir, out_dir):
        if method == "css":
            return {"survived": True, "method": method, "stage": "complete",
                    "combined_n": 1, "readback_n": 1, "readback_status": "ok",
                    "readback_validation": {}, "consistency_status": "PASS",
                    "consistency_report": {}, "admission": {"status": "ADMITTED"},
                    "source_sha256": {}, "output_sha256": {},
                    "fortran_executable_sha256": None, "match_report": {}}
        return {"survived": False, "method": method, "stage": method, "error": "x",
                "failure_category": "physical"}

    monkeypatch.setattr(pipeline, "run_pipeline", fake_pipeline)
    assert pipeline.main() == 0
    import json
    manifest = json.loads((out_root / "matching_manifest.json").read_text(encoding="utf-8"))
    assert manifest["source_duplicate_groups"]["crust"] == [["A", "B"]]
    assert manifest["source_duplicate_groups"]["core"] == []
    assert manifest["total_survivors"] == 3
    assert manifest["distinct_source_survivors"] == 2


def test_census_classifies_rejected_and_placeholder_reconstruction_markers(tmp_path, monkeypatch):
    import json
    crust_root, core_root, out_root = tmp_path / "crust", tmp_path / "core", tmp_path / "run"
    (crust_root / "REJ").mkdir(parents=True)
    (crust_root / "REJ" / "REJECTED.json").write_text(json.dumps(
        {"status": "REJECTED", "reason": "mechanical stability violated at rows [4]"}))
    (crust_root / "PH").mkdir()
    (crust_root / "PH" / "NOT_GENERATED.json").write_text(json.dumps(
        {"status": "NOT_GENERATED", "reason": "placeholder identical to M1"}))
    (core_root / "K").mkdir(parents=True)
    for name in ("eos.nb", "eos.t", "eos.yq", "eos.thermo"):
        (core_root / "K" / name).write_text(name, encoding="utf-8")
    monkeypatch.setattr(pipeline, "CRUST_MODELS", ["REJ", "PH"])
    monkeypatch.setattr(pipeline, "CORE_MODELS", ["K"])
    monkeypatch.setattr(pipeline, "CRUST_ROOT", crust_root)
    monkeypatch.setattr(pipeline, "CORE_ROOT", core_root)
    monkeypatch.setattr(pipeline, "OUT_ROOT", out_root)
    monkeypatch.setattr(pipeline, "executable", lambda: tmp_path / "exe")
    monkeypatch.setattr(pipeline, "run_pipeline", lambda *a: pytest.fail("must not run"))
    assert pipeline.main() == 0
    manifest = json.loads((out_root / "matching_manifest.json").read_text(encoding="utf-8"))
    assert manifest["crust_root"] == str(crust_root)
    assert manifest["evaluated_attempts"] == 0
    by = {r["crust_model"]: r for r in manifest["all_results"] if r["method"] == "css"}
    assert by["REJ"]["failure_category"] == "source_reconstruction_rejected"
    assert "mechanical stability" in by["REJ"]["error"]
    assert by["PH"]["failure_category"] == "data_provenance"
    assert "placeholder" in by["PH"]["error"]
