from __future__ import annotations

from pathlib import Path

from scripts import match_core_crust as pipeline


def _write_internal_csv(path: Path, *, q1_second: float = 2.0) -> None:
    path.write_text(
        "nb,Q1,Q2,Q3,Q4,Q5,Q6,Q7,q_extra,Iphase,Yn,Yp,Ye,Ymu,Abar,Zbar,Yh,"
        "row_status,yn_present,yp_present,ye_present,ymu_present,abar_present,"
        "zbar_present,yh_present,c3_applicable,c4_applicable\n"
        "0.1,1,0,0.1,0,0,0,0,0,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n"
        f"0.2,{q1_second},0,0.2,0,0,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n",
        encoding="utf-8",
    )


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
