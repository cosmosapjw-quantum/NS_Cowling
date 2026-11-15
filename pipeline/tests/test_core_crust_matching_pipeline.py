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
    assert result["output_sha256"]
    assert result["source_sha256"]["crust"]
    assert result["source_sha256"]["core"]
