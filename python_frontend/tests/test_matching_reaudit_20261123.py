"""Fortran regression controls for source-grid search and interval causality.

Python builds and runs the native test programs; it performs no EOS calculation.
"""
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("name", ["matching_search_exhaustive_20261123", "firstlaw_interval_20261123"])
@pytest.mark.parametrize("optimized", [False, True], ids=["checked", "optimized"])
def test_native_matching_reaudit(name, optimized, tmp_path):
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran unavailable")
    source = ROOT / "fortran_eos/src/compose_eos_matching.f90"
    text = source.read_text()
    marker = "public :: match_crust_core"
    exported = marker + ", pressure_work_integral, make_pressure_spline, spline_eval, construct_bridge, validate_bridge"
    probe = tmp_path / "matching_probe.f90"
    probe.write_text(text.replace(marker, exported, 1))
    assert probe.read_text().replace(exported, marker, 1) == text
    flags = ["-O2"] if optimized else ["-O1", "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow"]
    binary = tmp_path / name
    command = [compiler, *flags, "-ffree-line-length-none", "-J", str(tmp_path), "-I", str(tmp_path),
               str(ROOT / "fortran_eos/src/compose_eos_types.f90"), str(probe),
               str(ROOT / f"fortran_eos/tests/test_{name}.f90"), "-o", str(binary)]
    compiled = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert compiled.returncode == 0, compiled.stdout + compiled.stderr
    result = subprocess.run([str(binary)], cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert ("FAILURES=0" if name.startswith("matching_search") else
            "first-law interval counterexample and valid control: PASS") in result.stdout
