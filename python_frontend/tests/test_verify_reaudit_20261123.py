"""Native verifier input-domain and phase-interface regression (2026-11-23)."""
from __future__ import annotations

from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize(
    "flags",
    [
        ("-O0", "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow"),
        ("-O2", "-fcheck=all"),
    ],
    ids=("debug-fpe", "optimized"),
)
def test_native_verifier_domains_and_phase_boundaries(tmp_path, flags):
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran is required for the native verifier regression")
    sources = [
        ROOT / "fortran_eos/src/compose_eos_types.f90",
        ROOT / "fortran_eos/src/compose_eos_verify.f90",
        ROOT / "fortran_eos/tests/test_verify_domains_20261123.f90",
    ]
    executable = tmp_path / "verify_domains"
    compiled = subprocess.run(
        [compiler, "-std=f2008", "-ffree-line-length-none", *flags,
         *(str(path) for path in sources), "-o", str(executable)],
        cwd=tmp_path, capture_output=True, text=True, timeout=90,
    )
    assert compiled.returncode == 0, compiled.stdout + compiled.stderr
    completed = subprocess.run(
        [str(executable)], cwd=tmp_path, capture_output=True, text=True, timeout=30,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "FAILURES=0" in completed.stdout, completed.stdout
