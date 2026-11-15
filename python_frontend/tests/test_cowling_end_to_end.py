"""End-to-end Cowling spectrum tests (E390-06 / M6 remaining subscope:
compose_cowling runtime-compile path resolution).

These tests pin the M6 Cowling end-to-end physics path that was previously
blocked: the runtime-compile of the 67-file oscillation_source tree now
produces a working cell_solve, and a causal beta-equilibrated EOS yields a
physically reasonable f-mode (~2-3 kHz) and p1-mode (~7-8 kHz) spectrum.

Background on the two failure modes these tests guard against:
  1. Causality: the Fortran ns_table_eos requires g = eps - P to be strictly
     increasing in u = log P (c_s^2 = dP/deps <= 1). APR_APR has 23
     superluminal intervals at its high-density tail and is correctly rejected.
  2. Mass closure: solve_cowling_spectrum now defaults to target_mass_kg=0,
     computing the spectrum directly at pc_pa rather than bisecting the central
     pressure toward a fixed 1.4 M_sun target (which fails for EOS that do not
     reach that mass within the pc bracket).
"""
import os
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_ROOT", str(ROOT))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_EXE", str(ROOT / "build" / "compose_eos_fortran" / "compose_eos_fortran"))

from compose_eos import read_compose  # noqa: E402
from compose_cowling import solve_cowling_spectrum  # noqa: E402
from compose_cowling.solver import NativeExecutionError  # noqa: E402

RAW = ROOT / "inputs" / "raw_compose"
MEV_FM3_TO_PA = 1.602176634e32


def _read_or_skip(name):
    path = RAW / name
    if not path.is_dir():
        pytest.skip(f"inputs/raw_compose/{name} not staged (M0 not run)")
    return read_compose(str(path))


def test_cowling_spectrum_on_causal_eos_yields_f_and_p1_modes():
    """BBB_BHF2 is a causal (c_s^2 <= 1) beta-equilibrated EOS spanning the
    full crust-to-core density range. With the default frequency brackets
    (f in [2000,3000] Hz, p1 in [6500,8000] Hz) and target_mass_kg=0, the
    spectrum must return both modes in the physically expected bands."""
    eos = _read_or_skip("BBB_BHF2")
    pc_pa = 200.0 * MEV_FM3_TO_PA  # ~200 MeV/fm^3 central pressure
    result = solve_cowling_spectrum(eos, pc_pa)
    assert result["status"] == "ok", f"spectrum status: {result.get('status')}"
    assert result["backend"] == "Fortran"
    f_hz = result["f_hz"]
    p1_hz = result["p1_hz"]
    # f-mode: ~2-3 kHz for a typical neutron-star Cowling spectrum
    assert 2000.0 <= f_hz <= 3500.0, f"f-mode {f_hz} Hz outside [2000,3500]"
    # p1-mode: higher than f-mode, ~7-8 kHz
    assert 6500.0 <= p1_hz <= 9000.0, f"p1-mode {p1_hz} Hz outside [6500,9000]"
    assert p1_hz > f_hz, "p1-mode must be higher than f-mode"


def test_cowling_spectrum_frequency_brackets_are_independent():
    """A higher central pressure (denser, more compact star) should shift the
    Cowling f-mode upward. This pins the physical monotonicity of the mode
    frequency in central pressure (more compact => higher f-mode)."""
    eos = _read_or_skip("BBB_BHF2")
    f_low = solve_cowling_spectrum(eos, 100.0 * MEV_FM3_TO_PA)["f_hz"]
    f_high = solve_cowling_spectrum(eos, 500.0 * MEV_FM3_TO_PA)["f_hz"]
    assert f_high > f_low, (
        f"f-mode not monotone in pc: f(100)={f_low} >= f(500)={f_high}")


def test_cowling_rejects_superluminal_eos_with_causality_error():
    """APR_APR has 23 superluminal intervals at its high-density tail (max
    c_s^2 ~ 1.48). The Fortran ns_table_eos must reject it with a
    matching_inconsistent causality error, not silently produce a spectrum."""
    eos = _read_or_skip("APR_APR")
    pc_pa = 200.0 * MEV_FM3_TO_PA
    with pytest.raises(NativeExecutionError, match="not admissible for monotone causal g"):
        solve_cowling_spectrum(eos, pc_pa)
