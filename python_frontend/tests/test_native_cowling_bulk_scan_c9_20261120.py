"""C9: real-source Fortran bulk frequency scan, API and numerical parity."""
from __future__ import annotations

import ctypes as ct
import json
from pathlib import Path

import numpy as np
import pytest

from compose_eos.stellar_background import ContractEnthalpyEOS
from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_cowling_native import NativeCowlingDVODE
from compose_eos.stellar_cowling_tov_native import solve_background_native
from compose_eos import stellar_cowling_native as native

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def real_star():
    source = ROOT / "data/core_support_20261117/primary_v1/REPRESENTATION.json"
    mass = ROOT / "data/stellar_consumers_20261117/tov/CANONICAL_1P4.json"
    eos = ContractEnthalpyEOS.from_path(source)
    n_c = json.loads(mass.read_text())["central_density_fm3"]
    return solve_background_native(
        eos, n_c, surface_h=1e-4, rtol=1e-10,
        atol=1e-12, max_step_h=0.002, native_samples=4097)


@pytest.mark.parametrize("method", ["adams", "bdf"])
def test_bulk_scan_uses_fortran_and_matches_scalar_real_star(real_star, monkeypatch, method):
    settings = CowlingSettings(surface_h=1e-4, rtol=1e-9,
                               atol=1e-12, max_step_h=.002)
    solver = NativeCowlingDVODE(real_star, settings, samples=2049, method=method)
    frequencies = np.array([2200., 2350., 2475., 6200., 6550., 6750.])
    scalar = [solver.matching_residual(float(f)) for f in frequencies]

    def reject_python_shots(*args, **kwargs):
        raise AssertionError("bulk scan must not run a Python frequency loop")

    monkeypatch.setattr(solver, "matching_residual", reject_python_shots)
    brackets, scan = solver.scan_brackets(frequencies)
    np.testing.assert_allclose([p["residual"] for p in scan], scalar,
                               rtol=4e-13, atol=4e-13)
    expected = [(float(a), float(b)) for a, b, ra, rb in
                zip(frequencies[:-1], frequencies[1:], scalar[:-1], scalar[1:])
                if ra*rb < 0]
    assert brackets == expected
    assert all(p["frequency_hz"] == f for p, f in zip(scan, frequencies))


def test_native_bulk_scan_fails_closed_on_bad_frequency(real_star):
    solver = NativeCowlingDVODE(
        real_star, CowlingSettings(surface_h=1e-4), samples=1025)
    with pytest.raises(ValueError, match="strictly increasing"):
        solver.scan_brackets([2350., 2300.])
    with pytest.raises(ValueError, match="strictly increasing"):
        solver.scan_brackets([2350., float("nan")])
    with native._LOCK:
        solver._setup()
        lib = solver._lib
        bad = np.ascontiguousarray(np.array([2000., -1.], dtype=np.float64))
        residuals = np.zeros(2, dtype=np.float64)
        flags = np.zeros(1, dtype=np.int32)
        nfev, status = ct.c_int(-1), ct.c_int(-1)
        lib.cowling_native_scan(
            2, bad, 2, 1e-9, 1e-12, .0001, 10,
            solver.match_s, solver.surface_s, residuals, flags,
            ct.byref(nfev), ct.byref(status))
        assert status.value != 0
        assert nfev.value == 0
        assert np.array_equal(residuals, [0., 0.])
        assert flags[0] == 0


def test_bulk_smooth_scan_rejects_active_sharp_source_callback(real_star):
    """Public native ABI must never silently scan a sharp source phase."""
    solver = NativeCowlingDVODE(real_star, CowlingSettings(surface_h=1e-4),
                                samples=1025)
    calls = []
    prototype = ct.CFUNCTYPE(
        None, ct.c_double, ct.POINTER(ct.c_double), ct.POINTER(ct.c_int))

    @prototype
    def supplied_sharp_callback(s, values, status):
        calls.append(float(s))
        status[0] = 1

    with native._LOCK:
        solver._setup()
        lib = solver._lib
        register_status = ct.c_int(-1)
        lib.cowling_native_set_source_callback(
            ct.cast(supplied_sharp_callback, ct.c_void_p),
            ct.byref(register_status))
        assert register_status.value == 0
        try:
            frequencies = np.array([2200., 2400.], dtype=np.float64)
            residuals = np.ones(2, dtype=np.float64)
            flags = np.ones(1, dtype=np.int32)
            nfev, status = ct.c_int(-1), ct.c_int(-1)
            lib.cowling_native_scan(
                2, frequencies, 2, 1e-9, 1e-12, .0001, 10,
                solver.match_s, solver.surface_s, residuals, flags,
                ct.byref(nfev), ct.byref(status))
            assert status.value == 2
            assert calls == []
            assert nfev.value == 0
            assert np.array_equal(residuals, [0., 0.])
            assert flags[0] == 0
        finally:
            lib.cowling_native_set_source_callback(
                None, ct.byref(register_status))
            assert register_status.value == 0
