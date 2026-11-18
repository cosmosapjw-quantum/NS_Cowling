"""Adversarial Maxwell sharp-branch Fortran DVODE acceptance.

No Python ODE solver: test fixture is a physically causal analytic two-
branch P(epsilon) table. The original native Takatsy-Kovacs jump contract
is checked under Adams, BDF and corrected RK45.
"""
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture(scope="module")
def maxwell_files(tmp_path_factory):
    folder = tmp_path_factory.mktemp("maxwell_native_2phase")
    high = folder / "high.csv"
    low = folder / "low.csv"
    np.savetxt(high, np.array([[1.0, 153.0], [1.2, 153.6],
                               [1.5, 154.5], [2.0, 156.0]]),
               delimiter=",", fmt="%.17e")
    np.savetxt(low, np.array([[0.01, 100.03], [0.1, 100.3],
                              [0.5, 101.5], [1.0, 103.0]]),
               delimiter=",", fmt="%.17e")
    return high, low


def invoke(high, low, mode, delta=50.0, pt=1.0, pc=1.8):
    env = os.environ.copy()
    env["NS_COWLING_TOV_INTEGRATOR"] = mode
    env["NS_COWLING_TOV_RTOL"] = "2e-9"
    binary = Path(env["COMPOSE_EOS_FORTRAN_EXE"])
    return subprocess.run([str(binary), "tov-love", "maxwell",
                           str(high), str(low), str(pt), str(delta), str(pc)],
                          env=env, capture_output=True, text=True, timeout=45)


def test_three_fortran_methods_preserve_maxwell_jump(maxwell_files):
    high, low = maxwell_files
    results = {}
    for method in ("dvode_adams", "dvode_bdf", "rk45"):
        process = invoke(high, low, method)
        assert process.returncode == 0, (method, process.stdout, process.stderr)
        row = json.loads(process.stdout)
        assert row["status"] == "ok"
        assert row["converged"]
        assert row["has_transition"]
        assert 0.0 < row["r_transition_km"] < row["R_km"]
        assert 0.0 < row["m_transition_msun"] < row["M_msun"]
        assert row["delta_y"] < 0.0
        assert row["delta_eps_MeV_fm3"] == 50.0
        assert 0 < row["k2"] < 1
        results[method] = row

    ref = results["dvode_adams"]
    for method in ("dvode_bdf", "rk45"):
        for field, rel in (("M_msun", 8e-5), ("R_km", 8e-5),
                           ("k2", 5e-4), ("Lambda", 5e-4),
                           ("r_transition_km", 2e-4), ("delta_y", 2e-4)):
            p = results[method][field]
            assert abs(p-ref[field])/abs(ref[field]) < rel, (method, field, ref[field], p)


def test_unmatched_maxwell_latent_heat_fails_closed(maxwell_files):
    high, low = maxwell_files
    for method in ("dvode_adams", "dvode_bdf"):
        for delta in (49.0, -1.0, 0.0):
            p = invoke(high, low, method, delta=delta)
            assert p.returncode != 0, (method, delta, p.stdout)


def test_missing_bracket_and_out_of_range_fails_closed(maxwell_files, tmp_path):
    high, low = maxwell_files
    altered = tmp_path / "invalid_high.csv"
    altered.write_text("1.1,153.0\n1.3,153.6\n1.6,154.5\n2.0,156.0\n")
    for method in ("dvode_adams", "dvode_bdf"):
        assert invoke(altered, low, method).returncode != 0
        assert invoke(high, low, method, pt=.01).returncode != 0
        assert invoke(high, low, method, pc=2.01).returncode != 0
