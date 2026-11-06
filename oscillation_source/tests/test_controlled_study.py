"""T15 structure/acoustic diagnostics on analytic and admitted real backgrounds."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import pytest

from tools.runtime import CaseRunner

ROOT = Path(__file__).resolve().parents[1]
ANALYTIC = ROOT / "fixtures/analytic/t15_constant_background.csv"
REAL = ROOT / "fixtures/tables/compose69_l2_modes/background.csv"
REAL_CONTRACT = ROOT / "fixtures/tables/compose69_l2_modes/contract.json"
C = 299792458.0


def _runner(tmp_path):
    return CaseRunner(ROOT, tmp_path / "runs", workspace=os.environ.get("NS_WORKSPACE", "solutions"))


def test_constant_background_geometry_acoustic_and_dimensionless_frequency(tmp_path):
    r = _runner(tmp_path).run_case(
        "structure_diagnostics",
        {
            "background_path": str(ANALYTIC),
            "split_pressure_pa": 2.0e20,
            "f_hz": 100.0,
            "p1_hz": 200.0,
        },
    )
    assert r["status"] == "ok", r
    assert r["result_valid"] is True
    assert r["radius_m"] == pytest.approx(3.0, rel=0, abs=1e-14)
    assert r["proper_radius_m"] == pytest.approx(3.0, rel=3e-14)
    assert r["acoustic_time_s"] == pytest.approx(6.0 / C, rel=3e-14)
    assert r["r_split_m"] == pytest.approx(2.0, rel=0, abs=1e-14)
    assert r["outer_areal_thickness_m"] == pytest.approx(1.0, rel=3e-14)
    assert r["outer_proper_thickness_m"] == pytest.approx(1.0, rel=3e-14)
    assert r["outer_acoustic_time_s"] == pytest.approx(2.0 / C, rel=3e-14)
    dyn = math.sqrt(6.67430e-11 * 3.0 / 3.0**3)
    assert r["f_hat"] == pytest.approx(2 * math.pi * 100.0 / dyn, rel=5e-14)
    assert r["p1_hat"] == pytest.approx(2 * math.pi * 200.0 / dyn, rel=5e-14)


def test_compose69_structure_diagnostics_bind_to_t14_background(tmp_path):
    c = json.loads(REAL_CONTRACT.read_text())
    # Use an actual stored pressure so the split location itself is not an interpolation oracle.
    rows = REAL.read_text().splitlines()
    split_pressure = float(rows[len(rows)//2].split(",")[3])
    r = _runner(tmp_path).run_case(
        "structure_diagnostics",
        {
            "background_path": str(REAL),
            "split_pressure_pa": split_pressure,
            "f_hz": c["mode_calibration"]["f_hz_observed"],
            "p1_hz": c["mode_calibration"]["p1_hz_observed"],
        },
    )
    assert r["status"] == "ok", r
    assert r["result_valid"] is True
    assert r["radius_m"] == pytest.approx(c["star_observed"]["radius_m"], rel=2e-14)
    assert r["mass_kg"] > 0
    assert r["proper_radius_m"] > r["radius_m"]
    assert r["acoustic_time_s"] > 0
    assert 0 < r["sound_speed_min_m_s"] <= r["sound_speed_max_m_s"] <= C
    assert 0 < r["outer_areal_thickness_m"] < r["radius_m"]
    assert r["outer_proper_thickness_m"] > r["outer_areal_thickness_m"]
    assert r["outer_acoustic_time_s"] > 0
    assert r["f_hat"] > 0 and r["p1_hat"] > r["f_hat"]
