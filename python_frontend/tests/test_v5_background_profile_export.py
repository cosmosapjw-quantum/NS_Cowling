from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from pathlib import Path

import pytest

from compose_cowling import solver


HEADER = "cell,t,r_m,m_kg,p_pa,energy_j_m3,chi,nu,w_over_r3,v_over_r2\n"
ROWS = (
    "1,1.0,1.0e1,1.0e15,9.0e30,1.0e34,0.1,-0.3,0.0,0.0\n"
    "1,0.0,1.0e4,1.0e30,1.0e20,1.0e32,0.2,-0.1,0.0,0.0\n"
)


def test_v5_tov_exports_exact_native_background_profile(tmp_path, monkeypatch):
    eos = tmp_path / "eos.csv"
    eos.write_text("pressure_pa,energy_j_m3\n1.0e20,1.0e32\n", encoding="utf-8")
    destination = tmp_path / "profiles" / "profile.csv"

    monkeypatch.setattr(
        solver,
        "_compile_v5_cell_solve",
        lambda: (Path("/fake/cell_solve"), "test-build-key"),
    )

    def fake_run(command, capture_output, text):
        prefix = Path(command[-1])
        prefix.with_suffix(".json").write_text(
            json.dumps({"status": "ok", "mass_kg": 1.0e30, "radius_m": 1.0e4}),
            encoding="utf-8",
        )
        Path(str(prefix) + "_background.csv").write_text(HEADER + ROWS, encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(solver.subprocess, "run", fake_run)
    result = solver.solve_v5_tov_star(
        eos,
        1.0e34,
        p_cut_pa=1.0e20,
        background_output_path=destination,
    )

    assert destination.is_file()
    raw = destination.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == result["background_profile"]["sha256"]
    assert result["background_profile"]["row_count"] == 2
    assert result["backend_build_key"] == "test-build-key"
    with destination.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        assert tuple(reader.fieldnames) == (
            "cell", "t", "r_m", "m_kg", "p_pa", "energy_j_m3",
            "chi", "nu", "w_over_r3", "v_over_r2",
        )
        assert len(list(reader)) == 2


def test_v5_background_export_refuses_existing_destination(tmp_path, monkeypatch):
    eos = tmp_path / "eos.csv"
    eos.write_text("pressure_pa,energy_j_m3\n1.0e20,1.0e32\n", encoding="utf-8")
    destination = tmp_path / "profile.csv"
    destination.write_text("do not overwrite\n", encoding="utf-8")
    monkeypatch.setattr(
        solver,
        "_compile_v5_cell_solve",
        lambda: (Path("/fake/cell_solve"), "test-build-key"),
    )

    def fake_run(command, capture_output, text):
        prefix = Path(command[-1])
        prefix.with_suffix(".json").write_text(
            json.dumps({"status": "ok", "mass_kg": 1.0e30, "radius_m": 1.0e4}),
            encoding="utf-8",
        )
        Path(str(prefix) + "_background.csv").write_text(HEADER + ROWS, encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(solver.subprocess, "run", fake_run)
    with pytest.raises(FileExistsError):
        solver.solve_v5_tov_star(
            eos,
            1.0e34,
            p_cut_pa=1.0e20,
            background_output_path=destination,
        )
    assert destination.read_text(encoding="utf-8") == "do not overwrite\n"
