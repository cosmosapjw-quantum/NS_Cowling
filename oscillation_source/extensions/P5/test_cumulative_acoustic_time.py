from __future__ import annotations

import csv
import math
from pathlib import Path
import subprocess
import sys

C_SI = 299_792_458.0
HERE = Path(__file__).resolve().parent


def test_local_analytic_profile_has_tau_equal_two_r_over_c(tmp_path):
    out = tmp_path / "tau.csv"
    meta = tmp_path / "tau.json"
    proc = subprocess.run(
        [
            sys.executable,
            str(HERE / "cumulative_acoustic_time.py"),
            "--background",
            str(HERE / "analytic_background.csv"),
            "--output",
            str(out),
            "--metadata",
            str(meta),
        ],
        cwd=HERE.parent.parent,
        text=True,
        capture_output=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    with out.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 3
    for row in rows:
        r = float(row["r_m"])
        tau = float(row["cumulative_acoustic_time_s"])
        assert math.isclose(tau, 2.0 * r / C_SI, rel_tol=2e-15, abs_tol=0.0)


def _run_background(tmp_path, points):
    src = tmp_path/'background.csv'
    with src.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['r_m', 'chi', 'nu', 'lambda'])
        w.writeheader()
        w.writerows(points)
    out = tmp_path/'output.csv'
    proc = subprocess.run([sys.executable, str(HERE/'cumulative_acoustic_time.py'),
                           '--background', str(src), '--output', str(out),
                           '--metadata', str(tmp_path/'meta.json')],
                          capture_output=True, text=True, check=False)
    return proc, out


def test_nonuniform_grid_with_nonconstant_integrand(tmp_path):
    # r=1,3,7 m; g=(2,4,8)/c s/m. The center is first-value hold,
    # not a fictitious zero-valued central node. Trapezoid increments: 6/c,24/c.
    points = [dict(r_m=r, chi=0.25, nu=0.0, **{'lambda': math.log(a)})
              for r,a in [(1,1),(3,2),(7,4)]]
    p,out = _run_background(tmp_path,points)
    assert p.returncode == 0, p.stdout+p.stderr
    rows=list(csv.DictReader(out.open()))
    assert len(rows)==3
    for row,expected in zip(rows,[2.0,8.0,32.0]):
        assert math.isclose(float(row['cumulative_acoustic_time_s']),expected/C_SI,rel_tol=3e-15)


def test_invalid_profile_does_not_emit_success_output(tmp_path):
    points = [dict(r_m=r, chi=0.25, nu=0.0, **{'lambda':0.0}) for r in [1,1,3]]
    p,out = _run_background(tmp_path,points)
    assert p.returncode == 2, p.stdout+p.stderr
    assert not out.exists()


def test_lapse_shift_rescales_postprocessed_time(tmp_path):
    shift=math.log(3.0)
    points=[dict(r_m=r, chi=0.25, nu=shift, **{'lambda':0.0}) for r in [1,2,3]]
    p,out=_run_background(tmp_path,points)
    assert p.returncode==0,p.stdout+p.stderr
    for row in csv.DictReader(out.open()):
        assert math.isclose(float(row['cumulative_acoustic_time_s']),
                            2*float(row['r_m'])/(3*C_SI),rel_tol=3e-15)
