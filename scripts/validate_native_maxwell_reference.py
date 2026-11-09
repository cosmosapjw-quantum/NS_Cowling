#!/usr/bin/env python3
"""Independent SciPy cross-check and root-convergence study for the native fixture.

The EOS is reconstructed analytically from the same knots written by
maxwell_table_checks.f90. It is a generated integration fixture, not a
source-backed CompOSE EOS or an astrophysical model.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import EOSData
from compose_cowling import MaxwellEOS, scan_maxwell_cowling_modes, solve_maxwell_tov
from compose_cowling.maxwell import MEV_FM3_TO_J_M3

GL_X = np.array([
    0.095012509837637440185, 0.281603415044923588867,
    0.458016777657227386342, 0.617876244402643748447,
    0.755404408355003033895, 0.865631202387831743880,
    0.944575023073232576078, 0.989400934991649932596,
])
GL_W = np.array([
    0.189450610455068496285, 0.182603415044923588867,
    0.169156519395002538189, 0.149595988816576732081,
    0.124628971255533872052, 0.095158511682492784810,
    0.062253523938647892863, 0.027152459411754094852,
])


def log_density_integral(ua: float, ub: float, g_pt: float) -> float:
    mid, half = 0.5 * (ua + ub), 0.5 * (ub - ua)
    up, um = mid + half * GL_X, mid - half * GL_X
    slope = 2.0
    fp = (np.exp(up) + slope) / (2.0 * np.exp(up) + g_pt + slope * up)
    fm = (np.exp(um) + slope) / (2.0 * np.exp(um) + g_pt + slope * um)
    return float(half * np.dot(GL_W, fp + fm))


def phase_rows(pmin: float, pmax: float, g_pt: float, n_pt: float, phase: int) -> dict[str, np.ndarray]:
    p = np.geomspace(pmin, pmax, 24)
    u = np.log(p)
    logn = np.empty_like(p)
    logn[0] = math.log(n_pt)
    if u[0] < 0.0:
        logn[0] -= log_density_integral(float(u[0]), 0.0, g_pt)
    for i in range(1, p.size):
        logn[i] = logn[i - 1] + log_density_integral(float(u[i - 1]), float(u[i]), g_pt)
    n = np.exp(logn)
    g = g_pt + 2.0 * u
    eps = p + g
    mu = (eps + p) / n
    return {"p": p, "eps": eps, "n": n, "mu": mu, "phase": np.full(p.size, phase, dtype=int)}


def make_eos_data() -> EOSData:
    p_t, g_low, g_high, mu_t = 1.0, 939.0, 1139.0, 941.0
    low = phase_rows(0.1, p_t, g_low, 1.0, 1)
    high = phase_rows(p_t, 2.0, g_high, (2.0 * p_t + g_high) / mu_t, 2)
    fields = {k: np.concatenate((low[k], high[k])) for k in ("p", "eps", "n", "mu", "phase")}
    n, p, eps, mu = fields["n"], fields["p"], fields["eps"], fields["mu"]
    zeros = np.zeros_like(n)
    q7 = eps / (n * 939.0) - 1.0
    return EOSData(
        nb=n, mn=939.0, mp=939.0, Q1=p / n, Q2=zeros.copy(),
        Q3=mu / 939.0 - 1.0, Q4=zeros.copy(), Q5=zeros.copy(),
        Q6=q7.copy(), Q7=q7.copy(), q_extra=zeros.copy(), Iphase=fields["phase"],
        Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
        Aav=zeros.copy(), Zav=zeros.copy(), Yh=zeros.copy(),
        Il=1, iT=1, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_file_present=False, composition_available=False,
        beta_equilibrated=True,
    )


def parse_tsv(path: Path) -> dict[str, object]:
    out: dict[str, object] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        row = line.split("\t")
        if len(row) != 3:
            continue
        key, kind, raw = row
        if kind == "bool":
            out[key] = raw.strip().lower() == "true"
        elif kind == "int":
            out[key] = int(raw)
        elif kind == "real":
            out[key] = float(raw)
        else:
            out[key] = raw
    return out


def roots_from_log(data: dict[str, object]) -> list[float]:
    if data.get("spectrum_scan_valid") is not True or data.get("spectrum_all_trials_succeeded") is not True:
        raise RuntimeError("native Maxwell frequency scan is incomplete")
    count = int(data.get("spectrum_roots_found", -1))
    if count != 2:
        raise RuntimeError(f"expected two roots in this fixture scan, got {count}")
    roots = []
    for i in (1, 2):
        if data.get(f"spectrum_root_valid_{i:05d}") is not True:
            raise RuntimeError(f"native root {i} was not valid")
        roots.append(float(data[f"spectrum_root_hz_{i:05d}"]))
    return sorted(roots)


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), np.finfo(float).tiny)


def main() -> int:
    ap = argparse.ArgumentParser()
    for arg in (
        "slow", "rapid", "slow_ode_tight", "rapid_ode_tight",
        "slow_root_grid_tight", "rapid_root_grid_tight",
    ):
        ap.add_argument(f"--{arg.replace('_', '-')}-log", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    logs = {
        "slow": parse_tsv(args.slow_log),
        "rapid": parse_tsv(args.rapid_log),
        "slow_ode_tight": parse_tsv(args.slow_ode_tight_log),
        "rapid_ode_tight": parse_tsv(args.rapid_ode_tight_log),
        "slow_root_grid_tight": parse_tsv(args.slow_root_grid_tight_log),
        "rapid_root_grid_tight": parse_tsv(args.rapid_root_grid_tight_log),
    }
    convergence: dict[str, object] = {}
    max_ode_shift = 0.0
    max_grid_root_shift = 0.0
    for conv in ("slow", "rapid"):
        baseline = roots_from_log(logs[conv])
        conv_data: dict[str, object] = {"baseline_roots_hz": baseline}
        for variant in ("ode_tight", "root_grid_tight"):
            refined = roots_from_log(logs[f"{conv}_{variant}"])
            deltas = [rel(a, b) for a, b in zip(baseline, refined)]
            delta = max(deltas)
            if variant == "ode_tight":
                max_ode_shift = max(max_ode_shift, delta)
            else:
                max_grid_root_shift = max(max_grid_root_shift, delta)
            conv_data[variant] = {
                "roots_hz": refined,
                "relative_shifts": deltas,
                "max_relative_shift": delta,
            }
        convergence[conv] = conv_data
    if max_ode_shift > 2.0e-3:
        raise RuntimeError(f"ODE refinement root shift >2e-3: {max_ode_shift:.6e}")
    if max_grid_root_shift > 2.0e-5:
        raise RuntimeError(f"grid/root refinement shift >2e-5: {max_grid_root_shift:.6e}")

    eos = MaxwellEOS(make_eos_data())
    pc, pt, pcut = (1.4 * MEV_FM3_TO_J_M3, 1.0 * MEV_FM3_TO_J_M3, 0.2 * MEV_FM3_TO_J_M3)
    star = solve_maxwell_tov(
        eos, pc, pressure_cutoff_pa=pcut, r0_m=1.0, r_limit_m=2.0e6,
        rtol=2.0e-9, atol_mass_kg=1.0e13, atol_nu=1.0e-12,
        max_step_log_pressure=0.05,
    )
    frequency_grid = np.linspace(100.0, 20000.0, 41)
    python_roots: dict[str, list[float]] = {}
    root_comparison: dict[str, list[dict[str, float]]] = {}
    max_reference_root_delta = 0.0
    for conv in ("slow", "rapid"):
        refs = scan_maxwell_cowling_modes(
            star, frequency_grid, conversion=conv, rtol=2.0e-9, max_step_m=5.0
        )
        python_roots[conv] = sorted(float(row["frequency_hz"]) for row in refs)
        natives = roots_from_log(logs[conv])
        if len(python_roots[conv]) != 2:
            raise RuntimeError(f"Python/SciPy {conv} scan did not return two roots")
        pairs = [
            {"native_hz": n, "python_scipy_hz": p, "relative_difference": rel(n, p)}
            for n, p in zip(natives, python_roots[conv])
        ]
        root_comparison[conv] = pairs
        max_reference_root_delta = max(
            max_reference_root_delta,
            *(row["relative_difference"] for row in pairs),
        )
    native = logs["slow"]
    dm = rel(float(native["mass_kg"]), star.mass_kg)
    dr = rel(float(native["radius_m"]), star.radius_m)
    drt = rel(float(native["interface_radius_m"]), star.interface.radius_m)
    if max_reference_root_delta > 5.0e-3:
        raise RuntimeError(f"native/Python root discrepancy >5e-3: {max_reference_root_delta:.6e}")
    if max(dm, dr, drt) > 2.0e-4:
        raise RuntimeError(
            "native/Python TOV discrepancy >2e-4: "
            f"M={dm:.6e}, R={dr:.6e}, rt={drt:.6e}"
        )

    result = {
        "schema": "native-maxwell-convergence-and-scipy-crosscheck-v1",
        "status": "PASS_BOUNDED_SYNTHETIC_CONVERGENCE_AND_REFERENCE",
        "fixture": "generated analytic two-table EOS with g=epsilon-P linear in ln(P)",
        "frequency_grid_hz": [100.0, 20000.0],
        "frequency_grid_points": 41,
        "root_family_identity": "unresolved; root sets sorted only for numerical comparison",
        "root_convergence": {
            "max_relative_shift_ode_tight": max_ode_shift,
            "max_relative_shift_grid_root_tight": max_grid_root_shift,
            "slow": convergence["slow"],
            "rapid": convergence["rapid"],
            "acceptance_limits": {"ode_tight": 2.0e-3, "grid_root_tight": 2.0e-5},
        },
        "python_scipy_reference": {
            "tov": {
                "native_mass_kg": float(native["mass_kg"]), "python_mass_kg": star.mass_kg,
                "mass_relative_difference": dm,
                "native_radius_m": float(native["radius_m"]), "python_radius_m": star.radius_m,
                "radius_relative_difference": dr,
                "native_interface_radius_m": float(native["interface_radius_m"]),
                "python_interface_radius_m": star.interface.radius_m,
                "interface_radius_relative_difference": drt,
            },
            "cowling_root_comparison": root_comparison,
            "max_root_relative_difference": max_reference_root_delta,
            "settings": {
                "tov_rtol": 2.0e-9, "tov_max_step_log_pressure": 0.05,
                "cowling_rtol": 2.0e-9, "cowling_max_step_m": 5.0,
                "root_method": "SciPy brentq",
            },
            "acceptance_limits": {"tov_relative_difference": 2.0e-4, "root_relative_difference": 5.0e-3},
        },
        "interpretation": (
            "Native and independent Python/SciPy implementations are compared on the same generated analytic "
            "two-phase fixture. This is not a source-backed CompOSE model, an astrophysical prediction, "
            "or a mode-family identification."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
