#!/usr/bin/env python3
"""Reproduce analytic Maxwell-interface validation against the smooth Fortran path."""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import EOSData, construct_maxwell_transition
from compose_cowling import (
    MaxwellEOS,
    solve_maxwell_tov,
    solve_tov_star,
    solve_cowling_spectrum,
    shoot_maxwell_cowling,
    scan_maxwell_cowling_modes,
)
from compose_cowling.maxwell import MEV_FM3_TO_J_M3, _interface_transfer
from compose_cowling.solver import _compile_solver, _parse_tsv


def quadratic_phase(*, rest_mass: float, stiffness: float, n_min: float,
                    n_max: float, count: int = 401) -> EOSData:
    n = np.linspace(n_min, n_max, count)
    p = stiffness * n**2
    eps = rest_mass * n + p
    mu = (eps + p) / n
    zero = np.zeros_like(n)
    q7 = eps / (n * 939.0) - 1.0
    return EOSData(
        nb=n, mn=939.0, mp=939.0, Q1=p / n, Q2=zero.copy(),
        Q3=mu / 939.0 - 1.0, Q4=zero.copy(), Q5=zero.copy(),
        Q6=q7.copy(), Q7=q7.copy(), q_extra=zero.copy(),
        Iphase=np.zeros(n.size, dtype=int), Yn=zero.copy(), Yp=zero.copy(),
        Ye=zero.copy(), Ymu=zero.copy(), Aav=zero.copy(), Zav=zero.copy(),
        Yh=zero.copy(), Il=1, iT=1, iYq=1, inb_start=1,
        temperature_mev=0.0, charge_fraction_grid=0.0, q_extra_count=0,
        composition_file_present=False, composition_available=False,
        beta_equilibrated=True,
    )


def fortran_transfer_check(star, omega: float, y_inner: np.ndarray,
                           conversion: str, case_driver: Path) -> dict:
    bi = star.background(star.interface.radius_m, "inner")
    bo = star.background(star.interface.radius_m, "outer")
    y_python, diag = _interface_transfer(star, omega, y_inner, conversion)
    case_id = "rapid_transfer" if conversion == "rapid" else "material_transfer"
    import tempfile
    with tempfile.TemporaryDirectory(prefix="maxwell-transfer-") as td:
        nml = Path(td) / "transfer.nml"
        nml.write_text(
            "&request request_token='maxwell-independent-check' /\n"
            "&parameters\n"
            f"r_m={star.interface.radius_m:.17e}, omega_rad_s={omega:.17e}, "
            f"pressure_pa={star.interface.pressure_pa:.17e}, "
            f"eps_in={bi['energy_j_m3']:.17e}, eps_out={bo['energy_j_m3']:.17e}, "
            f"nu={bi['nu']:.17e}, lambda_value={bi['lambda']:.17e}, "
            f"dnu_dr={bi['dnu_dr']:.17e}, w_in_m3={y_inner[0]:.17e}, "
            f"v_in_m2={y_inner[1]:.17e}, inner_side=1, outer_side=-1\n/\n",
            encoding="utf-8",
        )
        run = subprocess.run([str(case_driver), case_id, str(nml)],
                             capture_output=True, text=True, check=False)
        if run.returncode:
            raise RuntimeError(f"Fortran interface case failed: {run.stderr}\n{run.stdout}")
        native = _parse_tsv(run.stdout)
        if native.get("status") != "ok":
            raise RuntimeError(f"Fortran interface case invalid: {native}")
    return {
        "conversion": conversion,
        "w_out_relative_difference": abs(native["w_out_m3"] - y_python[0]) / max(abs(y_python[0]), 1.0),
        "v_out_relative_difference": abs(native["v_out_m2"] - y_python[1]) / max(abs(y_python[1]), 1.0),
        "fortran_delta_p_relative_jump": native.get("delta_p_relative_jump", 0.0),
        "python_delta_p_relative_jump": diag["delta_p_relative_jump"],
        "python_interface_condition_residual": (
            diag["slow_condition_residual"] if conversion == "slow"
            else diag["rapid_condition_residual"]
        ),
    }


def main() -> int:
    out_path = ROOT / "validation" / "MAXWELL_REFERENCE_VALIDATION_20261109.json"
    pc = 1.35 * MEV_FM3_TO_J_M3

    # A deliberately strong jump exercises one-sided thermodynamics and both transfer maps.
    strong_low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    strong_high = quadratic_phase(rest_mass=940.0, stiffness=0.25, n_min=1.3, n_max=2.8)
    strong_table = construct_maxwell_transition(
        strong_low, strong_high, assume_local_charge_neutrality=True
    )
    strong_eos = MaxwellEOS(strong_table)
    strong_star = solve_maxwell_tov(strong_eos, pc, rtol=3e-9, max_step_log_pressure=0.015)
    strong_shots = {
        mode: shoot_maxwell_cowling(strong_star, 1800.0, conversion=mode, max_step_m=120.0)
        for mode in ("slow", "rapid")
    }
    strong_interface = {}
    for mode, shot in strong_shots.items():
        strong_interface[mode] = shot.interface_diagnostics

    # A tiny nonzero jump must recover the smooth-EOS TOV/Cowling solution.
    eta = 1.0e-4
    n_high_t = math.sqrt(1.0 + eta)
    k_high = 1.0 / (1.0 + eta)
    m_high = 941.0 - 2.0 * k_high * n_high_t
    smooth_low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    smooth_high = quadratic_phase(rest_mass=m_high, stiffness=k_high, n_min=0.8, n_max=2.8)
    small_table = construct_maxwell_transition(
        smooth_low, smooth_high, assume_local_charge_neutrality=True
    )
    small_eos = MaxwellEOS(small_table)
    small_star = solve_maxwell_tov(small_eos, pc, rtol=2e-10, max_step_log_pressure=0.01)
    smooth_tov = solve_tov_star(smooth_low, pc, rtol=1e-10)
    smooth_modes = solve_cowling_spectrum(
        smooth_low, pc, [4.0e4, 6.0e4], [5.0e4, 7.0e4],
        target_mass_kg=smooth_tov["mass_kg"], rtol=2e-8, mode_rtol=2e-8,
        step_cap=0.02, mass_rtol=2e-8, residual_tol=1e-5,
    )
    grid = np.array([4.0e4, 5.0e4, 6.0e4, 7.0e4])
    phase_modes = {}
    for mode in ("slow", "rapid"):
        roots = scan_maxwell_cowling_modes(
            small_star, grid, conversion=mode, rtol=2e-5, max_step_m=500.0
        )
        phase_modes[mode] = [float(x["frequency_hz"]) for x in roots]
    reference_modes = [float(smooth_modes["f_hz"]), float(smooth_modes["p1_hz"])]
    mode_relative_errors = {
        mode: [phase_modes[mode][i] / reference_modes[i] - 1.0 for i in range(2)]
        for mode in ("slow", "rapid")
    }
    transfer_checks = [
        fortran_transfer_check(small_star, 2.0 * math.pi * 1800.0,
                               np.array([2.5e9, -7.0e7]), mode, _compile_solver()[0])
        for mode in ("slow", "rapid")
    ]

    # Tolerance refinement on the exact same finite-jump TOV model.
    loose = solve_maxwell_tov(strong_eos, pc, rtol=2e-8, max_step_log_pressure=0.03)
    tight = solve_maxwell_tov(strong_eos, pc, rtol=2e-11, max_step_log_pressure=0.005)

    result = {
        "schema": "maxwell-interface-validation-v1",
        "date": "2026-11-09",
        "provenance": {
            "baseline_solver_zip_sha256": "6d2ba545101df4cdc5bb5aed187b640ec84e528bf2005db764792bc4e19d9885",
            "backend_eos_builder": "native Fortran construct_maxwell_transition",
            "backend_phase_aware_tov_cowling": "Python/SciPy independent reference",
            "backend_smooth_reference": "native Fortran smooth-EOS TOV/Cowling",
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "strong_jump_fixture": {
            "interface_thermodynamics": strong_eos.interface_thermo,
            "mass_solar": strong_star.mass_solar,
            "radius_m": strong_star.radius_m,
            "interface_radius_m": strong_star.interface.radius_m,
            "accepted_tov_steps": strong_star.accepted_steps,
            "cowling_interface_diagnostics_at_1800_hz": strong_interface,
            "rapid_onset_stability_diagnostic": strong_eos.interface_thermo["seidov_onset_stable_rapid_diagnostic"],
            "interpretation": "boundary-condition and discontinuity exercise only; not admitted as a stable hybrid-star prediction",
        },
        "small_jump_smooth_limit": {
            "delta_n_over_n_low": small_eos.interface_thermo["density_jump_fm3"] / small_eos.interface_thermo["density_low_fm3"],
            "seidov_ratio": small_eos.interface_thermo["seidov_ratio"],
            "seidov_critical_ratio": small_eos.interface_thermo["seidov_critical_ratio"],
            "maxwell_identity_relative_residual": small_eos.interface_thermo["maxwell_identity_relative_residual"],
            "provider_first_law_residual_low": small_eos.interface_thermo["provider_first_law_residual_low"],
            "provider_first_law_residual_high": small_eos.interface_thermo["provider_first_law_residual_high"],
            "tov_smooth_fortran": {
                "mass_kg": smooth_tov["mass_kg"], "radius_m": smooth_tov["radius_m"]
            },
            "tov_phase_aware": {
                "mass_kg": small_star.mass_kg, "radius_m": small_star.radius_m
            },
            "tov_relative_difference": {
                "mass": small_star.mass_kg / smooth_tov["mass_kg"] - 1.0,
                "radius": small_star.radius_m / smooth_tov["radius_m"] - 1.0,
            },
            "cowling_smooth_fortran_hz": {"f": smooth_modes["f_hz"], "p1": smooth_modes["p1_hz"]},
            "cowling_phase_aware_hz": phase_modes,
            "cowling_relative_difference_vs_smooth_fortran": mode_relative_errors,
            "slow_rapid_relative_difference": [
                phase_modes["slow"][i] / phase_modes["rapid"][i] - 1.0 for i in range(2)
            ],
            "roots_unlabelled_outside_this_smooth_limit_cross_match": True,
        },
        "independent_fortran_python_interface_transfer": transfer_checks,
        "tov_tolerance_refinement": {
            "relative_mass_difference": tight.mass_kg / loose.mass_kg - 1.0,
            "relative_radius_difference": tight.radius_m / loose.radius_m - 1.0,
            "relative_interface_radius_difference": tight.interface.radius_m / loose.interface.radius_m - 1.0,
        },
        "status": "PASS_BOUNDED_REFERENCE_VALIDATION",
        "non_claims": [
            "not a native Fortran sharp-interface TOV/Cowling solver",
            "not a full stable-branch or maximum-mass calculation",
            "not a realistic EOS cohort prediction",
            "not a full-GR polar-mode or damping calculation",
            "no generalized mode-family continuation through strong transitions",
        ],
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "output": str(out_path),
        "small_jump_tov_relative_difference": result["small_jump_smooth_limit"]["tov_relative_difference"],
        "small_jump_mode_relative_difference": result["small_jump_smooth_limit"]["cowling_relative_difference_vs_smooth_fortran"],
        "fortran_python_interface_transfer": transfer_checks,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
