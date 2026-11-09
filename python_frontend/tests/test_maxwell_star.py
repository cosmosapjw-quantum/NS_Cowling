from __future__ import annotations

import sys
import subprocess
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import EOSData, construct_maxwell_transition
from compose_cowling import (MaxwellEOS, solve_maxwell_tov, solve_tov_star, solve_cowling_spectrum,
                              shoot_maxwell_cowling, scan_maxwell_cowling_modes)
from compose_cowling.maxwell import MEV_FM3_TO_J_M3, _interface_transfer, _sampled_phase_nodes
from compose_cowling.solver import _compile_solver, _parse_tsv


def quadratic_phase(*, rest_mass: float, stiffness: float, n_min: float,
                    n_max: float, count: int = 121) -> EOSData:
    """Analytic cold phase P=K n^2, eps=m n+P in MeV/fm^3."""
    nb = np.linspace(n_min, n_max, count)
    pressure = stiffness * nb**2
    epsilon = rest_mass * nb + pressure
    mu = (epsilon + pressure) / nb
    zeros = np.zeros_like(nb)
    q7 = epsilon / (nb * 939.0) - 1.0
    return EOSData(
        nb=nb, mn=939.0, mp=939.0,
        Q1=pressure / nb, Q2=zeros.copy(), Q3=mu / 939.0 - 1.0,
        Q4=zeros.copy(), Q5=zeros.copy(), Q6=q7.copy(), Q7=q7.copy(),
        q_extra=zeros.copy(), Iphase=np.zeros(nb.size, dtype=int),
        Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
        Aav=zeros.copy(), Zav=zeros.copy(), Yh=zeros.copy(),
        Il=1, iT=1, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_file_present=False, composition_available=False,
        beta_equilibrated=True,
    )


def reference_transition():
    low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4)
    high = quadratic_phase(rest_mass=940.0, stiffness=0.25, n_min=1.3, n_max=2.8)
    matched = construct_maxwell_transition(
        low, high, assume_local_charge_neutrality=True,
    )
    return low, high, matched


def test_sampled_phase_wv_nodes_are_separate_and_ignore_interface_sign_jump():
    w_inner = np.array([1.0, -1.0, 1.0, 1.0])
    w_outer = np.array([-1.0, 1.0, -1.0])
    v_inner = np.array([1.0, 1.0, -1.0, 1.0])
    v_outer = np.array([1.0, -1.0, -1.0])
    diag = _sampled_phase_nodes(w_inner, w_outer, v_inner, v_outer)

    assert diag["w_nodes_inner_sampled"] == 2
    assert diag["w_nodes_outer_sampled"] == 2
    assert diag["w_nodes_total_sampled"] == 4
    assert diag["v_nodes_inner_sampled"] == 2
    assert diag["v_nodes_outer_sampled"] == 1
    assert diag["v_nodes_total_sampled"] == 3

    # Changing only the sign across the interface is not an interior node.
    no_jump_w = _sampled_phase_nodes(
        np.ones(4), -np.ones(3), np.ones(4), -np.ones(3)
    )
    assert no_jump_w["w_nodes_total_sampled"] == 0
    assert no_jump_w["v_nodes_total_sampled"] == 0

    # Counts must survive arbitrary finite amplitude normalization and sign.
    scaled = _sampled_phase_nodes(
        -1e-24 * w_inner, -1e-24 * w_outer,
        1e-24 * v_inner, 1e-24 * v_outer,
    )
    assert scaled["w_nodes_total_sampled"] == diag["w_nodes_total_sampled"]
    assert scaled["v_nodes_total_sampled"] == diag["v_nodes_total_sampled"]


def test_maxwell_eos_keeps_two_one_sided_states_and_exact_thermodynamics():
    _low, _high, table = reference_transition()
    eos = MaxwellEOS(table)
    lo = eos.evaluate(eos.transition_pressure_pa, "low")
    hi = eos.evaluate(eos.transition_pressure_pa, "high")
    assert lo["pressure_pa"] == hi["pressure_pa"]
    assert hi["number_density_m3"] > lo["number_density_m3"]
    assert hi["energy_j_m3"] > lo["energy_j_m3"]
    assert lo["chemical_potential_j"] == pytest.approx(hi["chemical_potential_j"], rel=2e-8)
    lhs = hi["energy_j_m3"] - lo["energy_j_m3"]
    rhs = lo["chemical_potential_j"] * (hi["number_density_m3"] - lo["number_density_m3"])
    assert lhs == pytest.approx(rhs, rel=2e-8)
    assert eos.interface_thermo["maxwell_identity_relative_residual"] < 2e-8


def test_maxwell_tov_restarts_exactly_at_transition_and_preserves_metric():
    _low, _high, table = reference_transition()
    eos = MaxwellEOS(table)
    star = solve_maxwell_tov(
        eos, 1.35 * MEV_FM3_TO_J_M3,
        rtol=3e-9, max_step_log_pressure=0.015,
    )
    assert star.interface.radius_m > star.r0_m
    assert star.radius_m > star.interface.radius_m
    assert star.mass_kg > 0
    assert star.interface.energy_jump_j_m3 > 0
    assert star.interface.density_jump_m3 > 0
    bi = star.background(star.interface.radius_m, "inner")
    bo = star.background(star.interface.radius_m, "outer")
    assert bi["pressure_pa"] == pytest.approx(eos.transition_pressure_pa, rel=1e-13)
    assert bo["pressure_pa"] == pytest.approx(eos.transition_pressure_pa, rel=1e-13)
    assert bi["r_m"] == bo["r_m"]
    assert bi["m_kg"] == pytest.approx(bo["m_kg"], rel=1e-13)
    assert bi["nu"] == pytest.approx(bo["nu"], abs=1e-13)
    assert bi["lambda"] == pytest.approx(bo["lambda"], abs=1e-13)
    assert bi["dnu_dr"] == pytest.approx(bo["dnu_dr"], rel=1e-13)
    assert bi["energy_j_m3"] > bo["energy_j_m3"]
    assert bi["number_density_m3"] > bo["number_density_m3"]


def test_slow_and_rapid_cowling_interface_conditions_are_not_conflated():
    _low, _high, table = reference_transition()
    star = solve_maxwell_tov(
        MaxwellEOS(table), 1.35 * MEV_FM3_TO_J_M3,
        rtol=3e-9, max_step_log_pressure=0.015,
    )
    slow = shoot_maxwell_cowling(star, 1800.0, conversion="slow", max_step_m=120.0)
    rapid = shoot_maxwell_cowling(star, 1800.0, conversion="rapid", max_step_m=120.0)
    assert slow.interface_diagnostics["delta_p_relative_jump"] < 2e-10
    assert rapid.interface_diagnostics["delta_p_relative_jump"] < 2e-10
    assert slow.interface_diagnostics["slow_condition_residual"] < 2e-10
    assert rapid.interface_diagnostics["rapid_condition_residual"] < 2e-10
    assert abs(rapid.interface_diagnostics["radial_displacement_jump_m"]) > 0
    assert abs(slow.interface_diagnostics["isobar_displacement_jump_m"]) > 0
    assert np.isfinite(slow.residual)
    assert np.isfinite(rapid.residual)


def test_maxwell_tov_rejects_a_central_pressure_below_transition():
    _low, _high, table = reference_transition()
    eos = MaxwellEOS(table)
    with pytest.raises(ValueError, match="above Pt"):
        solve_maxwell_tov(eos, 0.9 * eos.transition_pressure_pa)


def test_fortran_and_python_interface_transfers_agree_for_both_conversion_limits(tmp_path):
    _low, _high, table = reference_transition()
    star = solve_maxwell_tov(
        MaxwellEOS(table), 1.35 * MEV_FM3_TO_J_M3, rtol=3e-9, max_step_log_pressure=0.015
    )
    case_driver, _ = _compile_solver()
    bi = star.background(star.interface.radius_m, "inner")
    bo = star.background(star.interface.radius_m, "outer")
    omega = 2.0 * np.pi * 1800.0
    w_in, v_in = 2.5e9, -7.0e7
    for conversion, case_id in (("slow", "material_transfer"), ("rapid", "rapid_transfer")):
        y_out, diagnostics = _interface_transfer(
            star, omega, np.array([w_in, v_in], dtype=float), conversion
        )
        nml = tmp_path / f"{case_id}.nml"
        nml.write_text(
            "&request request_token='cross-check' /\n"
            "&parameters\n"
            f"r_m={star.interface.radius_m:.17e}, omega_rad_s={omega:.17e}, "
            f"pressure_pa={star.interface.pressure_pa:.17e}, "
            f"eps_in={bi['energy_j_m3']:.17e}, eps_out={bo['energy_j_m3']:.17e}, "
            f"nu={bi['nu']:.17e}, lambda_value={bi['lambda']:.17e}, "
            f"dnu_dr={bi['dnu_dr']:.17e}, w_in_m3={w_in:.17e}, v_in_m2={v_in:.17e}, "
            "inner_side=1, outer_side=-1\n/\n",
            encoding="utf-8",
        )
        run = subprocess.run([str(case_driver), case_id, str(nml)], capture_output=True, text=True)
        assert run.returncode == 0, run.stderr
        result = _parse_tsv(run.stdout)
        assert result.get("status") == "ok", run.stdout
        assert result["w_out_m3"] == pytest.approx(y_out[0], rel=3e-12, abs=1e-8)
        assert result["v_out_m2"] == pytest.approx(y_out[1], rel=3e-12, abs=1e-8)
        assert result["delta_p_out_pa"] == pytest.approx(result["delta_p_in_pa"], rel=3e-12, abs=1e-5)
        if conversion == "rapid":
            assert result["rapid_condition_residual"] < 3e-12
        else:
            assert diagnostics["slow_condition_residual"] < 3e-12


def test_maxwell_tov_converges_under_tighter_ode_tolerances():
    _low, _high, table = reference_transition()
    eos = MaxwellEOS(table)
    loose = solve_maxwell_tov(eos, 1.35 * MEV_FM3_TO_J_M3, rtol=2e-8, max_step_log_pressure=0.03)
    tight = solve_maxwell_tov(eos, 1.35 * MEV_FM3_TO_J_M3, rtol=2e-11, max_step_log_pressure=0.005)
    assert tight.radius_m == pytest.approx(loose.radius_m, rel=2e-7)
    assert tight.mass_kg == pytest.approx(loose.mass_kg, rel=2e-7)
    assert tight.interface.radius_m == pytest.approx(loose.interface.radius_m, rel=2e-7)
    assert tight.background(tight.interface.radius_m, "inner")["m_kg"] == pytest.approx(
        loose.background(loose.interface.radius_m, "inner")["m_kg"], rel=2e-7
    )


def test_vanishing_density_jump_recovers_smooth_fortran_tov_solution():
    eta = 1.0e-4
    n_high_t = np.sqrt(1.0 + eta)
    k_high = 1.0 / (1.0 + eta)
    m_high = 941.0 - 2.0 * k_high * n_high_t
    low = quadratic_phase(rest_mass=939.0, stiffness=1.0, n_min=0.3, n_max=1.4, count=401)
    high = quadratic_phase(rest_mass=m_high, stiffness=k_high, n_min=0.8, n_max=2.8, count=401)
    table = construct_maxwell_transition(
        low, high, assume_local_charge_neutrality=True
    )
    hybrid_eos = MaxwellEOS(table)
    assert hybrid_eos.interface_thermo["density_jump_fm3"] / hybrid_eos.interface_thermo["density_low_fm3"] < 6e-5
    assert hybrid_eos.interface_thermo["seidov_onset_stable_rapid_diagnostic"]
    pc = 1.35 * MEV_FM3_TO_J_M3
    hybrid = solve_maxwell_tov(hybrid_eos, pc, rtol=2e-10, max_step_log_pressure=0.01)
    smooth = solve_tov_star(low, pc, rtol=1e-10)
    assert smooth["backend"] == "Fortran"
    assert hybrid.radius_m / smooth["radius_m"] - 1.0 == pytest.approx(0.0, abs=2e-5)
    assert hybrid.mass_kg / smooth["mass_kg"] - 1.0 == pytest.approx(0.0, abs=6e-5)

    # Independent mode benchmark: smooth native Fortran vs phase-aware SciPy
    # with a tiny but nonzero density jump and both physical conversion limits.
    smooth_modes = solve_cowling_spectrum(
        low, pc, [4.0e4, 6.0e4], [5.0e4, 7.0e4],
        target_mass_kg=smooth["mass_kg"], rtol=2e-8, mode_rtol=2e-8,
        step_cap=0.02, mass_rtol=2e-8, residual_tol=1e-5,
    )
    assert smooth_modes["backend"] == "Fortran"
    grid = np.array([4.0e4, 5.0e4, 6.0e4, 7.0e4])
    phase_mode_roots = {}
    for conversion in ("slow", "rapid"):
        roots = scan_maxwell_cowling_modes(
            hybrid, grid, conversion=conversion, rtol=2e-5, max_step_m=500.0
        )
        frequencies = sorted(float(row["frequency_hz"]) for row in roots)
        assert len(frequencies) >= 2
        phase_mode_roots[conversion] = frequencies[:2]
    reference_f_p1 = [smooth_modes["f_hz"], smooth_modes["p1_hz"]]
    for conversion in ("slow", "rapid"):
        for actual, reference in zip(phase_mode_roots[conversion], reference_f_p1):
            tolerance = 2e-5 if conversion == "slow" else 1e-4
            assert actual / reference - 1.0 == pytest.approx(0.0, abs=tolerance)
    for i in range(2):
        assert phase_mode_roots["slow"][i] / phase_mode_roots["rapid"][i] - 1.0 == pytest.approx(
            0.0, abs=1e-4
        )
