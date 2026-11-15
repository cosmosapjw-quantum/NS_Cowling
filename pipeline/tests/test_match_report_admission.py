"""Independent method-specific admission of native core--crust match reports.

The admission layer must never infer success from ``status == "ok"``. It
re-derives every interface condition from the *serialized* rows (internal CSV
read back from the CompOSE artefact) and cross-checks the native report.
Each corruption below must be rejected.
"""
from __future__ import annotations

import copy
import csv
import math

import pytest

from scripts import match_report_admission as adm

MN = 939.56533
HEADER = [
    "nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "q_extra",
    "Iphase", "Yn", "Yp", "Ye", "Ymu", "Abar", "Zbar", "Yh",
    "row_status", "yn_present", "yp_present", "ye_present", "ymu_present",
    "abar_present", "zbar_present", "yh_present", "c3_applicable",
    "c4_applicable",
]


def _css_state(n0, mu0, p0, cs2, mu):
    """Exact CSS state P-P0 = cs2/(1+cs2) (n mu - n0 mu0), n = n0 (mu/mu0)^(1/cs2)."""
    n = n0 * (mu / mu0) ** (1.0 / cs2)
    p = p0 + cs2 / (1.0 + cs2) * (n * mu - n0 * mu0)
    return n, p, n * mu - p


def _crust_states():
    """Euler-consistent, causal, mechanically stable crust branch (polytrope)."""
    out = []
    for n in (0.01, 0.03, 0.06, 0.09, 0.12, 0.145):
        k, g = 2.0, 1.6
        p = k * n ** g
        e_over_n = MN - 8.0 + k * n ** (g - 1.0) / (g - 1.0)
        out.append((n, p, n * e_over_n))
    return out


def _row(n, p, eps, *, mu_q3=None):
    mu_e = (eps + p) / n
    mu = mu_e if mu_q3 is None else mu_q3
    f = eps / n
    return [n, p / n, 0.0, mu / MN - 1.0, 0.0, 0.0, f / MN - 1.0, f / MN - 1.0, 0.0,
            0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]


def _css_case():
    crust = _crust_states()
    n_t, p_t, eps_l = crust[-1]
    mu_t = (eps_l + p_t) / n_t
    cs2 = 0.25
    n_high = n_t * 1.04
    eps_high = mu_t * n_high - p_t
    mu_core = mu_t * 1.12
    bridge_mu = [mu_t * (mu_core / mu_t) ** (j / 9.0) for j in range(1, 9)]
    bridge = [_css_state(n_high, mu_t, p_t, cs2, mu) for mu in bridge_mu]
    n_core, p_core, eps_core = _css_state(n_high, mu_t, p_t, cs2, mu_core)
    # Core continues as a different CSS (cs2=0.5) from the join point.
    core = [(n_core, p_core, eps_core)] + [
        _css_state(n_core, mu_core, p_core, 0.5, mu_core * (1.0 + 0.03 * j))
        for j in range(1, 5)
    ]
    states = crust + [(n_high, p_t, eps_high)] + bridge + core
    rows = [_row(*s) for s in states]
    i_low = len(crust) - 1
    report = {
        "status": "ok", "operation": "match-core", "n": len(rows), "method": "css",
        "construction_kind": "css_single_speed_with_crust_density_jump",
        "spline_mode": "single_speed_css_analytic", "n_bridge": 8,
        "n1": n_t, "n2": n_core, "P1": p_t, "P2": p_core,
        "F1": eps_l / n_t, "F2": eps_core / n_core, "mu1": mu_t, "mu2": mu_core,
        "cs2_left": cs2, "cs2_right": cs2, "cs2_min": cs2, "cs2_max": cs2,
        "css_single_sound_speed": cs2,
        "css_transition_nb_low": n_t, "css_transition_nb_high": n_high,
        "css_transition_delta_nb": n_high - n_t,
        "css_transition_delta_eps": eps_high - eps_l,
        "css_transition_epsilon_low": eps_l, "css_transition_epsilon_high": eps_high,
        "css_transition_pressure": p_t, "css_transition_mu": mu_t,
        "css_transition_low_index": i_low, "css_transition_high_index": i_low + 1,
        "css_core_bridge_index": i_low + 1 + 8 + 1,
        "first_law_relative_residual": 1.0e-14,
        "composition_dropped": True, "is_exact_maxwell": True,
    }
    return report, rows


def _maxwell_case():
    crust = _crust_states()
    n_a, p_t, eps_a = crust[-1]
    mu_t = (eps_a + p_t) / n_a
    n_b = n_a * 1.3
    eps_b = mu_t * n_b - p_t
    high = [(n_b, p_t, eps_b)] + [
        _css_state(n_b, mu_t, p_t, 0.4, mu_t * (1.0 + 0.02 * j)) for j in range(1, 7)
    ]
    rows = [_row(*s) for s in crust + high]
    report = {
        "status": "ok", "operation": "maxwell-transition", "n": len(rows),
        "pressure_crossings": 1, "low_phase_rows": len(crust),
        "high_phase_rows": len(high), "mu_transition": mu_t,
        "pressure_transition": p_t, "n_low": n_a, "n_high": n_b,
        "epsilon_low": eps_a, "epsilon_high": eps_b, "density_jump": n_b - n_a,
        "latent_energy_density": eps_b - eps_a, "maxwell_identity_residual": 0.0,
        "max_euler_residual_low": 1.0e-12, "max_euler_residual_high": 1.0e-12,
        "min_cs2_low": 0.001, "max_cs2_low": 0.05, "min_cs2_high": 0.3,
        "max_cs2_high": 0.45,
    }
    return report, rows


META = {"mn": MN, "il": 1, "charge_fraction_grid": 0.0}


def _failed(result, name):
    assert result["status"] == "REJECTED", result
    assert name in result["failed_checks"], result["failed_checks"]


# ----------------------------------------------------------------------- CSS

def test_css_admits_exact_analytic_construction():
    report, rows = _css_case()
    result = adm.admit_css(report, adm.rows_from_table(rows, META), META)
    assert result["status"] == "ADMITTED", result["failed_checks"]
    assert result["schema"] == adm.SCHEMA
    assert not result["failed_checks"]


def test_css_rejects_wrong_construction_kind():
    report, rows = _css_case()
    report["construction_kind"] = "smooth_near_pressure_plateau"
    _failed(adm.admit_css(report, adm.rows_from_table(rows, META), META), "construction_kind")


@pytest.mark.parametrize("cs2", [0.0, -0.1, 1.2, float("nan")])
def test_css_rejects_noncausal_or_nonpositive_sound_speed(cs2):
    report, rows = _css_case()
    for key in ("css_single_sound_speed", "cs2_min", "cs2_max", "cs2_left", "cs2_right"):
        report[key] = cs2
    _failed(adm.admit_css(report, adm.rows_from_table(rows, META), META), "css_sound_speed_range")


def test_css_rejects_more_than_one_sound_speed_in_report():
    report, rows = _css_case()
    report["cs2_max"] = report["cs2_min"] * 1.01
    _failed(adm.admit_css(report, adm.rows_from_table(rows, META), META), "css_single_sound_speed")


def test_css_rejects_serialized_segment_with_second_sound_speed():
    report, rows = _css_case()
    k = report["css_transition_high_index"] + 4
    n, p = rows[k][0], rows[k][1] * rows[k][0]
    eps = n * MN * (rows[k][7] + 1.0)
    rows[k] = _row(n * 1.001, p, eps * 1.001)
    _failed(adm.admit_css(report, adm.rows_from_table(rows, META), META), "css_segment_single_speed")


def test_css_rejects_crust_side_pressure_discontinuity():
    report, rows = _css_case()
    k = report["css_transition_high_index"]
    rows[k][1] *= 1.0 + 1.0e-6
    _failed(adm.admit_css(report, adm.rows_from_table(rows, META), META), "crust_transition_pressure")


def test_css_rejects_corrupted_latent_energy_jump():
    report, rows = _css_case()
    k = report["css_transition_high_index"]
    # Shift the bridge-side energy per baryon: breaks Delta eps = mu_t Delta n
    rows[k][6] += 1.0e-5
    rows[k][7] += 1.0e-5
    result = adm.admit_css(report, adm.rows_from_table(rows, META), META)
    assert result["status"] == "REJECTED"
    assert {"crust_transition_latent_identity", "crust_transition_euler_mu"} & set(
        result["failed_checks"])


def test_css_rejects_nonpositive_density_jump():
    report, rows = _css_case()
    lo, hi = report["css_transition_low_index"], report["css_transition_high_index"]
    rows[hi] = list(rows[lo])
    _failed(adm.admit_css(report, adm.rows_from_table(rows, META), META), "crust_transition_positive_jump")


def test_css_rejects_core_endpoint_discontinuity():
    report, rows = _css_case()
    k = report["css_core_bridge_index"]
    rows[k][0] *= 1.0 + 1.0e-7
    _failed(adm.admit_css(report, adm.rows_from_table(rows, META), META), "css_core_continuity")


def test_css_rejects_first_law_residual_above_declared_gate():
    report, rows = _css_case()
    report["first_law_relative_residual"] = 6.0e-9
    _failed(adm.admit_css(report, adm.rows_from_table(rows, META), META), "first_law_relative_residual")


def test_css_rejects_report_rows_disagreement():
    report, rows = _css_case()
    report["css_transition_pressure"] *= 1.0 + 1.0e-6
    _failed(adm.admit_css(report, adm.rows_from_table(rows, META), META), "report_matches_rows")


def test_css_fails_closed_on_missing_report_field():
    report, rows = _css_case()
    del report["css_core_bridge_index"]
    result = adm.admit_css(report, adm.rows_from_table(rows, META), META)
    assert result["status"] == "REJECTED"
    assert "report_fields_present" in result["failed_checks"]


def test_status_ok_alone_is_never_admitted():
    result = adm.admit_css({"status": "ok"}, [], META)
    assert result["status"] == "REJECTED"
    result = adm.admit_maxwell({"status": "ok"}, [], META)
    assert result["status"] == "REJECTED"


# ------------------------------------------------------------------- Maxwell

def test_maxwell_admits_genuine_coexistence():
    report, rows = _maxwell_case()
    result = adm.admit_maxwell(report, adm.rows_from_table(rows, META), META)
    assert result["status"] == "ADMITTED", result["failed_checks"]


@pytest.mark.parametrize("crossings", [0, 2])
def test_maxwell_requires_exactly_one_pressure_crossing(crossings):
    report, rows = _maxwell_case()
    report["pressure_crossings"] = crossings
    _failed(adm.admit_maxwell(report, adm.rows_from_table(rows, META), META), "pressure_crossing")


def test_maxwell_rejects_chemical_potential_mismatch():
    report, rows = _maxwell_case()
    k = report["low_phase_rows"]
    rows[k][6] += 1.0e-4
    rows[k][7] += 1.0e-4
    result = adm.admit_maxwell(report, adm.rows_from_table(rows, META), META)
    assert result["status"] == "REJECTED"
    assert "coexistence_euler_mu" in result["failed_checks"]


def test_maxwell_rejects_nonpositive_jump():
    report, rows = _maxwell_case()
    k = report["low_phase_rows"]
    rows[k] = list(rows[k - 1])
    _failed(adm.admit_maxwell(report, adm.rows_from_table(rows, META), META), "coexistence_positive_jump")


def test_maxwell_rejects_acausal_branch():
    report, rows = _maxwell_case()
    # Steepen the last high-branch step so dP/deps > 1.
    n = rows[-1][0]
    n2, p2 = rows[-2][0], rows[-2][1] * rows[-2][0]
    eps2 = n2 * MN * (rows[-2][7] + 1.0)
    eps = n * MN * (rows[-1][7] + 1.0)
    rows[-1] = _row(n, p2 + 1.5 * (eps - eps2), eps)
    _failed(adm.admit_maxwell(report, adm.rows_from_table(rows, META), META), "high_branch_causal")


def test_maxwell_rejects_mechanically_unstable_branch():
    report, rows = _maxwell_case()
    n = rows[2][0]
    eps = n * MN * (rows[2][7] + 1.0)
    rows[2] = _row(n, rows[1][1] * rows[1][0] * 0.9, eps)
    _failed(adm.admit_maxwell(report, adm.rows_from_table(rows, META), META), "low_branch_stable")


def test_maxwell_rejects_report_euler_residual_above_native_gate():
    report, rows = _maxwell_case()
    report["max_euler_residual_high"] = 6.0e-4
    _failed(adm.admit_maxwell(report, adm.rows_from_table(rows, META), META), "report_euler_residuals")


def test_maxwell_rejects_serialized_q3_euler_violation():
    report, rows = _maxwell_case()
    rows[1][3] += 1.0e-2
    _failed(adm.admit_maxwell(report, adm.rows_from_table(rows, META), META), "serialized_euler_mu")


def test_maxwell_recomputes_identity_instead_of_trusting_report():
    report, rows = _maxwell_case()
    k = report["low_phase_rows"]
    # Report claims zero residual but the serialized jump violates it.
    n, p = rows[k][0], rows[k][1] * rows[k][0]
    eps = n * MN * (rows[k][7] + 1.0)
    rows[k] = _row(n, p, eps * (1.0 + 1.0e-6), mu_q3=(eps + p) / n)
    report["maxwell_identity_residual"] = 0.0
    result = adm.admit_maxwell(report, adm.rows_from_table(rows, META), META)
    assert result["status"] == "REJECTED"
    assert {"coexistence_latent_identity", "coexistence_euler_mu"} & set(result["failed_checks"])


# ----------------------------------------------------------- I/O + dispatch

def test_il0_charge_shift_is_applied_to_q3_chemical_potential():
    meta = {"mn": MN, "il": 0, "charge_fraction_grid": 0.3}
    rows = [_row(0.1, 1.0, 94.0)]
    mu_q = -50.0
    mu_e = (94.0 + 1.0) / 0.1
    rows[0][3] = (mu_e - 0.3 * mu_q) / MN - 1.0
    rows[0][4] = mu_q / MN
    parsed = adm.rows_from_table(rows, meta)
    assert math.isclose(parsed[0]["mu_q3"], mu_e, rel_tol=1e-14)


def test_admit_from_files_dispatches_by_method(tmp_path):
    report, rows = _css_case()
    csv_path = tmp_path / "eos_data.csv"
    nml_path = tmp_path / "eos_data.nml"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(HEADER)
        for r in rows:
            w.writerow([repr(v) if isinstance(v, float) else v for v in r])
    nml_path.write_text(
        "&metadata\n"
        f" mn={MN:.16E}\n il=1\n charge_fraction_grid=0.0000000000000000E+000\n/\n",
        encoding="utf-8",
    )
    result = adm.admit_match_report("css", report, csv_path, nml_path)
    assert result["status"] == "ADMITTED", result["failed_checks"]
    result = adm.admit_match_report("maxwell", report, csv_path, nml_path)
    assert result["status"] == "REJECTED"
    bad = adm.admit_match_report("cubic_F", report, csv_path, nml_path)
    assert bad["status"] == "REJECTED"
    assert "method_supported" in bad["failed_checks"]


def test_admission_does_not_mutate_report():
    report, rows = _css_case()
    snapshot = copy.deepcopy(report)
    adm.admit_css(report, adm.rows_from_table(rows, META), META)
    assert report == snapshot


def test_failure_taxonomy_separates_provenance_contract_physics_and_implementation():
    assert adm.classify_failure("source", "missing CompOSE source files: eos.nb", None) == "data_provenance"
    assert adm.classify_failure(
        "maxwell", "phase tables must use the same neutron and proton rest-mass conventions", None
    ) == "implementation_convention"
    assert adm.classify_failure(
        "css", "native matching requires beta-equilibrated cold barotropes", None
    ) == "input_contract"
    assert adm.classify_failure(
        "maxwell", "no pressure equality was found in the common mu_B interval", None
    ) == "physical"
    assert adm.classify_failure(
        "css", "single-speed CSS candidate rejected: no causal single-speed CSS "
        "solution (0<cs2<=1) matches the endpoints", None
    ) == "physical"
    assert adm.classify_failure(
        "css", "matching candidate violates the endpoint first-law integral tolerance", None
    ) == "numerical"
    assert adm.classify_failure("readback", "row 3 Q1 mismatch", None) == "implementation"
    assert adm.classify_failure("admission", "css_core_continuity", None) == "admission"
    assert adm.classify_failure("css", "something entirely new", None) == "unclassified"
