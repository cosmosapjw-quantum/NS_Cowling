import copy
import json
from types import SimpleNamespace

import numpy as np
import pytest

from compose_eos.physical_domain import physical_domain_report


def table(pressure=(0.1, 0.12, 0.15, 0.17), energy=(1.0, 1.1, 1.2, 1.3)):
    n = np.array([1.0, 1.1, 1.2, 1.3])
    p, e = np.array(pressure), np.array(energy)
    zeros = np.zeros(4)
    return SimpleNamespace(
        nb=n, p=p, eps=e, f=e, mub=(e+p)/n, muq=zeros.copy(), mul=zeros.copy(),
        mn=939.0, inb_start=17, Il=1, temperature_mev=0.0, charge_fraction_grid=0.0,
        Iphase=np.zeros(4, dtype=int), Ye=zeros.copy(), Ymu=zeros.copy(),
        Ye_present=np.zeros(4, dtype=bool), Ymu_present=np.zeros(4, dtype=bool),
        row_status=np.zeros(4, dtype=int), C3_applicable=np.zeros(4, dtype=bool),
        C4_applicable=np.zeros(4, dtype=bool), composition_available=False,
        composition_file_present=False,
    )


def test_acausal_tail_retains_lower_density_support():
    data = table(energy=(1, 1.12, 1.25, 1.27), pressure=(.1, .15, .2, .3))
    report = physical_domain_report(data)
    assert report["mechanically_causal_segments"] == [
        dict(first_row=1, last_row=3, n_min_fm3=1.0, n_max_fm3=1.2)]
    assert report["first_interval_findings"]["ACAUSAL_SECANT"]["left_row"] == 3
    assert len(report["rows"]) == 4
    assert report["rows"][0]["source_density_row_label"] == 17
    assert report["composition"]["status_counts"]["NOT_AVAILABLE"] == 4


def test_gap_and_phase_changes_split_segments():
    data = table(energy=(1, 1.12, 1.25, 1.4), pressure=(.1, .15, .2, .3))
    data.Iphase[:] = [0, 0, 1, 1]
    report = physical_domain_report(data)
    assert len(report["mechanically_causal_segments"]) == 2
    assert "PHASE_BOUNDARY" in report["intervals"][1]["reasons"]
    data.Iphase[:] = 0
    data.nb[:] = [1, 1.1, 2.0, 2.1]
    data.mub = (data.eps + data.p) / data.nb
    report = physical_domain_report(data)
    assert "COARSE_SAMPLING" in report["intervals"][1]["reasons"]
    assert report["intervals"][1]["continuity_status"] == "UNRESOLVED"


def test_missing_leptons_are_unavailable_and_do_not_reject_mechanical_domain():
    data = table(energy=(1, 1.12, 1.25, 1.4), pressure=(.1, .15, .2, .3))
    data.mul[:] = 20.0
    report = physical_domain_report(data)
    assert report["euler_unavailable_rows"] == [1, 2, 3, 4]
    assert report["mechanically_causal_segments"]
    assert not report["thermodynamically_supported_segments"]


def test_charge_term_is_used_in_il_zero_euler_identity():
    data = table(energy=(1, 1.12, 1.25, 1.4), pressure=(.1, .15, .2, .3))
    data.Il = 0
    data.charge_fraction_grid = .1
    data.muq[:] = .2
    data.mub -= .02
    report = physical_domain_report(data)
    assert not report["euler_failed_rows"]


def test_nonfinite_rows_remain_strict_json_and_split_support():
    data = table(energy=(1, 1.12, 1.25, 1.4), pressure=(.1, .15, .2, .3))
    data.p[2] = np.nan
    report = physical_domain_report(data)
    assert report["rows"][2]["P_MeV_fm3"] is None
    assert report["mechanically_causal_segments"][0]["last_row"] == 2
    json.dumps(report, allow_nan=False)


def test_domain_diagnostics_do_not_mutate_source():
    data = table(energy=(1, 1.12, 1.25, 1.4), pressure=(.1, .15, .2, .3))
    before = copy.deepcopy(data)
    physical_domain_report(data)
    for name, value in vars(data).items():
        if isinstance(value, np.ndarray):
            np.testing.assert_array_equal(value, getattr(before, name))


@pytest.mark.parametrize("density", [[1, 1, 2, 3], [1, 0, 2, 3], [1, np.nan, 2, 3]])
def test_invalid_density_fails_explicitly(density):
    data = table()
    data.nb = np.array(density)
    with pytest.raises(ValueError, match="density grid"):
        physical_domain_report(data)

def test_coarse_grid_alone_preserves_secant_diagnostic():
    data = table()
    data.nb = np.array([1., 2., 3., 4.])
    data.eps = 2 * data.nb
    data.f = data.eps.copy()
    data.p = .2 * data.nb
    data.mub = (data.eps + data.p) / data.nb
    report = physical_domain_report(data)
    assert all(row["secant_diagnostic_pass"] for row in report["intervals"])
    assert all(not row["physical_findings"] for row in report["intervals"])
    assert report["sampled_secant_pass_count"] == 3

def test_cold_energy_free_energy_mismatch_excludes_thermodynamic_support():
    data = table(energy=(1, 1.12, 1.25, 1.4), pressure=(.1, .15, .2, .3))
    data.eps = data.eps + .1
    report = physical_domain_report(data)
    assert report["mechanically_causal_segments"]
    assert not report["thermodynamically_supported_segments"]
    assert all(row["cold_free_energy_status"] == "FAIL" for row in report["rows"])


def test_charged_path_mu_b_decrease_is_separate_from_compressibility():
    data = table(energy=(1, 1.12, 1.25, 1.4), pressure=(.1, .15, .2, .3))
    data.Il = 0
    data.charge_fraction_grid = .5
    data.muq[:] = [.0, .2, .4, .6]
    data.mub = (data.eps + data.p) / data.nb - .5 * data.muq
    report = physical_domain_report(data)
    assert all(item["secant_diagnostic_pass"] for item in report["intervals"])
    assert any(item["baryon_chemical_potential_decreases"] for item in report["intervals"])
    assert not report["euler_failed_rows"]

def test_maxwell_pressure_roundoff_has_separate_machine_precision_allowance():
    data = table()
    data.p = np.array([1., 1., np.nextafter(1., 0.), 1.])
    data.eps = 3 * data.nb - data.p
    data.f = data.eps.copy()
    data.mub = (data.eps + data.p) / data.nb
    data.Iphase[:] = [6, 6, 7, 7]
    report = physical_domain_report(data)
    jump = report["intervals"][1]
    assert jump["physical_findings"] == []
    assert jump["secant_diagnostic_pass"]
    assert jump["pressure_monotonic_roundoff_MeV_fm3"] > abs(data.p[2]-data.p[1])
    assert jump["continuity_findings"] == ["PHASE_BOUNDARY"]


def test_small_pressure_real_decrease_is_not_swallowed_by_absolute_floor():
    data = table()
    data.p = np.array([1e-20, 1.1e-20, .9e-20, 1.2e-20])
    data.eps = 3 * data.nb
    data.f = data.eps.copy()
    data.mub = (data.eps + data.p) / data.nb
    report = physical_domain_report(data)
    assert "MECHANICAL_INSTABILITY" in report["intervals"][1]["physical_findings"]
