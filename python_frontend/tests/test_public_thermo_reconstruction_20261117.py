"""Native minimal-reconstruction controls; no generation or matching mocks."""
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from reconstruct_public_thermo_20261117 import (
    native_derivatives, reconstruct_data, hermite_interval_diagnostics,
)
from compose_eos.eos_data import EOSData


def analytic_source(pressure_factor=1.0, phase=None, q5=0.0):
    n = np.geomspace(.001, .1, 61)
    x = np.log(n/n[0])
    f = 930 + 2*x + .1*x*x + .01*x*x*x
    q = 2 + .2*x + .03*x*x
    z = np.zeros(len(n))
    data = EOSData(nb=n, mn=939.56533, mp=938.27209,
                   Q1=pressure_factor*q, Q2=z.copy(),
                   Q3=(f+pressure_factor*q)/939.56533-1,
                   Q4=z.copy(), Q5=z+q5, Q6=f/939.56533-1, Q7=f/939.56533-1,
                   q_extra=z.copy(), Iphase=np.ones(len(n), int) if phase is None else phase,
                   Yn=z.copy(), Yp=z.copy(), Ye=z.copy(), Ymu=z.copy(),
                   Aav=z.copy(), Zav=z.copy(), Yh=z.copy(), beta_equilibrated=q5 == 0.)
    return data, q


def test_native_cubic_derivative_independent_analytic_and_scipy():
    n = np.geomspace(.0002, .11, 71)
    x = np.log(n/n[0])
    binding = -8 + .3*x + .06*x*x + .01*x*x*x
    result = native_derivatives(n, binding)
    np.testing.assert_allclose(result[:, 0], .3+.12*x+.03*x*x, rtol=3e-12, atol=1e-12)
    np.testing.assert_allclose(result[:, 1], .12+.06*x, rtol=3e-11, atol=1e-11)
    np.testing.assert_allclose(result[:, 2], n*CubicSpline(n,binding)(n,1), rtol=5e-12)


def test_diagnosed_pressure_error_repaired_without_energy_or_gauge_change():
    source, exact_q = analytic_source(1.2)
    original_q7 = source.Q7.copy()
    rebuilt, report = reconstruct_data(source, model="analytic_control")
    np.testing.assert_allclose(rebuilt.Q1, exact_q, rtol=2e-10, atol=1e-10)
    np.testing.assert_array_equal(source.Q7, original_q7)
    np.testing.assert_array_equal(rebuilt.Q7, original_q7)
    assert report["changed_fields"] == ["Q1", "Q3"]
    assert report["physical_energy_offset_MeV_per_baryon"] == 0
    assert all(report["immutable_field_exact_preservation"].values())
    assert report["segments"][0]["pressure_inconsistency_diagnosed"]
    assert report["all_reconstructed_intervals_diagnostic_pass"]
    np.testing.assert_allclose(rebuilt.mn*(1+rebuilt.Q3), rebuilt.mn*(1+rebuilt.Q7)+rebuilt.Q1,
                               rtol=2e-16, atol=2e-13)


def test_consistent_public_table_is_exact_pass_through():
    source, _ = analytic_source()
    rebuilt, report = reconstruct_data(source)
    assert report["status"] == "UNCHANGED_THERMO"
    assert report["changed_fields"] == []
    np.testing.assert_array_equal(rebuilt.Q1, source.Q1)
    np.testing.assert_array_equal(rebuilt.Q3, source.Q3)


def test_algebraic_mu_error_only_does_not_change_pressure():
    source, _ = analytic_source()
    source.Q3 += .02/source.mn
    rebuilt, report = reconstruct_data(source)
    np.testing.assert_array_equal(rebuilt.Q1, source.Q1)
    np.testing.assert_array_equal(rebuilt.Q7, source.Q7)
    assert report["changed_fields"] == ["Q3"]
    assert report["segments"][0]["decision"] == "ALGEBRAIC_EULER_REPAIR_ONLY"


def test_phase_boundary_never_used_as_spline_interval():
    phase = np.r_[np.ones(30,int), np.full(31,7,int)]
    source, expected_q = analytic_source(1.2, phase)
    source.Q6[30:] += 1/source.mn
    source.Q7[30:] += 1/source.mn
    source.Q3[30:] += 1/source.mn
    rebuilt, report = reconstruct_data(source)
    np.testing.assert_allclose(rebuilt.Q1, expected_q, rtol=1e-9, atol=1e-9)
    np.testing.assert_array_equal(rebuilt.Q7, source.Q7)
    assert len(report["segments"]) == 2
    assert report["phase_boundaries_interpolated"] is False


def test_nonzero_q5_is_not_reset_or_reconstructed_under_unverified_premise():
    source, _ = analytic_source(1.2, q5=.001)
    rebuilt, report = reconstruct_data(source)
    for name in ("Q1", "Q3", "Q5", "Q7"):
        np.testing.assert_array_equal(getattr(rebuilt,name), getattr(source,name))
    assert report["status"] == "UNCHANGED_THERMO"
    assert not report["source_beta_convention_supported"]
    assert not report["ready_for_matching"]


def test_il_zero_charge_shift_preserved_in_repaired_q3():
    source, exact_q = analytic_source(1.2)
    source.Il = 0
    source.charge_fraction_grid = .1
    source.Q4[:] = -.1
    source.Q3 -= source.charge_fraction_grid*source.Q4
    rebuilt, report = reconstruct_data(source)
    mu = source.mn*(1+source.Q7)+exact_q
    np.testing.assert_allclose(rebuilt.mn*(1+rebuilt.Q3)+rebuilt.charge_fraction_grid*rebuilt.mn*rebuilt.Q4,
                               mu, rtol=2e-13)
    np.testing.assert_array_equal(rebuilt.Q4, source.Q4)


def test_bad_shape_is_reported_not_silently_clipped_or_removed():
    n=np.array([.01,.02,.03])
    mu=np.array([940.,941.,942.])
    p=np.array([.01,.011,.04])
    result=hermite_interval_diagnostics(n,mu-p/n,p/n)
    assert not result["all_intervals_pass"]
    assert result["mechanical_failure_intervals"]
    assert result["intervals"] == 2


def test_warm_data_cannot_be_cold_repaired():
    data,_=analytic_source(1.2)
    data.temperature_mev=.1
    with pytest.raises(ValueError, match="T=0"):
        reconstruct_data(data)


def test_named_inherited_outer_crust_is_preserved_until_internal_topology_resolved():
    source, _ = analytic_source(1.2, phase=np.zeros(61, int))
    rebuilt, report = reconstruct_data(source, model="SPG_M1")
    np.testing.assert_array_equal(rebuilt.Q1, source.Q1)
    np.testing.assert_array_equal(rebuilt.Q3, source.Q3)
    assert report["segments"][0]["decision"] == "UNCHANGED_INHERITED_OUTER_CRUST_TOPOLOGY"
