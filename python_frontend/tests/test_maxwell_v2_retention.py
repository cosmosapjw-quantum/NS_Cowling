"""Lossless source retention and truthful evidence for segmented Maxwell EOSs."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import EOSData, NativeExecutionError, construct_maxwell_transition, diagnose_maxwell_phase
from compose_eos.maxwell_v2 import (
    _local_neutrality_residual, _segment_eosdata, construct_maxwell_transition_v2,
)


def phase(rest_mass, stiffness, n_min, n_max, count=81):
    n = np.linspace(n_min, n_max, count)
    p = stiffness * n**2
    e = rest_mass * n + p
    z = np.zeros_like(n)
    return EOSData(
        nb=n, mn=939., mp=938., Q1=p/n, Q2=z.copy(),
        Q3=(e+p)/(n*939.)-1, Q4=z.copy(), Q5=z.copy(),
        Q6=e/(n*939.)-1, Q7=e/(n*939.)-1, q_extra=z.copy(),
        Iphase=np.zeros(count, int), Yn=z.copy(), Yp=z.copy(), Ye=z.copy(),
        Ymu=z.copy(), Aav=z.copy(), Zav=z.copy(), Yh=z.copy(),
        Il=1, iT=0, iYq=1, temperature_mev=0.,
        composition_available=False, beta_equilibrated=True,
    )


def pair():
    return phase(939., 1., .3, 1.4), phase(940., .25, 1.3, 2.8)


def test_segment_keeps_support_axes_and_beta_evidence():
    low, _ = pair()
    low.source_dir = "/retained/source/eos"
    low.beta_equilibrated = False
    low.inb_start = 19
    low.nT = low.nYq = 1
    low.t_first = low.t_last = low.selected_iT = 0
    low.yq_first = low.yq_last = low.selected_iYq = 1
    low.t_axis = np.array([0.])
    low.yq_axis = np.array([0.])
    segment = _segment_eosdata(low, (0, 20, 60))
    assert segment.source_dir == low.source_dir
    assert segment.beta_equilibrated is False
    assert segment.inb_start == 39
    assert segment.nT == 1 and segment.selected_iYq == 1
    assert segment.row_status.size == 40
    np.testing.assert_array_equal(segment.t_axis, low.t_axis)


def test_missing_composition_is_not_zero_neutrality_evidence():
    low, _ = pair()
    assert _local_neutrality_residual(low, slice(None)) is None


def test_v2_requires_explicit_local_neutrality_assumption():
    low, high = pair()
    with pytest.raises(NativeExecutionError) as caught:
        construct_maxwell_transition_v2(low, high)
    assert caught.value.category == "LOCAL_NEUTRALITY_NOT_CONFIRMED"


def test_native_error_keeps_category_and_unadmitted_report(monkeypatch):
    import compose_eos.maxwell_v2 as module
    low, high = pair()
    def reject(*args, **kwargs):
        raise NativeExecutionError("finite-domain crossing absent", category="NO_COEXISTENCE_CROSSING",
                                   operation="maxwell-transition")
    monkeypatch.setattr(module, "construct_maxwell_transition", reject)
    with pytest.raises(NativeExecutionError) as caught:
        construct_maxwell_transition_v2(low, high, assume_local_charge_neutrality=True)
    assert caught.value.category == "NO_COEXISTENCE_CROSSING"
    assert caught.value.report.admitted_pairs == []
    assert caught.value.report.attempted_pairs == [(0, 0)]


def test_v2_retains_both_source_tails_and_records_original_phase_labels():
    low, high = pair()
    low.Iphase[:20] = 4
    low.Iphase[20:] = 7
    high.Iphase[:60] = 3
    high.Iphase[60:] = 6
    # Header reference changes must not change any physical source column.
    old_mn = high.mn
    high.mn = 940.
    for name in ("Q3", "Q6", "Q7"):
        setattr(high, name, (getattr(high, name)+1)*old_mn/high.mn-1)
    for name in ("Q4", "Q5"):
        setattr(high, name, getattr(high, name)*old_mn/high.mn)
    out, report = construct_maxwell_transition_v2(
        low, high, assume_local_charge_neutrality=True)
    assert out.nb[0] == low.nb[0] and out.nb[-1] == high.nb[-1]
    np.testing.assert_array_equal(out.nb[:20], low.nb[:20])
    np.testing.assert_array_equal(out.nb[-21:], high.nb[60:])
    np.testing.assert_allclose(out.p[:20], low.p[:20], rtol=0, atol=0)
    np.testing.assert_allclose(out.eps[-21:], high.eps[60:], rtol=3e-15)
    np.testing.assert_allclose(out.mub[-21:], high.mub[60:], rtol=3e-15)
    assert report.admitted_pairs == [(1, 0)]
    assert report.native_report["pressure_transition"] == pytest.approx(1., rel=2e-10)
    a, b = report.transition_row_indices_zero_based
    assert b == a+1
    assert out.nb[a] == pytest.approx(1.) and out.nb[b] == pytest.approx(2.)
    assert out.p[a] == pytest.approx(out.p[b], rel=2e-11)
    assert [m["source_iphase"] for m in report.output_phase_map] == [4, 7, 3, 6]
    assert report.retained_source_rows == {"low_prefix": 20, "high_suffix": 21}
    assert report.inherited_interfaces_status == "PRESERVED_SOURCE_BOUNDARIES_NO_CROSS_PHASE_INTERPOLATION"
    assert report.local_neutrality_evidence_low["segment_0_iphase_4"]["charge_neutrality_residual"] is None


def test_v2_return_report_false_returns_only_eosdata():
    low, high = pair()
    out = construct_maxwell_transition_v2(
        low, high, assume_local_charge_neutrality=True, return_report=False)
    assert isinstance(out, EOSData)


def test_native_retained_source_thermodynamics_are_not_silently_reconstructed():
    low, high = pair()
    # A small source inconsistency lies inside the existing native acceptance
    # tolerance. It must stay visible on unchanged rows, rather than silently
    # turning the public-original cohort into a reconstruction cohort.
    low.Q3[4] += 1e-6
    low.Q2[4] = 1e-12
    low.Q5[4] = 1e-10
    low.Q6[4] += 1e-12
    out = construct_maxwell_transition(low, high, assume_local_charge_neutrality=True)
    hit = np.flatnonzero(out.nb == low.nb[4])
    assert hit.size == 1
    for name in ("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7"):
        assert getattr(out, name)[hit[0]] == pytest.approx(getattr(low, name)[4], rel=0, abs=2e-16)


def test_native_interval_diagnostic_is_local_and_does_not_bridge_phase_boundaries():
    low, _ = pair()
    low.Iphase[20:] = 7
    result = diagnose_maxwell_phase(low)
    assert result["global_validation_category"] == "COMPOSITION_UNSUPPORTED"
    assert result["source_mutated"] is False
    invalid = [r for r in result["intervals"] if not r["valid"]]
    assert len(invalid) == 1
    assert invalid[0]["left_row_zero_based"] == 19
    assert invalid[0]["category"] == "INHERITED_PHASE_BOUNDARY"
    assert all(r["max_cs2"] < 1 for r in result["intervals"] if r["valid"])


def test_native_interval_diagnostic_does_not_hide_local_beta_failure():
    low, _ = pair()
    low.Q5[30] = 1e-3
    low.beta_equilibrated = False
    result = diagnose_maxwell_phase(low)
    assert result["source_beta_equilibrated"] is False
    assert result["global_validation_category"] == "BETA_EQUILIBRIUM_FAIL"
    invalid = [r for r in result["intervals"] if not r["valid"]]
    assert [r["left_row_zero_based"] for r in invalid] == [29, 30]
    assert all(r["category"] == "BETA_EQUILIBRIUM_FAIL" for r in invalid)


def test_native_interval_diagnostic_reports_acausal_support_without_altering_it():
    low, _ = pair()
    n = low.nb.copy()
    low.Q1 = 1000*n**2
    low.Q7 = (939+500*n**2)/939-1
    low.Q6 = low.Q7.copy()
    low.Q3 = (939+1500*n**2)/939-1
    result = diagnose_maxwell_phase(low)
    valid = [r for r in result["intervals"] if r["valid"]]
    invalid = [r for r in result["intervals"] if not r["valid"]]
    assert valid and invalid
    assert all(r["category"] == "CAUSALITY_FAIL" for r in invalid)
    assert valid[-1]["n_right_fm3"] < (939/1500)**.5
    np.testing.assert_array_equal(low.nb, n)


def test_native_slice_convention_error_is_not_beta_failure():
    low, high = pair()
    high.Il = 0
    with pytest.raises(NativeExecutionError) as caught:
        construct_maxwell_transition(low, high, assume_local_charge_neutrality=True)
    assert caught.value.category == "SLICE_CONVENTION_MISMATCH"


def roundoff_pair():
    # Exact quadratic barotropes with a prescribed coexistence density. The
    # parameters reproduce a one-ulp decrease solely from (P_t/n)*n in the
    # two output rows; no source interval has decreasing pressure.
    ratio = 1.3087034113022535
    stiffness = 65.50878370948828
    n_low = 0.014567000677367797
    m_high = 939.+2*stiffness*n_low*(1-1/ratio)
    return (phase(939., stiffness, n_low*.5, n_low*1.5),
            phase(m_high, stiffness/ratio**2, n_low*ratio*.6, n_low*ratio*1.6))


def test_native_maxwell_allows_only_interface_pressure_reencoding_roundoff():
    low, high = roundoff_pair()
    out, report = construct_maxwell_transition(
        low, high, assume_local_charge_neutrality=True, return_report=True)
    i = int(report["low_phase_rows"])-1
    pressure_difference = out.p[i+1]-out.p[i]
    bound = 64*np.finfo(float).eps*max(abs(out.p[i]), abs(out.p[i+1]))
    assert abs(pressure_difference) <= bound
    assert report["interface_pressure_encoding_difference"] == pressure_difference
    assert report["interface_pressure_roundoff_bound"] == pytest.approx(bound, rel=1e-14)
    within_phase = out.Iphase[1:] == out.Iphase[:-1]
    assert np.all(np.diff(out.p)[within_phase] > 0)
    assert np.all(np.diff(out.nb) > 0) and np.all(np.diff(out.eps) > 0)
    for original in (low, high):
        for source_row in range(original.n):
            hits = np.flatnonzero(out.nb == original.nb[source_row])
            if hits.size:
                assert out.Q1[hits[0]] == original.Q1[source_row]


def test_native_maxwell_still_rejects_decreasing_source_pressure():
    low, high = roundoff_pair()
    previous_pressure = low.p[9]
    low.Q1[10] = np.nextafter(previous_pressure, -np.inf)/low.nb[10]
    assert low.p[10] < low.p[9]
    with pytest.raises(NativeExecutionError, match="strictly increasing n_B, P and epsilon"):
        construct_maxwell_transition(low, high, assume_local_charge_neutrality=True)
