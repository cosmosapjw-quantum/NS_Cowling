"""Certified native/DOP853 sharp Cowling: red-first regression.

The PDF mass-sweep's SPG_M3 CSS 0.7 Msun p1 is the minimum real
counterexample: native tabulation can fail while direct-EOS DOP853
computes a regular eigenfrequency. A failed solver is never relabelled PASS.
"""
import io,json,zipfile
from pathlib import Path
import pytest
import numpy as np

from compose_eos.stellar_sharp_certified import (
    classify_cowling_nodes, select_native_sample_convergence,
    count_reference_nodes_by_phase, solve_certified_sharp_mode
)
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_cowling import CowlingSettings

ROOT=Path(__file__).resolve().parents[2]
ZIP=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"


def test_strict_f_p1_not_confused_with_W1_V2():
    assert classify_cowling_nodes(0,0)=="STRICT_F"
    assert classify_cowling_nodes(1,1)=="STRICT_P1"
    assert classify_cowling_nodes(1,2)=="W1_V2_CLASSIFICATION_UNRESOLVED"
    assert classify_cowling_nodes(2,1)=="UNCLASSIFIED"
    with pytest.raises(ValueError):
        classify_cowling_nodes(-1,1)


def test_frequency_sample_convergence_is_independent_of_matching_residual():
    m1=dict(samples=4097,status="ADMITTED",frequency_hz=5329.201,
            W_nodes_per_phase=1,V_nodes_per_phase=2,
            relative_join_mismatch=2e-5,profile=[])
    m2=dict(samples=8193,status="ADMITTED",frequency_hz=5329.202,
            W_nodes_per_phase=1,V_nodes_per_phase=2,
            relative_join_mismatch=2e-5,profile=[])
    z=select_native_sample_convergence([m1,m2],
            absolute_tolerance_hz=.02,max_join=1e-5)
    assert z["status"]=="UNRESOLVED"
    assert z["reason"]=="INSUFFICIENT_ADMITTED_MODES"
    assert z["native_converged"] is False


def test_native_convergence_requires_both_node_identity_and_profile_overlap():
    rows=[
        dict(samples=4097,status="ADMITTED",frequency_hz=5200.,
             W_nodes_per_phase=1,V_nodes_per_phase=1,
             relative_join_mismatch=1e-7,profile=[]),
        dict(samples=8193,status="ADMITTED",frequency_hz=5200.001,
             W_nodes_per_phase=1,V_nodes_per_phase=2,
             relative_join_mismatch=1e-7,profile=[])]
    result=select_native_sample_convergence(rows,
             absolute_tolerance_hz=.02,max_join=1e-5)
    assert result["status"]=="UNRESOLVED"
    assert result["reason"]=="NODE_PAIR_CHANGED"
    with pytest.raises(ValueError):
        select_native_sample_convergence(rows,absolute_tolerance_hz=-.1)


def load_derived_eos(model,method):
    """Extract exactly source-pinned created EOS; no original raw sources."""
    with zipfile.ZipFile(ZIP) as outer:
        inner=io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))
        with zipfile.ZipFile(inner) as data:
            name=f"FRESH_5x2_WHOLE_STARS/{model}/{method}/REPRESENTATION.json"
            return SharpContractEnthalpyEOS(
                WholeCrustBarotrope.from_contract(json.loads(data.read(name))))


def test_reference_node_count_does_not_count_a_jump_as_radial_zero():
    transitions=[dict(h=.03),dict(h=.01)]
    profile=[
        dict(h=.035,side="high",W=1.,V_km_inverse=1.),
        dict(h=.03,side="high",W=.6,V_km_inverse=.4),
        dict(h=.03,side="low",W=-.8,V_km_inverse=-.6),
        dict(h=.02,side="low",W=-.5,V_km_inverse=-.2),
        dict(h=.01,side="high",W=-.2,V_km_inverse=-.2),
        dict(h=.01,side="low",W=.8,V_km_inverse=.7),
        dict(h=.005,side="low",W=.9,V_km_inverse=.9)]
    x=count_reference_nodes_by_phase(profile,transitions)
    assert x["W_nodes_per_phase"]==0
    assert x["V_nodes_per_phase"]==0
    assert x["n_phases"]==3


@pytest.mark.slow
def test_real_M3_CSS_0p7_native_is_rejected_but_direct_eos_oracle_recovers_mode():
    eos=load_derived_eos("SPG_M3","css")
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,
            atol=1e-12,max_step_h=.002)
    r=solve_certified_sharp_mode(
        eos,.192510767556007,(3710.,3760.),"slow",expected="p1",
        settings=settings,native_samples=(4097,8193),
        absolute_tolerance_hz=.05,profile_points=500)
    assert abs(r["frequency_hz"]-3736.20894)<.10
    assert r["status"]=="REFERENCE_ONLY_NATIVE_NOT_CERTIFIED"
    assert r["selected_backend"]=="PYTHON_DOP853_DIRECT_EOS"
    assert r["node_gate"]=="STRICT_P1"
    assert r["python_reference"]["relative_join_mismatch"]<1e-6
    assert r["python_reference"]["mass_deviation_from_native_Msun"]<3e-6
    assert len(r["native_attempts"])==2
    assert not r["native_convergence"]["native_converged"]
    assert r["source_eos_replaced"] is False


@pytest.mark.slow
def test_real_M3_CSS_1p4_W1_V2_not_unconditionally_called_p1():
    eos=load_derived_eos("SPG_M3","css")
    r=solve_certified_sharp_mode(
        eos,.304225985123765,(5290.,5360.),"slow",expected="p1",
        native_samples=(4097,8193),absolute_tolerance_hz=.05,
        profile_points=500)
    assert abs(r["frequency_hz"]-5329.20271)<.08
    assert r["node_gate"]=="W1_V2_CLASSIFICATION_UNRESOLVED"
    assert r["expected_mode_identification"]=="CONDITIONAL_NUMERIC_FAMILY_NOT_STRICT_P1"
    assert r["status"]=="REFERENCE_ONLY_NATIVE_NOT_CERTIFIED"


def test_native_sample_gate_does_not_skip_failed_intermediate_resolution():
    def profile():
        values=[]
        for phase in (0,1):
            for k in range(12):
                r=.2+phase*10+k*.1
                values.append(dict(phase_index=phase,r_km=r,
                    W=1.0/(1.+r),V_km_inverse=2./(1.+r)))
        return values
    a=dict(samples=4097,status="ADMITTED",frequency_hz=5148.0,
           W_nodes_per_phase=1,V_nodes_per_phase=1,
           relative_join_mismatch=1e-7,profile=profile())
    failed=dict(samples=8193,status="REJECTED",reason="DVODE interface mismatch")
    b=dict(samples=16385,status="ADMITTED",frequency_hz=5148.001,
           W_nodes_per_phase=1,V_nodes_per_phase=1,
           relative_join_mismatch=1e-7,profile=profile())
    result=select_native_sample_convergence([a,failed,b],
             absolute_tolerance_hz=.02,max_join=1e-5)
    assert result["native_converged"] is False
    assert result["reason"]=="NONCONSECUTIVE_RESOLUTION_VALIDATION"
