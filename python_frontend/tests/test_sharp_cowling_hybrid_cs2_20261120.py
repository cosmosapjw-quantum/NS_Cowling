"""Source-sound-speed / smoothly interpolated metric backend contract.

The high-accuracy physical source derivative c_s^2 is never approximated
by PCHIP across a discontinuity. Only (r,Phi,Lambda,Phi') is interpolated
within each physical phase. Never claim original microphysical EOS identity.
"""
from __future__ import annotations
import io,json,zipfile
from pathlib import Path
import numpy as np
import pytest
from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (
    NativeSharpCowlingDVODE, solve_sharp_background_native)
from compose_eos.whole_crust import WholeCrustBarotrope

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"


def load_eos(model,join):
    with zipfile.ZipFile(SOURCE) as a:
        with zipfile.ZipFile(io.BytesIO(a.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as b:
            key=f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json"
            contract=json.loads(b.read(key))
    return SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(contract))


def test_hybrid_accepts_new_backend_and_no_cross_phase_source_data():
    eos=load_eos("SPG_M3","css")
    bg=solve_sharp_background_native(eos,.192510767556007,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,
        max_step_h=.004,native_samples=1025)
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,
        atol=1e-12,max_step_h=.002)
    mode=NativeSharpCowlingDVODE(bg,settings,samples=1025,
        coefficient_backend="hybrid_cs2",
        source_restart_relative_jump=.01)
    assert mode.coefficient_backend=="hybrid_cs2"
    assert len(mode._phases)==2
    assert mode.source_restart_count>=30
    for phase in mode._phases:
        xs=np.linspace(phase["s"][0],phase["s"][-1],29)[1:-1]
        for sv in xs:
            vals=mode._source_coefficients(phase,float(sv))
            h=bg.h_start-float(sv)**2
            original=bg.state_at_enthalpy(h,eos.side_for_segment(h,phase["mid"]))
            assert vals[4]==original["cs2_over_c2"]  # exact physical EOS
            true=np.array([original[k] for k in
               ("r_km","Phi","lambda_metric","dPhi_dr_km")],float)
            relative=abs((np.array(vals[:4])-true)/
                       np.maximum(abs(true),[.05,1.,1.,1e-5]))
            assert max(relative)<.003


def test_hybrid_uses_different_coefficients_than_near_knot_pchip_smearing():
    eos=load_eos("SPG_M3","css")
    bg=solve_sharp_background_native(eos,.192510767556007,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,
        max_step_h=.004,native_samples=1025)
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,
        atol=1e-12,max_step_h=.002)
    hybrid=NativeSharpCowlingDVODE(bg,settings,samples=1025,
        coefficient_backend="hybrid_cs2")
    pchip=NativeSharpCowlingDVODE(bg,settings,samples=1025,
        coefficient_backend="pchip")
    errors=[]
    for phase,old in zip(hybrid._phases,pchip._phases):
        for s in phase["source_restart_s"]:
            # A strongly discontinuous c_s² may be smeared only inside
            # the deliberately narrow ~1e-10 s buffer around each knot.
            # Probe the one-sided limits at half that width.
            for dv in (-5e-11,5e-11):
                sv=float(np.clip(s+dv,phase["s"][0],phase["s"][-1]))
                if sv in (phase["s"][0],phase["s"][-1]):
                    continue
                physical=hybrid._source_coefficients(phase,sv)[4]
                interpolated=float(old["pchip"](sv)[4])
                errors.append(abs(physical-interpolated)/
                              max(physical,1e-6))
    assert max(errors)>.01


@pytest.mark.slow
def test_actual_M2_Maxwell_1p4_rapid_mode_hybrid_matches_direct_EOS_and_nodes():
    eos=load_eos("SPG_M2","maxwell")
    bg=solve_sharp_background_native(eos,.30385128926615257,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,
        max_step_h=.004,native_samples=2049)
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,
        atol=1e-12,max_step_h=.002)
    hybrid=NativeSharpCowlingDVODE(bg,settings,samples=4097,
        coefficient_backend="hybrid_cs2")
    mode=hybrid.solve_mode((5110.,5170.),conversion="rapid",
        profile_points=350)
    assert abs(mode["frequency_hz"]-5143.719721563306)<.05
    assert (mode["W_nodes_per_phase"],mode["V_nodes_per_phase"])==(1,1)
    assert mode["coefficient_backend"]=="hybrid_cs2"
    assert mode["relative_join_mismatch"]<1e-5


@pytest.mark.slow
def test_actual_M3_CSS_0p7_slow_hybrid_recovers_strict_p1():
    eos=load_eos("SPG_M3","css")
    bg=solve_sharp_background_native(eos,.192510767556007,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,
        max_step_h=.004,native_samples=2049)
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,
        atol=1e-12,max_step_h=.002)
    hybrid=NativeSharpCowlingDVODE(bg,settings,samples=4097,
        coefficient_backend="hybrid_cs2")
    mode=hybrid.solve_mode((3710.,3760.),conversion="slow",
        profile_points=350)
    assert abs(mode["frequency_hz"]-3736.2088951581254)<.05
    assert (mode["W_nodes_per_phase"],mode["V_nodes_per_phase"])==(1,1)
    assert mode["coefficient_backend"]=="hybrid_cs2"
    assert mode["relative_join_mismatch"]<1e-5


def test_hybrid_callback_failure_does_not_poison_default_native_mode(monkeypatch):
    eos=load_eos("SPG_M3","css")
    bg=solve_sharp_background_native(eos,.192510767556007,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,
        max_step_h=.004,native_samples=1025)
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,
        atol=1e-12,max_step_h=.002)
    hybrid=NativeSharpCowlingDVODE(bg,settings,samples=129,
        coefficient_backend="hybrid_cs2")
    original=eos.state_at_enthalpy
    def corrupted_source(h,side="high"):
        raise ValueError("INJECTED_HYBRID_SOURCE_FAILURE")
    monkeypatch.setattr(eos,"state_at_enthalpy",corrupted_source)
    with pytest.raises(RuntimeError,match="INJECTED_HYBRID_SOURCE_FAILURE"):
        hybrid.matching_residual(3736.2,conversion="slow")
    monkeypatch.setattr(eos,"state_at_enthalpy",original)
    # The F90 C callback pointer must be reset even on a physics source
    # exception: a new default PCHIP mode must not reuse stale callbacks.
    default=NativeSharpCowlingDVODE(bg,settings,samples=1025)
    assert default.coefficient_backend=="pchip"
    assert np.isfinite(default.matching_residual(3736.2,conversion="slow"))
    assert default._callback_rhs_calls==0
