"""Test source-evaluating Fortran Cowling against pinned direct-DOP853 oracle.

A tiny DVODE eigenvalue residual for interpolated coefficients does not imply
the physical EOS has been integrated accurately. The new optional source
callback may not smear internal cs2 jumps via global PCHIP.
"""
import io
import json
import zipfile
from pathlib import Path
import pytest

from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (
    solve_sharp_background_native,NativeSharpCowlingDVODE,
)
from compose_eos.whole_crust import WholeCrustBarotrope

ROOT=Path(__file__).resolve().parents[2]
ZIP=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"


def source_eos(model="SPG_M3",method="css"):
    with zipfile.ZipFile(ZIP) as archive:
        with zipfile.ZipFile(io.BytesIO(archive.read(
               "05_SPG5_PK1_10_WHOLE_STARS.zip"))) as sub:
            payload=sub.read(
                f"FRESH_5x2_WHOLE_STARS/{model}/{method}/REPRESENTATION.json")
    return SharpContractEnthalpyEOS(
        WholeCrustBarotrope.from_contract(json.loads(payload)))


def test_mode_source_callback_argument_is_opt_in_and_rejects_bad_input():
    eos=source_eos()
    bg=solve_sharp_background_native(eos,.304225985123765,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,
        max_step_h=.004,native_samples=1025)
    setting=CowlingSettings(surface_h=3e-8,rtol=1e-9,
        atol=1e-12,max_step_h=.002)
    with pytest.raises(ValueError,match="coefficient_backend"):
        NativeSharpCowlingDVODE(bg,setting,samples=129,
            coefficient_backend="invalid")
    fresh=NativeSharpCowlingDVODE(bg,setting,samples=129,
        coefficient_backend="direct_source")
    assert fresh.coefficient_backend=="direct_source"
    assert all(len(x["s"])<=80 for x in fresh._phases)
    assert fresh.source_break_count>=0


@pytest.mark.slow
def test_real_m3_css_0p7_direct_source_dvode_recovers_python_oracle():
    eos=source_eos()
    bg=solve_sharp_background_native(eos,.192510767556007,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,
        max_step_h=.004,native_samples=2049)
    setting=CowlingSettings(surface_h=3e-8,rtol=1e-9,
        atol=1e-12,max_step_h=.002)
    solver=NativeSharpCowlingDVODE(bg,setting,samples=129,
          coefficient_backend="direct_source")
    root=solver.solve_mode((3710.,3760.),conversion="slow",profile_points=350)
    assert root["backend"].startswith("Fortran DVODE sharp")
    assert root["coefficient_backend"]=="direct_source"
    assert abs(root["frequency_hz"]-3736.208895)<.05
    assert root["W_nodes_per_phase"]==root["V_nodes_per_phase"]==1
    assert root["relative_join_mismatch"]<1e-5
    assert root["source_callback_rhs_calls"]>10


@pytest.mark.slow
def test_real_m3_css_1p4_direct_source_retains_W1V2_identity():
    eos=source_eos()
    bg=solve_sharp_background_native(eos,.304225985123765,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,
        max_step_h=.004,native_samples=2049)
    setting=CowlingSettings(surface_h=3e-8,rtol=1e-9,
        atol=1e-12,max_step_h=.002)
    solver=NativeSharpCowlingDVODE(bg,setting,samples=129,
          coefficient_backend="direct_source")
    root=solver.solve_mode((5290.,5360.),conversion="slow",profile_points=350)
    assert abs(root["frequency_hz"]-5329.202674)<.05
    assert (root["W_nodes_per_phase"],root["V_nodes_per_phase"])==(1,2)
    assert root["relative_join_mismatch"]<1e-5


@pytest.mark.slow
def test_native_exact_callback_failure_does_not_poison_next_pchip_shot(monkeypatch):
    eos=source_eos()
    bg=solve_sharp_background_native(eos,.304225985123765,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,
        max_step_h=.004,native_samples=2049)
    setting=CowlingSettings(surface_h=3e-8,rtol=1e-9,
        atol=1e-12,max_step_h=.002)
    exact=NativeSharpCowlingDVODE(bg,setting,samples=129,
                                  coefficient_backend="direct_source")
    original=bg.state_at_enthalpy
    def inject_failure(h,side="high"):
        if .04 < h < .08:
            raise ValueError("INJECTED_SOURCE_COEFFICIENT_FAILURE")
        return original(h,side)
    monkeypatch.setattr(bg,"state_at_enthalpy",inject_failure)
    with pytest.raises(RuntimeError,match="INJECTED_SOURCE_COEFFICIENT_FAILURE"):
        exact.matching_residual(5329.2,conversion="slow")
    monkeypatch.setattr(bg,"state_at_enthalpy",original)
    # Default native calls must have the old PCHIP table, never a stale C
    # function pointer from the previous rejected callback shot.
    old=NativeSharpCowlingDVODE(bg,setting,samples=4097)
    old_mode=old.solve_mode((5290.,5360.),conversion="slow",profile_points=350)
    assert old_mode["coefficient_backend"]=="pchip"
    assert abs(old_mode["frequency_hz"]-5329.202674)>.20
    assert old_mode["source_callback_rhs_calls"]==0


@pytest.mark.slow
def test_native_direct_source_restarts_at_resolved_sound_speed_jumps():
    eos=source_eos()
    bg=solve_sharp_background_native(eos,.192510767556007,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,
        max_step_h=.004,native_samples=2049)
    settings=[CowlingSettings(surface_h=3e-8,rtol=r,
                 atol=r*1e-3,max_step_h=step)
              for r,step in ((1e-8,.002),(1e-9,.001),(1e-10,.0005))]
    frequencies=[]
    for setting in settings:
        native=NativeSharpCowlingDVODE(bg,setting,samples=129,
              coefficient_backend="direct_source",
              source_restart_relative_jump=.01)
        assert native.source_restart_count>=30
        mode=native.solve_mode((3710.,3760.),conversion="slow",
                               profile_points=350)
        assert (mode["W_nodes_per_phase"],mode["V_nodes_per_phase"])==(1,1)
        assert mode["relative_join_mismatch"]<1e-5
        assert abs(mode["frequency_hz"]-3736.208895)<.05
        assert mode["source_restart_count"]>=30
        frequencies.append(mode["frequency_hz"])
    assert max(frequencies)-min(frequencies)<.05


def test_source_restart_threshold_validates_non_numeric_and_limits():
    eos=source_eos()
    bg=solve_sharp_background_native(eos,.304225985123765,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,max_step_h=.004,
        native_samples=1025)
    setting=CowlingSettings(surface_h=3e-8,rtol=1e-9,
                            atol=1e-12,max_step_h=.002)
    for invalid in (0.,-0.1,1.,1.1,float("nan"),"0.01",True):
        with pytest.raises(ValueError,match="source_restart_relative_jump"):
            NativeSharpCowlingDVODE(bg,setting,samples=129,
                  coefficient_backend="direct_source",
                  source_restart_relative_jump=invalid)
    # None opts out of EOS knot restarts while continuing to evaluate
    # source coefficients exactly at every RHS callback.
    optout=NativeSharpCowlingDVODE(bg,setting,samples=129,
        coefficient_backend="direct_source",source_restart_relative_jump=None)
    assert optout.source_restart_count==0
