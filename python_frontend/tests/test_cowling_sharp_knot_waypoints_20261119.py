"""Exact EOS break-point waypoints must never cross a Maxwell interface."""
import json
from pathlib import Path
import numpy as np
import pytest

from compose_eos.stellar_sharp_native import (
    native_shoot_waypoints,NativeSharpCowlingDVODE,solve_sharp_background_native,
)
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_cowling import CowlingSettings


def test_ascending_waypoints_map_original_requested_targets():
    req=np.array([0,.5,1.],float)
    path,indices=native_shoot_waypoints(req,np.array([-.1,.2,.8,1.0,1.1]))
    np.testing.assert_allclose(path,[0.,.2,.5,.8,1.])
    np.testing.assert_allclose(path[indices],req)
    assert indices.tolist()==[0,2,4]


def test_descending_waypoints_retain_physical_order():
    req=np.array([1.,.7,0.],float)
    path,indices=native_shoot_waypoints(req,[.2,.8,.7,1.2])
    np.testing.assert_allclose(path,[1.,.8,.7,.2,0.])
    np.testing.assert_allclose(path[indices],req)
    assert indices.tolist()==[0,2,4]


def test_waypoints_do_not_silently_extrapolate_or_accept_duplicate_requested():
    with pytest.raises(ValueError):
        native_shoot_waypoints([0.,0.,1.],[.2])
    with pytest.raises(ValueError):
        native_shoot_waypoints([0.,1.],[float("nan")])


def test_native_maxwell_knot_aware_shooting_has_distinct_report():
    root=Path(__file__).resolve().parents[2]
    source=root/"data/cowling_maxwell_native_20261119/three_layer_original_fsu2/maxwell/REPRESENTATION.json"
    eos=SharpContractEnthalpyEOS(
        WholeCrustBarotrope.from_contract(json.loads(source.read_text())))
    bg=solve_sharp_background_native(eos,.30,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,max_step_h=.004,
        native_samples=2049)
    settings=CowlingSettings(surface_h=3e-8,rtol=3e-9,atol=3e-12,max_step_h=.001)
    sol=NativeSharpCowlingDVODE(bg,settings,samples=2049,knot_waypoints=True)
    assert sol.knot_waypoints
    assert sum(len(ph["source_break_s"]) for ph in sol._phases)>50
    freq=sol.solve_mode((1780,1840),conversion="slow",profile_points=120)
    assert freq["knot_waypoints_used"] is True
    assert abs(freq["frequency_hz"]-1809.53)<1.
    assert freq["W_nodes_per_phase"]==0
    assert freq["V_nodes_per_phase"]==0
