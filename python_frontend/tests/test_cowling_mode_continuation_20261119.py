"""Mode identity is a physical eigenfunction comparison, not root ordering."""
import numpy as np
import pytest

from compose_eos.stellar_mode_continuation import (
    eigenfunction_overlap, track_mode_family,
)


def mode(frequency,kind="inner",sign=1.,both=True):
    """Two-side synthetic profile with W and rV and distinct same-node shapes."""
    rows=[]
    for phase,r0,r1 in [(0,.001,10.),(1,10.,12.)]:
        for x in np.linspace(0,1,61):
            r=r0+(r1-r0)*x
            envelope=np.exp(-5*x) if kind=="inner" else np.exp(-5*(1-x))
            rows.append(dict(phase_index=phase,r_km=r,
                W=sign*envelope*(1+.08*phase),
                V_km_inverse=sign*envelope/(r+1e-4)))
    return dict(frequency_hz=frequency,profile=rows,
        W_nodes_per_phase=0,V_nodes_per_phase=0)


def test_overlap_is_sign_invariant_and_respects_two_sides_of_jump():
    original=mode(1900,kind="inner")
    same=mode(1950,kind="inner",sign=-1)
    dissimilar=mode(1902,kind="outer")
    x=eigenfunction_overlap(original,same)
    assert x["absolute_overlap"]>0.99999999999
    assert x["n_phases"]==2
    assert eigenfunction_overlap(original,dissimilar)["absolute_overlap"]<0.85


def test_frequency_nearest_can_be_wrong_when_eigenfunction_has_crossed():
    old=[mode(1800,"inner"),mode(2000,"outer")]
    current=[mode(1810,"outer"),mode(1980,"inner")]
    tracked=track_mode_family(old,current,min_overlap=.90,max_fractional_step=.25)
    assert tracked["status"]=="CONTINUED"
    assert tracked["pairs"]==[
         {"previous_index":0,"current_index":1,
          "absolute_overlap":pytest.approx(1.,abs=1e-12),
          "frequency_change_hz":pytest.approx(180.)},
         {"previous_index":1,"current_index":0,
          "absolute_overlap":pytest.approx(1.,abs=1e-12),
          "frequency_change_hz":pytest.approx(-190.)}]


def test_missing_continuity_returns_explicit_unresolved():
    old=[mode(1800,"inner")]
    next_cases=[mode(1802,"outer")]
    result=track_mode_family(old,next_cases,min_overlap=.90)
    assert result["status"]=="UNRESOLVED"
    assert result["pairs"]==[]


def test_missing_phase_or_wrong_root_node_identity_rejected():
    original=mode(1800,"inner")
    modified=mode(1820,"inner")
    modified["profile"]=[r for r in modified["profile"] if r["phase_index"]==0]
    with pytest.raises(ValueError,match="phase"):
        eigenfunction_overlap(original,modified)
    correct=mode(1820,"inner")
    correct["W_nodes_per_phase"]=1
    z=track_mode_family([original],[correct])
    assert z["status"]=="UNRESOLVED"


def test_scientific_profiles_must_have_finite_increasing_radius_in_each_phase():
    a=mode(1800)
    b=mode(1850)
    b["profile"][15]["r_km"]=float("nan")
    with pytest.raises(ValueError):
        eigenfunction_overlap(a,b)


def test_v_node_transition_can_be_conditionally_followed_by_shape():
    prev=mode(1800,"inner")
    curr=mode(1830,"inner")
    curr["V_nodes_per_phase"]=1
    strict=track_mode_family([prev],[curr])
    assert strict["status"]=="UNRESOLVED"
    continued=track_mode_family([prev],[curr],
        v_node_policy="shape_continuity")
    assert continued["status"]=="CONTINUED_WITH_V_NODE_CHANGE"
    assert continued["pairs"][0]["absolute_overlap"]>0.99999
    assert continued["pairs"][0]["V_node_transition"]==[0,1]
    assert "conditional" in continued["policy"]


def test_conditional_node_change_is_not_a_loose_gate():
    prev=mode(1800,"inner")
    curr=mode(1850,"inner")
    curr["V_nodes_per_phase"]=2
    assert track_mode_family([prev],[curr],
                 v_node_policy="shape_continuity")["status"]=="UNRESOLVED"
    curr=mode(1850,"outer")
    curr["V_nodes_per_phase"]=1
    assert track_mode_family([prev],[curr],
                 v_node_policy="shape_continuity")["status"]=="UNRESOLVED"
