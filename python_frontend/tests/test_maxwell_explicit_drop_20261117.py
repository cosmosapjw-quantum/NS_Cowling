"""Explicit thermo-only Maxwell transport retains physical and support data."""
from __future__ import annotations
import hashlib
from pathlib import Path
import numpy as np
import pytest
from compose_eos import NativeExecutionError, construct_maxwell_transition, write_compose, read_compose
from compose_eos import api
from test_maxwell_mass_gate import physical_pair
from scripts.match_core_crust import run_pipeline

ABUNDANCES = ("Yn","Yp","Ye","Ymu","Aav","Zav","Yh")
MASKS = ("Yn_present","Yp_present","Ye_present","Ymu_present",
         "Abar_present","Zbar_present","Yh_present","C3_applicable","C4_applicable")
THERMO = ("nb","Q1","Q2","Q3","Q4","Q5","Q6","Q7","q_extra","Iphase")

def composed_pair():
    pair=physical_pair()
    for phase, label in zip(pair, (6,7)):
        phase.Iphase[:]=label
        phase.Yn[:]=.8; phase.Yp[:]=.2; phase.Ye[:]=.2
        phase.Aav[:]=1.
        phase.composition_available=True
        phase.composition_file_present=True
        for name in MASKS:
            setattr(phase,name,np.ones(phase.n,dtype=bool))
        phase.row_status[:]=2
    return pair

def test_default_maxwell_still_rejects_composition():
    low,high=composed_pair()
    with pytest.raises(NativeExecutionError,match="composition-bearing"):
        construct_maxwell_transition(low,high,assume_local_charge_neutrality=True)

def test_thermo_only_copy_preserves_all_thermo_and_source_without_aliases():
    low,_=composed_pair();low.source_dir="/source/immutable-phase"
    low.q_extra_count=1;low.q_extra[:]=.3
    view=api._thermo_only_view(low)
    for name in THERMO:
        np.testing.assert_array_equal(getattr(view,name),getattr(low,name))
        assert not np.shares_memory(getattr(view,name),getattr(low,name))
    for name in ("mn","mp","Il","temperature_mev","charge_fraction_grid",
                 "q_extra_count","beta_equilibrated","source_dir"):
        assert getattr(view,name)==getattr(low,name)
    assert view.composition_file_present is False
    assert view.composition_available is False
    for name in ABUNDANCES:
        assert not np.any(getattr(view,name))
    for name in MASKS:
        assert not np.any(getattr(view,name))
    assert not np.any(view.row_status)
    assert low.composition_available and low.composition_file_present
    assert np.all(low.Yn==.8) and np.all(low.Yn_present)

def test_explicit_drop_matches_native_thermo_only_output_and_keeps_inputs(tmp_path):
    low,high=composed_pair()
    before=[{name:getattr(d,name).copy()for name in THERMO+ABUNDANCES+MASKS} for d in (low,high)]
    rawlow,rawhigh=physical_pair()
    ref,_=construct_maxwell_transition(rawlow,rawhigh,assume_local_charge_neutrality=True,return_report=True)
    out,report=construct_maxwell_transition(low,high,assume_local_charge_neutrality=True,
                                          drop_composition=True,return_report=True)
    for name in THERMO:
        np.testing.assert_allclose(getattr(out,name),getattr(ref,name),rtol=0.,atol=2e-10)
    assert report["composition_dropped"] is True
    assert report["drop_composition_requested"] is True
    assert report["composition_preserved"] is False
    assert out.composition_file_present and not out.composition_available
    for d,copy in zip((low,high),before):
        for name,original in copy.items():
            np.testing.assert_array_equal(getattr(d,name),original)
    write_compose(tmp_path/"out",out)
    back=read_compose(tmp_path/"out")
    np.testing.assert_array_equal(back.Iphase,out.Iphase)
    assert not back.composition_available
    assert not np.any(back.Yn_present) and not np.any(back.Abar_present)

@pytest.mark.parametrize("bad", [1,"true",None])
def test_drop_composition_requires_explicit_boolean(bad):
    low,high=physical_pair()
    with pytest.raises(TypeError,match="drop_composition must be bool"):
        construct_maxwell_transition(low,high,assume_local_charge_neutrality=True,drop_composition=bad)

def test_optin_does_not_remove_unsupported_extra_quantities():
    low,high=composed_pair();low.q_extra_count=1;low.q_extra[:]=.3
    with pytest.raises(NativeExecutionError,match="Nadd"):
        construct_maxwell_transition(low,high,assume_local_charge_neutrality=True,drop_composition=True)

def test_pipeline_declared_thermo_only_maxwell_actually_drops_composition(tmp_path):
    low,high=composed_pair()
    lo,hi=tmp_path/"low",tmp_path/"high"
    write_compose(lo,low);write_compose(hi,high)
    original={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
              for d in (lo,hi) for p in d.iterdir()if p.is_file()}
    result=run_pipeline("maxwell",lo,hi,tmp_path/"matched")
    assert result["survived"],result
    assert result["composition_dropped"] is True
    assert result["composition_transport"]["policy"]=="explicit_thermo_only_pipeline"
    assert original=={name:hashlib.sha256(Path(name).read_bytes()).hexdigest()for name in original}
    for side,raw in (("crust",lo),("core",hi)):
        work=tmp_path/"matched"/"work"
        source=api.EOSData.from_internal_directory(work/(side+"_internal"))
        view=api.EOSData.from_internal_directory(work/(side+"_thermo_internal"))
        assert source.composition_file_present
        assert not view.composition_file_present and not view.composition_available
        for name in THERMO:
            np.testing.assert_array_equal(getattr(source,name),getattr(view,name))
