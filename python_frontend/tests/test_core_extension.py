"""Behavioral tests for preserving A while appending an authenticated core tail."""
from copy import deepcopy
from dataclasses import replace
import importlib

import numpy as np
import pytest

from compose_eos.crust_assembly import ROW_FIELDS, select_rows
from compose_eos.eos_data import EOSData
from compose_eos.whole_crust import (
    HermiteInterval, PiecewisePolytrope, SurfacePolytrope, WholeCrustBarotrope,
)


def table(n, *, mass=939., phase=6):
    # Independent cold EOS: epsilon=4*n+n^2, P=n^2, mu=4+2*n.
    n = np.asarray(n, dtype=float)
    p, eps = n*n, 4*n+n*n
    z = np.zeros_like(n)
    return EOSData(nb=n, mn=mass, mp=mass, Q1=p/n, Q2=z.copy(),
        Q3=(4+2*n)/mass-1, Q4=z.copy(), Q5=z.copy(),
        Q6=eps/(n*mass)-1, Q7=eps/(n*mass)-1, q_extra=z.copy(),
        Iphase=np.full(n.size, phase), Yn=np.full(n.size,.9),
        Yp=np.full(n.size,.1), Ye=np.full(n.size,.1), Ymu=z.copy(),
        Aav=np.ones(n.size), Zav=z.copy(), Yh=z.copy(),
        composition_available=True, composition_file_present=True,
        beta_equilibrated=True)


def inputs():
    whole = table([.001,.01,.05,.1,.2,.3,.5], mass=1., phase=4)
    surface = SurfacePolytrope(.01,.0001,.0401,2.)
    outer = PiecewisePolytrope([.01,.05,.1],[.0001,.0025,.01],4.01)
    bridge = HermiteInterval(.1,.01,.41,.2,.04,.84)
    contract = dict(surface=surface.contract(), outer=outer.contract(),
        lower_bridge=bridge.contract(),
        upper=dict(n=[.2,.3,.5],p=[.04,.09,.25],
                   epsilon=[.84,1.29,2.25],sharp_interfaces=[]))
    core = table([.25,.6,.8,1.],mass=1.)
    anchor = table([.5],mass=1.)
    return whole, contract, core, anchor


def extension(*args):
    try:
        fn = importlib.import_module("compose_eos.core_extension").extend_whole_crust
    except ModuleNotFoundError:
        pytest.fail("core-tail extension is not implemented")
    return fn(*args)


def test_append_preserves_all_old_rows_lower_contract_and_original_inputs():
    whole, contract, core, anchor = inputs()
    before = deepcopy(contract)
    combined, record, report = extension(whole, contract, core, anchor)
    for field in ROW_FIELDS:
        np.testing.assert_array_equal(getattr(combined,field)[:whole.n],getattr(whole,field))
    for key in ("surface","outer","lower_bridge"):
        assert record[key] == before[key]
    assert contract == before
    assert combined.n == 10
    assert report["appended_rows"] == 3
    assert record["upper"]["n"][:3] == [.2,.3,.5]
    assert combined.Iphase[-1] == 4
    state = WholeCrustBarotrope.from_contract(record).state_at_pressure(.64)
    assert state["n"] == pytest.approx(.8,rel=2e-12)
    assert state["epsilon"] == pytest.approx(3.84,rel=2e-12)
    assert state["cs2_over_c2"] == pytest.approx(2/7,rel=2e-11)


def test_header_change_preserves_dimensional_core_energy():
    whole, contract, core, anchor = inputs()
    core,anchor = table([.25,.6,.8,1.],mass=1.1),table([.5],mass=1.1)
    combined, record, report = extension(whole,contract,core,anchor)
    np.testing.assert_allclose(combined.eps[-3:],[2.76,3.84,5.],rtol=2e-14)
    assert combined.mn == whole.mn
    assert report["energy_reference_shift"] == "NONE"


def test_rejects_an_energy_shift_at_existing_anchor():
    whole, contract, core, anchor = inputs()
    anchor = replace(anchor,Q7=anchor.Q7+.01/anchor.mn)
    with pytest.raises(ValueError,match="anchor"):
        extension(whole,contract,core,anchor)


def test_rejects_a_hidden_native_phase_change():
    whole, contract, core, anchor = inputs()
    phase = core.Iphase.copy(); phase[-1]=7
    with pytest.raises(ValueError,match="phase"):
        extension(whole,contract,replace(core,Iphase=phase),anchor)


def test_rejects_a_noncausal_tail_without_modifying_input():
    whole, contract, core, anchor = inputs()
    # Cold first-law EOS with asymptotic cs2=2: epsilon=n+10*n^3, P=20*n^3.
    n=np.asarray([.6,.8,1.])
    eps=n+10*n**3; p=20*n**3
    bad=table(n,mass=1.)
    bad=replace(bad,Q1=p/n,Q3=(eps+p)/n-1,Q6=eps/n-1,Q7=eps/n-1)
    saved=deepcopy(contract)
    with pytest.raises(ValueError,match="causal|stable|endpoints"):
        extension(whole,contract,bad,anchor)
    assert contract == saved


def test_refuses_metadata_only_support_with_no_high_density_rows():
    whole, contract, core, anchor = inputs()
    with pytest.raises(ValueError,match="extend|support"):
        extension(whole,contract,table([.1,.2,.4],mass=1.),anchor)


def test_consumer_refuses_pressure_beyond_actual_extended_support():
    whole, contract, core, anchor = inputs()
    _,record,_ = extension(whole,contract,core,anchor)
    model=WholeCrustBarotrope.from_contract(record)
    with pytest.raises(ValueError,match="domain|support"):
        model.state_at_pressure(np.nextafter(1.,np.inf))
