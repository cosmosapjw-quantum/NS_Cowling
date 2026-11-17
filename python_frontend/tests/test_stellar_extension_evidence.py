"""E physical artifact regressions, independent of the generation runner.

The immutable D prefix, actual native support and causal/turning-point evidence
are the scientific contracts; these do not test runner implementation details.
"""
from pathlib import Path
import json
import numpy as np
import pytest
from compose_eos import read_compose
from compose_eos.crust_assembly import ROW_FIELDS
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_background import ContractEnthalpyEOS

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"data/stellar_extension_20261117"
OLD=ROOT/"data/core_support_20261117/primary_v1"

@pytest.fixture(scope="module")
def evidence():
    if not (OUT/"eos/EOS_RESULT.json").exists():
        pytest.skip("E generated evidence is not present")
    return json.loads((OUT/"eos/EOS_RESULT.json").read_text())

@pytest.mark.parametrize("grid",["coarse","medium","fine"])
def test_d_prefix_and_exact_thermodynamic_consumer_are_preserved(evidence,grid):
    before=read_compose(OLD/"compose")
    after=read_compose(OUT/"eos"/grid/"compose")
    assert after.n>before.n
    for name in ROW_FIELDS:
        np.testing.assert_array_equal(getattr(after,name)[:before.n],getattr(before,name))
    old_record=json.loads((OLD/"REPRESENTATION.json").read_text())
    new_record=json.loads((OUT/"eos"/grid/"REPRESENTATION.json").read_text())
    for key in ("surface","outer","lower_bridge"):
        assert new_record[key]==old_record[key]
    for key in ("n","p","epsilon"):
        assert new_record["upper"][key][:len(old_record["upper"][key])]==old_record["upper"][key]
    old=WholeCrustBarotrope.from_contract(old_record)
    new=WholeCrustBarotrope.from_contract(new_record)
    pressures=np.unique(np.r_[0.,before.p,np.geomspace(before.p[0],before.p[-1],83)])
    a,b=old.state_at_pressure(pressures),new.state_at_pressure(pressures)
    for key in ("n","p","epsilon","mu","enthalpy"):
        np.testing.assert_array_equal(a[key],b[key])
    np.testing.assert_array_equal(a["cs2_over_c2"][:-1],b["cs2_over_c2"][:-1])

def test_native_crossing_is_equilibrium_derivative_and_not_completed(evidence):
    endpoint=evidence["causal_endpoint"]
    root=endpoint["native_causal_root_fm3"]
    state=endpoint["native_state"]
    assert endpoint["native_bracket_fm3"][0]<root<endpoint["native_bracket_fm3"][1]
    assert endpoint["native_bracket_chi"][0]<1<endpoint["native_bracket_chi"][1]
    assert abs(state["cs2"]-1)<1e-10
    assert abs(state["cs2_fd"]-state["cs2"])<1e-8
    assert state["composition_hessian"]>0
    assert state["cs2_frozen"]>=state["cs2"]
    assert max(abs(r["root_delta_fm3"]) for r in endpoint["fd_convergence"])<1e-8
    assert 1<endpoint["admissible_export_endpoint_fm3"]<root

@pytest.mark.parametrize("grid",["coarse","medium","fine"])
def test_export_is_continuously_causal_and_rejects_outside_support(evidence,grid):
    path=OUT/"eos"/grid
    model=WholeCrustBarotrope.from_contract(json.loads((path/"REPRESENTATION.json").read_text()))
    eos=ContractEnthalpyEOS(model)
    assert model.density_domain[1]==evidence["causal_endpoint"]["admissible_export_endpoint_fm3"]
    assert all(i.validate()["status"]=="PASS" for i in model.upper_intervals)
    with pytest.raises(ValueError,match="outside"):
        model.state_at_pressure(np.nextafter(model.pressure_domain[1],np.inf))
    with pytest.raises(ValueError,match="outside"):
        eos.enthalpy_at_density(np.nextafter(model.density_domain[1],np.inf))
    data=read_compose(path/"generated")
    assert data.nb[0]==1.
    assert data.nb[-1]==model.density_domain[1]
    # Chemical potential and physical energy are tied by the unchanged cold
    # Euler identity; a compensating metadata-only energy shift cannot pass.
    np.testing.assert_allclose(data.mub,(data.eps+data.p)/data.nb,rtol=1e-12,atol=0)
