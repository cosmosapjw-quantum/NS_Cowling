"""Cross-method physical identities on immutable H contracts."""
import importlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
PLAN = ROOT / "data/low_density_crust_20261117/stellar/causal_extension/TOV_PLAN.json"


def _source():
    assert importlib.util.find_spec("matching_interpolation_shapes_20261118") is not None, (
        "Expected new matching-interpolation crosscheck implementation")
    return importlib.import_module("matching_interpolation_shapes_20261118")


def _case_rows():
    return json.loads(PLAN.read_text())["cases"]


@pytest.mark.parametrize("entry", _case_rows(), ids=lambda c: c["id"])
@pytest.mark.parametrize("shape", ["power", "rational"])
def test_endpoint_first_law_and_causality(entry, shape):
    from compose_eos.whole_crust import WholeCrustBarotrope
    m = _source()
    contract = json.loads((ROOT / entry["contract"]).read_text())
    low = contract["lower_bridge"]["endpoint_states"]
    vals = [low[0][k] for k in ("n", "p", "epsilon")]
    vals += [low[1][k] for k in ("n", "p", "epsilon")]
    bridge = m.CausalShapeInterval(*vals, family=shape)
    assert bridge.validate()["status"] == "PASS"
    assert np.isclose(bridge.pressure(0), vals[1], atol=1e-14, rtol=0)
    assert np.isclose(bridge.pressure(1), vals[4], atol=1e-12, rtol=1e-12)
    assert np.isclose(bridge.number(0), vals[0], atol=1e-14, rtol=0)
    assert np.isclose(bridge.number(1), vals[3], atol=1e-12, rtol=1e-12)
    x = np.linspace(0, 1, 401)
    st = bridge.evaluate_mu(bridge.mu_left + bridge.width*x)
    assert np.min(st["cs2_over_c2"]) > 0
    assert np.max(st["cs2_over_c2"]) <= 1 + 4e-12
    assert np.min(np.diff(st["n"])) > 0
    assert np.max(np.abs(st["mu"]*st["n"]-st["p"]-st["epsilon"])) < 5e-11
    for t in (0.1, 0.3, 0.65, 0.9):
        dx = 1e-5
        dP = (bridge.pressure(t+dx)-bridge.pressure(t-dx))/(2*dx*bridge.width)
        dn = (bridge.number(t+dx)-bridge.number(t-dx))/(2*dx*bridge.width)
        assert abs(dP/bridge.number(t)-1) < 2e-6
        assert abs(dn/bridge.susceptibility(t)-1) < 2e-6
    b = WholeCrustBarotrope.from_contract(contract)
    old_sharp = set(b.sharp_interfaces)
    b.bridge = bridge
    assert old_sharp == b.sharp_interfaces
    assert b.upper_n[0] == vals[3]
    assert len(b.upper_intervals) == len(b.upper_n)-1


def test_negative_control_polytrope_preserves_pressure_but_not_energy():
    m = _source()
    entry = _case_rows()[0]
    low = json.loads((ROOT / entry["contract"]).read_text())["lower_bridge"]["endpoint_states"]
    v = [low[0][k] for k in ("n", "p", "epsilon")]
    v += [low[1][k] for k in ("n", "p", "epsilon")]
    c = m.naive_loglog_polytrope_endpoint(*v)
    assert abs(c["pressure_endpoint_residual"]) < 1e-12
    assert abs(c["energy_endpoint_residual_MeV_fm3"]) > 1e-9

@pytest.mark.parametrize("shape", ["power", "rational"])
def test_scalar_enthalpy_evaluator_contract(shape):
    m = _source()
    entry = _case_rows()[0]
    row = json.loads((ROOT / entry["contract"]).read_text())["lower_bridge"]["endpoint_states"]
    v = [row[0][k] for k in ("n", "p", "epsilon")]
    v += [row[1][k] for k in ("n", "p", "epsilon")]
    bridge = m.CausalShapeInterval(*v, family=shape)
    state = bridge.evaluate_mu(bridge.mu_left + 0.5*bridge.width)
    assert all(np.isfinite(vv) and np.isscalar(vv) for vv in state.values())

def test_numpy_validation_diagnostics_roundtrip_as_json():
    module = importlib.import_module("run_matching_interpolation_crosscheck_20261118")
    assert hasattr(module, "native_json_value"), "numpy diagnostic JSON adapter required"
    payload = dict(ok=np.bool_(True), cs2=np.float64(0.13), points=np.asarray([1, 2]))
    loaded = json.loads(json.dumps(payload, default=module.native_json_value))
    assert loaded == dict(ok=True, cs2=0.13, points=[1, 2])

@pytest.mark.parametrize("entry", _case_rows(), ids=lambda c: c["id"])
def test_bernstein_matches_neighbor_sound_speeds_and_causal_area(entry):
    from compose_eos.whole_crust import WholeCrustBarotrope
    m = _source()
    assert hasattr(m, "SlopeMatchedBernsteinInterval")
    contract=json.loads((ROOT/entry["contract"]).read_text())
    bar=WholeCrustBarotrope.from_contract(contract)
    end=bar.bridge
    cs_left=float(bar.outer.gamma[-1]*end.p_left/(end.epsilon_left+end.p_left))
    upper=bar.upper_intervals[0]
    cs_right=float(upper.number(0)/(upper.mu_left*upper.susceptibility(0)))
    v=[contract["lower_bridge"]["endpoint_states"][i][k] for i in (0,1)
       for k in ("n","p","epsilon")]
    b=m.SlopeMatchedBernsteinInterval(*v,left_cs2=cs_left,right_cs2=cs_right)
    assert b.validate()["status"]=="PASS"
    assert np.isclose(b.evaluate_mu(b.mu_left)["cs2_over_c2"],cs_left,rtol=1e-8)
    assert np.isclose(b.evaluate_mu(b.mu_right)["cs2_over_c2"],cs_right,rtol=1e-8)
    assert np.isclose(b.pressure(1),end.p_right,atol=1e-12)
    assert np.isclose(b.number(1),end.n_right,atol=1e-12)
    assert b.degree in (4,6,8,12,16,24,32)
    assert min(b.bernstein_coefficients)>=-1e-10

@pytest.mark.parametrize("entry", _case_rows(), ids=lambda c:c["id"])
def test_bernstein_pressure_matches_independent_quadrature(entry):
    from scipy.integrate import quad
    from compose_eos.whole_crust import WholeCrustBarotrope
    m=_source()
    c=json.loads((ROOT/entry["contract"]).read_text())
    bar=WholeCrustBarotrope.from_contract(c)
    q=bar.bridge;v=[c["lower_bridge"]["endpoint_states"][i][k] for i in (0,1)
                       for k in ("n","p","epsilon")]
    left=float(bar.outer.gamma[-1]*q.p_left/(q.epsilon_left+q.p_left))
    first=bar.upper_intervals[0]
    right=float(first.number(0)/(first.mu_left*first.susceptibility(0)))
    b=m.SlopeMatchedBernsteinInterval(*v,left_cs2=left,right_cs2=right)
    for x in (0.1,0.3,0.5,0.8,1.0):
        target_mu=b.mu_left+x*b.width
        independent=quad(lambda mu:b.evaluate_mu(mu)["n"],b.mu_left,target_mu,
                         epsabs=1e-13,epsrel=2e-12,limit=200)[0]
        analytical=b.pressure(x)-b.p_left
        assert abs(analytical-independent) < 1e-9*max(analytical,b.delta_p)
    for x in (0.15,0.55,0.85):
        h=1e-5
        numerical=(b.pressure(x+h)-b.pressure(x-h))/(2*h*b.width)
        assert abs(numerical/b.number(x)-1)<3e-6


def test_upper_power_law_fractional_exponent_has_endpoint_singularity():
    m=_source()
    c=json.loads((ROOT/_case_rows()[0]["contract"]).read_text())
    j=c["upper"]["sharp_interfaces"][0]
    v=[c["upper"][k][ii] for ii in (j+1,j+2) for k in ("n","p","epsilon")]
    b=m.CausalShapeInterval(*v,family="power")
    assert 0<b.shape_parameter<1
    assert b.validate()["status"]=="FAIL"
    assert b.evaluate_mu(b.mu_left)["cs2_over_c2"]==0
