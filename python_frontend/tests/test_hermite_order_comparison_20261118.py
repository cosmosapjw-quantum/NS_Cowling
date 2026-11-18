"""Cubic and quintic P(mu) Hermite contract tests for frozen native data."""
from pathlib import Path
import importlib
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]


def _mod():
    return importlib.import_module("compare_hermite_orders_20261118")


def test_quintic_exact_polynomial_and_two_derivatives():
    m=_mod()
    from scipy.interpolate import BPoly
    x=np.linspace(0.1,1.5,31)
    P=lambda z:8+5*z+z*z+.5*z**3+.02*z**4+.01*z**5
    n=lambda z:5+2*z+1.5*z*z+.08*z**3+.05*z**4
    chi=lambda z:2+3*z+.24*z*z+.20*z**3
    bp=m.quintic_hermite(x,P(x),n(x),chi(x))
    q=np.linspace(.13,1.45,117)
    np.testing.assert_allclose(bp(q),P(q),rtol=1e-11,atol=1e-11)
    np.testing.assert_allclose(bp(q,nu=1),n(q),rtol=1e-11,atol=1e-11)
    np.testing.assert_allclose(bp(q,nu=2),chi(q),rtol=1e-9,atol=1e-9)


def test_native_quintic_uses_only_training_knots():
    m=_mod()
    from benchmark_native_dense_eos_20261118 import load_native,train_holdout_indices
    source=load_native("rmf","PK1_CANONICAL")
    tr,hold=train_holdout_indices(len(source["n"]),16)
    old=m.predict("quintic",source,tr,hold[:60])
    modified=dict(source)
    for k in ("p","n","epsilon"):
        modified[k]=np.array(source[k]).copy()
        modified[k][hold]*=1.3
    new=m.predict("quintic",modified,tr,hold[:60])
    for k in ("p","n","epsilon","cs2"):
        np.testing.assert_array_equal(old[k],new[k])


def test_quintic_respects_cold_euler_identity_at_holdout():
    m=_mod()
    from benchmark_native_dense_eos_20261118 import load_native,train_holdout_indices
    data=load_native("core","RG_SLY4")
    tr,hold=train_holdout_indices(len(data["n"]),8)
    z=m.predict("quintic",data,tr,hold[:100])
    np.testing.assert_allclose(z["epsilon"]+z["p"],z["mu"]*z["n"],rtol=1e-12,atol=1e-10)


def test_cubic_and_quintic_are_distinct_functions():
    m=_mod()
    from benchmark_native_dense_eos_20261118 import load_native,train_holdout_indices
    data=load_native("rmf","NL3_PROVIDER_PRECISION")
    tr,hold=train_holdout_indices(len(data["n"]),16)
    a=m.predict("cubic",data,tr,hold[:100])
    b=m.predict("quintic",data,tr,hold[:100])
    assert np.max(np.abs(a["n"]-b["n"]))>1e-15

def test_run_receipt_serializes_numpy_vectors_and_booleans(tmp_path):
    import json
    import run_cubic_quintic_dense_cv_20261118 as m
    value=dict(n=np.asarray([.02,.03]),valid=np.bool_(True),a=np.float64(.12))
    f=tmp_path/"receipt.json"
    m.save(f,value)
    assert json.loads(f.read_text()) == dict(n=[.02,.03],valid=True,a=.12)

def test_matching_method_matrix_explicitly_contains_css_and_four_hermite_orders():
    import native_hermite_matching_matrix_20261118 as m
    assert m.METHODS == ("cubic_F","quintic_F","cubic_P","quintic_P","css")
    assert len(m.PAIRS)==4
    assert len({z[0] for z in m.PAIRS})==2
    assert len({z[1] for z in m.PAIRS})==2

def test_stellar_quintic_interval_exact_endpoints_and_causality():
    import stellar_quintic_contract_20261118 as m
    from compose_eos.whole_crust import WholeCrustBarotrope
    import json
    plan=json.loads((ROOT/"data/low_density_crust_20261117/stellar/causal_extension/TOV_PLAN.json").read_text())
    entry=next(v for v in plan["cases"] if v["label"]=="primary")
    data=json.loads((ROOT/entry["contract"]).read_text())
    bar=WholeCrustBarotrope.from_contract(data)
    changed,report=m.install_quintic_upper(bar)
    assert report["rejected"]==0
    assert report["installed"]==len(changed.upper_intervals)-len(changed.sharp_interfaces)
    assert changed.sharp_interfaces==bar.sharp_interfaces
    for i in (0,1,2,len(changed.upper_intervals)-1):
        if i in changed.sharp_interfaces:continue
        s=changed.upper_intervals[i]
        for t in (0.0,.25,.5,.75,1.0):
            mu=s.mu_left+t*s.width
            q=s.evaluate_mu(mu)
            assert 0 < q["cs2_over_c2"] <= 1+1e-12
            assert abs(q["epsilon"]+q["p"]-mu*q["n"])<1e-9
