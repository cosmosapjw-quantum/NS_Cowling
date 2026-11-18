"""Density holdout comparisons must use actual native states and immutable phase IDs."""
import importlib
import numpy as np
import pytest
cv=lambda:importlib.import_module("benchmark_native_dense_eos_20261118")
def test_loads_actual_rmf_uniform_and_core_tables():
    m=cv()
    a=m.load_native("rmf","NL3_PROVIDER_PRECISION")
    b=m.load_native("core","GM_GM1")
    assert len(a["n"])==4097 and len(b["n"])==1025
    assert np.all(np.diff(a["mu"])>0)
    assert len(a["sha256"])==64
def test_apr_phase_is_not_blended_across_6_8_7():
    m=cv()
    a=m.load_native("core","APR_APR")
    segments=m.phase_runs(a["phase"])
    assert len(segments)==3
    assert [int(a["phase"][i]) for i,j in segments]==[6,8,7]
    assert sum(j-i for i,j in segments)==len(a["n"])
def test_training_holdout_disjoint():
    m=cv()
    tr,hold=m.train_holdout_indices(103,stride=8)
    assert len(tr)>10 and len(hold)>0
    assert set(tr).isdisjoint(set(hold))
    assert len(set(tr)|set(hold))==103
    assert tr[0]==0 and tr[-1]==102
def test_hermite_predicts_finite_pressure_and_euler_identity():
    m=cv()
    a=m.load_native("rmf","PK1_CANONICAL")
    tr,hold=m.train_holdout_indices(len(a["n"]),stride=16)
    pred=m.predict_global("hermite",a,tr,hold[:25])
    assert np.all(np.isfinite(pred["p"]))
    assert np.max(np.abs(pred["mu"]*pred["n"]-pred["p"]-pred["epsilon"]))<1e-9
def test_linear_pressure_interpolant_has_no_finite_cs2():
    m=cv()
    a=m.load_native("core","RG_SLY4")
    tr,hold=m.train_holdout_indices(len(a["n"]),stride=8)
    pred=m.predict_global("linear",a,tr,hold[:20])
    assert not np.any(np.isfinite(pred["cs2"]))

def test_tf_receipt_serializes_numpy_boolean(tmp_path):
    module=importlib.import_module("benchmark_tf_native_holdout_20261118")
    output=tmp_path/"numpy_receipt.json"
    module.save(output,{"accepted":np.bool_(True),"value":np.float64(0.125)})
    import json
    assert json.loads(output.read_text())=={"accepted":True,"value":0.125}

def test_densified_rmf_holdout_provenance_and_data_leakage():
    m=cv()
    a=m.load_native("rmf","NL3_PROVIDER_PRECISION")
    train,hold=m.train_holdout_indices(len(a["n"]),stride=16)
    baseline=m.predict_global("hermite",a,train,hold[:100])
    altered=dict(a)
    import numpy as np
    for field in ("n","p","epsilon"):
        altered[field]=a[field].copy()
        altered[field][hold]*=1.3
    modified=m.predict_global("hermite",altered,train,hold[:100])
    for field in ("p","n","epsilon","cs2"):
        assert np.array_equal(baseline[field],modified[field]), field

def test_report_excludes_linear_in_every_figure():
    module=importlib.import_module("aggregate_native_dense_eos_20261118")
    assert "linear" not in module.PLOT_METHODS
    assert len(module.PLOT_METHODS)==8
    assert "linear" in module.DIAGNOSTIC_METHODS
