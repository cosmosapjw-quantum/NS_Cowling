from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"pipeline"))
from compose_cowling.adapter import build_adapter

def test_q_linear_exact_first_law_and_causality():
    mu=np.linspace(930.0,1500.0,80)
    q=0.0012+1.1e-6*(mu-930.0)
    n=mu*q
    p=np.ones_like(mu)
    for i,h in enumerate(np.diff(mu)):
        dq=q[i+1]-q[i]
        p[i+1]=p[i]+h*((mu[i]+.5*h)*q[i]+(.5*mu[i]+h/3)*dq)
    eps=mu*n-p
    eos,metrics=build_adapter(mu,n,p,eps)
    assert metrics["adapter_first_law_max_relative"] < 1e-13
    assert metrics["min_cs2"] > 0
    assert metrics["max_cs2"] <= 1+1e-13
    assert metrics["max_relative_hermite_pressure_residual"] < 2e-12

def test_pav_is_noop_for_monotone_q():
    mu=np.linspace(930.0,1600.0,90)
    q=0.001+np.arange(len(mu))*1e-6
    n=mu*q
    p=np.ones_like(mu)
    for i,h in enumerate(np.diff(mu)):
        dq=q[i+1]-q[i]
        p[i+1]=p[i]+h*((mu[i]+.5*h)*q[i]+(.5*mu[i]+h/3)*dq)
    eps=mu*n-p
    _,metrics=build_adapter(mu,n,p,eps)
    assert metrics["source_q_monotone"]
    assert metrics["pav_pooled_knots"] == 0
    assert metrics["max_abs_delta_n_over_n_source"] < 1e-13

def test_reject_non_first_law_source():
    mu=np.linspace(930.0,1300.0,12)
    q=np.full_like(mu,0.001)
    n=mu*q
    p=np.linspace(0.01,0.02,len(mu))
    eps=mu*n-p + 1.0
    try:
        build_adapter(mu,n,p,eps)
    except ValueError as exc:
        assert "first-law" in str(exc)
    else:
        raise AssertionError("bad source knots were accepted")

def test_existing_cohort_actual_table_matches_saved_adapter():
    original=ROOT/"data"/"eos_original_normalized"/"GPPVA_DD2__BBB_BHF2__css.csv"
    saved=ROOT/"data"/"eos_cowling_adapter_v1"/"GPPVA_DD2__BBB_BHF2__css.csv"
    assert original.is_file() and saved.is_file()
    from compose_cowling.adapter import build_adapter_csv
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        target=Path(tmp)/"out.csv"
        build_adapter_csv(original,target)
        got=np.genfromtxt(target,delimiter=",",names=True)
        ref=np.genfromtxt(saved,delimiter=",",names=True)
        for name in ("n_b_fm3","pressure_mev_fm3","energy_mev_fm3"):
            assert np.allclose(got[name],ref[name],rtol=3e-13,atol=1e-14)
