"""Checks of screened-electron elimination and composition semantics."""
import importlib.util
from pathlib import Path
import json
import numpy as np
import pytest
from compose_eos.spatial_tf import RadialGrid,initial_profile,E2
from compose_eos.spatial_tf_equilibrium import (
    screened_electrons,ReducedCellFunctional,solve_screened_cell)


@pytest.mark.parametrize("dimension",[1,2,3])
def test_screened_electron_pde_neutrality_and_chemical_equilibrium(dimension):
    g=RadialGrid(48,dimension)
    protons=.002*(1+.7*np.cos(np.pi*g.x))
    ne,v,mu,res=screened_electrons(g,protons,25.)
    assert res<2e-10
    assert abs(g.mean(protons-ne))<2e-12
    assert abs(g.mean(v))<2e-11
    assert np.max(abs(g.minus_laplacian@v/25**2-4*np.pi*E2*(protons-ne)))<2e-10
    from compose_eos.spatial_tf import fermions,ELECTRON_MASS
    assert np.max(abs(fermions(ne,ELECTRON_MASS)["ef"]-v-mu))<2e-10


def test_empty_electron_branch():
    g=RadialGrid(16,3);ne,v,mu,res=screened_electrons(g,np.zeros(16),25.)
    assert np.array_equal(ne,np.zeros(16))
    assert np.array_equal(v,np.zeros(16))
    assert res==0


@pytest.mark.parametrize("model",["GPPVA_DD2","GPPVA_FSU2R"])
@pytest.mark.parametrize("dimension",[1,2,3])
def test_reduced_functional_derivative(model,dimension):
    g=RadialGrid(20,dimension);f=ReducedCellFunctional(model,g,.03)
    nn,np_,_=initial_profile(g,.03);nn+=1e-4;np_+=1e-4
    x=np.r_[nn*g.sqrtw/.03,np_*g.sqrtw/.03,np.log(35.)]
    _,jac=f.objective(x)
    for i in (0,12,20,31,40):
        h=1e-6;a=x.copy();b=x.copy();a[i]+=h;b[i]-=h
        finite=(f.objective(a)[0]-f.objective(b)[0])/(2*h)
        assert finite==pytest.approx(jac[i],rel=3e-5,abs=4e-6)


def test_reduced_cell_pressure_is_independent_density_derivative():
    n=.03;h=2e-5
    states=[solve_screened_cell("GPPVA_DD2",v,points=32) for v in (n-h,n,n+h)]
    mid=states[1]
    assert mid["numerical_status"]=="STATIONARY"
    assert mid["max_chemical_KKT_residual_MeV"]<.02
    assert mid["electron_pde_residual"]<2e-10
    pressure=n*(states[2]["epsilon"]-states[0]["epsilon"])/(2*h)-mid["epsilon"]
    assert pressure==pytest.approx(mid["pressure"],rel=3e-4,abs=3e-6)


def test_tf_companion_roundtrip_and_mismatch_rejection(tmp_path):
    # Native writer/readback, then a deliberately altered physical composition.
    root=Path(__file__).resolve().parents[2]
    spec=importlib.util.spec_from_file_location("tf_runner",root/"scripts/build_gppva_tf.py")
    runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
    states=[solve_screened_cell("GPPVA_DD2",n,points=16) for n in (.025,.03,.035)]
    for state in states:state["geometry"]=runner.physical_geometry(state)
    runner.dump(tmp_path/"selected_states.json.gz",states)
    runner.export("GPPVA_DD2",states,tmp_path)
    from compose_eos.tf_composition import read_tf_composition
    result=read_tf_composition(tmp_path)
    assert np.max(abs(result["fractions"][:,0]+result["fractions"][:,1]-1))<1e-11
    assert not np.any(result["thermodynamics"].Yn_present)
    meta=json.loads((tmp_path/"composition.json").read_text())
    meta["rows"][0]["Y_p_total"]+=.01
    (tmp_path/"composition.json").write_text(json.dumps(meta))
    with pytest.raises(ValueError,match="fraction/profile mismatch"):read_tf_composition(tmp_path)


@pytest.mark.parametrize("dimension",[1,2,3])
def test_screened_electrons_stale_empty_gauge_shifted_guess(dimension):
    from compose_eos.spatial_tf import ELECTRON_MASS
    g=RadialGrid(32,dimension);protons=.003*(1+.8*np.cos(np.pi*g.x))
    reference=screened_electrons(g,protons,35.)
    for potential,mu in [(np.zeros(32),ELECTRON_MASS),
                         (np.full(32,100.),ELECTRON_MASS-100.)]:
        actual=screened_electrons(g,protons,35.,initial=(potential,mu))
        assert actual[0]==pytest.approx(reference[0],rel=2e-10,abs=1e-13)
        assert actual[1]==pytest.approx(reference[1],rel=2e-10,abs=1e-10)
        # A 2e-12 PDE tolerance permits chemical shifts of residual/(4*pi*e2*chi_e).
        # This bound is much tighter than the nucleon stationarity criterion.
        assert actual[2]==pytest.approx(reference[2],abs=1e-8)
