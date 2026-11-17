import math
import numpy as np
import pytest
from scipy.optimize import minimize_scalar
from compose_eos import gm1_uniform as gm1, bbb_uniform as bbb
from compose_eos.cold_leptons import free_fermion,lepton_state


def test_gm1_saturation_published_rounded_parameters():
    energy=lambda n:gm1.baryons(n,.5)["epsilon_baryons"]/n-gm1.MASS
    minimum=minimize_scalar(energy,bounds=(.1,.2),method="bounded")
    state=gm1.baryons(minimum.x,.5)
    h=1e-4
    incompressibility=9*minimum.x**2*(energy(minimum.x+h)-2*energy(minimum.x)+energy(minimum.x-h))/h**2
    assert minimum.x==pytest.approx(.153,abs=.0003)
    assert minimum.fun==pytest.approx(-16.3,abs=.07)
    assert state["effective_mass_mev"]/gm1.MASS==pytest.approx(.7,abs=.001)
    assert incompressibility==pytest.approx(300,abs=1)


@pytest.mark.parametrize("module,n",[(gm1,.2),(gm1,.8),(bbb,.3),(bbb,1.2)])
def test_independent_density_derivative_and_composition(module,n):
    state=module.beta_state(n)
    h=n*2e-5
    derivative=(module.beta_state(n+h)["epsilon"]-module.beta_state(n-h)["epsilon"])/(2*h)
    assert derivative==pytest.approx(state["mu_n"],rel=3e-6)
    assert state["Y_p"]==pytest.approx(state["Y_e"]+state["Y_mu"],abs=1e-11)
    assert state["Y_n"]+state["Y_p"]==1
    assert abs(state["euler_residual"])<2e-8


def test_bbb_published_composition_table():
    table=bbb.source_tables()[0]
    errors=[abs(bbb.beta_state(n)["Y_p"]-target)
            for n,target in zip(table["n"],table["table3_proton_fraction"])]
    assert max(errors)<.001
    for n,e0,en in zip(table["n"],table["symmetric_EA"],table["neutron_EA"]):
        assert bbb.baryons(n,.5)["epsilon_baryons"]/n-bbb.MASS==pytest.approx(e0,abs=1e-10)
        assert bbb.baryons(n,0.)["epsilon_baryons"]/n-bbb.MASS==pytest.approx(en,abs=1e-10)


@pytest.mark.parametrize("n",[0.,-1.,math.nan,math.inf])
def test_invalid_gm1_density(n):
    with pytest.raises(ValueError):
        gm1.baryons(n,.1)


@pytest.mark.parametrize("n",[.079,1.501,math.nan,math.inf])
def test_bbb_support(n):
    with pytest.raises(ValueError):
        bbb.baryons(n,.1)


def test_lepton_threshold_and_low_density():
    assert lepton_state(105.,105.658)==(0.,0.,0.)
    energy,pressure,chemical=free_fermion(1e-18,105.658)
    assert energy>0 and pressure>=0
    assert energy==pytest.approx(105.658e-18,rel=1e-9)


def test_gm1_species_contract():
    npe=gm1.beta_state(.3)
    npemu=gm1.beta_state(.3,include_muons=True)
    assert npe["Y_mu"]==0.
    assert npemu["Y_mu"]>0.
    assert npemu["Y_p"]>npe["Y_p"]


def test_generated_eos_internal_serialization(tmp_path):
    import importlib.util
    from pathlib import Path
    path=Path(__file__).resolve().parents[2]/"scripts/build_reference_cores.py"
    spec=importlib.util.spec_from_file_location("reference_core_builder",path)
    builder=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    rows=[gm1.beta_state(n) for n in (.1,.2,.3)]
    data=builder.eos_data(rows,gm1.MASS)
    data.write_internal_directory(tmp_path/"internal")
    assert data.q_extra.shape==(3,)
    assert (tmp_path/"internal/eos_data.csv").is_file()


@pytest.mark.parametrize("chemical,mass",[(math.nan,1.),(math.inf,1.),(1.,0.),(1.,math.nan)])
def test_invalid_lepton_inputs(chemical,mass):
    with pytest.raises(ValueError):
        lepton_state(chemical,mass)
