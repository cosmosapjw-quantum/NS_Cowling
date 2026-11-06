import os
import math
from pathlib import Path
import pytest
from tools.runtime import CaseRunner

ROOT=Path(__file__).resolve().parents[1]
H=6.62607015e-34
HBAR=H/(2*math.pi)
C=299792458.0
ME=9.1093837139e-31
MEC2=ME*C*C

@pytest.fixture(scope='module')
def eos(tmp_path_factory):
    r=CaseRunner(ROOT,tmp_path_factory.mktemp('beta'),workspace=os.environ.get('NS_WORKSPACE', 'learner' if (ROOT/'WORKSPACE.json').is_file() else 'solutions'))
    art=r.build_fpm()
    return lambda case,pars={}: r.run_artifact(art,case,pars)

def electron_oracle(n):
    if n==0.0:return MEC2,0.0,0.0
    pf=HBAR*(3*math.pi**2*n)**(1/3)
    x=pf/(ME*C)
    mu=MEC2*math.sqrt(1+x*x)
    # Direct exact-form subtraction loses ~1e-8 relative pressure accuracy at x~0.012.
    # Use the independently derived low-x expansion over the cancellation-prone range.
    if x < 5e-2:
        e=n*MEC2*(1+3*x*x/10-3*x**4/56+x**6/48-15*x**8/1408)
        p=n*MEC2*(x*x/5-x**4/14+x**6/24-5*x**8/176)
    else:
        pref=ME**4*C**5/(8*math.pi**2*HBAR**3)
        e=pref*(x*(1+2*x*x)*math.sqrt(1+x*x)-math.asinh(x))
        p=n*mu-e
    return mu,e,p

@pytest.mark.parametrize('n',[0.0,1e30,1e36,1e42,1e44,1e46])
def test_electron_gas_against_independent_oracle(eos,n):
    r=eos('electron_gas',{'n_e_m3':n})
    assert r['status']=='ok'
    mu,e,p=electron_oracle(n)
    assert r['mu_e_j']==pytest.approx(mu,rel=3e-12,abs=1e-30)
    assert r['energy_j_m3']==pytest.approx(e,rel=3e-10,abs=1e-18)
    assert r['pressure_pa']==pytest.approx(p,rel=3e-9,abs=1e-18)
    assert r['thermo_identity_closed'] is True

def test_electron_negative_density_is_rejected(eos):
    r=eos('electron_gas',{'n_e_m3':-1.0})
    assert r['status']=='invalid_input'

def test_beta_equilibrium_at_saturation(eos):
    r=eos('beta_equilibrium',{'n_b_fm3':0.16,'k0_mev':240.0,'s0_mev':30.0,'regulator_c':0.2})
    assert r['status']=='ok'
    assert 0.0 < r['proton_fraction'] < 0.5
    assert r['alpha']==pytest.approx(1-2*r['proton_fraction'],rel=2e-14)
    assert r['charge_neutrality_closed'] is True
    assert abs(r['beta_residual_mev']) < 1e-8
    assert r['muons_included'] is False
    assert r['energy_j_m3']>0 and r['pressure_pa']>0

def test_beta_pressure_matches_equilibrium_first_law(eos):
    n0=0.16
    h=1e-4
    m=eos('beta_equilibrium',{'n_b_fm3':n0*(1-h)})
    c=eos('beta_equilibrium',{'n_b_fm3':n0})
    p=eos('beta_equilibrium',{'n_b_fm3':n0*(1+h)})
    assert m['status']==c['status']==p['status']=='ok'
    nm=n0*(1-h)*1e45; nc=n0*1e45; np=n0*(1+h)*1e45
    em=m['energy_j_m3']/nm; ep=p['energy_j_m3']/np
    dedn=(ep-em)/(np-nm)
    rhs=c['pressure_pa']/nc**2
    assert dedn==pytest.approx(rhs,rel=3e-5)

def test_beta_fraction_changes_with_density_without_muon_threshold_logic(eos):
    a=eos('beta_equilibrium',{'n_b_fm3':0.08})
    b=eos('beta_equilibrium',{'n_b_fm3':0.32})
    assert a['status']==b['status']=='ok'
    assert a['proton_fraction'] != pytest.approx(b['proton_fraction'],rel=1e-5)
    assert a['muons_included'] is False and b['muons_included'] is False

def test_beta_equilibrium_reports_causal_local_response(eos):
    r=eos('beta_equilibrium',{'n_b_fm3':0.32})
    assert r['status']=='ok'
    assert 0.0 < r['chi_equilibrium'] < 1.0
    assert r['causal_local'] is True
