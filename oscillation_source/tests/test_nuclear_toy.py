import os
import math
from pathlib import Path
import pytest
from tools.runtime import CaseRunner

ROOT=Path(__file__).resolve().parents[1]

@pytest.fixture(scope='module')
def nuclear(tmp_path_factory):
    r=CaseRunner(ROOT,tmp_path_factory.mktemp('nuclear'),workspace=os.environ.get('NS_WORKSPACE', 'learner' if (ROOT/'WORKSPACE.json').is_file() else 'solutions'))
    art=r.build_fpm()
    return lambda case,pars={}: r.run_artifact(art,case,pars)

def test_regulated_saturation_constraints_close(nuclear):
    r=nuclear('nuclear_constraints',{'regulated':True,'k0_mev':240.0,'regulator_c':0.2})
    assert r['status']=='ok'
    assert r['model_ready'] is True
    assert abs(r['binding_residual_mev']) < 2e-10
    assert abs(r['pressure_residual_mev_fm3']) < 2e-10
    assert abs(r['k0_residual_mev']) < 2e-8
    assert r['regulator_fixed_externally'] is True
    assert r['regulator_c']==pytest.approx(0.2,rel=0,abs=1e-15)

def test_simple_model_matches_closed_form_calibration(nuclear):
    r=nuclear('nuclear_constraints',{'regulated':False,'k0_mev':400.0})
    assert r['status']=='ok'
    t=r['t0_mev']; be=-16.0; k0=400.0
    sigma=(k0+2*t)/(3*t-9*be)
    b=(sigma+1)/(sigma-1)*(t/3-be)
    a=be-5*t/3-b
    assert r['sigma']==pytest.approx(sigma,rel=3e-13)
    assert r['a_mev']==pytest.approx(a,rel=3e-13)
    assert r['b_mev']==pytest.approx(b,rel=3e-13)

def test_symmetry_energy_is_calibrated_at_saturation(nuclear):
    r=nuclear('nuclear_state',{'n_b_fm3':0.16,'alpha':0.7,'s0_mev':30.0})
    assert r['status']=='ok'
    assert r['symmetry_mev']==pytest.approx(30.0,rel=3e-13)
    assert r['energy_excess_mev']==pytest.approx(r['symmetric_energy_mev']+0.7**2*30.0,rel=4e-13)

def test_exact_kinetic_is_diagnostic_not_double_counted(nuclear):
    r=nuclear('nuclear_state',{'n_b_fm3':0.16,'alpha':1.0})
    assert r['status']=='ok'
    assert r['kinetic_double_count_detected'] is False
    assert r['exact_kinetic_excess_mev'] > 0.0
    # The parabolic model uses one alpha^2*S contribution, not exact kinetic + S again.
    assert r['energy_excess_mev']==pytest.approx(r['symmetric_energy_mev']+r['symmetry_mev'],rel=4e-13)

def test_isospin_evenness(nuclear):
    a=nuclear('nuclear_state',{'n_b_fm3':0.32,'alpha':0.4})
    b=nuclear('nuclear_state',{'n_b_fm3':0.32,'alpha':-0.4})
    assert a['status']==b['status']=='ok'
    assert a['energy_excess_mev']==pytest.approx(b['energy_excess_mev'],rel=4e-13)
    assert a['pressure_nuclear_pa']==pytest.approx(b['pressure_nuclear_pa'],rel=4e-13)

def test_unphysical_asymmetry_and_density_are_rejected(nuclear):
    for pars in ({'n_b_fm3':0.16,'alpha':1.01},{'n_b_fm3':0.0,'alpha':0.0},{'n_b_fm3':-0.1,'alpha':0.0}):
        r=nuclear('nuclear_state',pars)
        assert r['status'] in {'invalid_input','out_of_domain'}
        assert r.get('state_valid',False) is False

def test_regulator_is_not_inferred_from_three_saturation_constraints(nuclear):
    r1=nuclear('nuclear_constraints',{'regulated':True,'regulator_c':0.1,'k0_mev':240.0})
    r2=nuclear('nuclear_constraints',{'regulated':True,'regulator_c':0.2,'k0_mev':240.0})
    assert r1['status']==r2['status']=='ok'
    assert r1['regulator_fixed_externally'] and r2['regulator_fixed_externally']
    assert r1['sigma'] != pytest.approx(r2['sigma'],rel=1e-5)
    assert abs(r1['k0_residual_mev'])<2e-8 and abs(r2['k0_residual_mev'])<2e-8

def test_regulated_model_reports_causality_scan(nuclear):
    r=nuclear('nuclear_constraints',{'regulated':True,'regulator_c':0.2,'k0_mev':240.0,'scan_u_max':10.0})
    assert r['status']=='ok'
    assert r['scan_causal'] is True
    assert 0.0 < r['max_chi_scan'] <= 1.0
