"""T07 Cowling RHS, center regularity, and fixed-background access."""
import csv, math, os
from pathlib import Path
import pytest
from tools.runtime import CaseRunner

ROOT=Path(__file__).resolve().parents[1]
BG=ROOT/'fixtures/fixed_background/smooth_poly_t05/background.csv'

@pytest.fixture(scope='module')
def cow(tmp_path_factory):
    runner=CaseRunner(ROOT,tmp_path_factory.mktemp('cowling-rhs'),workspace=os.environ.get('NS_WORKSPACE','solutions'))
    art=runner.build_fpm()
    return lambda pars={},case='cowling_rhs_point':runner.run_artifact(art,case,pars)

def _row(i=200):
    with BG.open(newline='') as f:
        rows=list(csv.DictReader(f))
    return {k:float(v) for k,v in rows[i].items()}

def test_center_regular_l2(cow):
    r=cow({'r0_m':2.0,'amplitude':3.0,'ell':2},'cowling_center')
    assert r['status']=='ok'
    assert r['w_m3']==pytest.approx(24.0)
    assert r['v_m2']==pytest.approx(-6.0)
    assert r['w_over_r3']==pytest.approx(3.0)
    assert r['v_over_r2']==pytest.approx(-1.5)

def test_center_positive_scale_is_linear(cow):
    a=cow({'r0_m':1.7,'amplitude':1.0},'cowling_center')
    b=cow({'r0_m':1.7,'amplitude':1e20},'cowling_center')
    assert a['status']==b['status']=='ok'
    assert b['w_m3']/1e20==pytest.approx(a['w_m3'],rel=2e-14)
    assert b['v_m2']/1e20==pytest.approx(a['v_m2'],rel=2e-14)

@pytest.mark.parametrize('p',[{'r0_m':0.0},{'r0_m':-1.0},{'amplitude':0.0},{'ell':0},{'special_input':'nan_r'}])
def test_center_rejects_invalid_input(cow,p):
    r=cow(p,'cowling_center')
    assert r['status']!='ok' and not r['result_valid']
    assert 'w_m3' not in r

def test_rhs_point_matches_independent_formula(cow):
    x=_row(200); omega=2*math.pi*1800.0; W=2.5e10; V=-3.2e6
    pars={'r_m':x['r_m'],'omega_rad_s':omega,'w_m3':W,'v_m2':V,
          'chi':x['chi'],'nu':x['nu'],'lambda_value':x['lambda'],'dnu_dr':x['dnu_dr']}
    r=cow(pars)
    assert r['status']=='ok',r
    k=omega/299792458.0; L=6.0
    q=k*k*x['r_m']**2*math.exp(x['lambda']-2*x['nu'])*V+x['dnu_dr']*W
    dw=q/x['chi']-L*math.exp(x['lambda'])*V
    dv=2*x['dnu_dr']*V-math.exp(x['lambda'])*W/x['r_m']**2
    assert r['dw_dr_m2']==pytest.approx(dw,rel=3e-13)
    assert r['dv_dr_m']==pytest.approx(dv,rel=3e-13)
    assert r['q_m2']==pytest.approx(q,rel=3e-13)

def test_rhs_lapse_frequency_covariance(cow):
    x=_row(350); om=2*math.pi*1200; C=0.37
    base={'r_m':x['r_m'],'omega_rad_s':om,'w_m3':4e9,'v_m2':-2e6,'chi':x['chi'],
          'nu':x['nu'],'lambda_value':x['lambda'],'dnu_dr':x['dnu_dr']}
    a=cow(base); b=cow({**base,'omega_rad_s':math.exp(C)*om,'nu':x['nu']+C})
    assert a['status']==b['status']=='ok'
    assert b['dw_dr_m2']==pytest.approx(a['dw_dr_m2'],rel=5e-13)
    assert b['dv_dr_m']==pytest.approx(a['dv_dr_m'],rel=5e-13)

def test_hz_passed_as_omega_changes_physics_and_is_not_auto_converted(cow):
    x=_row(250); base={'r_m':x['r_m'],'w_m3':0.0,'v_m2':-2e6,'chi':x['chi'],
                      'nu':x['nu'],'lambda_value':x['lambda'],'dnu_dr':x['dnu_dr']}
    a=cow({**base,'omega_rad_s':1000.0}); b=cow({**base,'omega_rad_s':2*math.pi*1000.0})
    assert a['status']==b['status']=='ok'
    assert b['q_m2']/a['q_m2']==pytest.approx((2*math.pi)**2,rel=3e-13)

@pytest.mark.parametrize('p',[{'chi':0.0},{'chi':-1e-3},{'r_m':0.0},{'omega_rad_s':0.0},{'omega_rad_s':-1.0},{'special_input':'nan_w'}])
def test_rhs_rejects_invalid_physical_or_nonfinite_input(cow,p):
    r=cow(p)
    assert r['status']!='ok' and not r['result_valid']
    assert 'dw_dr_m2' not in r

def test_background_probe_exact_row(cow):
    x=_row(321)
    r=cow({'background_path':str(BG),'r_probe_m':x['r_m']},'background_probe')
    assert r['status']=='ok'
    for src,out in [('p_pa','pressure_pa'),('energy_j_m3','energy_j_m3'),('chi','chi'),('nu','nu'),('lambda','lambda_value'),('dnu_dr','dnu_dr')]:
        assert r[out]==pytest.approx(x[src],rel=3e-15,abs=1e-300)
    assert r['interpolation_id']=='linear_q_v1'

def test_background_probe_linear_midpoint(cow):
    a=_row(300);b=_row(301);r0=(a['r_m']+b['r_m'])/2
    r=cow({'background_path':str(BG),'r_probe_m':r0},'background_probe')
    assert r['status']=='ok'
    assert r['chi']==pytest.approx((a['chi']+b['chi'])/2,rel=5e-14)
    assert r['nu']==pytest.approx((a['nu']+b['nu'])/2,rel=5e-14)

def test_background_probe_rejects_extrapolation(cow):
    a=_row(0); b=_row(-1)
    for rr in (a['r_m']-1e-6,b['r_m']+1e-6):
        r=cow({'background_path':str(BG),'r_probe_m':rr},'background_probe')
        assert r['status']=='out_of_domain' and not r['result_valid']

def test_background_missing_is_explicit(cow,tmp_path):
    r=cow({'background_path':str(tmp_path/'missing.csv'),'r_probe_m':1.0},'background_probe')
    assert r['status']!='ok' and not r['result_valid']
