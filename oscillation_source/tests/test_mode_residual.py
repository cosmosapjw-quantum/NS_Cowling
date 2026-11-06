"""T07 amplitude-independent surface residual and single Cowling shots."""
import math, os
from pathlib import Path
import pytest
from tools.runtime import CaseRunner

ROOT=Path(__file__).resolve().parents[1]
BG=ROOT/'fixtures/fixed_background/smooth_poly_t05/background.csv'

@pytest.fixture(scope='module')
def cow(tmp_path_factory):
    runner=CaseRunner(ROOT,tmp_path_factory.mktemp('cowling-residual'),workspace=os.environ.get('NS_WORKSPACE','solutions'))
    art=runner.build_fpm()
    return lambda pars={},case='surface_residual':runner.run_artifact(art,case,pars)

def base_surface():
    return {'r_m':7800.0,'omega_rad_s':2*math.pi*1500.0,'nu':-0.11,'lambda_value':0.11,
            'dnu_dr':1.5e-5,'w_m3':2e11,'v_m2':-3e6}

def test_residual_scales_with_orientation_only(cow):
    p=base_surface(); a=cow(p); b=cow({**p,'w_m3':p['w_m3']*1e40,'v_m2':p['v_m2']*1e40})
    n=cow({**p,'w_m3':-p['w_m3'],'v_m2':-p['v_m2']})
    assert a['status']==b['status']==n['status']=='ok'
    assert b['residual']==pytest.approx(a['residual'],abs=3e-14)
    assert n['residual']==pytest.approx(-a['residual'],abs=3e-14)
    assert abs(a['residual'])<=1

def test_residual_has_no_single_component_pole(cow):
    p=base_surface()
    for w,v in ((0.0,p['v_m2']),(p['w_m3'],0.0)):
        r=cow({**p,'w_m3':w,'v_m2':v})
        assert r['status']=='ok' and math.isfinite(r['residual'])

def test_residual_rejects_zero_state_or_zero_b(cow):
    p=base_surface()
    a=cow({**p,'w_m3':0.0,'v_m2':0.0})
    b=cow({**p,'dnu_dr':0.0,'omega_rad_s':0.0})
    assert a['status']!='ok' and not a['result_valid']
    assert b['status']!='ok' and not b['result_valid']

def test_residual_matches_independent_cosine_definition(cow):
    p=base_surface();r=cow(p);assert r['status']=='ok'
    k=p['omega_rad_s']/299792458.0
    z=(p['w_m3']/p['r_m']**3,p['v_m2']/p['r_m']**2)
    b=(p['r_m']*p['dnu_dr'],(k*p['r_m'])**2*math.exp(p['lambda_value']-2*p['nu']))
    expected=(b[0]*z[0]+b[1]*z[1])/(math.hypot(*b)*math.hypot(*z))
    assert r['residual']==pytest.approx(expected,abs=3e-14)

def test_fixed_background_single_shot_is_finite(cow):
    r=cow({'background_path':str(BG),'omega_rad_s':2*math.pi*1500.0,'emit_profile':True},'cowling_shot')
    assert r['status']=='ok',r
    assert r['result_valid'] is True and math.isfinite(r['residual']) and abs(r['residual'])<=1
    assert r['profile_n']>50 and r['r_start_m']==pytest.approx(1.0)
    assert r['r_end_m']==pytest.approx(7848.961769861312,rel=2e-14)
    assert r['surface_policy']=='finite_cutoff_schwarzschild_time'
    assert r['mode_calibrated'] is False

def test_shot_lapse_frequency_covariance(cow):
    om=2*math.pi*1400.0; C=0.2
    a=cow({'background_path':str(BG),'omega_rad_s':om},'cowling_shot')
    b=cow({'background_path':str(BG),'omega_rad_s':math.exp(C)*om,'nu_shift':C},'cowling_shot')
    assert a['status']==b['status']=='ok'
    assert b['residual']==pytest.approx(a['residual'],abs=2e-9)

def test_shot_positive_amplitude_does_not_tune_root(cow):
    a=cow({'background_path':str(BG),'omega_rad_s':2*math.pi*1000.0,'amplitude':1.0},'cowling_shot')
    b=cow({'background_path':str(BG),'omega_rad_s':2*math.pi*1000.0,'amplitude':3.0},'cowling_shot')
    assert a['status']==b['status']=='ok'
    assert b['residual']==pytest.approx(a['residual'],abs=3e-10)
    assert b['surface_w_m3']/a['surface_w_m3']==pytest.approx(3.0,rel=5e-8)

def test_negative_amplitude_flips_residual_not_root_condition(cow):
    p={'background_path':str(BG),'omega_rad_s':2*math.pi*1100.0}
    a=cow({**p,'amplitude':1.0},'cowling_shot');b=cow({**p,'amplitude':-1.0},'cowling_shot')
    assert a['status']==b['status']=='ok'
    assert b['residual']==pytest.approx(-a['residual'],abs=3e-10)

def test_invalid_shot_has_no_physical_output(cow):
    for p in ({'background_path':str(BG),'omega_rad_s':0.0},{'background_path':str(BG),'amplitude':0.0},
              {'background_path':str(BG),'ode_rtol':-1.0},{'background_path':'/does/not/exist.csv'}):
        r=cow(p,'cowling_shot')
        assert r['status']!='ok' and not r['result_valid']
        assert 'residual' not in r

def test_t07_is_not_a_mode_calibration(cow):
    r=cow({'background_path':str(BG),'omega_rad_s':2*math.pi*1500.0},'cowling_shot')
    assert r['status']=='ok'
    assert r['mode_calibrated'] is False
    assert r['identity_status']=='not_attempted'

def test_shot_rejects_unrepresentable_output_amplitude_as_numerics(cow):
    r=cow({'background_path':str(BG),'omega_rad_s':2*math.pi*1000.0,'amplitude':1e300},'cowling_shot')
    assert r['status']=='range_error'
    assert r['category']=='numerics'
    assert not r['result_valid']
    assert 'residual' not in r
