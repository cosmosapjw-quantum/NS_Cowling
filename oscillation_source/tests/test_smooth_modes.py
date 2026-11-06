"""T08 actual Cowling spectra, distinct root and identity contracts."""
from pathlib import Path
import os
import pytest
from tools.runtime import CaseRunner
ROOT=Path(__file__).resolve().parents[1]
BG=ROOT/'fixtures/fixed_background/smooth_poly_t05/background.csv'

@pytest.fixture(scope='module')
def modes(tmp_path_factory):
    r=CaseRunner(ROOT,tmp_path_factory.mktemp('modes'),workspace=os.environ.get('NS_WORKSPACE','solutions'))
    a=r.build_fpm()
    return lambda p={},case='smooth_spectrum': r.run_artifact(a,case,p,timeout=90)

def test_spectrum_feature_is_real(modes):
    r=modes({'background_path':str(BG)},'smooth_spectrum')
    assert r['status']=='ok',r
    assert r['f_identity']=='identified'
    assert r['p1_identity']=='identified'

def test_normalization_and_node_contract(modes):
    r=modes({},'mode_diagnostics')
    assert r['status']=='ok',r
    assert r['positive_overlap']==pytest.approx(1.,abs=1e-12)
    assert r['negative_overlap']==pytest.approx(-1.,abs=1e-12)
    assert r['norm_squared']==pytest.approx(1.,abs=1e-12)

def test_successful_root_with_ambiguous_identity_is_separate(modes):
    r=modes({},'mode_identity_probe')
    assert r['status']=='ok',r
    assert r['ambiguous_identity']=='identity_unresolved'
    assert r['root_status']=='ok'

def test_identical_rounded_profile_duplicate_is_redundant(modes):
    r=modes({'diagnostic_action':'duplicate_identical'},'mode_diagnostics')
    assert r['status']=='ok',r
    assert r['w_nodes']==r['v_nodes']==1

@pytest.mark.parametrize('action',['zero','duplicate_r'])
def test_invalid_profile_is_not_a_mode(modes,action):
    r=modes({'diagnostic_action':action},'mode_diagnostics')
    assert r['status']!='ok' and not r['result_valid']

@pytest.mark.parametrize('bounds,label', [((2800.,3000.),'f'),((6300.,6550.),'p1')])
def test_refined_modes_have_independent_node_evidence(modes,bounds,label):
    r=modes({'background_path':str(BG),'frequency_min_hz':bounds[0],'frequency_max_hz':bounds[1]},'mode_refine')
    assert r['status']=='ok',r
    assert r['mode_root_valid'] and r['mode_identity']=='identified'
    assert r['mode_label']==label
    assert r['mode_w_nodes']==r['mode_v_nodes']==(0 if label=='f' else 1)
    assert abs(r['mode_residual'])<1e-5
    assert r['mode_bracket_width_hz']<1e-7

@pytest.mark.parametrize('p',[
    {'frequency_min_hz':100.,'frequency_max_hz':200.},
    {'diagnostic_action':'sample'},
    {'max_iter':1},
    {'frequency_min_hz':0.},
    {'background_path':'/nonexistent/background.csv'},
    {'method':'unknown'},
    {'amplitude':0.},
    {'residual_tol':0.},
    {'root_xtol':0.,'root_rtol':0.},
])
def test_bad_refinement_has_no_identified_root(modes,p):
    r=modes({'background_path':str(BG),'frequency_min_hz':2800.,'frequency_max_hz':3000.,**p},'mode_refine')
    assert r['status']!='ok'
    assert not r['result_valid']
    assert not r.get('mode_root_valid',False)
    assert r.get('mode_identity')!='identified'

def test_root_precision_below_representable_spacing_is_rejected(modes):
    r=modes({'background_path':str(BG),'frequency_min_hz':2800.,'frequency_max_hz':3000.,
             'root_xtol':1e-16,'root_rtol':0.},'mode_refine')
    assert r['status']=='stagnation'
    assert not r['mode_root_valid']

@pytest.mark.parametrize('n',[1,0,2001])
def test_invalid_scan_size_is_rejected(modes,n):
    r=modes({'background_path':str(BG),'scan_n':n})
    assert r['status']=='invalid_input'

def test_one_mode_window_is_not_complete_two_mode_result(modes):
    r=modes({'background_path':str(BG),'frequency_min_hz':2700.,'frequency_max_hz':3100.,'scan_n':9})
    assert r['status']=='identity_unresolved'
    assert r['f_identity']=='identified'
    assert 'p1_hz' not in r

def test_finer_scan_finds_same_scoped_mode_families(modes):
    p={'background_path':str(BG),'scan_n':33}
    a=modes(p);b=modes({**p,'scan_n':65})
    assert a['status']==b['status']=='ok'
    for label in ('f','p1'):
        assert a[label+'_hz']==pytest.approx(b[label+'_hz'],abs=3e-6)
        assert a[label+'_label']==b[label+'_label']==label

def test_phase_and_lapse_transform_do_not_relabel_modes(modes):
    import math
    base={'background_path':str(BG),'frequency_min_hz':2800.,'frequency_max_hz':3000.}
    a=modes(base,'mode_refine')
    b=modes({**base,'amplitude':-7.},'mode_refine')
    c=modes({**base,'frequency_min_hz':2800.*math.exp(.1),'frequency_max_hz':3000.*math.exp(.1),
             'nu_shift':.1},'mode_refine')
    assert a['status']==b['status']==c['status']=='ok'
    assert a['mode_label']==b['mode_label']==c['mode_label']=='f'
    assert b['mode_hz']==pytest.approx(a['mode_hz'],abs=3e-6)
    assert c['mode_hz']/math.exp(.1)==pytest.approx(a['mode_hz'],abs=3e-6)

def test_nonmatching_node_pattern_is_unresolved(modes):
    r=modes({'diagnostic_action':'extra_node'},'mode_diagnostics')
    assert r['status']=='ok'
    assert r['w_nodes']==2 and r['v_nodes']==1
