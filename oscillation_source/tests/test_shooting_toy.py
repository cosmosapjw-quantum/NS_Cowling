"""Analytic BVP oracle; the numerical solution is computed in Fortran."""
import math, os
from pathlib import Path
import pytest
from tools.runtime import CaseRunner

@pytest.fixture(scope='module')
def shoot(tmp_path_factory):
    runner=CaseRunner(Path(__file__).resolve().parents[1],tmp_path_factory.mktemp('shooting'),
                      workspace=os.environ.get('NS_WORKSPACE','solutions'))
    art=runner.build_fpm()
    return lambda pars={},case='string_spectrum':runner.run_artifact(art,case,pars)

def test_known_eigenvalues(shoot):
    r=shoot({'length_m':1.0,'k_min':0.2,'k_max':8.0,'scan_n':200})
    assert r['status']=='ok',r
    assert abs(r['k1_per_m']-math.pi)<1e-7
    assert abs(r['k2_per_m']-2*math.pi)<1e-7
    assert r['roots_found']==2 and not r['spectrum_complete_claim']

def test_zero_parameter_is_not_an_eigenvalue(shoot):
    r=shoot({'k_per_m':0.,'length_m':2.},'string_shot')
    assert r['status']=='ok',r
    assert r['boundary_y_m']==pytest.approx(2.,rel=2e-13)
    assert r['residual']==pytest.approx(1.,rel=2e-13)

@pytest.mark.parametrize('length,k',[(1.,0.7),(2.5,2.1),(0.2,20.),(1.,1e-10)])
def test_shot_against_independent_closed_form(shoot,length,k):
    r=shoot({'length_m':length,'k_per_m':k,'emit_profile':True},'string_shot')
    assert r['status']=='ok',r
    assert r['boundary_y_m']==pytest.approx(math.sin(k*length)/k,rel=2e-9,abs=2e-11)
    assert r['boundary_slope']==pytest.approx(math.cos(k*length),rel=2e-9,abs=2e-11)
    for i in range(1,r['profile_n']+1):
        x=r[f'x_{i:05d}']
        assert r[f'y_{i:05d}']==pytest.approx(math.sin(k*x)/k,rel=3e-8,abs=3e-10)
    assert r['x_end_m']==length

@pytest.mark.parametrize('amplitude',[1e-100,1e100,-1.,-1e-100,-1e100])
def test_global_amplitude_not_a_root_tuning_parameter(shoot,amplitude):
    a=shoot({'k_per_m':2.3},'string_shot')
    b=shoot({'k_per_m':2.3,'amplitude':amplitude},'string_shot')
    assert b['status']=='ok',b
    assert b['residual']==pytest.approx(math.copysign(1.,amplitude)*a['residual'],abs=2e-14)
    assert b['boundary_y_m']/amplitude==pytest.approx(a['boundary_y_m'],rel=2e-13)
    assert a['accepted']==b['accepted']

@pytest.mark.parametrize('length',[0.5,2.,3.5])
def test_length_scaling_for_first_two_roots(shoot,length):
    r=shoot({'length_m':length,'k_min':0.2/length,'k_max':8./length})
    assert r['status']=='ok',r
    assert r['k1_per_m']*length==pytest.approx(math.pi,abs=4e-8)
    assert r['k2_per_m']*length==pytest.approx(2*math.pi,abs=4e-8)
    assert r['max_root_residual']<=1e-8 and r['max_bracket_width']<3e-10

def test_negative_amplitude_preserves_all_roots(shoot):
    a=shoot();b=shoot({'amplitude':-7.})
    assert a['status']==b['status']=='ok'
    assert b['k1_per_m']==pytest.approx(a['k1_per_m'],abs=2e-11)
    assert b['k2_per_m']==pytest.approx(a['k2_per_m'],abs=2e-11)

def test_coarse_cell_can_miss_two_roots(shoot):
    coarse=shoot({'k_min':1.,'k_max':8.,'scan_n':2})
    fine=shoot({'k_min':1.,'k_max':8.,'scan_n':129})
    assert coarse['status']=='no_bracket' and coarse['roots_found']==0
    assert 'k1_per_m' not in coarse
    assert fine['status']=='ok' and fine['roots_found']==2
    assert not fine['spectrum_complete_claim']

def test_no_root_in_requested_interval(shoot):
    r=shoot({'k_min':0.,'k_max':2.})
    assert r['status']=='no_bracket' and r['roots_found']==0 and 'k1_per_m' not in r

def test_three_level_rk4_eigenvalue_convergence(shoot):
    v=[shoot({'method':'rk4','n_steps':n,'scan_n':41}) for n in (32,64,128)]
    assert all(x['status']=='ok' for x in v),v
    errors=[abs(x['k2_per_m']-2*math.pi) for x in v]
    assert errors[0]>errors[1]>errors[2]
    assert 12<errors[0]/errors[1]<20 and 12<errors[1]/errors[2]<20,errors

def test_scan_grid_refinement_preserves_roots(shoot):
    a=shoot({'scan_n':81});b=shoot({'scan_n':161})
    assert a['status']==b['status']=='ok'
    assert a['k1_per_m']==pytest.approx(b['k1_per_m'],abs=3e-10)
    assert a['k2_per_m']==pytest.approx(b['k2_per_m'],abs=3e-10)

def test_root_sized_samples_are_not_automatically_certified(shoot):
    r=shoot({'k_min':0.,'k_max':math.pi-1e-7,'zero_tol':1e-6,'scan_n':81})
    assert r['unresolved_candidates']>=1 and r['status']!='ok'
    assert 'k1_per_m' not in r

@pytest.mark.parametrize('p',[
    {'length_m':0.},{'length_m':-1.},{'amplitude':0.},{'method':'unknown'},
    {'k_per_m':-1.},{'ode_rtol':-1.},{'atol_y_m':0.},{'max_steps':1},
    {'special_input':'nan_k'},{'special_input':'failed_rhs'},
    {'special_input':'nan_rhs'}, {'method':'rk4','n_steps':0},
])
def test_invalid_or_failed_shot_has_no_physical_output(shoot,p):
    r=shoot(p,'string_shot')
    assert r['status']!='ok' and not r['result_valid'],r
    assert 'boundary_y_m' not in r and 'residual' not in r

@pytest.mark.parametrize('p',[
    {'scan_n':0},{'scan_n':1},{'k_min':-1.},{'k_min':4.,'k_max':2.},
    {'k_min':1.,'k_max':1.},{'root_xtol':-1.},{'root_rtol':-1.},
    {'max_iter':1},{'zero_tol':-1.},{'residual_tol':0.},
    {'special_input':'fail_scan_midpoint'},
])
def test_failed_spectrum_never_reuses_success(shoot,p):
    ok=shoot({'scan_n':31});assert ok['status']=='ok'
    r=shoot(p)
    assert r['status']!='ok' and not r['result_valid'],r
    assert 'k1_per_m' not in r

def test_failed_refinement_reason_is_preserved(shoot):
    r=shoot({'max_iter':1,'scan_n':21})
    assert r['status']=='max_iterations'
    assert r['refinement_status_00001']=='max_iterations'
    assert r['unresolved_candidates']==r['candidates']

def test_direct_refinement_reintegrates_actual_ivp(shoot):
    r=shoot({'k_min':3.,'k_max':4.,'emit_profile':True},'string_refine')
    assert r['status']=='ok',r
    assert r['root_per_m']==pytest.approx(math.pi,abs=2e-8)
    assert abs(r['residual'])<1e-8
    assert r['profile_n']>5

@pytest.mark.parametrize('special',['nan_residual_tol','positive_ftol','sample_candidate','unknown_candidate'])
def test_direct_refinement_rejects_unverifiable_input(shoot,special):
    r=shoot({'special_input':special,'k_min':3.,'k_max':4.},'string_refine')
    assert r['status']!='ok' and not r['result_valid'],r
    assert 'root_per_m' not in r

def test_late_interval_is_not_labeled_as_first_physical_mode(shoot):
    r=shoot({'k_min':4.,'k_max':10.})
    assert r['status']=='ok' and r['roots_found']==2
    assert r['k1_per_m']==pytest.approx(2*math.pi,abs=4e-8)
    assert r['k2_per_m']==pytest.approx(3*math.pi,abs=4e-8)
    assert 'no stellar mode labels' in r['labels']

def test_refine_recomputes_signs_instead_of_trusting_a_reported_bracket(shoot):
    # Fabricated signs in an external scan cannot override the actual IVP residual.
    r=shoot({'k_min':1.,'k_max':2.},'string_refine')
    assert r['status']=='no_bracket' and not r['result_valid']
    assert 'root_per_m' not in r

def test_root_profiles_match_two_distinct_eigenfunctions(shoot):
    for lo,hi,n in ((3.,4.,1),(6.,7.,2)):
        r=shoot({'k_min':lo,'k_max':hi,'emit_profile':True},'string_refine')
        assert r['status']=='ok'
        for i in range(1,r['profile_n']+1):
            x=r[f'x_{i:05d}']; y=r[f'y_{i:05d}']
            assert y==pytest.approx(math.sin(n*math.pi*x)/(n*math.pi),abs=3e-9)

def test_failed_trials_are_retained_in_scan_history(shoot):
    r=shoot({'scan_n':21,'special_input':'fail_scan_midpoint','emit_profile':True})
    assert r['status']=='out_of_domain' and not r['result_valid']
    assert r['trial_n']==21 and r['trial_status_00011']=='out_of_domain'
    assert 'trial_residual_00011' not in r and 'k1_per_m' not in r

def test_division_by_endpoint_slope_creates_pole_not_eigenvalue(shoot):
    a=shoot({'k_per_m':1.5},'string_shot');b=shoot({'k_per_m':1.6},'string_shot')
    assert a['status']==b['status']=='ok'
    assert a['boundary_y_m']/a['boundary_slope']>0 and b['boundary_y_m']/b['boundary_slope']<0
    assert a['residual']>0 and b['residual']>0
    r=shoot({'k_min':1.5,'k_max':1.6},'string_refine')
    assert r['status']=='no_bracket' and not r['result_valid']
