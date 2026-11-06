"""Exact solutions, order checks and failure-state contracts; no Python ODE solver."""
import math
import os
import pytest
from pathlib import Path
from tools.runtime import CaseRunner

@pytest.fixture(scope='module')
def numerics(tmp_path_factory):
    root=Path(__file__).resolve().parents[1]
    runner=CaseRunner(root,tmp_path_factory.mktemp('ode-numerics'),workspace=os.environ.get('NS_WORKSPACE','solutions'))
    artifact=runner.build_fpm()
    return lambda case,pars:runner.run_artifact(artifact,case,pars)

def test_rk4_first_success(numerics):
    r=numerics('rk4_exp',{})
    assert r['status']=='ok'
    assert abs(r['y1']-math.e)<3e-7

def test_adaptive_first_success(numerics):
    r=numerics('adaptive_exp',{})
    assert r['status']=='ok'
    assert abs(r['y1']-math.e)<2e-8

@pytest.mark.parametrize('method,order,min_order',[('euler',1,0.9),('rk4',4,3.8)])
def test_exponential_global_order(numerics,method,order,min_order):
    errors=[]
    for n in (16,32,64):
        r=numerics('ode',{'method':method,'n_steps':n})
        assert r['status']=='ok' and r['last_x']==1.
        errors.append(abs(r['y1']-math.e))
    observed=[math.log(errors[i]/errors[i+1],2) for i in range(2)]
    assert all(min_order<q<order+0.2 for q in observed)

@pytest.mark.parametrize('n',[16,32,64])
def test_oscillator_phase_amplitude(numerics,n):
    r=numerics('ode',{'method':'rk4','problem':'oscillator','n_steps':n,'y2_initial':0.})
    assert r['status']=='ok'
    assert abs(r['y1']-math.cos(1.))<2/n**4
    assert abs(r['y2']+math.sin(1.))<2/n**4
    assert abs(math.atan2(-r['y2'],r['y1'])-1.)<2/n**4

def test_nonautonomous_polynomial(numerics):
    r=numerics('ode',{'method':'rk4','problem':'cubic','y_initial':0.,'n_steps':3})
    assert r['status']=='ok' and abs(r['y1']-0.25)<2e-15

def test_rk4_one_step_polynomial(numerics):
    r=numerics('ode',{'method':'rk4','n_steps':1,'x1':0.2})
    p=sum(0.2**j/math.factorial(j) for j in range(5))
    assert r['y1']==pytest.approx(p,abs=2e-15)

@pytest.mark.parametrize('method',['euler','rk4','adaptive'])
@pytest.mark.parametrize('problem,status',[('nan','nonfinite_result'),('failure','out_of_domain')])
def test_failed_rhs_not_accepted(numerics,method,problem,status):
    r=numerics('ode',{'method':method,'problem':problem})
    assert r['status']==status and not r['result_valid']
    assert r['accepted']==0 and r['last_x']==0. and 'y1' not in r

@pytest.mark.parametrize('params',[{'h_initial':0.},{'h_min':0.},{'h_max':-1.},{'h_min':2.,'h_max':1.},
                                  {'rtol':-1.},{'atol1':0.},{'max_steps':0.}])
def test_invalid_adaptive_options(numerics,params):
    if 'max_steps' in params:params={'max_steps':0}
    r=numerics('ode',dict(method='adaptive',**params))
    assert r['status']=='invalid_input' and r['rhs_calls']==0

def test_zero_error_controller(numerics):
    r=numerics('control',{'action':'factor','error_value':0.})
    assert r['status']=='ok' and r['factor']==5.

@pytest.mark.parametrize('error',[1e-300,0.5,1.,2.,1e300])
def test_bounded_controller(numerics,error):
    r=numerics('control',{'action':'factor','error_value':error})
    assert r['status']=='ok' and 0.2<=r['factor']<=5.
    if error>1:assert r['factor']<1.
    r=numerics('control',{'action':'factor','error_value':error,'previous_rejected':True})
    assert r['factor']<=1.

def test_scaled_component_norm(numerics):
    r=numerics('control',{'action':'norm','old1':1.,'new1':2.,'error1':2e-3,'atol1':2e-3,
                         'old2':1e-20,'new2':1e-20,'error2':3e-24,'atol2':1e-24,'rtol':1e-4})
    # component 2 = 3e-24 / (1e-24 + 1e-4*1e-20) = 1.5
    assert r['error_norm']==pytest.approx(1.5,rel=3e-15)

def test_bad_norm_size_rejected(numerics):
    r=numerics('control',{'action':'bad_size'})
    assert r['status']=='invalid_input'

def test_dopri_pair_nonautonomous_order(numerics):
    errors=[]
    for h in (0.5,0.25,0.125):
        r=numerics('dopri_pair',{'problem':'quartic','h_initial':h,'y_initial':0.})
        assert r['status']=='ok'
        assert r['y1']==pytest.approx(h**5/5,rel=2e-13)
        assert r['rhs_calls']==7
        errors.append(abs(r['error1']))
    assert errors[0]/errors[1]==pytest.approx(32.,rel=1e-10)
    assert errors[1]/errors[2]==pytest.approx(32.,rel=1e-10)

def test_rejection_keeps_state_and_profile(numerics):
    r=numerics('ode',{'method':'adaptive','h_initial':1.,'h_max':1.,'rtol':1e-11,'atol1':1e-13})
    assert r['status']=='ok' and r['rejected']>=1
    assert r['rejection_rollback_ok'] is True
    assert r['profile_points']==r['accepted']+1
    assert r['attempts']==r['accepted']+r['rejected']
    assert r['all_accepted_errors_le_one'] is True
    assert r['y1']==pytest.approx(math.e,rel=1e-10)

def test_componentwise_adaptive_scales(numerics):
    r=numerics('ode',{'method':'adaptive','problem':'two_scale','y2_initial':1e-20,
                     'atol1':1e-12,'atol2':1e-32,'rtol':1e-9})
    assert r['status']=='ok'
    assert r['y1']==pytest.approx(math.e,rel=1e-8)
    assert r['y2']/1e-20==pytest.approx(math.exp(-3),rel=1e-8)

def test_endpoint_clipped_not_overshot(numerics):
    r=numerics('ode',{'method':'adaptive','problem':'constant','x1':1.03,'h_initial':0.7,'h_max':1.})
    assert r['status']=='ok' and r['last_x']==1.03
    assert r['max_stage_x']<=1.03
    assert r['y1']==pytest.approx(2.03,abs=1e-14)

def test_known_stops_are_hit_exactly(numerics):
    r=numerics('ode',{'method':'adaptive','problem':'constant','h_initial':0.7,'use_events':True,
                      'event1':0.37,'event2':0.8})
    assert r['status']=='ok' and r['stops_hit']==2
    assert r['crossed_stop'] is False

@pytest.mark.parametrize('params',[{'event1':0.8,'event2':0.37},{'event1':0.37,'event2':0.37},
                                   {'event1':-0.1,'event2':0.5},{'event1':0.5,'event2':1.1}])
def test_bad_stop_schedule(numerics,params):
    r=numerics('ode',dict(method='adaptive',use_events=True,**params))
    assert r['status']=='invalid_input' and r['rhs_calls']==0

def test_backwards_integration(numerics):
    r=numerics('ode',{'method':'adaptive','x0':1.,'x1':0.,'y_initial':math.e,
                     'use_events':True,'event1':0.8,'event2':0.37})
    assert r['status']=='ok' and r['last_x']==0. and r['stops_hit']==2
    assert r['y1']==pytest.approx(1.,rel=2e-8)

def test_zero_interval_does_not_call_rhs(numerics):
    r=numerics('ode',{'method':'adaptive','x0':0.25,'x1':0.25,'y_initial':3.,'problem':'failure'})
    assert r['status']=='ok' and r['rhs_calls']==0 and r['y1']==3.

def test_step_floor_returns_prefix(numerics):
    r=numerics('ode',{'method':'adaptive','h_initial':1.,'h_min':1.,'h_max':1.,'rtol':1e-14,'atol1':1e-15})
    assert r['status']=='step_too_small' and r['accepted']==0 and r['last_x']==0.
    assert r['rejected']>=1 and not r['result_valid']

def test_max_attempt_budget(numerics):
    r=numerics('ode',{'method':'adaptive','h_initial':0.01,'max_steps':1})
    assert r['status']=='max_steps' and r['attempts']==1 and r['last_x']<1.
    assert not r['result_valid']

def test_absolute_only_tolerances(numerics):
    r=numerics('ode',{'method':'adaptive','rtol':0.,'atol1':1e-10})
    assert r['status']=='ok' and abs(r['y1']-math.e)<1e-8

def test_tighter_tolerances_improve_analytic_fixture(numerics):
    e=[]
    for tol in (1e-4,1e-7,1e-10):
        r=numerics('ode',{'method':'adaptive','h_initial':1.,'rtol':tol,'atol1':tol*0.01})
        assert r['status']=='ok';e.append(abs(r['y1']-math.e))
    assert e[2]<e[1]<e[0]

def test_exact_stop_below_hmin_allowed_but_not_failed_trial(numerics):
    r=numerics('ode',{'method':'adaptive','problem':'constant','x1':0.03,
                     'h_initial':0.5,'h_min':0.1,'h_max':1.})
    assert r['status']=='ok' and r['last_x']==0.03

def test_step_stagnation_detected(numerics):
    r=numerics('ode',{'method':'adaptive','x0':1e20,'x1':1e20+1e6,'h_initial':1.,'h_max':10.})
    assert r['status']=='stagnation' and r['rhs_calls']==0

def test_endpoint_not_fabricated_by_roundoff(numerics):
    r=numerics('ode',{'method':'adaptive','problem':'constant','x0':-1e16,'x1':1.,
                     'h_initial':1e16,'h_max':1e16,'y_initial':0.})
    assert r['status']=='stagnation' and not r['result_valid']
    assert r['rhs_calls']==0

def test_unrecorded_trial_audit_is_not_claimed(numerics):
    r=numerics('ode',{'method':'adaptive','record_trials':False})
    assert r['status']=='ok' and not r['trial_audit_available']
    assert 'rejection_rollback_ok' not in r

def test_history_growth_beyond_initial_capacity(numerics):
    r=numerics('ode',{'method':'adaptive','h_max':0.005})
    assert r['status']=='ok' and r['accepted']>=200
    assert r['profile_points']==r['accepted']+1
    assert r['y1']==pytest.approx(math.e,abs=1e-12)

def test_failure_retains_only_accepted_prefix(numerics):
    r=numerics('ode',{'method':'adaptive','problem':'failure_after','h_max':0.1})
    assert r['status']=='out_of_domain' and r['accepted']>0
    assert r['last_x']<=0.3 and not r['result_valid'] and 'y1' not in r
    assert r['last_accepted_y1']==pytest.approx(math.exp(r['last_x']),rel=2e-9)

