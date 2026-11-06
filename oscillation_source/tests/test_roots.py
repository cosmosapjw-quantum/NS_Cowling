"""Analytic contracts for the locally implemented bisection kernel."""
import math
import os
import pytest
from pathlib import Path
from tools.runtime import CaseRunner

@pytest.fixture(scope='module')
def numerics(tmp_path_factory):
    root=Path(__file__).resolve().parents[1]
    runner=CaseRunner(root,tmp_path_factory.mktemp('root-numerics'),workspace=os.environ.get('NS_WORKSPACE','solutions'))
    artifact=runner.build_fpm()
    return lambda case,pars:runner.run_artifact(artifact,case,pars)

def test_root_sqrt2_first_success(numerics):
    r=numerics('root_sqrt2',{})
    assert r['status']=='ok'
    assert abs(r['root']-math.sqrt(2.0))<2e-12
    assert r['lo_final']<=math.sqrt(2.0)<=r['hi_final']

@pytest.mark.parametrize('a,b,expected',[(0.,2.,0.),(-2.,0.,0.),(0.,0.,0.)])
def test_endpoint_root(numerics,a,b,expected):
    r=numerics('root_probe',{'lo':a,'hi':b,'root_at':expected})
    assert r['status']=='ok' and r['root']==expected
    assert r['termination']=='exact'

def test_no_bracket_not_a_number(numerics):
    r=numerics('root_sqrt2',{'lo':2.,'hi':3.})
    assert r['status']=='no_bracket'
    assert r['result_valid'] is False and 'root' not in r

@pytest.mark.parametrize('params',[{'lo':3.,'hi':1.},{'xtol':-1.},{'rtol':-1.},{'ftol':-1.},{'max_iter':0}])
def test_invalid_options(numerics,params):
    r=numerics('root_sqrt2',params)
    assert r['status']=='invalid_input' and r['evaluations']==0

@pytest.mark.parametrize('scale',[1e300,1e-300])
def test_sign_without_multiplication(numerics,scale):
    r=numerics('root_probe',{'lo':0.,'hi':3.,'root_at':1.23456789,'scale':scale})
    assert r['status']=='ok'
    assert abs(r['root']-1.23456789)<2e-12

def test_midpoint_no_large_bracket_overflow(numerics):
    r=numerics('root_probe',{'problem':'wide','lo':-1e308,'hi':1e308,'xtol':1e292,'max_iter':180})
    assert r['status']=='ok'
    assert r['root']==pytest.approx(1e307,rel=4e-14)

def test_adjacent_bracket_stagnation(numerics):
    r=numerics('root_probe',{'problem':'adjacent','xtol':0.,'rtol':0.,'ftol':0.})
    assert r['status']=='stagnation' and r['result_valid'] is False
    assert r['lo_final']==1. and r['hi_final']==math.nextafter(1.,math.inf)

@pytest.mark.parametrize('problem,status',[('nan','nonfinite_result'),('failure','out_of_domain')])
def test_callback_failure(numerics,problem,status):
    r=numerics('root_probe',{'problem':problem})
    assert r['status']==status and not r['result_valid']

def test_midpoint_callback_failure_not_connected(numerics):
    r=numerics('root_probe',{'problem':'failure_mid','lo':0.,'hi':2.})
    assert r['status']=='out_of_domain' and r['evaluations']==3

def test_max_iter_preserves_real_bracket(numerics):
    r=numerics('root_sqrt2',{'max_iter':1,'xtol':0.,'rtol':0.,'ftol':0.})
    assert r['status']=='max_iterations' and not r['result_valid']
    assert r['lo_final']<=math.sqrt(2)<=r['hi_final']
    assert r['iterations']==1

def test_refinement_and_termination_are_separate(numerics):
    a=numerics('root_sqrt2',{'xtol':1e-5})
    b=numerics('root_sqrt2',{'xtol':1e-11})
    assert a['status']==b['status']=='ok'
    assert b['hi_final']-b['lo_final'] < a['hi_final']-a['lo_final']
    r=numerics('root_probe',{'lo':0.,'hi':2.,'root_at':0.3,'scale':1e-20,'xtol':0.,'ftol':1e-19})
    assert r['termination']=='residual'
    assert abs(r['root']-0.3)>0.1 # small residual is not a coordinate bound

def test_nested_bisection_has_independent_state(numerics):
    r=numerics('root_probe',{'problem':'nested','lo':1.,'hi':2.})
    assert r['status']=='ok'
    assert r['root']==pytest.approx(math.sqrt(2),abs=3e-12)
