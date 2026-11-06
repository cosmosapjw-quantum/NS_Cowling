from tests.test_tov import star

def test_target_mass_first(star):
    r=star({},'target_mass')
    assert r['status']=='ok',r
    assert r['target_converged'] is True
import math
import pytest

def test_same_mass_changes_central_pressure_not_target(star):
    a=star({},'target_mass');b=star({'p_anchor_pa':1.05e32},'target_mass')
    assert a['status']==b['status']=='ok'
    assert a['pc_pa']!=b['pc_pa']
    for r in (a,b):
        assert abs(r['mass_kg']-1.2e30)<=1e20+1.2e22
        assert r['sampled_branch_increasing'] and r['star_evaluations']>9
    assert a['radius_m']!=b['radius_m']

@pytest.mark.parametrize('pars',[
    {'m_target_kg':0.},{'m_target_kg':-1.},{'m_target_kg':1e31},
    {'pc_lo_pa':1e34,'pc_hi_pa':1e33},{'pc_lo_pa':-1.},{'pc_hi_pa':1e40},
    {'scan_n':2},{'scan_n':66},{'mass_atol_kg':-1.},{'mass_rtol':-1.},
    {'root_max_iter':1,'mass_rtol':1e-14},{'special_input':'fail_callback'},
    {'root_xtol':1.,'mass_rtol':1e-14},
])
def test_invalid_target_or_unresolved_mass_is_rejected(star,pars):
    r=star(pars,'target_mass');assert r['status']!='ok',r
    assert not r['target_converged'] and 'mass_kg' not in r

def test_turning_region_not_treated_as_increasing_branch(star):
    r=star({'pc_lo_pa':1e35,'pc_hi_pa':1e37,'n_max_m3':1e48},'target_mass')
    assert r['status']=='out_of_domain',r
    assert 'branch' in r['detail']
    assert not r['sampled_branch_increasing']
