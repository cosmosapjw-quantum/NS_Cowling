import math
import pytest


def test_material_interface_and_reverse(run_case):
    r=run_case('material_transfer',{})
    assert r['status']=='ok'
    assert r['w_out_m3']==pytest.approx(r['w_in_m3'],rel=2e-14)
    assert r['delta_p_out_pa']==pytest.approx(r['delta_p_in_pa'],rel=2e-13)
    assert r['v_reverse_m2']==pytest.approx(r['v_in_m2'],rel=2e-13)
    assert r['zero_jump_identity'] is True
    assert r['same_radius_two_sides'] is True
    assert r['v_continuous_delta_p_mismatch_pa'] != pytest.approx(0.0,abs=1e-6)


def test_material_interface_rejects_zero_frequency(run_case):
    r=run_case('material_transfer',{'omega_rad_s':0.0})
    assert r['status']=='invalid_input'
    assert r['result_valid'] is False


def test_material_interface_rejects_wrong_side_orientation(run_case):
    r=run_case('material_transfer',{'inner_side':-1,'outer_side':1})
    assert r['status']=='invalid_input'


def test_segment_events_are_pressure_sorted_and_exact_stops(run_case):
    r=run_case('segment_events',{})
    assert r['status']=='ok'
    assert r['events_valid'] is True
    assert r['q2'] > r['q1'] > 0.0
    assert r['cross_stage_allowed'] is False
    assert r['material_side_changes'] is True
    assert r['chi_knot_state_continuous'] is True


def test_unsorted_segment_events_are_rejected(run_case):
    r=run_case('segment_events',{'scenario':'unsorted'})
    assert r['status']=='invalid_input'


def test_material_event_cannot_reuse_same_side(run_case):
    r=run_case('segment_events',{'scenario':'same_side_material'})
    assert r['status']=='invalid_input'


def test_chi_knot_requires_same_state_side(run_case):
    r=run_case('segment_events',{'scenario':'chi_side_change'})
    assert r['status']=='invalid_input'


def test_segment_restart_hits_events_without_crossing_and_preserves_one_sided_state(run_case):
    r=run_case('segment_restart_probe',{})
    assert r['status']=='ok'
    assert r['segments_completed']==3
    assert r['material_event_hit'] is True
    assert r['chi_event_hit'] is True
    assert r['any_segment_crossed_event'] is False
    assert r['material_q_inner']==pytest.approx(r['material_q_outer'],abs=0.0)
    assert r['material_r_inner_m']==pytest.approx(r['material_r_outer_m'],abs=0.0)
    assert r['material_side_inner']==1
    assert r['material_side_outer']==-1
    assert r['material_w_outer_m3']==pytest.approx(r['material_w_inner_m3'],rel=2e-14)
    assert r['material_delta_p_outer_pa']==pytest.approx(r['material_delta_p_inner_pa'],rel=2e-13)
    assert abs(r['material_v_outer_m2']-r['material_v_inner_m2']) > 1e-6
    assert r['chi_q_before']==pytest.approx(r['chi_q_after'],abs=0.0)
    assert r['chi_w_after_m3']==pytest.approx(r['chi_w_before_m3'],rel=2e-14)
    assert r['chi_v_after_m2']==pytest.approx(r['chi_v_before_m2'],rel=2e-14)


def test_material_transfer_matches_independent_enthalpy_formula(run_case):
    pars={'r_m':8000.0,'omega_rad_s':2*math.pi*3000.0,'pressure_pa':1e33,
          'eps_in':1.1e34,'eps_out':7e33,'nu':-0.2,'lambda_value':0.15,
          'dnu_dr':1e-4,'w_in_m3':2e11,'v_in_m2':-3e7}
    r=run_case('material_transfer',pars)
    assert r['status']=='ok'
    c=299792458.0
    a=(pars['omega_rad_s']/c*pars['r_m'])**2*math.exp(pars['lambda_value']-2*pars['nu'])
    q=(pars['eps_in']+pars['pressure_pa'])/(pars['eps_out']+pars['pressure_pa'])
    expected=q*pars['v_in_m2']+(q-1)*pars['dnu_dr']*pars['w_in_m3']/a
    assert r['v_out_m2']==pytest.approx(expected,rel=3e-13)
