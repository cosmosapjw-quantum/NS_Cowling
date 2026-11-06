import math
import pytest

MU_T = 1.5e-10
P_T = 1.0e33
N_LOW = P_T/1.5e-11
N_HIGH = 2.0*P_T/1.5e-11


def test_continuous_two_polytrope_knot_preserves_all_state_variables(run_case):
    r=run_case('two_poly_knot',{
        'core_n_anchor_m3':1e44,'core_p_anchor_pa':1e32,'core_gamma':2.5,'core_c_j':1.5e-10,
        'core_n_min_m3':1e44,'core_n_max_m3':5e44,
        'transition_n_m3':1e44,'crust_n_min_m3':1e40,'gamma_crust':1.5,
    })
    assert r['status']=='ok'
    assert r['core_unchanged'] is True
    assert r['pressure_closed'] is True
    assert r['energy_closed'] is True
    assert r['chemical_closed'] is True
    assert r['density_closed'] is True
    assert r['chi_jump_nonzero'] is True
    assert r['core_pressure_pa']==pytest.approx(1e32,rel=3e-14)


def test_continuous_knot_bad_transition_is_rejected_without_refitting_core(run_case):
    r=run_case('two_poly_knot',{'transition_n_m3':5e43})
    assert r['status']=='out_of_domain'
    assert r['core_unchanged'] is True
    assert r['model_ready'] is False


def test_direct_join_pressure_match_can_expose_chemical_mismatch(run_case):
    r=run_case('direct_join',{'scenario':'maxwell_exact','pressure_pa':8e32})
    assert r['status']=='ok'
    assert r['pressure_residual_pa']==0.0
    assert abs(r['mu_mismatch_j']) > 1e-13
    assert r['chemical_closed'] is False
    assert r['thermodynamically_consistent'] is False


def test_direct_join_at_exact_crossing_is_not_declared_inconsistent(run_case):
    r=run_case('direct_join',{'scenario':'maxwell_exact','pressure_pa':P_T,'residual_tol':1e-12})
    assert r['status']=='ok'
    assert abs(r['mu_mismatch_j']) <= 1e-22
    assert r['chemical_closed'] is True
    assert r['thermodynamically_consistent'] is True


def test_maxwell_exact_fixture_recovers_known_crossing_and_jump(run_case):
    r=run_case('maxwell_match',{'scenario':'maxwell_exact','mu_left_j':1.46e-10,'mu_right_j':1.54e-10,
                                'root_xtol':1e-24,'root_rtol':0.0,'residual_tol':1e17,'max_iter':200})
    assert r['status']=='ok'
    assert r['mu_t_j']==pytest.approx(MU_T,rel=3e-12)
    assert r['p_t_pa']==pytest.approx(P_T,rel=3e-12)
    assert r['n_low_m3']==pytest.approx(N_LOW,rel=3e-12)
    assert r['n_high_m3']==pytest.approx(N_HIGH,rel=3e-12)
    assert r['n_high_m3'] > r['n_low_m3']
    assert r['branch_selection_valid'] is True
    assert r['thermodynamically_consistent'] is True
    assert r['delta_e_j_m3']==pytest.approx(r['mu_t_j']*r['delta_n_m3'],rel=5e-12)
    assert abs(r['pressure_residual_pa']) <= 1e18
    assert abs(r['mu_mismatch_j']) <= 1e-22


def test_maxwell_requires_a_common_mu_domain(run_case):
    r=run_case('maxwell_match',{'scenario':'maxwell_exact','mu_left_j':1.0e-10,'mu_right_j':1.1e-10})
    assert r['status']=='out_of_domain'
    assert r['transition_valid'] is False


def test_maxwell_missing_crossing_stays_no_bracket(run_case):
    r=run_case('maxwell_match',{'scenario':'no_crossing','mu_left_j':1.46e-10,'mu_right_j':1.54e-10})
    assert r['status']=='no_bracket'
    assert r['transition_valid'] is False


def test_maxwell_reversed_phase_roles_are_not_admitted(run_case):
    r=run_case('maxwell_match',{'scenario':'maxwell_reversed','mu_left_j':1.46e-10,'mu_right_j':1.54e-10})
    assert r['status']=='matching_inconsistent'
    assert r['transition_valid'] is False
    assert r['branch_selection_valid'] is False


def test_maxwell_nonfinite_input_is_rejected(run_case):
    r=run_case('maxwell_match',{'scenario':'maxwell_exact','special_input':'nan_mu_left'})
    assert r['status']=='nonfinite_input'


def test_matching_run_contract_records_thermodynamic_scope(run_case):
    import json
    from pathlib import Path
    r=run_case('maxwell_match',{'scenario':'maxwell_exact','mu_left_j':1.46e-10,'mu_right_j':1.54e-10})
    c=json.loads(Path(r['_evidence']['run_receipt']).read_text())['contract']
    assert c['units']=='SI'
    assert c['model']=='analytic_two_phase_polytropes'
    assert c['closure']=='cold_zero_temperature_thermodynamics'
    assert c['interface']=='equilibrium_matching_only_no_mode_transfer'


def test_maxwell_exact_fixture_documents_independent_analytic_oracle():
    import json
    from pathlib import Path
    path=Path('fixtures/analytic/maxwell_exact.json')
    assert path.is_file()
    f=json.loads(path.read_text())
    assert f['schema_version']==1
    assert f['gamma']==2.0
    assert float(f['mu_t_j'])==pytest.approx(MU_T,rel=1e-15)
    assert float(f['p_t_pa'])==pytest.approx(P_T,rel=1e-15)
    assert float(f['n_low_m3'])==pytest.approx(N_LOW,rel=1e-15)
    assert float(f['n_high_m3'])==pytest.approx(N_HIGH,rel=1e-15)
    # Recompute the two P(mu) branches without calling product code.
    mu=MU_T
    for side,n_key,k_key,c_key in [('low','n_low_m3','k_low_pa_m6','c_low_j'),
                                    ('high','n_high_m3','k_high_pa_m6','c_high_j')]:
        n=float(f[n_key]); k=float(f[k_key]); c=float(f[c_key])
        n_from_mu=(mu-c)/(2.0*k)
        p=k*n_from_mu*n_from_mu
        assert n_from_mu==pytest.approx(n,rel=3e-15)
        assert p==pytest.approx(P_T,rel=3e-15)
    assert f['derivation']['delta_e_identity']=='Delta epsilon = mu_t Delta n'
