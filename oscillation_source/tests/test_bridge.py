import math
from pathlib import Path
import pytest

MU_T=1.5e-10


def test_branch_based_hermite_bridge_matches_endpoint_pressure_and_density(run_case):
    r=run_case('mu_bridge',{'scenario':'maxwell_exact','mu_low_j':0.95*MU_T,'mu_high_j':1.05*MU_T,
                            'mu_query_j':MU_T})
    assert r['status']=='ok'
    assert r['bridge_ready'] is True
    assert r['endpoint_pressure_closed'] is True
    assert r['endpoint_density_closed'] is True
    assert r['stable_interval'] is True
    assert r['causal_interval'] is True
    assert r['admitted'] is True
    assert r['dn_dmu_m3_per_j'] > 0
    assert 0 < r['chi'] <= 1
    assert r['energy_j_m3']==pytest.approx(r['mu_query_j']*r['n_b_m3']-r['pressure_pa'],rel=3e-13)


def test_bridge_endpoint_evaluations_recover_branches(run_case):
    for which,mu in [('low',0.95*MU_T),('high',1.05*MU_T)]:
        r=run_case('mu_bridge',{'scenario':'maxwell_exact','mu_low_j':0.95*MU_T,'mu_high_j':1.05*MU_T,
                                'mu_query_j':mu})
        assert r['status']=='ok'
        pscale=r[f'bridge_{which}_pressure_pa']
        assert abs(r[f'{which}_pressure_error_pa']) <= 2e-15*pscale
        assert abs(r[f'{which}_density_relative_error']) <= 2e-15


def test_bridge_query_outside_interval_is_not_extrapolated(run_case):
    r=run_case('mu_bridge',{'scenario':'maxwell_exact','mu_low_j':0.95*MU_T,'mu_high_j':1.05*MU_T,
                            'mu_query_j':0.94*MU_T})
    assert r['status']=='out_of_domain'
    assert 'pressure_pa' not in r


def test_bridge_rejects_negative_compressibility_over_interval(run_case):
    r=run_case('mu_bridge',{'scenario':'negative_compressibility'})
    assert r['status']=='matching_inconsistent'
    assert r['bridge_ready'] is False
    assert r['stable_interval'] is False


def test_bridge_rejects_superluminal_interval(run_case):
    r=run_case('mu_bridge',{'scenario':'acausal'})
    assert r['status']=='matching_inconsistent'
    assert r['bridge_ready'] is False
    assert r['causal_interval'] is False


def test_bridge_checks_interval_not_only_midpoint_samples(run_case):
    r=run_case('mu_bridge',{'scenario':'endpoint_hidden_failure'})
    assert r['status']=='matching_inconsistent'
    assert r['bridge_ready'] is False
    assert r['interval_extrema_checked'] is True


def test_bridge_requires_positive_ordered_mu_interval(run_case):
    r=run_case('mu_bridge',{'scenario':'maxwell_exact','mu_low_j':1.05*MU_T,'mu_high_j':0.95*MU_T})
    assert r['status']=='invalid_input'


def test_bridge_nonfinite_query_is_rejected(run_case):
    r=run_case('mu_bridge',{'scenario':'maxwell_exact','special_input':'nan_query'})
    assert r['status']=='nonfinite_input'
