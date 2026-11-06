import os
from pathlib import Path
import pytest
from tools.runtime import CaseRunner

ROOT = Path(__file__).resolve().parents[1]

@pytest.fixture(scope='module')
def liquid_drop(tmp_path_factory):
    runner = CaseRunner(ROOT, tmp_path_factory.mktemp('liquid-drop'), workspace=os.environ.get('NS_WORKSPACE', 'learner' if (ROOT/'WORKSPACE.json').is_file() else 'solutions'))
    artifact = runner.build_fpm()
    return lambda case, pars={}: runner.run_artifact(artifact, case, pars)


def test_liquid_drop_scaling_and_domain(liquid_drop):
    r = liquid_drop('liquid_drop_scaling', {})
    assert r['status'] == 'ok'
    assert r['surface_a8_over_a'] == pytest.approx(4.0, rel=2e-14)
    assert r['coulomb_fixed_z_a8_over_a'] == pytest.approx(0.5, rel=2e-14)
    assert r['invalid_charge_rejected'] is True
    assert r['full_star_provider_enabled'] is False


def test_candidate_decomposition_and_charge_neutrality(liquid_drop):
    a, z, n_b = 56, 26, 1.0e-8
    r = liquid_drop('liquid_drop_candidate', {'a': a, 'z': z, 'n_b_fm3': n_b})
    assert r['status'] == 'ok'
    assert r['candidate_valid'] is True
    assert r['a'] == a and r['z'] == z
    assert r['volume_mev_per_baryon'] == pytest.approx(-15.75, rel=2e-14)
    assert r['surface_mev_per_baryon'] == pytest.approx(17.8*a**(-1/3), rel=3e-13)
    assert r['coulomb_mev_per_baryon'] == pytest.approx(0.711*z*z*a**(-4/3), rel=3e-13)
    assert r['asymmetry_mev_per_baryon'] == pytest.approx(23.7*(a-2*z)**2/a**2, rel=3e-13)
    assert r['n_e_m3'] == pytest.approx(n_b*1e45*z/a, rel=2e-14)
    assert r['charge_neutrality_closed'] is True
    nuclear = (r['volume_mev_per_baryon'] + r['surface_mev_per_baryon'] +
               r['coulomb_mev_per_baryon'] + r['asymmetry_mev_per_baryon'])
    assert r['nuclear_mev_per_baryon'] == pytest.approx(nuclear, rel=3e-13)
    assert r['electron_mev_per_baryon'] > 0.0
    assert r['total_mev_per_baryon'] == pytest.approx(nuclear+r['electron_mev_per_baryon'], rel=3e-13)
    assert r['pairing_included'] is False
    assert r['lattice_energy_included'] is False


@pytest.mark.parametrize('a,z', [(0, 0), (10, -1), (10, 11)])
def test_invalid_candidate_charge_or_mass_is_rejected(liquid_drop, a, z):
    r = liquid_drop('liquid_drop_candidate', {'a': a, 'z': z, 'n_b_fm3': 1e-8})
    assert r['status'] == 'invalid_input'
    assert r.get('candidate_valid', False) is False
    assert 'A>0' in r.get('detail', '') and '0<=Z<=A' in r.get('detail', '')


def test_toy_density_domain_stops_before_neutron_drip_physics(liquid_drop):
    high = liquid_drop('liquid_drop_candidate', {'a': 100, 'z': 35, 'n_b_fm3': 1e-2})
    low = liquid_drop('liquid_drop_candidate', {'a': 56, 'z': 26, 'n_b_fm3': 1e-12})
    assert high['status'] == 'out_of_domain'
    assert low['status'] == 'out_of_domain'
    assert high['full_star_provider_enabled'] is False
    assert high['dripped_neutrons_included'] is False


def test_discrete_scan_minimum_beats_adjacent_candidates(liquid_drop):
    pars = {'n_b_fm3': 1e-4, 'a_min': 20, 'a_max': 180, 'a_step': 1,
            'z_min': 1, 'z_max': 60, 'z_step': 1, 'objective': 'fixed_density_energy'}
    r = liquid_drop('liquid_drop_scan', pars)
    assert r['status'] == 'ok'
    assert r['minimum_kind'] == 'discrete_candidate_set'
    assert r['stationary_point_claimed'] is False
    assert r['candidate_count'] > 0
    assert pars['a_min'] <= r['best_a'] <= pars['a_max']
    assert pars['z_min'] <= r['best_z'] <= min(pars['z_max'], r['best_a'])
    assert r['charge_neutrality_closed'] is True
    for da, dz in ((-1,0),(1,0),(0,-1),(0,1)):
        aa, zz = r['best_a']+da, r['best_z']+dz
        if aa >= pars['a_min'] and aa <= pars['a_max'] and zz >= pars['z_min'] and zz <= min(pars['z_max'], aa):
            q = liquid_drop('liquid_drop_candidate', {'a': aa, 'z': zz, 'n_b_fm3': pars['n_b_fm3']})
            assert q['status'] == 'ok'
            assert r['best_total_mev_per_baryon'] <= q['total_mev_per_baryon'] + 2e-12


def test_finer_candidate_mesh_cannot_raise_minimum(liquid_drop):
    common = {'n_b_fm3': 1e-4, 'a_min': 20, 'a_max': 180, 'z_min': 2, 'z_max': 60,
              'objective': 'fixed_density_energy'}
    coarse = liquid_drop('liquid_drop_scan', {**common, 'a_step': 4, 'z_step': 2})
    fine = liquid_drop('liquid_drop_scan', {**common, 'a_step': 1, 'z_step': 1})
    assert coarse['status'] == fine['status'] == 'ok'
    assert fine['candidate_count'] > coarse['candidate_count']
    assert fine['best_total_mev_per_baryon'] <= coarse['best_total_mev_per_baryon'] + 2e-12
    assert fine['mesh_is_discrete'] is True and coarse['mesh_is_discrete'] is True


def test_boundary_minimum_is_not_called_stationary_point(liquid_drop):
    r = liquid_drop('liquid_drop_scan', {'n_b_fm3': 1e-8, 'a_min': 56, 'a_max': 56,
                                         'a_step': 1, 'z_min': 26, 'z_max': 26, 'z_step': 1,
                                         'objective': 'fixed_density_energy'})
    assert r['status'] == 'ok'
    assert r['minimum_on_boundary'] is True
    assert r['stationary_point_claimed'] is False


def test_fixed_pressure_gibbs_is_explicitly_out_of_scope(liquid_drop):
    r = liquid_drop('liquid_drop_scan', {'n_b_fm3': 1e-8, 'objective': 'fixed_pressure_gibbs'})
    assert r['status'] == 'out_of_domain'
    assert r['objective_supported'] is False
    assert r['fixed_density_only'] is True


def _electron_energy_per_baryon_mev(n_b_fm3, a, z):
    import math
    hbar = 6.62607015e-34/(2*math.pi)
    c = 299792458.0
    m_e = 9.1093837139e-31
    mev_j = 1.602176634e-13
    n_b = n_b_fm3*1e45
    n_e = n_b*z/a
    if n_e == 0.0:
        return 0.0
    pf = hbar*(3*math.pi**2*n_e)**(1/3)
    x = pf/(m_e*c)
    if x < 5e-2:
        e = n_e*m_e*c*c*(1+3*x*x/10-3*x**4/56+x**6/48-15*x**8/1408)
    else:
        pref = m_e**4*c**5/(8*math.pi**2*hbar**3)
        e = pref*(x*(1+2*x*x)*math.sqrt(1+x*x)-math.asinh(x))
    return e/n_b/mev_j


def _independent_total_mev_per_baryon(n_b_fm3, a, z):
    return (-15.75 + 17.8*a**(-1/3) + 0.711*z*z*a**(-4/3)
            + 23.7*(a-2*z)**2/a**2 + _electron_energy_per_baryon_mev(n_b_fm3, a, z))


def test_scan_matches_independent_small_grid_bruteforce(liquid_drop):
    pars = {'n_b_fm3': 1e-8, 'a_min': 40, 'a_max': 70, 'a_step': 1,
            'z_min': 15, 'z_max': 35, 'z_step': 1, 'objective': 'fixed_density_energy'}
    expected = min((_independent_total_mev_per_baryon(pars['n_b_fm3'], a, z), a, z)
                   for a in range(pars['a_min'], pars['a_max']+1)
                   for z in range(pars['z_min'], min(pars['z_max'], a)+1))
    r = liquid_drop('liquid_drop_scan', pars)
    assert r['status'] == 'ok'
    assert (r['best_a'], r['best_z']) == expected[1:]
    assert r['best_total_mev_per_baryon'] == pytest.approx(expected[0], rel=2e-12)


def test_liquid_drop_fixture_locks_empirical_inputs_and_scope(liquid_drop):
    import json
    fixture = json.loads((ROOT/'fixtures/nuclear/liquid_drop_toy.json').read_text())
    r = liquid_drop('liquid_drop_scaling', {})
    c = fixture['coefficient_set']
    assert r['a_volume_mev'] == pytest.approx(c['a_volume_mev'], rel=0, abs=0)
    assert r['a_surface_mev'] == pytest.approx(c['a_surface_mev'], rel=0, abs=0)
    assert r['a_coulomb_mev'] == pytest.approx(c['a_coulomb_mev'], rel=0, abs=0)
    assert r['a_asymmetry_mev'] == pytest.approx(c['a_asymmetry_mev'], rel=0, abs=0)
    assert r['n_b_min_fm3'] == pytest.approx(fixture['density_admission_fm3']['min'], rel=0, abs=0)
    assert r['n_b_max_fm3'] == pytest.approx(fixture['density_admission_fm3']['max'], rel=0, abs=0)
    assert fixture['conventions']['pairing_term_included'] is False
    assert 'realistic inner-crust EOS' in fixture['not_included']
