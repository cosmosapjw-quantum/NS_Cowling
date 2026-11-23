#!/usr/bin/env python3
"""Cross-check the native constrained profiles and the uniform beta control."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def rows(path):
    with path.open() as stream:
        records = list(csv.DictReader(stream))
    return {key: np.array([float(row[key]) for row in records]) for key in records[0]}


def check(condition, message):
    if not condition:
        raise ValueError(message)


def compose(directory, tab, has_beta):
    n = np.loadtxt(directory / 'eos.nb', skiprows=2)
    thermo = np.loadtxt(directory / 'eos.thermo', skiprows=1)
    mn = float((directory / 'eos.thermo').read_text().split()[0])
    composition = np.loadtxt(directory / 'eos.compo')
    check(np.allclose(n, tab['nb_fm3'], rtol=1e-14, atol=0), 'density mismatch')
    p = n * thermo[:, 3]
    energy = n * mn * (thermo[:, 9] + 1)
    check(np.allclose(p, tab['pressure_MeV_fm3'], rtol=1e-13), 'pressure mismatch')
    check(np.allclose(energy, tab['epsilon_MeV_fm3'], rtol=1e-13), 'energy mismatch')
    bsq = np.diff(p) / np.diff(energy)
    check(np.all((bsq > 0) & (bsq < 1)), 'secant stability/causality failure')
    check(np.all(np.diff(n) > 0) and np.all(np.diff(energy) > 0), 'nonmonotone input')
    beta_inferred = bool(np.max(abs(thermo[:, 7])) < 1e-8)
    check(beta_inferred == has_beta, 'composition-response beta mismatch')
    check(('beta_equilibrated=' + ('T' if has_beta else 'F')) in
          (directory / 'metadata.nml').read_text(), 'metadata beta mismatch')
    if has_beta:
        ye = tab['Ye']
        check(np.allclose(composition[:, 6] + composition[:, 8], 1), 'nucleon closure')
        check(np.allclose(composition[:, 8], composition[:, 10]), 'neutrality')
    else:
        ye = tab['Z'] / (tab['N'] + tab['Z'])
        check(np.allclose(composition[:, 6], ye), 'cell/electron charge mismatch')
        check(np.array_equal(composition[:, 3], tab['profile_sector']), 'sector label mismatch')
    euler = np.max(np.abs((energy + p) / n - mn * (1 + thermo[:, 5] + ye * thermo[:, 7])))
    check(euler < 1e-8, 'generalized Euler identity failed')
    return {'rows': int(n.size), 'n_min_fm-3': n[0], 'n_max_fm-3': n[-1],
            'cs2_over_c2_secant_range': [bsq.min(), bsq.max()],
            'generalized_euler_max_abs_MeV': euler}


def run(primary, refined):
    report = {'fine_radial_mesh': 256, 'published_control_radial_mesh': 128}
    for mode in ('printed', 'd03x100'):
        tab, fine = (rows(base / f'{mode}.csv') for base in (primary, refined))
        check(np.array_equal(tab['nb_fm3'], fine['nb_fm3']), 'mesh test density mismatch')
        stats = compose(primary / mode, tab, False)
        sector = tab['profile_sector']
        cs = tab['max_nucleon_density_fm3']
        check(np.all(cs <= .15), 'DME density-domain violation')
        check(set(sector) == {1, 2} and np.sum(np.diff(sector) != 0) == 1, 'shape partition')
        check(np.all((tab['N'] == 950) & (tab['Z'] == 32)), 'variable isotope on compression')
        anchors = rows(primary / f'{mode}_anchors.csv')
        reference_error = anchors['E_minus_mn_MeV'] - anchors['source_ground_E_minus_mn_MeV']
        fine_dp = np.max(abs(tab['pressure_MeV_fm3'] - fine['pressure_MeV_fm3']))
        check(fine_dp < 5e-5, 'radial refinement failed')
        check(max(tab['dilation_P_rel_error']) < 2e-7, 'fixed-shape pressure derivative failed')
        check(max(tab['mu_L_finite_difference_error_MeV']) < 2e-5, 'composition derivative failed')
        stats.update({
            'clustered_rows': int(sum(sector == 1)),
            'last_clustered_fm-3': max(tab['nb_fm3'][sector == 1]),
            'first_homogeneous_fm-3': min(tab['nb_fm3'][sector == 2]),
            'max_density_any_species_fm-3': max(cs),
            'mu_L_range_MeV': [min(tab['mu_L_MeV']), max(tab['mu_L_MeV'])],
            'max_particle_number_residual': max(tab['max_particle_error']),
            'max_dilation_pressure_relative_error': max(tab['dilation_P_rel_error']),
            'max_FD_charge_derivative_error_MeV': max(tab['mu_L_finite_difference_error_MeV']),
            'mesh_128_vs_256_max_abs_pressure_MeV_fm-3': fine_dp,
            'mesh_128_vs_256_max_abs_E_per_baryon_MeV': max(abs(tab['E_minus_mn_MeV'] - fine['E_minus_mn_MeV'])),
            'anchor_E_per_baryon_minus_table3_range_MeV': [min(reference_error), max(reference_error)],
            'negative_pressure_anchor_rows': int(sum(anchors['pressure_MeV_fm3'] < 0)),
        })
        report[mode] = stats
    uniform = rows(primary / 'uniform_beta.csv')
    stats = compose(primary / 'uniform_beta', uniform, True)
    stats.update({
        'Ye_range': [min(uniform['Ye']), max(uniform['Ye'])],
        'max_abs_beta_residual_MeV': max(abs(uniform['beta_residual_MeV'])),
        'max_first_law_P_rel_error': max(uniform['first_law_P_rel_error']),
        'analytic_charge_curvature_min_sampled_MeV_fm3': min(uniform['minimum_sampled_charge_curvature_MeV_fm3']),
        'cs2_over_c2_local_range': [min(uniform['cs2_over_c2_FD']), max(uniform['cs2_over_c2_FD'])],
    })
    report['uniform_beta'] = stats
    return report


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('primary', type=Path)
    cli.add_argument('refined', type=Path)
    args = cli.parse_args()
    print(json.dumps(run(args.primary, args.refined), indent=2))
