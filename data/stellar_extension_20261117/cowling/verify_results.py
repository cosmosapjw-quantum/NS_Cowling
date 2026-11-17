#!/usr/bin/env python3
"""Read back Stage E output identities, real-root receipts and serialized profiles."""
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
LIMIT = 2e-6


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def count_nodes(values):
    values = np.asarray(values, dtype=float)
    kept = values[np.abs(values) > 1e-9*np.max(np.abs(values))]
    return int(np.count_nonzero(kept[:-1]*kept[1:] < 0))


def main():
    result = json.loads((OUT/'SUMMARY.json').read_text())
    provenance = json.loads((OUT/'PROVENANCE.json').read_text())
    execution = json.loads((OUT/'EXECUTION.json').read_text())
    assert execution['status'] == 'FINISHED' and execution['exit_code'] == 0
    assert result['residual_acceptance'] == LIMIT
    identities = [provenance[key] for key in ['runner','primary_contract','comparison_contract','selection'] if provenance.get(key)]
    identities += provenance['reused_D_controls']
    for record in identities:
        assert sha(ROOT/record['path']) == record['sha256'], record['path']
    rows = list(result['configurations'])
    extra = json.loads((OUT/'TIGHT_BACKGROUND_COMPARISON.json').read_text())
    extra_execution = json.loads((OUT/'EXECUTION_TIGHT_BACKGROUND.json').read_text())
    assert extra_execution['status'] == 'FINISHED' and extra_execution['exit_code'] == 0
    assert sha(OUT/'tight_background_check.py') == extra['script_sha256']
    assert sha(ROOT/extra['reference_path']) == extra['reference_sha256']
    for record in extra['source_identities']:
        assert sha(ROOT/record['path']) == record['sha256']
    rows.append(extra['result'])
    profiles = []
    table = []
    for row in rows:
        for j, mode in enumerate(row['modes']):
            path = ROOT/mode['profile_path']
            assert sha(path) == mode['profile_sha256']
            with path.open() as stream:
                data = list(csv.DictReader(stream))
            r = np.array([float(p['r_km']) for p in data])
            w = np.array([float(p['W']) for p in data])
            v = np.array([float(p['V_km_inverse']) for p in data])
            assert np.all(np.isfinite(r)) and np.all(np.isfinite(w)) and np.all(np.isfinite(v))
            wn,vn = count_nodes(w),count_nodes(v)
            center_ratio = -mode['ell']*r[0]*v[0]/w[0]
            assert (wn,vn) == (j,j) == (mode['W_nodes'],mode['V_nodes'])
            assert mode['root_converged'] and mode['bracket_endpoint_residuals'][0]*mode['bracket_endpoint_residuals'][1] < 0
            assert mode['bracket_hz'][0] < mode['frequency_hz'] < mode['bracket_hz'][1]
            accepted = abs(mode['matching_residual']) <= LIMIT and mode['relative_join_mismatch'] <= LIMIT
            assert accepted == mode['residual_accepted']
            assert mode['accepted'] == bool(accepted and mode['root_and_nodes_accepted'])
            assert abs(center_ratio-1.) < 1e-12
            assert np.all(np.diff(r)>0) and r[-1] <= row['background']['radius_km']
            profiles.append(dict(model=row['model'],configuration=row['label'],mode=mode['name'],
                                 W_nodes=wn,V_nodes=vn,profile_points=len(data),center_regularity_ratio=float(center_ratio),
                                 matching_residual=mode['matching_residual'],relative_join_mismatch=mode['relative_join_mismatch'],
                                 accepted=mode['accepted'],path=mode['profile_path'],sha256=mode['profile_sha256']))
        table.append(dict(model=row['model'],configuration=row['label'],source_grid=row['source_grid'],
                          central_density_fm3=row['background']['central_density_fm3'],mass_msun=row['background']['mass_msun'],radius_km=row['background']['radius_km'],
                          f_hz=row['modes'][0]['frequency_hz'],p1_hz=row['modes'][1]['frequency_hz'],
                          delta_f_hz=row['frequency_shifts_hz'][0],delta_p1_hz=row['frequency_shifts_hz'][1],
                          shift_reference=row['shift_reference'],accepted=row['accepted']))
    background_checks = []
    for model in result['models']:
        path = OUT/model['label']/'background_profile.csv'
        with path.open() as stream:
            data = list(csv.DictReader(stream))
        assert float(data[0]['r_km']) == 0 and data[0]['location'] == 'analytic_center_limit'
        assert float(data[-1]['h']) == 0 and float(data[-1]['pressure_MeV_fm3']) == 0
        redshifted = [(float(p['energy_MeV_fm3'])+float(p['pressure_MeV_fm3']))/float(p['n_fm3'])*np.exp(float(p['Phi'])) for p in data if float(p['n_fm3']) > 0]
        spread = float(np.ptp(redshifted)/np.mean(redshifted))
        assert spread < 1e-12
        background_checks.append(dict(model=model['label'],path=str(path.relative_to(ROOT)),sha256=sha(path),
                                      rows=len(data),mu_exp_Phi_relative_spread=spread))
    spans = {}
    for model in result['models']:
        accepted = [r for r in rows if r['model']==model['label'] and r['source_grid']=='fine' and r['accepted']]
        spans[model['label']] = [float(np.ptp([r['modes'][j]['frequency_hz'] for r in accepted])) for j in [0,1]]
    validation = dict(status='OUTPUT_IDENTITY_AND_PROFILE_READBACK_VERIFIED',execution_exit_code=execution['exit_code'],
                      targeted_background_exit_code=extra_execution['exit_code'],
                      residual_acceptance=LIMIT,identity_count=len(identities),configuration_count=len(rows),
                      mode_count=len(profiles),accepted_mode_count=sum(p['accepted'] for p in profiles),
                      rejected_attempt_count=sum(not p['accepted'] for p in profiles),accepted_fine_frequency_spans_hz=spans,
                      profiles=profiles,background_profiles=background_checks,
                      command=[sys.executable,str(Path(__file__).relative_to(ROOT))],verifier_sha256=sha(__file__))
    (OUT/'VALIDATION.json').write_text(json.dumps(validation,indent=2,allow_nan=False)+'\n')
    with (OUT/'MODE_TABLE.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)
    print(json.dumps({key:validation[key] for key in ['status','configuration_count','mode_count','accepted_mode_count','rejected_attempt_count','accepted_fine_frequency_spans_hz']}),flush=True)


if __name__ == '__main__':
    main()
