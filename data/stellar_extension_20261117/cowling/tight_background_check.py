#!/usr/bin/env python3
"""Target the observed Stage E TOV integration-noise risk at the selected high star."""
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import time
from compose_eos.stellar_background import ContractEnthalpyEOS, solve_background
from compose_eos.stellar_cowling import CowlingSettings, matching_residual, solve_mode

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
LIMIT=2e-6


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    start=time.monotonic()
    summary=json.loads((OUT/'SUMMARY.json').read_text())
    provenance=json.loads((OUT/'PROVENANCE.json').read_text())
    reference=next(row for row in summary['configurations'] if row['model']=='endpoint' and row['label']=='ode_tight')
    source=provenance['primary_contract']
    assert sha(ROOT/source['path'])==source['sha256']
    eos=ContractEnthalpyEOS.from_path(ROOT/source['path'])
    bg=solve_background(eos,reference['background']['central_density_fm3'],r0_km=.001,rtol=1e-12,atol=1e-14,max_step_h=.0001)
    settings=CowlingSettings(**reference['settings'])
    modes=[]
    for j,old in enumerate(reference['modes']):
        bracket=old['bracket_hz']
        signs=[matching_residual(bg,f,settings) for f in bracket]
        mode=solve_mode(bg,bracket,settings)
        profile=mode.pop('profile')
        path=OUT/'endpoint'/('background_tight_'+old['name']+'_profile.csv')
        with path.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(profile[0]));writer.writeheader();writer.writerows(profile)
        residual_ok=abs(mode['matching_residual'])<=LIMIT and mode['relative_join_mismatch']<=LIMIT
        node_ok=mode['W_nodes']==j and mode['V_nodes']==j
        mode.update(name=old['name'],bracket_endpoint_residuals=signs,residual_accepted=bool(residual_ok),
                    root_and_nodes_accepted=bool(node_ok and mode['root_converged'] and signs[0]*signs[1]<0),
                    profile_path=str(path.relative_to(ROOT)),profile_sha256=sha(path))
        mode['accepted']=bool(mode['residual_accepted'] and mode['root_and_nodes_accepted'])
        modes.append(mode)
        print(old['name'],mode['frequency_hz'],mode['matching_residual'],mode['accepted'],flush=True)
    row=dict(model='endpoint',label='background_tight',source_grid='fine',settings=asdict(settings),
             background=bg.summary(),modes=modes,accepted=all(mode['accepted'] for mode in modes),
             frequency_shifts_hz=[m['frequency_hz']-old['frequency_hz'] for m,old in zip(modes,reference['modes'])],
             shift_reference='fine_ode_tight',wall_seconds=time.monotonic()-start)
    (OUT/'endpoint'/'background_tight.json').write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
    result=dict(status='TARGETED_BACKGROUND_COMPARISON_COMPLETE',reason='TOV refinement exposed local integration noise; test its effect on the fixed high-density mode observables',
                result=row,reference_path=str((OUT/'endpoint'/'ode_tight.json').relative_to(ROOT)),
                reference_sha256=sha(OUT/'endpoint'/'ode_tight.json'),primary_contract=source,
                residual_acceptance=LIMIT,script_sha256=sha(__file__),command=[sys.executable,str(Path(__file__).relative_to(ROOT))],
                source_identities=[dict(path=path,sha256=sha(ROOT/path)) for path in ['python_frontend/compose_eos/stellar_background.py','python_frontend/compose_eos/stellar_cowling.py']])
    (OUT/'TIGHT_BACKGROUND_COMPARISON.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return 0 if row['accepted'] else 2


if __name__=='__main__':
    raise SystemExit(main())
