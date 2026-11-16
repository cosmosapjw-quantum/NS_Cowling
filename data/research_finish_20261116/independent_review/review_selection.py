#!/usr/bin/env python3
"""Replay selection counterexamples without invoking an optimizer."""
from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.tf_selection import select_geometry_candidates
F=[(3,'droplet'),(2,'droplet'),(1,'droplet')]
def row(d,seed,energy,status='STATIONARY',**changes):
 r=dict(dimension=d,seed=seed,energy_per_baryon=energy,epsilon=energy*.05,
        n_b=.05,baryon_constraint_error=0.,charge_constraint_error=0.,
        scalar_pde_residual=0.,vector_pde_residual=0.,electron_pde_residual=0.,
        numerical_status=status,muon_threshold_crossed=False)
 r.update(changes);return r
base=[row(d,s,900+i) for i,(d,s) in enumerate(F)]
cases=[
 ('complete',base,'LOWEST_EVALUATED_STATIONARY_CANDIDATE',True),
 ('missing_rod',base[::2],'HOLD_UNRESOLVED_GEOMETRY',False),
 ('exception_and_restart',base+[dict(dimension=2,seed='droplet',execution='FAILED')],'LOWEST_EVALUATED_STATIONARY_CANDIDATE',True),
 ('lower_feasible_unstationary',base+[row(2,'droplet',899,'HOLD')],'LOWER_UNCONVERGED_ATTEMPT',False),
 ('lower_infeasible',base+[row(2,'droplet',899,'HOLD',baryon_constraint_error=.001)],'LOWEST_EVALUATED_STATIONARY_CANDIDATE',True),
 ('lower_unknown',base+[dict(dimension=2,seed='droplet',energy_per_baryon=899)],'LOWER_UNCONVERGED_ATTEMPT',False),
 ('false_stationary_field_failure',[base[0],row(2,'droplet',899,electron_pde_residual=1e-8),base[2]],'HOLD_UNRESOLVED_GEOMETRY',False),
 ('wrong_density',[base[0],row(2,'droplet',899,n_b=.06),base[2]],'HOLD_UNRESOLVED_GEOMETRY',False),
 ('higher_unstationary',base+[row(2,'droplet',904,'HOLD')],'LOWEST_EVALUATED_STATIONARY_CANDIDATE',True)]
out=[]
for label,rows,status,selected in cases:
 d=select_geometry_candidates(rows,F,n_b=.05)
 assert d['status']==status,(label,d)
 assert (d['selected'] is not None)==selected
 out.append(dict(case=label,status=d['status'],selected=d['selected'] is not None,
                 missing=d['unresolved_families'],lower=d['lower_unconverged_attempt_indices']))
module=ROOT/'python_frontend/compose_eos/tf_selection.py'
r=dict(scope='Nine synthetic selection counterexamples; solver stationarity status is an input contract, not independently inferred here',
       source_sha256=hashlib.sha256(module.read_bytes()).hexdigest(),cases=out)
Path(__file__).with_name('SELECTION_ORACLE.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
