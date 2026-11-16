"""Actual DD2 .0625 rod minimizations after screened-electron repair."""
from pathlib import Path
import json,gzip,hashlib,sys,time,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from compose_eos.spatial_tf_equilibrium import solve_screened_cell
from compose_eos.spatial_tf import RadialGrid
from tf_polish import profile_initial,polish_profile
OUT=Path(__file__).resolve().parent
seedfile=ROOT/'data/research_completion_20261116/gppva_tf/final/GPPVA_DD2/selected_states.json.gz'
old=next(r for r in json.loads(gzip.decompress(seedfile.read_bytes())) if r['n_b']==.065 and r['points']==192)
paths=['python_frontend/compose_eos/spatial_tf_equilibrium.py','python_frontend/compose_eos/spatial_tf.py','python_frontend/compose_eos/gppva_rmf.py']
hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
records=[];start=time.monotonic()
for label in ('warm','cold'):
    kwargs=dict(points=192,dimension=old['dimension'],seed=old['seed'],radius_bounds=(5.,180.),maxiter=900,ftol=1e-13)
    if label=='warm':
        initial=profile_initial(old,192)
        for k in ('nn','np','ne'):initial[k]*=.0625/old['n_b']
        kwargs['initial']=initial
    else:kwargs['radius']=30.
    try:
        row=solve_screened_cell('GPPVA_DD2',.0625,**kwargs)
        first=row
        if row['numerical_status']=='REQUIRES_REFINEMENT':row=polish_profile(row,192)
        record={'label':label,'first_attempt':first,'selected':row}
        print(label,row['numerical_status'],row['energy_per_baryon'],row['max_chemical_KKT_residual_MeV'],flush=True)
    except Exception as exc:
        record={'label':label,'execution':'FAILED','exception':str(exc),'traceback':traceback.format_exc()}
        print(label,'FAILED',str(exc),flush=True)
    records.append(record)
    (OUT/'REPAIRED_CELLS.json.gz').write_bytes(gzip.compress(json.dumps(records,allow_nan=False).encode(),mtime=0))
receipt={'source_sha256':hashes,'source_unchanged':all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in hashes.items()),'seed_sha256':hashlib.sha256(seedfile.read_bytes()).hexdigest(),'runtime_seconds':time.monotonic()-start,'cases':[{'label':r['label'],**{k:r['selected'][k] for k in ('numerical_status','energy_per_baryon','pressure','Y_p','max_chemical_KKT_residual_MeV','electron_pde_residual','charge_constraint_error','baryon_constraint_error','radius_stationarity_MeV')}} if 'selected'in r else r for r in records]}
(OUT/'GREEN_RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False))
