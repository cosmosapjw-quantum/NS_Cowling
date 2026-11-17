"""Actual N384 FSU2R minimization with isolated banded-field source copies."""
from pathlib import Path
import sys,hashlib,gzip,json,importlib.util,time,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from tf_polish import profile_initial

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec)
 sys.modules[name]=mod;spec.loader.exec_module(mod);return mod
base=load('compose_eos._vector_dense_prototype',OUT/'baseline_spatial_tf.py')
band=load('compose_eos._vector_band_prototype',OUT/'band_spatial_tf.py')
eq=load('compose_eos._vector_band_eq',OUT/'band_spatial_tf_equilibrium.py')
production={str(p.relative_to(ROOT)):digest(p) for p in [ROOT/'python_frontend/compose_eos/spatial_tf.py',ROOT/'python_frontend/compose_eos/spatial_tf_equilibrium.py']}
source=ROOT/'data/research_completion_20261116/gppva_tf/final/GPPVA_FSU2R/selected_states.json.gz'
old=min(json.loads(gzip.decompress(source.read_bytes())),key=lambda r:abs(r['n_b']-.05))
initial=profile_initial(old,384)
params=dict(points=384,dimension=old['dimension'],seed=old['seed'],radius_bounds=(5.,180.),maxiter=900,ftol=1e-13)
start=time.perf_counter()
try:
 row=eq.solve_screened_cell('GPPVA_FSU2R',.05,initial=initial,**params)
 elapsed=time.perf_counter()-start
 (OUT/'FSU2R_N384_CELL.json.gz').write_bytes(gzip.compress(json.dumps(row).encode()))
 grid=base.RadialGrid(384,row['dimension']);nn,np_,ne=[np.asarray(row['profile'][k]) for k in ('nn','np','ne')]
 states={};timings={}
 for label,mod in [('dense',base),('band',band)]:
  f=mod.CellFunctional(row['model_id'],grid,row['n_b']);t=time.perf_counter()
  states[label]=f.evaluate(nn,np_,ne,row['radius_fm']);timings[label]=time.perf_counter()-t
 diffs={k:float(np.max(abs(np.asarray(states['dense'][k])-np.asarray(states['band'][k])))) for k in states['dense']}
 pressure={}
 for label,s in states.items():
  active=[nn>max(1e-10,row['n_b']*1e-7),np_>max(1e-10,row['n_b']*1e-7)]
  weights=np.r_[grid.weights[active[0]],grid.weights[active[1]]]
  chemical=np.r_[s['mu_n'][active[0]],(s['mu_p']+row['mu_Q_lepton'])[active[1]]]
  mub=float(weights@chemical/sum(weights));pressure[label]=row['n_b']*mub-s['epsilon']
 receipt={'execution':'COMPLETED','wall_seconds':elapsed,'params':params,'source_sha256':production,
          'prototype_source_sha256':{f:digest(OUT/f) for f in ['band_spatial_tf.py','band_spatial_tf_equilibrium.py']},
          'seed_sha256':digest(source),'summary':{k:v for k,v in row.items() if k!='profile'},
          'final_profile_field_reevaluation_max_abs_differences':diffs,'final_profile_reevaluation_pressure':pressure,
          'final_profile_reevaluation_seconds':timings,'production_unchanged_during_execution':all(digest(ROOT/p)==v for p,v in production.items())}
 assert row['numerical_status']=='STATIONARY',row['numerical_status']
 assert max(diffs[k] for k in ('phi','omega','rho','epsilon'))<1e-10
 assert abs(pressure['dense']-pressure['band'])<1e-10
 print('STATIONARY',elapsed,'seconds',row['energy_per_baryon'],row['max_chemical_KKT_residual_MeV'],flush=True)
except Exception as exc:
 receipt=locals().get('receipt',{})
 receipt.update(execution='FAILED',exception=type(exc).__name__,message=str(exc),traceback=traceback.format_exc(),wall_seconds=time.perf_counter()-start)
 print(traceback.format_exc(),flush=True)
(OUT/'FSU2R_N384_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
