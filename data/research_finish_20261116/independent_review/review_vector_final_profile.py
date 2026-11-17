#!/usr/bin/env python3
"""Independent old/prototype operator reevaluation of new stationary N384 cell."""
from pathlib import Path
import sys,importlib.util,json,gzip,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'python_frontend'))
P=ROOT/'data/research_finish_20261116/vector_band_prototype'
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
dense=module('compose_eos._independent_vector_dense',P/'baseline_spatial_tf.py')
band=module('compose_eos._independent_vector_band',P/'band_spatial_tf.py')
path=P/'FSU2R_N384_CELL.json.gz';r=json.loads(gzip.decompress(path.read_bytes()))
# The owning runner may wrap the result; preserve explicit support for that schema.
if 'profile' not in r:r=r['cell']
n=r['n_b'];p=r['profile'];nn,np_,ne=[np.array(p[k]) for k in ['nn','np','ne']]
states={};physics={}
for label,mod in [('dense',dense),('band',band)]:
 grid=mod.RadialGrid(r['points'],r['dimension']);f=mod.CellFunctional(r['model_id'],grid,n)
 s=f.evaluate(nn,np_,ne,r['radius_fm']);states[label]=s
 active=[nn>max(1e-10,n*1e-7),np_>max(1e-10,n*1e-7)]
 weights=np.r_[grid.weights[active[0]],grid.weights[active[1]]]
 chemical=np.r_[s['mu_n'][active[0]],(s['mu_p']+r['mu_Q_lepton'])[active[1]]]
 mu=float(weights@chemical/sum(weights));kkt=[]
 for mask,q in zip(active,[s['mu_n'],s['mu_p']+r['mu_Q_lepton']]):
  kkt.extend(abs(q[mask]-mu));kkt.extend(np.maximum(mu-q[~mask],0.))
 kkt.extend(abs(s['mu_e']-r['mu_Q_lepton']))
 physics[label]=dict(epsilon=float(s['epsilon']),mu_B=mu,pressure=float(n*mu-s['epsilon']),
  KKT=float(max(kkt)),radius_stationarity=abs(float(s['d_epsilon_d_log_radius']/n)),
  scalar_PDE=float(s['scalar_pde_residual']),vector_PDE=float(s['vector_pde_residual']),
  baryon_constraint=abs(float(s['mean_baryon_density']-n)),charge_constraint=abs(float(s['mean_charge_density'])))
diffs={k:float(np.max(abs(np.asarray(states['dense'][k])-np.asarray(states['band'][k])))) for k in states['dense']}
assert r['numerical_status']=='STATIONARY'
assert max(diffs[k] for k in ['phi','omega','rho','epsilon'])<1e-10
assert abs(physics['dense']['pressure']-physics['band']['pressure'])<1e-10
assert all(s['KKT']<.02 and s['radius_stationarity']<2e-4 and s['scalar_PDE']<1e-8 and s['vector_PDE']<1e-8 for s in physics.values())
result=dict(scope='Actual stationary FSU2R N384 finalprofile independently reevaluated with cold dense and band field solvers; this is not a whole-optimization speed comparison',
 model=r['model_id'],n_B=n,points=r['points'],numerical_status=r['numerical_status'],
 fields=physics,max_state_differences=diffs,
 source_sha256={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in [P/'baseline_spatial_tf.py',P/'band_spatial_tf.py',path]},
 runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(__file__).with_name('VECTOR_FINAL_PROFILE_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
