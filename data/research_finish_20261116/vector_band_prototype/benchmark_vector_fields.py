"""Compare copied dense/banded full field and reduced-functional evaluations."""
from pathlib import Path
import sys,json,gzip,hashlib,importlib.util,time,statistics,platform
import numpy as np
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'python_frontend'))
import compose_eos.spatial_tf_equilibrium as eq
import compose_eos.spatial_tf as live

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path)
 mod=importlib.util.module_from_spec(spec);sys.modules[name]=mod;spec.loader.exec_module(mod);return mod
base=load('compose_eos._vector_dense_prototype',OUT/'baseline_spatial_tf.py')
band=load('compose_eos._vector_band_prototype',OUT/'band_spatial_tf.py')
class DenseReduced(eq.ReducedCellFunctional):fields=base.CellFunctional.fields
class BandReduced(eq.ReducedCellFunctional):fields=band.CellFunctional.fields
production=ROOT/'python_frontend/compose_eos/spatial_tf.py';before=digest(production)
records=[];profile_hashes={}
for model in ('GPPVA_FSU2R','GPPVA_FSU2H','GPPVA_TM1e'):
 src=ROOT/'data/research_completion_20261116/gppva_tf/final'/model/'selected_states.json.gz'
 profile_hashes[str(src.relative_to(ROOT))]=digest(src)
 old=min(json.loads(gzip.decompress(src.read_bytes())),key=lambda r:abs(r['n_b']-.05))
 for points in (192,384):
  grid=base.RadialGrid(points,old['dimension']);prof=old['profile'];oldx=np.asarray(prof['x'])
  nn,np_,ne=[np.interp(grid.x,oldx,prof[k]) for k in ('nn','np','ne')]
  radius=old['radius_fm'];nbar=grid.mean(nn+np_)
  retained=tuple(np.interp(grid.x,oldx,prof[k]) for k in ('phi','omega','rho'))
  for kind in ('COLD_FIELDS','PERTURBED_RETAINED_FIELDS'):
   guesses=None if kind=='COLD_FIELDS' else tuple(v*factor for v,factor in zip(retained,(.9,1.3,.7)))
   values={};times={};objects={}
   for label,mod in [('dense',base),('band',band)]:
    timings=[]
    for repeat in range(5):
     f=mod.CellFunctional(model,grid,nbar)
     f._fields=None if guesses is None else tuple(v.copy() for v in guesses)
     start=time.perf_counter();state=f.evaluate(nn,np_,ne,radius);timings.append(time.perf_counter()-start)
    values[label]=state;times[label]=timings;objects[label]=f
   a,b=values['dense'],values['band']
   diffs={k:float(np.max(abs(np.asarray(a[k])-np.asarray(b[k])))) for k in a}
   pressures={};physical={}
   for label,state in values.items():
    active=[nn>max(1e-10,nbar*1e-7),np_>max(1e-10,nbar*1e-7)]
    coeff=np.r_[grid.weights[active[0]],grid.weights[active[1]]]
    muq=grid.mean(state['mu_e'])
    chemical=np.r_[state['mu_n'][active[0]],(state['mu_p']+muq)[active[1]]]
    mub=float(coeff@chemical/sum(coeff));pressures[label]=nbar*mub-state['epsilon']
    physical[label]={k:float(state[k]) for k in ('epsilon','scalar_pde_residual','vector_pde_residual','d_epsilon_d_log_radius')}
   reduced={};rtimes={}
   x=np.r_[nn*grid.sqrtw/nbar,np_*grid.sqrtw/nbar,np.log(radius)]
   for label,cls in [('dense',DenseReduced),('band',BandReduced)]:
    timings=[]
    for repeat in range(5):
     f=cls(model,grid,nbar);f._fields=None if guesses is None else tuple(v.copy() for v in guesses)
     start=time.perf_counter();value,gradient=f.objective(x);timings.append(time.perf_counter()-start)
    reduced[label]=(value,gradient);rtimes[label]=timings
   row={'model_id':model,'source_n_b':old['n_b'],'points':points,'dimension':old['dimension'],'radius_fm':radius,
        'profile_kind':'RETAINED' if points==old['points'] else 'LINEAR_RESAMPLE_OF_RETAINED_N192',
        'initial_fields':kind,'max_absolute_state_differences':diffs,'fixed_profile_pressure':pressures,
        'pressure_difference':abs(pressures['dense']-pressures['band']),'physical_residuals':physical,
        'field_evaluation_seconds':times,'field_speedup_median':statistics.median(times['dense'])/statistics.median(times['band']),
        'reduced_objective_seconds':rtimes,'reduced_speedup_median':statistics.median(rtimes['dense'])/statistics.median(rtimes['band']),
        'reduced_value_difference':abs(float(reduced['dense'][0]-reduced['band'][0])),
        'reduced_gradient_max_abs_difference':float(np.max(abs(reduced['dense'][1]-reduced['band'][1])))}
   assert diffs['epsilon']<1e-10 and max(diffs[k] for k in ('phi','omega','rho'))<1e-10
   assert row['pressure_difference']<1e-10 and row['reduced_gradient_max_abs_difference']<1e-9
   assert all(v['scalar_pde_residual']<2e-8 and v['vector_pde_residual']<2e-8 for v in physical.values())
   records.append(row);print(model,points,kind,'field speedup',row['field_speedup_median'],'reduced',row['reduced_speedup_median'],'maxfield',max(diffs[k] for k in ('phi','omega','rho')),flush=True)
# Direct kernel adversaries: nonsymmetric radial-like blocks, indefinite local
# blocks, zero leading pivot, zero cross, and an exactly singular matrix.
rng=np.random.default_rng(97531);kernel=[]
for kind in ('nonsymmetric','indefinite','pivot','zero_cross','singular'):
 n=16
 a=np.diag(np.linspace(2,3,n))+np.diag(np.full(n-1,-.3),1)+np.diag(np.full(n-1,-.7),-1)
 b=a*.8;cross=np.full(n,.2)
 if kind=='indefinite':a-=5*np.eye(n)
 if kind=='pivot':a[0,0]=0.;cross[0]=1.
 if kind=='zero_cross':cross*=0
 if kind=='singular':a*=0;b*=0;cross*=0
 rhs=rng.normal(size=2*n);dense=np.block([[a,np.diag(cross)],[np.diag(cross),b]])
 results={}
 for label in ('dense','band'):
  try:
   result=np.linalg.solve(dense,rhs) if label=='dense' else np.concatenate(band._coupled_band_solve(a,b,cross,rhs[:n],rhs[n:]))
   results[label]={'status':'SOLVED','residual':float(np.max(abs(dense@result-rhs))),'solution':result.tolist()}
  except Exception as exc:results[label]={'status':'RAISED','type':type(exc).__name__,'message':str(exc)}
 assert results['dense']['status']==results['band']['status']
 if results['dense']['status']=='SOLVED':
  assert max(results['dense']['residual'],results['band']['residual'])<1e-12
  error=float(np.max(abs(np.asarray(results['dense'].pop('solution'))-np.asarray(results['band'].pop('solution')))))
 else:error=None;assert results['dense']['type']==results['band']['type']=='LinAlgError'
 kernel.append({'kind':kind,'results':results,'max_abs_difference':error})
assert before==digest(production)
receipt={'production_source_unchanged':True,'production_sha256':before,
         'prototype_sha256':digest(OUT/'band_spatial_tf.py'),'prototype_patch_sha256':digest(OUT/'PROPOSED.patch'),
         'profile_source_sha256':profile_hashes,'cases':records,'kernel_cases':kernel,
         'environment':{'python':platform.python_version(),'numpy':np.__version__,'threads':1,'timing_repeats':5},
         'scope':'Fixed physical profile field solves and reduced objective evaluations; N384 profiles are resampled N192 profiles. No full minimization timing or new stationary EOS result.'}
(OUT/'BENCHMARK.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS',len(records),'full-field/objective pairs;',len(kernel),'kernel adversaries; production source unchanged',flush=True)
