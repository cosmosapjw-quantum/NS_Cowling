#!/usr/bin/env python3
"""Independent replay of finite-family ledgers and corrected TF native exports."""
from pathlib import Path
from functools import lru_cache
import argparse, sys, json, gzip, hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.tf_composition import read_tf_composition
F={(3,'droplet'),(2,'droplet'),(1,'droplet'),(2,'hole'),(3,'hole'),(3,'uniform')}
PHASE={'uniform':3,'droplet':21,'rod':22,'slab':23,'tube':24,'bubble':25,
       'noncanonical_cylindrical_profile':26,'noncanonical_spherical_profile':27,
       'near_uniform_profile':28,'noncanonical_radial_profile':29}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
@lru_cache(None)
def load(name):
 p=Path(name)
 return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def number(r,k):
 try:
  v=float(r[k]);return v if np.isfinite(v) else None
 except (KeyError,TypeError,ValueError):return None
def feasible(r,n):
 names=['n_b','baryon_constraint_error','charge_constraint_error',
        'scalar_pde_residual','vector_pde_residual','electron_pde_residual']
 values=[number(r,k) for k in names]
 if any(v is None for v in values):return None
 density=max(1e-12,n*1e-8)
 limits=[density,density,density,1e-8,1e-8,2e-10]
 values[0]-=n
 return all(abs(v)<=tol for v,tol in zip(values,limits))
def decision(attempts,n):
 valid=[]
 for i,r in enumerate(attempts):
  if (r.get('numerical_status') in ('STATIONARY','UNIFORM_STATIONARY') and
      not r.get('muon_threshold_crossed',False) and feasible(r,n) is True and
      number(r,'epsilon') is not None and number(r,'energy_per_baryon') is not None):valid.append(i)
 best=min(valid,key=lambda i:attempts[i]['energy_per_baryon']) if valid else None
 missing=F-{(attempts[i]['dimension'],attempts[i]['seed']) for i in valid}
 lower=[] if best is None else [i for i,r in enumerate(attempts) if
   i not in valid and number(r,'energy_per_baryon') is not None and feasible(r,n) is not False
   and r['energy_per_baryon']<attempts[best]['energy_per_baryon']-1e-6]
 return best,valid,missing,lower
def occupancy(x):
 x=np.array(x,float)
 if np.ptp(x)==0:return False,False
 indices=np.flatnonzero(x>(x.min()+x.max())*.5)
 chunks=np.split(indices,np.flatnonzero(np.diff(indices)>1)+1)
 interior=len(chunks)==1 and len(indices)>0 and indices[0]!=0 and indices[-1]!=len(x)-1
 return interior,len(chunks)>1
def label(r):
 p=r['profile'];nn=np.array(p['nn']);np_=np.array(p['np']);ne=np.array(p['ne'])
 weights=np.array(p['weights']);weights/=weights.sum();total=nn+np_;scale=weights@total
 if np.ptp(total)/scale<1e-4:
  uniform=all(np.ptp(x)<=max(1e-18,scale*1e-12) or np.ptp(x)<=abs(weights@x)*1e-4 for x in (nn,np_,ne))
  return 'uniform' if uniform else 'near_uniform_profile'
 d=r['dimension'];b=occupancy(total);p=occupancy(np_)
 if d==1:return 'slab'
 if b[0] and p[0]:return {2:'noncanonical_cylindrical_profile',3:'noncanonical_spherical_profile'}[d]
 if (b[1] and p[1]) or b!=p:return 'noncanonical_radial_profile'
 return {2:('rod','tube'),3:('droplet','bubble')}[d][int(total[0]<=total[-1])]
def hashes(mapping,base):
 return {name:sha(base/name)==value for name,value in mapping.items()}
ap=argparse.ArgumentParser()
ap.add_argument('--base',default='data/research_finish_20261116/geometry/exported')
ap.add_argument('--output',default='GEOMETRY_EXPORT_REVIEW.json')
args=ap.parse_args();base=ROOT/args.base
nodes=[];tables=[];hash_records={};chosen={};attempt_counts=0;execution_failures=0
for ep in sorted(base.glob('*/N*/EVIDENCE.json')):
 parent=ep.parent;receipt=load(str(parent/'RECEIPT.json'));ledger=load(str(ep))
 checks=hashes(receipt['input_sha256'],ROOT);assert all(checks.values()),checks
 hash_records[str(ep.relative_to(ROOT))]={'evidence_sha256':sha(ep),'input_matches':checks,
  'export_source_matches':hashes(receipt['export_source_sha256'],ROOT)}
 assert all(hash_records[str(ep.relative_to(ROOT))]['export_source_matches'].values())
 rows=ledger['rows'];selected={};requested=[x['n_b'] for x in rows];missing_densities=[]
 for row in rows:
  raw=load(str(ROOT/row['base_attempt_file']))[row['base_record_index']]
  assert raw['n_b']==row['n_b'] and raw['points']==row['points']
  attempts=list(raw['attempts'])
  for extra in row['supplemental_attempts']:
   supplement=load(str(ROOT/extra['file']))
   attempts.append(supplement[extra['attempt_index']])
  best,valid,missing,lower=decision(attempts,row['n_b'])
  recorded=row['selection'];eligible=best is not None and not missing and not lower
  assert recorded['best_attempt_index']==best
  assert set(tuple(x) for x in recorded['unresolved_families'])==missing
  assert recorded['lower_unconverged_attempt_indices']==lower
  assert eligible==(recorded['status']=='LOWEST_EVALUATED_STATIONARY_CANDIDATE')
  attempt_counts+=len(attempts);execution_failures+=sum(x.get('execution')=='FAILED' for x in attempts)
  branch={}
  for j in valid:
   item=attempts[j];g=label(item)
   assert abs(item['epsilon']/row['n_b']-item['energy_per_baryon'])<1e-10
   assert abs(row['n_b']*item['mu_B']-item['epsilon']-item['pressure'])<1e-8
   branch[g]=min(branch.get(g,float('inf')),item['energy_per_baryon'])
  if eligible:selected[row['n_b']]=attempts[best]
  else:missing_densities.append(row['n_b'])
  nodes.append({'model':receipt['model_id'],'points':row['points'],'n_b':row['n_b'],
   'attempt_count':len(attempts),'stationary_count':len(valid),'selected':eligible,
   'best_attempt_index':best,'geometry':label(attempts[best]) if eligible else None,
   'branch_minima_MeV_per_baryon':branch,'unresolved_family_count':len(missing),
   'lower_unresolved_count':len(lower)})
 policy=load(str(parent/'INTERVAL_SUPPORT.json'))
 assert policy['unresolved_requested_densities_fm3']==missing_densities
 expected=[]
 for a,b in zip(requested,requested[1:]):
  if a in selected and b in selected and label(selected[a])!=label(selected[b]):
   expected.append({'density_bracket_fm3':[a,b],'phases':[label(selected[a]),label(selected[b])]})
 assert policy['phase_change_brackets']==expected
 exports=receipt['export'].get('segments',[receipt['export']])
 for e in exports:
  directory=parent/e.get('path','.')
  if 'output_sha256' in e:
   outcheck=hashes(e['output_sha256'],directory);assert all(outcheck.values())
   hash_records[str(directory.relative_to(ROOT))+'/native']=outcheck
  if not (directory/'compose/eos.thermo').exists():continue
  native=read_tf_composition(directory);t=native['thermodynamics']
  states=load(str(directory/'selected_states.json.gz'));assert t.n==len(states)
  maxerr=0.;maxfield=0.
  for i,state in enumerate(states):
   source=selected[state['n_b']]
   assert source['model_id']==receipt['model_id'] and source['points']==receipt['points']
   assert source['n_b']==state['n_b']
   for species in ['nn','np','ne']:
    array=np.asarray(source['profile'][species]);assert array.shape==(source['points'],)
    assert np.isfinite(array).all() and np.all(array>=0)
   weights=np.asarray(source['profile']['weights']);assert np.all(weights>0) and abs(weights.sum()-1)<1e-12
   assert state['profile']==source['profile']
   for field in ['epsilon','pressure','mu_B','mu_Q_lepton','radius_fm','energy_per_baryon']:
    assert state[field]==source[field]
   g=label(source);assert state['geometry']==g and int(t.Iphase[i])==PHASE[g]
   uniform=g=='uniform'
   for key in ['Yn_present','Yp_present','Yh_present','Abar_present','Zbar_present','C3_applicable','C4_applicable']:
    assert bool(getattr(t,key)[i])==uniform,(directory,i,key,g)
   assert t.Ye_present[i] and t.Ymu_present[i]
   p=source['profile'];w=np.array(p['weights']);n=source['n_b']
   fractions=np.array([w@np.array(p[k])/n for k in ['nn','np','ne']])
   maxerr=max(maxerr,abs(fractions[:2].sum()-1),abs(fractions[1]-fractions[2]),
              abs(t.Ye[i]-fractions[2]),float(np.max(abs(native['fractions'][i,:3]-fractions))))
   for key,target in [('nb',n),('eps',source['epsilon']),('p',source['pressure']),('mub',source['mu_B'])]:
    maxfield=max(maxfield,abs(float(getattr(t,key)[i])-target))
  assert maxerr<1e-8 and maxfield<1e-8
  tables.append({'directory':str(directory.relative_to(ROOT)),'rows':t.n,
    'max_composition_error':float(maxerr),'max_dimensional_field_error':float(maxfield),
    'phase_codes':t.Iphase.tolist(),'free_masks':t.Yn_present.tolist()})
 chosen[(receipt['model_id'],receipt['points'])]=len(selected)
result={'status':'PASS_SCOPED','scope':'Independent finite-family ledger reconstruction, native thermodynamics/integrated composition/masks, descriptive morphology, phase-gap manifest and byte provenance. No global minimum, common-tangent coexistence, continuum or full-crust certificate.',
 'base':str(base.relative_to(ROOT)),'model_count':len(set(x[0] for x in chosen)),
 'mesh_density_nodes':len(nodes),'selected_nodes':sum(x['selected'] for x in nodes),
 'attempt_count':attempt_counts,'preserved_execution_failures':execution_failures,
 'native_table_count':len(tables),'native_rows':sum(x['rows'] for x in tables),
 'records':nodes,'tables':tables,'hash_checks':hash_records,'runner_sha256':sha(Path(__file__))}
Path(__file__).with_name(args.output).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('records','tables','hash_checks')},indent=2))
