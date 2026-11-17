#!/usr/bin/env python3
"""Independent finite-family, radial-profile and native geometry review."""
from pathlib import Path
import sys,json,gzip,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.api import read_compose
BASE=ROOT/'data/research_finish_20261116/geometry/final'
F={(3,'droplet'),(2,'droplet'),(1,'droplet'),(2,'hole'),(3,'hole'),(3,'uniform')}
def load(p):return json.loads(p.read_text())
def zipped(p):return json.loads(gzip.decompress(p.read_bytes()))
def number(row,k):
 try:v=float(row[k]);return v if np.isfinite(v) else None
 except (KeyError,TypeError,ValueError):return None
def feasible(r,n):
 keys=['n_b','baryon_constraint_error','charge_constraint_error','scalar_pde_residual','vector_pde_residual','electron_pde_residual']
 v={k:number(r,k) for k in keys}
 if any(x is None for x in v.values()):return None
 tol=max(1e-12,n*1e-8)
 return abs(v['n_b']-n)<=tol and v['baryon_constraint_error']<=tol and v['charge_constraint_error']<=tol and v['scalar_pde_residual']<=1e-8 and v['vector_pde_residual']<=1e-8 and v['electron_pde_residual']<=2e-10
def geometry(r):
 if r['numerical_status']=='UNIFORM_STATIONARY' or r['baryon_density_contrast']<1e-4:return 'uniform'
 p=r['profile'];n=np.asarray(p['nn'])+np.asarray(p['np'])
 return 'slab' if r['dimension']==1 else {2:('rod','tube'),3:('droplet','bubble')}[r['dimension']][int(n[0]<n[-1])]
def shape(r):
 p=r['profile'];protons=np.array(p['np']);total=protons+np.array(p['nn']);scale=max(float(np.ptp(protons)),1e-15)
 total_scale=max(float(np.ptp(total)),1e-15)
 pp=float((max(protons)-max(protons[0],protons[-1]))/scale)
 bb=float((max(total)-max(total[0],total[-1]))/total_scale)
 delta=np.diff(protons);signs=np.sign(delta[abs(delta)>scale*1e-5]);turns=int(np.count_nonzero(np.diff(signs))) if len(signs)>1 else 0
 peaks=int(np.count_nonzero(np.diff(signs)<0)) if len(signs)>1 else 0
 return dict(geometry=geometry(r),proton_turns=turns,interior_proton_maxima=peaks,
   proton_total_edge_orientation_disagreement=bool((protons[0]-protons[-1])*(total[0]-total[-1])<0),
   significant_orientation_disagreement=bool(geometry(r) not in ['uniform','slab'] and abs(protons[0]-protons[-1])>scale*1e-3 and abs(total[0]-total[-1])>total_scale*1e-3 and (protons[0]-protons[-1])*(total[0]-total[-1])<0),
   shell_like_radial_shape=bool(geometry(r) not in ['uniform','slab'] and pp>.5 and bb>.5),
   baryon_normalized_interior_peak_above_both_edges=bb,
   normalized_interior_peak_above_both_edges=float(pp),
   relative_neutron_proton_sum_error=float(abs(np.array(p['weights'])@total/r['n_b']-1)))
def native_uniform(r):
 p=r['profile'];w=np.asarray(p['weights'],float);w=w/w.sum()
 total=np.asarray(p['nn'])+np.asarray(p['np']);scale=float(w@total)
 if np.ptp(total)/scale>=1e-4:return False
 return all(np.ptp(np.asarray(p[k]))<=max(1e-18,1e-12*scale) or np.ptp(np.asarray(p[k]))<=1e-4*abs(w@np.asarray(p[k])) for k in ['nn','np','ne'])
records=[];native=[];input_sha={}
for path in sorted(BASE.glob('*/attempts_N*.json.gz')):
 input_sha[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
 for record in zipped(path):
  n=record['n_b'];attempts=record['attempts']
  valid=[r for r in attempts if r.get('numerical_status') in ['STATIONARY','UNIFORM_STATIONARY'] and
   not r.get('muon_threshold_crossed',False) and feasible(r,n) is True and number(r,'epsilon') is not None and number(r,'energy_per_baryon') is not None]
  families={(r['dimension'],r['seed']) for r in valid};missing=F-families
  best=min(valid,key=lambda r:r['energy_per_baryon']) if valid else None
  lower=[r for r in attempts if r not in valid and number(r,'energy_per_baryon') is not None and
   feasible(r,n) is not False and best is not None and r['energy_per_baryon']<best['energy_per_baryon']-1e-6]
  eligible=best is not None and not missing and not lower
  assert eligible==(record['selected_geometry'] is not None),(path,n,missing,len(lower))
  if eligible:assert geometry(best)==record['selected_geometry']
  branch={}
  for r in valid:
   g=geometry(r)
   if g not in branch or r['energy_per_baryon']<branch[g]['energy_per_baryon']:branch[g]=r
  records.append(dict(model=path.parent.name,n_b=n,points=record['points'],status=record['status'],
    selected_geometry=record['selected_geometry'],missing=[list(x) for x in sorted(missing)],
    lower_unresolved_count=len(lower),attempt_count=len(attempts),stationary_count=len(valid),
    selected_shape=shape(best) if eligible else None,
    branch_minima={g:dict(energy_per_baryon=r['energy_per_baryon'],shape=shape(r)) for g,r in branch.items()}))
for p in sorted(BASE.glob('*/N*/**/compose/eos.thermo')):
 directory=p.parent;table=read_compose(directory);profile=directory.parent/'selected_states.json.gz'
 states=zipped(profile);assert table.n==len(states)
 err=[];mask_errors=0
 for i,s in enumerate(states):
  g=geometry(s);uniform=native_uniform(s);w=np.array(s['profile']['weights'])
  ye=w@np.array(s['profile']['ne'])/s['n_b'];yp=w@np.array(s['profile']['np'])/s['n_b']
  err.append(max(abs(ye-table.Ye[i]),abs(yp-ye)))
  assert all(bool(getattr(table,key)[i])==uniform for key in ['Yn_present','Yp_present','Yh_present','Abar_present','Zbar_present','C3_applicable','C4_applicable'])
  assert table.Ye_present[i] and table.Ymu_present[i]
 assert max(err)<1e-8
 native.append(dict(path=str(directory.relative_to(ROOT)),rows=table.n,max_charge_readback_error=max(err)))
campaign=load(BASE/'CAMPAIGN.json') if (BASE/'CAMPAIGN.json').exists() else None
result=dict(scope='Finite radial profiles/branch selection and native sparse composition; no arbitrary morphology or continuum/global-minimum certificate',
 campaign_complete=campaign is not None,source_unchanged=campaign.get('source_unchanged_during_execution') if campaign else None,
 input_sha256=input_sha,node_count=len(records),selected_nodes=sum(r['selected_geometry'] is not None for r in records),records=records,native_tables=native,
 runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(__file__).with_name('GEOMETRY_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['records','native_tables','input_sha256']},indent=2))
print('native_table_count',len(native))
print('orientation_disagreements',sum(bool(r['selected_shape'] and r['selected_shape']['proton_total_edge_orientation_disagreement'] and r['selected_geometry']!='uniform') for r in records))
print('multiple_proton_maxima',[(r['model'],r['n_b'],r['points'],r['selected_shape']) for r in records if r['selected_shape'] and r['selected_geometry']!='uniform' and r['selected_shape']['interior_proton_maxima']>1])

print('selected_shell_like',[(r['model'],r['n_b'],r['points']) for r in records if r['selected_shape'] and r['selected_shape']['shell_like_radial_shape']])
print('branch_shell_like',[(r['model'],r['n_b'],r['points'],g) for r in records for g,b in r['branch_minima'].items() if b['shape']['shell_like_radial_shape']])
