#!/usr/bin/env python3
"""Independent post-campaign audit of retained SPG states and native tables."""
from pathlib import Path
import sys,json,gzip,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.api import read_compose
BASE=ROOT/'data/research_finish_20261116/spg'
OUT=Path(__file__).parent
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
records=[];source_checks={}
for size in [129,257]:
 top=BASE/f'grid{size}';summary=load(top/'SUMMARY.json')
 for p,h in {**summary['source_sha256'],**summary['input_sha256']}.items():
  source_checks[p]=(sha(ROOT/p)==h)
 for recpath in sorted(top.glob('*/*/RECEIPT.json')):
  rec=load(recpath);root=recpath.parent;payload=json.load(gzip.open(root/'all_attempts.json.gz','rt'))
  states=payload['fine'];best=[r['best'] for r in states]
  native=read_compose(root/'compose')
  assert native.n==size==len(best)
  assert all(r['status']=='STATIONARY_CANDIDATE' and not r['violations'] for r in best)
  assert all(r['surface_parameters']==rec['surface_parameters'] for r in best)
  assert np.all(np.diff(native.nb)>0)
  mask_names=['Yn_present','Yp_present','Ye_present','Ymu_present','Yh_present','Abar_present','Zbar_present','C3_applicable','C4_applicable']
  assert all(getattr(native,k).all() for k in mask_names)
  def err(x,y):return float(np.max(abs(np.asarray(x)-np.asarray(y))))
  errors={'pressure':err(native.p,[r['pressure_MeV_fm3'] for r in best]),
    'epsilon':err(native.eps,[r['energy'] for r in best]),
    'mu_B':err(native.mub,[r['mu_B_MeV'] for r in best])}
  for k in ['Yn','Yp','Ye','Ymu','Yh','Aav','Zav']:
   errors[k]=err(getattr(native,k),[r['composition'][k] for r in best])
  errors['baryon_sum']=err(native.Yn+native.Yp+native.Yh,1)
  errors['charge_sum']=err(native.Yp+native.Yh*native.Zav/native.Aav,native.Ye+native.Ymu)
  assert max(errors.values())<3e-11,(root,errors)
  for p,h in rec['output_sha256'].items():assert sha(root/p)==h,(root,p)
  records.append(dict(grid=size,model=rec['model'],label=rec['label'],rows=size,
     max_projected_gradient=max(r['max_projected_gradient'] for r in best),errors=errors,
     hashes_checked=len(rec['output_sha256']),all_composition_masks_present=True,
     positive_gas_rows=sum(r.get('gas_protons')=='VARIATIONAL_POSITIVE' for r in best)))
comparison=load(BASE/'common_grid129/COMPARISON.json')
stencils=[]
for p in sorted((BASE/'common_grid129').glob('*/*/first_law_attempts.json.gz')):
 for c in json.load(gzip.open(p,'rt')):
  low,center,high=c['low'],c['center'],c['high']
  assert all(r['status']=='STATIONARY_CANDIDATE' for r in [low,center,high])
  # Recompute from the retained energies and actual input-density interval.
  nlow=low['n_B'];nhigh=high['n_B'];n=center['n_B']
  fd=(high['energy']-low['energy'])/(nhigh-nlow)
  stencils.append(dict(model=center['surface_parameters']['model_id'],label=center['surface_parameters']['label'],
    n_B=n,gas_protons=center.get('gas_protons'),relative_mu_error=abs(fd/center['mu_B_MeV']-1),
    absolute_pressure_error=abs(n*fd-center['energy']-center['pressure_MeV_fm3'])))
assert all(source_checks.values()),source_checks
assert len(records)==50 and len(stencils)==125
assert max(x['relative_mu_error'] for x in stencils)<1e-7
result=dict(scope='Native50tables, retained9650stationaryrows and125independently re-equilibrated stencils; no new optimizer solves in this audit',
 source_input_hash_matches=source_checks,table_count=len(records),row_count=sum(r['rows'] for r in records),
 max_projected_gradient=max(r['max_projected_gradient'] for r in records),tables=records,
 first_law_stencils=stencils,first_law_count=len(stencils),
 first_law_positive_gas_center_count=sum(x['gas_protons']=='VARIATIONAL_POSITIVE' for x in stencils),
 first_law_max_relative_mu_error=max(x['relative_mu_error'] for x in stencils),
 first_law_max_absolute_pressure_error=max(x['absolute_pressure_error'] for x in stencils),
 runner_sha256=sha(Path(__file__)))
(OUT/'SPG_TABLE_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['tables','first_law_stencils']},indent=2))
