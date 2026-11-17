#!/usr/bin/env python3
"""Narrow transition-report semantics: synthetic arithmetic and retained preview."""
from pathlib import Path
from copy import deepcopy
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[5];sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
import refine_gppva_geometry as rg
from compose_eos.tf_morphology import geometry_label,profile_morphology
OUT=Path(__file__).parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
paths=['scripts/refine_gppva_geometry.py','scripts/summarize_gppva_geometry.py','scripts/export_gppva_geometry.py','python_frontend/compose_eos/tf_morphology.py']
before={p:sha(ROOT/p) for p in paths}
def branches(delta):
 return {'rod':{'energy_per_baryon':950.,'pressure':.01,'mu_B':950.1},
         'slab':{'energy_per_baryon':950.+delta,'pressure':.0101,'mu_B':950.2}}
window=[{'left':.01,'right':.03,'geometries':['rod','slab']}]
def records(deltas):
 return [{'n_b':n,'points':N,'attempts':[branches(d) if d is not None else {}]}
         for n,pair in zip([.01,.02,.03],deltas) for N,d in zip([192,384],pair)]
old=rg.branch_minima
try:
 rg.branch_minima=lambda attempts:attempts[0]
 ulp=float(np.spacing(950.))
 near=rg.transition_evidence(window,records([[ulp,ulp],[-ulp,-ulp],[-ulp,-ulp]]),[192,384])[0]
 assert near['fine_mesh_sign_change_brackets'] and not near['fine_mesh_sign_change_brackets'][0]['both_endpoints_exceed_roundoff_scale']
 assert all(s['roundoff_near_degeneracy'] for s in near['samples'])
 gap=rg.transition_evidence(window,records([[.001,.001],[0.,None],[-.001,-.001]]),[192,384])[0]
 assert gap['fine_mesh_sign_change_brackets']==[]
 unresolved=rg.transition_evidence(window,records([[-.003,.001],[.003,-.001],[.003,-.001]]),[192,384])[0]
 assert not unresolved['fine_mesh_sign_change_brackets'][0]['both_signs_exceed_mesh_change']
finally:rg.branch_minima=old
def row(total,protons,dimension=2):
 total=np.asarray(total,float);protons=np.asarray(protons,float);N=len(total)
 return {'dimension':dimension,'profile':{'nn':(total-protons).tolist(),'np':protons.tolist(),
  'ne':np.full(N,float(np.mean(protons))).tolist(),'weights':np.diff(np.linspace(0,1,N+1)**dimension).tolist()}}
cases={'uniform':row([.1]*8,[.005]*8),
 'near_uniform':row([.1]*8,[.005-.000001,.005+.000001]*4),
 'cylindrical_shell':row([.05,.07,.1,.12,.12,.1,.07,.05],[.001,.005,.02,.03,.03,.02,.005,.001]),
 'spherical_shell':row([.05,.07,.1,.12,.12,.1,.07,.05],[.001,.005,.02,.03,.03,.02,.005,.001],3)}
expected={'uniform':'uniform','near_uniform':'near_uniform_profile','cylindrical_shell':'noncanonical_cylindrical_profile','spherical_shell':'noncanonical_spherical_profile'}
labels={k:geometry_label(r) for k,r in cases.items()};assert labels==expected
previews=[];inputs={}
for p in sorted((ROOT/'data/research_finish_20261116/geometry/exported_preview').glob('*/SUMMARY.json')):
 inputs[str(p.relative_to(ROOT))]=sha(p);report=json.loads(p.read_text())
 for w in report['transitions']:
  for sample in w['samples']:
   values=sample['mesh_values'];fine=values.get('384');coarse=values.get('192')
   if fine:
    delta=fine['energies_per_baryon'][1]-fine['energies_per_baryon'][0]
    assert delta==fine['delta_E_MeV_per_baryon']
    floor=64*np.finfo(float).eps*max(1.,*(abs(x) for x in fine['energies_per_baryon']))
    assert floor==sample['roundoff_scale_MeV_per_baryon']
    assert sample['roundoff_near_degeneracy']==(abs(delta)<=floor)
    if coarse:
     shift=abs(delta-coarse['delta_E_MeV_per_baryon'])
     assert shift==sample['adjacent_mesh_deltaE_change_MeV_per_baryon']
     assert sample['fine_sign_exceeds_mesh_change']==(abs(delta)>shift)
   assert 'neither an error bound' in sample['linearization_scope']
  assert 'not automatically a continuum transition' in w['interpretation']
  previews.append({'model':report['model_id'],'geometries':w['geometries'],'samples':len(w['samples']),
    'missing_fine_pair_count':sum('384' not in s['mesh_values'] for s in w['samples']),
    'brackets':w['fine_mesh_sign_change_brackets']})
result={'status':'PASS_SCOPED','scope':'Transition diagnostic arithmetic and morphology labels only; preview is not final ledger/native admission',
 'roundoff_counterexample':near,'missing_branch_counterexample':gap,
 'mesh_unresolved_counterexample':unresolved,'synthetic_shape_labels':labels,
 'preview_windows':previews,'preview_model_count':len({p['model'] for p in previews}),
 'preview_raw_sign_brackets':sum(len(p['brackets']) for p in previews),
 'input_sha256':inputs,'source_sha256':before,
 'source_unchanged':all(sha(ROOT/p)==h for p,h in before.items()),
 'runner_sha256':sha(Path(__file__))}
(OUT/'TRANSITION_SEMANTICS_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['roundoff_counterexample','missing_branch_counterexample','mesh_unresolved_counterexample','preview_windows','input_sha256']},indent=2))
assert result['source_unchanged']
