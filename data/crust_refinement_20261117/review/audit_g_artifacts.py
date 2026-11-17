"""Read-only independent evidence/provenance audit of G native production."""
from pathlib import Path
import ast,json,hashlib,sys,time,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'data/crust_refinement_20261117'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def fns(p):return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(Path(p).read_text()).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
before=BASE/'tf_regeneration/tf_candidate.interim_projection.py.txt'
after=ROOT/'python_frontend/compose_eos/tf_candidate.py'
old,new=fns(before),fns(after)
changed=[k for k in new if old.get(k)!=new[k]]
assert changed==['export_compose'],changed
assert sha(after)=='61ad299bfeab5ee5588758b6367f7b69e3911317257dcecea23665f08a3c4c18'
# These are byte checks, not another native solve or rematch.
outputs=[]
for folder in [BASE/'tf_regeneration/tables_v2',BASE/'support_repair/matching']:
 for p in sorted(folder.rglob('RESULT.json')):
  result=load(p)
  if result.get('status')!='GENERATED_VERIFIED_FINITE_DOMAIN':continue
  cp=Path(result.get('compose_path',result.get('output_path','')))
  if not cp.is_absolute():cp=ROOT/cp
  hashes=result.get('table_sha256',result.get('output_sha256'))
  assert hashes and all(sha(cp/name)==v for name,v in hashes.items())
  assert result['admission']['status']=='ADMITTED'
  verification=result.get('consistency',result.get('serialized_native_verification'))
  assert verification['overall']=='PASS'
  outputs.append(dict(result=str(p.relative_to(ROOT)),result_sha256=sha(p),table_sha256=hashes))
# Native state logs: support all attempts and TF states; failed solves remain present.
states=[];largest={'baryon_mean_error':0.,'charge_mean_error':0.,'euler_error':0.}
for folder in [BASE/'tf_regeneration/states',BASE/'support_repair/states']:
 for p in sorted(folder.rglob('STATE.json')):
  s=load(p)
  rec=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),status=s['status'],actual_exit_code=s.get('actual_exit_code'),accepted=bool(s.get('accepted_stationary')))
  if rec['accepted']:
   profile=np.array(s['profile']);w=np.diff(np.linspace(0.,1.,s['points']+1)**s['dimension'])
   baryon=abs(float(w@(profile[:,0]+profile[:,1]))-s['n_b'])
   charge=abs(float(w@(profile[:,1]-profile[:,2])))
   euler=abs(s['n_b']*s['mu_B']-s['epsilon']-s['pressure'])
   largest={k:max(largest[k],v) for k,v in zip(largest,[baryon,charge,euler])}
   assert np.all(profile[:,:3]>=0) and baryon<1e-10 and charge<1e-10 and euler<1e-12
   assert s['optimizer_success']==1. and s['muon_threshold_crossed']==0.
   assert s['kkt']<.02
   if s['status']=='STATIONARY':assert s['radius_stationarity']<2e-4 and s['radius_bound']==0 and s['density_bound']==0
  states.append(rec)
regen=load(BASE/'support_repair/REGENERATION.json')
point_count=sum(len(r['points_results']) for r in regen['runs'])
accepted=sum(p['accepted'] for r in regen['runs'] for p in r['points_results'])
native_calls=sum(len(p['attempts']) for r in regen['runs'] for p in r['points_results'])
assert (point_count,accepted,native_calls)==(462,453,674)
for r in regen['runs']:
 for p in r['points_results']:
  for a in p['attempts']:assert sha(a['state_path'])==a['state_sha256']
matching=load(BASE/'support_repair/MATCHING.json')
candidates=[c for cell in matching['cells'] for c in cell['candidates']]
assert sum(c.get('kernel_invoked',False) for c in candidates)==29
assert sum(c['status']=='GENERATED_VERIFIED_FINITE_DOMAIN' for c in candidates)==16
# F tracked input data stay untouched.
proc=subprocess.run(['git','diff','--name-only','d89d877d50b6ae3ad4d5cdea29d5ff67e3a7d09a','--','data/new_model_maxwell_20261117','data/research_native_20261116'],cwd=ROOT,text=True,capture_output=True)
assert proc.returncode==0 and proc.stdout.strip()==''
result=dict(status='PASS',scope='Read-only recorded-state physical invariants, manifests and output byte hashes; no solver calls',
  representation_repair={'before_sha256':sha(before),'after_sha256':sha(after),'AST_changed_functions':changed},
  output_tables_checked=len(outputs),outputs=outputs,native_states_checked=len(states),state_statuses={s:sum(x['status']==s for x in states) for s in sorted({x['status'] for x in states})},
  independently_recomputed_profile_residual_maxima=largest,
  support_counts=dict(requested_points=point_count,accepted_points=accepted,native_attempts=native_calls,maxwell_calls=29,admitted_outputs=16),
  unchanged_F_data_git_diff=proc.stdout,
  sources={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),BASE/'support_repair/REGENERATION.json',BASE/'support_repair/MATCHING.json',BASE/'tf_regeneration/PROJECTION_REPAIR.json']})
(BASE/'review/ARTIFACT_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['outputs','sources']},indent=2))
