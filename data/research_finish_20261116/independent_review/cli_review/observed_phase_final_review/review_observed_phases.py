#!/usr/bin/env python3
"""Read-only independent audit of an observed-phase finite-node report."""
from pathlib import Path
import contextlib,gzip,hashlib,importlib.util,io,json,math,sys,tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[5]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from compose_eos.tf_selection import select_geometry_candidates,feasibility
from compose_eos.tf_morphology import geometry_label
from refine_gppva_geometry import FAMILIES
OUT=Path(__file__).parent
REPORT=ROOT/(sys.argv[1] if len(sys.argv)>1 else 'data/research_finish_20261116/geometry/OBSERVED_PHASE_PREVIEW.json')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
d=json.loads(REPORT.read_text());report_sha=sha(REPORT)
inputs={p:sha(ROOT/p) for p in d['input_sha256']}
assert inputs==d['input_sha256']
sources={p:sha(ROOT/p) for p in d['source_sha256']}
assert sources==d['source_sha256']
assert 'scripts/refine_gppva_tf.py' in sources
loaded={p:json.loads(gzip.decompress((ROOT/p).read_bytes())) for p in inputs}
raw={};records_by_model={}
for p,rows in loaded.items():
 model=Path(p).parent.name
 for r in rows:
  key=(model,round(r['n_b'],12),r['points'])
  if key in raw:assert raw[key]['attempts']==r['attempts']
  raw[key]=r
assert len(raw)==d['unique_mesh_density_nodes']
checks={'sequences':0,'selected':0,'hold':0,'windows':0,'samples':0,'mesh_pairs':0,'raw_sign_brackets':0}
summary=[]
def independent_branches(record):
 out={}
 if record is None:return out
 for r in record['attempts']:
  if r.get('numerical_status') not in ('STATIONARY','UNIFORM_STATIONARY') or r.get('muon_threshold_crossed'):continue
  if feasibility(r,r['n_b'])!='FEASIBLE':continue
  assert math.isfinite(r['epsilon']) and math.isfinite(r['energy_per_baryon'])
  g=geometry_label(r)
  if g not in out or r['epsilon']<out[g]['epsilon']:out[g]=r
 return out
for model in d['models']:
 name=model['model_id'];expected_windows=set();meshes=sorted(map(int,model['sequences']));fine=str(meshes[-1]);coarse=str(meshes[-2]) if len(meshes)>1 else None
 for mesh,seq in model['sequences'].items():
  assert [x['n_b'] for x in seq]==sorted(x['n_b'] for x in seq)
  expected_nodes=sorted(r['n_b'] for k,r in raw.items() if k[0]==name and k[2]==int(mesh))
  assert [x['n_b'] for x in seq]==expected_nodes
  for item in seq:
   checks['sequences']+=1;r=raw[(name,round(item['n_b'],12),int(mesh))]
   for origin in item['evidence']:
    origin_row=loaded[origin['file']][origin['record_index']]
    # An extension anchor adds a provenance pointer but retains identical science.
    assert {k:v for k,v in origin_row.items() if k!='extension_anchor_provenance'}=={k:v for k,v in r.items() if k!='extension_anchor_provenance'}
    for anchor in (origin_row,r):
     if 'extension_anchor_provenance' in anchor:
      pointer=anchor['extension_anchor_provenance']
      assert sha(ROOT/pointer['file'])==pointer['sha256']
      original=loaded[pointer['file']][pointer['record_index']]
      assert {k:v for k,v in anchor.items() if k!='extension_anchor_provenance'}==original
   decision=select_geometry_candidates(r['attempts'],FAMILIES.values(),n_b=r['n_b']);q=decision['selected']
   assert item['status']==decision['status'] and item['best_attempt_index']==decision['best_attempt_index']
   assert item['geometry']==(geometry_label(q) if q else None)
   if q:
    checks['selected']+=1
    for k in ['energy_per_baryon','pressure','mu_B','Y_n','Y_p','Y_e','Y_mu','max_chemical_KKT_residual_MeV','radius_fm']:assert item[k]==q.get(k)
    valid=[x['attempt'] for x in decision['attempt_diagnostics'] if x['stationary_supported']]
    assert q['energy_per_baryon']==min(r['attempts'][i]['energy_per_baryon'] for i in valid)
   else:checks['hold']+=1
  for a,b in zip(seq,seq[1:]):
   if a['geometry'] and b['geometry'] and a['geometry']!=b['geometry']:expected_windows.add((a['n_b'],b['n_b'],tuple((a['geometry'],b['geometry']))))
 actual={(w['left'],w['right'],tuple(w['geometries'])) for w in model['observed_changes']}
 assert actual==expected_windows and len(actual)==len(model['observed_changes'])
 model_brackets=[]
 for w in model['observed_changes']:
  checks['windows']+=1
  assert 'no bound on residual energy error' in w['optimizer_uncertainty'] and 'do not certify an energy-order sign' in w['optimizer_uncertainty']
  expected_ns=sorted({r['n_b'] for k,r in raw.items() if k[0]==name and w['left']<=r['n_b']<=w['right']})
  assert [s['n_b'] for s in w['samples']]==expected_ns
  assert w['simple_physical_morphologies']==all(g in ('droplet','rod','slab','tube','bubble','uniform') for g in w['geometries'])
  for s in w['samples']:
   checks['samples']+=1;expected_values={}
   for mesh in meshes:
    r=raw.get((name,round(s['n_b'],12),mesh));branches=independent_branches(r)
    if all(g in branches for g in w['geometries']):
     pair=[branches[g] for g in w['geometries']]
     expected_values[str(mesh)]={'delta_E_MeV_per_baryon':pair[1]['energy_per_baryon']-pair[0]['energy_per_baryon'],'geometries':w['geometries'],'energies_per_baryon':[x['energy_per_baryon'] for x in pair],'pressures':[x['pressure'] for x in pair],'mu_B':[x['mu_B'] for x in pair]}
   for mesh,value in expected_values.items():
    pair=independent_branches(raw[(name,round(s['n_b'],12),int(mesh))])
    value['stationarity_diagnostics']=[{k:pair[g].get(k) for k in ['max_chemical_KKT_residual_MeV','radius_stationarity_MeV','scalar_pde_residual','vector_pde_residual','electron_pde_residual','baryon_constraint_error','charge_constraint_error','radius_fm','dimension']} for g in w['geometries']]
   assert s['mesh_values']==expected_values;checks['mesh_pairs']+=len(expected_values)
   f=expected_values.get(fine);c=expected_values.get(coarse)
   delta=f['delta_E_MeV_per_baryon'] if f else None
   shift=abs(delta-c['delta_E_MeV_per_baryon']) if f and c else None
   floor=64*np.finfo(float).eps*max(1.,*(abs(x) for x in f['energies_per_baryon'])) if f else None
   assert s['adjacent_mesh_deltaE_change_MeV_per_baryon']==shift
   assert s['roundoff_scale_MeV_per_baryon']==floor
   assert s['roundoff_near_degeneracy']==bool(f and abs(delta)<=floor)
   assert s['fine_sign_exceeds_mesh_change']==bool(f and shift is not None and abs(delta)>shift)
  brackets=[]
  for a,b in zip(w['samples'],w['samples'][1:]):
   av=a['mesh_values'].get(fine);bv=b['mesh_values'].get(fine)
   if av and bv and av['delta_E_MeV_per_baryon']*bv['delta_E_MeV_per_baryon']<0:
    brackets.append({'n_left':a['n_b'],'n_right':b['n_b'],'both_signs_exceed_mesh_change':a['fine_sign_exceeds_mesh_change'] and b['fine_sign_exceeds_mesh_change'],'both_endpoints_exceed_roundoff_scale':not a['roundoff_near_degeneracy'] and not b['roundoff_near_degeneracy']})
  assert brackets==w['fine_mesh_sign_change_brackets'];checks['raw_sign_brackets']+=len(brackets)
  if brackets:model_brackets.append({'geometries':w['geometries'],'brackets':brackets})
  assert 'not automatically a continuum transition' in w['interpretation']
 summary.append({'model':name,'node_count':sum(len(x) for x in model['sequences'].values()),'changes':len(model['observed_changes']),'sign_brackets':model_brackets})
# Isolate the new report's adjacency/union logic with synthetic selectors; actual
# scientific selector and morphology are independently replayed above on all rows.
script=ROOT/'data/research_finish_20261116/geometry/observed_phase_report.py'
spec=importlib.util.spec_from_file_location('observed_phase_audit_target',script);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
mod.select_geometry_candidates=lambda attempts,*a,**kw:{'selected':attempts[0] if attempts else None,'status':'SYNTHETIC_OK' if attempts else 'SYNTHETIC_HOLD','best_attempt_index':0 if attempts else None}
mod.geometry_label=lambda q:q['label']
mod.transition_evidence=lambda windows,records,meshes:windows
synthetic={}
with tempfile.TemporaryDirectory(prefix='phase_adjacency_',dir=OUT) as tmp:
 base=Path(tmp);p=base/'primary'/'SYNTHETIC';p.mkdir(parents=True);dest=base/'report.json'
 records=[{'n_b':.01,'points':192,'attempts':[{'label':'rod'}]},{'n_b':.02,'points':192,'attempts':[]},{'n_b':.03,'points':192,'attempts':[{'label':'slab'}]}]
 def replay(rows,extension=False):
  (p/'attempts_N192.json.gz').write_bytes(gzip.compress(json.dumps(rows).encode()))
  argv=sys.argv;sys.argv=[str(script),str(p.parent),str(dest)]+(['--extension',str(p.parent)] if extension else [])
  try:
   with contextlib.redirect_stdout(io.StringIO()):mod.main()
   return json.loads(dest.read_text())
  finally:sys.argv=argv
 no_bridge=replay(records);assert not no_bridge['models'][0]['observed_changes'];synthetic['hold_node_blocks_adjacent_change']=True
 duplicate=replay(records,True);assert duplicate['unique_mesh_density_nodes']==3;synthetic['identical_union_nodes_deduplicate']=True
 conflicting=records+[{'n_b':.01,'points':192,'attempts':[{'label':'slab'}]}]
 try:replay(conflicting)
 except ValueError as e:assert 'conflicting repeated anchor attempts' in str(e);synthetic['conflicting_union_nodes_rejected']=True
 else:raise AssertionError('conflicting duplicate accepted')
assert all(sha(ROOT/p)==h for p,h in inputs.items())
assert all(sha(ROOT/p)==h for p,h in sources.items()) and sha(REPORT)==report_sha
result={'status':'PASS_SCOPED','scope':'Actual descriptive phase discovery, canonical selection replay, finite two-mesh branch arithmetic and adjacency only; no continuum, common-tangent, original calibration or optimizer energy-error certification.','report':str(REPORT.relative_to(ROOT)),'report_sha256':report_sha,'input_sha256':inputs,'source_sha256':sources,'checks':checks,'models':summary,'synthetic':synthetic,'source_and_input_unchanged':True,'runner_sha256':sha(Path(__file__))}
(OUT/'OBSERVED_PHASE_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'checks':checks,'synthetic':synthetic,'report_sha256':report_sha},indent=2))
