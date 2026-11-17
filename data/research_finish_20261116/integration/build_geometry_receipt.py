#!/usr/bin/env python3
"""Bind completed finite geometry campaign; no physical solves are performed."""
from pathlib import Path
import hashlib, json
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'data/research_finish_20261116'
OUT=DATA/'geometry/RECEIPT.json'
def read(p): return json.loads((ROOT/p).read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
review=read('data/research_finish_20261116/independent_review/GEOMETRY_UNION_RESUME_REVIEW.json')
assert review['status']=='PASS_SCOPED'
assert (review['unique_nodes'],review['selected'],review['native_tables'],review['native_rows'])==(156,156,18,158)
assert review['unresolved_families']==review['lower_unresolved']==0
reg=read('data/research_finish_20261116/integration/resume_verification/REGRESSION.json')
assert reg['status']=='PASS' and reg['tested_sources_unchanged']
for name,digest in reg['tested_source_sha256_after'].items():
 assert sha(ROOT/name)==digest,name
evidence=[
'geometry/exported/RECEIPT.json','geometry/nl3_upper_exported/RECEIPT.json',
'geometry/retained_revalidation/RECEIPT.json','geometry/OBSERVED_PHASES.json',
'independent_review/REVIEW_RECEIPT.json','independent_review/GEOMETRY_UNION_RESUME_REVIEW.json',
'integration/resume_verification/REGRESSION.json','integration/ENVIRONMENT.json']
bindings={str((DATA/p).relative_to(ROOT)):sha(DATA/p) for p in evidence}
# Raw ledgers, source snapshots and native/profile products remain immutable.
geometry_files=[p for p in (DATA/'geometry').rglob('*') if p.is_file()
 and '__pycache__' not in p.parts and p!=OUT and p.suffix not in ('.pyc','.pyo')]
result={
 'schema':'ns-cowling.finite-geometry-closeout.v1','status':'PASS_SCOPED',
 'base_commit':'573a31485bc9af30915e01fe74f259e4fa5ac54f',
 'primary_mesh_density_nodes':152,'upper_extension_new_nodes':4,
 'summary':review,'model_count':8,'meshes':[192,384],
 'retained_dense_nodes':112,'new_band_nodes':40,
 'retained_profile_unique_count':1241,'retained_profile_attempt_count':1254,
 'raw_failure_history_preserved':True,'result_validity':'finite-selection/readback checked',
 'provenance_validity':'source/input/profile/native identities checked',
 'global_minimum_proven':False,'continuum_transition_proven':False,
 'scope':'Completed declared finite campaign only; no original-model full-crust equivalence implied.',
 'binary_sha256':reg['binary_sha256'],'source_sha256':{
 n:d for n,d in reg['tested_source_sha256_after'].items()
 if n.startswith('python_frontend/compose_eos/') or n in [
 'scripts/complete_gppva_geometry.py','scripts/refine_gppva_geometry.py',
 'scripts/export_gppva_geometry.py','scripts/summarize_gppva_geometry.py']},
 'evidence_sha256':bindings,
 'geometry_artifact_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(geometry_files)}}
OUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({'status':result['status'],'unique_nodes':156,'bound_geometry_files':len(geometry_files),'receipt_sha256':sha(OUT)}))
