#!/usr/bin/env python3
"""Seal the completed independent review; refuses incomplete geometry evidence."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[3];BASE=Path(__file__).parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
audits=[json.loads((BASE/name).read_text()) for name in ['GEOMETRY_EXPORT_REVIEW.json','GEOMETRY_EXTENSION_REVIEW.json']]
assert all(a['status']=='PASS_SCOPED' for a in audits)
unique={};duplicate_rows=0
for audit in audits:
 assert 0<=audit['native_rows']<=audit['selected_nodes']
 for row in audit['records']:
  key=(row['model'],row['points'],row['n_b'])
  if key in unique:
   assert unique[key]==row
   duplicate_rows+=1
  else:unique[key]=row
assert len(unique)==156 and duplicate_rows==2
assert len({key[0] for key in unique})==8
geometry_summary={
 'model_count':8,'mesh_density_nodes':len(unique),
 'selected_nodes':sum(row['selected'] for row in unique.values()),
 'attempt_count':sum(row['attempt_count'] for row in unique.values()),
 'native_table_count':sum(a['native_table_count'] for a in audits),
 'native_rows_including_repeated_anchors':sum(a['native_rows'] for a in audits),
 'repeated_anchor_rows':duplicate_rows}
# Correctly preserved HOLDs do not invalidate a scoped evidence review.
report=ROOT/'docs/research_finish_20261116/INDEPENDENT_REVIEW.md'
assert 'Status: **PASS_SCOPED**' in report.read_text()
source_paths=[
 'python_frontend/compose_eos/spatial_tf.py',
 'python_frontend/compose_eos/spatial_tf_equilibrium.py',
 'python_frontend/compose_eos/candidate_gdfm_cld.py',
 'python_frontend/compose_eos/spg_crust_generator.py',
 'python_frontend/compose_eos/restored_spg_candidate.py',
 'python_frontend/compose_eos/crust_assembly.py',
 'python_frontend/compose_eos/tf_selection.py',
 'python_frontend/compose_eos/tf_morphology.py',
 'python_frontend/compose_eos/generation_dispatch.py',
 'scripts/build_gppva_tf.py']
files=[p for p in BASE.rglob('*') if p.is_file() and '__pycache__' not in p.parts
       and p.name not in ('REVIEW_RECEIPT.json','SHA256SUMS.txt')]
receipt={
 'status':'PASS_SCOPED','project_status':'PARTIALLY_COMPLETE',
 'base_commit':'573a31485bc9af30915e01fe74f259e4fa5ac54f',
 'independent_production_edits':False,
 'scope':'Electron constrained mathematics and independent numerical oracles; explicit SPG derivatives and native composition; finite geometry selection/ledger/native exports; partial crust energy references and critical-cone/cache defects; public CLI evidence; vector-band algebra and actual-profile equivalence.',
 'excluded_claims':['original-author SPG calibration','full source-equivalent crust','global minimum','continuum phase transition','common-tangent interpretation of fixed-density geometry crossings','complete composition thermochemistry of the BSk representation','core-crust joining','stellar calculations'],
 'final_geometry':geometry_summary,
 'source_sha256':{p:sha(ROOT/p) for p in source_paths},
 'report_sha256':{str(report.relative_to(ROOT)):sha(report)},
 'evidence_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(files)}}
(BASE/'REVIEW_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
files.append(BASE/'REVIEW_RECEIPT.json');files.append(report)
manifest=BASE/'SHA256SUMS.txt'
manifest.write_text(''.join(sha(p)+'  '+str(p.relative_to(ROOT))+'\n' for p in sorted(files)))
for line in manifest.read_text().splitlines():
 digest,name=line.split('  ',1);assert sha(ROOT/name)==digest
print(json.dumps({'status':'PASS_SCOPED','verified_checksum_entries':len(files),
 'source_files':len(source_paths),'final_geometry':receipt['final_geometry']},indent=2))
