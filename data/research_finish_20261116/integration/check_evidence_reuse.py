from pathlib import Path
import json,subprocess,hashlib
root=Path(__file__).resolve().parents[3]
out=root/'data/research_finish_20261116/integration'
base='573a31485bc9af30915e01fe74f259e4fa5ac54f'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
changed=set(subprocess.check_output(['git','diff','--name-only',base],cwd=root,text=True).splitlines())
results={}
for file in ('SHA256SUMS_v4.txt','eos_generator/SHA256SUMS_v3_3.txt'):
 raw=subprocess.check_output(['git','show',base+':'+file],cwd=root,text=True)
 good=[];intentional=[];unexpected=[]
 for line in raw.splitlines():
  if not line.strip():continue
  expected,name=line.split(maxsplit=1);name=name.lstrip('*')
  parent=root/Path(file).parent
  p=parent/name
  if not p.exists():p=root/name
  actual=digest(p) if p.is_file() else None
  rel=str(p.relative_to(root))
  if actual==expected:good.append(rel)
  elif rel in changed:intentional.append({'path':rel,'before':expected,'after':actual})
  else:unexpected.append({'path':rel,'expected':expected,'actual':actual})
 results[file]={'unchanged_verified':len(good),'intentional_modified_sources':intentional,'unexpected_mismatches':unexpected}
cores={}
for p in sorted((root/'data/research_completion_20261116/core/generated').glob('*/RECEIPT.json')):
 r=json.loads(p.read_text()); prov=r['provenance']
 deps={n:digest(root/'python_frontend/compose_eos'/n)==v for n,v in prov['code_sha256'].items()}
 outputs={n:digest(p.parent/n)==v for n,v in r['output_sha256'].items()}
 cores[p.parent.name]={'source_dependencies_match_current':all(deps.values()),'dependency_checks':deps,
  'prior_native_output_hashes_match':all(outputs.values()),'rows':len(r['rows']),
  'evidence':str(p.relative_to(root))}
result={'base_commit':base,'historical_manifests':results,'reused_core_evidence':cores,
 'scope':'Historical unaffected core science is reused; this does not describe a new full core campaign.',
 'status':'PASS' if all(not r['unexpected_mismatches'] for r in results.values()) and all(r['source_dependencies_match_current'] and r['prior_native_output_hashes_match'] for r in cores.values()) else 'FAIL'}
(out/'EVIDENCE_REUSE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
