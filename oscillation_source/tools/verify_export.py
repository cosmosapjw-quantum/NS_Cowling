"""Check local dependencies and optionally build and run the exported learner."""
from pathlib import Path
import argparse,json,sys,tempfile,tomllib,hashlib
if __package__ in (None,''):sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from tools.runtime import CaseRunner

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--export',type=Path,required=True);ap.add_argument('--build',action='store_true');ap.add_argument('--output',type=Path,required=True)
 a=ap.parse_args();root=a.export.resolve()
 if a.output.exists():ap.error('Result path exists')
 result={'root':str(root),'build_executed':False}
 try:
  if (root/'solutions').exists():raise ValueError('An independent learner directory must not include solutions')
  for name in ['learner','support']:
   f=root/name/'fpm.toml';d=tomllib.loads(f.read_text())
   for dep in d.get('dependencies',{}).values():
    p=(f.parent/dep['path']).resolve()
    if not p.is_relative_to(root) or not p.is_dir():raise ValueError('Dependency escapes the directory')
  manifest=json.loads((root/'SOURCE_FILES.json').read_text())
  for package,rows in manifest['source_files'].items():
   for rel,sha in rows.items():
    if hashlib.sha256((root/package/rel).read_bytes()).hexdigest()!=sha:raise ValueError('Changed source: '+rel)
  if a.build:
   with tempfile.TemporaryDirectory(prefix='ns-check-') as tmp:
    runner=CaseRunner(root,Path(tmp),'learner');art=runner.build_fpm()
    run=runner.run_artifact(art,'environment',{})
    if run['status']!='ok' or run['workspace']!='learner':raise ValueError('Unexpected environment result')
    result.update(build_executed=True,binary_sha256=art.binary_sha256,source_sha256=art.source_sha256,environment_status=run['status'])
  result['status']='ok';code=0
 except Exception as e:
  result.update(status='failed',error=type(e).__name__+': '+str(e));code=1
 a.output.parent.mkdir(parents=True,exist_ok=True)
 with a.output.open('x') as f:json.dump(result,f,indent=2)
 print(json.dumps(result,indent=2));return code
if __name__=='__main__':raise SystemExit(main())
