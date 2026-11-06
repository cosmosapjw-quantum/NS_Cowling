"""Copy the current learner and its inputs into a self-contained source directory."""
from pathlib import Path
import argparse,json,hashlib,sys
if __package__ in (None,''):sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from tools.prepare_workspace import create,ROOT
def fingerprint(path):
 return {str(f.relative_to(path)):hashlib.sha256(f.read_bytes()).hexdigest()
  for f in sorted(path.rglob('*')) if f.is_file() and (f.suffix=='.f90' or f.name=='fpm.toml')}
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 try:
  before={k:fingerprint(ROOT/k) for k in ['learner','support']}
  create(a.output,'foundations')
  after={k:fingerprint(a.output/k) for k in ['learner','support']}
  if before!=after:raise ValueError('Source changed during copy')
  result={'source_files':after,'learner_unchanged':True,'support_unchanged':True,'solutions_included':False}
  (a.output/'SOURCE_FILES.json').write_text(json.dumps(result,indent=2)+'\n')
  print(json.dumps({'output':str(a.output),'files':sum(map(len,after.values()))}));return 0
 except (OSError,ValueError) as e:print(str(e),file=sys.stderr);return 2
if __name__=='__main__':raise SystemExit(main())
