"""Fortran execution CLI and byte-bound reader; no pytest runtime dependency."""
import argparse,json,sys
from pathlib import Path
if __package__ in (None,''):sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from tools.atomic_io import read_json,digest
from tools.runtime import CaseRunner,BootstrapError,parse_records
class RunNotComplete(RuntimeError):pass
class RunRecordMismatch(RuntimeError):pass

def read_completed_run(directory,expected_token=None,expected_source_sha256=None):
    p=Path(directory)
    if p.is_symlink():raise RunRecordMismatch('Run symlink')
    try:s=read_json(p/'run.json')
    except (OSError,ValueError) as e:raise RunRecordMismatch('Missing/malformed run manifest') from e
    if not isinstance(s,dict):raise RunRecordMismatch('Invalid manifest')
    if s.get('status')!='complete':raise RunNotComplete('Only complete runs provide success results')
    try:
        if s.get('schema_version')!=2 or s.get('domain_status')!='ok':raise ValueError('Inconsistent completion')
        t=s['request_token']
        if s.get('run_id')!=t or (expected_token is not None and t!=expected_token):raise ValueError('Wrong request identity')
        if expected_source_sha256 is not None and s['source_sha256']!=expected_source_sha256:raise ValueError('Wrong source identity')
        for name,key in [('input.nml','input_sha256'),('result.json','result_sha256'),('process/stdout.txt','stdout_sha256'),('process/command.json','command_sha256')]:
            parts=Path(name).parts
            if any(p.joinpath(*parts[:i]).is_symlink() for i in range(1,len(parts)+1)):raise ValueError('Run artifact symlink')
            if digest(p/name)!=s[key]:raise ValueError('Changed run artifact: '+name)
        cmd=read_json(p/'process/command.json')
        if cmd['returncode']!=0 or cmd['timed_out']:raise ValueError('Incomplete process')
        r=read_json(p/'result.json')
        if r!=parse_records((p/'process/stdout.txt').read_text()):raise ValueError('Result not actual stdout')
        if r.get('request_token')!=t or r.get('workspace')!=s['workspace'] or r.get('case_id')!=s['case_id'] or r.get('status')!='ok':raise ValueError('Result identity/status')
        return r
    except (OSError,ValueError,KeyError,TypeError,BootstrapError) as e:raise RunRecordMismatch(str(e)) from e

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--work-root',type=Path,required=True)
    p.add_argument('--workspace',choices=['solutions','learner'],default='solutions');p.add_argument('--case',required=True);p.add_argument('--parameters',type=Path)
    a=p.parse_args(argv)
    try:
        r=CaseRunner(a.root,a.work_root,a.workspace).run_case(a.case,read_json(a.parameters) if a.parameters else {})
        print(json.dumps(r,ensure_ascii=False,allow_nan=False));return 0 if r['status']=='ok' else 1
    except (OSError,ValueError,BootstrapError) as e:print(type(e).__name__+': '+str(e),file=sys.stderr);return 2
if __name__=='__main__':raise SystemExit(main())
