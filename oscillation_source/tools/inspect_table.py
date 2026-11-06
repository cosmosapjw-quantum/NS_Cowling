#!/usr/bin/env python3
"""Call the supplied Fortran table reader and native-Hermite admissibility method."""
from __future__ import annotations
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def inspect_table(table:Path,output:Path)->dict:
    table=Path(table).resolve();output=Path(output).resolve()
    if output.exists():raise FileExistsError('Output exists: '+str(output))
    if not table.is_file():raise FileNotFoundError(table)
    if output==ROOT or output.is_relative_to(ROOT):raise ValueError('Use an output outside the source tree')
    fpm=os.environ.get('NS_FPM') or shutil.which('fpm')
    if not fpm or not Path(fpm).is_file():raise FileNotFoundError('fpm not found; set NS_FPM')
    output.mkdir(parents=True);project=output/'call';(project/'app').mkdir(parents=True)
    driver=ROOT/'examples/method_calls/table_query.f90'
    shutil.copy2(driver,project/'app/table_query.f90')
    (project/'fpm.toml').write_text('\n'.join([
        'name = "table-call"','[build]','auto-executables = false',
        'auto-tests = false','auto-examples = false','[dependencies]',
        'ns-methods = { path = '+json.dumps(os.path.relpath(ROOT/'methods',project))+' }',
        'ns-support = { path = '+json.dumps(os.path.relpath(ROOT/'support',project))+' }',
        '[[executable]]','name = "table_query"','source-dir = "app"',
        'main = "table_query.f90"','']))
    identity=hashlib.sha256()
    for base in ('methods','support'):
        for p in sorted([ROOT/base/'fpm.toml',*(ROOT/base/'src').rglob('*.f90')]):
            identity.update(str(p.relative_to(ROOT)).encode());identity.update(p.read_bytes())
    receipt={'input_sha256':hashlib.sha256(table.read_bytes()).hexdigest(),
             'methods_sha256':identity.hexdigest(),
             'driver_sha256':hashlib.sha256(driver.read_bytes()).hexdigest(),
             'fpm_sha256':hashlib.sha256(Path(fpm).read_bytes()).hexdigest(),'status':'running'}
    rp=output/'execution.json';rp.write_text(json.dumps(receipt,indent=2)+'\n')
    cmd=[fpm,'run','table_query','--profile','release','--',str(table),str(output/'table.json')]
    receipt['command']=cmd
    try:
        run=subprocess.run(cmd,cwd=project,capture_output=True,text=True,timeout=120)
        (output/'stdout.log').write_text(run.stdout);(output/'stderr.log').write_text(run.stderr)
        receipt['returncode']=run.returncode
        if run.returncode:raise RuntimeError('Fortran table inspection failed; see '+str(output))
        result=json.loads((output/'table.json').read_text())
        receipt['result_sha256']=hashlib.sha256((output/'table.json').read_bytes()).hexdigest()
        receipt['status']='complete'
    except (OSError,subprocess.SubprocessError,ValueError,RuntimeError) as exc:
        receipt.update(status='failed',error=str(exc));raise
    finally:
        rp.write_text(json.dumps(receipt,indent=2)+'\n')
    return result

def main()->int:
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('table',type=Path)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:
        result=inspect_table(a.table,a.output);print(json.dumps(result,indent=2))
        return 0 if result['model_ready'] else 2
    except (OSError,ValueError,RuntimeError,subprocess.SubprocessError) as exc:
        print(str(exc),file=sys.stderr);return 1
if __name__=='__main__':raise SystemExit(main())
