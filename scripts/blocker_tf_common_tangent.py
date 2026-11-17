#!/usr/bin/env python3
"""Run two bounded native common-tangent candidates from actual N128 profiles."""
import csv,hashlib,json,os,subprocess,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'data/research_blockers_20261116/tf'
DRIVER=ROOT/'build/blocker_tf_native/blocker_tf_common_tangent'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');t.replace(p)
def run(model):
 dest=BASE/'common_tangent'/model;dest.mkdir(parents=True,exist_ok=False)
 source=BASE/'ambiguous_N128'/model/'output/generated';rows=list(csv.DictReader((source/'selected.csv').open()))
 if model=='TM1e':rows=rows[-3:]
 selected=(rows[0],rows[-1]);attempts=list(csv.DictReader((source/'attempts.csv').open()));profiles=list(csv.DictReader((source/'profiles.csv').open()))
 text=''
 for row in selected:
  a=next(x for x in attempts if (x['node'],x['mesh'],x['attempt'])==(row['node'],row['mesh'],row['attempt']))
  prof=sorted([x for x in profiles if (x['node'],x['mesh'],x['attempt'])==(row['node'],row['mesh'],row['attempt'])],key=lambda x:int(x['cell']))
  assert len(prof)==128
  text+=f"{a['dimension']} {a['seed']} {row['n_b']} {a['radius']} {row['geometry']}\n"
  text+='\n'.join(' '.join(x[key] for x in prof) for key in ('nn','np'))+'\n'
 (dest/'input.txt').write_text(text)
 cmd=[str(DRIVER),'GPPVA_'+model,'128']
 meta={'status':'RUNNING','model':model,'command':cmd,'input_sha256':sha(dest/'input.txt'),'source_files':{str((source/name).relative_to(ROOT)):sha(source/name) for name in ('selected.csv','attempts.csv','profiles.csv')},'driver_sha256':sha(DRIVER),'driver_source_sha256':sha(ROOT/'fortran_eos/tests/blocker_tf_common_tangent.f90'),'scope':'Two separately neutral TF branches; common tangent only after actual native solves, not table extrapolation, not global phase certification'}
 write(dest/'RECEIPT.json',meta);start=time.monotonic()
 try:
  proc=subprocess.run(cmd,input=text.encode(),cwd=ROOT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'},stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=360)
  (dest/'stdout.txt').write_bytes(proc.stdout);(dest/'stderr.txt').write_bytes(proc.stderr);meta.update(status='COMPLETED' if proc.returncode==0 else 'HOLD_NATIVE_SOLVE',returncode=proc.returncode)
 except subprocess.TimeoutExpired as e:
  (dest/'stdout.txt').write_bytes(e.stdout or b'');(dest/'stderr.txt').write_bytes(e.stderr or b'');meta['status']='TIMEOUT'
 meta.update(seconds=time.monotonic()-start,stdout_sha256=sha(dest/'stdout.txt'),stderr_sha256=sha(dest/'stderr.txt'));write(dest/'RECEIPT.json',meta)
 print(json.dumps({k:meta[k] for k in ('model','status','seconds')}),flush=True)
def main():
 with ThreadPoolExecutor(max_workers=2) as p:
  for f in as_completed([p.submit(run,m) for m in ('FSU2R','TM1e')]):f.result()
if __name__=='__main__':main()
