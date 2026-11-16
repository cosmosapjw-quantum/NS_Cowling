#!/usr/bin/env python3
"""One bounded TM1e half-step precision repair; retain initial tangent evidence."""
import hashlib,json,os,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'data/research_blockers_20261116/tf/common_tangent'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 src=BASE/'TM1e/stdout.txt';lines=[x.split() for x in src.read_text().splitlines()]
 co=next(x for x in lines if x[0]=='COEXISTENCE');dest=BASE/'TM1e_polish';dest.mkdir(exist_ok=False)
 inp=''
 for j,(d,seed,geom) in enumerate(((1,'droplet','slab'),(2,'hole','tube')),start=1):
  n=co[j];states=[x for x in lines if x[0]=='STATE' and int(x[2])==j and float(x[5])==float(n)]
  state=states[0];radius=state[11]
  prof=sorted([x for x in lines if x[0]=='PROFILE' and int(x[1])==j],key=lambda x:int(x[2]))
  assert len(prof)==128
  inp+=f'{d} {seed} {n} {radius} {geom}\n'+' '.join(x[3] for x in prof)+'\n'+' '.join(x[4] for x in prof)+'\n'
 (dest/'input.txt').write_text(inp);exe=ROOT/'build/blocker_tf_native/blocker_tf_common_tangent_polish'
 cmd=[str(exe),'GPPVA_TM1e','128'];meta={'status':'RUNNING','command':cmd,'source_stdout_sha256':sha(src),'driver_sha256':sha(exe),'source_sha256':sha(ROOT/'fortran_eos/tests/blocker_tf_common_tangent_polish.f90'),'input_sha256':sha(dest/'input.txt'),'change':'h 2e-5 -> 1e-5; ftol 1e-13 -> 1e-14; identical physical gates and root targets'}
 (dest/'RECEIPT.json').write_text(json.dumps(meta,indent=2)+'\n');start=time.monotonic()
 with (dest/'stdout.txt').open('w') as out,(dest/'stderr.txt').open('w') as err:
  try:
   p=subprocess.run(cmd,cwd=ROOT,input=inp,text=True,stdout=out,stderr=err,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'},timeout=180)
   meta.update(status='COMPLETED' if p.returncode==0 else 'HOLD_NATIVE_SOLVE',returncode=p.returncode)
  except subprocess.TimeoutExpired:meta['status']='TIMEOUT'
 meta.update(seconds=time.monotonic()-start,stdout_sha256=sha(dest/'stdout.txt'),stderr_sha256=sha(dest/'stderr.txt'));(dest/'RECEIPT.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta))
if __name__=='__main__':main()
