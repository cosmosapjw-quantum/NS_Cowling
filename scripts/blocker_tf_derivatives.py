#!/usr/bin/env python3
"""Transport selected native profiles to native energy derivative continuation."""
import argparse,csv,hashlib,json,os,subprocess,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'data/research_blockers_20261116/tf'
DRIVER=ROOT/'build/blocker_tf_native/blocker_tf_branch_derivative'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,obj):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');t.replace(p)
def read(path):return list(csv.DictReader(path.open()))
def task(job,out):
 m,r,a,p,src=job;n=float(r['n_b']);mesh=int(r['mesh'])
 ident=f'{m}_N{mesh}_n{n:.12g}';dest=out/ident;dest.mkdir(exist_ok=True)
 receipt=dest/'RECEIPT.json'
 if receipt.exists():
  old=json.loads(receipt.read_text())
  if old['status']!='RUNNING':return old
 inp=' '.join(r[k] for k in ('epsilon','pressure','mu_b'))+'\n'
 inp+='\n'.join(' '.join(row[key] for row in p) for key in ('nn','np'))+'\n'
 (dest/'input.txt').write_text(inp)
 cmd=[str(DRIVER),'GPPVA_'+m,str(mesh),a['dimension'],r['n_b'],a['radius'],a['seed']]
 result={'id':ident,'model':m,'density_fm3':n,'mesh':mesh,'source_selected':str(src.relative_to(ROOT)),'source_selected_sha256':sha(src),'selected':r,'attempt':a,'input_sha256':sha(dest/'input.txt'),'command':cmd,'status':'RUNNING'}
 write(receipt,result);start=time.monotonic()
 try:
  run=subprocess.run(cmd,input=inp.encode(),cwd=ROOT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'},stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180)
  (dest/'stdout.txt').write_bytes(run.stdout);(dest/'stderr.txt').write_bytes(run.stderr)
  result.update(status='COMPLETED' if run.returncode==0 else 'PROCESS_FAILURE',returncode=run.returncode)
 except subprocess.TimeoutExpired as exc:
  (dest/'stdout.txt').write_bytes(exc.stdout or b'');(dest/'stderr.txt').write_bytes(exc.stderr or b'');result['status']='TIMEOUT'
 result['seconds']=time.monotonic()-start
 result['stdout_sha256']=sha(dest/'stdout.txt');result['stderr_sha256']=sha(dest/'stderr.txt');write(receipt,result)
 print(json.dumps({k:result[k] for k in ('id','status','seconds')}),flush=True);return result
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);ap.add_argument('--mesh',type=int,default=64);ap.add_argument('--ambiguous',action='store_true');args=ap.parse_args()
 assert 1<=args.workers<=4
 out=BASE/f'branch_derivative_N{args.mesh}';out.mkdir(exist_ok=True)
 plan=json.loads((BASE/('ambiguous_N128/PLAN.json' if args.ambiguous else 'PLAN.json')).read_text());jobs=[]
 for m,v in plan['models'].items():
  wanted=set(map(float,v['densities'] if args.ambiguous else v['endpoint_union']));seen=set()
  sources=(BASE/'ambiguous_N128'/m/'output/generated/selected.csv',) if args.ambiguous else (ROOT/v['reuse_selected_path'],BASE/m/'output/generated/selected.csv')
  for src in sources:
   rows=read(src);attempts=read(src.parent/'attempts.csv');profiles=read(src.parent/'profiles.csv')
   amap={(r['node'],r['mesh'],r['attempt']):r for r in attempts}
   pmap={}
   for row in profiles:pmap.setdefault((row['node'],row['mesh'],row['attempt']),[]).append(row)
   for r in rows:
    n=float(r['n_b'])
    if int(r['mesh'])!=args.mesh or n not in wanted or n in seen:continue
    key=(r['node'],r['mesh'],r['attempt']);p=sorted(pmap[key],key=lambda x:int(x['cell']))
    assert len(p)==args.mesh
    jobs.append((m,r,amap[key],p,src));seen.add(n)
  assert seen==wanted,(m,len(seen),len(wanted))
 summary={'status':'RUNNING','driver_sha256':sha(DRIVER),'driver_source_sha256':sha(ROOT/'fortran_eos/tests/blocker_tf_branch_derivative.f90'),'jobs':len(jobs),'workers':args.workers,'threads_per_child':1,'h_fm3':'min(1e-5,1e-3*n)','ftol':1e-13,'physical_tolerances':'unchanged native KKT and stationarity gates; raw derivative residuals reported','results':[]}
 write(out/'SUMMARY.json',summary)
 with ThreadPoolExecutor(max_workers=args.workers) as pool:
  for f in as_completed([pool.submit(task,j,out) for j in jobs]):
   r=f.result();summary['results'].append({k:r[k] for k in ('id','status','seconds')});write(out/'SUMMARY.json',summary)
 summary['status']='COMPLETED' if all(r['status']=='COMPLETED' for r in summary['results']) else 'COMPLETED_WITH_FAILURES';write(out/'SUMMARY.json',summary)
if __name__=='__main__':main()
