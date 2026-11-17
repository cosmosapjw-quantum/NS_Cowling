#!/usr/bin/env python3
"""Select only unresolved N64 pressure-drop intervals for native N128 midpoint runs."""
import csv,json,os,subprocess,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'data/research_blockers_20261116/tf'
def write(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2)+'\n');t.replace(p)
def main():
 original=json.loads((BASE/'PLAN.json').read_text());plan={'status':'PREPARED','meshes':[128],'models':{}}
 for m,v in original['models'].items():
  rows=[]
  for src in (ROOT/v['reuse_selected_path'],BASE/m/'output/generated/selected.csv'):rows+=list(csv.DictReader(src.open()))
  selected={float(r['n_b']):r for r in rows if int(r['mesh'])==64}
  intervals=[]
  for k in v['remaining_intervals']:
   left,right=float(k['left']),float(k['right'])
   delta=float(selected[right]['pressure'])-float(selected[left]['pressure'])
   if delta<0:intervals.append({**k,'delta_P_N64':delta,'midpoint':(left+right)/2})
  if not intervals:continue
  densities=sorted({float(k[q]) for k in intervals for q in ('left','right','midpoint')})
  target=BASE/'ambiguous_N128'/m
  cmd=['python3','-m','compose_eos.cli','generate','GPPVA_'+m,str(target/'output'),'--repository-root',str(ROOT),'--densities',*[format(n,'.17g') for n in densities],'--resolutions','128']
  plan['models'][m]={'intervals':intervals,'densities':densities,'command':cmd}
 dest=BASE/'ambiguous_N128/PLAN.json'
 if dest.exists():raise FileExistsError(dest)
 write(dest,plan)
 env={**os.environ,'PYTHONPATH':str(ROOT/'python_frontend'),'COMPOSE_EOS_NATIVE_BUILD_DIR':str(ROOT/'build/blocker_tf_native'),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
 def run(m,p):
  d=BASE/'ambiguous_N128'/m;d.mkdir(exist_ok=True);start=time.monotonic()
  with (d/'stdout.log').open('w') as out,(d/'stderr.log').open('w') as err:
   r=subprocess.run(p['command'],cwd=ROOT,env=env,stdout=out,stderr=err)
  result={'status':'COMPLETED' if r.returncode==0 else 'PROCESS_FAILURE','returncode':r.returncode,'seconds':time.monotonic()-start,'command':p['command']}
  write(d/'RECEIPT.json',result);return m,result
 with ThreadPoolExecutor(max_workers=2) as pool:
  for f in as_completed([pool.submit(run,m,p) for m,p in plan['models'].items()]):print(f.result(),flush=True)
 plan['status']='EXECUTIONS_FINISHED';write(dest,plan)
if __name__=='__main__':main()
