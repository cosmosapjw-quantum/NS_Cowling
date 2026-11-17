#!/usr/bin/env python3
"""Run untouched native generator on exact remaining negative-interval endpoints."""
import argparse,csv,hashlib,json,os,subprocess,tarfile,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/research_blockers_20261116/tf'
REUSED=ROOT/'build/blocker_tf_reused'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');t.replace(p)
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 if (OUT/'PLAN.json').exists():raise FileExistsError('Preserved campaign exists; use its checkpoint, do not overwrite')
 old=ROOT/'data/research_native_20261116/source_grid_cohort/SUMMARY.json'
 deps=json.loads(old.read_text())['source_sha256']; mismatches=[p for p,h in deps.items() if sha(ROOT/p)!=h]
 # Historical runner/source-grid copy can differ independently; every producer dependency remains bound.
 producer={p:h for p,h in deps.items() if p.startswith(('fortran_eos/','python_frontend/'))}
 bad=[p for p,h in producer.items() if sha(ROOT/p)!=h]
 write(OUT/'SOURCE_IDENTITY.json',{'baseline':'a49aaf9f82ddfda484e0f965aefdfe73aa7918f2','prior_receipt':str(old.relative_to(ROOT)),'prior_dependency_count':len(deps),'all_dependency_mismatches':mismatches,'producer_dependency_count':len(producer),'producer_mismatches':bad,'status':'PASS' if not bad else 'HOLD'})
 if bad: raise RuntimeError(bad)
 archive=ROOT/'data/research_native_refinement_20261116/native_mesh_runs.tar.gz'
 if not REUSED.exists():
  REUSED.mkdir(parents=True)
  with tarfile.open(archive) as t:t.extractall(REUSED,filter='data')
 plan={'baseline':'a49aaf9f82ddfda484e0f965aefdfe73aa7918f2','status':'PREPARED','threads_per_child':1,'workers':3,'existing_failure_history':'data/research_native_20261116/source_grid_cohort','reused_archive_sha256':sha(archive),'models':{}}
 for m in ('FSU2R','FSU2H','TM1e'):
  source=ROOT/f'data/research_native_20261116/source_grid_cohort/GPPVA_{m}/attempt_001/output/generated/selected.csv'
  rows=list(csv.DictReader(source.open()))
  intervals=[{'left':a['n_b'].strip(),'right':b['n_b'].strip(),'delta_P_N16':float(b['pressure'])-float(a['pressure'])} for a,b in zip(rows,rows[1:]) if float(b['pressure'])<float(a['pressure'])]
  worst=min(intervals,key=lambda x:x['delta_P_N16']);remaining=[x for x in intervals if x!=worst]
  endpoint=sorted(set([x[k] for x in remaining for k in ('left','right')]),key=float)
  prior=list(REUSED.glob(f'**/native_refine_{m}/generated/selected.csv'))
  assert len(prior)==1,prior
  oldrows=list(csv.DictReader(prior[0].open()))
  reused={r['n_b'].strip() for r in oldrows if int(r['mesh']) in (32,64)}
  todo=[n for n in endpoint if float(n) not in {float(v) for v in reused}]
  target=OUT/m
  cmd=['python3','-m','compose_eos.cli','generate','GPPVA_'+m,str(target/'output'),'--repository-root',str(ROOT),'--densities',*todo,'--resolutions','32','64']
  plan['models'][m]={'original_negative_intervals':len(intervals),'remaining_intervals':remaining,'worst_reused':worst,'endpoint_union':endpoint,'reused_endpoint_values':[n for n in endpoint if float(n) in {float(v) for v in reused}],'new_density_values':todo,'new_mesh_density_nodes':2*len(todo),'original_selected_sha256':sha(source),'reuse_selected_path':str(prior[0].relative_to(ROOT)),'reuse_selected_sha256':sha(prior[0]),'command':cmd}
 write(OUT/'PLAN.json',plan)
 env={**os.environ,'PYTHONPATH':str(ROOT/'python_frontend'),'COMPOSE_EOS_FORTRAN_ROOT':str(ROOT),'COMPOSE_EOS_NATIVE_BUILD_DIR':str(ROOT/'build/blocker_tf_native'),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
 # Replays use the supported launcher build-directory setting.
 build=['bash','fortran_eos/build.sh',str(ROOT/'build/blocker_tf_native')]
 with (OUT/'build.log').open('w') as log: r=subprocess.run(build,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError('build failed')
 def run(m,p):
  target=OUT/m;target.mkdir(exist_ok=True);started=time.monotonic()
  write(target/'RECEIPT.json',{'status':'RUNNING','command':p['command']})
  with (target/'stdout.log').open('w') as log,(target/'stderr.log').open('w') as err:
   proc=subprocess.run(p['command'],cwd=ROOT,env=env,stdout=log,stderr=err)
  receipt={'status':'COMPLETED' if proc.returncode==0 else 'PROCESS_FAILURE','returncode':proc.returncode,'seconds':time.monotonic()-started,'command':p['command']}
  write(target/'RECEIPT.json',receipt);return m,receipt
 plan['status']='RUNNING';write(OUT/'PLAN.json',plan)
 with ThreadPoolExecutor(max_workers=3) as pool:
  for future in as_completed([pool.submit(run,m,p) for m,p in plan['models'].items()]):
   m,r=future.result();print(m,json.dumps(r),flush=True)
 plan['status']='EXECUTIONS_FINISHED';write(OUT/'PLAN.json',plan)
if __name__=='__main__':main()
