import argparse,concurrent.futures,hashlib,json,os,subprocess,time
from pathlib import Path
parser=argparse.ArgumentParser();parser.add_argument('--driver',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args();root=Path(__file__).resolve().parents[3];out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
exe=args.driver.resolve()
models=('DD2','DDME2','TW','FSU2','FSU2R','FSU2H','NL3wrL55','TM1e')
geometry=((3,'droplet'),(2,'droplet'),(1,'droplet'),(2,'hole'),(3,'hole'),(3,'uniform'))
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
def run(task):
 model,n,d,seed=task;name=f'{model}_{n:.4f}_d{d}_{seed}';t=time.time()
 command=[str(exe),'family',model,'16',str(d),str(n),'30',seed]
 r=subprocess.run(command,env=env,text=True,capture_output=True)
 (out/(name+'.out')).write_text(r.stdout);(out/(name+'.log')).write_text(r.stderr)
 item={'model':model,'n_b':n,'points':16,'dimension':d,'seed':seed,'returncode':r.returncode,'seconds':time.time()-t,'command':command}
 lines=r.stdout.splitlines();item['status']=lines[0] if lines else 'NO_OUTPUT'
 if r.returncode==0:
  vals=list(map(float,lines[1].split()));item.update(dict(zip(('epsilon','pressure','mu_B','mu_e','scalar_residual','vector_residual','electron_residual','d_epsilon_d_log_radius','kkt','radius_stationarity','Y_p','radius','iterations','evaluations','optimizer_mode','stationary'),vals)))
 item['stdout_sha256']=hashlib.sha256(r.stdout.encode()).hexdigest();item['stderr_sha256']=hashlib.sha256(r.stderr.encode()).hexdigest()
 (out/(name+'.json')).write_text(json.dumps(item,indent=2)+'\n')
 print(name,item['status'],f'{item["seconds"]:.1f}s',flush=True)
 return item
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 results=list(pool.map(run,[(m,n,d,s) for m in models for n in (.003,.03,.075) for d,s in geometry]))
summary={'status':'COMPLETE_WITH_RETAINED_FAILURES' if any(x['returncode'] for x in results) else 'ALL_FAMILIES_STATIONARY','backend':'fortran','n_jobs':len(results),'success':sum(x['returncode']==0 for x in results),'results':results}
(out/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
