from pathlib import Path
import os, sys, subprocess, json, time, hashlib
from datetime import datetime, timezone
root=Path('/home/cosmosapjw/ns_cowling_wt/audit_eos_matching_20261117')
py='/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python'
env=os.environ.copy()
env.update(PYTHONPATH=str(root/'python_frontend')+':/tmp/ns_cowling_audit_20261117_pydeps',COMPOSE_EOS_FORTRAN_ROOT=str(root),COMPOSE_EOS_FORTRAN_EXE=str(root/'build/compose_eos_fortran/compose_eos_fortran'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
steps=[('build',['bash','fortran_eos/build.sh','build/compose_eos_fortran'],{}),('frontend_pipeline',[py,'-m','pytest','-q','-ra','python_frontend/tests','pipeline/tests'],{}),('generator',[py,'-m','pytest','-q','-ra','eos_generator/tests'],{'PYTHONPATH':str(root/'eos_generator')+':'+str(root/'eos_generator/compose_eos_pkg')}),('theory',[py,'-m','pytest','-q','fortran_eos/tests/matching_theory_soft_barotrope.py','fortran_eos/tests/matching_theory_maxwell_topology.py'],{})]
source_names=subprocess.check_output(['git','diff','--name-only'],cwd=root,text=True).splitlines()
source_names += subprocess.check_output(['git','ls-files','--others','--exclude-standard','fortran_eos','python_frontend'],cwd=root,text=True).splitlines()
source_names=sorted(n for n in set(source_names) if n.endswith(('.py','.f90')))
def hashes():return {n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in source_names}
record={'base_sha':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'started_utc':datetime.now(timezone.utc).isoformat(),'cwd':str(root),'python_executable':py,'threads':1,'source_sha256_before':hashes(),'steps':[]}
record_path=Path('/tmp/ns_cowling_audit_20261117_final_execution.json')
for name,command,extra in steps:
 log=Path('/tmp/ns_cowling_audit_20261117_final_'+name+'.log')
 started=time.monotonic()
 with log.open('w') as stream: result=subprocess.run(command,cwd=root,env={**env,**extra},stdout=stream,stderr=subprocess.STDOUT,check=False)
 item={'name':name,'command':command,'environment':{k:({**env,**extra})[k] for k in ['PYTHONPATH','COMPOSE_EOS_FORTRAN_ROOT','COMPOSE_EOS_FORTRAN_EXE','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']},'returncode':result.returncode,'seconds':time.monotonic()-started,'log':str(log),'tail':log.read_text().splitlines()[-15:]}
 record['steps'].append(item); record_path.write_text(json.dumps(record,indent=2)+'\n')
 print(name,'returncode',result.returncode,'seconds',round(item['seconds'],2),flush=True)
 if result.returncode:break
record['source_sha256_after']=hashes()
record['source_unchanged_during_run']=record['source_sha256_before']==record['source_sha256_after']
record['finished_utc']=datetime.now(timezone.utc).isoformat()
record['binary_sha256']={n:hashlib.sha256((root/'build/compose_eos_fortran'/n).read_bytes()).hexdigest() for n in ['compose_eos_fortran','compose_eos_generate']}
record['status']='PASS' if len(record['steps'])==len(steps) and all(s['returncode']==0 for s in record['steps']) and record['source_unchanged_during_run'] else 'FAILED'
record_path.write_text(json.dumps(record,indent=2)+'\n')
print(record['status'],flush=True)
sys.exit(0 if record['status']=='PASS' else 1)
