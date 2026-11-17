"""Final native build and regression evidence; no production mutation."""
import datetime,hashlib,json,os,platform,re,subprocess,time
from pathlib import Path
root=Path(__file__).resolve().parents[3]
out=Path(__file__).resolve().parent
os.chdir(root)
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
         PYTHONPATH='python_frontend:scripts',COMPOSE_EOS_FORTRAN_ROOT=str(root),
         COMPOSE_EOS_FORTRAN_EXE=str(root/'build/compose_eos_fortran/compose_eos_fortran'))
folders=('python_frontend','pipeline','eos_generator','scripts','fortran_eos','.github/workflows')
def snapshot():
 paths=[p for d in folders for p in (root/d).rglob('*') if p.is_file() and p.suffix in ('.py','.f90','.f','.sh','.toml','.yaml','.yml') and '__pycache__' not in p.parts]
 paths+=[root/'pyproject.toml']
 return {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(paths))}
receipt={'status':'RUNNING','base_head':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
         'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'environment':{k:env[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','PYTHONPATH','COMPOSE_EOS_FORTRAN_ROOT','COMPOSE_EOS_FORTRAN_EXE')},
         'python':platform.python_version(),'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
         'source_sha256_before':snapshot(),'steps':[]}
def save():
 (out/'REGRESSION.json').write_text(json.dumps(receipt,indent=2)+'\n')
def run(name,cmd,extra=None):
 t=time.monotonic();start=datetime.datetime.now(datetime.timezone.utc).isoformat()
 with (out/(name+'.log')).open('w') as log:
  r=subprocess.run(cmd,env=dict(env,**(extra or {})),stdout=log,stderr=subprocess.STDOUT)
 text=(out/(name+'.log')).read_text(errors='replace')
 summary=[s for s in text.splitlines() if re.search(r'\d+ (passed|failed|skipped)|PASS|FAIL',s)]
 step={'name':name,'command':cmd,'environment_override':extra or {},'started_utc':start,'elapsed_seconds':time.monotonic()-t,'exit_code':r.returncode,'summary':summary[-12:],'log_sha256':hashlib.sha256((out/(name+'.log')).read_bytes()).hexdigest()}
 receipt['steps'].append(step);save();print(name,r.returncode,step['elapsed_seconds'],step['summary'],flush=True)
 return r.returncode
save()
run('native_build',['bash','fortran_eos/build.sh','build/compose_eos_fortran'])
run('frontend_pipeline',['python3','-m','pytest','-q','python_frontend/tests','pipeline/tests'])
run('generator',['python3','-m','pytest','-q','eos_generator/tests'],{'PYTHONPATH':'eos_generator:eos_generator/compose_eos_pkg'})
run('optimizer_build',['gfortran','-std=f2018','-O2','-ffree-line-length-none','-fcheck=all','-Ibuild/compose_eos_fortran/mod','fortran_eos/tests/test_native_optimizer.f90','build/compose_eos_fortran/obj/compose_eos_types.o','build/compose_eos_fortran/obj/compose_eos_optimize_native.o','build/compose_eos_fortran/obj/slsqp_optmz.o','-o',str(out/'test_native_optimizer')])
run('optimizer',[str(out/'test_native_optimizer')])
receipt['source_sha256_after']=snapshot();receipt['source_identity_unchanged']=receipt['source_sha256_before']==receipt['source_sha256_after']
receipt['binaries_sha256']={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'build/compose_eos_fortran/compose_eos_fortran',root/'build/compose_eos_fortran/compose_eos_generate') if p.exists()}
receipt['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
receipt['status']='PASS' if receipt['source_identity_unchanged'] and all(x['exit_code']==0 for x in receipt['steps']) else 'FAIL_OR_SOURCE_CHANGED'
save()
# Keep the executable in build, not in committed evidence.
exe=out/'test_native_optimizer'
if exe.exists():
 receipt['optimizer_binary_sha256']=hashlib.sha256(exe.read_bytes()).hexdigest();save();exe.unlink()
print(receipt['status'],flush=True)
