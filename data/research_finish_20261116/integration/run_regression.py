from pathlib import Path
import hashlib,json,os,subprocess,time,sys,argparse
root=Path(__file__).resolve().parents[3]
os.chdir(root)
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=Path('build/replay_finish_regression'))
args=parser.parse_args();out=(root/args.output).resolve()
out.mkdir(parents=True,exist_ok=True)
if any((out/name).exists() for name in ('REGRESSION.json','frontend_pipeline.log','generator.log')):
 raise FileExistsError('use a fresh regression output; historical logs are not overwritten')
def snapshot():
 files=[]
 for folder in ('python_frontend','pipeline','eos_generator','scripts','fortran_eos'):
  files.extend(p for p in (root/folder).rglob('*') if p.is_file() and p.suffix in ('.py','.f90','.F90','.c','.h') and '__pycache__' not in p.parts and '.venv' not in p.parts)
 return {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(files))}
before=snapshot()
env=os.environ.copy()
env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONPATH='python_frontend',
 COMPOSE_EOS_FORTRAN_ROOT=str(root),COMPOSE_EOS_FORTRAN_EXE=str(root/'build/compose_eos_fortran/compose_eos_fortran'))
runs=[]
for name,command,extra in [
 ('frontend_pipeline',[sys.executable,'-m','pytest','-q','python_frontend/tests','pipeline/tests'],{}),
 ('generator',[sys.executable,'-m','pytest','-q','eos_generator/tests'],{'PYTHONPATH':'eos_generator:eos_generator/compose_eos_pkg'})]:
 start=time.monotonic()
 with (out/(name+'.log')).open('w') as f:
  run=subprocess.run(command,env={**env,**extra},stdout=f,stderr=subprocess.STDOUT)
 record={'name':name,'command':command,'returncode':run.returncode,'seconds':time.monotonic()-start,'log':name+'.log'}
 runs.append(record);print(json.dumps(record),flush=True)
after=snapshot()
receipt={'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
 'runs':runs,'tested_source_sha256_before':before,'tested_source_sha256_after':after,
 'tested_sources_unchanged':before==after,'status':'PASS' if before==after and all(r['returncode']==0 for r in runs) else 'FAIL',
 'binary_sha256':hashlib.sha256((root/'build/compose_eos_fortran/compose_eos_fortran').read_bytes()).hexdigest(),
 'runtime':{'python':sys.version},'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(out/'REGRESSION.json').write_text(json.dumps(receipt,indent=2)+'\n')
for name in ['frontend_pipeline','generator']:
 print(name,(out/(name+'.log')).read_text()[-1800:])
print(receipt['status'],flush=True)
raise SystemExit(0 if receipt['status']=='PASS' else 1)
