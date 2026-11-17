from pathlib import Path
import json,os,subprocess,time,sys,platform
root=Path('/home/cosmosapjw/ns_cowling_wt/matched_eos_reference_20261117')
base=root/'data/matched_eos_reference_20261117'
python='/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python'
env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','COMPOSE_EOS_FORTRAN_ROOT':str(root),'COMPOSE_EOS_FORTRAN_EXE':str(root/'build/compose_eos_fortran/compose_eos_fortran'),'PYTHONPATH':str(root/'python_frontend')+':/tmp/ns_cowling_audit_20261117_pydeps'}
commands=[('comparison',[python,'scripts/compare_generated_reference.py','--reference-root','data/matched_eos_reference_20261117/references','--retained-reference-root','eos_generator/eos_data','--output','data/matched_eos_reference_20261117/comparison','--plots']),('diagnosis',[python,'scripts/diagnose_reference_differences_20261117.py','--source-root','data/matched_eos_reference_20261117/references','--output','data/matched_eos_reference_20261117/diagnosis']),('reference_support',[python,'scripts/build_reference_seed_catalog_20261117.py','--references','data/matched_eos_reference_20261117/references','--out','data/matched_eos_reference_20261117/reference_support']),('frontend_pipeline',[python,'-m','pytest','-q','-ra','python_frontend/tests','pipeline/tests'])]
runs=[]
for name,command in commands:
 start=time.time()
 with (base/'validation'/f'{name}.log').open('w') as log:
  run=subprocess.run(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT)
 record={'name':name,'command':command,'returncode':run.returncode,'elapsed_seconds':time.time()-start}
 runs.append(record)
 (base/'validation/EXECUTION.json').write_text(json.dumps({'runtime':python,'runs':runs,'all_passed':all(r['returncode']==0 for r in runs)},indent=2)+'\n')
 print(json.dumps(record),flush=True)
 if run.returncode: sys.exit(run.returncode)
