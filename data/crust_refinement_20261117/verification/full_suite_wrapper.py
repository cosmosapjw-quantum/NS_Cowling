from pathlib import Path
import os,subprocess,time,json,hashlib,sys
r=Path("/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117")
o=r/"data/crust_refinement_20261117/verification";o.mkdir(exist_ok=True)
def snapshot():
 paths=[]
 for d in ["python_frontend","pipeline","fortran_eos/src","fortran_eos/app","fortran_eos/tools","scripts"]:
  paths += [p for p in (r/d).rglob("*") if p.is_file() and p.suffix in [".py",".f90",".sh"]]
 return {str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
before=snapshot()
(o/"CODE_BEFORE_FULL_SUITE.json").write_text(json.dumps(before,indent=2)+"\n")
env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
 PYTHONPATH=str(r/"python_frontend")+":"+str(r/"pipeline")+":/tmp/ns_cowling_audit_20261117_pydeps",
 COMPOSE_EOS_FORTRAN_ROOT=str(r),COMPOSE_EOS_FORTRAN_EXE=str(r/"build/G_verified_general/compose_eos_fortran"),
 COMPOSE_EOS_GENERATOR_EXE=str(r/"build/G_verified_general/compose_eos_generate"))
argv=["/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python","-m","pytest","-q","python_frontend/tests","pipeline/tests"]
t=time.time()
with (o/"FULL_SUITE.log").open("w") as f:p=subprocess.run(argv,cwd=r,env=env,stdout=f,stderr=subprocess.STDOUT)
after=snapshot()
(o/"CODE_AFTER_FULL_SUITE.json").write_text(json.dumps(after,indent=2)+"\n")
diff=subprocess.run(["git","diff","--check"],cwd=r,capture_output=True,text=True)
rec=dict(argv=argv,actual_exit_code=p.returncode,wall_time_seconds=time.time()-t,
 source_count=len(before),source_bytes_unchanged=before==after,
 changed_paths=[k for k in sorted(set(before)|set(after)) if before.get(k)!=after.get(k)],
 git_diff_check_exit=diff.returncode,git_diff_check_stdout=diff.stdout,git_diff_check_stderr=diff.stderr,
 environment={k:env[k] for k in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","PYTHONPATH","COMPOSE_EOS_FORTRAN_ROOT","COMPOSE_EOS_FORTRAN_EXE","COMPOSE_EOS_GENERATOR_EXE"]})
(o/"FULL_SUITE_RECEIPT.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(rec));print((o/"FULL_SUITE.log").read_text()[-2000:])
raise SystemExit(p.returncode or diff.returncode or (0 if before==after else 9))
