"""Run complete matching worktree frontend and pipeline regression with explicit native executable provenance."""
from pathlib import Path
import hashlib,json,os,subprocess,time
root=Path(__file__).resolve().parents[1]
out=root/"data/matching_interpolation_20261118"
out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
paths=(list((root/"python_frontend"/"compose_eos").rglob("*.py"))+
       list((root/"python_frontend"/"tests").rglob("*.py"))+
       list((root/"pipeline").rglob("*.py"))+
       list((root/"scripts").glob("*matching*20261118.py"))+
       list((root/"scripts").glob("*bernstein*20261118.py")))
paths=sorted(set(p for p in paths if p.is_file()))
before={str(p.relative_to(root)):sha(p) for p in paths}
py="/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python"
g="/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117/build/G_verified_general"
env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
         PYTHONPATH=f"{root}/scripts:{root}/python_frontend:{root}/pipeline:/tmp/ns_cowling_audit_20261117_pydeps",
         COMPOSE_EOS_FORTRAN_ROOT=str(root),
         COMPOSE_EOS_FORTRAN_EXE=str(Path(g)/"compose_eos_fortran"),
         COMPOSE_EOS_GENERATOR_EXE=str(Path(g)/"compose_eos_generate"))
args=[py,"-m","pytest","-q","python_frontend/tests","pipeline/tests"]
binary_hashes={k:sha(Path(env[k])) for k in ("COMPOSE_EOS_FORTRAN_EXE","COMPOSE_EOS_GENERATOR_EXE")}
t=time.time()
with (out/"FULL_SUITE.log").open("w") as f:
    p=subprocess.run(args,cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT)
after={str(q.relative_to(root)):sha(q) for q in paths}
r=dict(args=args,actual_exit_code=p.returncode,elapsed_seconds=time.time()-t,
       source_snapshot_count=len(before),source_bytes_unchanged=before==after,
       modified_source_files=[k for k in before if before[k]!=after.get(k)],
       native_binaries=binary_hashes,environment={k:env[k] for k in
       ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
        "PYTHONPATH","COMPOSE_EOS_FORTRAN_EXE","COMPOSE_EOS_GENERATOR_EXE")},
       log_tail=(out/"FULL_SUITE.log").read_text()[-2500:])
(out/"FULL_SUITE_RECEIPT.json").write_text(json.dumps(r,indent=2)+"\n")
print(json.dumps(r,indent=2))
raise SystemExit(p.returncode or (0 if before==after else 42))
