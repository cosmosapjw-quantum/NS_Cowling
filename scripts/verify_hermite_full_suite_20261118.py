"""Full frontend/pipeline regression and immutable source checksum receipt."""
from pathlib import Path
import os,sys,subprocess,json,hashlib,time
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/cubic_quintic_matching_20261118"
OUT.mkdir(exist_ok=True,parents=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files=sorted(set([*list((ROOT/"python_frontend").rglob("*.py")),
                  *list((ROOT/"pipeline").rglob("*.py")),
                  *list((ROOT/"scripts").glob("*20261118.py"))]))
files=[f for f in files if f.is_file()]
snapshot={str(f.relative_to(ROOT)):sha(f) for f in files}
native="/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117/build/G_verified_general"
env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
         PYTHONPATH=f"{ROOT}/scripts:{ROOT}/python_frontend:{ROOT}/pipeline:/tmp/ns_cowling_audit_20261117_pydeps",
         COMPOSE_EOS_FORTRAN_ROOT=str(ROOT),
         COMPOSE_EOS_FORTRAN_EXE=native+"/compose_eos_fortran",
         COMPOSE_EOS_GENERATOR_EXE=native+"/compose_eos_generate")
args=["/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python","-m","pytest","-q",
      "python_frontend/tests","pipeline/tests"]
t=time.monotonic()
with (OUT/"FULL_SUITE.log").open("w") as f:
    result=subprocess.run(args,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
after={str(f.relative_to(ROOT)):sha(f) for f in files}
receipt=dict(status="PASS" if result.returncode==0 and after==snapshot else "FAILED",
    actual_exit_code=result.returncode,
    elapsed_seconds=time.monotonic()-t,
    source_file_count=len(snapshot),
    source_bytes_unchanged=after==snapshot,
    source_snapshot_before=snapshot,
    fixed_native_binary_sha256={k:sha(Path(env[k])) for k in ("COMPOSE_EOS_FORTRAN_EXE","COMPOSE_EOS_GENERATOR_EXE")},
    exact_command=args,
    numerical_environment={k:env[k] for k in ("PYTHONPATH","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS")},
    log_tail=(OUT/"FULL_SUITE.log").read_text()[-2000:])
(OUT/"FULL_SUITE_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(receipt["status"],"exit",result.returncode,"source",len(snapshot),"seconds",receipt["elapsed_seconds"])
print(receipt["log_tail"][-550:])
raise SystemExit(result.returncode or (0 if after==snapshot else 4))
