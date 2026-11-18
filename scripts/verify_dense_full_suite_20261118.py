"""Prove native-dense interpolation new tests and existing pipeline keep passing."""
from pathlib import Path
import os,sys,json,hashlib,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/native_dense_interpolation_20261118/aggregated"
OUT.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
paths=sorted(set([*list((ROOT/"python_frontend").rglob("*.py")),
                  *list((ROOT/"pipeline").rglob("*.py")),
                  *list((ROOT/"scripts").glob("*20261118.py"))]))
paths=[p for p in paths if p.is_file()]
before={str(p.relative_to(ROOT)):sha(p) for p in paths}
native="/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117/build/G_verified_general"
env=os.environ.copy()
env.update(PYTHONPATH=f"{ROOT}/scripts:{ROOT}/python_frontend:{ROOT}/pipeline:/tmp/ns_cowling_audit_20261117_pydeps",
           OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
           COMPOSE_EOS_FORTRAN_ROOT=str(ROOT),
           COMPOSE_EOS_FORTRAN_EXE=native+"/compose_eos_fortran",
           COMPOSE_EOS_GENERATOR_EXE=native+"/compose_eos_generate")
cmd=["/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python","-m","pytest","-q",
     "python_frontend/tests","pipeline/tests"]
t=time.monotonic()
with (OUT/"FULL_SUITE.log").open("w") as f:
    result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
after={str(p.relative_to(ROOT)):sha(p) for p in paths}
receipt=dict(schema="ns_cowling.dense_interp_full_regression.v1",
             command=cmd,actual_exit_code=result.returncode,elapsed_seconds=time.monotonic()-t,
             source_snapshot_count=len(before),source_bytes_unchanged=before==after,
             source_hashes=before,
             fixed_native_exe_sha256={k:sha(Path(env[k])) for k in ("COMPOSE_EOS_FORTRAN_EXE","COMPOSE_EOS_GENERATOR_EXE")},
             environment={k:env[k] for k in ("PYTHONPATH","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS")},
             log_tail=(OUT/"FULL_SUITE.log").read_text()[-1800:])
(OUT/"FULL_SUITE_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({k:v for k,v in receipt.items() if k!="source_hashes"},indent=2))
raise SystemExit(result.returncode or (0 if before==after else 5))
