"""Full NS_Cowling frontend/pipeline regression with byte-preserved sources.

Does not relax scientific acceptance gates to turn missing external crust
thermodynamics into matching success.
"""
from pathlib import Path
import os,sys,subprocess,json,hashlib,time
R=Path(__file__).resolve().parents[1]
D=R/"data/all_origin_maxwell_css_20261118"
D.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
files=sorted(set([*list((R/"python_frontend").rglob("*.py")),
                  *list((R/"pipeline").rglob("*.py")),
                  *list((R/"scripts").glob("*20261118.py"))]))
files=[f for f in files if f.is_file()]
before={str(f.relative_to(R)):sha(f) for f in files}
py="/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python"
g="/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117/build/G_verified_general"
env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
         PYTHONPATH=f"{R}/scripts:{R}/python_frontend:{R}/pipeline:/tmp/ns_cowling_audit_20261117_pydeps",
         COMPOSE_EOS_FORTRAN_ROOT=str(R),
         COMPOSE_EOS_FORTRAN_EXE=g+"/compose_eos_fortran",
         COMPOSE_EOS_GENERATOR_EXE=g+"/compose_eos_generate")
cmd=[py,"-m","pytest","-q","python_frontend/tests","pipeline/tests"]
t=time.monotonic()
with (D/"FULL_SUITE.log").open("w") as f:
    result=subprocess.run(cmd,cwd=R,env=env,stdout=f,stderr=subprocess.STDOUT)
after={str(f.relative_to(R)):sha(f) for f in files}
receipt=dict(status="PASS" if result.returncode==0 and before==after else "FAILED",
       actual_exit_code=result.returncode,elapsed_seconds=time.monotonic()-t,
       source_snapshot_count=len(before),source_bytes_unchanged=before==after,
       source_drift=[k for k in before if before[k]!=after.get(k)],
       fixed_native_binaries={k:sha(Path(env[k])) for k in ("COMPOSE_EOS_FORTRAN_EXE","COMPOSE_EOS_GENERATOR_EXE")},
       exact_command=cmd,environment={k:env[k] for k in ("PYTHONPATH","OMP_NUM_THREADS","OPENBLAS_NUM_THREADS")},
       log_tail=(D/"FULL_SUITE.log").read_text()[-2000:])
(D/"FULL_SUITE_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
print("FULL_SUITE",receipt["status"],result.returncode,"checks",len(files),"seconds",receipt["elapsed_seconds"])
print(receipt["log_tail"][-620:])
raise SystemExit(result.returncode or (0 if before==after else 3))
