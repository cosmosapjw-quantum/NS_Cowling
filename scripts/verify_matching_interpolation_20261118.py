"""Reproducible selected legacy plus exploratory matching regression receipt."""
from pathlib import Path
import subprocess,os,time,hashlib,json,sys
R=Path(__file__).resolve().parents[1];O=R/"data/matching_interpolation_20261118";O.mkdir(exist_ok=True)
tests=["test_matching_interpolation_20261118.py","test_causal_lower_bridge.py",
       "test_whole_crust.py","test_stellar_sharp.py","test_stellar_background.py"]
sources=[R/"python_frontend/tests"/n for n in tests]
sources += [R/"scripts"/n for n in
            ["matching_interpolation_shapes_20261118.py",
             "run_matching_interpolation_crosscheck_20261118.py",
             "run_slope_matched_bernstein_20261118.py",
             "compare_transition_adjacent_interpolation_20261118.py",
             "verify_matching_interpolation_20261118.py"]]
def snap():return {str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
before=snap()
env=os.environ.copy()
env.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
           PYTHONPATH=f"{R}/scripts:{R}/python_frontend:{R}/pipeline:/tmp/ns_cowling_audit_20261117_pydeps",
           COMPOSE_EOS_FORTRAN_ROOT=str(R),
           COMPOSE_EOS_FORTRAN_EXE="/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117/build/G_verified_general/compose_eos_fortran",
           COMPOSE_EOS_GENERATOR_EXE="/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117/build/G_verified_general/compose_eos_generate")
cmd=["/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python","-m","pytest","-q"]+[str(p) for p in sources[:5]]
t=time.time()
with (O/"REGRESSION.log").open("w") as out:
    result=subprocess.run(cmd,cwd=R,env=env,stdout=out,stderr=subprocess.STDOUT)
after=snap()
receipt=dict(command=cmd,actual_exit_code=result.returncode,elapsed_seconds=time.time()-t,
             source_snapshot_count=len(before),source_bytes_unchanged=before==after,
             source_hashes=before,stdout_tail=(O/"REGRESSION.log").read_text()[-1400:])
(O/"REGRESSION_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({k:v for k,v in receipt.items() if k!="source_hashes"},indent=2))
raise SystemExit(result.returncode or (0 if before==after else 5))
