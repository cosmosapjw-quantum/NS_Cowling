from pathlib import Path
import subprocess,os,json,time
r=Path("/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117")
o=r/"data/crust_refinement_20261117/stellar"
env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
argv=["/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python",str(r/"scripts/run_sharp_stellar_20261117.py")]
t=time.time()
with (o/"PRODUCTION.log").open("w") as log:
 p=subprocess.run(argv,cwd=r,env=env,stdout=log,stderr=subprocess.STDOUT)
rec=dict(argv=argv,actual_exit_code=p.returncode,wall_time_seconds=time.time()-t)
(o/"PRODUCTION_RECEIPT.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(rec))
raise SystemExit(p.returncode)
