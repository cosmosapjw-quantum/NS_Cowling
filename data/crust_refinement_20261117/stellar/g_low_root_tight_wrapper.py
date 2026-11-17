from pathlib import Path
import os,subprocess,time,json
r=Path("/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117")
o=r/"data/crust_refinement_20261117/stellar"
argv=["/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python","/home/cosmosapjw/ns_cowling_wt/scratch/g_low_root_tight_control.py"]
env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
t=time.time()
with (o/"LOW_ROOT_TIGHT.log").open("w") as f:p=subprocess.run(argv,env=env,stdout=f,stderr=subprocess.STDOUT)
rec=dict(argv=argv,actual_exit_code=p.returncode,wall_time_seconds=time.time()-t)
(o/"LOW_ROOT_TIGHT_RECEIPT.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(rec));raise SystemExit(p.returncode)
