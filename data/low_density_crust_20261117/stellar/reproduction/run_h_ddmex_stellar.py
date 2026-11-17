from pathlib import Path
import os,subprocess,time,json,traceback
root=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117")
out=root/"data/low_density_crust_20261117/stellar/ddmex_comparison"
env=os.environ.copy();env.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
argv=["/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python",str(root/"scripts/run_low_density_stellar_20261117.py"),str(out/"PLAN.json"),str(out)]
t=time.time()
with (out/"RUN.stdout.log").open("w") as a,(out/"RUN.stderr.log").open("w") as b:
 r=subprocess.run(argv,cwd=root,env=env,stdout=a,stderr=b)
receipt=dict(argv=argv,actual_exit_code=r.returncode,seconds=time.time()-t)
(out/"RUN.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(receipt))
print((out/"RUN.stdout.log").read_text()[-6000:])
print((out/"RUN.stderr.log").read_text()[-6000:])
raise SystemExit(r.returncode)
