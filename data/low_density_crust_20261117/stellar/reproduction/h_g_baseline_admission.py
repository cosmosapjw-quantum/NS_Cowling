from pathlib import Path
import subprocess,os,json,time,hashlib
r=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117")
out=r/"data/low_density_crust_20261117/stellar/g_baseline";out.mkdir(parents=True,exist_ok=True)
g=json.loads((r/"data/crust_refinement_20261117/G_SUMMARY.json").read_text())
inputs=[dict(id="G_"+x["crust"],crust=x["crust"],core=x["core"],mesh=x["mesh"],compose_path=x["compose_path"]) for x in g["finest_candidates"]]
(out/"INPUTS.json").write_text(json.dumps(inputs,indent=2)+"\n")
env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
 COMPOSE_EOS_FORTRAN_ROOT=str(r),COMPOSE_EOS_FORTRAN_EXE="/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117/build/G_verified_general/compose_eos_fortran")
argv=["/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python",str(r/"scripts/prepare_low_density_stellar_20261117.py"),str(out/"INPUTS.json"),str(out)]
t=time.time();p=subprocess.run(argv,cwd=r,env=env,capture_output=True,text=True)
(out/"RUN.stdout.log").write_text(p.stdout);(out/"RUN.stderr.log").write_text(p.stderr)
(out/"RUN.json").write_text(json.dumps(dict(argv=argv,actual_exit_code=p.returncode,seconds=time.time()-t,
 script_sha256=hashlib.sha256((r/"scripts/prepare_low_density_stellar_20261117.py").read_bytes()).hexdigest()),indent=2)+"\n")
print(p.stdout);print(p.stderr);print("ACTUAL_EXIT",p.returncode);raise SystemExit(p.returncode)
