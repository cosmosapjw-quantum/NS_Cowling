from pathlib import Path
import subprocess,os,json,time
r=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117")
o=r/"data/low_density_crust_20261117/stellar/h_initial";o.mkdir(parents=True,exist_ok=True)
env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
COMPOSE_EOS_FORTRAN_ROOT=str(r),COMPOSE_EOS_FORTRAN_EXE="/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117/build/G_verified_general/compose_eos_fortran")
argv=["/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python",str(r/"scripts/prepare_low_density_stellar_20261117.py"),str(o/"INPUTS.json"),str(o)]
t=time.time();p=subprocess.run(argv,env=env,cwd=r,capture_output=True,text=True)
(o/"RUN.stdout.log").write_text(p.stdout);(o/"RUN.stderr.log").write_text(p.stderr)
(o/"RUN.json").write_text(json.dumps(dict(argv=argv,actual_exit_code=p.returncode,seconds=time.time()-t),indent=2)+"\n")
print(p.stdout);print(p.stderr);print("ACTUAL_EXIT",p.returncode);raise SystemExit(p.returncode)
