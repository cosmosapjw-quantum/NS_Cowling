from pathlib import Path
import json,os,subprocess,time,hashlib
root=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117")
out=root/"data/low_density_crust_20261117/stellar/h_final"
out.mkdir(parents=True,exist_ok=True)
entries=[]
for model,core in [("DDME2","PK1_CANONICAL"),("DDMEX_ORIGINAL_TABLE_RHO0","GM1"),("NL3_PROVIDER_PRECISION","GM1")]:
 entries.append(dict(id="H_"+model+"_N128",crust=model,core=core,mesh=128,compose_path=f"data/low_density_crust_20261117/low_extension/repaired/tables_v2/{model}/N128/block_01/match_00_00/compose"))
entries.append(dict(id="H_GPPVA_NL3wrL55_N128",crust="GPPVA_NL3wrL55",core="PK1_CANONICAL",mesh=128,compose_path="data/low_density_crust_20261117/support_refinement/matching/GPPVA_NL3wrL55/continued_G_lower_adjacent_branch/N128/run0_component_000__component_000/compose"))
inp=out/"INPUTS.json";inp.write_text(json.dumps(entries,indent=2)+"\n")
env=os.environ.copy()
env.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",PYTHONPATH=str(root/"python_frontend")+":"+str(root/"pipeline")+":/tmp/ns_cowling_audit_20261117_pydeps",COMPOSE_EOS_FORTRAN_ROOT=str(root),COMPOSE_EOS_FORTRAN_EXE="/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117/build/G_verified_general/compose_eos_fortran")
argv=["/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python",str(root/"scripts/prepare_low_density_stellar_20261117.py"),str(inp),str(out)]
t=time.time()
with (out/"RUN.stdout.log").open("w") as a,(out/"RUN.stderr.log").open("w") as b:r=subprocess.run(argv,cwd=root,env=env,stdout=a,stderr=b)
rec=dict(argv=argv,actual_exit_code=r.returncode,seconds=time.time()-t)
(out/"RUN.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(rec))
print((out/"RUN.stdout.log").read_text())
print((out/"RUN.stderr.log").read_text())
raise SystemExit(r.returncode)
