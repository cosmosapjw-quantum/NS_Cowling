from pathlib import Path
import subprocess,os,json,time,hashlib
r=Path("/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117")
o=r/"data/crust_refinement_20261117/verification";o.mkdir(exist_ok=True)
env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
argv=["bash",str(r/"fortran_eos/build.sh"),str(r/"build/G_verified_general")]
before=subprocess.check_output(["bash",str(r/"fortran_eos/build.sh"),"--signature"],text=True).strip()
t=time.time()
with (o/"NATIVE_BUILD.log").open("w") as log:p=subprocess.run(argv,cwd=r,env=env,stdout=log,stderr=subprocess.STDOUT)
rec=dict(argv=argv,actual_exit_code=p.returncode,wall_time_seconds=time.time()-t,
source_signature_before=before,source_signature_after=subprocess.check_output(["bash",str(r/"fortran_eos/build.sh"),"--signature"],text=True).strip())
rec["binaries"]={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (r/"build/G_verified_general").glob("compose_eos_*") if p.is_file()}
(o/"NATIVE_BUILD_RECEIPT.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(rec));raise SystemExit(p.returncode)
