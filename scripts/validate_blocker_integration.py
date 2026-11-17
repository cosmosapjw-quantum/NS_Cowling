import os,sys,json,subprocess,hashlib,time,datetime,csv
from pathlib import Path
root=Path(__file__).resolve().parents[1];os.chdir(root)
out=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else root/"build/blocker_integration_validation";out.mkdir(parents=True,exist_ok=True)
env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",PYTHONPATH="python_frontend:scripts",COMPOSE_EOS_FORTRAN_ROOT=str(root),COMPOSE_EOS_FORTRAN_EXE=str(root/"build/compose_eos_fortran/compose_eos_fortran"),COMPOSE_EOS_NATIVE_BUILD_DIR=str(root/"build/compose_eos_fortran"))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
tracked=[*root.glob("fortran_eos/src/*.f90"),*root.glob("fortran_eos/app/*.f90"),root/"fortran_eos/build.sh",*root.glob("python_frontend/compose_eos/*.py"),root/"python_frontend/tests/test_native_generation_receipts.py"]
receipt={"base_head":subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),"started":datetime.datetime.now(datetime.timezone.utc).isoformat(),"mode":sys.argv[1],"source_sha256_before":{str(p.relative_to(root)):sha(p) for p in tracked},"steps":[]}
def run(name,cmd,extra=None):
 e=env.copy();e.update(extra or {});t=time.monotonic();p=out/(name+".log")
 with p.open("w") as log:r=subprocess.run(cmd,env=e,cwd=root,stdout=log,stderr=subprocess.STDOUT)
 entry={"name":name,"command":cmd,"exit_code":r.returncode,"seconds":time.monotonic()-t,"log_sha256":sha(p),"last_lines":p.read_text().splitlines()[-3:]}
 receipt["steps"].append(entry);print(json.dumps(entry),flush=True)
 if r.returncode:raise RuntimeError(name+" failed")
try:
 run("native_build",["bash","fortran_eos/build.sh","--ensure",str(root/"build/compose_eos_fortran")])
 if sys.argv[1]=="regression":
  run("frontend_pipeline",["python3","-m","pytest","-q","python_frontend/tests","pipeline/tests"])
  run("generator",["python3","-m","pytest","-q","eos_generator/tests"],{"PYTHONPATH":"eos_generator:eos_generator/compose_eos_pkg"})
 else:
  models=["APR_APR","BBB_BHF2","GM_GM1","RG_KDE0v1","RG_SLY4"]
  for model in models:run(model,["python3","-m","compose_eos.cli","generate",model,str(out/model),"--points","17"])
  run("lower_smoke",["python3","-m","compose_eos.cli","generate","GPPVA_DD2",str(out/"lower_smoke"),"--densities",".003",".004",".005","--resolutions","16","--include-lower-crust"])
  for model in models:
   d=out/model/"generated";rows=list(csv.DictReader((d/"sound_speed.csv").open()))
   assert len(rows)==json.loads((d/"RESULT.json").read_text())["rows"]
   assert all(x["cs2_over_c2_fd"]=="" for x in rows if x["fd_status"]!="ok")
  a=json.loads((out/"lower_smoke/generated/attachment_diagnostics.json").read_text())
  assert a["outer_to_lower"]["density_overlap"] and a["outer_to_lower"]["source_prescription_valid"]
  assert a["lower_to_inner"]["density_gap"] and a["lower_to_inner"]["source_prescription_valid"]
  assert all(not a[k]["cross_interval_supported"] and not a[k]["equilibrium_established"] for k in ("outer_to_lower","lower_to_inner"))
  receipt["public_diagnostic_assertions"]="PASS"
 receipt["status"]="PASS"
except Exception as exc:receipt.update(status="FAILED",error=repr(exc))
receipt["source_sha256_after"]={str(p.relative_to(root)):sha(p) for p in tracked}
receipt["sources_unchanged"]=receipt["source_sha256_before"]==receipt["source_sha256_after"]
receipt["finished"]=datetime.datetime.now(datetime.timezone.utc).isoformat()
(out/(sys.argv[1].upper()+".json")).write_text(json.dumps(receipt,indent=2)+"\n")
print(receipt["status"],receipt["sources_unchanged"],flush=True)
sys.exit(0 if receipt["status"]=="PASS" and receipt["sources_unchanged"] else 1)
