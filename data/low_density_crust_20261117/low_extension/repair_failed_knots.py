#!/usr/bin/env python3
"""One bounded opposite-direction continuation repair for observed optimizer failures."""
from pathlib import Path
from datetime import datetime, timezone
import argparse,json,os,sys,time,traceback
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend"),str(ROOT/"pipeline")]
import extend_low_density_crust_20261117 as production
from compose_eos.tf_candidate import sha,save,run_state,interpolate_initial
import run_crust_refinement_20261117 as reference
BASE=production.BASE/"repaired"
def load(p):return json.loads(Path(p).read_text())
def stamp():return datetime.now(timezone.utc).isoformat()
def run(model,execute):
 source=production.BASE/f"{model}_N128_STATES.json";rows=load(source)
 original_points=load(production.BASE/f"{model}_POINTS.json")
 tasks=[]
 for i,s in enumerate(rows):
  if s["accepted_stationary"]:continue
  first=load(original_points[i]["attempts"][0]["path"])
  parent=first["parent"]
  parent_path=parent["state_path"] if isinstance(parent,dict) else parent
  parent_n=load(parent_path)["n_b"]
  opposite=[t for t in rows if t["accepted_stationary"] and (t["n_b"]-s["n_b"])*(parent_n-s["n_b"])<0]
  if not opposite:
   tasks.append({"index":i,"n_b":s["n_b"],"eligible":False,"reason":"no accepted opposite-direction source"});continue
  other=min(opposite,key=lambda t:abs(t["n_b"]-s["n_b"]))
  eligible=(s.get("optimizer_mode") in (8.,9.) and s.get("kkt",1)<.02 and s.get("radius_stationarity",1)<2e-4 and s.get("radius_bound",1)==0 and s.get("density_bound",1)==0)
  tasks.append({"index":i,"n_b":s["n_b"],"eligible":eligible,"reason":"observed mode 8/9 with small physical residuals; change continuation direction using an independently accepted profile",
    "original_failure":{"path":s["state_path"],"sha256":sha(s["state_path"]),"mode":s["optimizer_mode"],"kkt":s["kkt"],"radius_stationarity":s["radius_stationarity"]},
    "original_parent_n":parent_n,"repair_parent":{"path":other["state_path"],"sha256":sha(other["state_path"]),"n_b":other["n_b"]}})
 plan={"schema":"ns_cowling.H_opposite_direction_repair_plan.v1","model":model,"production_plan_sha256":sha(production.BASE/"PLAN.json"),
  "source_index_sha256":sha(source),"script_sha256":sha(__file__),"driver_sha256":sha(production.EXE),
  "mesh_points":128,"ftol":1e-14,"maxiter":1800,"maximum_extra_runs_per_failed_state":1,
  "policy":"Distinct starting profile and continuation direction; no additional returned-profile restart, no tolerance/gate relaxation, no spatial-mesh mixing.",
  "tasks":tasks}
 dest=BASE/model;dest.mkdir(parents=True,exist_ok=True)
 pp=dest/"PLAN.json"
 if pp.exists():assert load(pp)==plan
 else:save(pp,plan)
 print(json.dumps({"event":"REPAIR_PLAN","model":model,"tasks":len(tasks),"plan_sha256":sha(pp)}),flush=True)
 if not execute:return
 receipt={"status":"RUNNING","started_at_utc":stamp(),"pid":os.getpid(),"actual_exit_code":None,"plan_sha256":sha(pp)}
 save(dest/"PROCESS.json",receipt)
 start=time.monotonic();results=[]
 try:
  for task in tasks:
   if not task["eligible"]:
    results.append(dict(task,status="NOT_ATTEMPTED"));continue
   i=task["index"];parent=load(task["repair_parent"]["path"])
   work=dest/"states"/f"n{task['n_b']:.9f}";p=work/"STATE.json"
   if p.exists():
    s=load(p);assert s["status"]!="RUNNING" and s["executable_sha256"]==production.EXE_SHA and s["points"]==128
   else:
    nn,np_,r=interpolate_initial(parent,128,task["n_b"])
    s=run_state(production.EXE,model,task["n_b"],128,3,"droplet",r,work_dir=work,
      initial_nn=nn,initial_np=np_,parent={"operation":"one opposite-direction accepted-profile continuation repair","state_path":parent["state_path"],"state_sha256":sha(parent["state_path"])},
      maxiter=1800,ftol=1e-14,timeout=1200)
   rec=dict(task,status=s["status"],accepted=s["accepted_stationary"],actual_exit_code=s.get("actual_exit_code"),
      repaired_state_path=s["state_path"],repaired_state_sha256=sha(s["state_path"]),wall_seconds=s.get("wall_seconds"))
   results.append(rec)
   if s["accepted_stationary"]:rows[i]=s
   print(json.dumps({"event":"REPAIR_STATE","model":model,"n_b":s["n_b"],"status":s["status"],"kkt":s.get("kkt"),"optimizer_mode":s.get("optimizer_mode")}),flush=True)
  save(dest/"REPAIRS.json",results);save(dest/"FINAL_N128_STATES.json",rows)
  reference.BASE=BASE;reference.CONFIG=production.CONFIG;reference.EXE=production.EXE
  match=reference.export_and_match(model,128,rows)
  summary={"model":model,"source":str(source),"source_sha256":sha(source),"final_source":str(dest/"FINAL_N128_STATES.json"),
    "final_source_sha256":sha(dest/"FINAL_N128_STATES.json"),"requested":len(rows),"accepted":sum(s["accepted_stationary"] for s in rows),
    "repairs":results,"failed_densities":[s["n_b"] for s in rows if not s["accepted_stationary"]],
    "match_summary":str(BASE/"tables_v2"/model/"N128"/"SUMMARY.json")}
  save(dest/"SUMMARY.json",summary)
  receipt.update(status="COMPLETED",actual_exit_code=0,finished_at_utc=stamp(),wall_seconds=time.monotonic()-start,
    accepted=summary["accepted"],requested=len(rows),native_invocations=sum(t["eligible"] for t in tasks))
  save(dest/"PROCESS.json",receipt)
 except BaseException as exc:
  receipt.update(status="FAILED",actual_exit_code=1,finished_at_utc=stamp(),wall_seconds=time.monotonic()-start,
    exception_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
  save(dest/"PROCESS.json",receipt);raise
def main():
 a=argparse.ArgumentParser();a.add_argument("phase",choices=("prepare","run"));a.add_argument("--models",nargs="+",required=True,choices=list(production.CONFIG));args=a.parse_args()
 os.environ["COMPOSE_EOS_FORTRAN_ROOT"]=str(ROOT);os.environ["COMPOSE_EOS_FORTRAN_EXE"]=str(production.MATCH)
 for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):os.environ[k]="1"
 assert sha(production.EXE)==production.EXE_SHA and sha(production.MATCH)==production.MATCH_SHA
 for model in args.models:run(model,args.phase=="run")
if __name__=="__main__":main()
