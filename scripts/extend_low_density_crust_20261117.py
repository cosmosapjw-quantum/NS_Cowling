#!/usr/bin/env python3
"""H: bounded, actual N128 TF continuation toward the frozen outer-crust reference."""
from __future__ import annotations
import argparse, json, os, sys, time, traceback
from pathlib import Path
from datetime import datetime, timezone
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/"python_frontend"), str(ROOT/"scripts"), str(ROOT/"pipeline")]
from compose_eos.tf_candidate import run_state, interpolate_initial, sha, save
import run_crust_refinement_20261117 as reference

BASE = ROOT/"data/low_density_crust_20261117/low_extension"
GROOT = ROOT.parent/"crust_refinement_20261117"
EXE = GROOT/"build/tf_candidate_G_final/tf_candidate_driver"
MATCH = GROOT/"build/G_verified_general/compose_eos_fortran"
EXE_SHA = "21f00ecd9cd492e37329e5690d342842c19025d4b649f97f99e84e592522c189"
MATCH_SHA = "533d55ba2b5264e3e5248395baefc9463b4732b818e3e8af14beb774779f5356"
BASE_COMMIT = "905b845ad0f14d2f1e7f67f9eb802691144850cd"
POINTS = 128
LOW = [.0003,.0004,.0005,.00075,.001,.0015,.002,.003,.004,.005]
CONFIG = {
 "DDME2": {"densities": LOW + np.linspace(.006,.014,9).tolist(),
           "anchor": .006, "mn":939., "mp":939., "core":"new_native:PK1_CANONICAL:fine",
           "reuse":"data/crust_refinement_20261117/tf_regeneration/DDME2_N128_STATES.json"},
 "DDMEX_ORIGINAL_TABLE_RHO0": {"densities": LOW[:-1] + np.linspace(.005,.012,29).tolist(),
           "anchor":.008, "mn":938.5, "mp":938.5, "core":"native_generated:GM_GM1"},
 "NL3_PROVIDER_PRECISION": {"densities": LOW + [.006,.008,.010,.012,.014,.016,.018,.020,.022,.024] + np.linspace(.026,.038,25).tolist(),
           "anchor":.030, "mn":939., "mp":939., "core":"native_generated:GM_GM1"},
}
G_INDEX = ROOT/"data/crust_refinement_20261117/support_repair/REGENERATION.json"

def load(p):
 return json.loads(Path(p).read_text())

def source_state(path):
 p=Path(path)
 # Existing G source identity is frozen; use its exact bytes at the original path.
 assert p.is_file(), p
 s=load(p)
 assert s["state_path"]==str(p) and s["accepted_stationary"]
 return s

def stamp():
 return datetime.now(timezone.utc).isoformat()

def source_records(model):
 if "reuse" in CONFIG[model]:
  rows=load(ROOT/CONFIG[model]["reuse"])
  return [source_state(s["state_path"]) for s in rows]
 run=next(r for r in load(G_INDEX)["runs"] if r["model"]==model and r["points"]==64 and r["anchor"]=="provider_radius_seed")
 return [source_state(r["state_path"]) for r in run["points_results"] if r["accepted"]]

def plan():
 assert sha(EXE)==EXE_SHA and sha(MATCH)==MATCH_SHA
 sources={}
 for model in CONFIG:
  states=source_records(model)
  sources[model]=[{"n_b":s["n_b"],"points":s["points"],"path":s["state_path"],"sha256":sha(s["state_path"])} for s in states]
 files=[Path(__file__),ROOT/"python_frontend/compose_eos/tf_candidate.py",
        ROOT/"fortran_eos/src/compose_eos_tf_models.f90",ROOT/"fortran_eos/tools/tf_candidate_driver.f90",
        ROOT/"scripts/run_crust_refinement_20261117.py",G_INDEX,
        ROOT/CONFIG["DDME2"]["reuse"],ROOT/"data/new_model_maxwell_20261117/prepared/MANIFEST.json"]
 p={"schema":"ns_cowling.H_low_extension_plan.v1","created_at_utc":stamp(),"base_commit":BASE_COMMIT,
    "configurations":CONFIG,"mesh_points":POINTS,"source_states":sources,
    "source_sha256":{str(f.relative_to(ROOT)):sha(f) for f in files},
    "executables":{str(EXE):EXE_SHA,str(MATCH):MATCH_SHA},
    "units":{"n_b":"fm^-3","pressure":"MeV fm^-3","epsilon":"MeV fm^-3","chemical_potential":"MeV","radius":"fm"},
    "stationarity":{"ftol":1e-14,"maxiter":1800,"native_acceptance":"unchanged","warm_restarts_maximum":1},
    "scope":"Fixed-energy, fixed-parameter, spherical droplet continuation on N128; finite branch evidence only.",
    "policy":"Keep all failed requested nodes; export only contiguous accepted same-mesh runs. No energy shifts, parameter fitting, failed-gap filling, or global-minimum claim.",
    "warm_restart_condition":"One returned-profile restart only for mode 8/9, KKT < .02, radius stationarity < 2e-4, and inactive bounds.",
    "outer_reference":{"n_b":.000269021,"pressure":.0004981910045986,"epsilon":.2520862865914808,"mu_B":938.9024559275276},
    "future_repair":"Only after observed failure; targeted spacing/radius seed or N256 with separate checkpoint and preserved N128 evidence."}
 existing=BASE/"PLAN.json"
 if existing.exists():
  previous=load(existing)
  for k in ("configurations","mesh_points","source_states","source_sha256","executables","stationarity"):
   assert previous[k]==p[k], ("PLAN_MISMATCH",k)
  return previous
 save(existing,p)
 return p

def attempt(parent,model,n,dest,operation):
 p=dest/"STATE.json"
 if p.exists():
  s=load(p)
  assert s["model"]==model and s["points"]==POINTS and abs(s["n_b"]-n)<1e-15
  assert s["executable_sha256"]==EXE_SHA
  assert s["status"]!="RUNNING", "Inspect interrupted native state before resuming"
  return s
 nn,pp,r=interpolate_initial(parent,POINTS,n)
 provenance={"operation":operation,"state_path":parent["state_path"],"state_sha256":sha(parent["state_path"])}
 return run_state(EXE,model,n,POINTS,3,"droplet",r,work_dir=dest,initial_nn=nn,initial_np=pp,
                  parent=provenance,maxiter=1800,ftol=1e-14,timeout=1200)

def solve(parent,model,n,index):
 dest=BASE/"states"/model/f"N128_i{index:03d}"
 first=attempt(parent,model,n,dest/"primary","N128 spatial transfer / adjacent-density continuation; initial guess only")
 attempts=[first]
 s=first
 if (not first.get("accepted_stationary") and "profile" in first and
     first.get("optimizer_mode") in (8.,9.) and first.get("kkt",1)<.02 and
     first.get("radius_stationarity",1)<2e-4 and first.get("radius_bound",1)==0 and
     first.get("density_bound",1)==0):
  s=attempt(first,model,n,dest/"warm_restart","one returned-profile optimizer restart; unchanged physical gates and tolerance")
  attempts.append(s)
 rec={"model":model,"n_b":n,"points":POINTS,"accepted":bool(s.get("accepted_stationary")),
      "state_path":s["state_path"],"state_sha256":sha(s["state_path"]),
      "attempts":[{"path":a["state_path"],"sha256":sha(a["state_path"]),"status":a["status"],
                   "accepted_stationary":a["accepted_stationary"],"actual_exit_code":a.get("actual_exit_code"),
                   "wall_seconds":a.get("wall_seconds"),"optimizer_mode":a.get("optimizer_mode"),
                   "kkt":a.get("kkt"),"radius_stationarity":a.get("radius_stationarity"),
                   "radius_bound":a.get("radius_bound"),"density_bound":a.get("density_bound")} for a in attempts]}
 save(dest/"POINT.json",rec)
 print(json.dumps({"event":"STATE",**{k:s.get(k) for k in ("model","n_b","points","status","pressure","mu_B","radius","wall_seconds")}}),flush=True)
 return s,rec

def build_model(model):
 cfg=CONFIG[model];dens=cfg["densities"];source=source_records(model)
 rows=[None]*len(dens);records=[None]*len(dens)
 if "reuse" in cfg:
  for old in source:
   i=min(range(len(dens)),key=lambda j:abs(dens[j]-old["n_b"]))
   assert abs(dens[i]-old["n_b"])<1e-15 and old["points"]==POINTS
   rows[i]=old;records[i]={"model":model,"n_b":old["n_b"],"points":POINTS,"accepted":True,"reused_exact_G":True,
      "state_path":old["state_path"],"state_sha256":sha(old["state_path"]),"attempts":[]}
  anchor=min(i for i,s in enumerate(rows) if s is not None)
 else:
  anchor=min(range(len(dens)),key=lambda j:abs(dens[j]-cfg["anchor"]))
  old=min(source,key=lambda s:abs(s["n_b"]-dens[anchor]))
  rows[anchor],records[anchor]=solve(old,model,dens[anchor],anchor)
 def checkpoint():
  save(BASE/f"{model}_PROGRESS.json",{"updated_at_utc":stamp(),"requested":len(dens),
    "finished":sum(s is not None for s in rows),"accepted":sum(bool(s and s.get("accepted_stationary")) for s in rows),"records":records})
  # Completed, ordered snapshots can be inspected while later source knots run.
  save(BASE/f"{model}_COMPLETED_STATES.json",[s for s in rows if s is not None])
 checkpoint()
 # Downward production is first: make lower-boundary evidence available promptly.
 for order in (range(anchor-1,-1,-1),range(anchor+1,len(dens))):
  last=rows[anchor]
  if not last.get("accepted_stationary"):
   last=min(source,key=lambda s:abs(s["n_b"]-dens[anchor]))
  for i in order:
   if rows[i] is None:
    rows[i],records[i]=solve(last,model,dens[i],i)
    checkpoint()
   if rows[i].get("accepted_stationary"):last=rows[i]
 assert all(s is not None for s in rows)
 save(BASE/f"{model}_N128_STATES.json",rows)
 save(BASE/f"{model}_POINTS.json",records)
 reference.BASE=BASE;reference.CONFIG=CONFIG;reference.EXE=EXE
 result=reference.export_and_match(model,POINTS,rows)
 summary={"model":model,"requested":len(rows),"accepted":sum(s["accepted_stationary"] for s in rows),
  "failures":[{"n_b":s["n_b"],"state_path":s["state_path"],"status":s["status"]} for s in rows if not s["accepted_stationary"]],
  "accepted_density_bounds":[min(s["n_b"] for s in rows if s["accepted_stationary"]),max(s["n_b"] for s in rows if s["accepted_stationary"])],
  "actual_native_invocations":sum(len(r["attempts"]) for r in records),
  "source_index_sha256":sha(BASE/f"{model}_N128_STATES.json"),
  "match_summary":str(BASE/"tables_v2"/model/"N128"/"SUMMARY.json"),
  "matches":[{"status":m["status"],"compose_path":m.get("compose_path"),"error":m.get("error"),
    "native_report":m.get("native_report")} for b in result["blocks"] for m in b["matches"]]}
 save(BASE/f"{model}_SUMMARY.json",summary)
 print(json.dumps({"event":"MODEL_DONE","model":model,"accepted":summary["accepted"],"requested":len(rows),"native_calls":summary["actual_native_invocations"],"match_statuses":[m["status"] for m in summary["matches"]]}),flush=True)
 return summary

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("phase",choices=("prepare","run"))
 ap.add_argument("--models",nargs="+",choices=list(CONFIG),default=list(CONFIG))
 args=ap.parse_args()
 BASE.mkdir(parents=True,exist_ok=True)
 os.environ["COMPOSE_EOS_FORTRAN_ROOT"]=str(ROOT)
 os.environ["COMPOSE_EOS_FORTRAN_EXE"]=str(MATCH)
 for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):os.environ[k]="1"
 p=plan()
 print(json.dumps({"event":"PLAN_VERIFIED","path":str(BASE/"PLAN.json"),"sha256":sha(BASE/"PLAN.json"),"new_grid_points":{m:len(c["densities"]) for m,c in CONFIG.items()}}),flush=True)
 if args.phase=="prepare":return
 started=time.monotonic()
 receipt={"started_at_utc":stamp(),"models":args.models,"plan_sha256":sha(BASE/"PLAN.json"),
    "pid":os.getpid(),"status":"RUNNING","actual_exit_code":None}
 receipt_path=BASE/("PROCESS_"+"_".join(args.models)+".json")
 save(receipt_path,receipt)
 try:
  results=[build_model(m) for m in args.models]
  receipt.update(status="COMPLETED",actual_exit_code=0,finished_at_utc=stamp(),wall_seconds=time.monotonic()-started,
    results=[str(BASE/f"{m}_SUMMARY.json") for m in args.models],source_bytes_unchanged=all(sha(ROOT/f)==v for f,v in p["source_sha256"].items()))
  assert receipt["source_bytes_unchanged"], "Source code/data identity changed during native production"
  save(receipt_path,receipt)
 except BaseException as exc:
  receipt.update(status="FAILED",actual_exit_code=1,finished_at_utc=stamp(),wall_seconds=time.monotonic()-started,
    exception_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
  save(receipt_path,receipt)
  raise
if __name__=="__main__":main()
