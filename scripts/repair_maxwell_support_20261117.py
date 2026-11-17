#!/usr/bin/env python3
"""G targeted support diagnosis and actual source regeneration orchestration."""
from __future__ import annotations
import argparse, json, sys, hashlib
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from numpy.polynomial import Polynomial
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from maxwell_three_family_campaign import read_prepared
from compose_eos.maxwell_v2 import _segment_eosdata
BASE=ROOT/"data/crust_refinement_20261117/support_repair"
FBASE=ROOT/"data/new_model_maxwell_20261117"
TARGETS=[
 ("old_native_crust__new_native_core","GPPVA_NL3wrL55__PK1_CANONICAL",1),
 ("old_native_crust__new_native_core","GPPVA_TM1e__PK1_CANONICAL",1),
 ("reference_assisted_hybrid","XMLSLZ_NL3__GM_GM1",8),
 ("reference_assisted_hybrid","XMLSLZ_DDME-X__GM_GM1",4),
]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def save(p,d):
 Path(p).parent.mkdir(parents=True,exist_ok=True)
 Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2)+"\n")
def poly_on_unit(spline,a,b):
 i=min(max(np.searchsorted(spline.x,a,side="right")-1,0),len(spline.x)-2)
 c=spline.c[:,i];x=Polynomial([a-spline.x[i],b-a])
 return ((c[0]*x+c[1])*x+c[2])*x+c[3]
def diagnose():
 records={r["id"]:r for r in load(FBASE/"prepared/MANIFEST.json")["source_records"]}
 campaign=load(FBASE/"campaign/CAMPAIGN.json")
 cells={(c["category"],c["id"]):c for c in campaign["cells"]}
 results=[]
 for category,name,index in TARGETS:
  c=cells[category,name];original=c["candidate_results"][index]
  arrays=[]
  for key,skey in [("crust_component","crust_source_id"),("core_component","core_source_id")]:
   comp=original[key];data=read_prepared(records[c[skey]])
   arrays.append(_segment_eosdata(data,(comp["phase_id"],comp["start_zero_based"],comp["stop_exclusive"])))
  low,high=arrays
  ml=(low.eps+low.p)/low.nb;mh=(high.eps+high.p)/high.nb
  sl=CubicHermiteSpline(ml,low.p,low.nb);sh=CubicHermiteSpline(mh,high.p,high.nb)
  lo,hi=max(ml[0],mh[0]),min(ml[-1],mh[-1])
  union=np.unique(np.r_[lo,ml[(ml>lo)&(ml<hi)],mh[(mh>lo)&(mh<hi)],hi]);roots=[]
  for a,b in zip(union[:-1],union[1:]):
   for z in (poly_on_unit(sh,a,b)-poly_on_unit(sl,a,b)).roots():
    if abs(z.imag)<=1e-9 and -1e-10<=z.real<=1+1e-10:
     mu=float(a+(b-a)*z.real)
     if not roots or abs(mu-roots[-1])>1e-8:roots.append(mu)
  root_records=[]
  for mu in roots:
   tol=1e-10*max(1,abs(mu));nlow=int(np.sum(ml<mu-tol));nhigh=int(np.sum(mh>mu+tol))
   root_records.append({"mu_transition_MeV":mu,"pressure_MeV_fm3":float(sl(mu)),
    "n_low_fm3":float(sl(mu,1)),"n_high_fm3":float(sh(mu,1)),
    "retained_low_source_knots":nlow,"retained_high_source_knots":nhigh,
    "native_mu_exclusion_tolerance":tol,
    "low_source_bracket_nb":low.nb[max(0,np.searchsorted(ml,mu)-1):min(low.n,np.searchsorted(ml,mu)+1)].tolist(),
    "low_source_bracket_mu":ml[max(0,np.searchsorted(ml,mu)-1):min(low.n,np.searchsorted(ml,mu)+1)].tolist()})
  results.append({"category":category,"cell":name,"candidate_index":index,
   "original_status":original["status"],"original_error":original["error"],
   "crust_source_id":c["crust_source_id"],"core_source_id":c["core_source_id"],
   "crust_component":original["crust_component"],"core_component":original["core_component"],
   "roots":root_records,"source_identity":{k:{"internal":records[c[k]]["internal_sha256"],"raw":records[c[k]]["source_table_sha256"]} for k in ["crust_source_id","core_source_id"]},
   "interpretation":"Independent root diagnostics only. No output admitted or source row replaced."})
 d={"schema":"ns_cowling.G_concrete_support_diagnosis.v1","F_campaign_sha256":sha(FBASE/"campaign/CAMPAIGN.json"),
  "F_prepared_manifest_sha256":sha(FBASE/"prepared/MANIFEST.json"),"script_sha256":sha(__file__),"results":results}
 save(BASE/"DIAGNOSIS.json",d)
 print(json.dumps({"diagnosed":len(results),"roots":[{"cell":r["cell"],"roots":r["roots"]} for r in results]},indent=2))

from collections import Counter
from datetime import datetime, timezone
import csv, os, subprocess, time
from compose_eos.tf_candidate import build_driver,run_state,interpolate_initial,export_compose
from compose_eos.api import read_compose,write_compose
from compose_eos.native import executable,run_native
from maxwell_three_family_campaign import plain,rel,observe_kernel
from run_new_model_maxwell_20261117 import prepare_record
from match_report_admission import admit_match_report
from match_core_crust import verify_internal_readback,verify_internal_metadata_readback

CONFIGS=[
 {"id":"GPPVA_NL3wrL55","model":"GPPVA_NL3wrL55","core_id":"new_native:PK1_CANONICAL:fine",
  "bounds":[.010,.015],"step":.00025,"extra_density":[.01201,.01202,.01204,.01208],"mn":939.,"mp":939.,
  "anchors":[{"name":"upper_F_component","kind":"native","node":11},
             {"name":"lower_adjacent_branch","kind":"native","node":10}]},
 {"id":"GPPVA_TM1e","model":"GPPVA_TM1e","core_id":"new_native:PK1_CANONICAL:fine",
  "bounds":[.006,.012],"step":.00025,"mn":938.,"mp":938.,
  "anchors":[{"name":"upper_F_component","kind":"native","node":7},
             {"name":"lower_adjacent_branch","kind":"native","node":6}]},
 {"id":"NL3_PROVIDER_PRECISION","model":"NL3_PROVIDER_PRECISION","core_id":"native_generated:GM_GM1",
  "bounds":[.026,.038],"step":.0005,"mn":939.,"mp":939.,
  "anchors":[{"name":"provider_radius_seed","kind":"public","source_model":"XMLSLZ_NL3","density":.03}]},
 {"id":"DDMEX_ORIGINAL_TABLE_RHO0","model":"DDMEX_ORIGINAL_TABLE_RHO0","core_id":"native_generated:GM_GM1",
  "bounds":[.005,.012],"step":.00025,"mn":938.5,"mp":938.5,
  "anchors":[{"name":"provider_radius_seed","kind":"public","source_model":"XMLSLZ_DDME-X","density":.008}]},
]
MESHES=[16,32,64]
MATCH_EXE_SHA="3893fbe061c2564ed4138211f7ea379a973851eb6e785fad22edcec0875152b1"

def seed_for(config,anchor):
 if anchor["kind"]=="native":
  src=ROOT/"data/research_native_20261116/source_grid_cohort"/config["model"]/"attempt_001/output/generated"
  sel=next(x for x in csv.DictReader((src/"selected.csv").open()) if int(x["node"])==anchor["node"])
  at=next(x for x in csv.DictReader((src/"attempts.csv").open()) if int(x["node"])==anchor["node"] and int(x["attempt"])==int(sel["attempt"]))
  prof=[x for x in csv.DictReader((src/"profiles.csv").open()) if int(x["node"])==anchor["node"] and int(x["attempt"])==int(sel["attempt"])]
  old={"dimension":int(at["dimension"]),"radius":float(at["radius"]),"n_b":float(sel["n_b"]),
       "profile":[[float(x[k]) for k in ["nn","np"]] for x in prof]}
  meta={"kind":"retained_native_profile","source":rel(src),"source_node":anchor["node"],"source_selected":sel,"source_attempt":at,
   "source_sha256":{f:sha(src/f) for f in ["selected.csv","attempts.csv","profiles.csv"]}}
  return old,meta
 src=FBASE/"sources"/anchor["source_model"]
 n=list(map(float,(src/"eos.nb").read_text().split()[2:]))
 rows=[line.split() for line in (src/"eos.thermo").read_text().splitlines()[1:] if line.strip()]
 idx=min(range(len(n)),key=lambda i:abs(n[i]-anchor["density"]))
 row=rows[idx];assert int(row[10])==7 and int(row[13])==3
 radius=float(row[17])
 return None,{"kind":"public_radius_start_only","source":rel(src),"source_row_zero_based":idx,
  "source_density":n[idx],"radius":radius,"target_anchor_density":anchor["density"],
  "source_Q_tokens":row[3:10],"auxiliary_tokens":row[11:],
  "source_sha256":{f:sha(src/f) for f in ["eos.nb","eos.thermo","eos.compo"]},
  "role":"Only the initial radius is informed by source; all thermodynamic states are new native field solutions."}

def cached_run(exe,config,n,points,old,seed_meta,work,parent):
 p=work/"STATE.json"
 if p.exists():
  rec=load(p)
  assert rec["executable_sha256"]==sha(exe) and rec["model"]==config["model"] and rec["points"]==points and rec["n_b"]==float(n)
  if rec["status"]=="RUNNING": raise RuntimeError("Unfinished native run must be inspected before resume: "+str(p))
  return rec
 kwargs={}
 if old is not None:
  nn,pp,radius=interpolate_initial(old,points,float(n))
  kwargs.update(initial_nn=nn,initial_np=pp)
 else: radius=seed_meta["radius"]*(seed_meta["target_anchor_density"]/float(n))**(1/3)
 return run_state(exe,config["model"],float(n),points,3,"droplet",radius,work_dir=work,parent=parent,
                  maxiter=1800,ftol=1e-14,**kwargs)

def solve_point(exe,config,n,points,initial,seed_meta,destination,parent):
 first=cached_run(exe,config,n,points,initial,seed_meta,destination/"primary",parent)
 attempts=[first]
 chosen=first if first.get("accepted_stationary") else None
 if chosen is None and first.get("status")=="REQUIRES_REFINEMENT" and first.get("optimizer_mode") in (8.,9.) and "profile" in first:
  if first["kkt"]<.02 and first["radius_stationarity"]<2e-4 and first["radius_bound"]==0 and first["density_bound"]==0:
   second=cached_run(exe,config,n,points,first,seed_meta,destination/"warm_restart",
       {"operation":"single optimizer warm restart from returned finite profile; unchanged ftol/physical acceptance",
        "state_path":first["state_path"],"state_sha256":sha(first["state_path"])})
   attempts.append(second)
   if second.get("accepted_stationary"):chosen=second
 rec={"n_b":float(n),"points":points,"model":config["model"],"accepted":chosen is not None,
      "state_path":chosen["state_path"] if chosen else None,
      "attempts":[{"state_path":a["state_path"],"state_sha256":sha(a["state_path"]),"status":a["status"],
                   "accepted_stationary":a["accepted_stationary"],"actual_exit_code":a.get("actual_exit_code")} for a in attempts]}
 save(destination/"POINT.json",rec)
 return chosen,rec

def regenerate():
 BASE.mkdir(parents=True,exist_ok=True)
 exe=build_driver(ROOT,BASE/"build/production")
 registry=ROOT/"fortran_eos/src/compose_eos_tf_models.f90"
 source_files=[Path(__file__),ROOT/"python_frontend/compose_eos/tf_candidate.py",registry,ROOT/"fortran_eos/tools/tf_candidate_driver.f90"]
 identity={rel(p):sha(p) for p in source_files}
 plan={"schema":"ns_cowling.G_support_regeneration_plan.v1","created_at_utc":datetime.now(timezone.utc).isoformat(),
  "configurations":CONFIGS,"meshes":MESHES,"source_sha256":identity,"executable_sha256":sha(exe),
  "branch_policy":"Frozen source anchors and full predetermined density windows; no partner-based trimming. All accepted contiguous runs are diagnosed afterwards.",
  "repair_policy":"At most one optimizer reset using its returned finite profile for mode8/9, with unchanged ftol1e-14 and native acceptance.",
  "public_sources":"Original public tables remain unchanged; their native counterparts are explicitly relabelled as regenerated sources."}
 if (BASE/"PLAN.json").exists():
  old=load(BASE/"PLAN.json")
  assert old["source_sha256"]==identity and old["configurations"]==CONFIGS and old["meshes"]==MESHES
 else:save(BASE/"PLAN.json",plan)
 runs=[];start=time.monotonic()
 for config in CONFIGS:
  count=round((config["bounds"][1]-config["bounds"][0])/config["step"])+1
  grid=np.unique(np.r_[np.linspace(*config["bounds"],count),config.get("extra_density",[])])
  for anchor in config["anchors"]:
   original,meta=seed_for(config,anchor)
   a_density=original["n_b"] if original else meta["target_anchor_density"]
   anchor_index=int(np.argmin(abs(grid-a_density)))
   branch=BASE/"states"/config["id"]/anchor["name"]
   save(branch/"SEED.json",meta)
   prior_anchor=original
   for points in MESHES:
    output=branch/f"N{points}";records=[None]*len(grid)
    order=[anchor_index]+list(range(anchor_index+1,len(grid)))+list(range(anchor_index-1,-1,-1))
    last=prior_anchor;anchor_state=None
    for index in order:
     if index==anchor_index-1:last=anchor_state or prior_anchor
     parent={"seed":rel(branch/"SEED.json"),"seed_sha256":sha(branch/"SEED.json"),
       "previous_state":last.get("state_path") if isinstance(last,dict) else None,
       "initial_source":"same declared continuation branch; interpolated profiles are initial guesses only"}
     chosen,point=solve_point(exe,config,grid[index],points,last,meta,output/f"point_{index:03d}",parent)
     point["grid_index"]=index;records[index]=point
     if chosen is not None:last=chosen
     if index==anchor_index:anchor_state=chosen
    if anchor_state is not None:prior_anchor=anchor_state
    record={"model":config["model"],"config_id":config["id"],"anchor":anchor["name"],"points":points,
      "density_grid_fm3":grid.tolist(),"points_results":records,"seed":meta,
      "accepted_points":sum(p["accepted"] for p in records),"requested_points":len(grid)}
    save(output/"RUN.json",record);runs.append(record)
    save(BASE/"PROGRESS.json",{"completed_mesh_branches":len(runs),"last_model":config["model"],"last_anchor":anchor["name"],"last_mesh":points,
      "accepted_states":sum(r["accepted_points"] for r in runs),"requested_states_so_far":sum(r["requested_points"] for r in runs),"wall_seconds":time.monotonic()-start})
    print(json.dumps({"model":config["model"],"anchor":anchor["name"],"mesh":points,"accepted":record["accepted_points"],"requested":len(grid)},ensure_ascii=False),flush=True)
 assert {rel(p):sha(p) for p in source_files}==identity,"Production sources changed during run"
 result={"schema":"ns_cowling.G_support_regeneration.v1","plan_sha256":sha(BASE/"PLAN.json"),"source_sha256":identity,"executable_sha256":sha(exe),
  "runs":runs,"wall_seconds":time.monotonic()-start,"actual_native_invocations":sum(len(p["attempts"]) for r in runs for p in r["points_results"]),
  "accepted_grid_states":sum(r["accepted_points"] for r in runs),"requested_grid_states":sum(r["requested_points"] for r in runs),
  "source_bytes_unchanged_during_production":True}
 save(BASE/"REGENERATION.json",result)
 print(json.dumps({k:v for k,v in result.items() if k not in ["runs","source_sha256"]},ensure_ascii=False),flush=True)

def accepted_runs(points):
 runs=[];start=None
 for i,p in enumerate(points+[{"accepted":False}]):
  if p["accepted"] and start is None:start=i
  if not p["accepted"] and start is not None:
   if i-start>=3:runs.append((start,i))
   start=None
 return runs

def match_candidate(model,core_record,low,high,lc,hc,destination,source_record):
 destination.mkdir(parents=True,exist_ok=False)
 assert any(c["model"]==model and c["core_id"]==core_record["id"] for c in CONFIGS)
 result={"status":"STARTED","kernel_invoked":False,"model":model,"core_id":core_record["id"],
         "crust_component":lc,"core_component":hc,"source_record":source_record["id"]}
 try:
  from compose_eos.maxwell_v2 import construct_maxwell_transition_v2
  ld=_segment_eosdata(low,(lc["phase_id"],lc["start_zero_based"],lc["stop_exclusive"]))
  hd=_segment_eosdata(high,(hc["phase_id"],hc["start_zero_based"],hc["stop_exclusive"]))
  with observe_kernel(result):
   matched,report=construct_maxwell_transition_v2(ld,hd,assume_local_charge_neutrality=True,return_report=True)
  native=report.native_report;result["v2_report"]=plain(report)
  matched.write_internal_directory(destination/"expected_internal")
  write_compose(destination/"compose",matched);back=read_compose(destination/"compose")
  back.write_internal_directory(destination/"readback_internal")
  numerical=verify_internal_readback(destination/"expected_internal/eos_data.csv",destination/"readback_internal/eos_data.csv",expected_n=matched.n,actual_n=back.n)
  metadata=verify_internal_metadata_readback(destination/"expected_internal/eos_data.nml",destination/"readback_internal/eos_data.nml")
  result["readback"]={"numerical":numerical,"metadata":metadata}
  if numerical["status"]!="PASS" or metadata["status"]!="PASS":raise ValueError("Actual serialized readback failed")
  result["admission"]=admit_match_report("maxwell",native,destination/"readback_internal/eos_data.csv",destination/"readback_internal/eos_data.nml")
  if result["admission"]["status"]!="ADMITTED":raise ValueError("Independent interface admission rejected")
  result["serialized_native_verification"]=run_native("verify",destination/"compose",destination/"CONSISTENCY.json")
  if result["serialized_native_verification"].get("overall")!="PASS":raise ValueError("Actual output consistency rejected")
  result.update(status="GENERATED_VERIFIED_FINITE_DOMAIN",output_path=rel(destination/"compose"),
    output_sha256={p.name:sha(p) for p in (destination/"compose").glob("eos.*")},
    transition={k:native[k] for k in ["mu_transition","pressure_transition","n_low","n_high","epsilon_low","epsilon_high","density_jump","latent_energy_density","maxwell_identity_residual"]},
    density_fm3=[float(matched.nb[0]),float(matched.nb[-1])],
    scope="Explicitly regenerated native crust; finite connected thermodynamic source support; no public table replacement or global microphysical minimum claim")
 except Exception as exc:
  result.update(status="REJECTED",exception_type=type(exc).__name__,error=str(exc),category=getattr(exc,"category",None) or "EXECUTION_OR_VALIDATION_FAILURE")
  if hasattr(exc,"report"):result["v2_report"]=plain(exc.report)
 save(destination/"RESULT.json",plain(result))
 return plain(result)

def match():
 assert sha(executable())==MATCH_EXE_SHA
 records={r["id"]:r for r in load(FBASE/"prepared/MANIFEST.json")["source_records"]}
 generated=load(BASE/"REGENERATION.json")
 cells=[];prepared=[];configs={c["model"]:c for c in CONFIGS}
 for run in generated["runs"]:
  config=configs[run["model"]];core=records[config["core_id"]];high=read_prepared(core)
  cell={"model":run["model"],"anchor":run["anchor"],"mesh":run["points"],
        "requested_grid_points":run["requested_points"],"accepted_grid_points":run["accepted_points"],
        "core_id":core["id"],"generated_native_crust":True,"candidates":[]}
  for j,(a,b) in enumerate(accepted_runs(run["points_results"])):
   states=[load(p["state_path"]) for p in run["points_results"][a:b]]
   target=BASE/"prepared"/run["model"]/run["anchor"]/f"N{run['points']}"/f"accepted_run_{j:03d}"
   if (target/"GENERATION.json").exists():
    receipt=load(target/"GENERATION.json")
    for s in receipt["source_states"]:assert sha(s["path"])==s["sha256"]
    for f,h in receipt["table_sha256"].items():assert sha(target/"compose"/f)==h
    data=read_compose(target/"compose")
   else:data,receipt=export_compose(target,states,mn=config["mn"],mp=config["mp"],model=config["model"])
   rp=target/"diagnostic/SOURCE.json"
   if rp.exists():record=load(rp)
   else:
    record={"id":f"native_regenerated:{run['model']}:{run['anchor']}:N{run['points']}:run{j}",
     "origin":"native_regenerated","family":"native_regenerated","model":run["model"],"role":"crust",
     "resolution":run["points"],"source_path":rel(target/"compose"),"source_table_sha256":receipt["table_sha256"],
     "selection":{"source_start_zero_based":0,"source_stop_exclusive":len(states),"basis":"Every consecutive accepted node of the predetermined density grid; failures split support"},
     "grid_indices":[a,b],"generation_receipt":rel(target/"GENERATION.json"),"generation_receipt_sha256":sha(target/"GENERATION.json")}
    record=prepare_record(data,record,target/"diagnostic")
   prepared.append(record)
   for lc in record["components"]:
    for hc in core["components"]:
     dest=BASE/"matching"/run["model"]/run["anchor"]/f"N{run['points']}"/f"run{j}_{lc['id']}__{hc['id']}"
     result=load(dest/"RESULT.json") if (dest/"RESULT.json").exists() else match_candidate(run["model"],core,data,high,lc,hc,dest,record)
     cell["candidates"].append(result)
  cell["successes"]=sum(x["status"]=="GENERATED_VERIFIED_FINITE_DOMAIN" for x in cell["candidates"])
  cells.append(cell)
  print(json.dumps({k:v for k,v in cell.items() if k!="candidates"}),flush=True)
 result={"schema":"ns_cowling.G_support_matching.v1","regeneration_sha256":sha(BASE/"REGENERATION.json"),
  "F_prepared_manifest_sha256":sha(FBASE/"prepared/MANIFEST.json"),"matching_executable_sha256":sha(executable()),
  "cells":cells,"actual_kernel_invocations":sum(x["kernel_invoked"] for c in cells for x in c["candidates"]),
  "verified_outputs":sum(c["successes"] for c in cells),"physical_pair_status":{m:any(c["successes"] for c in cells if c["model"]==m) for m in configs},
  "rejections":dict(Counter(x.get("category") for c in cells for x in c["candidates"] if x["status"]=="REJECTED"))}
 save(BASE/"PREPARED_MANIFEST.json",{"source_records":prepared})
 save(BASE/"MATCHING.json",result)
 print(json.dumps({k:v for k,v in result.items() if k!="cells"},ensure_ascii=False),flush=True)

if __name__=="__main__":
 parser=argparse.ArgumentParser();parser.add_argument("action",choices=["diagnose","regenerate","match","all"])
 args=parser.parse_args()
 if args.action=="diagnose":diagnose()
 if args.action in ["regenerate","all"]:regenerate()
 if args.action in ["match","all"]:match()
