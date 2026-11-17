#!/usr/bin/env python3
"""H: bounded 1.4-Msun TOV inventory for frozen opt-in causal closures."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,os,sys,time,traceback
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"data/low_density_crust_20261117/stellar/causal_extension"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def rel(p):return str(Path(p).relative_to(ROOT))
def save(p,x):
 Path(p).parent.mkdir(parents=True,exist_ok=True)
 Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def stamp():return datetime.now(timezone.utc).isoformat()
def main(plan_path):
 plan_path=Path(plan_path);plan=load(plan_path)
 if len(plan["cases"])>8:raise ValueError("at most eight frozen causal contracts")
 if (BASE/"TOV_RESULTS.json").exists():raise ValueError("TOV results already exist; no evidence overwrite")
 start=time.monotonic();receipt={"status":"RUNNING","started_at_utc":stamp(),"pid":os.getpid(),"actual_exit_code":None,
   "plan":rel(plan_path),"plan_sha256":sha(plan_path)}
 save(BASE/"TOV_PROCESS.json",receipt)
 try:
  expected=plan["source_and_code_sha256"]
  for p,h in expected.items():
   if sha(ROOT/p)!=h:raise ValueError("frozen source/code changed: "+p)
  sys.path.insert(0,str(ROOT/"python_frontend"))
  from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background
  files=[ROOT/p for p in expected]+[plan_path,ROOT/plan["admission"]]+[ROOT/c["contract"] for c in plan["cases"]]
  before={rel(p):sha(p) for p in files}
  for c in plan["cases"]:
   if sha(ROOT/c["contract"])!=c["contract_sha256"]:raise ValueError("frozen contract changed")
  result=dict(schema="ns_cowling.H_causal_background_inventory.v1",status="RUNNING",scope=plan["scope"],
    plan=rel(plan_path),plan_sha256=sha(plan_path),target_mass_msun=plan["mass_msun"],
    settings=plan["settings"],source_hashes_before=before,cases=[],numerical_sensitivities=[],
    radial_stability_claim=False,complete_mode_spectrum_computed=False)
  save(BASE/"TOV_RESULTS.json",result)
  primaries=[];all_calls=0
  def structure(bg):
   eos=bg.eos;b=eos.barotrope
   hl=float(np.log(b.bridge.mu_left/eos.a));hr=float(np.log(b.bridge.mu_right/eos.a))
   rl,ml=map(float,bg.solution.sol(hl));rr,mr=map(float,bg.solution.sol(hr))
   continuity=[]
   for interface in bg.interfaces:
    hi=bg.solution.sol(interface["h"],side="high")
    lo=bg.solution.sol(interface["h"],side="low")
    continuity.append(dict(h=interface["h"],r_jump_km=float(lo[0]-hi[0]),m_jump_geom_km=float(lo[1]-hi[1]),
      two_sided_dense_output=True,source_density_jump_fm3=interface["n_high_fm3"]-interface["n_low_fm3"]))
   return dict(background=bg.summary(),lower_bridge_radial_width_km=rl-rr,
      outer_and_surface_radial_width_km=bg.radius_km-rl,
      sharp_crust_radial_width_km=bg.radius_km-bg.interfaces[0]["r_km"] if bg.interfaces else None,
      lower_bridge_boundaries=dict(outer_radius_km=rl,upper_radius_km=rr,outer_mass_geom_km=ml,upper_mass_geom_km=mr),
      interface_continuity=continuity)
  for case in plan["cases"]:
   case_start=time.monotonic();dest=BASE/"calculations"/case["id"]
   record=dict(**case,status="RUNNING",sequence=[],solver_calls=[],selection=plan["selection"])
   result["cases"].append(record);save(dest/"BACKGROUND.json",record);save(BASE/"TOV_RESULTS.json",result)
   try:
    eos=SharpContractEnthalpyEOS.from_path(ROOT/case["contract"])
    def solve(n,settings,purpose):
     nonlocal all_calls
     bg=solve_sharp_background(eos,float(n),**settings);all_calls+=1
     record["solver_calls"].append(dict(purpose=purpose,central_density_fm3=float(n),mass_msun=bg.mass_msun,
       radius_km=bg.radius_km,nfev=bg.solution.nfev))
     return bg
    for n in plan["density_scan_fm3"]:
     if n>eos.density_domain[1]:
      record["sequence"].append(dict(central_density_fm3=n,status="OUTSIDE_SOURCE_SUPPORT"));continue
     try:
      bg=solve(n,plan["settings"],"fixed density bracket inventory")
      record["sequence"].append(dict(status="SOLVED",**bg.summary()))
     except (ValueError,RuntimeError) as exc:
      record["sequence"].append(dict(central_density_fm3=n,status="FAILED",exception_type=type(exc).__name__,error=str(exc)))
    brackets=[]
    for a,b in zip(record["sequence"][:-1],record["sequence"][1:]):
     if a["status"]=="SOLVED" and b["status"]=="SOLVED" and a["mass_msun"]<plan["mass_msun"]<b["mass_msun"]:
      brackets.append((a["central_density_fm3"],b["central_density_fm3"]))
    record["all_sampled_ascending_target_brackets"]=brackets
    if not brackets:raise ValueError("No sampled ascending target-mass bracket in source support")
    record["selected_mass_bracket"]=brackets[0]
    nc,root=brentq(lambda n:solve(n,plan["settings"],"target-mass root").mass_msun-plan["mass_msun"],
       *brackets[0],xtol=plan["mass_root_xtol"],full_output=True)
    bg=solve(nc,plan["settings"],"final target-mass background")
    record.update(structure(bg));record["mass_root_converged"]=bool(root.converged)
    record["mass_target_residual_msun"]=bg.mass_msun-plan["mass_msun"]
    if not root.converged or abs(record["mass_target_residual_msun"])>plan["mass_target_absolute_tolerance_msun"]:
     raise RuntimeError("Target mass root failed fixed acceptance")
    if any(abs(x["r_jump_km"])>1e-12 or abs(x["m_jump_geom_km"])>1e-12 for x in record["interface_continuity"]):
     raise RuntimeError("Two-sided dense TOV interface geometry is discontinuous")
    record.update(status="PASS_SCOPED",wall_seconds=time.monotonic()-case_start)
    if case["label"]=="primary":primaries.append((case,eos,bg,record))
    print(json.dumps({"event":"BACKGROUND","id":case["id"],"mass":bg.mass_msun,"radius":bg.radius_km,
      "nc":nc,"bridge_km":record["lower_bridge_radial_width_km"],"crust_km":record["sharp_crust_radial_width_km"]}),flush=True)
   except (ValueError,RuntimeError) as exc:
    record.update(status="FAILED",exception_type=type(exc).__name__,error=str(exc),wall_seconds=time.monotonic()-case_start)
    print(json.dumps({"event":"CASE_FAILED","id":case["id"],"error":str(exc)}),flush=True)
   save(dest/"BACKGROUND.json",record);save(BASE/"TOV_RESULTS.json",result)
  for case,eos,bg,record in primaries:
   tight_start=time.monotonic()
   control=dict(source_id=case["source_id"],primary_case_id=case["id"],fixed_central_density=True,
       central_density_fm3=bg.central_density_fm3,settings=plan["tighter_settings"])
   try:
    tight=solve_sharp_background(eos,bg.central_density_fm3,**plan["tighter_settings"]);all_calls+=1
    control.update(structure(tight))
    control.update(status="PASS_SCOPED",mass_difference_msun=tight.mass_msun-bg.mass_msun,
      radius_difference_km=tight.radius_km-bg.radius_km,
      lower_bridge_width_difference_km=control["lower_bridge_radial_width_km"]-record["lower_bridge_radial_width_km"],
      sharp_crust_width_difference_km=control["sharp_crust_radial_width_km"]-record["sharp_crust_radial_width_km"],
      wall_seconds=time.monotonic()-tight_start)
   except (ValueError,RuntimeError) as exc:
    control.update(status="FAILED",exception_type=type(exc).__name__,error=str(exc),wall_seconds=time.monotonic()-tight_start)
   result["numerical_sensitivities"].append(control)
   save(BASE/"calculations"/case["source_id"]/"NUMERICAL_SENSITIVITY.json",control)
   save(BASE/"TOV_RESULTS.json",result)
   print(json.dumps({"event":"TIGHTER_CONTROL","source":case["source_id"],"status":control["status"],
      "delta_R_km":control.get("radius_difference_km")}),flush=True)
  variations=[]
  for source in sorted({c["source_id"] for c in result["cases"]}):
   group=[c for c in result["cases"] if c["source_id"]==source and c["status"]=="PASS_SCOPED"]
   if len(group)==2:
    a=next(c for c in group if c["label"]=="primary");b=next(c for c in group if c["label"]!="primary")
    variations.append(dict(source_id=source,interpretation="same target mass; different endpoint-selected phenomenological lower closure",
      delta_radius_km=b["background"]["radius_km"]-a["background"]["radius_km"],
      delta_central_density_fm3=b["background"]["central_density_fm3"]-a["background"]["central_density_fm3"],
      delta_bridge_width_km=b["lower_bridge_radial_width_km"]-a["lower_bridge_radial_width_km"],
      delta_sharp_crust_width_km=b["sharp_crust_radial_width_km"]-a["sharp_crust_radial_width_km"]))
  result["next_join_sensitivities"]=variations
  result["source_hashes_after"]={rel(p):sha(p) for p in files}
  result["source_bytes_unchanged"]=before==result["source_hashes_after"]
  result["actual_TOV_solver_calls_completed"]=all_calls
  result["status"]="PASS_SCOPED" if result["cases"] and all(c["status"]=="PASS_SCOPED" for c in result["cases"]+result["numerical_sensitivities"]) and result["source_bytes_unchanged"] else "PARTIAL_OR_REJECTED"
  result["wall_seconds"]=time.monotonic()-start
  save(BASE/"TOV_RESULTS.json",result)
  exitcode=0 if result["status"]=="PASS_SCOPED" else 1
  receipt.update(status="COMPLETED" if exitcode==0 else "COMPLETED_WITH_REJECTIONS",actual_exit_code=exitcode,
     finished_at_utc=stamp(),wall_seconds=result["wall_seconds"],results_sha256=sha(BASE/"TOV_RESULTS.json"),
     source_bytes_unchanged=result["source_bytes_unchanged"],actual_TOV_solver_calls_completed=all_calls)
  save(BASE/"TOV_PROCESS.json",receipt)
  print(json.dumps({"event":"FINISHED","status":result["status"],"calls":all_calls,"wall_seconds":result["wall_seconds"]}),flush=True)
  return exitcode
 except BaseException as exc:
  receipt.update(status="FAILED",actual_exit_code=1,finished_at_utc=stamp(),wall_seconds=time.monotonic()-start,
    exception_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
  save(BASE/"TOV_PROCESS.json",receipt);raise
if __name__=="__main__":
 a=argparse.ArgumentParser();a.add_argument("plan",nargs="?",default=str(BASE/"TOV_PLAN.json"));args=a.parse_args();raise SystemExit(main(args.plan))
