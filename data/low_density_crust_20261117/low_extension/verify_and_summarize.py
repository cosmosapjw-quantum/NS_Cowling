#!/usr/bin/env python3
"""Finite evidence manifest for the H native low-density lane; no new solver calls."""
from pathlib import Path
import hashlib,json
from collections import Counter
ROOT=Path(__file__).resolve().parents[3]
BASE=Path(__file__).resolve().parent
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(p):return str(Path(p).relative_to(ROOT))
def save(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2)+"\n")
plan=load(BASE/"PLAN.json")
assert all(sha(ROOT/p)==h for p,h in plan["source_sha256"].items())
assert all(sha(p)==h for p,h in plan["executables"].items())
assert all(sha(s["path"])==s["sha256"] for states in plan["source_states"].values() for s in states)
process=list(BASE.glob("PROCESS_*.json"))+list((BASE/"repaired").glob("*/PROCESS.json"))
assert all(load(p)["status"]=="COMPLETED" and load(p)["actual_exit_code"]==0 for p in process)
models=[];all_attempts={}
fields=("mu_transition","pressure_transition","n_low","n_high","epsilon_low","epsilon_high","density_jump","latent_energy_density","maxwell_identity_residual")
for model,cfg in plan["configurations"].items():
 original=load(BASE/f"{model}_N128_STATES.json")
 points=load(BASE/f"{model}_POINTS.json")
 repair_base=BASE/"repaired"/model
 final=load(repair_base/"FINAL_N128_STATES.json")
 repair=load(repair_base/"REPAIRS.json")
 for p in points:
  for a in p["attempts"]:
   assert sha(a["path"])==a["sha256"]
   all_attempts[a["path"]]=load(a["path"])
 for a in repair:
  if "repaired_state_path" in a:
   assert sha(a["repaired_state_path"])==a["repaired_state_sha256"]
   all_attempts[a["repaired_state_path"]]=load(a["repaired_state_path"])
 smpath=BASE/"repaired/tables_v2"/model/"N128/SUMMARY.json"
 sm=load(smpath)
 successes=[]
 for bi,b in enumerate(sm["blocks"]):
  assert not b["interval_rejections"]
  for m in b["matches"]:
   if m["status"]!="GENERATED_VERIFIED_FINITE_DOMAIN":continue
   assert m["admission"]["status"]=="ADMITTED" and m["readback"]["status"]=="PASS" and m["metadata"]["status"]=="PASS"
   assert m["consistency"]["overall"]=="PASS"
   cp=Path(m["compose_path"])
   assert all(sha(cp/f)==h for f,h in m["table_sha256"].items())
   successes.append({"compose_path":rel(cp),"result_path":rel(cp.parent/"RESULT.json"),
     "result_sha256":sha(cp.parent/"RESULT.json"),"table_sha256":m["table_sha256"],
     "connected_source_density_fm3":m["crust_component"]["density"],
     "source_component_knots":m["crust_component"]["stop"]-m["crust_component"]["start"],
     "transition":{k:m["native_report"][k] for k in fields},"serialized_native_verification":"PASS",
     "independent_admission":"ADMITTED","readback":"PASS","physical_energy_shift_mev":m["native_report"]["physical_energy_shift_mev"]})
 assert len(successes)==1
 baseline=None
 if model!="DDME2":
  for p in (ROOT/"data/crust_refinement_20261117/support_repair/matching"/model/"provider_radius_seed/N64").glob("*/RESULT.json"):
   d=load(p)
   if d["status"]=="GENERATED_VERIFIED_FINITE_DOMAIN":
    assert baseline is None
    pt64=d["transition"]["pressure_transition"];pt128=successes[0]["transition"]["pressure_transition"]
    baseline={"path":rel(p),"sha256":sha(p),"mesh":64,"transition":d["transition"],
      "N64_to_N128_abs_relative_pressure_change":abs(pt128-pt64)/abs(pt64)}
 sources=[]
 for s in final:
  if s["accepted_stationary"]:
   assert s["points"]==128
   sources.append({"n_b":s["n_b"],"path":s["state_path"],"sha256":sha(s["state_path"]),
      "epsilon_MeV_fm3":s["epsilon"],"pressure_MeV_fm3":s["pressure"],"mu_B_MeV":s["mu_B"]})
 models.append({"model":model,"requested_knots":len(final),"reused_G_N128_knots":sum(p.get("reused_exact_G",False) for p in points),
   "first_pass_accepted":sum(s["accepted_stationary"] for s in original),"final_accepted":len(sources),
   "failed_source_densities":[s["n_b"] for s in final if not s["accepted_stationary"]],
   "state_index":rel(repair_base/"FINAL_N128_STATES.json"),"state_index_sha256":sha(repair_base/"FINAL_N128_STATES.json"),
   "accepted_sources":sources,"maxwell":successes[0],"N64_baseline":baseline,
   "production_native_invocations":sum(len(p["attempts"]) for p in points),
   "repair_native_invocations":sum("repaired_state_path" in p for p in repair)})
assert all(s["actual_exit_code"]==0 for s in all_attempts.values())
out={"schema":"ns_cowling.H_low_extension_finite_evidence.v1","plan_path":rel(BASE/"PLAN.json"),"plan_sha256":sha(BASE/"PLAN.json"),
 "verifier_sha256":sha(__file__),"source_identity_verification":"PASS","process_receipts":[{"path":rel(p),"sha256":sha(p),**load(p)} for p in process],
 "native_invocations":len(all_attempts),"native_status_counts":dict(Counter(s["status"] for s in all_attempts.values())),
 "native_exit_code_counts":dict(Counter(str(s["actual_exit_code"]) for s in all_attempts.values())),
 "first_pass_requested_knots":sum(m["requested_knots"] for m in models),"reused_G_knots":sum(m["reused_G_N128_knots"] for m in models),
 "final_accepted_knots":sum(m["final_accepted"] for m in models),"models":models,
 "limitations":["Three numerical optimizer failures remain despite one returned-profile warm restart and one accepted opposite-direction seed.",
  "Their density gaps split source support; no failed point or alternate spatial mesh is substituted.",
  "Finite spherical N128 branch evidence is not a continuum or global-ground-state proof.",
  "Full-EOS outer-crust completion and stellar admission are separate downstream root-owned results."]}
save(BASE/"MANIFEST.json",out)
print(json.dumps({"native_invocations":out["native_invocations"],"native_status_counts":out["native_status_counts"],
 "final_accepted_knots":out["final_accepted_knots"],"requested_knots":out["first_pass_requested_knots"],
 "manifest_sha256":sha(BASE/"MANIFEST.json"),
 "models":[{"model":m["model"],"accepted":m["final_accepted"],"failed":m["failed_source_densities"],
 "production_calls":m["production_native_invocations"],"repair_calls":m["repair_native_invocations"],
 "Pt":m["maxwell"]["transition"]["pressure_transition"],"N64_relative":(m["N64_baseline"] or {}).get("N64_to_N128_abs_relative_pressure_change")} for m in models]},indent=2))
