from pathlib import Path
import json,hashlib,datetime
import numpy as np
ROOT=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117");D=ROOT/"data/low_density_crust_20261117";B=D/"stellar/causal_extension";O=D/"review"
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rawtable(path):
 path=Path(path);n=np.loadtxt(path/"eos.nb",skiprows=2);lines=(path/"eos.thermo").read_text().splitlines()
 mn=float(lines[0].split()[0]);q=np.array([[float(v) for v in l.split()] for l in lines[1:] if l.strip()])
 return dict(n=n,p=n*q[:,3],epsilon=n*mn*(q[:,9]+1))
errors=[]
def check(ok,msg):
 if not ok:errors.append(msg)
admission=load(B/"FULL_CONTRACT_ADMISSION.json");plan=load(B/"TOV_PLAN.json");run=load(B/"TOV_PROCESS.json");res=load(B/"TOV_RESULTS.json");base=load(ROOT/admission["reference"])
check(admission["source_bytes_unchanged"] and admission["stellar_observables_computed"] is False,"admission state")
check(sha(B/"FULL_CONTRACT_ADMISSION.json")==plan["admission_sha256"],"admission hash")
check(sha(B/"TOV_RESULTS.json")==run["results_sha256"],"result hash")
check(run["actual_exit_code"]==0 and res["status"]=="PASS_SCOPED","TOV completion")
check(datetime.datetime.fromisoformat(plan["created_at_utc"])<datetime.datetime.fromisoformat(run["started_at_utc"]),"plan after observables")
check(len(plan["cases"])==len(res["cases"])==8 and len(res["numerical_sensitivities"])==4,"case count")
for path,h in res["source_hashes_after"].items():check(sha(ROOT/path)==h,"source identity "+path)
check(res["source_hashes_before"]==res["source_hashes_after"],"source bytes changed")
rows=[];upper_residual=0.
for entry in admission["inputs"]:
 source=rawtable(ROOT/entry["compose_path"]);passed=[t["source_row"] for t in entry["lower_join_trials"] if t["status"]=="ADMITTED"]
 check(passed==[c["source_row_start"] for c in entry["admitted"]],entry["id"]+" prechoice mismatch")
 for selected in entry["admitted"]:
  p=ROOT/selected["contract"];c=load(p);check(sha(p)==selected["contract_sha256"],entry["id"]+" contract hash")
  check(c["surface"]==base["surface"] and c["outer"]==base["outer"],entry["id"]+" E source law changed")
  prov=c["provenance"];check(prov["disconnected_low_states_imported"] is False and prov["microscopic_support_claim_for_bridge"] is False and prov["closure_does_not_repair_native_failures"] is True,entry["id"]+" wrong closure semantics")
  check(c["lower_bridge"]["kind"]=="causal_exponential_n_over_mu" and c["lower_bridge"]["validation"]["status"]=="PASS",entry["id"]+" bridge law")
  j=selected["source_row_start"]
  for k in ["n","p","epsilon"]:
   arr=np.asarray(c["upper"][k]);original=source[k][j:];check(arr.shape==original.shape,entry["id"]+" upper shape")
   err=float(np.max(abs(arr-original)/np.maximum(1.,abs(original))));upper_residual=max(upper_residual,err);check(err<5e-13,entry["id"]+" upper "+k+" changed")
  check(c["upper"]["sharp_interfaces"]==[entry["original_interface_index"]-j],entry["id"]+" interface index")
for c in res["cases"]:
 check(c["status"]=="PASS_SCOPED" and c["mass_root_converged"],c["id"]+" status")
 check(abs(c["mass_target_residual_msun"])<=plan["mass_target_absolute_tolerance_msun"],c["id"]+" mass error")
 check(len(c["interface_continuity"])==1 and all(x["two_sided_dense_output"] and x["r_jump_km"]==0 and x["m_jump_geom_km"]==0 for x in c["interface_continuity"]),c["id"]+" interface geometry")
 lo,hi=c["selected_mass_bracket"];check(lo<c["background"]["central_density_fm3"]<hi,c["id"]+" density bracket")
 rows.append(dict(id=c["id"],mass_msun=c["background"]["mass_msun"],radius_km=c["background"]["radius_km"],lower_join=c["upper_start_density"],kappa=c["kappa"],mass_target_residual_msun=c["mass_target_residual_msun"]))
for c in res["numerical_sensitivities"]:
 check(c["status"]=="PASS_SCOPED" and c["fixed_central_density"],c["source_id"]+" numerical control")
 check(all(x["r_jump_km"]==0 and x["m_jump_geom_km"]==0 for x in c["interface_continuity"]),c["source_id"]+" tighter interface")
report=dict(status="PASS" if not errors else "FAIL",scope="Readback of frozen phenomenological whole-EOS contracts, exact source identities and actual TOV receipts; no additional stellar solve",
 admitted=8,backgrounds=8,tighter_controls=4,actual_TOV_solver_calls_completed=res["actual_TOV_solver_calls_completed"],actual_exit_code=run["actual_exit_code"],
 source_upper_max_relative_difference=upper_residual,max_abs_target_mass_error_msun=max(abs(x["mass_target_residual_msun"]) for x in rows),
 largest_tighter_radius_difference_km=max(abs(x["radius_difference_km"]) for x in res["numerical_sensitivities"]),cases=rows,
 next_join_sensitivities=res["next_join_sensitivities"],numerical_radius_differences={x["source_id"]:x["radius_difference_km"] for x in res["numerical_sensitivities"]},
 evidence_sha256={str(p.relative_to(ROOT)):sha(p) for p in [B/"FULL_CONTRACT_ADMISSION.json",B/"TOV_PLAN.json",B/"TOV_PROCESS.json",B/"TOV_RESULTS.json"]},errors=errors)
(O/"CAUSAL_TOV_REVIEW.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2));raise SystemExit(bool(errors))
