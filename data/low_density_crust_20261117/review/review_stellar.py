"""Final readback within the single H independent review."""
from pathlib import Path
import json,hashlib,math
R=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117");D=R/"data/low_density_crust_20261117";O=D/"review";B=D/"stellar/ddmex_comparison"
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
errors=[]
def check(ok,msg):
 if not ok:errors.append(msg)
raw=load(B/"STELLAR_RESULTS.json");plan=load(B/"PLAN.json");run=load(B/"RUN.json");mode=load(B/"MODE_ACCEPTANCE.json");fix=load(B/"RUNNER_SETTINGS_CORRECTION.json");summary=load(D/"H_SUMMARY.json")
finalize=load(B/"FINALIZE_PROCESS.json")
check(finalize.get("actual_exit_code",finalize.get("returncode"))==0,"finalizer execution exit")
check(raw["source_bytes_unchanged"] and raw["source_hashes_before"]==raw["source_hashes_after"],"historical execution source changed")
check(raw["plan_sha256"]==sha(B/"PLAN.json"),"plan identity")
check(mode["raw_results_sha256"]==sha(B/"STELLAR_RESULTS.json"),"raw results modified")
check(mode["runner_correction_receipt_sha256"]==sha(B/"RUNNER_SETTINGS_CORRECTION.json"),"correction receipt changed")
check(sha(R/fix["executed_source"])==fix["executed_source_sha256"],"archive snapshot changed")
check(sha(R/fix["final_runner"])==fix["final_runner_sha256"],"corrected runner changed")
old_runner=str(R/fix["final_runner"])
check(raw["source_hashes_before"][old_runner]==fix["executed_source_sha256"],"executed source not original runner")
for path,h in raw["source_hashes_before"].items():
 actual=R/fix["executed_source"] if path==old_runner else Path(path)
 check(sha(actual)==h,"execution input bytes changed: "+path)
new=(R/fix["final_runner"]).read_text()
check("repair_cfg=bounded_repair_settings(tight,cfg)" in new and "solve_mode(bg,bracket,repair_cfg" in new,"corrected surface helper not applied at solve")
check('mode["effective_settings"]=asdict(repair_cfg)' in new and 'mode["effective_settings"]=asdict(cfg)' in new,"effective configuration absent")
check(len(fix["configuration_regressions"])==3 and all(c["status"]=="PASS" and c["original_immutable"] and c["requested_surface_h"]==c["effective"]["surface_h"] for c in fix["configuration_regressions"]),"configuration regression receipt")
check(len(raw["cases"])==4 and len(raw["numerical_refinement"]["modes"]["slow"])==2 and len(raw["numerical_refinement"]["modes"]["rapid"])==1,"planned matrix incomplete")
records=mode["all_selected_root_records"];check(len(records)==15,"selected record count")
accepted=0;failed=[];mismatches=[];profile_count=0
for row in records:
 p=R/row["selected_file"];m=load(p);orig=R/row["original_file"]
 check(sha(p)==row["selected_file_sha256"],row["id"]+" selected file hash")
 check(sha(orig)==row["original_file_sha256"],row["id"]+" original file hash")
 check(m["frequency_hz"]==row["frequency_hz"] and m["matching_residual"]==row["matching_residual"] and m["relative_join_mismatch"]==row["relative_join_mismatch"],row["id"]+" result field mismatch")
 a=plan["acceptance"]
 passed=bool(m["root_converged"] and abs(m["matching_residual"])<a["determinant_abs"] and m["relative_join_mismatch"]<a["join_relative"] and m["surface_BC_relative"]<a["surface_relative"] and all(j["traction_relative_residual"]<a["traction_relative"] and abs(j["kinematic_residual"])<a["kinematic_abs"] for j in m["interface_junctions"]))
 check(len(m["interface_junctions"])==1,row["id"]+" missing interface transfer")
 check(passed==(row["acceptance"]=="PASS"),row["id"]+" admission differs from raw residual")
 check(all(math.isfinite(v) for q in m["profile"] for v in [q["h"],q["r_km"],q["W"],q["V_km_inverse"]]),row["id"]+" nonfinite profile")
 # Check distinct one-sided rows at the actual sharp interface.
 h=m["interface_junctions"][0]["h"];sides=[q["side"] for q in m["profile"] if q["h"]==h]
 check(sorted(sides)==["high","low"],row["id"]+" profile failed to retain both sides")
 profile_count+=len(m["profile"])
 if not passed:failed.append({k:row[k] for k in ["id","conversion","root_index","frequency_hz","matching_residual","relative_join_mismatch"]})
 if not row["planned_refinement_surface_preserved"]:
  mismatches.append(row)
  check(row["id"]=="numerical_refinement" and not row["used_for_accepted_frequency_comparison"],"unpaired refinement included")
 if row["used_for_accepted_frequency_comparison"]:
  accepted+=1;check(passed and row["planned_refinement_surface_preserved"],"unapproved root used for comparison")
check(accepted==11 and len(failed)==4 and len(mismatches)==1,"admission counts")
check(run["actual_exit_code"]==1 and raw["status"]=="REJECTED_OR_EMPTY" and mode["status"]=="COMPLETED_WITH_REJECTED_CANDIDATES","raw rejection exit misrepresented")
for c in mode["finite_comparisons"]:
 case=next(x for x in raw["cases"] if x["id"]==c["id"]);base=raw["cases"][0]
 for conversion,entries in c["frequencies"].items():
  for x in entries:
   i=x["root_index"];a,b=case["modes"][conversion][i],base["modes"][conversion][i]
   check(a["acceptance"]==b["acceptance"]=="PASS","frequency delta from rejected root")
   check(x["delta_frequency_hz"]==a["frequency_hz"]-b["frequency_hz"],"frequency delta arithmetic")
for path,h in summary["key_files"].items():check(sha(R/path)==h,"summary pointer identity "+path)
check(summary["Cowling"]["accepted"]==accepted and summary["Cowling"]["rejected"]==len(failed) and summary["Cowling"]["raw_process_exit"]==run["actual_exit_code"],"top-level Cowling summary")
check(summary["Cowling"]["new_causal_models_mode_spectrum_computed"] is False and summary["causal_TOV"]["microscopic_lower_gap_repaired"] is False,"scope metadata overclaim")
report=dict(status="PASS_SCOPED_REJECTIONS_PRESERVED" if not errors else "FAIL",ordinary_cases=4,selected_records=15,accepted_records=accepted,rejected_records=len(failed),unpaired_refinement_records=len(mismatches),raw_process_actual_exit_code=run["actual_exit_code"],raw_status=raw["status"],mode_admission_status=mode["status"],finite_profile_points_checked=profile_count,unaccepted_low_candidates=failed,original_code_snapshot_sha256=fix["executed_source_sha256"],final_runner_sha256=fix["final_runner_sha256"],current_core_sha256=sha(R/"python_frontend/compose_eos/whole_crust.py"),evidence_sha256={str(p.relative_to(R)):sha(p) for p in [B/"STELLAR_RESULTS.json",B/"MODE_ACCEPTANCE.json",B/"RUNNER_SETTINGS_CORRECTION.json",B/"FINALIZE_PROCESS.json",B/"RUN.json",D/"H_SUMMARY.json",R/"docs/low_density_crust_20261117/RESULTS_KO.md"]},errors=errors)
(O/"COWLING_REVIEW.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2));raise SystemExit(bool(errors))
