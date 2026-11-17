"""Finalize H mode receipts without rerunning preserved failed calculations."""
from pathlib import Path
import json,hashlib,sys,subprocess
ROOT=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117")
OUT=ROOT/"data/low_density_crust_20261117/stellar/ddmex_comparison"
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
run=json.loads((OUT/"RUN.json").read_text())
raw=json.loads((OUT/"STELLAR_RESULTS.json").read_text())
assert raw["status"]!="RUNNING" and raw["source_bytes_unchanged"]
runner=ROOT/"scripts/run_low_density_stellar_20261117.py"
before=runner.read_text()
expected=raw["source_hashes_before"][str(runner)]
assert sha(runner)==expected
snap=OUT/"source_snapshot";snap.mkdir(exist_ok=True)
executed=snap/"run_low_density_stellar_executed.py";executed.write_bytes(runner.read_bytes())
after=before.replace("from dataclasses import asdict","from dataclasses import asdict, replace")
after=after.replace('def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()',
'''def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def bounded_repair_settings(tight_settings, requested_settings):
    """Tighten accuracy while preserving the requested surface location."""
    return replace(tight_settings, surface_h=requested_settings.surface_h)''')
after=after.replace('mode=solve_mode(bg,bracket,cfg,conversion=conversion)\n        original=',
                    'mode=solve_mode(bg,bracket,cfg,conversion=conversion)\n        mode["effective_settings"]=asdict(cfg)\n        original=')
after=after.replace('mode=solve_mode(bg,bracket,tight,conversion=conversion)\n            save(repair,mode);used=repair',
                    'repair_cfg=bounded_repair_settings(tight,cfg)\n            mode=solve_mode(bg,bracket,repair_cfg,conversion=conversion)\n            mode["effective_settings"]=asdict(repair_cfg)\n            save(repair,mode);used=repair')
assert before!=after and 'bounded_repair_settings(tight,cfg)' in after
runner.write_text(after)
sys.path.insert(0,str(ROOT/"python_frontend"))
sys.path.insert(0,str(ROOT/"scripts"))
import run_low_density_stellar_20261117 as fixed
from compose_eos.stellar_cowling import CowlingSettings
from dataclasses import asdict
checks=[]
tight=CowlingSettings(rtol=3e-10,atol=3e-13,max_step_h=.0005,root_xtol_hz=1e-10)
for surface in [1e-7,3e-8,7e-8]:
 requested=CowlingSettings(rtol=3e-9,atol=3e-12,surface_h=surface)
 original=asdict(requested); repaired=fixed.bounded_repair_settings(tight,requested)
 assert repaired.surface_h==surface and repaired.rtol==tight.rtol and asdict(requested)==original
 assert all(asdict(repaired)[k]==v for k,v in asdict(tight).items() if k!="surface_h")
 checks.append(dict(requested_surface_h=surface,effective=asdict(repaired),original_immutable=True,status="PASS"))
ordinary=raw["cases"];refinement=raw["numerical_refinement"]
rows=[]
for c in ordinary+[dict(refinement,id="numerical_refinement")]:
 for conv,modes in c["modes"].items():
  for i,m in enumerate(modes):
   isfine=c["id"]=="numerical_refinement"
   settings=asdict(tight) if m["bounded_tighter_repair"] else (refinement["settings"] if isfine else raw["settings"])
   planned_matches=not(isfine and m["bounded_tighter_repair"] and settings["surface_h"]!=refinement["settings"]["surface_h"])
   original_path=OUT/"calculations"/c["id"]/(conv+f"_ROOT_{i:02d}.json")
   original_mode=json.loads(original_path.read_text())
   rows.append(dict(id=c["id"],conversion=conv,root_index=i,frequency_hz=m["frequency_hz"],
      acceptance=m["acceptance"],matching_residual=m["matching_residual"],relative_join_mismatch=m["relative_join_mismatch"],
      actual_settings=settings,bounded_repair=m["bounded_tighter_repair"],
      planned_refinement_surface_preserved=planned_matches,
      selected_file=m["file"],selected_file_sha256=m["file_sha256"],
      original_file=str(original_path.relative_to(ROOT)),original_file_sha256=sha(original_path),
      original_matching_residual=original_mode["matching_residual"],
      used_for_accepted_frequency_comparison=m["acceptance"]=="PASS" and planned_matches))
recipe=ROOT/"data/low_density_crust_20261117/stellar/reproduction"
recipe.mkdir(exist_ok=True)
for src in ["/home/cosmosapjw/ns_cowling_wt/scratch/prepare_h_stellar_plan.py",
            "/home/cosmosapjw/ns_cowling_wt/scratch/run_h_ddmex_stellar.py",
            "/home/cosmosapjw/ns_cowling_wt/scratch/h_g_baseline_admission.py",
            "/home/cosmosapjw/ns_cowling_wt/scratch/h_initial_admission.py",
            "/home/cosmosapjw/ns_cowling_wt/scratch/h_final_admission.py"]:
 p=Path(src)
 if p.exists():(recipe/p.name).write_bytes(p.read_bytes())
rec=dict(reason="Fresh H review found that a rejected refinement could fall back to the ordinary surface_h; preserve requested surface on future bounded repairs",
 executed_source=str(executed.relative_to(ROOT)),executed_source_sha256=sha(executed),
 final_runner=str(runner.relative_to(ROOT)),final_runner_sha256=sha(runner),
 production_original_receipts_modified=False,configuration_regressions=checks,
 core_physics_changed=False,full_suite_source_scope_unaffected=True,
 refinement_repairs_with_different_surface=[r for r in rows if not r["planned_refinement_surface_preserved"]])
save(OUT/"RUNNER_SETTINGS_CORRECTION.json",rec)
accepted=[r for r in rows if r["used_for_accepted_frequency_comparison"]]
rejected=[r for r in rows if r["acceptance"]!="PASS"]
summary=dict(status="COMPLETED_WITH_REJECTED_CANDIDATES" if rejected else "PASS_SCOPED",
 raw_run_actual_exit_code=run["actual_exit_code"],raw_result_status=raw["status"],
 scope=raw["scope"],all_selected_root_records=rows,
 accepted_selected_records=len(accepted),rejected_selected_records=len(rejected),
 unpaired_refinement_settings_records=len(rec["refinement_repairs_with_different_surface"]),
 complete_spectrum_claim=False,raw_results_sha256=sha(OUT/"STELLAR_RESULTS.json"),
 runner_correction_receipt_sha256=sha(OUT/"RUNNER_SETTINGS_CORRECTION.json"),
 low_root_policy="Rejected candidates retained; do not infer mode absence or use their frequency differences as validated model effects",
 finite_comparisons=[])
base=ordinary[0]
for c in ordinary[1:]:
 item=dict(id=c["id"],delta_R_km=c["background"]["radius_km"]-base["background"]["radius_km"],frequencies={})
 for conv in ["slow","rapid"]:
  item["frequencies"][conv]=[dict(root_index=i,delta_frequency_hz=m["frequency_hz"]-base["modes"][conv][i]["frequency_hz"])
     for i,m in enumerate(c["modes"][conv]) if m["acceptance"]=="PASS" and base["modes"][conv][i]["acceptance"]=="PASS"]
 summary["finite_comparisons"].append(item)
save(OUT/"MODE_ACCEPTANCE.json",summary)
print(json.dumps({k:v for k,v in summary.items() if k!="all_selected_root_records"},indent=2))
