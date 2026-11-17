#!/usr/bin/env python3
"""Summarize actual H N128 support states and immutable G-to-H mesh comparison."""
from pathlib import Path
import sys,json,hashlib,math
ROOT=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117")
BASE=Path(__file__).resolve().parent
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
def rel(p):
    try:return str(Path(p).relative_to(ROOT))
    except ValueError:return str(p)
regen=load(BASE/"REGENERATION.json")
matching=load(BASE/"MATCHING.json")
original=load(ROOT/"data/crust_refinement_20261117/support_repair/MATCHING.json")
sources=[];failed=[];models=[]
for run in regen["runs"]:
    accepted=[]
    for p in run["points_results"]:
        if p["accepted"]:
            st=load(p["state_path"])
            record={k:st[k] for k in ["model","n_b","points","dimension","geometry","epsilon","pressure","mu_B","mu_e","Y_p","radius","kkt","radius_stationarity","optimizer_mode","baryon_error","charge_error","actual_exit_code","density_contrast"]}
            record.update(state_path=p["state_path"],state_sha256=sha(p["state_path"]))
            accepted.append(record);sources.append(record)
        else:
            attempts=[]
            for a in p["attempts"]:
                st=load(a["state_path"])
                attempts.append({**a,**{k:st[k] for k in ["optimizer_mode","kkt","radius_stationarity","baryon_error","charge_error","radius_bound","density_bound"]}})
            failed.append(dict(model=run["model"],n_b=p["n_b"],classification="optimizer_mode8_at_finite_low_physical_residual",attempts=attempts))
    low=[a for a in accepted if a["n_b"]<=.006]
    models.append(dict(model=run["model"],mesh=128,selected_G_anchor=run["seed"]["selected_G_anchor"],
        requested=run["requested_points"],accepted=run["accepted_points"],lowest_accepted_density=min(a["n_b"] for a in accepted),
        geometries=sorted(set(a["geometry"] for a in accepted)),max_kkt=max(a["kkt"] for a in accepted),
        max_radius_stationarity=max(a["radius_stationarity"] for a in accepted),
        max_baryon_error=max(abs(a["baryon_error"]) for a in accepted),max_charge_error=max(abs(a["charge_error"]) for a in accepted),
        lowest_state=low[0]))
successes=[]
for cell in matching["cells"]:
    for r in cell["candidates"]:
        if r["status"]=="GENERATED_VERIFIED_FINITE_DOMAIN":
            sr=dict(model=cell["model"],mesh=cell["mesh"],transition=r["transition"],
                input_support_fm3=r["crust_component"]["density_fm3"],matched_support_fm3=r["density_fm3"],
                output_path=r["output_path"],output_sha256=r["output_sha256"],
                readback_numerical=r["readback"]["numerical"]["status"],readback_metadata=r["readback"]["metadata"]["status"],
                independent_serialized_admission=r["admission"]["status"],serialized_native_verification=r["serialized_native_verification"],
                source_rows_preserved=r["v2_report"]["source_thermodynamic_rows_preserved"],
                physical_energy_shift_MeV=r["v2_report"]["native_report"]["physical_energy_shift_mev"])
            gsel=next(c for c in original["cells"] if c["model"]==cell["model"] and c["anchor"]=="lower_adjacent_branch" and c["mesh"]==64)
            g64=next(q for q in gsel["candidates"] if q["status"]=="GENERATED_VERIFIED_FINITE_DOMAIN")
            sr["G_N64_transition"]=g64["transition"]
            sr["N64_to_N128_relative_change_percent"]={k:100*(r["transition"][k]-g64["transition"][k])/g64["transition"][k] for k in ["pressure_transition","n_low","n_high","latent_energy_density"]}
            g32c=next(c for c in original["cells"] if c["model"]==cell["model"] and c["anchor"]=="lower_adjacent_branch" and c["mesh"]==32)
            g32=next(q for q in g32c["candidates"] if q["status"]=="GENERATED_VERIFIED_FINITE_DOMAIN")
            sr["G_N32_transition"]=g32["transition"]
            successes.append(sr)
control=load(BASE/"bounded_anchor_control/CONTROL_RESULT.json")
st=load(control["state_path"]);control["optimizer_mode"]=st["optimizer_mode"];control["status"]=st["status"]
summary=dict(schema="ns_cowling.H_support_refinement_summary.v1",
  base_G_commit="905b845ad0f14d2f1e7f67f9eb802691144850cd",
  execution_receipt=load(BASE/"PRODUCTION_RECEIPT.json"),
  actual_native_invocations=regen["actual_native_invocations"]+1,initial_native_invocations=regen["actual_native_invocations"],
  bounded_new_invocations=1,initial_maxwell_invocations=matching["actual_kernel_invocations"],
  accepted_requested_states=dict(accepted=regen["accepted_grid_states"],requested=regen["requested_grid_states"]),
  models=models,verified_matches=successes,retained_failed_points=failed,bounded_control=control,
  source_sha256=regen["source_sha256"],
  evidence_sha256={str(p.relative_to(BASE)):sha(p) for p in [BASE/"PLAN.json",BASE/"REGENERATION.json",BASE/"MATCHING.json",BASE/"PREPARED_MANIFEST.json",BASE/"PRODUCTION_RECEIPT.json",BASE/"bounded_anchor_control/PLAN.json",BASE/"bounded_anchor_control/CONTROL_RESULT.json",BASE/"bounded_anchor_control/EXECUTION_RECEIPT.json"]},
  status_by_model={"GPPVA_NL3wrL55":"IMPLEMENTATION_VERIFIED_FINITE_DOMAIN_MATCH_N128","GPPVA_TM1e":"NUMERICALLY_CHECKED_LOW_DENSITY_N128_STATES_MATCH_BLOCKED_BY_OPTIMIZER_SOURCE_GAP"},
  scope=dict(no_global_minimum_claim=True,no_continuum_limit_claim=True,no_stellar_admission_claim=True,no_cross_mesh_splice=True,no_physical_energy_shift=True,
       thermodynamic_composition_scope="Thermodynamics plus phase labels only; native local neutrality is checked in source states. CompOSE composition abundances are unavailable."))
save(BASE/"ACCEPTED_SOURCE_STATES.json",dict(schema="ns_cowling.H_accepted_source_states.v1",states=sources))
save(BASE/"SUMMARY.json",summary)
print(json.dumps(dict(summary_path=str(BASE/"SUMMARY.json"),summary_sha256=sha(BASE/"SUMMARY.json"),successes=successes,model_summary=models),ensure_ascii=False,indent=2))
