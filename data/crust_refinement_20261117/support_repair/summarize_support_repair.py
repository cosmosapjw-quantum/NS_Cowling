#!/usr/bin/env python3
"""Bounded postprocessing of G support production; no native solver calls."""
from pathlib import Path
import ast
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
import numpy as np
from numpy.polynomial import Polynomial
from scipy.interpolate import CubicHermiteSpline

ROOT=Path(__file__).resolve().parents[3]
BASE=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from maxwell_three_family_campaign import read_prepared
from compose_eos.maxwell_v2 import _segment_eosdata

def load(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name, value): (BASE/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+"\n")
def rel(p): return str(Path(p).relative_to(ROOT))
def segment(record, comp):
    return _segment_eosdata(read_prepared(record),(comp["phase_id"],comp["start_zero_based"],comp["stop_exclusive"]))

# Reuse exactly the three pure F audit functions, without executing its F-writing top level.
old=ROOT/"data/new_model_maxwell_20261117/matching_validation/audit_success_envelopes.py"
tree=ast.parse(old.read_text())
names={"interval_poly","unit_roots","audit"}
nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
assert len(nodes)==3
scope={"np":np,"Polynomial":Polynomial,"CubicHermiteSpline":CubicHermiteSpline}
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(old),"exec"),scope)
assert np.allclose(scope["unit_roots"](Polynomial.fromroots([.1,.5,.9])),[.1,.5,.9],rtol=0,atol=1e-12)

matching=load(BASE/"MATCHING.json")
regeneration=load(BASE/"REGENERATION.json")
prepared={r["id"]:r for r in load(BASE/"PREPARED_MANIFEST.json")["source_records"]}
fprepared={r["id"]:r for r in load(ROOT/"data/new_model_maxwell_20261117/prepared/MANIFEST.json")["source_records"]}
sources={**fprepared,**prepared}
audits=[];successes=[];rejections=[];cells=[]
for cell in matching["cells"]:
    cells.append({k:cell[k] for k in ["model","anchor","mesh","requested_grid_points","accepted_grid_points","successes"]})
    for index, result in enumerate(cell["candidates"]):
        ident={k:cell[k] for k in ["model","anchor","mesh","core_id"]}
        ident.update(candidate_index=index,source_record=result["source_record"])
        if result["status"]=="GENERATED_VERIFIED_FINITE_DOMAIN":
            low=segment(sources[result["source_record"]],result["crust_component"])
            high=segment(sources[result["core_id"]],result["core_component"])
            a=scope["audit"](low,high,result["transition"])
            a.update(ident);audits.append(a)
            for file,digest in result["output_sha256"].items():
                assert sha(ROOT/result["output_path"]/file)==digest
            assert result["readback"]["numerical"]["status"]=="PASS"
            assert result["readback"]["metadata"]["status"]=="PASS"
            assert result["admission"]["status"]=="ADMITTED"
            assert result["serialized_native_verification"]["overall"]=="PASS"
            successes.append({**ident,**result["transition"],"output_path":result["output_path"],
                "result_path":str(Path(result["output_path"]).parent/"RESULT.json"),
                "density_min":result["density_fm3"][0],"density_max":result["density_fm3"][1],
                "source_crust_component_density":result["crust_component"]["density_fm3"],
                "envelope_status":a["status"]})
        else:
            r={**ident,"reported_category":result["category"],"error":result["error"],
               "interpreted_category":result["category"]}
            if "fewer than two retained source knots" in result["error"]:
                r["interpreted_category"]="INSUFFICIENT_RETAINED_SOURCE_KNOTS"
                low=segment(sources[result["source_record"]],result["crust_component"])
                high=segment(sources[result["core_id"]],result["core_component"])
                ml=(low.eps+low.p)/low.nb;mh=(high.eps+high.p)/high.nb
                sl=CubicHermiteSpline(ml,low.p,low.nb);sh=CubicHermiteSpline(mh,high.p,high.nb)
                lo,hi=max(ml[0],mh[0]),min(ml[-1],mh[-1])
                union=np.unique(np.r_[lo,ml[(ml>lo)&(ml<hi)],mh[(mh>lo)&(mh<hi)],hi])
                roots=[]
                for left,right in zip(union[:-1],union[1:]):
                    poly=scope["interval_poly"](sh,left,right)-scope["interval_poly"](sl,left,right)
                    for x in scope["unit_roots"](poly):
                        mu=float(left+(right-left)*x)
                        if roots and abs(mu-roots[-1]["mu_transition"])<=1e-8: continue
                        tol=1e-10*max(1,abs(mu))
                        roots.append({"mu_transition":mu,"pressure_transition":float(sl(mu)),
                          "n_low":float(sl(mu,1)),"n_high":float(sh(mu,1)),
                          "retained_low_knots":int(np.sum(ml<mu-tol)),
                          "retained_high_knots":int(np.sum(mh>mu+tol)),
                          "low_component_density":low.nb.tolist(),"low_component_mu":ml.tolist()})
                r["independent_root_diagnosis"]=roots
            rejections.append(r)
assert len(audits)==matching["verified_outputs"]
envelope={"schema":"ns_cowling.G_support_envelope_audit.v1",
    "reused_F_function_source":rel(old),"reused_F_function_source_sha256":sha(old),
    "functions":["interval_poly","unit_roots","audit"],"hidden_three_root_control":"PASS",
    "matching_sha256":sha(BASE/"MATCHING.json"),"prepared_sha256":sha(BASE/"PREPARED_MANIFEST.json"),
    "postprocessor_sha256":sha(__file__),"results":audits,
    "status":"PASS" if audits and all(r["status"]=="PASS" for r in audits) else "FAIL"}
save("ENVELOPE_AUDIT.json",envelope)
assert envelope["status"]=="PASS"
with (BASE/"SUCCESS_INDEX.csv").open("w",newline="") as f:
    fields=["model","anchor","mesh","core_id","source_record","mu_transition","pressure_transition",
       "n_low","n_high","epsilon_low","epsilon_high","density_jump","latent_energy_density",
       "maxwell_identity_residual","density_min","density_max","envelope_status","output_path","result_path"]
    writer=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore");writer.writeheader();writer.writerows(successes)

bybranch=defaultdict(dict)
for s in successes: bybranch[(s["model"],s["anchor"])][s["mesh"]]=s
mesh=[]
metrics=["pressure_transition","n_low","n_high","mu_transition","latent_energy_density"]
for (model,anchor),values in bybranch.items():
    for low_n,high_n in [(16,32),(32,64)]:
        if low_n in values and high_n in values:
            lo,hi=values[low_n],values[high_n]
            mesh.append({"model":model,"anchor":anchor,"mesh_pair":[low_n,high_n],
               "relative_change_reference":"absolute difference divided by absolute higher-mesh value",
               "absolute_changes":{k:abs(hi[k]-lo[k]) for k in metrics},
               "relative_changes":{k:abs(hi[k]-lo[k])/abs(hi[k]) for k in metrics}})
states={};failed=[]
for run in regeneration["runs"]:
    key=(run["model"],run["anchor"],run["points"]);states[key]={}
    for p in run["points_results"]:
        for attempt in p["attempts"]:
            assert sha(attempt["state_path"])==attempt["state_sha256"]
        if p["accepted"]:
            state=load(p["state_path"])
            assert state["accepted_stationary"] and state["actual_exit_code"]==0
            states[key][round(p["n_b"],14)]=state
        else:
            last=load(p["attempts"][-1]["state_path"])
            failed.append({"model":run["model"],"anchor":run["anchor"],"mesh":run["points"],"n_b":p["n_b"],
                "last_state_path":p["attempts"][-1]["state_path"],"attempt_count":len(p["attempts"]),
                **{k:last.get(k) for k in ["status","optimizer_mode","optimizer_success","kkt","radius_stationarity","radius_bound","density_bound"]}})
branches=[]
for model in ["GPPVA_NL3wrL55","GPPVA_TM1e"]:
    for n in [16,32,64]:
        upper=states[model,"upper_F_component",n];lower=states[model,"lower_adjacent_branch",n]
        rows=[]
        for density in sorted(set(upper)&set(lower)):
            u,l=upper[density],lower[density]
            rows.append({"n_b":density,"lower_minus_upper_epsilon":l["epsilon"]-u["epsilon"],
               "lower_minus_upper_energy_per_baryon":(l["epsilon"]-u["epsilon"])/density,
               "lower_minus_upper_mu":l["mu_B"]-u["mu_B"],"lower_minus_upper_pressure":l["pressure"]-u["pressure"],
               "lower_radius":l["radius"],"upper_radius":u["radius"]})
        de=np.array([x["lower_minus_upper_energy_per_baryon"] for x in rows])
        branches.append({"model":model,"mesh":n,"common_accepted_grid_nodes":len(rows),
          "energy_per_baryon_difference_min_MeV":float(de.min()),
          "energy_per_baryon_difference_max_MeV":float(de.max()),
          "energy_per_baryon_difference_max_abs_MeV":float(np.max(abs(de))),
          "lower_anchor_lower_energy_nodes_at_1e-9_MeV":int(np.sum(de < -1e-9)),
          "upper_anchor_lower_energy_nodes_at_1e-9_MeV":int(np.sum(de > 1e-9)),
          "near_equal_nodes_at_1e-9_MeV":int(np.sum(abs(de)<=1e-9)),
          "rows":rows})
# Re-read all four original F source identities; no public data replacement is inferred.
original_ids=[]
for r in load(BASE/"DIAGNOSIS.json")["results"]:
    for key in ["crust_source_id","core_source_id"]:
        sid=r[key]
        if sid not in original_ids: read_prepared(fprepared[sid]);original_ids.append(sid)
summary={"schema":"ns_cowling.G_support_results.v1",
    "postprocessor_sha256":sha(__file__),
    "input_sha256":{name:sha(BASE/name) for name in ["DIAGNOSIS.json","PLAN.json","REGENERATION.json","MATCHING.json","PREPARED_MANIFEST.json","TRANSPORT_REPAIR.json","PRODUCTION_RUN.json","MATCHING_RUN.json"]},
    "status":"FINITE_SOURCE_MATCHING_VERIFIED_NUMERICAL_LIMITS_RETAINED",
    "requested_grid_states":regeneration["requested_grid_states"],"accepted_grid_states":regeneration["accepted_grid_states"],
    "native_TF_invocations":regeneration["actual_native_invocations"],"native_TF_repeats_in_matching":0,
    "native_Maxwell_invocations":matching["actual_kernel_invocations"],"verified_outputs":len(successes),
    "pair_types":4,"N64_pair_types_with_success":len({s["model"] for s in successes if s["mesh"]==64}),
    "envelope_audit_status":envelope["status"],"envelope_outputs":len(audits),
    "max_envelope_root_difference_MeV":max(a["root_difference_from_native_MeV"] for a in audits),
    "max_envelope_pressure_relative_residual":max(a["pressure_relative_residual"] for a in audits),
    "original_F_sources_hash_verified_unchanged":original_ids,
    "unmodified_original_public_failures":"NL3 retained-knot and DDME-X C2 failures remain. New native counterparts are separate source identities.",
    "cells":cells,"successes":successes,"mesh_sensitivity":mesh,"branch_comparisons":branches,
    "failed_grid_states":failed,"rejections":rejections,
    "rejection_counts":dict(Counter(r["interpreted_category"] for r in rejections)),
    "limits":["N16/N32/N64 sensitivity is measured, with no new arbitrary convergence tolerance or N128 result.",
      "The 16 outputs are four pair types at different meshes/anchors, not 16 independent physical combinations.",
      "Native finite droplet branches are stationary; no all-geometry global-minimum claim.",
      "Failure gaps are preserved; no interpolation across a failed point or disconnected component.",
      "Public-source supplemental precision was not acquired; DDME-X saturation-density source identity remains a declared convention.",
      "C3/C4 composition checks are inapplicable to thermo-only exports; actual field baryon/charge/beta checks are retained in the states."]}
save("SUMMARY.json",summary)
print(json.dumps({k:summary[k] for k in ["status","requested_grid_states","accepted_grid_states","native_TF_invocations","native_Maxwell_invocations","verified_outputs","pair_types","N64_pair_types_with_success","envelope_audit_status","max_envelope_root_difference_MeV","max_envelope_pressure_relative_residual","rejection_counts"]}))
for m in mesh:
    if m["mesh_pair"]==[32,64]:print("MESH",m["model"],m["anchor"],json.dumps(m["relative_changes"]))
for b in branches:
    print("BRANCH",b["model"],b["mesh"],b["energy_per_baryon_difference_min_MeV"],b["energy_per_baryon_difference_max_MeV"])
