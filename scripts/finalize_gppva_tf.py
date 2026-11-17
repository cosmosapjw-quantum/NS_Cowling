#!/usr/bin/env python3
"""Repair observed numerical holes, retain originals, export and audit TF tables."""
from __future__ import annotations
import argparse,collections,gzip,hashlib,json,shutil
from pathlib import Path
import numpy as np
from build_gppva_tf import dump,export,physical_geometry,MODELS
from refine_gppva_tf import good,select
from tf_polish import polish_profile
from compose_eos.tf_composition import read_tf_composition

def read(path):
    return json.loads(gzip.decompress(path.read_bytes())) if path.name.endswith(".gz") else json.loads(path.read_text())

def comparison(fine,coarse):
    return {"coarse_points":coarse["points"],"fine_points":fine["points"],
        "delta_energy_per_baryon_MeV":fine["energy_per_baryon"]-coarse["energy_per_baryon"],
        "relative_pressure_change":(fine["pressure"]-coarse["pressure"])/max(abs(fine["pressure"]),1e-30),
        "delta_Yp":fine["Y_p"]-coarse["Y_p"],"geometry":physical_geometry(fine)}

def blocks(requested,accepted):
    result=[];current=[]
    for n in requested:
        if n in accepted:current.append(n)
        elif current:result.append([current[0],current[-1]]);current=[]
    if current:result.append([current[0],current[-1]])
    return result

def run(model,source,coarse,outroot):
    out=outroot/model;shutil.copytree(source/model,out)
    before=read(out/"RECEIPT.json");states=read(out/"selected_states.json.gz")
    attempts=read(out/"attempts.json.gz");coarse_rows=read(coarse/model/"attempts_N096.json.gz")
    repairs=[]
    for item in before["refinements"]:
        n=item["n_b"]
        if n<.002 or item["status"]=="LOWEST_EVALUATED_STATIONARY_CANDIDATE":continue
        relevant=[r for r in attempts if r["n_b"]==n];points=max(r["points"] for r in relevant)
        old=[a for r in relevant if r["points"]==points for a in r["attempts"]]
        trials=[a for a in old if good(a)];low=min(a["energy_per_baryon"] for a in old if "epsilon" in a)
        new=[]
        for candidate in old:
            if "profile" not in candidate or good(candidate) or candidate["energy_per_baryon"]>low+.01:continue
            repaired=polish_profile(candidate);repaired["geometry"]=physical_geometry(repaired);new.append(repaired);trials.append(repaired)
        best,status=select(trials);mesh=None
        if best:
            if n<=.002 and best["points"]<384:
                finer=polish_profile(best,384);finer["geometry"]=physical_geometry(finer);new.append(finer)
                if good(finer):mesh=comparison(finer,best);best=finer
                else:best=None;status="N384_REQUIRES_REFINEMENT"
            else:
                base=[a for r in coarse_rows if r["n_b"]==n for a in r["attempts"] if "epsilon" in a and physical_geometry(a)==physical_geometry(best)]
                mesh=comparison(best,min(base,key=lambda a:a["epsilon"]))
        repairs.append({"n_b":n,"prior_status":item["status"],"status":status,"attempts":new})
        if best:
            states=[s for s in states if s["n_b"]!=n]+[best]
            item.update(status=status,selected_geometry=best["geometry"],mesh_comparison=mesh,
                        mesh_pressure_within_1percent=abs(mesh["relative_pressure_change"])<.01,
                        repair="Reoptimized exact functional; same KKT/PDE/bounds acceptance")
    states.sort(key=lambda s:s["n_b"]);dump(out/"repair_attempts.json.gz",repairs);dump(out/"selected_states.json.gz",states)
    shutil.rmtree(out/"compose");exported=export(model,states,out)
    read_tf_composition(out)
    requested=[r["n_b"] for r in before["refinements"]]
    byn={s["n_b"]:s for s in states};bycheck={r["n_b"]:r for r in before["refinements"]}
    accepted={n for n in requested if n in byn and bycheck[n]["mesh_pressure_within_1percent"]}
    midpoint=[]
    for j in range(1,len(requested)-1,2):
        l,n,r=requested[j-1:j+2]
        if not all(v in byn for v in (l,n,r)):continue
        left,mid,right=(byn[v] for v in (l,n,r))
        interp=(left["pressure"]+right["pressure"])/2
        thermo=n*n*(right["energy_per_baryon"]-left["energy_per_baryon"])/(r-l)
        midpoint.append({"n_b":n,"interval":[l,r],"relative_linear_pressure_interpolation_error":(interp-mid["pressure"])/mid["pressure"],
                         "relative_secant_energy_pressure_error":(thermo-mid["pressure"])/mid["pressure"],
                         "geometries":[physical_geometry(v) for v in (left,mid,right)]})
    summary={"model_id":model,"rows":len(states),"requested_rows":len(requested),
        "source_interval_requested_rows":sum(n>=.002 for n in requested),
        "source_interval_selected_rows":sum(s["n_b"]>=.002 for s in states),
        "held_densities":[n for n in requested if n not in byn],
        "sampled_mesh_validated_runs_fm3":blocks(requested,accepted),
        "source_interval_mesh_validated_runs_fm3":blocks([n for n in requested if n>=.002],accepted),
        "geometry_counts":dict(collections.Counter(physical_geometry(s) for s in states)),
        "last_clustered_density_fm3":max((s["n_b"] for s in states if physical_geometry(s)!="uniform"),default=None),
        "first_uniform_density_fm3":min((s["n_b"] for s in states if physical_geometry(s)=="uniform"),default=None),
        "maximum_source_mesh_pressure_change":max(abs(bycheck[n]["mesh_comparison"]["relative_pressure_change"]) for n in accepted if n>=.002),
        "maximum_chemical_KKT_residual_MeV":max(s["max_chemical_KKT_residual_MeV"] for s in states),
        "maximum_baryon_error_fm3":max(s["baryon_constraint_error"] for s in states),
        "maximum_charge_error_fm3":max(s["charge_constraint_error"] for s in states),
        "midpoint_density_diagnostics":midpoint,
        "support_semantics":"Runs contain consecutive requested sampled densities with stationary cells and <1% adjacent-mesh pressure change. They are not a proof over every density, nor a global-minimum certificate. Interpolation diagnostics are independent.",
        "composition_representation":"spatial_TF_species_densities"}
    before.update(rows=len(states),export=exported,repair_attempts="repair_attempts.json.gz",summary=summary)
    dump(out/"RECEIPT.json",before);dump(out/"SUPPORT.json",summary)
    print(json.dumps({k:summary[k] for k in ["model_id","rows","held_densities","source_interval_selected_rows"]}),flush=True)
    return summary

def main():
    p=argparse.ArgumentParser();p.add_argument("source",type=Path);p.add_argument("coarse",type=Path);p.add_argument("output",type=Path);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    root=Path(__file__).resolve().parents[1]
    paths=["scripts/finalize_gppva_tf.py","scripts/tf_polish.py","scripts/refine_gppva_tf.py","scripts/build_gppva_tf.py","python_frontend/compose_eos/spatial_tf.py","python_frontend/compose_eos/spatial_tf_equilibrium.py","python_frontend/compose_eos/gppva_rmf.py","python_frontend/compose_eos/tf_composition.py"]
    ids={s:hashlib.sha256((root/s).read_bytes()).hexdigest() for s in paths}
    dump(a.output/"INPUT.json",{"initial_refinement":str(a.source),"coarse_campaign":str(a.coarse),"source_sha256":ids})
    results=[run(m,a.source,a.coarse,a.output) for m in MODELS]
    assert ids=={s:hashlib.sha256((root/s).read_bytes()).hexdigest() for s in paths}
    dump(a.output/"SUMMARY.json",{"models":results,"source_sha256":ids,"source_unchanged_during_execution":True})
if __name__=="__main__":main()
