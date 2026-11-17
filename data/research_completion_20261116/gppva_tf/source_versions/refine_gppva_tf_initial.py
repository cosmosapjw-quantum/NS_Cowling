#!/usr/bin/env python3
"""Refine observed competing TF branches; retain numerical holes explicitly."""
from __future__ import annotations
import argparse,gzip,hashlib,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from build_gppva_tf import dump,export,physical_geometry,MODELS
from compose_eos.spatial_tf_equilibrium import solve_screened_cell


def interpolate_initial(old,points):
    p=old["profile"];x=(np.arange(points)+.5)/points
    vals={k:np.interp(x,p["x"],p[k]) for k in ("nn","np","ne")}
    w=np.diff(np.linspace(0,1,points+1)**old["dimension"])
    fac=old["n_b"]/(w@(vals["nn"]+vals["np"]))
    vals["nn"]*=fac;vals["np"]*=fac
    vals["ne"]*=float(w@vals["np"])/max(float(w@vals["ne"]),1e-30)
    return {**vals,"radius":old["radius_fm"]}


def refine(old,points):
    options=dict(points=points,dimension=old["dimension"],seed=old["seed"],
                 radius_bounds=(5.,180.),initial=interpolate_initial(old,points),ftol=1e-13)
    row=solve_screened_cell(old["model_id"],old["n_b"],**options)
    row["geometry"]=physical_geometry(row)
    if row["numerical_status"]=="REQUIRES_REFINEMENT":
        retry=solve_screened_cell(old["model_id"],old["n_b"],
            **{**options,"initial":interpolate_initial(row,points),"ftol":1e-14})
        retry["geometry"]=physical_geometry(retry)
        row["polish_attempt"]=retry
        if retry["numerical_status"]!="REQUIRES_REFINEMENT" and retry["epsilon"]<=row["epsilon"]+1e-11:
            retry["unpolished_result"]={k:v for k,v in row.items() if k not in ("profile","polish_attempt")}
            row=retry
    return row


def good(a):
    return a.get("numerical_status") in ("STATIONARY","UNIFORM_STATIONARY") and not a["muon_threshold_crossed"]


def select(attempts):
    valid=[a for a in attempts if good(a)]
    if not valid:return None,"NO_STATIONARY_CANDIDATE"
    best=min(valid,key=lambda a:a["epsilon"])
    if any("epsilon" in a and not good(a) and a["energy_per_baryon"]<best["energy_per_baryon"]-1e-6 for a in attempts):
        return None,"LOWER_UNCONVERGED_ATTEMPT"
    return best,"LOWEST_EVALUATED_STATIONARY_CANDIDATE"


def run(model,source,output):
    started=time.monotonic();out=Path(output)/model;out.mkdir()
    records=json.loads(gzip.decompress((Path(source)/model/"attempts_N096.json.gz").read_bytes()))
    report=[];selected=[];all_attempts=[]
    for record in records:
        candidates=[a for a in record["attempts"] if "epsilon" in a]
        if not candidates:
            report.append({"n_b":record["n_b"],"status":"NO_INITIAL_BRANCH"});continue
        low=min(a["energy_per_baryon"] for a in candidates)
        # Every observed geometry within 10 keV/baryon of the coarse minimum,
        # plus the uniform reference, is refined; this is a finite search.
        by_geometry={}
        for a in candidates:
            geom=physical_geometry(a)
            if a["energy_per_baryon"]>low+.01 and geom!="uniform":continue
            if geom not in by_geometry or a["epsilon"]<by_geometry[geom]["epsilon"]:by_geometry[geom]=a
        attempts=[]
        for old in by_geometry.values():
            try:attempts.append(refine(old,192))
            except Exception as exc:attempts.append({"execution":"FAILED","geometry":physical_geometry(old),"message":str(exc),"exception":type(exc).__name__})
        best,status=select(attempts);comp=None
        if best:
            same=[a for a in candidates if physical_geometry(a)==best["geometry"]]
            coarse=min(same or candidates,key=lambda a:a["epsilon"])
            comp={"coarse_points":coarse["points"],"fine_points":192,
                  "delta_energy_per_baryon_MeV":best["energy_per_baryon"]-coarse["energy_per_baryon"],
                  "relative_pressure_change":(best["pressure"]-coarse["pressure"])/max(abs(best["pressure"]),1e-30),
                  "delta_Yp":best["Y_p"]-coarse["Y_p"],"geometry":best["geometry"]}
            # Very dilute cells need a finer physical radial spacing.
            if record["n_b"]<=.002:
                refined=[]
                for old in [a for a in attempts if "epsilon" in a and a["energy_per_baryon"]<best["energy_per_baryon"]+.01]:
                    try:refined.append(refine(old,384))
                    except Exception as exc:refined.append({"execution":"FAILED","geometry":physical_geometry(old),"message":str(exc),"exception":type(exc).__name__})
                finer,stat384=select(refined)
                all_attempts.append({"n_b":record["n_b"],"points":384,"attempts":refined})
                if finer:
                    comp={"coarse_points":192,"fine_points":384,
                          "delta_energy_per_baryon_MeV":finer["energy_per_baryon"]-best["energy_per_baryon"],
                          "relative_pressure_change":(finer["pressure"]-best["pressure"])/max(abs(finer["pressure"]),1e-30),
                          "delta_Yp":finer["Y_p"]-best["Y_p"],"geometry":finer["geometry"],"prior_N96_N192":comp}
                    best=finer
                else:best=None;status=stat384
        all_attempts.append({"n_b":record["n_b"],"points":192,"attempts":attempts})
        if best:selected.append(best)
        report.append({"n_b":record["n_b"],"status":status,"selected_geometry":best["geometry"] if best else None,
                       "mesh_comparison":comp,"mesh_pressure_within_1percent":bool(comp and abs(comp["relative_pressure_change"])<.01)})
        dump(out/"progress.json",{"completed":len(report),"total":len(records),"selected":len(selected),"last_density":record["n_b"]})
    dump(out/"attempts.json.gz",all_attempts);dump(out/"selected_states.json.gz",selected)
    exported=export(model,selected,out)
    result={"model_id":model,"seconds":time.monotonic()-started,"rows":len(selected),"requested_rows":len(records),
            "refinements":report,"export":exported,"selection_scope":"N192 all observed geometries within 0.01 MeV/baryon of N96 minimum plus uniform; N384 these competing branches at n<=0.002 fm^-3. Finite search, not proof of global minimum."}
    dump(out/"RECEIPT.json",result)
    return {"model":model,"rows":len(selected),"seconds":result["seconds"],"mesh_1percent":sum(r["mesh_pressure_within_1percent"] for r in report)}


def main():
    p=argparse.ArgumentParser();p.add_argument("source",type=Path);p.add_argument("output",type=Path);p.add_argument("--workers",type=int,default=6);p.add_argument("--models",nargs="+",default=list(MODELS));a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[1]
    paths=["scripts/refine_gppva_tf.py","scripts/build_gppva_tf.py","python_frontend/compose_eos/spatial_tf.py","python_frontend/compose_eos/spatial_tf_equilibrium.py","python_frontend/compose_eos/gppva_rmf.py"]
    identities={s:hashlib.sha256((root/s).read_bytes()).hexdigest() for s in paths}
    dump(a.output/"INPUT.json",{"source":str(a.source),"models":a.models,"source_sha256":identities})
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        jobs={pool.submit(run,m,a.source,a.output):m for m in a.models}
        for job in as_completed(jobs):
            try:r=job.result()
            except Exception as e:r={"model":jobs[job],"status":"FAILED","message":str(e)}
            results.append(r);print(json.dumps(r),flush=True)
    if identities!={s:hashlib.sha256((root/s).read_bytes()).hexdigest() for s in paths}:raise RuntimeError("source changed during refinement")
    dump(a.output/"CAMPAIGN.json",{"models":results,"source_sha256":identities,"source_unchanged_during_execution":True})
if __name__=="__main__":main()
