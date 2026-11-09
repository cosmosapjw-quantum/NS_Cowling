#!/usr/bin/env python3
"""Recompute v5 Cowling f/p1 spectra over the initial rising-mass TOV branch."""
from __future__ import annotations
import argparse, csv, hashlib, json, math, sys
from collections import Counter
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_cowling.solver import solve_v5_tov_star,solve_v5_cowling_spectrum
from compose_cowling.mass_sequence import (extract_candidate_mass_segments, extract_stable_branch,
    find_mass_bracket,paired_sensitivity)
MASS_GRID=(.8,1.,1.2,1.33,1.4,1.6,1.8,2.)
METHODS=("css","cubic_F","cubic_P","maxwell","quintic_F","quintic_P")
REF=ROOT/"oscillation_release_v1/data/spectrum_summary_234.csv"
EOS_DIR=ROOT/"oscillation_release_v1/data/eos_projected"
RELEASE=ROOT/"oscillation_release_v1/provenance/release_summary.json"
V5_COMPILE_MANIFEST=ROOT/"oscillation_release_v1/workflow/compile_manifest.json"


def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()


def read_csv(path):
    with path.open(newline="",encoding="utf-8") as f: return list(csv.DictReader(f))


def write_csv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore");w.writeheader()
        for r in rows:
            def encode(v):
                if isinstance(v,(dict,list,tuple)):
                    return json.dumps(v,sort_keys=True)
                if isinstance(v,bool):
                    return str(v).lower()
                return v
            w.writerow({k:encode(v) for k,v in r.items()})


def pressure_bounds(eos):
    p=[float(r["pressure_pa"]) for r in read_csv(eos)]
    if len(p)<3 or any(not math.isfinite(x) or x<=0 for x in p) or any(b<=a for a,b in zip(p,p[1:])):
        raise ValueError(f"invalid normalized EOS pressure domain: {eos}")
    return p[0],p[-1]


def node_identity(result):
    try:
        n=result["node_counts"]
        f=[n["f"][r][v] for r in ("inner","outer") for v in ("W","V")]
        p=[n["p1"][r][v] for r in ("inner","outer") for v in ("W","V")]
    except (KeyError,TypeError):
        return False,"missing_node_diagnostics"
    if any(x is None for x in f+p): return False,"incomplete_node_diagnostics"
    if f != [0,0,0,0]: return False,f"f_node_identity_mismatch:{f}"
    if p != [1,1,0,0]: return False,f"p1_node_identity_mismatch:{p}"
    return True,"smooth_l2_WV_nodes_v1"


def solve_model(ref,eos_dir,masses,solar_mass,*,nscan,tov_rtol,spectrum_rtol,mode_rtol,mass_rtol,residual_tol,radius_limit_m):
    name=ref["name"];eos=eos_dir/f"{name}.csv";pcut,pmax=pressure_bounds(eos);pc_ref=float(ref["pc_pa"])
    grid=np.geomspace(max(1.05*pcut,5e-4*pc_ref),min(.95*pmax,100*pc_ref),nscan)
    scan=[]
    for i,pc in enumerate(grid):
        r={"scan_index":i,"pc_pa":float(pc),"status":"failed"}
        try:
            star=solve_v5_tov_star(eos,float(pc),p_cut_pa=pcut,rtol=tov_rtol,step_cap=.05,
                r0_m=10.0,r_limit_m=radius_limit_m)
            r.update(status="ok",mass_kg=float(star["mass_kg"]),radius_m=float(star["radius_m"]),
                     surface_residual_pa=star.get("surface_residual_pa"))
        except Exception as exc: r["failure"]=f"{type(exc).__name__}: {exc}"[:1000]
        scan.append(r)
    branch=extract_stable_branch(scan)
    candidate_segments=extract_candidate_mass_segments(scan)
    candidate_by_scan={}
    for segment in candidate_segments:
        for scan_index in range(segment["start_scan_index"],segment["end_scan_index"]+1):
            candidate_by_scan.setdefault(scan_index,[]).append(segment)
    initial_rise_ids={int(p["scan_index"]) for p in branch.points}
    tov_rows=[]
    for pnt in scan:
        segments=candidate_by_scan.get(int(pnt["scan_index"]),[])
        tov_rows.append({**pnt,"name":name,"crust":ref["crust"],"core":ref["core"],
            "method":ref["matching_method"],
            "initial_rising_branch_proxy_member":int(pnt["scan_index"]) in initial_rise_ids,
            "candidate_segment_ids":";".join(segment["candidate_segment_id"] for segment in segments),
            "candidate_segment_directions":";".join(segment["direction"] for segment in segments),
            "radial_stability_status":"unresolved_requires_radial_mode_spectrum" if segments else "unresolved_invalid_tov_sample"})
    candidate_rows=[{**segment,"name":name,"crust":ref["crust"],"core":ref["core"],
        "method":ref["matching_method"],"kinematic_only":True} for segment in candidate_segments]
    baseline=solve_v5_tov_star(eos,pc_ref,p_cut_pa=pcut,rtol=tov_rtol,step_cap=.05,
        r0_m=10.0,r_limit_m=radius_limit_m)
    bm=float(baseline["mass_kg"])/float(ref["mass_kg"])-1
    br=float(baseline["radius_m"])/float(ref["radius_m"])-1
    rows=[]
    for m in masses:
        target=m*solar_mass;bracket=find_mass_bracket(branch,target)
        base={"name":name,"crust":ref["crust"],"core":ref["core"],"method":ref["matching_method"],
              "mass_msun":m,"turning_point_found":branch.turning_point_found,"branch_status":branch.status}
        if bracket is None:
            rows.append({**base,"valid":False,"status":"target_outside_initial_rising_branch_proxy",
                "first_turning_point_mass_msun":branch.maximum_mass_kg/solar_mass if branch.maximum_mass_kg else None,
                "branch_endpoint_mass_msun":branch.branch_endpoint_mass_kg/solar_mass if branch.branch_endpoint_mass_kg else None})
            continue
        scale=math.sqrt((target/bracket["radius_guess_m"]**3)/(float(ref["mass_kg"])/float(ref["radius_m"])**3))
        f0=float(ref["f_hz"])*scale;p10=float(ref["p1_hz"])*scale
        flo=[.72*f0,.72*p10];fhi=[1.28*f0,1.28*p10]
        try:
            out=solve_v5_cowling_spectrum(eos,bracket["pc_seed_pa"],flo,fhi,target_mass_kg=target,p_cut_pa=pcut,
                rtol=spectrum_rtol,mode_rtol=mode_rtol,mass_rtol=mass_rtol,residual_tol=residual_tol,
                pc_bracket_pa=(bracket["pc_lo_pa"],bracket["pc_hi_pa"]),r_limit_m=radius_limit_m)
            identity_ok,identity=node_identity(out);mrel=float(out["mass_kg"])/target-1
            valid=(identity_ok and abs(mrel)<=5*mass_rtol and abs(float(out["f_residual"]))<=residual_tol
                   and abs(float(out["p1_residual"]))<=residual_tol)
            rows.append({**base,"valid":bool(valid),"status":"ok" if valid else identity,
                "target_mass_kg":target,"mass_kg":float(out["mass_kg"]),"mass_relative_residual":mrel,
                "radius_m":float(out["radius_m"]),"pc_pa":float(out["pc_pa"]),
                "f_hz":float(out["f_hz"]),"p1_hz":float(out["p1_hz"]),
                "f_residual":float(out["f_residual"]),"p1_residual":float(out["p1_residual"]),
                "f_bracket_hz":float(out["f_bracket_hz"]),"p1_bracket_hz":float(out["p1_bracket_hz"]),
                "frequency_scale_mean_density":scale,"node_counts_json":out.get("node_counts"),
                "mode_identity_id":identity,"pc_lo_pa":bracket["pc_lo_pa"],"pc_hi_pa":bracket["pc_hi_pa"]})
        except Exception as exc:
            rows.append({**base,"valid":False,"status":"solver_failure","failure":f"{type(exc).__name__}: {exc}"[:2000]})
    model={"name":name,"crust":ref["crust"],"core":ref["core"],"method":ref["matching_method"],
        "projected_eos_path":str(eos.relative_to(ROOT)),"projected_eos_sha256":sha256(eos),
        "reference_pc_pa":pc_ref,"baseline_mass_kg":float(baseline["mass_kg"]),
        "baseline_radius_m":float(baseline["radius_m"]),
        "baseline_reproduction_mass_relative":bm,"baseline_reproduction_radius_relative":br,
        "tov_backend":"v5 release compile_manifest source family; calculation/methods/support cell_solve task=star",
        "v5_backend_build_key":baseline.get("backend_build_key"),
        "v5_compile_manifest_sha256":sha256(V5_COMPILE_MANIFEST),
        "v5_reference_source_zip_sha256":json.loads(RELEASE.read_text()).get("source_zip_sha256"),
        "v5_source_byte_identity_vs_reference_archive":"NOT_VERIFIED_IN_THIS_SESSION",
        "r0_m":10.0,"r_limit_m":radius_limit_m,"tov_step_cap":0.05,"tov_rtol":tov_rtol,
        "baseline_reproduction_status":"pass" if abs(bm)<=2e-5 and abs(br)<=2e-5 else "mismatch",
        "branch_status":branch.status,"turning_point_found":branch.turning_point_found,
        "branch_stability_classification":"unresolved_requires_radial_mode_spectrum",
        "candidate_mass_segment_count":len(candidate_segments),
        "candidate_mass_segment_directions":[row["direction"] for row in candidate_segments],
        "initial_rising_branch_proxy_points":len(branch.points),
        "first_turning_point_mass_kg":branch.maximum_mass_kg,
        "first_turning_point_pc_pa":branch.maximum_mass_pc_pa,"branch_endpoint_mass_kg":branch.branch_endpoint_mass_kg,
        "first_branch_pc_pa":float(branch.points[0]["pc_pa"]) if branch.points else None,
        "last_branch_pc_pa":float(branch.points[-1]["pc_pa"]) if branch.points else None,
        "low_pressure_truncated":branch.low_pressure_truncated,
        "target_mass_valid_count":sum(bool(r["valid"]) for r in rows)}
    return tov_rows,rows,model,candidate_rows


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",type=Path,default=ROOT/"dist/mass-dependent-cowling")
    p.add_argument("--mass-grid",default=",".join(map(str,MASS_GRID)))
    p.add_argument("--crusts",default="");p.add_argument("--cores",default="");p.add_argument("--limit-models",type=int,default=0)
    p.add_argument("--tov-samples",type=int,default=36);p.add_argument("--tov-rtol",type=float,default=1e-9)
    p.add_argument("--spectrum-rtol",type=float,default=1e-9);p.add_argument("--mode-rtol",type=float,default=1e-9)
    p.add_argument("--mass-rtol",type=float,default=1e-10);p.add_argument("--residual-tol",type=float,default=1e-6)
    p.add_argument("--radius-limit-m",type=float,default=1e8)
    a=p.parse_args(argv);masses=tuple(float(x) for x in a.mass_grid.split(",") if x)
    if not masses or tuple(sorted(set(masses)))!=masses or any(x<=0 or not math.isfinite(x) for x in masses):
        raise SystemExit("mass grid must be strictly increasing finite positive solar masses")
    if a.tov_samples<8: raise SystemExit("--tov-samples must be >= 8")
    if not math.isfinite(a.radius_limit_m) or a.radius_limit_m<=10.0:
        raise SystemExit("--radius-limit-m must be finite and > 10 m")
    refs=read_csv(REF)
    if len(refs)!=234 or len({r["name"] for r in refs})!=234: raise SystemExit("frozen 234-row cohort is not intact")
    if {r["matching_method"] for r in refs}!=set(METHODS): raise SystemExit("six-method cohort contract mismatch")
    if Counter(r["core"] for r in refs)!=Counter({"BBB_BHF2":78,"GM_GM1":78,"RG_SLY4":78}):
        raise SystemExit("frozen core-family cohort balance changed")
    if Counter(r["crust"] for r in refs)!=Counter({name:18 for name in sorted({r["crust"] for r in refs})}):
        raise SystemExit("frozen crust-family cohort balance changed")
    pair_methods={}
    for row in refs:
        pair_methods.setdefault((row["crust"],row["core"]),[]).append(row["matching_method"])
    if len(pair_methods)!=39 or any(len(v)!=6 or set(v)!=set(METHODS) for v in pair_methods.values()):
        raise SystemExit("frozen core/crust paired six-method design is not balanced")
    cs=set(x for x in a.crusts.split(",") if x);ks=set(x for x in a.cores.split(",") if x)
    selected=[r for r in refs if (not cs or r["crust"] in cs) and (not ks or r["core"] in ks)]
    if a.limit_models: selected=selected[:a.limit_models]
    if not selected: raise SystemExit("no EOS rows selected")
    if any(not (EOS_DIR/f"{r['name']}.csv").is_file() for r in selected): raise SystemExit("projected EOS table missing")
    solar=float(json.loads(RELEASE.read_text())["solar_mass_nominal_kg"])
    a.output_dir.mkdir(parents=True,exist_ok=True);tov_all=[];modes=[];models=[];candidate_all=[];failures=[]
    for i,ref in enumerate(selected,1):
        print(f"[{i}/{len(selected)}] {ref['name']}",flush=True)
        try:
            tr,sr,model,candidates=solve_model(ref,EOS_DIR,masses,solar,nscan=a.tov_samples,tov_rtol=a.tov_rtol,
                spectrum_rtol=a.spectrum_rtol,mode_rtol=a.mode_rtol,mass_rtol=a.mass_rtol,residual_tol=a.residual_tol,radius_limit_m=a.radius_limit_m)
            tov_all.extend(tr);modes.extend(sr);models.append(model);candidate_all.extend(candidates)
            write_csv(a.output_dir/"tov_sequence_partial.csv",tov_all);write_csv(a.output_dir/"spectrum_partial.csv",modes)
            cap=(f"first-turning-point candidate={model['first_turning_point_mass_kg']/solar:.4f} Msun" if model["turning_point_found"]
                 else f"rising-proxy-end={model['branch_endpoint_mass_kg']/solar:.4f} Msun")
            print(f"  {model['branch_status']}; {model['target_mass_valid_count']}/{len(masses)} spectra; {cap}",flush=True)
        except Exception as exc:
            failures.append({"name":ref["name"],"error":f"{type(exc).__name__}: {exc}"})
            print(f"  MODEL FAILURE: {failures[-1]['error']}",flush=True);break
    paired,deltas=paired_sensitivity(modes,METHODS)
    for filename,rows in (("tov_sequence.csv",tov_all),("candidate_mass_segments.csv",candidate_all),("spectrum_mass_grid.csv",modes),
        ("paired_mass_sensitivity.csv",paired),("paired_method_deviations.csv",deltas),("model_branch_summary.csv",models)):
        write_csv(a.output_dir/filename,rows)
    manifest={"schema":"ns-cowling-mass-dependent-spectrum-v1","reference_cohort_sha256":sha256(REF),
        "release_summary_sha256":sha256(RELEASE),
        "v5_compile_manifest_sha256":sha256(V5_COMPILE_MANIFEST),
        "v5_reference_source_zip_sha256":json.loads(RELEASE.read_text()).get("source_zip_sha256"),
        "v5_source_byte_identity_vs_reference_archive":"NOT_VERIFIED_IN_THIS_SESSION",
        "v5_backend_build_keys":sorted({m.get("v5_backend_build_key") for m in models if m.get("v5_backend_build_key")}),"selected_model_count":len(selected),"completed_model_count":len(models),
        "failed_model_count":len(failures),"mass_grid_msun":list(masses),"expected_methods":list(METHODS),
        "expected_spectrum_count":len(selected)*len(masses),"valid_spectrum_count":sum(bool(r.get("valid")) for r in modes),
        "complete_paired_mass_cells":sum(bool(r["pair_complete"]) for r in paired),
        "baseline_reproduction_failures":[m["name"] for m in models if m["baseline_reproduction_status"]!="pass"],
        "models":models,"model_failures":failures,
        "scope":{"cohort":"finite 13-crust x 3-core x 6-smooth-method crossed design; not 234 independent EOS samples",
        "maxwell_label":"historical smooth bridge, not a sharp first-order Maxwell transition",
        "stability":"the initial rising-M(Pc) branch is only a turning-point proxy; every monotonic candidate segment remains radially unresolved until a radial pulsation spectrum/mode count is computed",
        "mode_identity":"ell=2 smooth Cowling W/V nodes: f=(0,0)/(0,0), p1=(1,1)/(0,0), inner/outer",
        "frequency_seed":"0.72-1.28 mean-density scaling about frozen 1.33-Msun row, then node gate",
        "background_solver":"frozen v5 compile_manifest cell_solve; center_order=2 for TOV and Cowling",
        "tov_r0_m":10.0,"tov_step_cap":0.05,"tov_r_limit_m":a.radius_limit_m,"center_order":2,
        "candidate_segments":"kinematic M(Pc) monotonic pieces only; not radial-stability labels",
        "surface":"historical finite positive pressure cutoff; no zero-pressure extrapolation",
        "gravity":"Cowling frequencies, not full-GR complex frequencies or damping times"}}
    (a.output_dir/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:manifest[k] for k in ("selected_model_count","completed_model_count","valid_spectrum_count",
        "expected_spectrum_count","complete_paired_mass_cells","baseline_reproduction_failures","model_failures")},indent=2),flush=True)
    return 0 if not failures and len(models)==len(selected) else 2

if __name__=="__main__": raise SystemExit(main())
