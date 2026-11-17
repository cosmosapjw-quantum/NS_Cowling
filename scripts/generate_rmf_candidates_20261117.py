#!/usr/bin/env python3
"""Frozen-source native RMF uniform regeneration and bounded scientific controls."""
from __future__ import annotations
import argparse
from pathlib import Path
import json
import os
import platform
import sys
import time
import numpy as np
from scipy.optimize import brentq
from compose_eos.rmf_candidate import (
    PRIMARY_MODELS, DIAGNOSTIC_MODELS, load_parameters, build_driver,
    evaluate_states, export_compose, residual_report, sha256)

def dump(path, obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")

def first_law(exe,p,out):
    ns=np.array([.015,.03,.06,.1,.2,.4,.8,1.])
    results=[]
    for step in [2e-4,1e-4]:
        dens=np.array([[n*(1-step),n,n*(1+step)] for n in ns]).flatten()
        s=evaluate_states(exe,p,dens,work_dir=out/f"firstlaw-{step}")
        rows=[]
        for i,n in enumerate(ns):
            j=3*i;h=n*step
            de=(s["epsilon"][j+2]-s["epsilon"][j])/(2*h)
            dp=(s["pressure"][j+2]-s["pressure"][j])/(2*h)
            dm=(s["mu_B"][j+2]-s["mu_B"][j])/(2*h)
            rows.append({"density":float(n),"dedn_minus_mu_scaled":float((de-s["mu_B"][j+1])/s["mu_B"][j+1]),
                         "dPdn_minus_n_dmudn_scaled":float((dp-n*dm)/max(1.,abs(dp))),
                         "cs2_dP_de":float(dp/de),"cs2_n_dmu_mu":float(n*dm/s["mu_B"][j+1])})
        results.append({"relative_step":step,"rows":rows,
            "max_firstlaw_scaled":max(abs(r["dedn_minus_mu_scaled"]) for r in rows),
            "max_gibbs_duhem_scaled":max(abs(r["dPdn_minus_n_dmudn_scaled"]) for r in rows)})
    return {"status":"NUMERICALLY_CHECKED","two_steps":results,
            "passed":all(r["max_firstlaw_scaled"]<2e-6 and r["max_gibbs_duhem_scaled"]<2e-6 for r in results)}

def saturation(exe,p,out):
    count=0
    def state(n,yp=.5):
        nonlocal count
        count+=1
        return evaluate_states(exe,p,[n],mode="baryons",yp=yp,work_dir=out/f"eval-{count:03d}")
    def pressure(n):return float(state(n)["pressure"][0])
    nsat=brentq(pressure,.12,.18,xtol=2e-12)
    s=state(nsat)
    h=nsat*1e-4
    dp=(pressure(nsat+h)-pressure(nsat-h))/(2*h)
    def symmetry(n):
        delta=.001
        v=[]
        for d in [-2*delta,-delta,delta,2*delta]:
            st=state(n,(1-d)/2);v.append(float(st["mu_n"][0]-st["mu_p"][0]))
        return (v[0]-8*v[1]+8*v[2]-v[3])/(48*delta)
    sym=symmetry(nsat)
    slope=3*nsat*(symmetry(nsat+h)-symmetry(nsat-h))/(2*h)
    result={"n_sat_fm3":nsat,"binding_mev":float(s["epsilon"][0]/nsat-(p["mn"]+p["mp"])/2),
        "K_mev":9*dp,"J_mev":sym,"L_mev":slope,
        "effective_mass_ratio_n":float(s["effective_mass_n"][0]/p["mn"]),
        "effective_mass_ratio_p":float(s["effective_mass_p"][0]/p["mp"]),
        "pressure_at_root":float(s["pressure"][0]),
        "source_controls":{k:v for k,v in p.items() if k.startswith("nmp_controls")},
        "status":"NUMERICALLY_CHECKED_SOURCE_CONTROLS_REPORTED_SEPARATELY"}
    return result

def provider_source(directory, densities):
    directory=Path(directory)
    nb=np.loadtxt(directory/"eos.nb",skiprows=2)
    lines=(directory/"eos.thermo").read_text().splitlines()
    header=np.array([float(x) for x in lines[0].split()])
    thermo=np.array([[float(x) for x in ln.split()] for ln in lines[1:] if ln.strip()])
    compo=[ln.split() for ln in (directory/"eos.compo").read_text().splitlines() if ln.strip()]
    output={"n_b":[],"pressure":[],"epsilon":[],"mu_B":[],"mu_Q":[],"Y_p":[],"Y_e":[],"Y_mu":[]}
    indices=[]
    for n in densities:
        i=int(np.argmin(abs(nb-n)))
        if abs(nb[i]-n)>2e-12:raise RuntimeError("comparison requires exact common source density")
        row=compo[i];phase=int(row[3])
        if phase!=0:raise RuntimeError("comparison row is not raw uniform phase 0")
        fractions={int(row[5+2*k]):float(row[6+2*k]) for k in range(int(row[4]))}
        q=thermo[i,3:10];pressure=nb[i]*q[0];eps=nb[i]*header[0]*(1+q[6])
        values=[nb[i],pressure,eps,(eps+pressure)/nb[i],header[0]*(1+q[2]),
                fractions.get(11,0),fractions.get(0,0),fractions.get(1,0)]
        for key,value in zip(output,values):output[key].append(float(value))
        indices.append(i)
    output.update(header_mn=float(header[0]),header_mp=float(header[1]),source_row_indices_zero_based=indices,
       raw_sha256={f:sha256(directory/f) for f in ("eos.nb","eos.thermo","eos.compo")})
    return output

def compare_providers(exe,params,source_root,out):
    densities=[.1,.2,.4,.8]
    providers={name:provider_source(source_root/name,densities) for name in
       ["XMLSLZ_NL3","XMLSLZ_PK1","XMLSLZ_TM1","XMLSLZ_DDME-X"]}
    models={name:evaluate_states(exe,p,densities,work_dir=out/"comparison_native"/name) for name,p in params.items()}
    matrix=[]
    for model,s in models.items():
        for provider,r in providers.items():
            fields={}
            for key in ("pressure","epsilon","mu_B","Y_p","Y_e","Y_mu"):
                delta=s[key]-np.array(r[key])
                fields[key]={"signed_difference":delta.tolist(),"max_abs":float(max(abs(delta))),
                    "max_scaled":float(max(abs(delta)/np.maximum(1.,abs(np.array(r[key])))))}
            matrix.append({"model":model,"provider":provider,"residuals":fields,
                "sum_thermo_scaled":sum(fields[k]["max_scaled"] for k in ("pressure","epsilon","mu_B")),
                "native_mu_n_minus_provider_mu_Q":(s["mu_B"]-r["mu_Q"]).tolist()})
    return {"densities":densities,"raw_provider_data":providers,"matrix":matrix,
        "status":"NUMERICAL_COMPARISON_NO_PARAMETER_FIT",
        "reference_policy":"Physical epsilon reconstructed with original provider header; no physical energy rebasing.",
        "closest_thermo_provider":{model:min((x for x in matrix if x["model"]==model),key=lambda x:x["sum_thermo_scaled"])["provider"] for model in models}}

def dd2_control(exe,out):
    from compose_eos.gppva_rmf import model
    from compose_eos.native import evaluate_ddrmf
    p=model("DD2");p.update(id="DD2_CONTROL",density_dependent=True)
    ns=[.1,.4]
    s=evaluate_states(exe,p,ns,work_dir=out/"DD2_native")
    old=[evaluate_ddrmf("GPPVA_DD2",n) for n in ns]
    mapping={"epsilon":"eps_total_mev_fm3","pressure":"pressure_total_mev_fm3",
             "mu_n":"mu_n_mev","mu_p":"mu_p_mev","mu_e":"mu_e_mev","Y_p":"yp","Y_e":"ye","Y_mu":"ymu",
             "rearrangement":"sigma_rearrangement_mev"}
    delta={key:(s[key]-np.array([r[other] for r in old])).tolist() for key,other in mapping.items()}
    scaled={key:float(max(abs(np.array(v))/np.maximum(1,abs(s[key])))) for key,v in delta.items()}
    return {"densities":ns,"independent_backend":"compose_eos_rmf.f90:evaluate_ddrmf",
            "existing_native_results":old,"differences":delta,"max_scaled":scaled,
            "passed":max(scaled.values())<2e-8}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--repo",required=True,type=Path)
    parser.add_argument("--source-root",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    validation=out/"uniform_validation";validation.mkdir(exist_ok=True)
    start=time.time()
    source=args.source_root/"PARAMETER_INPUTS.json"
    params={model:load_parameters(source,model) for model in PRIMARY_MODELS+DIAGNOSTIC_MODELS}
    dump(validation/"CHECKPOINT.json",{"status":"PRE_PRODUCTION_CHECKPOINT","models":list(params),
      "source_parameter_sha256":sha256(source),"physical_energy_adjustment":"NONE",
      "command":sys.argv,"python":sys.version,"platform":platform.platform(),
      "stage_scope":"C5-C8 uniform regeneration; E publication prerequisite confirmed by parent"})
    exe=build_driver(args.repo,validation/"build")
    generation=[]
    for model,p in params.items():
        for level,m in [("coarse",1),("fine",2)]:
            n=np.r_[np.linspace(.015,.16,256*m+1),np.geomspace(.16,1.,128*m+1)[1:]]
            folder=out/"generated"/model/level
            s=evaluate_states(exe,p,n,work_dir=folder)
            receipt=export_compose(folder,s,p)
            generation.append({"model":model,"level":level,"receipt":str(folder/"GENERATION.json"),
                               "compose_directory":receipt["compose_directory"]})
        dump(out/"UNIFORM_MANIFEST.json",{"status":"TABLES_GENERATED_VALIDATION_PENDING","tables":generation,
                    "executable":str(exe)})
        print(f"TABLES_READY {model}",flush=True)
    validation_results={}
    for model,p in params.items():
        mdir=validation/model;mdir.mkdir(exist_ok=True)
        result={"first_law":first_law(exe,p,mdir),"saturation":saturation(exe,p,mdir/"saturation")}
        dump(mdir/"VALIDATION.json",result);validation_results[model]=result
        print(f"VALIDATION_READY {model}",flush=True)
    controls=dd2_control(exe,validation)
    dump(validation/"DD2_CONTROL.json",controls)
    comparison=compare_providers(exe,params,args.source_root,validation)
    dump(validation/"PROVIDER_COMPARISON.json",comparison)
    passed=all(v["first_law"]["passed"] for v in validation_results.values()) and controls["passed"]
    dump(validation/"VALIDATION_SUMMARY.json",{
        "status":"IMPLEMENTATION_VERIFIED" if passed else "VALIDATION_FAILED",
        "two_step_firstlaw_passed":all(v["first_law"]["passed"] for v in validation_results.values()),
        "independent_DD2_passed":controls["passed"],"provider_identity":"NUMERICAL_RESIDUALS_REPORTED",
        "models":validation_results,"elapsed_seconds":time.time()-start})
    dump(out/"UNIFORM_MANIFEST.json",{"status":"IMPLEMENTATION_VERIFIED" if passed else "VALIDATION_FAILED",
          "tables":generation,"executable":str(exe),"validation_summary":str(validation/"VALIDATION_SUMMARY.json"),
          "provider_comparison":str(validation/"PROVIDER_COMPARISON.json")})
    if not passed:raise RuntimeError("scientific controls failed")
    print("UNIFORM_GENERATION_COMPLETE",flush=True)
if __name__=="__main__":main()
