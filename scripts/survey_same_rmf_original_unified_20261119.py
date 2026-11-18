#!/usr/bin/env python3
"""Conservative original-crust × same-RMF homogeneous core survey.

Finite native Maxwell/CSS construction from source identity, with measured
cubic/quintic source-interpolant sensitivity. Preserve original source rows.
This is a bounded newly authorized original+generated same-functional pairing,
NOT an attribution of generated homogeneous rows to the public original EOS.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator, BPoly
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.eos_data import EOSData
from compose_eos import read_compose,write_compose,match_crust_core
from compose_eos import gppva_rmf
from compose_eos.native import evaluate_ddrmf,executable
from compose_eos.api import construct_maxwell_transition, _thermo_only_view, diagnose_maxwell_phase
from compose_eos.maxwell_v2 import _segment_eosdata

OUT=ROOT/"data/cowling_maxwell_native_20261119/same_rmf"
MODELS=("DDME2","TM1e","FSU2")
CORPUS=ROOT/"data/all_origin_maxwell_css_20261118/SOURCE_CATALOG.json"

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n")

def slim(e):
    return {"exception":type(e).__name__,"message":str(e)[:600]}

def prepared(d):
    return (float(d.nb[0]),float(d.nb[-1]),int(d.n))

def core_model(name):
    name="GPPVA_"+name
    model=gppva_rmf.model(name)
    mn,mp=model["mn"],model["mp"]
    # Controlled broad grid; low homogeneous states outside physically
    # admissible support are not silently extrapolated into phase coexistence.
    ngrid=np.unique(np.r_[np.linspace(.006,.13,162),
                          np.geomspace(.13,1.2,130)])
    rows=[];invalid=[]
    t=time.perf_counter()
    for n in ngrid:
        try:
            if name=="GPPVA_DDME2":
                r=evaluate_ddrmf(name,float(n))
                if r["status"]!="ok":raise RuntimeError("Fortran DDME2 status")
                row=dict(n=float(n),e=float(r["eps_total_mev_fm3"]),
                         p=float(r["pressure_total_mev_fm3"]),
                         mu=float(r["mu_n_mev"]),mue=float(r["mu_e_mev"]),
                         yp=float(r["yp"]),ye=float(r["ye"]),ymu=float(r["ymu"]),
                         beta_res=float(r["beta_residual_mev"]),
                         charge_res=float(r["charge_residual_fm3"]),
                         euler_res=float(r["eps_total_mev_fm3"])+float(r["pressure_total_mev_fm3"])-n*float(r["mu_n_mev"]))
            else:
                r=gppva_rmf.beta_state(name,float(n))
                row=dict(n=float(n),e=float(r["epsilon"]),p=float(r["pressure"]),
                         mu=float(r["mu_B"]),mue=float(r["mu_e"]),
                         yp=float(r["Y_p"]),ye=float(r["Y_e"]),ymu=float(r["Y_mu"]),
                         beta_res=float(r["beta_residual_mev"]),
                         charge_res=float(r["charge_residual"]),
                         euler_res=float(r["total_euler_residual"]))
            if not all(np.isfinite(z) for z in row.values()):
                raise ValueError("nonfinite RMF state")
            if not row["p"]>0 or not row["mu"]>0:
                raise ValueError("nonpositive homogeneous pressure or mu")
            if abs(row["euler_res"])>2e-6:
                raise ValueError("Euler first-law residual too large")
            if abs(row["beta_res"])>2e-7 or abs(row["charge_res"])>2e-11:
                raise ValueError("RMF equilibrium residual too large")
            rows.append(row)
        except Exception as e:
            invalid.append({"n":float(n),**slim(e)})
    if len(rows)<8:raise RuntimeError("no substantial positive-pressure same-RMF core")
    arr={k:np.array([r[k] for r in rows]) for k in rows[0]}
    n=arr["n"];z=np.zeros(len(n));rest=mn
    muE=(arr["e"]+arr["p"])/n
    assert np.max(np.abs(muE-arr["mu"]))<2e-5
    # The original RMF physical function is unchanged; Q3 encodes mu_n.
    # Phase=0 is the CompOSE homogeneous numeric default; we do NOT call
    # this a provider-original core. The physical branch remains homogeneous.
    d=EOSData(nb=n,mn=mn,mp=mp,Q1=arr["p"]/n,Q2=z.copy(),
        Q3=arr["mu"]/rest-1,Q4=-arr["mue"]/rest,Q5=z.copy(),
        Q6=arr["e"]/(n*rest)-1,Q7=arr["e"]/(n*rest)-1,
        q_extra=z.copy(),Iphase=np.full(len(n),0,dtype=int),
        Yn=z.copy(),Yp=z.copy(),Ye=z.copy(),Ymu=z.copy(),
        Aav=z.copy(),Zav=z.copy(),Yh=z.copy(),beta_equilibrated=True)
    out=OUT/"generated_core"/name
    out.mkdir(parents=True,exist_ok=True)
    write_compose(out,d)
    back=read_compose(out)
    if not np.allclose(back.p,arr["p"],rtol=1e-12,atol=1e-13):
        raise RuntimeError("native core CompOSE pressure readback differs")
    rep=dict(model=name,functional_family=name.removeprefix("GPPVA_"),
       source_registry="python_frontend/compose_eos/gppva_rmf.py",
       evaluator="native_fortran_DD_RMF" if name=="GPPVA_DDME2" else "python_stationary_RMF_oracle",
       parameter_family="exact current source parameter registry; no fit/energy offset",
       original_provider_homogeneous_core=False,
       output_compose_path=str(out.relative_to(ROOT)),
       output_eos_thermo_sha256=sha(out/"eos.thermo"),
       n_points=len(rows),invalid_points=invalid,
       density_range=prepared(back),seconds=time.perf_counter()-t,
       max_abs_beta_residual=max(abs(x["beta_res"]) for x in rows),
       max_abs_euler_residual=max(abs(x["euler_res"]) for x in rows))
    dump(OUT/"generated_core"/(name+".json"),rep)
    print("SAME_RMF_CORE",name,len(rows),"invalid",len(invalid),flush=True)
    return d,rep

def numerical_poly(v,order):
    P=np.array(v.p,float);n=np.array(v.nb,float)
    # Use conserved-baryon Euler mu_E, not provider Q3 which can disagree
    # on original tabulated crusts; the native Maxwell kernel does likewise.
    mu=(np.array(v.eps,float)+P)/n
    if len(mu)<4 or min(np.diff(mu))<=0 or min(np.diff(P))<=0 or min(np.diff(n))<=0:
        raise ValueError("unadmitted polynomial source interval")
    if order=="cubic":
        return CubicHermiteSpline(mu,P,n)
    if order=="quintic":
        slope=PchipInterpolator(mu,n).derivative()(mu)
        if min(slope)<=0:raise ValueError("unphysical quintic susceptibility")
        return BPoly.from_derivatives(mu,[[P[i],n[i],slope[i]] for i in range(len(mu))])
    raise ValueError(order)

def poly_roots(l,h,method):
    plo=numerical_poly(l,method);phi=numerical_poly(h,method)
    mul=(np.asarray(l.eps)+np.asarray(l.p))/np.asarray(l.nb)
    muh=(np.asarray(h.eps)+np.asarray(h.p))/np.asarray(h.nb)
    x0=max(float(mul[0]),float(muh[0]))
    x1=min(float(mul[-1]),float(muh[-1]))
    if x1<=x0:return {"status":"NO_MU_OVERLAP"}
    nodes=np.unique(np.r_[np.asarray(l.mub),np.asarray(h.mub),np.linspace(x0,x1,257)])
    nodes=nodes[(nodes>=x0)&(nodes<=x1)]
    d=phi(nodes)-plo(nodes)
    roots=[]
    for i in range(len(nodes)-1):
        if d[i]*d[i+1]<0:
            mu=float(brentq(lambda v:float(phi(v)-plo(v)),float(nodes[i]),float(nodes[i+1]),xtol=1e-11))
            nl=float(plo(mu,nu=1));nh=float(phi(mu,nu=1))
            if nh>nl>0 and phi(mu)>0:
                roots.append(dict(mu_MeV=mu,P_MeVfm3=float(phi(mu)),n_low=nl,n_high=nh))
    return dict(status="POSITIVE_JUMP_ROOT" if roots else "SCREEN_NO_POSITIVE_ROOT_NOT_PROOF",
                mu_overlap=[x0,x1],roots=roots,checked_knots=int(len(nodes)))

def eligible_components(d):
    # Partner-independent derivative/thermodynamic/phase constraints. CompOSE
    # phase changes never stitched into one polynomial.
    diag=diagnose_maxwell_phase(d)
    from maxwell_three_family_campaign import valid_components
    groups=valid_components(diag["intervals"],d.n)
    return [dict(start=a,stop=b) for a,b in groups if b-a>=4],diag

def survey_pair(label,source_path,core,model,tier,selection=None):
    low=read_compose(source_path)
    if selection:
        a,b=selection
        low=_segment_eosdata(low,(0,a,b))
    low=_thermo_only_view(low)
    high=_thermo_only_view(core)
    try:
        groups,_=eligible_components(low)
        core_groups,_=eligible_components(high)
    except Exception as e:
        return dict(label=label,tier=tier,status="SOURCE_COMPONENT_DIAG_FAILED",error=slim(e))
    row=dict(label=label,tier=tier,model=model,
             original_source_Q3_minus_Euler_mu_max_abs_MeV=float(np.max(np.abs(low.mub-(low.eps+low.p)/low.nb))),
             core_Q3_minus_Euler_mu_max_abs_MeV=float(np.max(np.abs(high.mub-(high.eps+high.p)/high.nb))),
             numerical_interpolation_mu="Euler mu_E=(epsilon+P)/n; original Q3 preserved separately",
             source_path=str(source_path.relative_to(ROOT)),
             source_sha256=sha(source_path/"eos.thermo"),
             low_n=low.n,high_n=high.n,
             components_low=groups,components_high=core_groups,
             matches=[],tier_contract="explicit original or thermodynamically corrected crust + newly computed same-RMF homogeneous core")
    for il,cl in enumerate(groups):
        l=_segment_eosdata(low,(0,cl["start"],cl["stop"]))
        for ih,ch in enumerate(core_groups):
            h=_segment_eosdata(high,(0,ch["start"],ch["stop"]))
            pair=dict(low_component=il,high_component=ih,
                      source_support_low=prepared(l),source_support_high=prepared(h),
                      interpolation={})
            for order in ("cubic","quintic"):
                try:pair["interpolation"][order]=poly_roots(l,h,order)
                except Exception as e:pair["interpolation"][order]=dict(status="REJECTED_INPUT_INTERPOLATION",error=slim(e))
            if not any(v.get("status")=="POSITIVE_JUMP_ROOT" for v in pair["interpolation"].values()):
                row["matches"].append(pair);continue
            try:
                out,rep=construct_maxwell_transition(l,h,
                     assume_local_charge_neutrality=True,
                     drop_composition=True,return_report=True)
                rep_comp=rep["consistency"]["overall"]["pass_"]
                if not rep_comp:raise RuntimeError("native Maxwell consistency FAIL")
                q=dict(status="PASS_NATIVE_FINITE",mu_t=rep["mu_transition"],
                       P_t=rep["pressure_transition"],n_low=rep["n_low"],
                       n_high=rep["n_high"],delta_eps=rep["latent_energy_density"],
                       input_rows=[l.n,h.n])
                pair["maxwell"]=q
            except Exception as e:
                pair["maxwell"]=dict(status="REJECTED_NATIVE",error=slim(e))
            try:
                out,rep=match_crust_core(l,h,label,model,"css",40,
                                        drop_composition=True)
                if not rep["consistency"]["overall"]["pass_"]:
                    raise RuntimeError("native CSS consistency FAIL")
                pair["css"]=dict(status="PASS_NATIVE_FINITE",
                     P_t=rep["css_transition_pressure"],
                     n_low=rep["css_transition_nb_low"],
                     n_high=rep["css_transition_nb_high"],
                     cs2=rep["css_single_sound_speed"])
            except Exception as e:
                pair["css"]=dict(status="REJECTED_NATIVE",error=slim(e))
            row["matches"].append(pair)
    row["dual_pass_count"]=sum(
        p.get("maxwell",{}).get("status")=="PASS_NATIVE_FINITE" and
        p.get("css",{}).get("status")=="PASS_NATIVE_FINITE" for p in row["matches"])
    row["maxwell_pass_count"]=sum(p.get("maxwell",{}).get("status")=="PASS_NATIVE_FINITE" for p in row["matches"])
    row["css_pass_count"]=sum(p.get("css",{}).get("status")=="PASS_NATIVE_FINITE" for p in row["matches"])
    row["status"]="SCOPED_FINITE_SURVEY"
    return row

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    source=json.load(CORPUS.open())["sources"]
    original=[x for x in source if x.get("role")=="crust" and x.get("tier") in ("original","corrected")]
    index={(s["model"],s["tier"]):s for s in original}
    record=dict(schema="ns_cowling.same_rmf_original_crust_survey.v1",
      source_catalog_sha256=sha(CORPUS),model_status={},BBP_BPS_status={
        "status":"BLOCKED_ORIGINAL_STANDALONE_SOURCE_NOT_AVAILABLE",
        "BPS":"outer-only original historical EOS; thermodynamically complete standalone path absent from current sealed source catalog",
        "BBP":"inner-only original historical EOS; thermodynamically complete standalone path absent from current sealed source catalog",
        "no_false_claim":"not evaluated as a BBP+BPS full crust; original embedded BPS/SLy4 alternatives are not BBP"},
      cases=[])
    for name in MODELS:
        model="GPPVA_"+name
        try:
            core,core_record=core_model(name)
            record["model_status"][model]=core_record
        except Exception as e:
            record["model_status"][model]=dict(status="CORE_GENERATION_BLOCKED",error=slim(e));continue
        for tier in ("original","corrected"):
            s=index[(model,tier)]
            src=ROOT/s["compose_path"]
            row=survey_pair(model,src,core,model,tier)
            row["source_catalog_id"]=s["source_id"]
            record["cases"].append(row)
            print("SAME_RMF_PAIR",model,tier,"components",len(row.get("matches",[])),
                  "Maxwell",row.get("maxwell_pass_count"),
                  "CSS",row.get("css_pass_count"),
                  "DUAL",row.get("dual_pass_count"),flush=True)
            dump(OUT/"PARTIAL.json",record)
    # Original unified 1D tables: actually split published table on Iphase.
    for label in ("XMLSLZ_TM1","XMLSLZ_NL3","XMLSLZ_PK1"):
        p=ROOT/f"data/new_model_maxwell_20261117/sources/{label}"
        try:
            original=read_compose(p)
            outer_idx=np.where(original.Iphase!=0)[0]
            inner_idx=np.where(original.Iphase==0)[0]
            if not(len(outer_idx)>4 and len(inner_idx)>4):
                raise ValueError("missing original crust/core phases")
            if int(outer_idx[-1])+1!=int(inner_idx[0]):
                raise ValueError("original crust/core not contiguous in raw file")
            core=_segment_eosdata(original,(0,int(inner_idx[0]),int(inner_idx[-1])+1))
            row=survey_pair(label,p,core,label,"original_unified_split",
                   (int(outer_idx[0]),int(outer_idx[-1])+1))
            row["original_unified"]=dict(total_rows=int(original.n),
               phases={str(x):int(sum(original.Iphase==x)) for x in np.unique(original.Iphase)},
               crust_rows=len(outer_idx),core_rows=len(inner_idx),
               boundary_original=dict(n_low=float(original.nb[outer_idx[-1]]),
                  n_high=float(original.nb[inner_idx[0]]),
                  P_low=float(original.p[outer_idx[-1]]),
                  P_high=float(original.p[inner_idx[0]]),
                  mu_low=float(original.mub[outer_idx[-1]]),
                  mu_high=float(original.mub[inner_idx[0]])),
               original_table_sha256=sha(p/"eos.thermo"))
            record["cases"].append(row)
            print("ORIGINAL_UNIFIED_SPLIT",label,
                  "Maxwell",row.get("maxwell_pass_count"),
                  "CSS",row.get("css_pass_count"),
                  "DUAL",row.get("dual_pass_count"),flush=True)
        except Exception as e:
            record["cases"].append(dict(model=label,status="ORIGINAL_SPLIT_BLOCKED",error=slim(e)))
            print("ORIGINAL_SPLIT_BLOCKED",label,repr(e),flush=True)
        dump(OUT/"PARTIAL.json",record)
    record["status"]="COMPLETED_BOUNDED_SOURCE_RESPECTFUL_SURVEY"
    dump(OUT/"SURVEY.json",record)
    print("SURVEY_DONE",len(record["cases"]),flush=True)

if __name__=="__main__":
    main()
