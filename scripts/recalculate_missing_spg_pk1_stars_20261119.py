#!/usr/bin/env python3
"""Recalculate 5 SPG(Mi)+PK1 crust/core alternatives and anchored whole stars.

Original native CLD 129-point source remains a calibrated reconstructed
candidate, not the original SPG paper coefficients. A *new* native 4097 PK1
core is used. BSk22 outer/surface are copied EXACTLY from the frozen E
reconstruction, while every new lower join, Maxwell/CSS interface and
1.4M native TOV/Cowling solution is recalculated.

This is explicitly phenomenological, not source-unified microphysics.
"""
from __future__ import annotations
import argparse,hashlib,json,sys,traceback,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos import read_compose,match_crust_core
from compose_eos.api import construct_maxwell_transition
from compose_eos.whole_crust import WholeCrustBarotrope,CausalExponentialInterval
from run_fixed_mass_maxwell_css_20261119 import run_one as modes_1p4

BASE=ROOT/"data/eos_missing_data_recalc_20261119"
SOURCE=ROOT/"data/research_matching_20261117/spg/native_129_repaired"
CORE=BASE/"pk1_core_fresh/table/compose"
OUTER_SOURCE=ROOT/"data/stellar_extension_20261117/eos/fine/REPRESENTATION.json"
OUT=BASE/"spg_pk1_complete_stars"
SOURCE_REP=json.loads(OUTER_SOURCE.read_text())

def sh(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")

def model(name,method):
    first=SOURCE/name/"generated"
    low=read_compose(first); high=read_compose(CORE)
    if method=="maxwell":
        matched,match=construct_maxwell_transition(
            low,high,assume_local_charge_neutrality=True,
            drop_composition=True,return_report=True)
        j=int(match["low_phase_rows"])-1
        pt=float(match["pressure_transition"])
    elif method=="css":
        matched,match=match_crust_core(low,high,name,"PK1_CANONICAL",
                          "css",40,drop_composition=True)
        j=int(match["css_transition_low_index"])
        pt=float(match["css_transition_pressure"])
    else:raise ValueError(method)
    if not match["consistency"]["overall"]["pass_"]:
        raise ValueError("native matching consistency rejected")
    n=np.asarray(matched.nb,dtype=float)
    p=np.asarray(matched.p,dtype=float)
    e=np.asarray(matched.eps,dtype=float)
    assert 0<=j<len(n)-1 and n[j+1]>n[j] and abs(p[j]-p[j+1])<2e-11
    outer=SOURCE_REP["outer"];surface=SOURCE_REP["surface"]
    lower_end=(float(outer["n_fm3"][-1]),float(outer["p_MeV_fm3"][-1]),
               float(outer["n_fm3"][-1]*outer["F_MeV"][-1]))
    failures=[];accepted=None
    for k in range(min(55,j-2)):
        if n[k]<=lower_end[0]:continue
        right=(float(n[k]),float(p[k]),float(e[k]))
        try:
            joint=CausalExponentialInterval(*lower_end,*right)
            diag=joint.validate()
            if diag["status"]!="PASS":raise ValueError("causal lower join rejected")
            contract=dict(surface=surface,outer=outer,
                lower_bridge=joint.contract(),upper=dict(
                    n=n[k:].tolist(),p=p[k:].tolist(),epsilon=e[k:].tolist(),
                    sharp_interfaces=[j-k]))
            eos=WholeCrustBarotrope.from_contract(contract)
            filename=OUT/name/method/"REPRESENTATION.json"
            save(filename,contract)
            # Native fixed-mass wrapper performs actual DVODE TOV and slow
            # and rapid Cowling f/p1 with strict one-sided interface gates.
            star=modes_1p4(name+"__PK1_FRESH4097__"+method,filename)
            if not star.get("status","").startswith("FIXED_1P4_"):
                raise RuntimeError("native 1.4M execution not admitted "+star.get("status","unknown"))
            accepted=dict(closure_source="E_BSk22_SOURCE_INHERITED_NOT_REGENERATED",
                original_outer_source_sha256=sh(OUTER_SOURCE),
                outer_original_217_knots_unchanged=True,
                generator_spg_grid="129 original native reconstructed candidate",
                pk1_4097_source_sha256=sh(CORE/"eos.thermo"),
                original_spg129_source_sha256=sh(first/"eos.thermo"),
                matching=method,transition_pressure_MeV_fm3=pt,
                density_jump=float(n[j+1]-n[j]),join_index=k,
                join_n_fm3=float(n[k]),join_pressure=float(p[k]),
                bridge=diag,source_compose_rows=int(matched.n),
                contract_path=str(filename.relative_to(ROOT)),
                contract_sha256=sh(filename),
                native_star=star,
                no_pretence_original_BSk22_author_ETFSI=True,
                no_global_unified_microphysical_claim=True)
            break
        except Exception as ex:
            failures.append(dict(index=k,type=type(ex).__name__,error=str(ex)[:210]))
    if accepted is None:raise RuntimeError("no lower join and valid 1.4M star: "+str(failures[:5]))
    return accepted

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--models",nargs="+",default=[f"SPG_M{i}" for i in range(1,6)])
    args=ap.parse_args()
    print("SPG_PK1_START",args.models,flush=True)
    report=dict(schema="ns_cowling.missing_eos.spg_pk1_5_complete_stars.v1",
         source="native reconstructed SPG 129-point + freshly regenerated PK1 4097-point",
         fixed_E_outer_source_sha256=sh(OUTER_SOURCE),
         same_source_original_outer=True,
         baseline_source_calibration_unresolved=True,
         models={},status="RUNNING")
    for name in args.models:
        rec={}
        for method in ("maxwell","css"):
            t=time.monotonic()
            try:
                z=model(name,method);rec[method]=dict(status="NATIVE_1P4_EXECUTED",elapsed_seconds=time.monotonic()-t,**z)
                print("SPG_PK1_WHOLE",name,method,"R",z["native_star"]["native_TOV"]["R_km"],
                      "status",z["native_star"]["status"],flush=True)
            except Exception as ex:
                rec[method]=dict(status="BLOCKED",error_type=type(ex).__name__,
                                  error=str(ex),traceback=traceback.format_exc()[-1400:])
                print("SPG_PK1_BLOCKED",name,method,type(ex).__name__,str(ex)[:260],flush=True)
        report["models"][name]=rec
        save(OUT/"PARTIAL.json",report)
    report["status"]="ALL_REQUESTED_ATTEMPTED"
    save(OUT/"RESULT.json",report)
    print("SPG_PK1_EXIT",len(report["models"]),flush=True)

if __name__=="__main__":main()
