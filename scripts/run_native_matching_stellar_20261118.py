"""Finite 1.4 M_sun TOV/Cowling on actual native Hermite/CSS/Maxwell constructions.

All methods share exactly the same source crust/core pair, H/E outer crust and
phenomenological lower connector, and stellar closure. CSS differs physically
from Maxwell; it is not an interpolation approximation to Maxwell coexistence.
"""
from __future__ import annotations
from pathlib import Path
import sys,json,hashlib,time,argparse,traceback,copy
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from native_hermite_matching_matrix_20261118 import source
from compose_eos import read_compose,write_compose
from compose_eos.api import construct_maxwell_transition,_thermo_only_view
from compose_eos.whole_crust import WholeCrustBarotrope,CausalExponentialInterval
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background,matching_residual,solve_mode
from compose_eos.stellar_cowling import CowlingSettings
OUT=ROOT/"data/cubic_quintic_matching_20261118/native_stellar"
NATIVE=ROOT/"data/cubic_quintic_matching_20261118/native_matching"
SOURCE_J={"SPG_M1":15,"SPG_M3":13}
REFERENCE=ROOT/"data/stellar_extension_20261117/eos/fine/REPRESENTATION.json"
TOV=dict(r0_km=.001,rtol=1e-9,atol=1e-11,max_step_h=.005)
COARSE=CowlingSettings(ell=2,rtol=2e-7,atol=2e-10,max_step_h=.004,surface_h=3e-8,root_xtol_hz=1e-6)
FINE=CowlingSettings(ell=2,rtol=3e-9,atol=3e-12,max_step_h=.001,surface_h=3e-8,root_xtol_hz=1e-8)
WINDOWS={"f":(1550.,2900.),"p1":(4600.,7800.)}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def whole_from_match(crust,matched,transition_index,method):
    j=SOURCE_J[crust]
    if j>=len(matched.nb)-3:raise ValueError("source row start outside matching support")
    source_ref=json.loads(REFERENCE.read_text())
    base=WholeCrustBarotrope.from_contract(source_ref)
    left=(float(base.outer.nb[-1]),float(base.outer.p[-1]),
          float(base.outer.nb[-1]*base.outer.F[-1]))
    right=tuple(float(s[j]) for s in (matched.nb,matched.p,matched.eps))
    bridge=CausalExponentialInterval(*left,*right)
    if bridge.validate()["status"]!="PASS":raise ValueError("frozen H/E outer connection is not causal")
    issharp=transition_index is not None
    if issharp and transition_index<=j:raise ValueError("sharp transition lies below selected source start")
    c=dict(surface=copy.deepcopy(source_ref["surface"]),outer=copy.deepcopy(source_ref["outer"]),
           lower_bridge=bridge.contract(),
           upper=dict(n=matched.nb[j:].tolist(),p=matched.p[j:].tolist(),epsilon=matched.eps[j:].tolist(),
                      sharp_interfaces=[int(transition_index-j)] if issharp else []))
    physics=WholeCrustBarotrope.from_contract(c)
    return physics,c
def star(eos):
    def solve(n):return solve_sharp_background(eos,n,**TOV)
    l,h=.22,.5
    if not solve(l).mass_msun<1.4<solve(h).mass_msun:
        raise ValueError("No first ascending M=1.4 bracket on fixed EOS support")
    nc=brentq(lambda n:solve(n).mass_msun-1.4,l,h,xtol=1e-9)
    bg=solve(nc)
    if abs(bg.mass_msun-1.4)>2e-7:raise ValueError("target mass residual exceeds tolerance")
    return bg
def root_admitted(x):
    return (x["root_converged"] and abs(x["matching_residual"])<1e-6
            and x["relative_join_mismatch"]<1e-5
            and x["surface_BC_relative"]<1e-5
            and all(j["traction_relative_residual"]<1e-5 and abs(j["kinematic_residual"])<1e-5
                    for j in x["interface_junctions"]))
def mode_calc(bg,method,target):
    sharp=bool(bg.eos.transitions)
    conversions=("slow","rapid") if sharp else ("slow",)
    result={}
    for conversion in conversions:
        for label,(a,b) in WINDOWS.items():
            key=conversion+"_"+label
            try:
                xa=matching_residual(bg,a,COARSE,conversion=conversion)
                xb=matching_residual(bg,b,COARSE,conversion=conversion)
                if xa*xb>=0:
                    pts=np.linspace(a,b,7)
                    vals=[matching_residual(bg,f,COARSE,conversion=conversion) for f in pts]
                    brackets=[(float(z),float(w)) for z,w,q,s in zip(pts[:-1],pts[1:],vals[:-1],vals[1:]) if q*s<0]
                    if len(brackets)!=1:raise ValueError(f"frequency root not unique in prescribed {label} window")
                    a,b=brackets[0]
                q=solve_mode(bg,(a,b),COARSE,conversion=conversion,profile_points=140)
                save(target/f"{method}_{key}_COARSE.json",q)
                t=solve_mode(bg,(q["frequency_hz"]*.996,q["frequency_hz"]*1.004),FINE,
                             conversion=conversion,profile_points=140)
                save(target/f"{method}_{key}_TIGHT.json",t)
                record={k:v for k,v in t.items() if k!="profile"}
                record["status"]="PASS_SCOPED" if root_admitted(t) else "REJECTED_RESIDUAL"
                record["method"]=method;record["frequency_window"]=[a,b]
                record["frequency_control_delta_hz"]=t["frequency_hz"]-q["frequency_hz"]
                record["requested_mode_family"]=label
                result[key]=record
                print("NATIVE_STAR_MODE",method,key,record["status"],record["frequency_hz"],flush=True)
            except Exception as exc:
                result[key]=dict(status="BLOCKED",type=type(exc).__name__,error=str(exc))
                print("NATIVE_STAR_MODE_FAIL",method,key,str(exc)[:140],flush=True)
    return result
def run_one(crust):
    core="PK1_CANONICAL"
    pair=f"{crust}__{core}"
    dest=OUT/pair
    if dest.exists():raise FileExistsError(dest)
    dest.mkdir(parents=True)
    provenance=json.loads((NATIVE/pair/"RESULT.json").read_text())
    raw_crust=read_compose(source("crust",crust))
    raw_core=read_compose(source("core",core))
    result=dict(pair=pair,status="RUNNING",
        fixed_source_sha256=provenance["input_hashes"],
        actual_native_methods=["cubic_F","quintic_F","cubic_P","quintic_P","css","native_maxwell"],
        outer_reference_sha256=sha(REFERENCE),lower_join_source_index=SOURCE_J[crust],
        smooth_vs_sharp_not_same_physics=True,
        mass_target_msun=1.4,closure="equilibrium perfect fluid Gamma1_eq, l2 Cowling",
        method_results={})
    save(dest/"RESULT.json",result)
    for method in result["actual_native_methods"]:
        entry=dict(status="RUNNING");result["method_results"][method]=entry
        folder=dest/method;folder.mkdir()
        try:
            if method=="native_maxwell":
                matched,rep=construct_maxwell_transition(_thermo_only_view(raw_crust),
                    _thermo_only_view(raw_core),
                    assume_local_charge_neutrality=True,return_report=True)
                transition=int(rep["low_phase_rows"])-1
                info=dict(native_maxwell_report=rep,transition_index=transition)
                write_compose(folder/"compose",matched)
            else:
                info=provenance["methods"][method]["native_report"]
                matched=read_compose(NATIVE/pair/method/"compose")
                transition=int(info["css_transition_low_index"]) if method=="css" else None
            physics,contract=whole_from_match(crust,matched,transition,method)
            entry["constructed_full_eos"]=dict(source_rows=matched.n,sharp_transition_index=transition,
                n_join=physics.bridge.n_right,p_join=physics.bridge.p_right,
                source_whole_upper_n=physics.upper_n.tolist()[:2],
                contract_hash=hashlib.sha256(json.dumps(contract,sort_keys=True).encode()).hexdigest(),
                native_method=method)
            save(folder/"REPRESENTATION.json",contract)
            bg=star(SharpContractEnthalpyEOS(physics))
            entry["background"]=bg.summary()
            entry["status"]="TOV_COMPLETE"
            print("NATIVE_STAR",pair,method,bg.mass_msun,bg.radius_km,len(bg.interfaces),flush=True)
            entry["modes"]=mode_calc(bg,method,folder)
            entry["status"]="COMPLETED_WITH_MODE_STATUSES"
        except Exception as exc:
            entry["status"]="REJECTED_OR_BLOCKED";entry["type"]=type(exc).__name__
            entry["error"]=str(exc);entry["traceback"]=traceback.format_exc()[-2200:]
        save(dest/"RESULT.json",result)
    result["status"]="COMPLETED"
    save(dest/"RESULT.json",result)
    return result
if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--crust",choices=list(SOURCE_J),required=True)
    a=ap.parse_args();run_one(a.crust)
