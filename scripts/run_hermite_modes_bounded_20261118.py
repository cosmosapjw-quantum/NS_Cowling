"""Finite 1.4-Msun paired cubic/quintic Cowling f/p1 reference.

This is NOT native CSS matching: fixed Maxwell H strict whole-EOS source,
unchanged lower crust/surface, same gamma1=gamma_eq and interface conversion.
Frequencies found in predeclared broad windows, with tighter refinement.
"""
from __future__ import annotations
from pathlib import Path
from dataclasses import asdict
import sys,json,hashlib,time,traceback,argparse
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background,matching_residual,solve_mode
from compose_eos.stellar_cowling import CowlingSettings
from stellar_quintic_contract_20261118 import install_quintic_upper
H=ROOT/"data/low_density_crust_20261117/stellar/causal_extension"
OUT=ROOT/"data/cubic_quintic_matching_20261118/stellar_bounded"
TOV=dict(r0_km=.001,rtol=1e-9,atol=1e-11,max_step_h=.005)
COARSE=CowlingSettings(ell=2,rtol=2e-7,atol=2e-10,max_step_h=.004,surface_h=3e-8,root_xtol_hz=1e-6)
FINE=CowlingSettings(ell=2,rtol=3e-9,atol=3e-12,max_step_h=.001,surface_h=3e-8,root_xtol_hz=1e-8)
WINDOWS={"f":[1550.,2900.],"p1":[4600.,7800.]}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def mode_pass(m):
    return bool(m["root_converged"] and abs(m["matching_residual"])<1e-6
                and m["relative_join_mismatch"]<1e-5 and m["surface_BC_relative"]<1e-5
                and all(x["traction_relative_residual"]<1e-5
                        and abs(x["kinematic_residual"])<1e-5 for x in m["interface_junctions"]))
def mass_star(bar,nc):
    e=SharpContractEnthalpyEOS(bar)
    def calc(x):return solve_sharp_background(e,x,**TOV)
    l,h=nc*.9,nc*1.1
    if not calc(l).mass_msun<1.4<calc(h).mass_msun:
        raise ValueError("same target mass bracket missing")
    root=brentq(lambda n:calc(n).mass_msun-1.4,l,h,xtol=1e-9)
    bg=calc(root)
    if abs(bg.mass_msun-1.4)>2e-7:raise ValueError("fixed mass root not accurate")
    return bg
def single(case):
    cid=case["id"]
    target=OUT/cid
    if target.exists():raise FileExistsError("prior bounded mode evidence exists")
    target.mkdir(parents=True)
    hcase=next(x for x in json.loads((H/"TOV_RESULTS.json").read_text())["cases"] if x["id"]==cid)
    nc=float(hcase["background"]["central_density_fm3"])
    physical=ROOT/case["contract"]
    if sha(physical)!=case["contract_sha256"]:raise ValueError("H source contract changed")
    cubic=WholeCrustBarotrope.from_contract(json.loads(physical.read_text()))
    quintic,upgrade=install_quintic_upper(cubic)
    if upgrade["rejected"]:raise ValueError("Cannot admit noncausal quintic source")
    res=dict(status="RUNNING",case=cid,source_contract=case["contract"],
        source_sha256=case["contract_sha256"],source_basis="H actual source Maxwell sharp 1.4-Msun",
        upper_quintic_installation=upgrade,cubic_vs_quintic_only_upper_smooth_intervals=True,
        fixed_outer_lower_surface=True,source_transition_unchanged=True,
        TOV_settings=TOV,Cowling_coarse=asdict(COARSE),Cowling_fine=asdict(FINE),
        physical_closure="Gamma1_eq perfect fluid, l=2, slow/rapid conversion, Cowling delta_g=0",
        cases={})
    save(target/"RESULT.json",res)
    for method,bar in (("cubic",cubic),("quintic",quintic)):
        row=dict(status="RUNNING",modes={});res["cases"][method]=row
        try:
            bg=mass_star(bar,nc)
            row["background"]=bg.summary()
            row["status"]="BACKGROUND_PASS"
            print("MSTAR",cid,method,bg.mass_msun,bg.radius_km,flush=True)
            for conv in ("slow","rapid"):
                for family,window in WINDOWS.items():
                    key=conv+"_"+family
                    try:
                        a,b=window
                        fa=matching_residual(bg,a,COARSE,conversion=conv)
                        fb=matching_residual(bg,b,COARSE,conversion=conv)
                        if fa*fb>=0:
                            # bounded diagnostic 5 subsegments, not a blind exhaustive spectrum search
                            points=np.linspace(a,b,6)
                            vals=[matching_residual(bg,f,COARSE,conversion=conv) for f in points]
                            choices=[(float(x),float(y)) for x,y,v,w in zip(points[:-1],points[1:],vals[:-1],vals[1:]) if v*w<0]
                            if len(choices)!=1:raise ValueError(f"no unique root in fixed {family} window {choices}")
                            a,b=choices[0]
                        mode=solve_mode(bg,(a,b),COARSE,conversion=conv,profile_points=180)
                        save(target/f"{method}_{key}_COARSE.json",mode)
                        tight=solve_mode(bg,(mode["frequency_hz"]*.997,mode["frequency_hz"]*1.003),
                                          FINE,conversion=conv,profile_points=180)
                        save(target/f"{method}_{key}_FINE.json",tight)
                        m={k:v for k,v in tight.items() if k!="profile"}
                        m.update(status="PASS" if mode_pass(tight) else "REJECTED_PHYSICAL_OR_NUMERIC",
                                 coarse_frequency_hz=mode["frequency_hz"],
                                 frequency_control_delta_hz=tight["frequency_hz"]-mode["frequency_hz"],
                                 mode_family_requested=family,
                                 source_H_reference="fixed P(mu) source support; no new Maxwell matching algorithm",
                                 original_window_hz=window)
                        row["modes"][key]=m
                        print("MSTAR_MODE",cid,method,key,m["status"],m["frequency_hz"],
                              m.get("W_sign_changes_including_interfaces"),m.get("V_sign_changes_including_interfaces"),flush=True)
                    except Exception as exc:
                        row["modes"][key]=dict(status="BLOCKED",error_type=type(exc).__name__,error=str(exc))
                        print("MSTAR_MODE_REJECT",cid,method,key,type(exc).__name__,str(exc)[:150],flush=True)
                    save(target/"RESULT.json",res)
            row["status"]="COMPUTED_WITH_REJECTED_ROOTS"
        except Exception as exc:
            row["status"]="BACKGROUND_REJECTED";row["error"]=str(exc)
        save(target/"RESULT.json",res)
    res["status"]="COMPLETED_WITH_EXPLICIT_MODE_STATUSES"
    save(target/"RESULT.json",res)
    return res
if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--case",required=True)
    args=ap.parse_args()
    plans=json.loads((H/"TOV_PLAN.json").read_text())
    c=next((x for x in plans["cases"] if x["id"]==args.case and x["label"]=="primary"),None)
    if c is None:raise SystemExit("unapproved H case")
    single(c)
