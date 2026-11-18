"""Matched-mass Cowling f/p1 sensitivity on four real H sharp Maxwell EOSs.

Physics held identical: H lower bridge/outer/surface, source Gibbs/Maxwell
doublet, bulk Gamma1=Gamma_eq, l=2, perfect fluid. Only cubic vs
quintic P(mu) evaluation of within-phase *upper* EOS intervals changes.
Slow and rapid phase conversion computed separately. Native CSS is
not relabeled as this H Maxwell background.
"""
from __future__ import annotations
from pathlib import Path
from dataclasses import asdict
import sys,json,hashlib,time,traceback
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from compose_eos.whole_crust import WholeCrustBarotrope
from stellar_quintic_contract_20261118 import install_quintic_upper
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background,scan_brackets,solve_mode
from compose_eos.stellar_cowling import CowlingSettings
H=ROOT/"data/low_density_crust_20261117/stellar/causal_extension"
OUT=ROOT/"data/cubic_quintic_matching_20261118/stellar"
SETTINGS=dict(r0_km=.001,rtol=1e-10,atol=1e-12,max_step_h=.002)
COWLING=CowlingSettings(ell=2,rtol=3e-9,atol=3e-12,max_step_h=.001,
                        surface_h=3e-8,root_xtol_hz=1e-7)
REFINE=CowlingSettings(ell=2,rtol=3e-10,atol=3e-13,max_step_h=.0005,
                       surface_h=3e-8,root_xtol_hz=1e-8)
FREQ={"f":np.linspace(1550,2900,17),"p1":np.linspace(4500,7700,22)}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,rec):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(rec,indent=2,allow_nan=False)+"\n")
def fixed_mass(bg_eos,ref_center):
    lower=max(.15,ref_center*.94);upper=min(.8,ref_center*1.06)
    def solve(n):return solve_sharp_background(bg_eos,n,**SETTINGS)
    ml=solve(lower).mass_msun;mr=solve(upper).mass_msun
    if not ml<1.4<mr:raise ValueError(f"fixed target mass not bracketed in original supported region {ml},{mr}")
    n=float(brentq(lambda n:solve(n).mass_msun-1.4,lower,upper,xtol=1e-10))
    result=solve(n)
    if abs(result.mass_msun-1.4)>1e-7:raise ValueError("fixed mass root miss")
    return result
def admit(root):
    return bool(root["root_converged"] and
                abs(root["matching_residual"])<1e-6 and
                root["relative_join_mismatch"]<1e-5 and
                root["surface_BC_relative"]<1e-5 and
                all(z["traction_relative_residual"]<1e-5 and
                    abs(z["kinematic_residual"])<1e-5 for z in root["interface_junctions"]))
def run_one(entry):
    cid=entry["id"]
    target=OUT/cid
    if target.exists():raise FileExistsError(target)
    target.mkdir(parents=True)
    rec=dict(case=cid,contract_path=entry["contract"],source_sha256=entry["contract_sha256"],
             baseline_H_commit="634a586fdd0c214bf7fa9718b02ebafa3e927008",
             source_sharp_preserved=True,scope="Upper crust/core P(mu) interpolation only, fixed Maxwell transition and H lower closure",
             settings=dict(tov=SETTINGS,cowling=asdict(COWLING),refine=asdict(REFINE)),
             status="RUNNING",orders={})
    save(target/"RESULT.json",rec)
    original=json.loads((ROOT/entry["contract"]).read_text())
    if sha(ROOT/entry["contract"])!=entry["contract_sha256"]:raise ValueError("source identity changed")
    original_bar=WholeCrustBarotrope.from_contract(original)
    quintic_bar,upgrade=install_quintic_upper(original_bar)
    if upgrade["rejected"]:raise ValueError("quintic upper source phase contains noncausal intervals; do not mix laws")
    ref=json.loads((H/"TOV_RESULTS.json").read_text())
    center=next(q["background"]["central_density_fm3"] for q in ref["cases"] if q["id"]==cid)
    for order,b in (("cubic",original_bar),("quintic",quintic_bar)):
        st=time.monotonic()
        item=dict(status="RUNNING",upgrade=upgrade if order=="quintic" else None,modes={})
        rec["orders"][order]=item
        save(target/"RESULT.json",rec)
        try:
            eos=SharpContractEnthalpyEOS(b)
            bg=fixed_mass(eos,center)
            item["background"]=bg.summary()
            item["center_n_fm3"]=bg.central_density_fm3
            item["status"]="BACKGROUND_PASS"
            print("TOV",cid,order,bg.mass_msun,bg.radius_km,flush=True)
            for conversion in ("slow","rapid"):
                for family,grid in FREQ.items():
                    key=f"{conversion}_{family}"
                    try:
                        brackets,scan=scan_brackets(bg,grid,COWLING,conversion=conversion)
                        if len(brackets)!=1:raise ValueError(f"require exactly one {family} root in predeclared window; got {brackets}")
                        mode=solve_mode(bg,brackets[0],COWLING,conversion=conversion,profile_points=400)
                        accepted=admit(mode)
                        compact={k:v for k,v in mode.items() if k!="profile"}
                        compact.update(status="PASS" if accepted else "REJECTED_RESIDUAL",
                                     selected_window_hz=[float(grid[0]),float(grid[-1])],
                                     scan_brackets=brackets,
                                     expected_mode_family=family)
                        # Always retain raw mode profile and scan, even for rejected attempts.
                        save(target/f"{order}_{key}_RAW.json",mode)
                        save(target/f"{order}_{key}_SCAN.json",dict(brackets=brackets,values=scan))
                        if accepted:
                            ctl=solve_mode(bg,[mode["frequency_hz"]*.998,mode["frequency_hz"]*1.002],
                                  REFINE,conversion=conversion,profile_points=320)
                            control={k:v for k,v in ctl.items() if k!="profile"}
                            control["status"]="PASS" if admit(ctl) else "REJECTED_RESIDUAL"
                            save(target/f"{order}_{key}_TIGHT_RAW.json",ctl)
                            compact["tighter_control"]=control
                            if control["status"]!="PASS":
                                compact["status"]="REJECTED_REFINEMENT"
                        item["modes"][key]=compact
                        print("MODE",cid,order,key,compact["status"],
                              compact.get("frequency_hz"),flush=True)
                    except Exception as exc:
                        item["modes"][key]=dict(status="BLOCKED",error_type=type(exc).__name__,error=str(exc))
                        print("MODE_BLOCKED",cid,order,key,repr(exc),flush=True)
                    save(target/"RESULT.json",rec)
            item["status"]="MODES_COMPUTED_WITH_REJECTIONS"
        except Exception as exc:
            item.update(status="BACKGROUND_FAILED",error_type=type(exc).__name__,error=str(exc),
                        traceback=traceback.format_exc()[-2800:])
        item["wall_seconds"]=time.monotonic()-st
        save(target/"RESULT.json",rec)
    rec["status"]="COMPLETED_SCOPED"
    save(target/"RESULT.json",rec)
    return rec
def main():
    from concurrent.futures import ThreadPoolExecutor,as_completed
    plan=json.loads((H/"TOV_PLAN.json").read_text())
    cases=[c for c in plan["cases"] if c["label"]=="primary"]
    with ThreadPoolExecutor(max_workers=4) as pool:
        ret=[]
        for f in as_completed([pool.submit(run_one,c) for c in cases]):
            x=f.result()
            ret.append(dict(case=x["case"],statuses={k:v["status"] for k,v in x["orders"].items()}))
    save(OUT/"CAMPAIGN.json",dict(status="COMPLETED",results=ret))
    print("STELLAR_DONE",len(ret))
if __name__=="__main__":main()
