"""Classify ambiguous p1 frequency roots by eigenfunction nodes, not root order.

Repair is separately versioned so original search failures and raw diagnostics
remain immutable. Only actual node-classified, residual-admitted (W,V)=(1,1)
roots can enter the f/p1 comparison table. No spectrum-completeness claim.
"""
from __future__ import annotations
from pathlib import Path
import sys,json,hashlib,time,argparse
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background,matching_residual,solve_mode
from stellar_quintic_contract_20261118 import install_quintic_upper
from run_native_matching_stellar_20261118 import COARSE,FINE,TOV,root_admitted
from run_hermite_modes_bounded_20261118 import mode_pass
H=ROOT/"data/low_density_crust_20261117/stellar/causal_extension"
OUT=ROOT/"data/cubic_quintic_matching_20261118/p1_node_refinement"
STELLAR=ROOT/"data/cubic_quintic_matching_20261118/stellar_bounded"
NATIVE=ROOT/"data/cubic_quintic_matching_20261118/native_stellar"
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def load_task(group,case,method):
    if group=="H":
        p=STELLAR/case/"RESULT.json"
        x=json.loads(p.read_text());z=x["cases"][method]
        root=json.loads((H/"TOV_PLAN.json").read_text())
        contract=next(i["contract"] for i in root["cases"] if i["id"]==case)
        bar=WholeCrustBarotrope.from_contract(json.loads((ROOT/contract).read_text()))
        if method=="quintic":
            bar,admission=install_quintic_upper(bar)
            if admission["rejected"]:raise ValueError("quintic source not fully admitted")
        accept=mode_pass
    else:
        p=NATIVE/case/"RESULT.json"
        x=json.loads(p.read_text());z=x["method_results"][method]
        bar=WholeCrustBarotrope.from_contract(json.loads((NATIVE/case/method/"REPRESENTATION.json").read_text()))
        accept=root_admitted
    if z.get("status") not in ("COMPUTED_WITH_REJECTED_ROOTS","COMPLETED_WITH_MODE_STATUSES","TOV_COMPLETE"):
        raise ValueError(f"background not admitted: {z.get('status')}")
    nc=float(z["background"]["central_density_fm3"])
    bg=solve_sharp_background(SharpContractEnthalpyEOS(bar),nc,**TOV)
    return bg,z,accept
def one(group,case,method,conversion):
    dest=OUT/group/case/method/conversion
    if dest.exists():raise FileExistsError(dest)
    dest.mkdir(parents=True)
    result=dict(group=group,case=case,method=method,conversion=conversion,
        status="RUNNING",classification="p1 requires exactly one (W_nodes,V_nodes)=(1,1)",
        native_source_not_changed=True,cs2_closure="Gamma1_eq",candidates=[])
    save(dest/"RESULT.json",result)
    try:
        bg,earlier,pass_fn=load_task(group,case,method)
        source_modes=earlier.get("modes",{})
        key=conversion+"_p1"
        if key in source_modes and source_modes[key]["status"] in ("PASS","PASS_SCOPED"):
            result.update(status="ALREADY_ADMITTED_PREVIOUS_RUN",selected=source_modes[key])
            save(dest/"RESULT.json",result);return result
        grid=np.linspace(4450,7850,12)
        values=[matching_residual(bg,float(f),COARSE,conversion=conversion) for f in grid]
        brackets=[(float(a),float(b)) for a,b,fa,fb in zip(grid[:-1],grid[1:],values[:-1],values[1:]) if fa*fb<0]
        result["scan_brackets"]=brackets;result["fixed_scan_grid"]=grid.tolist()
        if len(brackets)>3:raise ValueError("too many roots for bounded classification")
        for i,bracket in enumerate(brackets):
            row=dict(candidate_bracket=bracket)
            try:
                m=solve_mode(bg,bracket,COARSE,conversion=conversion,profile_points=200)
                save(dest/f"CANDIDATE_{i:02d}_COARSE.json",m)
                nW=int(m["W_sign_changes_including_interfaces"])
                nV=int(m["V_sign_changes_including_interfaces"])
                row.update(frequency_hz=m["frequency_hz"],W_nodes=nW,V_nodes=nV,
                           coarse_admitted=bool(pass_fn(m)))
                if (nW,nV)==(1,1) and pass_fn(m):
                    tight=solve_mode(bg,(m["frequency_hz"]*.998,m["frequency_hz"]*1.002),FINE,conversion=conversion,profile_points=200)
                    save(dest/f"CANDIDATE_{i:02d}_TIGHT.json",tight)
                    row["tighter_admitted"]=bool(pass_fn(tight))
                    row["tighter_frequency_hz"]=tight["frequency_hz"]
                    row["delta_numerical_hz"]=tight["frequency_hz"]-m["frequency_hz"]
                    row["tighter_nodes"]=[tight["W_sign_changes_including_interfaces"],tight["V_sign_changes_including_interfaces"]]
                    if row["tighter_nodes"]!=[1,1]:row["tighter_admitted"]=False
            except Exception as exc:
                row["error"]=type(exc).__name__+": "+str(exc)
            result["candidates"].append(row)
            save(dest/"RESULT.json",result)
        accepted=[r for r in result["candidates"] if r.get("tighter_admitted")]
        if len(accepted)==1:
            result.update(status="PASS_SCOPED_NODE_CLASSIFIED_P1",selected=accepted[0])
        else:
            result.update(status="UNRESOLVED_OR_REJECTED",accepted_candidate_count=len(accepted))
    except Exception as exc:
        result.update(status="BLOCKED",error=str(exc),error_type=type(exc).__name__)
    save(dest/"RESULT.json",result)
    print("P1_CLASS",group,case,method,conversion,result["status"],
          result.get("selected",{}).get("tighter_frequency_hz"),flush=True)
    return result
def main():
    from concurrent.futures import ProcessPoolExecutor,as_completed
    jobs=[]
    hcases=json.loads((H/"TOV_PLAN.json").read_text())
    for case in (x["id"] for x in hcases["cases"] if x["label"]=="primary"):
        for method in ("cubic","quintic"):
            for conv in ("slow","rapid"):
                jobs.append(("H",case,method,conv))
    for case in ("SPG_M1__PK1_CANONICAL","SPG_M3__PK1_CANONICAL"):
        for method in ("cubic_P","quintic_P","cubic_F","quintic_F","css","native_maxwell"):
            for conv in (("slow","rapid") if method in ("css","native_maxwell") else ("slow",)):
                jobs.append(("native",case,method,conv))
    results=[]
    with ProcessPoolExecutor(max_workers=5) as pool:
        for f in as_completed([pool.submit(one,*job) for job in jobs]):
            x=f.result();results.append(dict(group=x["group"],case=x["case"],method=x["method"],
               conversion=x["conversion"],status=x["status"]))
    save(OUT/"CAMPAIGN.json",dict(status="FINISHED_WITH_PRESERVED_REJECTIONS",jobs=len(jobs),results=results))
    print("P1_NODE_CLASSIFICATION_DONE",len(results),sum(x["status"]=="PASS_SCOPED_NODE_CLASSIFIED_P1" for x in results))
if __name__=="__main__":main()
