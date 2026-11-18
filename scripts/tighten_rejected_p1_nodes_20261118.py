"""Bounded targeted refinement of 1-node p1 candidates rejected by coarse residuals.

No gate relaxation. Retains original scan/rejection and new tight raw JSON.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import json,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from rescue_p1_node_classified_20261118 import load_task,FINE,OUT as PREVIOUS,save
from compose_eos.stellar_sharp import solve_mode
OUT=ROOT/"data/cubic_quintic_matching_20261118/p1_tight_controls"
def work(item):
    x=json.loads(Path(item).read_text())
    target=OUT/x["group"]/x["case"]/x["method"]/x["conversion"]
    if target.exists():raise FileExistsError(target)
    target.mkdir(parents=True)
    rec=dict(group=x["group"],case=x["case"],method=x["method"],conversion=x["conversion"],
             prior_status=x["status"],prior_source=str(Path(item).relative_to(ROOT)),
             status="RUNNING",candidates=[])
    save(target/"RESULT.json",rec)
    try:
        bg,previous,passed=load_task(x["group"],x["case"],x["method"])
        old=[z for z in x["candidates"] if
             (z.get("W_nodes"),z.get("V_nodes"))==(1,1) and
             not z.get("tighter_admitted") and z.get("frequency_hz")]
        for i,candidate in enumerate(old):
            q=dict(original_bracket=candidate["candidate_bracket"],
                   earlier_frequency=candidate["frequency_hz"],raw_original_coarse_admission=False)
            try:
                f=candidate["frequency_hz"]
                trial=solve_mode(bg,(f*.997,f*1.003),FINE,conversion=x["conversion"],profile_points=220)
                save(target/f"TIGHT_{i:02d}_RAW.json",trial)
                q.update(frequency_hz=trial["frequency_hz"],
                    W_nodes=trial["W_sign_changes_including_interfaces"],
                    V_nodes=trial["V_sign_changes_including_interfaces"],
                    root_converged=trial["root_converged"],
                    relative_join_mismatch=trial["relative_join_mismatch"],
                    surface_BC_relative=trial["surface_BC_relative"],
                    status="PASS" if passed(trial) and
                        (trial["W_sign_changes_including_interfaces"],trial["V_sign_changes_including_interfaces"])==(1,1)
                         else "REJECTED_REFINED")
            except Exception as e:
                q.update(status="BLOCKED",exception_type=type(e).__name__,error=str(e))
            rec["candidates"].append(q)
        approved=[v for v in rec["candidates"] if v["status"]=="PASS"]
        if len(approved)==1:rec.update(status="PASS_SCOPED_NODE_CLASSIFIED_P1",selected=approved[0])
        else:rec.update(status="UNRESOLVED",accepted_count=len(approved))
    except Exception as exc:
        rec.update(status="BLOCKED",error=str(exc),error_type=type(exc).__name__)
    save(target/"RESULT.json",rec)
    print("P1_TIGHT",rec["group"],rec["case"],rec["method"],rec["conversion"],rec["status"],flush=True)
    return rec
def main():
    paths=[]
    for p in sorted(PREVIOUS.glob("**/RESULT.json")):
        x=json.loads(p.read_text())
        if x["status"]!="UNRESOLVED_OR_REJECTED":continue
        if any((v.get("W_nodes"),v.get("V_nodes"))==(1,1)
               and not v.get("tighter_admitted") for v in x["candidates"]):
            paths.append(p)
    with ProcessPoolExecutor(max_workers=5) as pool:
        outcomes=[p.result() for p in as_completed([pool.submit(work,x) for x in paths])]
    save(OUT/"CAMPAIGN.json",dict(status="COMPLETED_WITH_PRESERVED_REJECTIONS",requested=len(paths),
          accepted=sum(z["status"]=="PASS_SCOPED_NODE_CLASSIFIED_P1" for z in outcomes),
          results=[dict(case=z["case"],method=z["method"],status=z["status"]) for z in outcomes]))
    print("P1_TIGHT_DONE",len(paths),sum(z["status"]=="PASS_SCOPED_NODE_CLASSIFIED_P1" for z in outcomes))
if __name__=="__main__":main()
