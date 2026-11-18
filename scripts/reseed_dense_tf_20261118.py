"""One independent physics-preserving optimizer seed per rejected native N128 TF point.

Uses two *previously accepted* neighboring profiles at same source branch/mesh
instead of nearest-neighbor-only. Physical stationarity tolerances unchanged.
No false acceptance for optimizer mode-8 or interpolated profiles.
"""
from pathlib import Path
import argparse,json,time,traceback,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from run_dense_tf_20261118 import CASES,OUT,EXE,existing_anchors,sha,save
from compose_eos.tf_candidate import run_state,weights

def bracket(anchors,n):
    left=[s for s in anchors if s["n_b"]<n-1e-12]
    right=[s for s in anchors if s["n_b"]>n+1e-12]
    if not left or not right:raise ValueError("no native source profile bracket")
    l=max(left,key=lambda a:a["n_b"]);r=min(right,key=lambda a:a["n_b"])
    if r["n_b"]-l["n_b"]>.00201:
        raise ValueError("bracket too wide to serve as bounded initialization")
    return l,r
def calculate(case):
    spec=CASES[case];base=json.loads((OUT/case/"RESULT.json").read_text())
    anchors=existing_anchors(spec)
    folder=OUT/"refined_source_bracket"/case
    if folder.exists():raise FileExistsError(folder)
    folder.mkdir(parents=True)
    outcomes=[]
    result=dict(schema="ns_cowling.dense_tf_two_anchor_refinement.v1",case=case,
        original_result_sha256=sha(OUT/case/"RESULT.json"),
        original_rejections=sum(q["status"]=="REJECTED_NATIVE_STATE" for q in base["records"]),
        source="Only original accepted same-branch N128 H/G profiles",
        attempts_per_rejected_state=1,
        original_optimizer_tolerance=1e-14,
        status="RUNNING",outcomes=outcomes)
    save(folder/"REFINEMENT.json",result)
    for i,item in enumerate(base["records"]):
        if item["status"]!="REJECTED_NATIVE_STATE":continue
        n=item["n_fm3"]
        record=dict(n_fm3=n,original_status=item["status"],status="NOT_RUN")
        path=folder/f"index_{i:03d}"
        try:
            lo,hi=bracket(anchors,n)
            t=(n-lo["n_b"])/(hi["n_b"]-lo["n_b"])
            ll=np.asarray(lo["profile"],float)
            rr=np.asarray(hi["profile"],float)
            both=(1-t)*ll[:,:2]+t*rr[:,:2]
            w=weights(128,3)
            both*=n/(w@(both[:,0]+both[:,1]))
            radius=(1-t)*float(lo["radius"])+t*float(hi["radius"])
            parent=dict(seed_method="interpolation_between_two_frozen_native_stationary_profiles",
                left_sha256=lo["_source_sha256"],right_sha256=hi["_source_sha256"],
                left_relative_path=lo["_source_path"],right_relative_path=hi["_source_path"],
                density_bracket_fm3=[lo["n_b"],hi["n_b"]],weight=t,
                no_physical_tolerance_change=True)
            s=run_state(EXE,spec["model"],n,128,3,"droplet",radius,
              work_dir=path,initial_nn=both[:,0],initial_np=both[:,1],
              parent=parent,maxiter=1800,ftol=1e-14,timeout=180)
            okay=bool(s.get("accepted_stationary") and s.get("geometry")=="droplet")
            record.update(status="ACCEPTED_NATIVE_STATIONARY" if okay else "REJECTED_NATIVE_STATE",
                 native_status=s.get("status"),optimizer_mode=s.get("optimizer_mode"),
                 kkt=s.get("kkt"),parent=parent,
                 state_path=str(Path(s["state_path"]).relative_to(ROOT)),
                 state_sha256=sha(s["state_path"]),native_exit=s.get("actual_exit_code"))
        except Exception as exc:
            record.update(status="NOT_APPLICABLE_OR_ENVIRONMENT_ERROR",error=str(exc),
                          error_type=type(exc).__name__,traceback=traceback.format_exc()[-1200:])
        outcomes.append(record)
        save(folder/"REFINEMENT.json",result)
        print("RESEED",case,round(n,9),record["status"],flush=True)
    result.update(status="SCOPED_COMPLETED_REJECTIONS_RETAINED",
        recovered_count=sum(x["status"]=="ACCEPTED_NATIVE_STATIONARY" for x in outcomes),
        still_rejected_count=sum(x["status"]=="REJECTED_NATIVE_STATE" for x in outcomes),
        not_applicable_count=sum(x["status"]=="NOT_APPLICABLE_OR_ENVIRONMENT_ERROR" for x in outcomes))
    save(folder/"REFINEMENT.json",result)
    print("RESEED_SUMMARY",case,result["recovered_count"],result["still_rejected_count"],result["not_applicable_count"],flush=True)
    return result
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--case",choices=list(CASES),required=True)
    args=ap.parse_args()
    calculate(args.case)
