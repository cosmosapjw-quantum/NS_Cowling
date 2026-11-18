"""Blind interpolation validation against new spatial-TF stationary N128 states.

Frozen H/G branch-matched source stationary states are training knots; newly
computed same-branch TF states at previously uncomputed densities are holdouts.
Rejected optimizer outputs NEVER become reference thermodynamic states.
"""
import json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from run_dense_tf_20261118 import CASES,existing_anchors,OUT,sha
from benchmark_native_dense_eos_20261118 import (save,GLOBAL_METHODS,SHAPE_METHODS,predict_global,
            prediction_from_shape,native_truth_cs2,summarize)
OUTPUT=ROOT/"data/native_dense_interpolation_20261118/tf_holdout"
def collect(case):
    spec=CASES[case]
    start=spec["n_low"]-spec["step"]*2
    end=spec["n_high"]+spec["step"]*2
    orig=existing_anchors(spec)
    anchors={}
    for v in orig:
        n=round(float(v["n_b"]),9)
        if not start<=n<=end:continue
        if n not in anchors or v["epsilon"]<anchors[n]["epsilon"]:
            anchors[n]=v
    base=json.loads((OUT/case/"RESULT.json").read_text())
    refpath=OUT/"refined_source_bracket"/case/"REFINEMENT.json"
    recovered=json.loads(refpath.read_text())["outcomes"] if refpath.is_file() else []
    recovered={round(x["n_fm3"],9):x for x in recovered if x["status"]=="ACCEPTED_NATIVE_STATIONARY"}
    values=[]
    for item in base["records"]:
        n=round(item["n_fm3"],9)
        src=None
        if item["status"]=="ACCEPTED_NATIVE_STATIONARY":
            src=ROOT/item["accepted_path"]
        elif n in recovered:
            src=ROOT/recovered[n]["state_path"]
        if src is not None:
            s=json.loads(src.read_text())
            if not s.get("accepted_stationary") or s.get("points")!=128 or s.get("geometry")!="droplet":
                raise ValueError("not accepted same-mesh native TF holdout")
            values.append(dict(n=n,state=s,path=str(src.relative_to(ROOT)),sha=sha(src)))
    return anchors,values,base,recovered
def run(case):
    start_time=time.monotonic()
    anchors,values,base,recovered=collect(case)
    srcrows=sorted(anchors.items())
    src_n=np.array([x for x,s in srcrows])
    train={float(k):dict(n=float(k),p=float(s["pressure"]),epsilon=float(s["epsilon"]),
                         mu=float(s["mu_B"]),geometry=s["geometry"],origin="H/G_NATIVE_SOURCE",
                         sha=s["_source_sha256"],path=s["_source_path"]) for k,s in srcrows}
    seen={};repeat=[];excluded=[]
    for item in values:
        n=item["n"];s=item["state"]
        if n in train:
            ref=train[n]
            repeat.append(dict(n=n,native_ref_sha=item["sha"],
                delta_P=float(s["pressure"]-ref["p"]),delta_epsilon=float(s["epsilon"]-ref["epsilon"]),
                delta_mu=float(s["mu_B"]-ref["mu"])))
            continue
        if n not in seen:
            seen[n]=dict(n=n,p=float(s["pressure"]),epsilon=float(s["epsilon"]),
                         mu=float(s["mu_B"]),geometry=s["geometry"],origin="NEW_NATIVE_N128",
                         sha=item["sha"],path=item["path"])
    table=[(n,record,"train") for n,record in train.items()]
    for n,record in seen.items():
        if not min(src_n)<n<max(src_n):
            excluded.append(dict(n=n,reason="OUTSIDE_FROZEN_TRAIN_SUPPORT"))
            continue
        near=np.searchsorted(src_n,n)
        if not(train[src_n[near-1]]["mu"]<record["mu"]<train[src_n[near]]["mu"]):
            excluded.append(dict(n=n,reason="NATIVE_CHEMICAL_POTENTIAL_NOT_BETWEEN_H_KNOTS",
                      left_train_mu=train[src_n[near-1]]["mu"],right_train_mu=train[src_n[near]]["mu"],
                      direct_native_mu=record["mu"]))
            continue
        table.append((n,record,"holdout"))
    table=sorted(table)
    data=dict(n=np.asarray([t[1]["n"] for t in table]),mu=np.asarray([t[1]["mu"] for t in table]),
        p=np.asarray([t[1]["p"] for t in table]),epsilon=np.asarray([t[1]["epsilon"] for t in table]))
    idx_train=np.array([i for i,x in enumerate(table) if x[2]=="train"],int)
    idx_hold=np.array([i for i,x in enumerate(table) if x[2]=="holdout"],int)
    if len(idx_train)<5 or np.any(np.diff(data["mu"])<=0) or np.any(np.diff(data["p"])<=0):
        raise RuntimeError("native source is not a consistent monotone single-branch reference")
    detail=dict(schema="ns_cowling.tf_native_holdout.v1",case=case,model=CASES[case]["model"],
        known_original_accepted=len(srcrows),original_targets=len(base["records"]),
        original_success=base["accepted_count"],refined_success=len(recovered),
        heldout_distinct=len(idx_hold),independent_repeat_at_knot=repeat,
        excluded_direct_states=excluded,
        original_rejected=[x["n_fm3"] for x in base["records"] if x["status"]=="REJECTED_NATIVE_STATE" and round(x["n_fm3"],9) not in recovered],
        training_shas={str(k):v["sha"] for k,v in train.items()},
        new_direct_shas={str(k):v["sha"] for k,v in seen.items()},model_qualifier="fixed droplet N128 single source branch, no minimum over geometry",
        test_density_node_count=len(data["n"]),
        benchmark_coordinate="at native holdout chemical potential, withheld n,P,epsilon",
        status="NO_NEW_WITHHELD_SUPPORT" if len(idx_hold)<2 else "RUNNING",
        metrics={})
    if len(idx_hold)>=2:
        truth=native_truth_cs2(data)
        detail["direct_native_fd_diagnostics"]=truth
        for method in GLOBAL_METHODS+SHAPE_METHODS:
            try:
                if method in GLOBAL_METHODS:pred=predict_global(method,data,idx_train,idx_hold);aux=None
                else:pred,aux=prediction_from_shape(method,data,idx_train,idx_hold)
                met=summarize(method,data,idx_hold,pred,truth,aux)
                met["status"]="EVALUATED_AGAINST_NEW_NATIVE_DROPLET_N128"
                detail["metrics"][method]=met
            except Exception as e:
                detail["metrics"][method]=dict(status="METHOD_FAILED",reason=str(e))
        detail["status"]="SCOPED_FINITE_NATIVE_VALIDATION"
    detail["wall_seconds"]=time.monotonic()-start_time
    OUTPUT.mkdir(parents=True,exist_ok=True)
    save(OUTPUT/f"{case}.json",detail)
    print("TF_HOLDOUT",case,"new",len(idx_hold),"train",len(idx_train),
          "repeat",len(repeat),"excluded",len(excluded),
          "status",detail["status"],flush=True)
    return detail
if __name__=="__main__":
    from concurrent.futures import ThreadPoolExecutor,as_completed
    cases=list(CASES)
    with ThreadPoolExecutor(max_workers=5) as p:
        future=[p.submit(run,case) for case in cases]
        rows=[]
        for f in as_completed(future):
            rows.append(f.result())
    out=dict(status="COMPLETED_ALL_FIVE",cases=[
        dict(case=x["case"],new=x["heldout_distinct"],repeat=len(x["independent_repeat_at_knot"]),
             status=x["status"]) for x in rows])
    save(OUTPUT/"SUMMARY.json",out)
    print(json.dumps(out,indent=2))
