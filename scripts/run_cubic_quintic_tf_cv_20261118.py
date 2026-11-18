"""Quintic/cubic sparse crust TF holdouts, same immutable H source/30 new native nodes.

No TF rejected trial is ever used as truth, no cross-geometry interpolation,
no newly generated model parameters and no original dataset mutation.
"""
from pathlib import Path
import json,sys,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from benchmark_tf_native_holdout_20261118 import collect
from benchmark_native_dense_eos_20261118 import summarize,native_truth_cs2
from compare_hermite_orders_20261118 import predict
from run_dense_tf_20261118 import CASES
OUT=ROOT/"data/cubic_quintic_matching_20261118/tf_cv"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def one(case):
    anchors,newvals,base,recovered=collect(case)
    source=ROOT/"data/native_dense_interpolation_20261118/tf_holdout"/f"{case}.json"
    evidence=json.loads(source.read_text())
    # Only the identical H accepted profile files and entirely new state hashes
    train_by_sha={s["_source_sha256"]:s for s in anchors.values()}
    held_by_sha={entry["sha"]:entry["state"] for entry in newvals}
    records=[]
    for nstring,digest in evidence["training_shas"].items():
        row=train_by_sha[digest]
        records.append(dict(n=float(nstring),p=float(row["pressure"]),epsilon=float(row["epsilon"]),
             mu=float(row["mu_B"]),kind="train",state_sha256=digest))
    for nstring,digest in evidence["new_direct_shas"].items():
        row=held_by_sha[digest]
        if not row.get("accepted_stationary") or row.get("geometry")!="droplet" or int(row.get("points",0))!=128:
            raise ValueError("TF holdout is not native N128 droplet stationary")
        records.append(dict(n=float(nstring),p=float(row["pressure"]),epsilon=float(row["epsilon"]),
             mu=float(row["mu_B"]),kind="held",state_sha256=digest))
    records.sort(key=lambda r:r["n"])
    data=dict(n=np.array([i["n"] for i in records]),
              p=np.array([i["p"] for i in records]),
              epsilon=np.array([i["epsilon"] for i in records]),
              mu=np.array([i["mu"] for i in records]))
    if np.any(np.diff(data["n"])<=0) or np.any(np.diff(data["mu"])<=0):
        raise ValueError("source-and-holdout combined phase not monotone")
    tr=np.array([i for i,v in enumerate(records) if v["kind"]=="train"],int)
    hold=np.array([i for i,v in enumerate(records) if v["kind"]=="held"],int)
    assert len(hold)==evidence["heldout_distinct"]
    truth=native_truth_cs2(data)
    result=dict(case=case,source_holdout_sha256=sha(source),training_ids=[records[i]["state_sha256"] for i in tr],
                heldout_ids=[records[i]["state_sha256"] for i in hold],
                known_native_failures=evidence["original_rejected"],
                n_training=len(tr),n_holdout=len(hold),
                cs2_native_fd_max_cross_form=truth["max_abs_cross_form"],
                direct_native_reference="N128 same droplet branch; no interpolated optimizer failure as truth",
                status="RUNNING",methods={})
    for method in ("cubic","quintic"):
        try:
            z=predict(method,data,tr,hold)
            m=summarize(method,data,hold,z,truth)
            m["status"]="CHECKED_NATIVE_TF_HOLDOUT"
            m["max_abs_euler_identity"]=float(np.max(abs(z["epsilon"]+z["p"]-z["mu"]*z["n"])))
            result["methods"][method]=m
        except Exception as ex:
            result["methods"][method]=dict(status="BLOCKED",message=str(ex),exception_type=type(ex).__name__)
    result["status"]="COMPLETED_WITH_PRESERVED_REJECTIONS"
    save(OUT/f"{case}.json",result)
    print("TF_CV",case,result["n_holdout"],
          [(k,v["status"],v.get("number_density_relative",{}).get("rms")) for k,v in result["methods"].items()],flush=True)
    return result
def main():
    from concurrent.futures import ThreadPoolExecutor,as_completed
    with ThreadPoolExecutor(max_workers=5) as pool:
        result=[v.result() for v in as_completed([pool.submit(one,case) for case in CASES])]
    save(OUT/"CAMPAIGN.json",dict(status="COMPLETED",n_cases=len(result),
          n_native_heldout=sum(x["n_holdout"] for x in result),
          methods=["cubic","quintic"]))
    print("TF_CV_DONE",sum(x["n_holdout"] for x in result),len(result))
if __name__=="__main__":main()
