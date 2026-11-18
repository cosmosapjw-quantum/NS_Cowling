"""Paired sparse-node cubic/quintic Hermite validation on native dense models.

The 25,612 native states belong to 10 fixed EOS input IDs and include
published-table reconstructions (BBB), not global nuclear matter truth.
Phase 6|8|7 APR segments are never blended.
"""
from __future__ import annotations
from pathlib import Path
import json,hashlib,time,sys,os
import numpy as np
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from benchmark_native_dense_eos_20261118 import load_native,train_holdout_indices,phase_runs,native_truth_cs2,summarize
from compare_hermite_orders_20261118 import predict
OUT=ROOT/"data/cubic_quintic_matching_20261118/dense_cv_v2"
MODELS={"rmf":("NL3_PROVIDER_PRECISION","PK1_CANONICAL","TM1_CANONICAL","DDMEX_PROVIDER_ROUNDED_N0","DDMEX_ORIGINAL_TABLE_RHO0"),
        "core":("APR_APR","BBB_BHF2","GM_GM1","RG_KDE0v1","RG_SLY4")}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def json_numpy(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    raise TypeError(str(type(v)))
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,allow_nan=False,default=json_numpy)+"\n")
def run_one(family,model):
    folder=OUT/family/model
    if folder.exists():raise FileExistsError(folder)
    folder.mkdir(parents=True)
    src=load_native(family,model)
    x=dict(model=model,family=family,source_sha256=src["sha256"],
           baseline_native_model="direct Fortran fixed parameter/function/table provider",
           status="RUNNING",total_native_rows=len(src["mu"]),
           method_definitions={"cubic":"scipy cubic Hermite P(mu) with exact native n at retained knots",
                               "quintic":"piecewise quintic Hermite P(mu), exact native P,n and PCHIP one-sided dn/dmu from retained knots only"},
           runs=[])
    save(folder/"RESULT.json",x)
    for a,b in phase_runs(src["phase"]):
        d={k:(v[a:b] if isinstance(v,np.ndarray) else v) for k,v in src.items()}
        if len(d["mu"])<16:
            x["runs"].append(dict(phase=int(d["phase"][0]),status="TOO_FEW_ROWS"))
            continue
        truth=native_truth_cs2(d)
        for stride in ([16,8] if family=="rmf" else [8,4]):
            tr,hold=train_holdout_indices(len(d["mu"]),stride)
            rec=dict(phase=int(d["phase"][0]),native_row_slice=[a,b],
                     training_stride=stride,training_rows=len(tr),holdout_rows=len(hold),
                     training_indices=tr.tolist(),holdout_indices_sha256=hashlib.sha256(np.asarray(hold,np.int64).tobytes()).hexdigest(),
                     native_reference_diagnostics=truth,methods={})
            for name in ("cubic","quintic"):
                try:
                    pred=predict(name,d,tr,hold)
                    metrics=summarize(name,d,hold,pred,truth)
                    metrics.update(status="CHECKED_NATIVE_HOLDOUT",sampled_positive_n=bool(np.all(pred["n"]>0)),
                        no_dense_truth_leakage=True,
                        source_and_training_identity="same exact source P,n,mu training indices",
                        max_abs_euler_residual_mev_fm3=float(np.max(np.abs(pred["epsilon"]+pred["p"]-pred["mu"]*pred["n"]))))
                    rec["methods"][name]=metrics
                except Exception as exc:
                    rec["methods"][name]=dict(status="BLOCKED_METHOD",
                                             error_type=type(exc).__name__,message=str(exc))
            x["runs"].append(rec)
            save(folder/"RESULT.json",x)
            print("DENSE",family,model,rec["phase"],"stride",stride,
                  [(k,v["status"],v.get("number_density_relative",{}).get("rms")) for k,v in rec["methods"].items()],flush=True)
    x["status"]="COMPLETED_WITH_SCOPED_REJECTIONS"
    save(folder/"RESULT.json",x)
    return x
def main():
    jobs=[(f,m) for f,mm in MODELS.items() for m in mm]
    start=time.monotonic()
    rows=[]
    with ThreadPoolExecutor(max_workers=5) as pool:
        for fut in as_completed([pool.submit(run_one,*arg) for arg in jobs]):
            z=fut.result()
            rows.append(dict(family=z["family"],model=z["model"],runs=len(z["runs"]),status=z["status"]))
    result=dict(status="FINISHED",parallel_workers=5,inputs=rows,number_models=len(rows),
                elapsed_seconds=time.monotonic()-start)
    save(OUT/"PARALLEL_SUMMARY.json",result)
    print("ALL_FINISHED",len(rows),len(jobs),result["elapsed_seconds"])
if __name__=="__main__":main()
