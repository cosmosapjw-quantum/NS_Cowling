"""Native dense RMF reference generation for blind interpolation holdout.

Each model is an independent Fortran solver process, model parameters frozen.
No interpolation may be used as ground truth. Grid: 2049 linear nodes
on [0.015,0.16] and 2049 geometric on [0.16,1.0] fm^-3, overlap one.
"""
from __future__ import annotations
import argparse,hashlib,json,sys,time,traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.rmf_candidate import (PRIMARY_MODELS,DIAGNOSTIC_MODELS,load_parameters,
                                       evaluate_states,export_compose,residual_report,sha256)
MODELS=tuple(PRIMARY_MODELS)+tuple(DIAGNOSTIC_MODELS)
PARAM=ROOT/"data/new_model_maxwell_20261117/sources/PARAMETER_INPUTS.json"
OUT=ROOT/"data/native_dense_interpolation_20261118/rmf"
EXE=Path("/home/cosmosapjw/ns_cowling_wt/scratch/native_dense_interp_20261118/rmf_build/rmf_candidate_driver")
def grid():
    n=np.r_[np.linspace(.015,.16,2049),np.geomspace(.16,1.,2049)[1:]]
    assert len(n)==4097 and np.all(np.diff(n)>0)
    return n
def run(model):
    if model not in MODELS:raise ValueError(f"unsupported model {model}")
    out=OUT/model
    if out.exists():raise FileExistsError(f"refusing to overwrite native evidence: {out}")
    out.mkdir(parents=True)
    n=grid()
    files={"parameters_sha256":sha256(PARAM),"executable_sha256":sha256(EXE),
           "native_driver_generator":sha256(ROOT/"python_frontend/compose_eos/rmf_candidate.py"),
           "density_sha256":hashlib.sha256(n.tobytes()).hexdigest()}
    receipt=dict(schema="ns_cowling.native_dense_rmf_reference.v1",model=model,
                 status="RUNNING",scope="direct native uniform beta-equilibrated state; finite range not global ground state",
                 grid_count=len(n),range_fm3=[float(n[0]),float(n[-1])],provenance=files)
    def save(): (out/"RESULT.json").write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    save()
    t=time.monotonic()
    try:
        params=load_parameters(PARAM,model)
        states=evaluate_states(EXE,params,n,work_dir=out/"native_raw")
        direct=residual_report(states)
        receipt["direct_native_residual"]=direct
        if not direct["physical_state_checks"]:raise ValueError("direct Fortran model physical residual failed")
        io=export_compose(out/"table",states,params)
        receipt.update(status="GENERATED_VERIFIED_DIRECT_NATIVE",
           elapsed_seconds=time.monotonic()-t,
           readback=io["native_readback_scaled_residuals"],composition_scope="homogeneous npe(mu) beta equilibrium",
           table_sha256=io["table_sha256"],
           native_state_csv_sha256=sha256(out/"native_raw/states.csv"),
           numerical_samplings=dict(holdout_all=4097,training_every_16=257,training_every_8=513,
                                    excluded_truth_points_for_stride16=3840))
        save()
        print("RMF_SUCCESS",model,4097,f"{receipt['elapsed_seconds']:.3f}",flush=True)
        return receipt
    except Exception as exc:
        receipt.update(status="FAILED",elapsed_seconds=time.monotonic()-t,
                       exception_type=type(exc).__name__,exception_text=str(exc),
                       exception_traceback=traceback.format_exc()[-4000:])
        save()
        print("RMF_FAILED",model,type(exc).__name__,str(exc),flush=True)
        return receipt
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--model",required=True,choices=MODELS)
    args=ap.parse_args()
    x=run(args.model)
    raise SystemExit(0 if x["status"]=="GENERATED_VERIFIED_DIRECT_NATIVE" else 2)
