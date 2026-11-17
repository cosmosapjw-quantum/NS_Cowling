#!/usr/bin/env python3
"""Exercise the public native generator for the complete declared model cohort.

This script orchestrates CLI requests, hashes files and checks native readback.
It performs no nuclear physics, phase selection, interpolation or table repair.
Each retry receives a new directory; prior failures are never overwritten.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
CORE=("APR_APR","BBB_BHF2","GM_GM1","RG_KDE0v1","RG_SLY4")
SPG=tuple(f"SPG_M{i}" for i in range(1,6))
TF=("GPPVA_DD2","GPPVA_DDME2","GPPVA_TW","GPPVA_FSU2","GPPVA_FSU2R",
    "GPPVA_FSU2H","GPPVA_NL3wrL55","GPPVA_TM1e")
MODELS=CORE+SPG+TF
NATIVE_FILES=("eos.nb","eos.t","eos.yq","eos.thermo","eos.compo")

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save(path,value):
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")
    temporary.replace(path)

def utc():
    return datetime.now(timezone.utc).isoformat()

def code_files():
    native=ROOT/"fortran_eos"
    paths=[native/"build.sh",native/"run-generator.sh",Path(__file__)]
    paths+=list((native/"src").glob("*.f90"))+list((native/"app").glob("*.f90"))
    paths+=list((native/"vendor").rglob("*.f"))
    paths+=[ROOT/"python_frontend/compose_eos"/name for name in
            ("cli.py","generation_dispatch.py","native_generation.py","native_model_contracts.py",
             "api.py","native.py","eos_data.py")]
    return sorted(set(paths))

def identities(paths):
    return {str(p.relative_to(ROOT)):sha(p) for p in paths}

def attempt_directory(out,model):
    parent=out/model
    parent.mkdir(parents=True,exist_ok=True)
    for i in range(1,10000):
        attempt=parent/f"attempt_{i:03d}"
        try:
            attempt.mkdir()
            return attempt
        except FileExistsError:
            pass
    raise RuntimeError("too many preserved attempts")

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"data/research_native_20261116/public_cohort")
    parser.add_argument("--workers",type=int,default=4)
    parser.add_argument("--models",nargs="+",choices=MODELS,default=MODELS)
    parser.add_argument("--timeout-seconds",type=int,default=3600)
    args=parser.parse_args(argv)
    if not 1<=args.workers<=4:
        parser.error("workers must be between 1 and 4")
    if len(set(args.models))!=len(args.models):
        parser.error("duplicate models are not accepted")
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    paths=code_files();frozen=identities(paths)
    env=dict(os.environ,PYTHONPATH=str(ROOT/"python_frontend"),
             OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",NUMEXPR_NUM_THREADS="1")
    command_prefix=[sys.executable,"-m","compose_eos.cli"]
    commit=subprocess.run(["git","rev-parse","HEAD"],cwd=ROOT,capture_output=True,text=True,check=True).stdout.strip()
    summary={"status":"RUNNING","started_utc":utc(),"baseline_commit":commit,
             "backend":"fortran","requested_models":list(args.models),"concurrency":args.workers,
             "threads_per_child":1,"source_sha256":frozen,"models":{},
             "scope":"finite native public-generation cohort; source fidelity, physical support, and CompOSE comparison are separate",
             "core_grid_points":161,"SPG_grid_points":17,"SPG_variant":"AME2408_candidate",
             "TF_densities_fm3":[.003,.03,.075],"TF_meshes":[16],"crust_lower_segments":True,
             "source_equivalence":"SPG surface calibration remains candidate; BSk is analytic; BBB is published-energy reconstruction"}
    summary_path=out/"SUMMARY.json"
    # Preserve earlier aggregate receipts too.
    if summary_path.exists():
        archived=out/("SUMMARY_before_"+str(time.time_ns())+".json")
        archived.write_bytes(summary_path.read_bytes())
    save(summary_path,summary)

    def execute(model):
        attempt=attempt_directory(out,model)
        generated=attempt/"output"
        options=["--points","161"] if model in CORE else (
            ["--points","17","--variant","AME2408_candidate","--include-lower-crust"] if model in SPG else
            ["--densities","0.003","0.03","0.075","--resolutions","16","--include-lower-crust"])
        command=[*command_prefix,"generate",model,str(generated),"--repository-root",str(ROOT),*options]
        result={"model":model,"attempt_directory":str(attempt.relative_to(ROOT)),
                "command":command,"started_utc":utc(),"status":"RUNNING"}
        start=time.monotonic()
        try:
            with (attempt/"public_cli.log").open("w") as log:
                proc=subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,
                                    timeout=args.timeout_seconds,check=False)
            result["returncode"]=proc.returncode
            dispatch_path=generated/"DISPATCH.json"
            if dispatch_path.is_file():
                dispatch=json.loads(dispatch_path.read_text())
                result["dispatch_status"]=dispatch.get("status")
                result["dispatch_sha256"]=sha(dispatch_path)
                result["request"]=dispatch.get("request")
                result["native_evidence"]=dispatch.get("evidence")
            else:
                dispatch={}
            if proc.returncode!=0 or dispatch.get("status")!="REQUESTED_ROWS_NATIVE_VERIFIED":
                result["status"]="PUBLIC_GENERATION_FAILED"
            else:
                directories=[generated/"generated"]
                if model not in CORE:
                    directories += [directories[0]/"outer",directories[0]/"lower_source_fit"]
                readbacks=[]
                for directory in directories:
                    read_command=[*command_prefix,"read",str(directory)]
                    check=subprocess.run(read_command,cwd=ROOT,env=env,capture_output=True,
                                         text=True,timeout=120,check=False)
                    item={"directory":str(directory.relative_to(ROOT)),"returncode":check.returncode}
                    log_name="read_"+("inner" if directory==directories[0] else directory.name)+".log"
                    (attempt/log_name).write_text(check.stdout+check.stderr)
                    if check.returncode:
                        raise ValueError("public native readback failed: "+str(directory))
                    parsed=json.loads(check.stdout)
                    if parsed.get("status")!="ok":
                        raise ValueError("native readback did not report ok")
                    item["result"]=parsed
                    item["output_sha256"]={name:sha(directory/name) for name in NATIVE_FILES}
                    readbacks.append(item)
                result["readbacks"]=readbacks
                result["status"]="REQUESTED_ROWS_NATIVE_VERIFIED"
        except Exception as exc:
            result.update(status="EXECUTION_OR_READBACK_FAILED",error_type=type(exc).__name__,error=str(exc))
        result["elapsed_seconds"]=time.monotonic()-start
        result["finished_utc"]=utc()
        result["preserved_files"]={str(p.relative_to(attempt)):{"size":p.stat().st_size,"sha256":sha(p)}
                                   for p in sorted(attempt.rglob("*")) if p.is_file()}
        save(attempt/"COHORT_RESULT.json",result)
        return result

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        future={pool.submit(execute,model):model for model in args.models}
        for job in as_completed(future):
            model=future[job]
            try:
                item=job.result()
            except Exception as exc:
                item={"model":model,"status":"ORCHESTRATOR_FAILED","error_type":type(exc).__name__,"error":str(exc)}
            summary["models"][model]=item
            summary["finished_models"]=len(summary["models"])
            summary["verified_models"]=sum(v["status"]=="REQUESTED_ROWS_NATIVE_VERIFIED" for v in summary["models"].values())
            save(summary_path,summary)
            print(json.dumps({"model":model,"status":item["status"],
                              "elapsed_seconds":item.get("elapsed_seconds"),
                              "completed":len(summary["models"]),"total":len(args.models)}),flush=True)
    current=identities(paths)
    summary["source_unchanged"]=current==frozen
    if current!=frozen:
        summary["changed_source_files"]=[key for key in frozen if frozen[key]!=current.get(key)]
    summary["finished_utc"]=utc()
    summary["verified_readback_tables"]=sum(len(item.get("readbacks",[])) for item in summary["models"].values())
    summary["main_rows"]=sum(item["readbacks"][0]["result"]["n"]
                             for item in summary["models"].values() if item.get("readbacks"))
    summary["all_segment_rows"]=sum(read["result"]["n"] for item in summary["models"].values()
                                   for read in item.get("readbacks",[]))
    success=summary["verified_models"]==len(args.models) and summary["source_unchanged"]
    summary["status"]="FINITE_COHORT_NATIVE_VERIFIED" if success else "INCOMPLETE_OR_SOURCE_CHANGED"
    save(summary_path,summary)
    print(json.dumps({"status":summary["status"],"verified_models":summary["verified_models"],
                      "requested_models":len(args.models),"readback_tables":summary["verified_readback_tables"],
                      "main_rows":summary["main_rows"],"all_segment_rows":summary["all_segment_rows"]}),flush=True)
    return 0 if success else 1

if __name__=="__main__":
    raise SystemExit(main())
