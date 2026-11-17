#!/usr/bin/env python3
"""Prepare or execute the prospective final-API plan in a separate directory."""
from pathlib import Path
import argparse
import copy
import hashlib
import importlib.util
import json
import os

ROOT=Path(__file__).resolve().parents[3]
HISTORICAL=Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=["prepare","run"])
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--plan",type=Path,default=HISTORICAL/"PLAN_FINAL_API.json")
    args=parser.parse_args()
    plan_path=args.plan.resolve()
    plan=json.loads(plan_path.read_text())
    original=json.loads((HISTORICAL/"PLAN.json").read_text())
    prospective=plan["prospective_reproduction"]
    assert prospective["status"]=="PROSPECTIVE_NOT_EXECUTED"
    assert sha(HISTORICAL/"PLAN.json")==prospective["original_production_plan_sha256"]
    key="python_frontend/compose_eos/tf_candidate.py"
    comparison=copy.deepcopy(plan)
    comparison.pop("prospective_reproduction")
    comparison["source_sha256"][key]=original["source_sha256"][key]
    assert comparison==original,"Unexpected physical-grid/gate/source-plan change"
    assert plan["source_sha256"][key]==prospective["final_API_sha256"]
    for path,digest in plan["source_sha256"].items():
        assert sha(ROOT/path)==digest,"Frozen source hash mismatch: "+path
    for variable in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"]:
        os.environ[variable]="1"
    runner_path=ROOT/"scripts/repair_maxwell_support_20261117.py"
    spec=importlib.util.spec_from_file_location("support_reproduction_driver",runner_path)
    runner=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    assert runner.CONFIGS==plan["configurations"] and runner.MESHES==plan["meshes"]
    output=args.output.resolve()
    output.relative_to(ROOT)  # runner's source/output receipts use repository-relative paths
    assert output!=HISTORICAL,"Historical production directory cannot be used"
    if output.exists():
        assert output.is_dir() and {p.name for p in output.iterdir()}=={"PLAN.json"}, "Output must be new or contain only the same prospective PLAN.json"
        assert sha(output/"PLAN.json")==sha(plan_path),"Output plan differs"
    else:
        output.mkdir(parents=True)
        (output/"PLAN.json").write_bytes(plan_path.read_bytes())
    receipt={"status":"PROSPECTIVE_PLAN_PREPARED","action":args.action,
        "plan_path":str(plan_path),"plan_sha256":sha(plan_path),
        "output":str(output),"historical_plan_sha256":sha(HISTORICAL/"PLAN.json"),
        "launcher_sha256":sha(__file__),"source_hashes_verified":plan["source_sha256"],
        "physical_plan_unchanged":True,"native_calls":0,
        "execution_evidence":"Preparation only; no native calculations have run in this launcher yet."}
    print(json.dumps(receipt,ensure_ascii=False),flush=True)
    if args.action=="prepare":
        return
    # Retain the original production runner byte-for-byte; only its output BASE changes.
    runner.BASE=output
    original_builder=runner.build_driver
    def checked_builder(repo,build_dir):
        executable=original_builder(repo,build_dir)
        assert sha(executable)==plan["executable_sha256"],"Rebuilt native executable identity differs from the declared plan"
        return executable
    runner.build_driver=checked_builder
    runner.regenerate()
    runner.match()

if __name__=="__main__":
    main()
