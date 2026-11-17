#!/usr/bin/env python3
"""Actual native-generated EOS to public matching integration, not a new solver.

The core interval ends at 0.5 fm^-3 deliberately: this bounded matching campaign
does not repair known high-density model acausality or claim complete stellar EOS.
Python performs transport, hashing and serialization comparisons only.
"""
from __future__ import annotations
import argparse, hashlib, json, os, subprocess, sys, time
from collections import Counter
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "python_frontend"), str(ROOT / "scripts")]
from compose_eos import read_compose, write_compose, match_crust_core
from compose_eos.generation_dispatch import CORE_IDS, SPG_IDS
from compose_eos.native_generation import verify_native_generation
from match_core_crust import run_pipeline
from match_report_admission import classify_failure

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False, default=str) + "\n")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--build", type=Path, required=True)
    p.add_argument("--spg-root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--reuse-core-root", type=Path)
    args = p.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    env = {**os.environ, "OMP_NUM_THREADS":"1", "OPENBLAS_NUM_THREADS":"1",
           "MKL_NUM_THREADS":"1", "PYTHONPATH":str(ROOT/"python_frontend")+":"+str(ROOT/"scripts"),
           "COMPOSE_EOS_FORTRAN_ROOT":str(ROOT),
           "COMPOSE_EOS_FORTRAN_EXE":str(args.build.resolve()/"compose_eos_fortran"),
           "COMPOSE_EOS_NATIVE_BUILD_DIR":str(args.build.resolve())}
    os.environ.update(env)
    deps = sorted(list((ROOT/"fortran_eos/src").glob("*.f90")) +
                  list((ROOT/"fortran_eos/app").glob("*.f90")) +
                  list((ROOT/"python_frontend/compose_eos").glob("*.py")) +
                  [ROOT/"scripts/match_core_crust.py", ROOT/"scripts/match_report_admission.py",
                   Path(__file__)])
    initial = {str(x.relative_to(ROOT)):sha(x) for x in deps}
    manifest = {"baseline":"478c7b5de1514ac512f91156f36551478f36f688",
                "started_unix":time.time(), "core_n_max_fm3":0.5,
                "scope":"finite regenerated tables, conditional inferred SPG; thermo-only matching",
                "production_physics":"Fortran", "workers":1, "threads":1,
                "source_sha256_before":initial,
                "binary_sha256":{p:sha(args.build/p) for p in ("compose_eos_fortran","compose_eos_generate")},
                "core_generation":[], "pairs":[]}
    save(out/"RECEIPT.json",manifest)
    cores = {}
    for model in CORE_IDS:
        if args.reuse_core_root is not None:
            directory = args.reuse_core_root.resolve()/model
            receipt_path = directory/"DISPATCH.json"
            receipt = json.loads(receipt_path.read_text())
            evidence = verify_native_generation(directory/"generated",receipt["request"])
            manifest["core_generation"].append({"model":model,"mode":"REUSED_NATIVE_VERIFIED",
                "source_dispatch":str(receipt_path),"dispatch_sha256":sha(receipt_path),
                "output_sha256":evidence["output_sha256"],"rows":evidence["rows"]})
            cores[model] = directory/"generated"
            continue
        directory = out/"core"/model
        command = [sys.executable,"-m","compose_eos.cli","generate",model,str(directory),
                   "--points","129","--n-max","0.5"]
        with (out/(model+".log")).open("w") as log:
            r = subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        entry = {"model":model,"command":command,"returncode":r.returncode}
        manifest["core_generation"].append(entry)
        save(out/"RECEIPT.json",manifest)
        if r.returncode:
            raise RuntimeError("native public core generation failed: "+model)
        cores[model] = directory / "generated"
    for crust in SPG_IDS:
        source = args.spg_root.resolve()/crust/"generated"
        if not (source/"RESULT.json").is_file():
            raise RuntimeError("missing successful SPG source: "+str(source))
        cr = read_compose(source)
        for core, directory in cores.items():
            co = read_compose(directory)
            for method in ("css","cubic_F","cubic_P","quintic_F","quintic_P","maxwell"):
                dest = out/"matches"/(crust+"__"+core+"__"+method)
                item = {"crust":crust,"core":core,"method":method,
                        "source_sha256":{side:{f.name:sha(f) for f in d.glob("*")
                            if f.is_file() and (f.name.startswith("eos.") or f.name=="SUPPORT.json")}
                            for side,d in (("crust",source),("core",directory))}}
                try:
                    if method == "maxwell":
                        result = run_pipeline(method,source,directory,dest)
                        item["pipeline"] = result
                        item["outcome"] = "ADMITTED" if result.get("survived") else result.get("outcome","NOT_ADMITTED")
                        item["failure_category"] = result.get("failure_category")
                    else:
                        matched, info = match_crust_core(cr,co,cr_name=crust,co_name=core,
                            method=method,n_bridge=40,drop_composition=True)
                        write_compose(dest/"compose",matched)
                        back = read_compose(dest/"compose")
                        comparisons = {}
                        for column in ("nb","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Iphase"):
                            a,b = np.asarray(getattr(matched,column)),np.asarray(getattr(back,column))
                            ok = a.shape==b.shape and bool(np.allclose(a,b,rtol=5e-13,atol=1e-14,equal_nan=False))
                            comparisons[column] = ok
                        if not all(comparisons.values()):
                            raise RuntimeError("native serialized matching readback mismatch: "+str(comparisons))
                        item.update(report=info,readback=comparisons,rows=matched.n,
                                    outcome="CONSTRUCTED_DIAGNOSTIC_ONLY" if info["diagnostic_only"] else "CONSTRUCTED_VERIFIED")
                except Exception as exc:
                    item.update(outcome="REJECTED",error_type=type(exc).__name__,error=str(exc),
                                failure_category=classify_failure("matching",str(exc),None))
                save(dest/"PAIR.json",item)
                manifest["pairs"].append({"crust":crust,"core":core,"method":method,
                    "outcome":item["outcome"],"failure_category":item.get("failure_category"),
                    "receipt":str((dest/"PAIR.json").relative_to(out))})
                save(out/"RECEIPT.json",manifest)
                print(crust,core,method,item["outcome"],flush=True)
    final = {str(x.relative_to(ROOT)):sha(x) for x in deps}
    manifest.update(finished_unix=time.time(),source_sha256_after=final,
                    sources_unchanged=initial==final,
                    outcomes=dict(Counter(x["outcome"] for x in manifest["pairs"])),
                    by_method={m:dict(Counter(x["outcome"] for x in manifest["pairs"] if x["method"]==m))
                        for m in ("css","cubic_F","cubic_P","quintic_F","quintic_P","maxwell")})
    manifest["execution_status"] = "COMPLETE" if initial==final and len(manifest["pairs"])==150 else "IDENTITY_CHANGED"
    save(out/"RECEIPT.json",manifest)
    return 0 if manifest["execution_status"]=="COMPLETE" else 1
if __name__=="__main__":
    raise SystemExit(main())
