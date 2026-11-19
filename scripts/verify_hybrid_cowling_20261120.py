#!/usr/bin/env python3
"""Fresh regression gate for the optional original-source cs2 hybrid Cowling.

Existing default pchip and explicit direct_source stay unchanged.
Two nonoverlapping test groups: 172 legacy + 31 new/scoped real-star tests.
Checks actual command status and immutable code and scientific result identity.
"""
from __future__ import annotations
import hashlib,json,os,re,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from verify_direct_source_dvode_20261119 import (
    LEGACY_TESTS,CERT_TESTS,NEW_TESTS as DIRECT_TESTS)
NEW_TESTS=[
    "test_sharp_cowling_hybrid_cs2_20261120.py",
    "test_sharp_cowling_hybrid_evidence_20261120.py",
]
D=ROOT/"data/cowling_hybrid_geometry_20261120"
OUT=D/"verification"
ENV=dict(os.environ,PYTHONPATH="python_frontend",OPENBLAS_NUM_THREADS="1",
         OMP_NUM_THREADS="1")
SOURCE=[
    "python_frontend/compose_eos/stellar_sharp_native.py",
    "python_frontend/compose_eos/stellar_sharp_certified.py",
    "python_frontend/compose_eos/stellar_mode_continuation.py",
    "python_frontend/compose_eos/stellar_cowling_native.py",
    "python_frontend/tests/test_sharp_cowling_hybrid_cs2_20261120.py",
    "python_frontend/tests/test_sharp_cowling_hybrid_evidence_20261120.py",
    "fortran_eos/src/compose_eos_cowling_dvode.f90",
    "fortran_eos/src/compose_eos_cowling_tov_dvode.f90",
    "scripts/benchmark_cowling_hybrid_cs2_20261120.py",
    "scripts/summarize_hybrid_accuracy_20261120.py",
    "scripts/verify_hybrid_cowling_20261120.py",
    "docs/COWLING_HYBRID_GEOMETRY_EXACT_CS2_20261120.md",
    "data/cowling_hybrid_geometry_20261120/RESULT.json",
    "data/cowling_hybrid_geometry_20261120/DECISIONS.json",
    "data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip",
]
def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def command(cmd,log):
    v=subprocess.run(cmd,cwd=ROOT,env=ENV,text=True,
                     stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    log.write_text(v.stdout)
    if v.returncode:
        print("COMMAND_FAILURE",cmd[0],"return",v.returncode,
              v.stdout[-3000:],flush=True)
        raise RuntimeError("verification failure "+str(cmd))
    print("RUN_PASS",cmd[0],v.stdout.splitlines()[-1:] if v.stdout else [],flush=True)
    return v.stdout
def test_summary(txt):
    p=re.search(r"(\d+) passed, (\d+) skipped",txt)
    if p:return tuple(map(int,p.groups()))
    p=re.search(r"(\d+) passed in",txt)
    if p:return int(p.group(1)),0
    raise RuntimeError("passing pytest summary not observed")
def main():
    started=time.monotonic()
    OUT.mkdir(parents=True,exist_ok=True)
    before={v:sha(ROOT/v) for v in SOURCE}
    pyfiles=[x for x in SOURCE if x.endswith(".py")]
    command([sys.executable,"-m","compileall","-q",*pyfiles],
            OUT/"FINAL_COMPILEALL.log")
    compilelogs={}
    for src in [x for x in SOURCE if x.endswith(".f90")]:
        path=OUT/(Path(src).stem+".o")
        try:
            log=OUT/(Path(src).stem+".FINAL_BUILD.log")
            command(["gfortran","-c","-O0","-Wall","-Wextra","-fcheck=all",
                     "-ffree-line-length-none",src,"-o",str(path)],log)
            compilelogs[src]=sha(log)
        finally:
            path.unlink(missing_ok=True)
    command(["bash","fortran_eos/build.sh","--ensure"],
            OUT/"FINAL_NATIVE_MATCHER_BUILD.log")
    matcher=ROOT/"build/compose_eos_fortran/compose_eos_fortran"
    if not matcher.is_file():
        raise RuntimeError("original native matcher not compiled")
    a=command([sys.executable,"-m","pytest","-q",
           *(f"python_frontend/tests/{n}" for n in LEGACY_TESTS),
           "--disable-warnings"],OUT/"FINAL_LEGACY_TEST_RUN.log")
    b=command([sys.executable,"-m","pytest","-q",
           *(f"python_frontend/tests/{n}" for n in
             (CERT_TESTS+DIRECT_TESTS+NEW_TESTS)),
           "--run-slow","--disable-warnings"],
           OUT/"FINAL_SOURCE_HYBRID_TEST_RUN.log")
    pass_a,skip_a=test_summary(a)
    pass_b,skip_b=test_summary(b)
    if (pass_a,skip_a)!=(172,3) or (pass_b,skip_b)!=(31,0):
        raise RuntimeError("test scope/count drift")
    after={v:sha(ROOT/v) for v in SOURCE}
    if before!=after:
        raise RuntimeError("tested source or science receipt mutated during test")
    if subprocess.run(["git","diff","--check"],cwd=ROOT).returncode:
        raise RuntimeError("git staged/working whitespace check failed")
    derived=json.loads((D/"DECISIONS.json").read_text())
    observed=json.loads((D/"RESULT.json").read_text())
    if derived["input_benchmark_sha256"]!=sha(D/"RESULT.json"):
        raise RuntimeError("decision no longer pinned to exact benchmark result")
    if observed["frozen_EOS_bundle_sha256"]!=before[SOURCE[-1]]:
        raise RuntimeError("physical source-digest mismatch")
    for c in observed["cases"].values():
        for variant in c["backend_variants"].values():
            if variant.get("status")=="MATCHED_REFERENCE_PASS":
                if (variant["join_relative"]>1e-5 or
                    variant["nodes"]!=[1,1] or
                    abs(variant["reference_delta_hz"])>.05):
                    raise RuntimeError("scientific numerical acceptance invalid")
    summary=dict(schema="ns_cowling.cowling_hybrid_cs2_source_locked_verification.v1",
         scientific_parent="6253dc40eee188ea207450ae13a5d4588aa6d398",
         exact_source_and_science_input_sha256=before,
         native_matcher_sha256=sha(matcher),
         fortran_bounds_checked=compilelogs,
         original_EOS_source_unchanged=True,
         status="SCOPED_HYBRID_SOURCE_AND_IMPLEMENTATION_VERIFIED",
         suites={"legacy":dict(passed=pass_a,skipped=skip_a,
             log_sha256=sha(OUT/"FINAL_LEGACY_TEST_RUN.log")),
             "source_hybrid":dict(passed=pass_b,skipped=skip_b,
             log_sha256=sha(OUT/"FINAL_SOURCE_HYBRID_TEST_RUN.log"),
             actual_slow_and_rapid_stars=True)},
         unique_tests=dict(passed=pass_a+pass_b,skipped=skip_a+skip_b),
         benchmark_result_sha256=sha(D/"RESULT.json"),
         empirical_decision_sha256=sha(D/"DECISIONS.json"),
         python_compileall_passed=True,
         git_diff_check_passed=True,
         wall_seconds=time.monotonic()-started)
    (OUT/"VERIFICATION_RECEIPT.json").write_text(
        json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print("HYBRID_STAGE_VERIFIED",pass_a+pass_b,"passed",
          skip_a+skip_b,"skipped",round(summary["wall_seconds"],3),
          flush=True)
if __name__=="__main__":
    main()
