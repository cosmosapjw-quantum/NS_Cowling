#!/usr/bin/env python3
"""Fresh scoped verification of source-coefficient-aware sharp Cowling.

Legacy suite is intentionally run with its unchanged normal pytest marker
policy, while expensive reference-matched star tests run explicitly with
--run-slow. The source and evidence SHA-256 are frozen around those calls.
"""
from __future__ import annotations
import hashlib,json,os,re,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data/cowling_certified_modes_20261119"
OUT=DATA/"verification"
sys.path.insert(0,str(ROOT/"scripts"))
from verify_cowling_rmf_onepoint_20261119 import TESTS as LEGACY_TESTS

NEW_TESTS=["test_cowling_sharp_certified_20261119.py",
           "test_cowling_sharp_certified_evidence_20261119.py"]
SOURCE_FILES=[
    "python_frontend/compose_eos/stellar_sharp_certified.py",
    "python_frontend/tests/test_cowling_sharp_certified_20261119.py",
    "python_frontend/tests/test_cowling_sharp_certified_evidence_20261119.py",
    "scripts/run_certified_sharp_reaudit_20261119.py",
    "scripts/verify_certified_sharp_20261119.py",
    "python_frontend/compose_eos/stellar_sharp_native.py",
    "python_frontend/compose_eos/stellar_sharp.py",
    "python_frontend/compose_eos/stellar_mode_continuation.py",
    "fortran_eos/src/compose_eos_cowling_tov_dvode.f90",
    "fortran_eos/src/compose_eos_cowling_dvode.f90",
    "data/cowling_certified_modes_20261119/RESULT.json",
    "data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip",
    "docs/COWLING_SHARP_SOURCE_REFERENCE_GATE_20261119.md",
]
ENV=dict(os.environ,PYTHONPATH="python_frontend",OPENBLAS_NUM_THREADS="1",
         OMP_NUM_THREADS="1")

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run(command,log):
    p=subprocess.run(command,cwd=ROOT,env=ENV,
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    log.write_text(p.stdout)
    if p.returncode:
        raise RuntimeError("actual command failed "+str(command)+": "+p.stdout[-2500:])
    print("VERIFIED_COMMAND",command[0],"exit",p.returncode,
          "last",p.stdout.splitlines()[-1:] if p.stdout else [],flush=True)
    return p.stdout

def counts(raw):
    p=re.search(r"(\d+) passed, (\d+) skipped",raw)
    if p:
        return list(map(int,p.groups()))
    p=re.search(r"(\d+) passed in",raw)
    if p:
        return [int(p.group(1)),0]
    raise ValueError("pytest output not a verified passing test summary")

def main():
    OUT.mkdir(exist_ok=True,parents=True)
    since=time.monotonic()
    old_sha={name:sha(ROOT/name) for name in SOURCE_FILES}
    pyfiles=[v for v in SOURCE_FILES if v.endswith(".py")]
    run([sys.executable,"-m","compileall","-q",*pyfiles],
        OUT/"PYTHON_COMPILEALL.log")
    for src in SOURCE_FILES:
        if src.endswith(".f90"):
            obj=OUT/(Path(src).stem+".o")
            try:
                run(["gfortran","-c","-O0","-fcheck=all","-Wall","-Wextra",
                     "-ffree-line-length-none",src,"-o",str(obj)],
                    OUT/(Path(src).stem+".build.log"))
            finally:
                obj.unlink(missing_ok=True)
    if subprocess.run(["git","diff","--check"],cwd=ROOT).returncode:
        raise ValueError("git diff --check rejected source edits")
    # The older native matcher is a separate executable. An isolated Git
    # worktree does not inherit old build/ files; compile it from the
    # exact currently pinned Fortran sources before running legacy tests.
    build_result=run(["bash","fortran_eos/build.sh","--ensure"],
                     OUT/"BUILD_NATIVE_MATCHER.log")
    native_exe=ROOT/"build/compose_eos_fortran/compose_eos_fortran"
    if not native_exe.is_file():
        raise RuntimeError("native Fortran matcher missing after source build")
    baseline=[f"python_frontend/tests/{name}" for name in LEGACY_TESTS]
    a=run([sys.executable,"-m","pytest","-q",*baseline,"--disable-warnings"],
          OUT/"BASELINE_TEST_RUN.log")
    b=run([sys.executable,"-m","pytest","-q",
           *[f"python_frontend/tests/{name}" for name in NEW_TESTS],
           "--run-slow","--disable-warnings"],OUT/"NEW_TEST_RUN.log")
    pa,sa=counts(a)
    pb,sb=counts(b)
    if pa<172 or sa!=3 or pb<11 or sb!=0:
        raise RuntimeError("scientific regression count lower than required")
    new_sha={name:sha(ROOT/name) for name in SOURCE_FILES}
    if old_sha!=new_sha:
        raise RuntimeError("source/derived evidence changed during full verification")
    out=dict(schema="ns_cowling.sharp_certified_source_accuracy.v1",
        original_science_parent="1e7823ecd62a2914d40c0948d3811410b1ea16da",
        input_source_evidence_sha256=old_sha,
        baseline=dict(passed=pa,skipped=sa,
                 log_sha256=sha(OUT/"BASELINE_TEST_RUN.log")),
        upgrade=dict(passed=pb,skipped=sb,
                 log_sha256=sha(OUT/"NEW_TEST_RUN.log"),
                 actual_real_star_tests_executed=True),
        total_passed=pa+pb,total_skipped=sa+sb,
        python_compileall_passed=True,
        fortran_cowling_tov_bounds_checked_compilation=True,
        native_matcher_compiled_fresh_in_isolated_worktree=True,
        native_matcher_binary_sha256=sha(native_exe),
        native_matcher_builder_log_sha256=sha(OUT/"BUILD_NATIVE_MATCHER.log"),
        source_identity_unchanged_over_test=True,
        verified_numerical_results_sha256=sha(DATA/"RESULT.json"),
        runtime_seconds=time.monotonic()-since,
        status="SCOPED_CERTIFIED_SHARP_IMPLEMENTATION_VERIFIED",
        no_global_p1_family_or_full96_star_claim=True)
    (OUT/"VERIFICATION_RECEIPT.json").write_text(json.dumps(out,indent=2)+"\n")
    print("FULL_SHARP_CERTIFIED_PASS",pa+pb,"passed",sa+sb,"skipped",
         "run_seconds",out["runtime_seconds"],flush=True)

if __name__=="__main__":
    main()
