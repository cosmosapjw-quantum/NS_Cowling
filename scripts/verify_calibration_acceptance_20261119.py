#!/usr/bin/env python3
"""Run fresh scoped native TOV/Cowling tests, compiler checks and source hash gate.

No test is silently reused from previous commit. The receipt binds the
results to exact Python/Fortran source bytes and to the running executable.
"""
from __future__ import annotations
import hashlib,json,os,re,subprocess,sys,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/"data/cowling_eos_source_calibration_20261119/verification"
TESTS=[
"test_cowling_mode_continuation_20261119.py",
"test_cowling_sharp_knot_waypoints_20261119.py",
"test_cowling_source_calibration_20261119.py",
"test_cowling_calibration_evidence_20261119.py",
"test_cowling_maxwell_native_20261119.py",
"test_cowling_maxwell_sources_20261119.py",
"test_native_cowling_tov_20261119.py",
"test_native_polytrope_cowling_20261119.py",
"test_native_cowling_dvode_20261118.py",
"test_stellar_background.py",
"test_stellar_cowling.py",
"test_cowling_reaudit_contracts.py",
"test_css_cowling_adapter.py",
"test_cowling_end_to_end.py",
"test_maxwell_star.py",
"test_native_maxwell_transition.py",
]
SOURCES=[
"fortran_eos/src/compose_eos_cowling_dvode.f90",
"fortran_eos/src/compose_eos_cowling_tov_dvode.f90",
"python_frontend/compose_eos/stellar_sharp_native.py",
"python_frontend/compose_eos/stellar_cowling_native.py",
"python_frontend/compose_eos/source_calibration.py",
"python_frontend/compose_eos/stellar_mode_continuation.py",
]+[f"python_frontend/tests/{v}" for v in TESTS]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def call(command,log):
    run=subprocess.run(command,cwd=ROOT,stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT,text=True,env=ENV)
    log.write_text(run.stdout)
    print("RUN",command[0],"exit",run.returncode,"stdout_last",
          run.stdout.splitlines()[-1:] if run.stdout else [],flush=True)
    if run.returncode:
        raise RuntimeError(f"{command[0]} failed: {run.stdout[-1600:]}")
    return run

ENV=dict(os.environ)
ENV["OPENBLAS_NUM_THREADS"]="1"
ENV["OMP_NUM_THREADS"]="1"
ENV["PYTHONPATH"]="python_frontend"

def main():
    DEST.mkdir(exist_ok=True,parents=True)
    start=time.monotonic()
    sha_before={name:sha(ROOT/name) for name in SOURCES}
    check=subprocess.run(["git","diff","--check"],cwd=ROOT,
            text=True,capture_output=True)
    if check.returncode:raise RuntimeError("git diff --check: "+check.stderr)
    pyfiles=[v for v in SOURCES if v.endswith(".py")]
    call([sys.executable,"-m","compileall","-q",*pyfiles],DEST/"PYTHON_COMPILEALL.log")
    build=[]
    for name in SOURCES[:2]:
        output=DEST/(Path(name).stem+".o")
        try:
            cmd=["gfortran","-c","-O0","-Wall","-Wextra","-fcheck=all",
                 "-ffree-line-length-none",name,"-o",str(output)]
            r=call(cmd,DEST/(Path(name).stem+".build.log"))
            build.append(dict(source=name,exit_status=r.returncode,
                       compiler_output_sha256=sha(DEST/(Path(name).stem+".build.log"))))
        finally:
            output.unlink(missing_ok=True)
    test_command=[sys.executable,"-m","pytest","-q",*[
        f"python_frontend/tests/{name}" for name in TESTS],"--disable-warnings"]
    tested=call(test_command,DEST/"TEST_RUN.log")
    pattern=re.search(r"(\d+) passed, (\d+) skipped",tested.stdout)
    if not pattern:raise RuntimeError("scoped test result not summarized with pass/skip")
    passing,skipping=map(int,pattern.groups())
    if passing<120:
        raise RuntimeError("fewer observed passing tests than intended scoped regression")
    sha_after={name:sha(ROOT/name) for name in SOURCES}
    if sha_before!=sha_after:raise RuntimeError("tested source changed during execution")
    receipts=dict(schema="ns_cowling.source_calibration_acceptance.v1",
        parent_03_commit="0d802d09fff57df1cc5557916be824cbc457d2b7",
        git_head_during_verification=subprocess.check_output(
            ["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        exact_tested_source_sha256=sha_before,
        runtime_python=sys.version,
        runtime_gfortran=subprocess.check_output(["gfortran","--version"],
            text=True).splitlines()[0],
        test_command=test_command,
        pytest=dict(exit_status=tested.returncode,
             passed=passing,skipped=skipping,
             test_log_sha256=sha(DEST/"TEST_RUN.log"),
             test_log_relative_path=str((DEST/"TEST_RUN.log").relative_to(ROOT))),
        Fortran_compile_with_bounds_checks=build,
        python_compileall_passed=True,git_diff_check_passed=True,
        elapsed_wall_seconds=time.monotonic()-start,
        scope="stage04 + M03 complete relevant F90 Cowling/TOV regression, no global EOS proof",
        status="SCOPED_IMPLEMENTATION_VERIFIED")
    (DEST/"VERIFICATION_RECEIPT.json").write_text(json.dumps(receipts,indent=2)+"\n")
    print("FRESH_STAGE_04_ACCEPTANCE",passing,"passed",skipping,"skipped",
          "in",round(receipts["elapsed_wall_seconds"],2),"seconds",flush=True)

if __name__=="__main__":main()
