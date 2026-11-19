#!/usr/bin/env python3
"""Fresh source-locked verification of new direct-source DVODE backend.

Three separate scopes: previous 172-pass native EOS suite, 11 prior
certified-source physical tests, and 11 new direct-source tests including
slow/rapid numerical results and sensitivity under original source knots.
No existing source or historical evidence receipt is overwritten.
"""
from __future__ import annotations
import hashlib,json,os,re,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from verify_certified_sharp_20261119 import LEGACY_TESTS,NEW_TESTS as CERT_TESTS

DATA=ROOT/"data/cowling_direct_source_callback_20261119"
OUT=DATA/"verification"
NEW_TESTS=[
    "test_sharp_cowling_exact_source_callback_20261119.py",
    "test_sharp_cowling_direct_source_evidence_20261119.py",
]
SOURCES=[
    "fortran_eos/src/compose_eos_cowling_dvode.f90",
    "fortran_eos/src/compose_eos_cowling_tov_dvode.f90",
    "python_frontend/compose_eos/stellar_cowling_native.py",
    "python_frontend/compose_eos/stellar_sharp_native.py",
    "python_frontend/tests/test_sharp_cowling_exact_source_callback_20261119.py",
    "python_frontend/tests/test_sharp_cowling_direct_source_evidence_20261119.py",
    "scripts/benchmark_sharp_direct_source_callback_20261119.py",
    "scripts/verify_direct_source_rapid_20261119.py",
    "scripts/check_direct_source_tolerance_20261119.py",
    "scripts/audit_sharp_sound_speed_knots_20261119.py",
    "scripts/verify_direct_source_dvode_20261119.py",
    "docs/COWLING_NATIVE_DIRECT_SOURCE_DVODE_20261119.md",
    "data/cowling_direct_source_callback_20261119/RESULT.json",
    "data/cowling_direct_source_callback_20261119/RESULT_BEFORE_RESTART.json",
    "data/cowling_direct_source_callback_20261119/RAPID.json",
    "data/cowling_direct_source_callback_20261119/RAPID_BEFORE_RESTART.json",
    "data/cowling_direct_source_callback_20261119/TOLERANCE.json",
    "data/cowling_direct_source_callback_20261119/TOLERANCE_BEFORE_RESTART.json",
    "data/cowling_direct_source_callback_20261119/KNOT_AUDIT.json",
    "data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip",
    "data/cowling_certified_modes_20261119/RESULT.json",
]
ENV=dict(os.environ,PYTHONPATH="python_frontend",OPENBLAS_NUM_THREADS="1",
         OMP_NUM_THREADS="1")
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run(args,log):
    result=subprocess.run(args,cwd=ROOT,env=ENV,
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    log.write_text(result.stdout)
    print("SCOPE_VERIFY",args[0],"exit",result.returncode,
         "tail",result.stdout.splitlines()[-2:] if result.stdout else [],flush=True)
    if result.returncode:
        raise RuntimeError("failed source-locked command: "+
                           str(args)+": "+result.stdout[-1900:])
    return result.stdout

def testcounts(raw):
    match=re.search(r"(\d+) passed, (\d+) skipped",raw)
    if match:return [int(match.group(1)),int(match.group(2))]
    match=re.search(r"(\d+) passed in",raw)
    if match:return [int(match.group(1)),0]
    raise ValueError("pytest output missing passing summary")

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    before={name:sha(ROOT/name) for name in SOURCES}
    if subprocess.run(["git","diff","--check"],cwd=ROOT).returncode:
        raise ValueError("git diff whitespace conflict")
    pyfiles=[name for name in SOURCES if name.endswith(".py")]
    run([sys.executable,"-m","compileall","-q",*pyfiles],
         OUT/"PYTHON_COMPILEALL.log")
    builds=[]
    for src in SOURCES[:2]:
        obj=OUT/(Path(src).stem+".o")
        try:
            run(["gfortran","-c","-O0","-fcheck=all","-Wall","-Wextra",
                 "-ffree-line-length-none",src,"-o",str(obj)],
                 OUT/(Path(src).stem+".build.log"))
            builds.append(dict(source=src,bounds_checked=True,
                compiler_log_sha256=sha(OUT/(Path(src).stem+".build.log"))))
        finally:
            obj.unlink(missing_ok=True)
    run(["bash","fortran_eos/build.sh","--ensure"],
        OUT/"NATIVE_MATCHER_BUILD.log")
    exe=ROOT/"build/compose_eos_fortran/compose_eos_fortran"
    if not exe.is_file():
        raise RuntimeError("original native matcher missing after source build")
    jobs=[
        ("previous_native",LEGACY_TESTS,False,172,3),
        ("previous_certified",CERT_TESTS,True,11,0),
        ("direct_source",NEW_TESTS,True,11,0),
    ]
    receipts={}
    for name,tests,slow,expected_pass,expected_skip in jobs:
        command=[sys.executable,"-m","pytest","-q",
                 *[f"python_frontend/tests/{p}" for p in tests],
                 "--disable-warnings"]
        if slow:command.append("--run-slow")
        raw=run(command,OUT/(name.upper()+"_TEST_RUN.log"))
        observed=testcounts(raw)
        if observed!=[expected_pass,expected_skip]:
            raise RuntimeError("unanticipated test coverage in "+name+
                               ": "+repr(observed))
        receipts[name]=dict(passed=observed[0],skipped=observed[1],
             log_sha256=sha(OUT/(name.upper()+"_TEST_RUN.log")),command=command)
    after={name:sha(ROOT/name) for name in SOURCES}
    if before!=after:
        raise ValueError("tested files or original-archive source changed")
    num=sum(x["passed"] for x in receipts.values())
    skip=sum(x["skipped"] for x in receipts.values())
    if num!=194 or skip!=3:
        raise ValueError("full source-bound test gate must be exactly 194/3")
    receipt=dict(schema="ns_cowling.fortran_direct_eos_source_callback.v1",
        exact_scientific_parent="daf5de300b06009f8654f795aa99fe40d512666e",
        input_source_SHA256=before,
        tests=receipts,total_passed=num,total_skipped=skip,
        native_matcher_build_sha256=sha(exe),
        native_source_builds=builds,
        new_direct_source_fortran_callback_has_physical_restarts=True,
        original_GPPVA_Pais_sources_untouched=True,
        original_source_derived_EOS_unchanged=True,
        no_full96_stars_or_global_frequency_error_claim=True,
        wall_seconds=time.monotonic()-started,
        status="SCOPED_DIRECT_SOURCE_NUMERICAL_AND_IMPLEMENTATION_VERIFIED")
    (OUT/"VERIFICATION_RECEIPT.json").write_text(
        json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    print("DIRECT_SOURCE_FULL_SCOPED",num,"passed",skip,
          "skipped; source_hash_count",len(before),flush=True)
if __name__=="__main__":
    main()
