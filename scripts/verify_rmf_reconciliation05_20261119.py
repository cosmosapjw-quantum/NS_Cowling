#!/usr/bin/env python3
"""Fresh source-locked verification for COWLING-EOS-SOURCE-CALIBRATION-05.

Build native modules with Fortran array bounds checking, run inherited
03/04 tests plus new source, charge and same-crust stellar checks,
and verify SHA inputs from the actual authenticated provider files.
"""
from __future__ import annotations
import hashlib,json,os,subprocess,sys,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data/cowling_rmf_reconciliation_20261119"
OUT=DATA/"verification"
sys.path.insert(0,str(ROOT/"scripts"))
from verify_calibration_acceptance_20261119 import TESTS as STAGE_04_TESTS

TESTS=["test_cowling_rmf_reconciliation_20261119.py",
       "test_cowling_rmf_reconciliation_evidence_20261119.py",*STAGE_04_TESTS]
SOURCE=[
"fortran_eos/src/compose_eos_cowling_dvode.f90",
"fortran_eos/src/compose_eos_cowling_tov_dvode.f90",
"python_frontend/compose_eos/gppva_rmf.py",
"python_frontend/compose_eos/rmf_source_reconciliation.py",
"python_frontend/compose_eos/source_calibration.py",
"python_frontend/compose_eos/stellar_mode_continuation.py",
"python_frontend/compose_eos/stellar_sharp_native.py",
"python_frontend/compose_eos/stellar_cowling_native.py",
]+[f"python_frontend/tests/{v}" for v in TESTS]+[
"scripts/audit_rmf_core_provenance_20261119.py",
"scripts/audit_compose_charge_thermo_20261119.py",
"scripts/audit_fixed_fraction_rmf_20261119.py",
"scripts/run_paired_core_source_stars_20261119.py",
"scripts/audit_compose_mr_provenance_20261119.py",
"scripts/build_rmf_source_reconciliation_outcome_20261119.py"]

SOURCE_EOS=ROOT/"data/cowling_maxwell_native_20261119/provider_original_unified"
S04=ROOT/"data/cowling_eos_source_calibration_20261119/source_audit/RESULT.json"
env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",
         PYTHONPATH="python_frontend")

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def cmd(args,log):
    out=subprocess.run(args,cwd=ROOT,env=env,capture_output=True,
                       text=True)
    log.write_text(out.stdout+out.stderr)
    print("STAGE05_VERIFY",args[0],"exit",out.returncode,
          "last",((out.stdout+out.stderr).splitlines()[-1:] or [""])[0],flush=True)
    if out.returncode:raise RuntimeError(str(args[:4])+" failed: "+(out.stdout+out.stderr)[-1500:])
    return out

def main():
    OUT.mkdir(exist_ok=True,parents=True)
    start=time.monotonic()
    input_hashes={p:digest(ROOT/p) for p in SOURCE}
    source_prev=json.loads(S04.read_text())["models"]
    source_hashes={}
    for name in ("FSU2","TM1e","DDME2"):
        source=SOURCE_EOS/name/"raw/eos.thermo"
        data_sha=digest(source)
        original_expected=source_prev[name]["raw_eos_source_sha256"]["eos.thermo"]
        if data_sha!=original_expected:
            raise ValueError("original GPPVA provider EOS changed from stage04 "+name)
        source_hashes[name]=data_sha
    code=[s for s in SOURCE if s.endswith(".py")]
    cmd([sys.executable,"-m","compileall","-q",*code],OUT/"PYTHON_COMPILEALL.log")
    builds=[]
    for rel in SOURCE[:2]:
        o=OUT/(Path(rel).stem+".o")
        try:
            args=["gfortran","-c","-O0","-fcheck=all","-Wall","-Wextra",
                  "-ffree-line-length-none",rel,"-o",str(o)]
            out=cmd(args,OUT/(Path(rel).stem+".build.log"))
            builds.append(dict(source=rel,status="PASS",
                     compiler_log_sha256=digest(OUT/(Path(rel).stem+".build.log"))))
        finally:o.unlink(missing_ok=True)
    tests=[f"python_frontend/tests/{t}" for t in TESTS]
    result=cmd([sys.executable,"-m","pytest","-q",*tests,
                 "--disable-warnings"],OUT/"TEST_RUN.log")
    import re
    match=re.search(r"(\d+) passed, (\d+) skipped",result.stdout)
    if match is None:raise RuntimeError("pytest result summary absent")
    npassed,nskipped=map(int,match.groups())
    if npassed<135 or nskipped!=3:
        raise RuntimeError("unexpected regression scope or failure count")
    if input_hashes!={p:digest(ROOT/p) for p in SOURCE}:
        raise RuntimeError("source changes during regression invalidated test")
    staged=subprocess.run(["git","diff","--cached","--check"],cwd=ROOT,
                          capture_output=True,text=True)
    if staged.returncode:raise RuntimeError("staged whitespace check failed: "+staged.stderr)
    receipt=dict(schema="ns_cowling.cowling_source_calibration05_verified.v1",
       parent_04_commit="f191abe843d81c6a1e5797dd9edc9a1558f095fa",
       working_head_at_execution=subprocess.check_output(
          ["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
       python=sys.version,
       gfortran=subprocess.check_output(["gfortran","--version"],
                        text=True).splitlines()[0],
       tested_files_sha256=input_hashes,
       original_provider_raw_eos_thermo_sha256=source_hashes,
       original_provider_compositions_unmodified=True,
       source_audit_root_receipt_sha256=digest(DATA/"SOURCE_RECONCILIATION_OUTCOME.json"),
       checkout_staged_diff_check_exit=0,
       native_fortran_bounds_checked_compile=builds,
       pytest=dict(passed=npassed,skipped=nskipped,exit_code=0,
          log_sha256=digest(OUT/"TEST_RUN.log"),
          command=[sys.executable,"-m","pytest","-q",*tests,"--disable-warnings"]),
       passed_syntax=True,
       elapsed_seconds=time.monotonic()-start,
       native_subagent="NOT_RUN",
       certification_scope="Cold locally neutral 3-layer TOV and slow/rapid Cowling under original source identity and nonfitted RMF calibration only",
       status="PASS_SCOPED_IMPLEMENTATION_SOURCE_LOCKED")
    (OUT/"VERIFICATION_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print("SOURCE_CALIBRATION05_SCOPED_VERIFIED",npassed,nskipped,
          "elapsed",round(receipt["elapsed_seconds"],2),flush=True)

if __name__=="__main__":main()
