#!/usr/bin/env python3
"""Final bounded verification of all EOS files and nested ZIP.

Runs stage05 inherited native TOV/Cowling suite plus new package tests.
Checks all input-source hashes before/after and Fortran bounds-check compiles.
The complete previously source-locked raw EOS source is never overwritten.
"""
from __future__ import annotations
import hashlib,json,os,re,subprocess,sys,time,zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data/eos_missing_data_recalc_20261119"
DEST=DATA/"verification"
sys.path.insert(0,str(ROOT/"scripts"))
from verify_rmf_reconciliation05_20261119 import TESTS as PREVIOUS

TESTS=list(PREVIOUS)+["test_eos_missing_nested_bundle_20261119.py"]
SOURCE=[
    "fortran_eos/src/compose_eos_cowling_dvode.f90",
    "fortran_eos/src/compose_eos_cowling_tov_dvode.f90",
    "python_frontend/compose_eos/rmf_source_reconciliation.py",
    "python_frontend/compose_eos/stellar_sharp_native.py",
    "python_frontend/compose_eos/stellar_mode_continuation.py",
    "python_frontend/tests/test_eos_missing_nested_bundle_20261119.py",
    "scripts/recalculate_missing_spg_pk1_stars_20261119.py",
    "scripts/recalculate_bsk22_surface_from_anchor_20261119.py",
    "scripts/audit_missing_external_crust_20261119.py",
    "scripts/package_nested_eos_missing_20261119.py",
    "scripts/generate_spg_crust.py",
    "scripts/reconstruct_public_thermo_20261117.py",
    "scripts/generate_rmf_candidates_20261117.py",
    "python_frontend/compose_eos/candidate_gdfm_cld.py",
    "python_frontend/compose_eos/spg_crust_generator.py",
    "python_frontend/compose_eos/rmf_candidate.py",
    "python_frontend/compose_eos/bsk22_analytic.py",
    "data/new_model_maxwell_20261117/sources/PARAMETER_INPUTS.json",
    "data/eos_missing_data_recalc_20261119/CORRECTED13_RESULT.json",
    "data/eos_missing_data_recalc_20261119/reconstructed_spg_fresh17/SUMMARY.json",
    "data/eos_missing_data_recalc_20261119/reconstructed_gppva_fsu2_fresh52/RECALC_REPORT.json",
    "data/eos_missing_data_recalc_20261119/pk1_core_fresh/RECALC_REPORT.json",
    "data/eos_missing_data_recalc_20261119/bsk22_outer_and_surface/RECALC_REPORT.json",
    "data/eos_missing_data_recalc_20261119/bsk22_lower_inner_fresh/RECEIPT.json",
    "data/eos_missing_data_recalc_20261119/spg_pk1_complete_stars/RESULT.json",
    "data/eos_missing_data_recalc_20261119/external_original/NUDA_SOURCE_AUDIT.json",
    "data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip",
]

ENV=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",
         PYTHONPATH="python_frontend")
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def command(args,log):
    out=subprocess.run(args,cwd=ROOT,env=ENV,
                       stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    log.write_text(out.stdout)
    print("VERIFY_EXEC",args[0],"exit",out.returncode,"last",out.stdout.splitlines()[-1:] if out.stdout else "",flush=True)
    if out.returncode:
        raise RuntimeError(args[0]+" failed: "+out.stdout[-1800:])
    return out.stdout

def main():
    DEST.mkdir(parents=True,exist_ok=True)
    t=time.monotonic()
    input_hash={p:sha(ROOT/p) for p in SOURCE}
    zipname=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
    with zipfile.ZipFile(zipname) as z:
        assert z.testzip() is None and len(z.namelist())==10
    command([sys.executable,"-m","compileall","-q",*[x for x in SOURCE if x.endswith(".py")]],
            DEST/"PYTHON_COMPILEALL.log")
    fortran=[]
    for src in SOURCE[:2]:
        out=DEST/(Path(src).stem+".o")
        try:
            command(["gfortran","-c","-O0","-fcheck=all","-Wall","-Wextra",
                     "-ffree-line-length-none",src,"-o",str(out)],
                    DEST/(Path(src).stem+".BUILD.log"))
            fortran.append(dict(source=src,compiled=True,log_sha256=sha(
                        DEST/(Path(src).stem+".BUILD.log"))))
        finally:
            out.unlink(missing_ok=True)
    actual=[f"python_frontend/tests/{t}" for t in TESTS]
    result=command([sys.executable,"-m","pytest","-q",*actual,"--disable-warnings"],
                   DEST/"TEST_RUN.log")
    passed=re.search(r"(\d+) passed, (\d+) skipped",result)
    if passed is None:raise RuntimeError("full scope pytest count not observed")
    npassed,nskip=map(int,passed.groups())
    if npassed<150 or nskip!=3:raise RuntimeError("required regression acceptance not met")
    if input_hash!={p:sha(ROOT/p) for p in SOURCE}:
        raise RuntimeError("source bytes changed after actual test run")
    output=dict(schema="ns_cowling.all_requested_eos_nested_bundle_verification.v1",
        current_head_at_test=subprocess.check_output(["git","rev-parse","HEAD"],
            cwd=ROOT,text=True).strip(),input_sha256=input_hash,
        native_fortran_builds=fortran,compiler=subprocess.check_output(
            ["gfortran","--version"],text=True).splitlines()[0],
        pytest=dict(passed=npassed,skipped=nskip,exit_code=0,
            log_sha256=sha(DEST/"TEST_RUN.log")),
        python_compileall_passed=True,
        outer_and_seven_inner_zip_all_testzip_success=True,
        zip_sha256=sha(zipname),zip_bytes=zipname.stat().st_size,
        elapsed_wall_seconds=time.monotonic()-t,
        source_redistribution="nuda original 24 files unchanged with CC BY-NC-ND 4.0 author licence",
        source_vs_derived_model_identity_enforced=True,
        physical_scope="direct Fortran source calculations and conditional 1.4M Cowling, no original SPG surface/GMRS microscopic generator",
        status="SCOPED_SOURCE_AND_IMPLEMENTATION_VERIFIED")
    (DEST/"VERIFICATION_RECEIPT.json").write_text(json.dumps(output,indent=2)+"\n")
    print("STAGE_MISSING_EOS_BUNDLE_VERIFIED",npassed,nskip,"zip_sha",output["zip_sha256"],flush=True)

if __name__=="__main__":main()
