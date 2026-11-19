#!/usr/bin/env python3
"""Fresh scoped end-to-end scientific verification for RMF onepoint/figures.

Tests all inherited stage05 + EOS bundle regressions and 3 new onepoint test
modules, bounds-checked Fortran and original source bytes. Never call a pass
unless a fresh test exit=0 with explicit source SHAs before and after.
"""
from __future__ import annotations
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from verify_missing_eos_nested_20261119 import TESTS as ORIGINAL_TESTS

DATA=ROOT/"data/cowling_original_rmf_onepoint_20261119"
OUT=DATA/"verification"
TESTS=ORIGINAL_TESTS+[
    "test_cowling_rmf_onepoint_20261119.py",
    "test_cowling_rmf_onepoint_evidence_20261119.py",
    "test_cowling_rmf_onepoint_figure_20261119.py",
]
PYTHON_FILES=[
"python_frontend/compose_eos/rmf_onepoint_components.py",
"python_frontend/tests/test_cowling_rmf_onepoint_20261119.py",
"python_frontend/tests/test_cowling_rmf_onepoint_evidence_20261119.py",
"python_frontend/tests/test_cowling_rmf_onepoint_figure_20261119.py",
"scripts/fetch_shen_eos4_t00_20261119.py",
"scripts/run_cowling_original_rmf_onepoint_20261119.py",
"scripts/check_rmf_onepoint_fortran_parity_20261119.py",
"scripts/plot_cowling_rmf_compose_onepoint_20261119.py",
"scripts/plot_cowling_compose_density_trend_20261119.py",
"scripts/package_cowling_rmf_onepoint_20261119.py",
"scripts/verify_cowling_rmf_onepoint_20261119.py",
]
FORTRAN_FILES=[
"fortran_eos/src/compose_eos_types.f90",
"fortran_eos/src/compose_eos_tf_models.f90",
"fortran_eos/app/rmf_original_onepoint_components.f90",
"fortran_eos/src/compose_eos_cowling_dvode.f90",
"fortran_eos/src/compose_eos_cowling_tov_dvode.f90",
]
EVIDENCE_FILES=[
"data/cowling_original_rmf_onepoint_20261119/RESULT.json",
"data/cowling_original_rmf_onepoint_20261119/FORTRAN_COMPONENT_PARITY.json",
"data/cowling_original_rmf_onepoint_20261119/figures/COMPOSE_vs_LOCAL_RMF_exact_onepoint.png",
"data/cowling_original_rmf_onepoint_20261119/figures/COMPOSE_vs_LOCAL_RMF_exact_onepoint.svg",
"data/cowling_original_rmf_onepoint_20261119/figures/COMPOSE_vs_LOCAL_RMF_exact_onepoint.csv",
"data/cowling_original_rmf_onepoint_20261119/figures/COMPOSE_vs_LOCAL_RMF_original_density_trend.png",
"data/cowling_original_rmf_onepoint_20261119/figures/COMPOSE_vs_LOCAL_RMF_original_density_trend.svg",
"data/cowling_original_rmf_onepoint_20261119/figures/COMPOSE_vs_LOCAL_RMF_original_density_trend.csv",
"data/cowling_original_rmf_onepoint_20261119/artifacts/COWLING_RMF_ONEPOINT_COMPOSE_AUTHOR_COMPARE_20261119.zip",
"data/cowling_original_rmf_onepoint_20261119/external_reference/shen_EOS4_original_eos4.t00.gz",
"data/cowling_rmf_reconciliation_20261119/EXACT_DENSITY_AUDIT.json",
"docs/COWLING_RMF_ONEPOINT_RESULT_20261119.md",
"docs/COWLING_RMF_ONEPOINT_PROTOCOL_20261119.md",
]
ROOT_SOURCES=PYTHON_FILES+FORTRAN_FILES+EVIDENCE_FILES
ENV=dict(os.environ)
ENV["OPENBLAS_NUM_THREADS"]="1"
ENV["OMP_NUM_THREADS"]="1"
ENV["PYTHONPATH"]="python_frontend"

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run(command,log):
    p=subprocess.run(command,cwd=ROOT,env=ENV,stdout=subprocess.PIPE,
                     stderr=subprocess.STDOUT,text=True,check=False)
    log.write_text(p.stdout)
    print("RUN",command[0],"EXIT",p.returncode,"LAST",
          p.stdout.splitlines()[-1:] if p.stdout else "",flush=True)
    if p.returncode:
        raise RuntimeError(command[0]+" failed "+p.stdout[-1600:])
    return p.stdout

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    before={name:sha(ROOT/name) for name in ROOT_SOURCES}
    source=DATA/"external_reference/shen_EOS4_original_eos4.t00.gz"
    if sha(source)!="346b2701f7372817c516c3c945fdee0b3da5af9f08e2fe83f743951d1b5ebaef":
        raise ValueError("original author source SHA mismatch")
    pkg=DATA/"artifacts/COWLING_RMF_ONEPOINT_COMPOSE_AUTHOR_COMPARE_20261119.zip"
    with zipfile.ZipFile(pkg) as z:
        if z.testzip() is not None:
            raise ValueError("report package CRC invalid")
        manifest=json.loads(z.read("INTERNAL_SOURCE_IDENTITY_SHA256.json"))
        for name,spec in manifest["members"].items():
            if hashlib.sha256(z.read(name)).hexdigest()!=spec["sha256"]:
                raise ValueError("report member SHA mismatch")
    if subprocess.run(["git","diff","--check"],cwd=ROOT).returncode:
        raise RuntimeError("git diff --check rejected working source")
    run([sys.executable,"-m","compileall","-q",*PYTHON_FILES],
        OUT/"PYTHON_COMPILEALL.log")
    build_dir=OUT/"f90_compile"
    build_dir.mkdir(exist_ok=True)
    native=run(["gfortran","-std=f2018","-O0","-fcheck=all","-Wall","-Wextra",
                "-ffree-line-length-none","-J",str(build_dir),"-I",str(build_dir),
                *FORTRAN_FILES[:3],"-o",str(build_dir/"rmf_onepoint")],
                OUT/"COMPONENT_FORTRAN_BUILD.log")
    for name in FORTRAN_FILES[3:]:
        out=build_dir/(Path(name).stem+".o")
        try:
            run(["gfortran","-c","-O0","-fcheck=all","-Wall","-Wextra",
                 "-ffree-line-length-none",name,"-o",str(out)],
                 OUT/(Path(name).stem+".BUILD.log"))
        finally:
            out.unlink(missing_ok=True)
    paths=[f"python_frontend/tests/{n}" for n in TESTS]
    stdout=run([sys.executable,"-m","pytest","-q",*paths,"--disable-warnings"],
               OUT/"TEST_RUN.log")
    m=re.search(r"(\d+) passed, (\d+) skipped",stdout)
    if not m:
        raise RuntimeError("no successful pytest pass/skip summary")
    passes,skips=map(int,m.groups())
    if passes<171 or skips!=3:
        raise RuntimeError("unexpected regression scope counts")
    after={name:sha(ROOT/name) for name in ROOT_SOURCES}
    if before!=after:
        raise RuntimeError("test source or original input changed during verification")
    report=dict(schema="ns_cowling.rmf_original_author_onepoint_with_compose_figures_verified.v1",
        frozen_git_parent=subprocess.check_output(
            ["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        exact_tested_files_sha256=before,
        original_author_shen_compressed_source_sha256=sha(source),
        original_author_shen_source_binary_redistributed=False,
        GPPVA_provider_micro_fields_or_code_obtained=False,
        original_GPPVA_source_points_compared=213,
        exact_onepoint_source_models=["FSU2","TM1e","DDME2"],
        report_package=dict(path=str(pkg.relative_to(ROOT)),
            sha256=sha(pkg),bytes=pkg.stat().st_size,n_members=len(manifest["members"])),
        pytest=dict(passed=passes,skipped=skips,
           log_sha256=sha(OUT/"TEST_RUN.log"),exit_code=0),
        compiler=dict(GNU_Fortran_version=subprocess.check_output(
            ["gfortran","--version"],text=True).splitlines()[0],
            bounds_checked_native_components=True,
            native_source_SHA256={k:before[k] for k in FORTRAN_FILES},
            build_log_sha256=sha(OUT/"COMPONENT_FORTRAN_BUILD.log")),
        python_compileall_passed=True,
        git_diff_check_passed=True,
        original_data_and_plot_SHA_unchanged_during_test=True,
        conclusion="Author TM1e bulk agrees at one exact zero-T fixed-Yp point; Pais GPPVA original field data still missing",
        status="SCOPED_SOURCE_AUTHOR_RMF_AND_COMPOSE_FIGURES_IMPLEMENTATION_VERIFIED",
        wall_seconds=time.monotonic()-started)
    (OUT/"VERIFICATION_RECEIPT.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print("RMF_ONEPOINT_FULL_SCOPED_VERIFIED",passes,skips,
          "input_sha_count",len(before),"ZIP",report["report_package"]["sha256"],flush=True)

if __name__=="__main__":main()
