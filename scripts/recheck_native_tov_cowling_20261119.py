#!/usr/bin/env python3
"""Final current-code provenance check; reuse unchanged exact Python reference.

The expensive exact DOP853 data have original source SHA binding. Current
Fortran native computations are rerun at current source SHA for every case.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.stellar_background import ContractEnthalpyEOS,solve_background
from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_cowling_tov_native import solve_background_native
from compose_eos.stellar_cowling_native import NativeCowlingDVODE

OUT=ROOT/"data/cowling_native_tov_20261119"
FIELDS=("frequency_hz","W_nodes","V_nodes","nfev",
        "matching_residual","relative_join_mismatch","surface_BC_relative")
OPTIONS=dict(surface_h=1e-4,rtol=1e-10,atol=1e-12,max_step_h=.002)
MODE=CowlingSettings(surface_h=1e-4,rtol=1e-9,atol=1e-12)


def hashfile(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    reference_path=OUT/"BENCHMARK.json"
    old=json.loads(reference_path.read_text())
    assert len(old["cases"])==5
    # Only Python exact-contract references may be reused; if either physical
    # reference implementation changed, fail rather than silently reattribute.
    for source in ("python_frontend/compose_eos/stellar_background.py",
                   "python_frontend/compose_eos/stellar_cowling.py"):
        assert old["source_sha256"][source]==hashfile(ROOT/source),source
    meta=dict(schema="ns_cowling.native_tov_cowling_final_recheck.v1",
       original_reference_sha256=hashfile(reference_path),
       reference_type="Python DOP853 exact-contract baseline copied from immutable original run by matching source+inputs",
       current_commit_at_run=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
       compiler=subprocess.check_output(["gfortran","--version"],text=True).splitlines()[0],
       hardware=platform.platform(),
       current_source_sha256={
          p:hashfile(ROOT/p) for p in (
             "fortran_eos/src/compose_eos_cowling_tov_dvode.f90",
             "fortran_eos/src/compose_eos_cowling_dvode.f90",
             "python_frontend/compose_eos/stellar_cowling_tov_native.py",
             "python_frontend/compose_eos/stellar_cowling_native.py"
          )},
       controls=dict(background=OPTIONS,mode=vars(MODE),
                     native_samples=8193,
                     native_method="DVODE Adams",
                     source_contract="Original exact piece laws evaluated once in Python; no sharp interfaces",
                     timing="No JIT/recompile cost; background EOS sampling and native ODE included in measured background time"),
       cases=[])
    for row in old["cases"]:
        eosfile=ROOT/row["source_path"]
        assert hashfile(eosfile)==row["source_sha256"]
        eos=ContractEnthalpyEOS.from_path(eosfile)
        nc=row["central_density_fm3"]
        start=time.perf_counter()
        bg=solve_background_native(eos,nc,native_samples=8193,**OPTIONS)
        setup_seconds=time.perf_counter()-start
        start=time.perf_counter()
        solver=NativeCowlingDVODE(bg,MODE,samples=8193)
        mode_init_seconds=time.perf_counter()-start
        rec=dict(case_id=row["case_id"],
                 eos_sha256=row["source_sha256"],
                 central_density_fm3=nc,
                 native_background=dict(
                     M_msun=bg.mass_msun,R_km=bg.radius_km,
                     source_knots=len(bg.s),tov_nfev=bg.nfev),
                 exact_python_background=row["background"]["python"],
                 native_background_setup_seconds=setup_seconds,
                 native_mode_setup_seconds=mode_init_seconds,
                 modes={})
        medians=[]
        for label in ("f","p1"):
            src=row["modes"][label]
            bracket=src["root_bracket_hz"]
            out=[]
            samples=[]
            for _ in range(7):
                start=time.perf_counter()
                result=solver.solve_mode(bracket,profile_points=240)
                samples.append(time.perf_counter()-start)
                out.append(result)
            x=out[-1]
            py=src["python_exact_dop853"]
            delta=x["frequency_hz"]-py["frequency_hz"]
            rec["modes"][label]={
               "bracket_hz":bracket,
               "python_original":py,
               "native_current":{k:x[k] for k in FIELDS},
               "delta_hz":delta,
               "relative_frequency_error":abs(delta)/py["frequency_hz"],
               "native_seven_wall_seconds":samples,
               "native_warm_median_seconds":statistics.median(samples),
               "python_original_wall_seconds":src["python_solve_seconds"],
               "speedup_warm_vs_python":src["python_solve_seconds"]/statistics.median(samples)
            }
            assert x["W_nodes"]==py["W_nodes"] and x["V_nodes"]==py["V_nodes"],(row["case_id"],label)
            medians.append(statistics.median(samples))
            print("FINAL_CASE",row["case_id"],label,
                  "native",x["frequency_hz"],"py",py["frequency_hz"],
                  "delta_hz",delta,"seconds",statistics.median(samples),flush=True)
        native_total=setup_seconds+mode_init_seconds+sum(medians)
        python_total=row["python_background_seconds"]+sum(
            row["modes"][label]["python_solve_seconds"] for label in ("f","p1"))
        rec["full_pipeline"]=dict(
            python_original_seconds=python_total,
            native_current_seconds= native_total,
            speedup_python_over_native=python_total/native_total,
            note="Two modes per source with preidentified brackets; native setup includes Python EOS sampling, not compilation")
        rec["sample_convergence"]={}
        for resolution in (4097,16385):
            nbg=solve_background_native(eos,nc,native_samples=resolution,**OPTIONS)
            ns=NativeCowlingDVODE(nbg,MODE,samples=resolution)
            nrow=dict(M_msun=nbg.mass_msun,R_km=nbg.radius_km,modes={})
            for label in ("f","p1"):
                br=rec["modes"][label]["bracket_hz"]
                result=ns.solve_mode(br,profile_points=240)
                nrow["modes"][label]=result["frequency_hz"]
            rec["sample_convergence"][str(resolution)]=nrow
        meta["cases"].append(rec)
        (OUT/"FINAL_RECHECK_PARTIAL.json").write_text(json.dumps(meta,indent=2)+"\n")
    final=OUT/"FINAL_RECHECK.json"
    final.write_text(json.dumps(meta,indent=2)+"\n")
    print("FINAL_RECHECK_COMPLETE",str(final),len(meta["cases"]),flush=True)

if __name__=="__main__":
    main()
