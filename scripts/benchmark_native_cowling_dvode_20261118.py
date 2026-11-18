#!/usr/bin/env python3
"""Reproduce smooth 1.4 Msun exact-contract Cowling f/p1 accuracy and timing.

Python TOV background is common to both paths. The baseline calls the exact
source EOS at every DOP853 RHS evaluation; native shots use a one-time PCHIP
coefficient table with explicit one-sided sound-speed breakpoints.
"""
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.stellar_background import ContractEnthalpyEOS,solve_background
from compose_eos.stellar_cowling import CowlingSettings,solve_mode
from compose_eos.stellar_cowling_native import NativeCowlingDVODE

EOS=ROOT/"data/core_support_20261117/primary_v1/REPRESENTATION.json"
CENTRAL=ROOT/"data/stellar_consumers_20261117/tov/CANONICAL_1P4.json"
OUT=ROOT/"data/cowling_native_dvode_20261118"
CASES=[("f",[2250.,2500.]),("p1",[6500.,6750.])]
FIELDS=["frequency_hz","W_nodes","V_nodes","nfev","matching_residual",
        "relative_join_mismatch","surface_BC_relative"]
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run():
    OUT.mkdir(exist_ok=True,parents=True)
    central=json.loads(CENTRAL.read_text())
    es=ContractEnthalpyEOS.from_path(EOS)
    start=time.perf_counter()
    bg=solve_background(es,central["central_density_fm3"],
                        surface_h=1e-4,rtol=1e-10,atol=1e-12,max_step_h=.002)
    setup_background=time.perf_counter()-start
    setting=CowlingSettings(surface_h=1e-4,rtol=1e-9,atol=1e-12)
    record=dict(schema="ns_cowling.cowling_dvode_smooth_parity.v1",
        git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        dirty_status=subprocess.check_output(["git","status","--short"],cwd=ROOT,text=True).splitlines(),
        machine=platform.platform(),processor=platform.processor(),
        compiler=subprocess.check_output(["gfortran","--version"],text=True).splitlines()[0],
        python=platform.python_version(),numpy_thread_env={
            k:os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS")},
        eos_contract_sha256=sha(EOS),central_configuration_sha256=sha(CENTRAL),
        fortran_source_sha256=sha(ROOT/"fortran_eos/src/compose_eos_cowling_dvode.f90"),
        background=dict(mass_msun=bg.mass_msun,radius_km=bg.radius_km,
                        h_start=bg.h_start,h_c=bg.h_c,
                        seconds=setup_background),
        settings=vars(setting),run_condition="single CPU process; identical physical target and mode brackets",
        scope="smooth equilibrium barotrope; Maxwell interfaces unsupported",
        cases={})
    for label,br in CASES:
        start=time.perf_counter()
        baseline=solve_mode(bg,br,setting,profile_points=400)
        seconds=time.perf_counter()-start
        record["cases"][label]={"bracket_hz":br,
             "python_exact_dop853":{"seconds_single_run":seconds,
                     "output":{k:baseline[k] for k in FIELDS}}}
        print("REFERENCE",label,"seconds",round(seconds,3),
              "freq",baseline["frequency_hz"],flush=True)
        (OUT/"BENCHMARK_PARTIAL.json").write_text(json.dumps(record,indent=2)+"\n")
    start=time.perf_counter()
    native=NativeCowlingDVODE(bg,setting,samples=8193)
    sampled_seconds=time.perf_counter()-start
    first=native.solve_mode(CASES[0][1],profile_points=400)
    library=Path(native._lib._name)
    record["native_setup"]=dict(sample_construction_seconds=sampled_seconds,
                                compile_and_first_solve_included_separately=True,
                                compiled_lib_sha256=sha(library),
                                requested_samples=8193,actual_samples=len(native.s),
                                source_breakpoint_count=len(native.break_s),
                                library_path=str(library))
    for label,br in CASES:
        samples=[]; output=None
        for i in range(7):
            start=time.perf_counter()
            output=native.solve_mode(br,profile_points=400)
            samples.append(time.perf_counter()-start)
        comparison=record["cases"][label]
        reference=comparison["python_exact_dop853"]
        record["cases"][label]["fortran_dvode_adams"]={
            "seconds_7_runs":samples,
            "median_seconds":statistics.median(samples),
            "output":{k:output[k] for k in FIELDS},
            "frequency_minus_reference_hz":output["frequency_hz"]-
                reference["output"]["frequency_hz"],
            "speedup_warm_shooting_vs_python_single_run":
                reference["seconds_single_run"]/statistics.median(samples)}
        print("NATIVE",label,record["cases"][label]["fortran_dvode_adams"],
              flush=True)
        (OUT/"BENCHMARK_PARTIAL.json").write_text(json.dumps(record,indent=2)+"\n")
    (OUT/"BENCHMARK.json").write_text(json.dumps(record,indent=2)+"\n")
    print("DELIVERED",OUT/"BENCHMARK.json",flush=True)

if __name__=="__main__":
    run()
