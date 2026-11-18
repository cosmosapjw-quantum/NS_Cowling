#!/usr/bin/env python3
"""Reproduce native TOV+Cowling f/p1 parity across original & extended EOS.
Reference source: original exact-contract Python EOS + DOP853 TOV/mode.
Native: exact EOS sampled once, Fortran DVODE TOV and Cowling shots.
"""
import hashlib
import json
import os
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
from compose_eos.stellar_cowling import CowlingSettings,solve_mode
from compose_eos.stellar_cowling_tov_native import solve_background_native
from compose_eos.stellar_cowling_native import NativeCowlingDVODE

OUT=ROOT/"data/cowling_native_tov_20261119"
CASES=[
    ("primary_n035","data/core_support_20261117/primary_v1/REPRESENTATION.json",0.35),
    ("primary_n1p4","data/core_support_20261117/primary_v1/REPRESENTATION.json",0.5377339776363226),
    ("primary_n070","data/core_support_20261117/primary_v1/REPRESENTATION.json",0.7),
    ("extension_medium_n115","data/stellar_extension_20261117/eos/medium/REPRESENTATION.json",1.15),
    ("extension_fine_n115","data/stellar_extension_20261117/eos/fine/REPRESENTATION.json",1.15),
]
BG=dict(surface_h=1e-4,rtol=1e-10,atol=1e-12,max_step_h=.002)
SETTING=CowlingSettings(surface_h=1e-4,rtol=1e-9,atol=1e-12)
NATIVE_SAMPLES=8193
MODE_ROWS=("frequency_hz","W_nodes","V_nodes","matching_residual",
           "relative_join_mismatch","surface_BC_relative","nfev")


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def result_slim(x):
    return {k:x[k] for k in MODE_ROWS}


def execute():
    OUT.mkdir(parents=True,exist_ok=True)
    record={
        "schema":"ns_cowling.native_tov_cowling_parity.v1",
        "parent_commit":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "source_sha256":{
            str(x.relative_to(ROOT)):digest(x) for x in (
                ROOT/"python_frontend/compose_eos/stellar_background.py",
                ROOT/"python_frontend/compose_eos/stellar_cowling.py",
                ROOT/"python_frontend/compose_eos/stellar_cowling_tov_native.py",
                ROOT/"python_frontend/compose_eos/stellar_cowling_native.py",
                ROOT/"fortran_eos/src/compose_eos_cowling_tov_dvode.f90",
                ROOT/"fortran_eos/src/compose_eos_cowling_dvode.f90")},
        "platform":platform.platform(),"compiler":subprocess.check_output(
            ["gfortran","--version"],text=True).splitlines()[0],
        "python":platform.python_version(),
        "threads":{k:os.getenv(k) for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS")},
        "definitions":{
            "metric":"(-,+,+,+) with r,m in geometric km",
            "background":"TOV h-sqrt(hstart-h) center, surface_h=1e-4",
            "mode":"equilibrium barotropic ell=2, matched determinant roots",
            "python_ref":"DOP853 with exact EOS at every RHS evaluation",
            "native":"PCHIP-sampled exact EOS; both ODEs with native Netlib DVODE Adams",
            "scope":"no Maxwell sharp interfaces, no full-native EOS generator",
            "TOV_controls":BG,
            "mode_controls":vars(SETTING),
            "reference_timing":"one Python invocation per mode",
            "native_timing":"seven warm invocations per mode after native TOV sampling",
        },
        "cases":[],
    }
    part=OUT/"BENCHMARK_PARTIAL.json"
    for tag,relative,nc in CASES:
        p=ROOT/relative;eos=ContractEnthalpyEOS.from_path(p)
        row=dict(case_id=tag,source_path=relative,source_sha256=digest(p),
                 central_density_fm3=nc,eos_max_density=eos.density_domain[1])
        before=time.perf_counter()
        bgp=solve_background(eos,nc,**BG)
        row["python_background_seconds"]=time.perf_counter()-before
        before=time.perf_counter()
        bgn=solve_background_native(eos,nc,native_samples=NATIVE_SAMPLES,**BG)
        row["native_background_seconds_including_python_eos_sampling"]=time.perf_counter()-before
        row["background"]={
            "python":dict(mass_msun=bgp.mass_msun,radius_km=bgp.radius_km),
            "native":dict(mass_msun=bgn.mass_msun,radius_km=bgn.radius_km,
                          ode_rhs_calls=bgn.nfev,profile_points=len(bgn.s)),
            "delta_mass_msun":bgn.mass_msun-bgp.mass_msun,
            "delta_radius_km":bgn.radius_km-bgp.radius_km
        }
        native=NativeCowlingDVODE(bgn,SETTING,samples=NATIVE_SAMPLES)
        brackets,_=native.scan_brackets(np.linspace(1000.,13000.,85))
        matches=[]
        for br in brackets:
            out=native.solve_mode(br,profile_points=240)
            if out["W_nodes"]==len(matches) and out["V_nodes"]==len(matches):
                matches.append((br,out))
            if len(matches)>=2:
                break
        if len(matches)!=2:
            raise RuntimeError(f"{tag}: insufficient correctly labelled f,p1 roots: {len(matches)}")
        row["modes"]={}
        for label,(br,initial) in zip(("f","p1"),matches):
            before=time.perf_counter()
            python_result=solve_mode(bgp,br,SETTING,profile_points=240)
            python_seconds=time.perf_counter()-before
            warm=[]
            for _ in range(7):
                before=time.perf_counter()
                native_result=native.solve_mode(br,profile_points=240)
                warm.append(time.perf_counter()-before)
            row["modes"][label]={
                "root_bracket_hz":[float(z) for z in br],
                "python_exact_dop853":result_slim(python_result),
                "native_tov_native_cowling":result_slim(native_result),
                "frequency_difference_hz":native_result["frequency_hz"]-python_result["frequency_hz"],
                "python_solve_seconds":python_seconds,
                "native_warm_solve_median_seconds":statistics.median(warm),
                "native_warm_seconds":warm,
                "speedup_warm_shooting":python_seconds/statistics.median(warm)
            }
            if python_result["W_nodes"]!=native_result["W_nodes"]:
                raise RuntimeError(f"{tag}/{label}: mismatched W nodes")
            print("CASE_MODE",tag,label,
                  "native_frequency",native_result["frequency_hz"],
                  "python_frequency",python_result["frequency_hz"],
                  "delta",row["modes"][label]["frequency_difference_hz"],
                  "warm_speedup",row["modes"][label]["speedup_warm_shooting"],flush=True)
            part.write_text(json.dumps(record,indent=2)+"\n")
        # resolution study: same EOS/center inputs, repeated native TOV + mode shots
        row["resolution_convergence"]={}
        for ns in (4097,16385):
            bg=solve_background_native(eos,nc,native_samples=ns,**BG)
            solver=NativeCowlingDVODE(bg,SETTING,samples=ns)
            found={}
            for mode_name,(br,_) in zip(("f","p1"),matches):
                result=solver.solve_mode(br,profile_points=240)
                found[mode_name]=dict(frequency_hz=result["frequency_hz"],
                                      W_nodes=result["W_nodes"],
                                      V_nodes=result["V_nodes"])
            row["resolution_convergence"][str(ns)]={"M_msun":bg.mass_msun,
                  "R_km":bg.radius_km,"modes":found}
        record["cases"].append(row)
        part.write_text(json.dumps(record,indent=2)+"\n")
    (OUT/"BENCHMARK.json").write_text(json.dumps(record,indent=2)+"\n")
    print("FINAL_NATIVE_TOV_COWLING_BENCHMARK",len(record["cases"]),
          str(OUT/"BENCHMARK.json"),flush=True)


if __name__=="__main__":
    execute()
