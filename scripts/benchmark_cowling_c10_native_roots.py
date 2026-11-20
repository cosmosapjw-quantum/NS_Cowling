#!/usr/bin/env python3
"""Timing frontend only: Fortran root polishing vs legacy SciPy dispatch.

Existing EOS/PCHIP construction and native Fortran ODE physics unchanged.
"""
from __future__ import annotations

import hashlib
import json
import os
import statistics
import sys
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))

from compose_eos.stellar_background import ContractEnthalpyEOS
from compose_eos.stellar_cowling_tov_native import solve_background_native
from compose_eos.stellar_cowling_native import NativeCowlingDVODE
from compose_eos.stellar_cowling import CowlingSettings

SOURCE=ROOT/"data/core_support_20261117/primary_v1/REPRESENTATION.json"
CENTRAL=ROOT/"data/stellar_consumers_20261117/tov/CANONICAL_1P4.json"
OUT=ROOT/"data/cowling_fortran_c10_native_root_20261120/TIMING.json"
SOURCE_SHA="591ac27aa99448512e291684c98f8273efc729d17c63be88e331c20e83bd2653"
BRACKETS=(("f",(2250.,2500.),0),("p1",(6500.,6750.),1))


def run():
    if hashlib.sha256(SOURCE.read_bytes()).hexdigest()!=SOURCE_SHA:
        raise RuntimeError("source EOS drifted")
    eos=ContractEnthalpyEOS.from_path(SOURCE)
    nc=json.loads(CENTRAL.read_text())["central_density_fm3"]
    t0=time.perf_counter()
    star=solve_background_native(eos,nc,surface_h=1e-4,rtol=1e-10,
               atol=1e-12,max_step_h=.002,native_samples=4097)
    bg_s=time.perf_counter()-t0
    settings=CowlingSettings(surface_h=1e-4,rtol=1e-9,atol=1e-12,
                            max_step_h=.002,root_xtol_hz=1e-6)
    record=dict(schema="ns_cowling.c10.fortran_brent_real_source_benchmark.v1",
        scientific_parent="0952a5c8ca34c5d94136a2c89c705e000c21ef3f",
        source_sha256=SOURCE_SHA,
        compiler="GNU Fortran native DVODE from checked-out source",
        constraints="smooth equilibrium-barotropic, no frozen-composition/g-mode claims",
        threading={k:os.environ.get(k) for k in
                   ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS")},
        background=dict(central_density_fm3=nc,mass_msun=star.mass_msun,
                        radius_km=star.radius_km,construction_seconds=bg_s),
        settings=dict(ell=settings.ell,rtol=settings.rtol,atol=settings.atol,
                      max_step_h=settings.max_step_h,root_xtol_hz=settings.root_xtol_hz),
        cases={})
    for method in ("adams","bdf"):
        setup=time.perf_counter()
        solver=NativeCowlingDVODE(star,settings,method=method,samples=8193)
        setup_s=time.perf_counter()-setup
        for name,br,nodes in BRACKETS:
            # Warm each implementation; do not include source sampling or compilation.
            baseline=solver.solve_mode(br,profile_points=240,root_backend="python_brentq")
            native=solver.solve_mode(br,profile_points=240,root_backend="fortran")
            if abs(baseline["frequency_hz"]-native["frequency_hz"])>2.e-3:
                raise AssertionError((name,method,"frequency mismatch",baseline["frequency_hz"],native["frequency_hz"]))
            if (native["W_nodes"],native["V_nodes"])!=(nodes,nodes):
                raise AssertionError((name,method,"wrong node identity",native["W_nodes"],native["V_nodes"]))
            if native["relative_join_mismatch"]>=2.e-5:
                raise AssertionError((name,method,"join",native["relative_join_mismatch"]))
            samples={"python_brentq":[],"fortran":[]}
            for i in range(7):
                order=("python_brentq","fortran") if i%2==0 else ("fortran","python_brentq")
                for backend in order:
                    begin=time.perf_counter()
                    result=solver.solve_mode(br,profile_points=240,root_backend=backend)
                    samples[backend].append(time.perf_counter()-begin)
                    if abs(result["frequency_hz"]-baseline["frequency_hz"])>2.e-3:
                        raise AssertionError((name,method,i,backend))
            median={k:statistics.median(v) for k,v in samples.items()}
            entry=dict(method=method,mode=name,bracket_hz=list(br),
                frontend_setup_seconds=setup_s,
                old_root=dict(frequency_hz=baseline["frequency_hz"],
                              residual=baseline["matching_residual"],
                              iterations=baseline["root_iterations"]),
                native_root=dict(frequency_hz=native["frequency_hz"],
                                 residual=native["matching_residual"],
                                 iterations=native["root_iterations"],
                                 join=native["relative_join_mismatch"],
                                 W_nodes=native["W_nodes"],V_nodes=native["V_nodes"]),
                absolute_frequency_difference_hz=abs(baseline["frequency_hz"]-native["frequency_hz"]),
                full_mode_wall_runs_seconds=samples,
                medians_seconds=median,
                warm_speedup_ratio=median["python_brentq"]/median["fortran"])
            record["cases"][f"{method}_{name}"]=entry
            print(f"{method} {name}: native_f={native['frequency_hz']:.9f} Hz,"
                  f" |delta|={entry['absolute_frequency_difference_hz']:.6g} Hz,"
                  f" speedup={entry['warm_speedup_ratio']:.3f}x",flush=True)
            OUT.parent.mkdir(parents=True,exist_ok=True)
            OUT.write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    record["status"]="PASS_SCOPED_REAL_STAR_ROOT_NUMERIC_PARITY"
    OUT.write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    print("C10_BENCHMARK_WRITTEN",OUT)


if __name__=="__main__":
    run()
