#!/usr/bin/env python3
"""Bounded direct-source DVODE numerical sensitivity, same EOS and TOV.

Tests whether tighter DVODE tolerances/maximum source step change f in the
challenging M3 CSS 0.7 Msun slow p1. These are local solver sensitivity
checks, not a formal global numerical error bound.
"""
from __future__ import annotations
import io,zipfile,json,hashlib,time,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (
    solve_sharp_background_native,NativeSharpCowlingDVODE)
from compose_eos.stellar_cowling import CowlingSettings
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
OUT=ROOT/"data/cowling_direct_source_callback_20261119/TOLERANCE.json"

def main():
    with zipfile.ZipFile(SOURCE) as z:
        with zipfile.ZipFile(io.BytesIO(z.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as a:
            raw=a.read("FRESH_5x2_WHOLE_STARS/SPG_M3/css/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
    bg=solve_sharp_background_native(eos,.192510767556007,
                  surface_h=3e-8,rtol=3e-10,atol=3e-11,
                  max_step_h=.004,native_samples=2049)
    report=dict(schema="ns_cowling.native_direct_source_rhs_tolerance.v1",
        original_contract_sha256=hashlib.sha256(raw).hexdigest(),
        original_archive_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        case="SPG_M3_CSS_0p7Msun_slow_p1",
        fixed_native_TOV_mass_msun=bg.mass_msun,
        independent_Python_DOP853_source_reference_hz=3736.2088951581254,
        parameter_variation="DVODE mode integrator rtol, atol and maximum step only",
        trials=[],status="RUNNING")
    for rtol,maxstep in ((1e-8,.002),(1e-9,.001),(1e-10,.0005)):
        setting=CowlingSettings(surface_h=3e-8,rtol=rtol,
              atol=rtol*1e-3,max_step_h=maxstep)
        t=time.monotonic()
        try:
            mode=NativeSharpCowlingDVODE(bg,setting,samples=129,
                  coefficient_backend="direct_source").solve_mode(
                    (3710.,3760.),conversion="slow",profile_points=450)
            record=dict(status="ADMITTED",rtol=rtol,max_step_h=maxstep,
                frequency_hz=float(mode["frequency_hz"]),
                DOP853_delta_hz=float(mode["frequency_hz"]-
                                 report["independent_Python_DOP853_source_reference_hz"]),
                join=float(mode["relative_join_mismatch"]),
                W_nodes=int(mode["W_nodes_per_phase"]),
                V_nodes=int(mode["V_nodes_per_phase"]),
                source_rhs_evaluations=int(mode["source_callback_rhs_calls"]),
                elapsed_seconds=time.monotonic()-t)
        except Exception as exc:
            record=dict(status="REJECTED",rtol=rtol,max_step_h=maxstep,
                reason=str(exc)[:310],elapsed_seconds=time.monotonic()-t)
        report["trials"].append(record)
        print("DIRECT_SOURCE_TOL",record,flush=True)
    report["status"]="ALL_THREE_SETTINGS_ATTEMPTED"
    OUT.write_text(json.dumps(report,indent=2)+"\n")
    print("DIRECT_SOURCE_TOL_DONE",len(report["trials"]),flush=True)

if __name__=="__main__":
    main()
