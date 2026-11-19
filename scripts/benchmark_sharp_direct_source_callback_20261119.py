#!/usr/bin/env python3
"""Native exact-EOS callback versus old PCHIP and pinned DOP853 anchor modes.

Same source-derived EOS, same TOV native background, the same mode shooting
brackets, rtol/atol, interface transfer, and root admission. Vary ONLY
coefficient evaluation backend, keeping the original phase-boundary policy.
"""
from __future__ import annotations
import hashlib,io,json,time,zipfile,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp_native import (
    NativeSharpCowlingDVODE,solve_sharp_background_native)
from compose_eos.stellar_cowling import CowlingSettings

SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
REFERENCE=ROOT/"data/cowling_certified_modes_20261119/RESULT.json"
OUT=ROOT/"data/cowling_direct_source_callback_20261119"
CASES=[
    dict(model="SPG_M3",method="css",nc=.192510767556007,
         mass=.7,conversion="slow",name="p1",
         bracket=[3710.,3760.],ref_key="SPG_M3_CSS_0p7_slow_p1"),
    dict(model="SPG_M3",method="css",nc=.304225985123765,
         mass=1.4,conversion="slow",name="W1V2",
         bracket=[5290.,5360.],ref_key="SPG_M3_CSS_1p4_slow_W1V2"),
    dict(model="SPG_M2",method="maxwell",nc=.30385128926615257,
         mass=1.4,conversion="slow",name="p1",
         bracket=[5130.,5160.],ref_key="SPG_M2_MAXWELL_1p4_slow_p1"),
]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")

def eos_data():
    lookup={}
    with zipfile.ZipFile(SOURCE) as outer:
        b=outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip")
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            for case in CASES:
                pair=(case["model"],case["method"])
                if pair in lookup:continue
                key=f"FRESH_5x2_WHOLE_STARS/{pair[0]}/{pair[1]}/REPRESENTATION.json"
                raw=z.read(key)
                lookup[pair]=dict(sha=hashlib.sha256(raw).hexdigest(),
                    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(
                        json.loads(raw))))
    return lookup

def main():
    lookup=eos_data()
    ref=json.loads(REFERENCE.read_text())["modes"]
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,max_step_h=.002)
    out=dict(schema="ns_cowling.direct_source_vs_PCHIP_vs_DOP853.v1",
        science_parent_commit="daf5de300b06009f8654f795aa99fe40d512666e",
        original_derived_eos_zip_sha256=sha(SOURCE),
        original_reference_sha256=sha(REFERENCE),
        phasewise_sharp_interfaces_unchanged=True,
        no_new_density_interpolation_between_phases=True,
        test_settings=vars(settings),
        cases={},status="RUNNING")
    for case in CASES:
        eos=lookup[(case["model"],case["method"])]["eos"]
        t=time.monotonic()
        bg=solve_sharp_background_native(eos,case["nc"],native_samples=2049,
            surface_h=settings.surface_h,rtol=3e-10,atol=3e-11,max_step_h=.004)
        bg_elapsed=time.monotonic()-t
        expected_freq=ref[case["ref_key"]]["python_reference"]["frequency_hz"]
        entry=dict(reference_case=case,source_contract_sha256=
                   lookup[(case["model"],case["method"])]["sha"],
                   oracle_Python_direct_EOS_frequency_hz=expected_freq,
                   background_native_tov=dict(mass_msun=bg.mass_msun,
                       radius_km=bg.radius_km,elapsed_seconds=bg_elapsed),
                   variants={})
        for backend,samples in [("pchip",8193),("direct_source",129),
                                ("pchip",16385)]:
            key=f"{backend}_{samples}"
            start=time.monotonic()
            try:
                t0=time.monotonic()
                solver=NativeSharpCowlingDVODE(bg,settings,samples=samples,
                                             coefficient_backend=backend)
                setup=time.monotonic()-t0
                root=solver.solve_mode(case["bracket"],conversion=case["conversion"],
                                       profile_points=450)
                trial=dict(status="NATIVE_ROOT_ADMITTED",
                    frequency_hz=root["frequency_hz"],
                    delta_against_DOP853_hz=root["frequency_hz"]-expected_freq,
                    relative_join_mismatch=root["relative_join_mismatch"],
                    nodes_W=root["W_nodes_per_phase"],nodes_V=root["V_nodes_per_phase"],
                    source_callback_rhs_calls=root["source_callback_rhs_calls"],
                    tabulation_points_per_phase=[len(ph["s"]) for ph in solver._phases],
                    setup_elapsed_seconds=setup,
                    total_elapsed_seconds=time.monotonic()-start)
            except Exception as exc:
                trial=dict(status="NATIVE_ROOT_REJECTED",
                     reason=type(exc).__name__+": "+str(exc)[:360],
                     total_elapsed_seconds=time.monotonic()-start)
            entry["variants"][key]=trial
            print("COMPARE",case["ref_key"],key,
                "outcome",trial["status"],"f",trial.get("frequency_hz"),
                "deltaHz",trial.get("delta_against_DOP853_hz"),
                "seconds",round(trial["total_elapsed_seconds"],3),flush=True)
        out["cases"][case["ref_key"]]=entry
        save(OUT/"PARTIAL.json",out)
    out["status"]="BOUNDED_SOURCE_CALLBACK_COMPARISON_RECORDED"
    save(OUT/"RESULT.json",out)
    print("SOURCE_CALLBACK_COMPARISON_DONE",len(out["cases"]),flush=True)

if __name__=="__main__":main()
