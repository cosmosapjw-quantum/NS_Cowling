#!/usr/bin/env python3
"""Native Cowling PCHIP vs source callback vs hybrid on same star and frequency.

All source-EOS contracts are extracted from the pinned 05 nested archive.
Independent original-EOS DOP853 frequencies are already checked and pinned
in the parent scientific stage, not recomputed per variant.
Runtimes are end-to-end native setup + root shot, not matched single RHS.
"""
from __future__ import annotations
import hashlib,io,json,sys,time,traceback,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (
    NativeSharpCowlingDVODE,solve_sharp_background_native)
from compose_eos.stellar_cowling import CowlingSettings

BUNDLE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
DATA=ROOT/"data/cowling_hybrid_geometry_20261120"
PARENT=ROOT/"data/cowling_direct_source_callback_20261119"
CASES=[
 dict(id="M3_CSS_0p7_SLOW_P1",model="SPG_M3",join="css",
      n_c_fm3=.192510767556007,conversion="slow",
      bracket=[3710.,3760.],reference_hz=3736.2088951581254,
      expected_nodes=[1,1]),
 dict(id="M2_MAXWELL_1p4_RAPID_P1",model="SPG_M2",join="maxwell",
      n_c_fm3=.30385128926615257,conversion="rapid",
      bracket=[5110.,5170.],reference_hz=5143.719721563306,
      expected_nodes=[1,1]),
]
VARIANTS=[
 ("PCHIP_8193","pchip",8193),
 ("PCHIP_16385","pchip",16385),
 ("PCHIP_32769","pchip",32769),
 ("DIRECT_SOURCE_129","direct_source",129),
 ("HYBRID_CS2_2049","hybrid_cs2",2049),
 ("HYBRID_CS2_4097","hybrid_cs2",4097),
]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,obj):
    path.parent.mkdir(exist_ok=True,parents=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")
def load_eos(bundle_bytes):
    out={}
    with zipfile.ZipFile(io.BytesIO(bundle_bytes)) as outer:
        nested=outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip")
        with zipfile.ZipFile(io.BytesIO(nested)) as inside:
            hashlist=json.loads(inside.read("INTERNAL_PROVENANCE_AND_SHA256.json"))
            for c in CASES:
                source=f"FRESH_5x2_WHOLE_STARS/{c['model']}/{c['join']}/REPRESENTATION.json"
                b=inside.read(source)
                if hashlib.sha256(b).hexdigest()!=hashlist["file_sha256"][source]:
                    raise RuntimeError("EOS contract SHA-256 changed: "+source)
                out[c["id"]]=dict(eos=SharpContractEnthalpyEOS(
                     WholeCrustBarotrope.from_contract(json.loads(b))),
                     hash=hashlib.sha256(b).hexdigest(),source_path=source)
    return out

def main():
    original=BUNDLE.read_bytes()
    models=load_eos(original)
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,
             atol=1e-12,max_step_h=.002)
    report=dict(schema="ns_cowling.source_hybrid_cs2_native_benchmark.v1",
       baseline_commit="6253dc40eee188ea207450ae13a5d4588aa6d398",
       frozen_EOS_bundle_sha256=sha(BUNDLE),test_settings=vars(settings),
       reference_method="SOURCE_EOS_PYTHON_DOP853_PARITY_PREVIOUS_STAGES",
       oracle_from="data/cowling_direct_source_callback_20261119/RESULT.json plus RAPID.json",
       measurement="single process, OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1",
       tolerance_hz=.05,maximum_join=1e-5,
       cases={},status="RUNNING")
    for case in CASES:
        t0=time.perf_counter()
        base=models[case["id"]]
        background=solve_sharp_background_native(
            base["eos"],case["n_c_fm3"],
            surface_h=settings.surface_h,rtol=3e-10,atol=3e-11,
            max_step_h=.004,native_samples=2049)
        item=dict(model=case["model"],join=case["join"],
            n_c_fm3=case["n_c_fm3"],conversion=case["conversion"],
            reference_hz=case["reference_hz"],
            source_contract_SHA256=base["hash"],
            backend_variants={},
            native_TOV=dict(M_msun=background.mass_msun,
                R_km=background.radius_km,setup_seconds=time.perf_counter()-t0))
        for label,backend,samples in VARIANTS:
            start=time.perf_counter()
            rec=dict(samples=samples,backend=backend,status="RUNNING")
            try:
                solver=NativeSharpCowlingDVODE(
                    background,settings,samples=samples,
                    coefficient_backend=backend,
                    source_restart_relative_jump=.01)
                setup_s=time.perf_counter()-start
                mode=solver.solve_mode(case["bracket"],
                    conversion=case["conversion"],profile_points=350)
                total=time.perf_counter()-start
                df=float(mode["frequency_hz"]-case["reference_hz"])
                nodes=[int(mode["W_nodes_per_phase"]),int(mode["V_nodes_per_phase"])]
                accepted=(abs(df)<=.05 and
                    float(mode["relative_join_mismatch"])<=1e-5 and
                    nodes==case["expected_nodes"])
                rec.update(status="MATCHED_REFERENCE_PASS" if accepted else "NUMERICALLY_UNRESOLVED",
                    frequency_hz=float(mode["frequency_hz"]),
                    reference_delta_hz=df,join_relative=float(mode["relative_join_mismatch"]),
                    nodes=nodes,source_restart_count=solver.source_restart_count,
                    source_callback_calls=mode.get("source_callback_rhs_calls",0),
                    setup_seconds=setup_s,total_seconds=total)
            except Exception as ex:
                rec.update(status="NATIVE_ROOT_REJECTED",
                    error_type=type(ex).__name__,error=str(ex)[:350],
                    total_seconds=time.perf_counter()-start)
            item["backend_variants"][label]=rec
            print("HYBRID_BENCHMARK",case["id"],label,
                rec["status"],"frequency",rec.get("frequency_hz"),
                "df",rec.get("reference_delta_hz"),
                "time",round(rec["total_seconds"],3),flush=True)
            save(DATA/"PARTIAL.json",report)
        direct=item["backend_variants"]["DIRECT_SOURCE_129"]
        hybrids=[(label,item["backend_variants"][label]) for label in
            ("HYBRID_CS2_2049","HYBRID_CS2_4097")
            if item["backend_variants"][label]["status"]=="MATCHED_REFERENCE_PASS"]
        if direct["status"]=="MATCHED_REFERENCE_PASS" and hybrids:
            winner=min(hybrids,key=lambda v:v[1]["total_seconds"])
            item["fastest_hybrid_accepted"]=dict(
                label=winner[0],
                speedup_source_callback_vs_hybrid=
                  direct["total_seconds"]/winner[1]["total_seconds"],
                frequency_error_hz=winner[1]["reference_delta_hz"],
                reference_tolerance_hz=.05)
        else:
            item["fastest_hybrid_accepted"]=None
        accurate=[(k,v) for k,v in item["backend_variants"].items()
             if v["status"]=="MATCHED_REFERENCE_PASS"
             and abs(v["reference_delta_hz"])<=.005]
        if accurate:
            choice=min(accurate,key=lambda pair:pair[1]["total_seconds"])
            item["fastest_matched_5mHz"]=dict(
                label=choice[0],backend=choice[1]["backend"],
                frequency_hz=choice[1]["frequency_hz"],
                signed_reference_error_hz=choice[1]["reference_delta_hz"],
                elapsed_seconds=choice[1]["total_seconds"])
        else:
            item["fastest_matched_5mHz"]=None
        report["cases"][case["id"]]=item
        save(DATA/"PARTIAL.json",report)
    report["status"]="BOUNDED_BENCHMARK_RECORDED_NOT_GENERAL_PERFORMANCE_PROOF"
    save(DATA/"RESULT.json",report)
    print("HYBRID_BENCHMARK_END",report["status"],flush=True)
if __name__=="__main__":
    main()
