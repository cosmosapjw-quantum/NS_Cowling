#!/usr/bin/env python3
"""Real source-locked sharp-star comparison of Python-callback vs Fortran c_s².

Both backends execute identical native Fortran DVODE shooting, geometry PCHIP,
strong-knot restarts, Maxwell transfer and SciPy Brent root search. Only the
original-source c_s² evaluator changes from Python callback to Fortran.
"""
from __future__ import annotations
import hashlib,io,json,statistics,sys,time,zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.stellar_cowling import CowlingSettings

SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SOURCE_SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
OUT=ROOT/"data/cowling_fortran_source_cs2_20261120/REAL_STAR_PERFORMANCE.json"
CASES=(
    ("M2_MAXWELL_RAPID_P1","SPG_M2","maxwell",.30385128926615257,
     "rapid",(5110.,5170.),5143.719721563306),
    ("M3_CSS_SLOW_P1","SPG_M3","css",.192510767556007,
     "slow",(3710.,3760.),3736.2088951581254))

def run():
    raw=SOURCE.read_bytes()
    assert hashlib.sha256(raw).hexdigest()==SOURCE_SHA
    out=dict(schema="ns_cowling.native_exact_source_cs2_callback_replacement.v1",
        parent_commit="19e5b16cd0c5542e8fcb91f2fab14f36b1e663df",
        source_bundle_sha256=SOURCE_SHA,
        closure="GAMMA1_EQUALS_GAMMA_EQ",
        frozen_composition_supported=False,
        numerical_comparison="same source EOS and geometry, Fortran DVODE, knot restarts; Python SciPy Brent root remains",
        samples_requested=4097,
        cases={})
    with zipfile.ZipFile(io.BytesIO(raw)) as outer:
        with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as nested:
            for name,model,join,density,conversion,bracket,ref in CASES:
                path=f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json"
                blob=nested.read(path)
                eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(blob)))
                star=solve_sharp_background_native(eos,density,surface_h=3e-8,
                    rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
                settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,max_step_h=.002)
                r=dict(model=model,join=join,conversion=conversion,bracket_hz=list(bracket),
                       original_contract_sha256=hashlib.sha256(blob).hexdigest(),
                       density_fm3=density,mass_msun=float(star.mass_msun),
                       reference_exact_source_DOP853_hz=ref,variants={})
                for backend in ("hybrid_cs2","native_source_cs2"):
                    start=time.perf_counter()
                    solver=NativeSharpCowlingDVODE(star,settings,samples=4097,
                      coefficient_backend=backend,source_restart_relative_jump=.01)
                    setup_s=time.perf_counter()-start
                    # Warm mode once, then three bounded timings.
                    mode=solver.solve_mode(bracket,conversion=conversion,profile_points=300)
                    runs=[]
                    for _ in range(3):
                        t=time.perf_counter()
                        mode=solver.solve_mode(bracket,conversion=conversion,profile_points=300)
                        runs.append(time.perf_counter()-t)
                    assert (mode["W_nodes_per_phase"],mode["V_nodes_per_phase"])==(1,1)
                    assert mode["relative_join_mismatch"]<1e-5
                    assert abs(mode["frequency_hz"]-ref)<.05
                    r["variants"][backend]=dict(
                        frequency_hz=float(mode["frequency_hz"]),
                        frequency_minus_external_reference_hz=float(mode["frequency_hz"]-ref),
                        nfev=int(mode["nfev"]),
                        node_counts=[int(mode["W_nodes_per_phase"]),int(mode["V_nodes_per_phase"])],
                        root_iterations=int(mode["root_iterations"]),
                        join_relative=float(mode["relative_join_mismatch"]),
                        interface_traction_residual_max=max(
                          j["traction_relative_residual"] for j in mode["interface_junctions"]),
                        source_callback_rhs_calls=int(mode["source_callback_rhs_calls"]),
                        restart_count=int(mode["source_restart_count"]),
                        setup_seconds=setup_s,
                        warm_mode_runs_seconds=runs,
                        median_seconds=statistics.median(runs))
                    print("MODE",name,backend,"F_HZ",mode["frequency_hz"],
                      "WARM_MEDIAN_S",statistics.median(runs),
                      "PYTHON_RHS_CALLS",mode["source_callback_rhs_calls"],flush=True)
                    OUT.parent.mkdir(exist_ok=True,parents=True)
                    out["cases"][name]=r
                    OUT.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
                old=r["variants"]["hybrid_cs2"]; new=r["variants"]["native_source_cs2"]
                delta=new["frequency_hz"]-old["frequency_hz"]
                assert abs(delta)<.005,(name,delta)
                assert new["source_callback_rhs_calls"]==0
                assert new["restart_count"]==old["restart_count"]
                r["frequency_difference_native_minus_python_callback_hz"]=delta
                r["warm_speedup_python_callback_over_fortran_source"]=old["median_seconds"]/new["median_seconds"]
                print("COMPARISON",name,"DELTA_HZ",delta,"SPEEDUP",
                    r["warm_speedup_python_callback_over_fortran_source"],flush=True)
    out["status"]="PASS_SCOPED_SOURCE_EXACT_CS2_AND_REAL_STAR_MODE_PARITY"
    OUT.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print("SOURCE_CS2_PERFORMANCE_COMPARISON_SAVED",OUT,flush=True)

if __name__=="__main__":
    run()
