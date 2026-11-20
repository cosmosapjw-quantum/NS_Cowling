#!/usr/bin/env python3
"""C12: matched-accuracy real-source M2/M3 Fortran sharp-root benchmark.

Root iterations, original source c_s² and strong-knot restarts execute
inside Fortran for the new backend; Python is source/frontend only.
SciPy Brent on same native DVODE supplies a reference comparison.
"""
from __future__ import annotations
import hashlib,io,json,statistics,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (
  NativeSharpCowlingDVODE,solve_sharp_background_native)
from compose_eos.stellar_cowling import CowlingSettings

SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
OUT=ROOT/"data/cowling_fortran_source_exact_root_c12_20261121/REAL_SOURCE_BENCHMARK.json"
CASES=(
 ("M2_MAXWELL_1P4_RAPID_P1","SPG_M2","maxwell",.30385128926615257,"rapid",(5110.,5170.),5143.719721563306),
 ("M3_CSS_0P7_SLOW_P1","SPG_M3","css",.192510767556007,"slow",(3710.,3760.),3736.2088951581254),
)
BACKENDS=("python_brentq","fortran_source_exact")

def run():
    raw=SOURCE.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SHA:
        raise RuntimeError("source archive drift")
    record=dict(schema="ns_cowling.c12.source_exact_native_brent_benchmark.v1",
       baseline_commit="7b6bebee997e2a3cad069bddcd8580ea36f63839",
       source_bundle_sha256=SHA,
       scientific_scope="one sharp Maxwell interface, Gamma1=Gamma_eq, source-exact cs2",
       backend_comparison="same original EOS, Fortran TOV and DVODE, same phase geometry, strong-knot set",
       method="adams",tolerances=dict(rtol=1e-9,atol=1e-12,
             max_step_h=.002,root_xtol_hz=1e-5),
       repeated_run_count=7,cases={})
    with zipfile.ZipFile(io.BytesIO(raw)) as outer:
      with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as source:
        for name,model,join,nc,conversion,br,ref in CASES:
          source_path=f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json"
          original=source.read(source_path)
          eos=SharpContractEnthalpyEOS(
            WholeCrustBarotrope.from_contract(json.loads(original)))
          t=time.perf_counter()
          star=solve_sharp_background_native(eos,nc,surface_h=3e-8,
            rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
          background_seconds=time.perf_counter()-t
          settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,
            max_step_h=.002,root_xtol_hz=1e-5)
          t=time.perf_counter()
          solver=NativeSharpCowlingDVODE(star,settings,samples=4097,
             coefficient_backend="native_source_cs2",source_restart_relative_jump=.01)
          setup_seconds=time.perf_counter()-t
          samples={name:[] for name in BACKENDS}
          examples={name:solver.solve_mode(br,conversion=conversion,
                    profile_points=300,root_backend=name) for name in BACKENDS}
          baseline=examples["python_brentq"]
          native=examples["fortran_source_exact"]
          assert abs(native["frequency_hz"]-baseline["frequency_hz"])<.005
          assert abs(native["frequency_hz"]-ref)<.05
          assert native["root_fortran_restart_count"]==solver.source_restart_count
          assert native["source_callback_rhs_calls"]==0
          assert (native["W_nodes_per_phase"],native["V_nodes_per_phase"])==(1,1)
          assert native["relative_join_mismatch"]<1e-5
          for j in range(7):
            order=BACKENDS if j%2==0 else BACKENDS[::-1]
            for backend in order:
              start=time.perf_counter()
              m=solver.solve_mode(br,conversion=conversion,
                      profile_points=300,root_backend=backend)
              elapsed=time.perf_counter()-start
              samples[backend].append(elapsed)
              assert abs(m["frequency_hz"]-baseline["frequency_hz"])<.005
              assert (m["W_nodes_per_phase"],m["V_nodes_per_phase"])==(1,1)
          median={k:statistics.median(v) for k,v in samples.items()}
          record["cases"][name]=dict(model=model,join=join,conversion=conversion,
             source_contract_sha256=hashlib.sha256(original).hexdigest(),
             density_fm3=nc,mass_msun=star.mass_msun,
             radius_km=star.radius_km,background_seconds=background_seconds,
             coefficient_setup_seconds=setup_seconds,
             bracket_hz=list(br),external_reference_hz=ref,
             original_py_brent_freq_hz=baseline["frequency_hz"],
             fortran_brent_freq_hz=native["frequency_hz"],
             frequency_difference_hz=native["frequency_hz"]-baseline["frequency_hz"],
             fortran_residual=native["matching_residual"],
             fortran_join_relative=native["relative_join_mismatch"],
             fortran_root_iterations=native["root_iterations"],
             fortran_root_nfev=native["root_native_rhs_calls"],
             root_restart_cuts_in_fortran=native["root_fortran_restart_count"],
             root_restart_cuts_previous_python=solver.source_restart_count,
             python_rhs_callbacks=native["source_callback_rhs_calls"],
             nodes=[native["W_nodes_per_phase"],native["V_nodes_per_phase"]],
             wall_seconds=samples,median_wall_seconds=median,
             median_speedup=median["python_brentq"]/median["fortran_source_exact"])
          OUT.parent.mkdir(parents=True,exist_ok=True)
          OUT.write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
          print("BENCH",name,"PY",baseline["frequency_hz"],
                "NATIVE",native["frequency_hz"],
                "DIFF",native["frequency_hz"]-baseline["frequency_hz"],
                "SPEEDUP",record["cases"][name]["median_speedup"],
                "FORTRAN_RHS",native["root_native_rhs_calls"],flush=True)
    record["status"]="PASS_SCOPED_MATCHED_REAL_SOURCE_MODES"
    OUT.write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    print("RECORD",OUT,flush=True)

if __name__=="__main__":run()
