#!/usr/bin/env python3
"""Warm original-EOS sharp Maxwell Brent vs new all-Fortran root timing.

Python creates the immutable source contract, TOV and phasewise PCHIP arrays.
No new Python root computation; old SciPy brentq is the explicit comparison.
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
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
OUT=ROOT/"data/cowling_sharp_native_root_20261120/PCHIP_ROOT_BENCHMARK.json"
CASES=(("slow_f","slow",(1800.,1900.),(0,0)),
       ("rapid_f","rapid",(1800.,1900.),(0,0)),
       ("slow_p1","slow",(5100.,5200.),(1,1)),
       ("rapid_p1","rapid",(5110.,5170.),(1,1)))

def main():
    raw=SOURCE.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SHA:raise RuntimeError("source bundle hash mismatch")
    with zipfile.ZipFile(io.BytesIO(raw)) as outer:
        with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inner:
            payload=inner.read("FRESH_5x2_WHOLE_STARS/SPG_M2/maxwell/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(payload)))
    t=time.perf_counter()
    bg=solve_sharp_background_native(eos,.30385128926615257,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
    background_s=time.perf_counter()-t
    st=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,max_step_h=.002,root_xtol_hz=1e-6)
    t=time.perf_counter()
    solver=NativeSharpCowlingDVODE(bg,st,samples=8193,coefficient_backend="pchip")
    coeff_s=time.perf_counter()-t
    rec=dict(schema="ns_cowling.sharp_onejump_native_root_benchmark.v1",
        parent_commit="436d7e4213d94696cf07defe075fed778587d12a",
        source_bundle_sha256=SHA,source_contract_sha256=hashlib.sha256(payload).hexdigest(),
        model="SPG_M2_maxwell",central_density_fm3=.30385128926615257,
        background=dict(mass_msun=bg.mass_msun,radius_km=bg.radius_km,
                        phase_count=len(solver._phases),source_tov_seconds=background_s,
                        coefficient_preparation_seconds=coeff_s,
                        coefficient_backend="PCHIP_ONE_SIDED_NOT_EXACT_SOURCE"),
        settings=dict(rtol=st.rtol,atol=st.atol,max_step_h=st.max_step_h,
                      root_xtol_hz=st.root_xtol_hz),
        cases={})
    for label,conversion,bracket,nodes in CASES:
        samples={"python_brentq":[],"fortran_phase_pchip":[]}
        example={}
        # Same star, phase tables, physical controls and profile size.
        for backend in samples:
            example[backend]=solver.solve_mode(bracket,conversion=conversion,
                profile_points=200,root_backend=backend)
        for run in range(7):
            order=("fortran_phase_pchip","python_brentq") if run%2 else ("python_brentq","fortran_phase_pchip")
            for backend in order:
                start=time.perf_counter()
                result=solver.solve_mode(bracket,conversion=conversion,
                    profile_points=200,root_backend=backend)
                samples[backend].append(time.perf_counter()-start)
                if abs(result["frequency_hz"]-example["python_brentq"]["frequency_hz"])>.005:
                    raise RuntimeError("native root frequency drift: "+label)
                if (result["W_nodes_per_phase"],result["V_nodes_per_phase"])!=nodes:
                    raise RuntimeError("p/f physical node-class mismatch: "+label)
        baseline=example["python_brentq"]; native=example["fortran_phase_pchip"]
        ms={k:statistics.median(v) for k,v in samples.items()}
        rec["cases"][label]=dict(
            conversion=conversion,bracket_hz=list(bracket),expected_nodes=list(nodes),
            baseline_frequency_hz=baseline["frequency_hz"],
            native_frequency_hz=native["frequency_hz"],
            delta_frequency_hz=native["frequency_hz"]-baseline["frequency_hz"],
            baseline_join=baseline["relative_join_mismatch"],
            native_join=native["relative_join_mismatch"],
            native_interface_traction=[j["traction_relative_residual"] for j in native["interface_junctions"]],
            native_native_rhs_calls=native["root_native_rhs_calls"],
            native_root_iterations=native["root_iterations"],
            elapsed_s=samples,median_elapsed_s=ms,
            speedup_python_brent_vs_fortran=ms["python_brentq"]/ms["fortran_phase_pchip"])
        print(label,"delta_Hz",rec["cases"][label]["delta_frequency_hz"],
          "medians",ms,"speedup",rec["cases"][label]["speedup_python_brent_vs_fortran"],flush=True)
        OUT.write_text(json.dumps(rec,indent=2,allow_nan=False)+"\n")
    rec["status"]="SCOPED_SOURCE_LOCKED_ONEJUMP_NUMERICAL_PARITY"
    OUT.write_text(json.dumps(rec,indent=2,allow_nan=False)+"\n")
    print("PCHIP_SHARP_BENCHMARK_DONE",OUT,flush=True)
if __name__=="__main__":main()
