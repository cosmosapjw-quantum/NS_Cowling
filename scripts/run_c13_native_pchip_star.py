#!/usr/bin/env python3
"""C13 frontend-only independent real-source star process.

EOS identity, one-sided TOV sample preparation and report serialization
are Python frontend. Native Fortran owns geometry PCHIP, source c_s²,
strong-knot scheduler, DVODE, Maxwell junction and Brent polishing.
"""
from __future__ import annotations
import argparse,hashlib,io,json,os,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.stellar_cowling import CowlingSettings

SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
CASES={
"M2":("SPG_M2","maxwell",.30385128926615257,"rapid",(5110.,5170.)),
"M3":("SPG_M3","css",.192510767556007,"slow",(3710.,3760.))
}

def run(case,out,samples):
    original=SOURCE.read_bytes()
    if hashlib.sha256(original).hexdigest()!=SHA:raise RuntimeError("source archive hash changed")
    model,join,nc,conversion,br=CASES[case]
    with zipfile.ZipFile(io.BytesIO(original)) as full:
        with zipfile.ZipFile(io.BytesIO(full.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inside:
            raw=inside.read(f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
    t=time.perf_counter()
    star=solve_sharp_background_native(eos,nc,surface_h=3e-8,
       rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
    background_s=time.perf_counter()-t
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,
       atol=1e-12,max_step_h=.002,root_xtol_hz=1e-5)
    t=time.perf_counter()
    solver=NativeSharpCowlingDVODE(star,settings,samples=samples,
       coefficient_backend="native_source_cs2",source_restart_relative_jump=.01)
    geometry_s=time.perf_counter()-t
    t=time.perf_counter()
    mode=solver.solve_mode(br,conversion=conversion,profile_points=300,
       root_backend="fortran_source_exact")
    mode_s=time.perf_counter()-t
    assert mode["geometry_interpolation_backend"]=="fortran_pchip"
    assert mode["strong_knot_selector_backend"]=="fortran_source_exact"
    assert mode["source_callback_rhs_calls"]==0
    assert (mode["W_nodes_per_phase"],mode["V_nodes_per_phase"])==(1,1)
    assert mode["relative_join_mismatch"]<1e-5
    assert mode["root_fortran_restart_count"]==solver.source_restart_count
    r=dict(schema="ns_cowling.c13.fortran_pchip_original_star.v1",
      source_sha256=SHA,source_contract_sha256=hashlib.sha256(raw).hexdigest(),
      original_model=model,join=join,central_density_fm3=nc,
      conversion=conversion,mode="p1",samples_requested=samples,
      mass_msun=float(star.mass_msun),radius_km=float(star.radius_km),
      frequency_hz=float(mode["frequency_hz"]),
      original_ode_backend=mode["backend"],geometry_backend="fortran_pchip",
      source_sound_speed_backend="fortran_original_contract",
      root_backend=mode["root_backend"],source_restart_backend="fortran",
      native_rhs_evals=int(mode["root_native_rhs_calls"]),
      native_source_restarts=int(mode["root_fortran_restart_count"]),
      W_nodes=int(mode["W_nodes_per_phase"]),V_nodes=int(mode["V_nodes_per_phase"]),
      relative_join_mismatch=float(mode["relative_join_mismatch"]),
      normalized_determinant=float(mode["matching_residual"]),
      python_rhs_eos_callback_count=int(mode["source_callback_rhs_calls"]),
      validation="PASS_SCOPED_SOURCE_LOCKED_P1")
    out.parent.mkdir(parents=True,exist_ok=True)
    partial=out.with_name(out.name+".part")
    partial.write_text(json.dumps(r,sort_keys=True,indent=2,allow_nan=False)+"\n")
    os.replace(partial,out)
    print("C13_FORTRAN_STAR",case,"P1_HZ",r["frequency_hz"],
      "source_calls",r["python_rhs_eos_callback_count"],
      "FRONTEND_BG_SECONDS",background_s,
      "FRONTEND_GEOMETRY_SECONDS",geometry_s,
      "NATIVE_MODE_SECONDS",mode_s,flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--case",choices=CASES,required=True)
    parser.add_argument("--samples",type=int,default=4097)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    run(args.case,args.output,args.samples)
