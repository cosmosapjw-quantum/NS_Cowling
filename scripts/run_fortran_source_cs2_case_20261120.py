#!/usr/bin/env python3
"""One independent original-source star and Fortran c_s² mode calculation.

Python: source parameter marshaling and result formatting.
Fortran: TOV and Cowling DVODE, phase transfer, source c_s² every RHS.
The existing sharp-mode Brent outer iteration is still Python.
"""
from __future__ import annotations
import argparse,hashlib,io,json,os,sys,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.stellar_cowling import CowlingSettings
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
CONFIG={
  "M2":("SPG_M2","maxwell",.30385128926615257,"rapid",(5110.,5170.),5143.719721563306),
  "M3":("SPG_M3","css",.192510767556007,"slow",(3710.,3760.),3736.2088951581254)
}

def run(case,out):
    src=SOURCE.read_bytes()
    if hashlib.sha256(src).hexdigest()!=SHA: raise RuntimeError("source EOS identity mismatch")
    model,join,nc,conversion,bracket,reference=CONFIG[case]
    with zipfile.ZipFile(io.BytesIO(src)) as archive:
        with zipfile.ZipFile(io.BytesIO(archive.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as nested:
            raw=nested.read(f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
    star=solve_sharp_background_native(eos,nc,surface_h=3e-8,
       rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,max_step_h=.002)
    native=NativeSharpCowlingDVODE(star,settings,samples=4097,
       coefficient_backend="native_source_cs2",source_restart_relative_jump=.01)
    m=native.solve_mode(bracket,conversion=conversion,profile_points=300)
    if (m["W_nodes_per_phase"],m["V_nodes_per_phase"])!=(1,1):
        raise RuntimeError("physical p1 mode node contract failed")
    if m["relative_join_mismatch"]>=1e-5 or abs(m["frequency_hz"]-reference)>=.05:
        raise RuntimeError("physical p1 source parity not satisfied")
    if m["source_callback_rhs_calls"]!=0:
        raise RuntimeError("Python callback re-entered RHS")
    row=dict(schema="ns_cowling.real_source_native_cs2_process_case.v1",
       eos_archive_sha256=SHA,model=model,join=join,
       source_contract_sha256=hashlib.sha256(raw).hexdigest(),
       n_c_fm3=nc,mass_msun=star.mass_msun,radius_km=star.radius_km,
       conversion=conversion,
       frequency_hz=m["frequency_hz"],reference_hz=reference,
       residual=m["matching_residual"],join_relative=m["relative_join_mismatch"],
       nodes=[m["W_nodes_per_phase"],m["V_nodes_per_phase"]],
       source_restart_count=m["source_restart_count"],
       source_callback_rhs_calls=m["source_callback_rhs_calls"],
       native_dvode_nfev=m["nfev"],
       code_fortran_eos_source_eval=True,
       root_polishing_frontend="python_brentq_legacy",
       scientific_status="PASS_SCOPED_ORIGINAL_EOS_ONEJUMP")
    out.parent.mkdir(parents=True,exist_ok=True)
    tmp=out.with_suffix(".json.part")
    tmp.write_text(json.dumps(row,indent=2,allow_nan=False,sort_keys=True)+"\n")
    os.replace(tmp,out)
    print("SOURCE_NATIVE_CASE",case,conversion,
          "frequency",m["frequency_hz"],"calls",m["source_callback_rhs_calls"],flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--case",required=True,choices=tuple(CONFIG))
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    run(args.case,args.output)
