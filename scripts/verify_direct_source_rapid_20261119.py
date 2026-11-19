#!/usr/bin/env python3
"""Verify rapid-conversion exact-source native shooting against DOP853.

A separate physical continuity condition [xi^r - DeltaP/P'] = 0 is used
by the existing unchanged native transfer. New code only replaces the
coefficient interpolation with direct one-sided EOS state evaluation.
"""
from __future__ import annotations
import io,zipfile,json,time,hashlib,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background,solve_mode
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_cowling import CowlingSettings
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
OUT=ROOT/"data/cowling_direct_source_callback_20261119/RAPID.json"

def main():
    with zipfile.ZipFile(SOURCE) as z:
        with zipfile.ZipFile(io.BytesIO(z.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as a:
            contract=json.loads(a.read("FRESH_5x2_WHOLE_STARS/SPG_M2/maxwell/REPRESENTATION.json"))
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(contract))
    n=.30385128926615257
    setting=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,max_step_h=.002)
    native_bg=solve_sharp_background_native(eos,n,surface_h=setting.surface_h,
                  rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
    bracket=(5110.,5170.)
    result=dict(schema="ns_cowling.direct_source_rapid_conversion.v1",
       derived_EOS_bundle_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
       source_model="SPG_M2_PK1_MAXWELL",n_c_fm3=n,
       conversion="rapid",junction_conditions_unchanged=True,
       modes={},oracle={})
    for backend,samples in [("pchip",8193),("direct_source",129)]:
        t=time.monotonic()
        try:
            sol=NativeSharpCowlingDVODE(native_bg,setting,
                                        samples=samples,coefficient_backend=backend)
            mode=sol.solve_mode(bracket,conversion="rapid",profile_points=350)
            result["modes"][backend]=dict(status="ADMITTED",
                samples=samples,frequency_hz=mode["frequency_hz"],
                W_nodes=mode["W_nodes_per_phase"],V_nodes=mode["V_nodes_per_phase"],
                join=mode["relative_join_mismatch"],
                traction_residual=max(x["traction_relative_residual"]
                                     for x in mode["interface_junctions"]),
                kinematic_residual=max(x["kinematic_relative_residual"]
                                      for x in mode["interface_junctions"]),
                callback_rhs_count=mode["source_callback_rhs_calls"],
                elapsed_seconds=time.monotonic()-t)
        except Exception as ex:
            result["modes"][backend]=dict(status="REJECTED",
                error=str(ex)[:250],elapsed_seconds=time.monotonic()-t)
        print("RAPID_NATIVE",backend,result["modes"][backend],flush=True)
    t=time.monotonic()
    py_bg=solve_sharp_background(eos,n,surface_h=setting.surface_h,
                                  rtol=3e-10,atol=3e-12,max_step_h=.003)
    oracle=solve_mode(py_bg,bracket,setting,conversion="rapid",profile_points=350)
    result["oracle"]=dict(frequency_hz=oracle["frequency_hz"],
            join=oracle["relative_join_mismatch"],
            mass_Msun=py_bg.mass_msun,
            elapsed_seconds=time.monotonic()-t)
    print("RAPID_ORACLE",result["oracle"],flush=True)
    result["status"]="ALL_RAPID_CONVERSION_REFERENCE_CHECKS_RECORDED"
    OUT.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("RAPID_RESULT",OUT,flush=True)
if __name__=="__main__":main()
