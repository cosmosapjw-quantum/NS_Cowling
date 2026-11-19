#!/usr/bin/env python3
"""Quantify PCHIP cs2 smearing at exact source EOS knot neighborhoods.

This *does not* assert physical EOS discontinuities arise from a first
order phase conversion: these are inner spline/polytrope derivative
jumps WITHIN one physical phase. It records enough data to falsify the
source restart criterion used by the optional new Fortran backend.
"""
from __future__ import annotations
import hashlib,io,json,sys,zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (
    solve_sharp_background_native,NativeSharpCowlingDVODE)
from compose_eos.stellar_cowling import CowlingSettings

SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
OUT=ROOT/"data/cowling_direct_source_callback_20261119/KNOT_AUDIT.json"

def main():
    with zipfile.ZipFile(SOURCE) as z:
        with zipfile.ZipFile(io.BytesIO(z.read(
                 "05_SPG5_PK1_10_WHOLE_STARS.zip"))) as d:
            raw=d.read(
              "FRESH_5x2_WHOLE_STARS/SPG_M3/css/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
    bg=solve_sharp_background_native(eos,.192510767556007,
        native_samples=2049,surface_h=3e-8,rtol=3e-10,
        atol=3e-11,max_step_h=.004)
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,
        atol=1e-12,max_step_h=.002)
    pchip=NativeSharpCowlingDVODE(bg,settings,samples=4097)
    direct=NativeSharpCowlingDVODE(bg,settings,samples=129,
                  coefficient_backend="direct_source",
                  source_restart_relative_jump=.01)
    result=dict(schema="ns_cowling.sharp_original_knot_coefficient_audit.v1",
       exact_original_derived_EOS_zip_SHA256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
       source_contract_SHA256=hashlib.sha256(raw).hexdigest(),
       modeled_case="SPG_M3_CSS_0p7Msun_native_TOV",
       source_sound_speed_knots=direct.source_break_count,
       strong_restart_knots=direct.source_restart_count,
       source_restart_relative_jump=.01,
       relative_denominator_floors=dict(
           r=.05,Phi=1.,lambda_metric=1.,dPhi_dr_km=1e-5,cs2=1e-3),
       physically_sharp_interfaces_unchanged=True,
       phases=[])
    floors=np.array([.05,1.,1.,1e-5,1e-3])
    for k,ph in enumerate(pchip._phases):
        s=ph["s"]; mids=.5*(s[:-1]+s[1:])
        original=[]
        for sx in mids:
            h=float(bg.h_start-sx*sx)
            st=bg.state_at_enthalpy(h,
                bg.eos.side_for_segment(h,ph["mid"]))
            original.append([st["r_km"],st["Phi"],
                            st["lambda_metric"],st["dPhi_dr_km"],
                            st["cs2_over_c2"]])
        original=np.asarray(original)
        interpolated=ph["pchip"](mids)
        rel=np.abs(original-interpolated)/np.maximum(np.abs(original),floors)
        cs=rel[:,4]
        worst=np.argsort(cs)[-6:][::-1]
        row=dict(phase_index=k,
           input_grid_n=len(s),source_knots_in_phase=len(direct._phases[k]["source_break_s"]),
           selected_restart_knots=len(direct._phases[k]["source_restart_s"]),
           max_abs_cs2_normalized_error=float(np.max(cs)),
           median_abs_cs2_normalized_error=float(np.median(cs)),
           cs2_error_99th_percentile=float(np.quantile(cs,.99)),
           max_all_channel_normalized_error=float(np.max(rel)),
           source_left_and_right_of_knot=True,
           worst_cs2_pchip_mismatch=[dict(
               h=float(bg.h_start-mids[i]**2),
               exact_cs2=float(original[i,4]),
               interpolated_cs2=float(interpolated[i,4]),
               relative_error=float(cs[i]),
               s_spacing=float(s[i+1]-s[i]))
               for i in worst])
        result["phases"].append(row)
        print("KNOT_ERROR",k,"grid",len(s),"strong",
          row["selected_restart_knots"],"max",row["max_abs_cs2_normalized_error"],flush=True)
    result["status"]="SOURCE_KNOTS_AND_ORIGINAL_ONE_SIDED_SOUND_SPEED_AUDITED"
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("SOUND_SPEED_KNOT_REPORT",OUT,flush=True)
if __name__=="__main__":
    main()
