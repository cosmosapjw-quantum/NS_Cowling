#!/usr/bin/env python3
"""C14 physical mode admissibility at selected ODE tolerance extremes.

Python only marshals source EOS/stellar inputs, Fortran performs TOV,
source-exact c_s², native geometry PCHIP, DVODE and Brent; the
frontend checks final eigenfunction node/junction classification.
"""
from __future__ import annotations
import hashlib,io,json,sys,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.stellar_cowling import CowlingSettings

SRC=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
OUT=ROOT/"data/cowling_c14_source_error_budget_20261121/MODE_ADMISSIBILITY.json"
CASES=(("M2","SPG_M2","maxwell",.30385128926615257,"rapid",(5110.,5170.)),
       ("M3","SPG_M3","css",.192510767556007,"slow",(3710.,3760.)))
RTOLS=(3e-8,1e-9,1e-10)

def run():
    raw=SRC.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SHA:raise RuntimeError("original EOS ZIP SHA mismatch")
    cases={}
    with zipfile.ZipFile(io.BytesIO(raw)) as a:
      with zipfile.ZipFile(io.BytesIO(a.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as z:
       for short,model,join,nc,conv,br in CASES:
        rep=z.read(f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json")
        eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(rep)))
        bg=solve_sharp_background_native(eos,nc,surface_h=3e-8,
               rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
        rows=[]
        for mf in ("adams","bdf"):
         for rtol in RTOLS:
          st=CowlingSettings(surface_h=3e-8,rtol=rtol,
             atol=1e-12,max_step_h=.002,root_xtol_hz=1e-5)
          obj=NativeSharpCowlingDVODE(bg,st,samples=4097,
             coefficient_backend="native_source_cs2",
             source_restart_relative_jump=.01,method=mf)
          result=dict(method=mf,rtol=rtol,method_fortran_root="fortran_source_exact")
          try:
            mode=obj.solve_mode(br,conversion=conv,profile_points=300,
                             root_backend="fortran_source_exact")
            result.update(status="PASS_MODE",
              frequency_hz=mode["frequency_hz"],
              W_nodes=mode["W_nodes_per_phase"],
              V_nodes=mode["V_nodes_per_phase"],
              join_relative=mode["relative_join_mismatch"],
              det=mode["matching_residual"],
              source_knots=mode["root_fortran_restart_count"],
              source_callback_count=mode["source_callback_rhs_calls"])
            if (result["W_nodes"],result["V_nodes"])!=(1,1):
              raise RuntimeError("non-p1 node identity")
            if result["join_relative"]>=1e-5:
              raise RuntimeError("join outside original fixed threshold")
          except Exception as e:
            result["status"]="PHYSICAL_OR_NUMERICAL_REJECTED"
            result["reason"]=str(e)
          rows.append(result)
          print("C14_MODE_ADMISSIBILITY",short,mf,rtol,result["status"],
                result.get("frequency_hz"),"join",result.get("join_relative"),flush=True)
        cases[short]=dict(contract_sha256=hashlib.sha256(rep).hexdigest(),rows=rows)
    OUT.write_text(json.dumps(dict(schema="ns_cowling.c14.mode_admissibility.v1",
       source_zip_sha256=SHA,closure="Gamma1=Gamma_eq",
       controls="3 rtols x 2 DVODE methods x 2 source-locked stars",
       cases=cases),sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("C14_PHYSICAL_PROFILE_CHECK_WRITTEN",OUT)
if __name__=="__main__":run()
