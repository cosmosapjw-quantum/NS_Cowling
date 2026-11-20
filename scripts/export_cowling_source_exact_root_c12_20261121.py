#!/usr/bin/env python3
"""Export source-locked EOS laws and PCHIP-only geometry to independent Fortran.

This frontend produces no eigenfrequency and performs no mode shooting.
The Fortran executable reads the numeric contract and solves all roots
using its own original-source c_s², strong-knot scheduler and DVODE.
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
  NativeSharpCowlingDVODE,solve_sharp_background_native)
from compose_eos.stellar_cowling import CowlingSettings
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
OUT=ROOT/"data/cowling_fortran_source_exact_root_c12_20261121"
CASES=(
 ("M2","SPG_M2","maxwell",.30385128926615257,"rapid",(5110.,5170.)),
 ("M3","SPG_M3","css",.192510767556007,"slow",(3710.,3760.)),
)

def run():
    archive=SOURCE.read_bytes()
    if hashlib.sha256(archive).hexdigest()!=SHA:raise RuntimeError("source lock drift")
    recorded=json.loads((OUT/"REAL_SOURCE_BENCHMARK.json").read_text())["cases"]
    with zipfile.ZipFile(io.BytesIO(archive)) as outer:
      with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inside:
        for key,model,join,nc,conversion,br in CASES:
          blob=inside.read(f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json")
          sha=hashlib.sha256(blob).hexdigest()
          eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(blob)))
          bg=solve_sharp_background_native(eos,nc,surface_h=3e-8,
                rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
          st=CowlingSettings(surface_h=3e-8,rtol=1e-9,
                atol=1e-12,max_step_h=.002,root_xtol_hz=1e-5)
          solver=NativeSharpCowlingDVODE(bg,st,samples=4097,
                 coefficient_backend="native_source_cs2",
                 source_restart_relative_jump=.01)
          outer_eos,upper_eos,surface,bridge,jump=solver._native_source_eos
          high,low=solver._phases
          ht=bg.interfaces[0]["h"]
          physical_low=bg.state_at_enthalpy(ht,"low")
          physical_high=bg.state_at_enthalpy(ht,"high")
          geom=np.array([
              physical_low["r_km"],physical_low["Phi"],
              physical_low["lambda_metric"],physical_low["dPhi_dr_km"],
              physical_low["p"],physical_low["epsilon"],
              physical_high["epsilon"]],float)
          bench_key=next(k for k in recorded if k.startswith(key+"_"))
          ref=recorded[bench_key]
          assert sha==ref["source_contract_sha256"]
          expectedcuts=ref["root_restart_cuts_in_fortran"]
          assert expectedcuts==solver.source_restart_count
          inp=OUT/f"ORIGINAL_FORTRAN_ROOT_{key}.dat"
          with inp.open("w") as f:
            f.write(f"{outer_eos.shape[1]} {upper_eos.shape[1]} {jump}\n")
            f.write(" ".join(f"{v:.17e}" for v in surface)+"\n")
            f.write(" ".join(f"{v:.17e}" for v in bridge)+"\n")
            np.savetxt(f,outer_eos.T,fmt="%.17e")
            np.savetxt(f,upper_eos.T,fmt="%.17e")
            f.write(f"{len(high['s'])} {len(low['s'])}\n")
            max_ds=min(solver._smax/200,st.max_step_h/(2*solver._smax))
            config=[bg.h_start,.5*high["s"][-1],
                    st.rtol,st.atol,max_ds,st.root_xtol_hz,.01]
            f.write(" ".join(f"{v:.17e}" for v in config)+"\n")
            f.write(" ".join(f"{v:.17e}" for v in geom)+"\n")
            for phase in (high,low):
              matrix=np.column_stack((phase["s"],phase["vals"].T,phase["slopes"].T))
              np.savetxt(f,matrix,fmt="%.17e")
            f.write(f"{1 if conversion=='slow' else 2} "+
              f"{br[0]:.17e} {br[1]:.17e} "+
              f"{ref['original_py_brent_freq_hz']:.17e} {expectedcuts}\n")
          print("EXPORTED_FORTRAN_ROOT",key,"model",model,"sha_source",sha,
                "sha_input",hashlib.sha256(inp.read_bytes()).hexdigest(),
                "phase_samples",len(high["s"]),len(low["s"]),
                "restart_cuts",expectedcuts,flush=True)

if __name__=="__main__":run()
