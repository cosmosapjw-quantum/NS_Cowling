#!/usr/bin/env python3
"""C16 matched-physics warm native vs Python-loop source determinant scan.

Both backends integrate the ORIGINAL analytic EOS source c_s^2 and
the same Fortran DVODE with the same physical Maxwell junction.
Only the *frequency scan orchestration* differs:
- native: one Fortran ABI, source-knot schedule prepared once;
- python loop: one Python call per Fortran shot (reference only).
This is not a new EOS model nor full-pipeline speedup.
"""
from __future__ import annotations
import hashlib,io,json,statistics,sys,time,zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.stellar_cowling import CowlingSettings
SRC=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
OUT=ROOT/"data/cowling_c16_source_continuation_20261121/SCANNER_PERFORMANCE.json"
CASES={
 "M2":("SPG_M2","maxwell",.30385128926615257,"rapid",(1710.,1980.),(4920.,5320.)),
 "M3":("SPG_M3","css",.192510767556007,"slow",(1430.,1710.),(3490.,3980.))
}
def run():
    raw=SRC.read_bytes()
    assert hashlib.sha256(raw).hexdigest()==SHA
    payload=dict(schema="ns_cowling.c16.original_source_scan_performance.v1",
       source_archive_sha256=SHA,scope="32 warm determinant samples, same original EOS physical source and Fortran DVODE",
       repeats=5,setup_included=False,cases={})
    with zipfile.ZipFile(io.BytesIO(raw)) as outer:
      with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as z:
       for short,(model,kind,nc,conv,fb,pb) in CASES.items():
        rep=z.read(f"FRESH_5x2_WHOLE_STARS/{model}/{kind}/REPRESENTATION.json")
        eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(rep)))
        bg=solve_sharp_background_native(eos,nc,surface_h=3e-8,
                  rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
        opt=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,
                  max_step_h=.002,root_xtol_hz=1e-5)
        solver=NativeSharpCowlingDVODE(bg,opt,samples=4097,
                  coefficient_backend="native_source_cs2",
                  source_restart_relative_jump=.01)
        frequencies=np.ascontiguousarray(np.r_[np.linspace(*fb,15),np.linspace(*pb,17)],dtype=float)
        scores={"native":[],"python_loop":[]}
        max_abs_residual_error=0.
        for it in range(5):
            for name in (("native","python_loop") if it%2==0 else ("python_loop","native")):
                t=time.perf_counter()
                if name=="native":
                    brackets,rows=solver.scan_source_brackets_native(frequencies,conversion=conv)
                    vals=np.asarray([r["determinant"] for r in rows],float)
                    assert all(r["status"]==0 for r in rows)
                else:
                    brackets,rows=solver.scan_brackets(frequencies,conversion=conv)
                    vals=np.asarray([r["residual"] for r in rows],float)
                elapsed=time.perf_counter()-t
                scores[name].append(elapsed)
                if name=="python_loop":
                    old=vals
                if name=="native":
                    new=vals
            max_abs_residual_error=max(max_abs_residual_error,float(np.max(np.abs(old-new))))
        med={k:statistics.median(v) for k,v in scores.items()}
        payload["cases"][short]=dict(original_contract_sha256=hashlib.sha256(rep).hexdigest(),
            original_model=model,method="Adams",
            num_frequencies=len(frequencies),
            serial_times_seconds=scores,
            median_seconds=med,
            speedup_same_physics=med["python_loop"]/med["native"],
            maximum_absolute_normalized_determinant_difference=max_abs_residual_error,
            callback_physics_changed=False)
        print("C16_SCANNER_BENCH",short,"PYTHON_S",med["python_loop"],
              "FORTRAN_S",med["native"],
              "SPEEDUP",payload["cases"][short]["speedup_same_physics"],
              "MAX_NORMALIZED_DET_DELTA",max_abs_residual_error,flush=True)
    OUT.write_text(json.dumps(payload,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("SCANNER_BENCH_RECORD",OUT,flush=True)
if __name__=="__main__":
    run()
