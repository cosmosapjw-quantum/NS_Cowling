#!/usr/bin/env python3
"""C16 original-source pilot: Fortran determinant scan + Fortran eigenmodes.

A dense fixed frequency grid is a *candidate screen*, not a physical
eigenfrequency proof; every candidate needs native Brent and profile
node/interface admission. Preserve all failed status codes.
"""
import hashlib,io,json,sys,time,zipfile
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
OUT=ROOT/"data/cowling_c16_source_continuation_20261121/PILOT_SCAN.json"
CASES=[("M2","SPG_M2","maxwell",.30385128926615257,"rapid"),("M3","SPG_M3","css",.192510767556007,"slow")]

def run():
 data=SRC.read_bytes()
 assert hashlib.sha256(data).hexdigest()==SHA
 items={}
 with zipfile.ZipFile(io.BytesIO(data)) as outer:
  with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inner:
   for label,model,kind,nc,conv in CASES:
    raw=inner.read(f"FRESH_5x2_WHOLE_STARS/{model}/{kind}/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
    bg=solve_sharp_background_native(eos,nc,surface_h=3e-8,
             rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
    opts=CowlingSettings(surface_h=3e-8,rtol=1e-9,
             atol=1e-12,max_step_h=.002,root_xtol_hz=1e-5)
    solver=NativeSharpCowlingDVODE(bg,opts,samples=4097,
             coefficient_backend="native_source_cs2",source_restart_relative_jump=.01)
    freq=np.ascontiguousarray(np.r_[np.linspace(700,3400,55),np.linspace(3450,6850,69)])
    t=time.perf_counter()
    brackets,samples=solver.scan_source_brackets_native(freq,conversion=conv)
    scan_s=time.perf_counter()-t
    print("SCAN",label,"samples",len(freq),"brackets",brackets,"scan_s",scan_s,flush=True)
    selected=[]
    for lo,hi in brackets[:15]:
     rr=dict(bracket_hz=[lo,hi])
     try:
      mode=solver.solve_mode((lo,hi),conversion=conv,
                profile_points=360,root_backend="fortran_source_exact")
      rr.update(status="PASS_ROOT_AND_JOIN",
           frequency_hz=mode["frequency_hz"],
           W_nodes=mode["W_nodes_per_phase"],V_nodes=mode["V_nodes_per_phase"],
           determinant=mode["matching_residual"],
           join=mode["relative_join_mismatch"],
           native_rhs_calls=mode["root_native_rhs_calls"],
           node_count_backend=mode["node_count_backend"])
     except (ValueError,RuntimeError) as e:
      rr.update(status="REJECTED",reason=str(e))
     selected.append(rr)
     print("CANDIDATE",label,lo,hi,rr["status"],rr.get("frequency_hz"),
         rr.get("W_nodes"),rr.get("V_nodes"),rr.get("reason"),flush=True)
    items[label]=dict(source_sha=hashlib.sha256(raw).hexdigest(),
        central_density=nc,conversion=conv,
        M_solar=bg.mass_msun,R_km=bg.radius_km,
        scan_points=len(freq),scan_s=scan_s,
        failed_scan_points=sum(p["status"]!=0 for p in samples),
        scanned_values=samples,candidate_brackets=brackets,
        roots=selected)
    OUT.write_text(json.dumps(dict(schema="c16.pilot.original_source.v1",archive_sha=SHA,models=items),
                         indent=2,sort_keys=True,allow_nan=False)+"\n")
if __name__=="__main__":run()
