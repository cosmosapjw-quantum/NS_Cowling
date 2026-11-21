#!/usr/bin/env python3
"""C13 source-locked phase geometry PCHIP accuracy and runtime study.

Python is a source frontend and comparison oracle only. Actual production
PCHIP derivatives, phase geometry interpolation, source EOS c_s², Maxwell
junction, Brent eigenroot and DVODE are all Fortran for native_source_cs2.
"""
from __future__ import annotations
import hashlib,io,json,statistics,sys,time,zipfile
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos import stellar_cowling_native as native
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (
  NativeSharpCowlingDVODE,solve_sharp_background_native)
from compose_eos.stellar_cowling import CowlingSettings
import ctypes as ct
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SOURCE_SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
OUT=ROOT/"data/cowling_native_pchip_c13_20261121/REAL_SOURCE_GEOMETRY_CONVERGENCE.json"
CASES=(
  ("M2","SPG_M2","maxwell",.30385128926615257,"rapid",(5110.,5170.),5143.719608336112),
  ("M3","SPG_M3","css",.192510767556007,"slow",(3710.,3760.),3736.210698739710),
)
SAMPLES=(1025,2049,4097,8193)

def run():
    raw=SOURCE.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:raise RuntimeError("original EOS byte drift")
    record=dict(schema="ns_cowling.c13.original_source_phase_pchip_convergence.v1",
       baseline="a709208588367032cf34b04f502238f773625041",
       source_sha256=SOURCE_SHA,
       studied_samples=SAMPLES,
       geometry_slope_algorithm="Fortran Fritsch-Butland PCHIP",
       slope_comparison="SciPy PchipInterpolator, reference only",
       mode_algorithm="Fortran source-exact c_s² + Fortran Brent and Maxwell + DVODE",
       cases={})
    with zipfile.ZipFile(io.BytesIO(raw)) as outer:
      with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as archive:
       for label,model,join,nc,conv,br,contract_reference in CASES:
        rawrep=archive.read(f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json")
        eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(rawrep)))
        bg=solve_sharp_background_native(eos,nc,surface_h=3e-8,rtol=3e-10,
              atol=3e-11,max_step_h=.004,native_samples=2049)
        opts=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,
              max_step_h=.002,root_xtol_hz=1e-5)
        case=dict(source_contract_sha256=hashlib.sha256(rawrep).hexdigest(),
             model=model,join=join,conversion=conv,central_density_fm3=nc,
             mass_msun=bg.mass_msun,radius_km=bg.radius_km,
             original_c12_4097_frequency_hz=contract_reference,runs={})
        for n in SAMPLES:
          t=time.perf_counter()
          solver=NativeSharpCowlingDVODE(bg,opts,samples=n,
               coefficient_backend="native_source_cs2",
               source_restart_relative_jump=.01)
          construction=time.perf_counter()-t
          max_rel_slope=0.
          max_abs_slope=0.
          medians={}
          for phase_i,phase in enumerate(solver._phases):
            sci=PchipInterpolator(phase["s"],phase["vals"].T,axis=0,
                                 extrapolate=False)
            original=np.asfortranarray(sci(phase["s"],nu=1).T)
            diff=np.abs(original-phase["slopes"])
            max_rel_slope=max(max_rel_slope,float(np.max(diff / np.maximum(np.abs(original),1e-8))))
            max_abs_slope=max(max_abs_slope,float(np.max(diff)))
            assert np.allclose(original,phase["slopes"],rtol=2e-10,atol=5e-9)
            # Five warm samples each; the source PCHIP is an oracle only.
            old=[];fortran=[]
            v=np.asfortranarray(phase["vals"])
            slopes=np.empty_like(v,order="F")
            dll=native._native_library()
            status=ct.c_int(-1)
            for j in range(7):
                t1=time.perf_counter()
                _=PchipInterpolator(phase["s"],v.T,axis=0,extrapolate=False)
                old.append(time.perf_counter()-t1)
                t1=time.perf_counter()
                dll.cowling_native_pchip_slopes(len(phase["s"]),
                    phase["s"],v,slopes,ct.byref(status))
                fortran.append(time.perf_counter()-t1)
                assert status.value==0
            medians[str(phase_i)]=dict(python_scipy=statistics.median(old),
                 native_fortran=statistics.median(fortran),
                 native_over_scipy=statistics.median(old)/statistics.median(fortran),
                 knot_count=len(phase["s"]))
          row=dict(construction_seconds=construction,
             max_abs_slope_delta=max_abs_slope,
             max_scaled_slope_delta=max_rel_slope,
             source_restart_count=solver.source_restart_count,
             profile_geometry_backend=solver.geometry_interpolation_backend,
             per_phase_pchip_build_medians=medians)
          try:
            t=time.perf_counter()
            m=solver.solve_mode(br,conversion=conv,profile_points=300,
                                root_backend="fortran_source_exact")
            row.update(status="PASS_MODE",frequency_hz=m["frequency_hz"],
                delta_previous_4097_hz=m["frequency_hz"]-contract_reference,
                matched_nodes=[m["W_nodes_per_phase"],m["V_nodes_per_phase"]],
                join_relative=m["relative_join_mismatch"],
                det=m["matching_residual"],
                python_rhs_eos_calls=m["source_callback_rhs_calls"],
                root_cut_count=m["root_fortran_restart_count"],
                mode_seconds=time.perf_counter()-t)
            assert m["root_fortran_restart_count"]==solver.source_restart_count
            assert m["source_callback_rhs_calls"]==0
            assert row["matched_nodes"]==[1,1]
            assert m["relative_join_mismatch"]<1e-5
            if n==4097:assert abs(m["frequency_hz"]-contract_reference)<.005
          except (RuntimeError,ValueError) as exc:
            row.update(status="ROOT_OR_MODE_REJECTED",reason=str(exc))
            if n==4097:raise
          case["runs"][str(n)]=row
          record["cases"][label]=case
          OUT.parent.mkdir(exist_ok=True,parents=True)
          OUT.write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
          print("C13_CONVERGENCE",label,n,row["status"],
                "hz",row.get("frequency_hz"),
                "slope_max_scaled",max_rel_slope,
                "construct_s",round(construction,4),flush=True)
    record["status"]="SCOPED_CONVERGENCE_COMPLETED_WITH_REJECTIONS_PRESERVED"
    OUT.write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":run()
