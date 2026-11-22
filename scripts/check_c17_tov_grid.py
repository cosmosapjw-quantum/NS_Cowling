#!/usr/bin/env python3
"""C17 independent TOV source-grid sensitivity at four extended branch stars.

The prior C17 original 2049-source-sample roots are baseline. Rebuild
original EOS source/TOV background at 4097 and 8193 samples, holding
cold closure, density, Maxwell conversion, Cowling controls and root
bracket fixed. All differential equations, original EOS c_s2 and root
polishing use Fortran. Python only marshals original source and writes
a status-preserving finite-control comparison.

Differences are NOT global frequency error bounds or Richardson rates.
"""
from __future__ import annotations
import hashlib,io,json,sys,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (
    solve_sharp_background_native,NativeSharpCowlingDVODE)
from compose_eos.stellar_cowling import CowlingSettings
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SOURCE_SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
D=ROOT/"data/cowling_c17_extended_branch_20261122"
SAMPLE_CONTROLS=(4097,8193)
SCOPES=(("M2",.6),("M2",1.),("M3",.9),("M3",1.))

def run():
  source=SOURCE.read_bytes()
  if hashlib.sha256(source).hexdigest()!=SOURCE_SHA:
      raise RuntimeError("source EOS lock drift")
  cases={}
  out=dict(schema="ns_cowling.c17.original_source_TOV_grid_crosscheck.v1",
      original_source_sha256=SOURCE_SHA,
      c17_sweeps="REAL_P2_INITIAL original 2049 TOV grid anchors",
      new_background_requested_source_samples=SAMPLE_CONTROLS,
      cowling_requested_samples=4097,
      numerical_differences_not_absolute_error_bounds=True,
      cases=cases)
  with zipfile.ZipFile(io.BytesIO(source)) as full:
    with zipfile.ZipFile(io.BytesIO(full.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as nested:
      for model,factor in SCOPES:
        config={"M2":("SPG_M2","maxwell",.30385128926615257,"rapid"),
                "M3":("SPG_M3","css",.192510767556007,"slow")}[model]
        name,kind,nc,conversion=config
        raw=nested.read(f"FRESH_5x2_WHOLE_STARS/{name}/{kind}/REPRESENTATION.json")
        eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
        canonical=json.loads((D/"REAL_P2_INITIAL"/f"{model}.json").read_text())
        row=next(x for x in canonical["densities"] if x["density_factor"]==factor)
        optics=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,max_step_h=.002,
                               root_xtol_hz=1e-5)
        baseline={}
        for label,required_nodes in (("f",0),("p1",1)):
            selected=row["local_branches"][label]
            if selected["frequency_hz"] is None:
                continue
            baseline[label]=dict(legacy_frequency_hz=selected["frequency_hz"],
                 W_nodes=required_nodes,V_nodes=required_nodes)
        results=[]
        for samples in SAMPLE_CONTROLS:
          rec=dict(TOV_source_samples_requested=samples)
          try:
            bg=solve_sharp_background_native(eos,nc*factor,
              surface_h=3e-8,rtol=3e-10,atol=3e-11,
              max_step_h=.004,native_samples=samples)
            rec.update(star_status="NATIVE_TOV_SUCCESS",
                       mass_msun=bg.mass_msun,radius_km=bg.radius_km)
          except (ValueError,RuntimeError) as exc:
            rec.update(star_status="TOV_REJECTED",reason=str(exc))
            results.append(rec)
            continue
          try:
            solver=NativeSharpCowlingDVODE(bg,optics,samples=4097,
                coefficient_backend="native_source_cs2",
                source_restart_relative_jump=.01)
            rec["geometry_status"]="FORTRAN_SOURCE_GEOMETRY_SUCCESS"
          except (ValueError,RuntimeError) as exc:
            rec.update(geometry_status="GEOMETRY_REJECTED",reason=str(exc))
            results.append(rec)
            continue
          modes={}
          for label,ref in baseline.items():
            f=ref["legacy_frequency_hz"]
            mode_rec=dict(reference_2049_hz=f,expected_W_V_nodes=ref["W_nodes"])
            try:
              m=solver.solve_mode((f-40.,f+40.),conversion=conversion,
                   profile_points=360,root_backend="fortran_source_exact")
              if m["W_nodes_per_phase"]!=ref["W_nodes"] or m["V_nodes_per_phase"]!=ref["V_nodes"]:
                mode_rec.update(status="NODE_IDENTITY_REJECTED",
                     observed_nodes=[m["W_nodes_per_phase"],m["V_nodes_per_phase"]])
              elif m["relative_join_mismatch"]>=1e-5:
                mode_rec.update(status="JOIN_REJECTED",
                     original_join=m["relative_join_mismatch"])
              else:
                mode_rec.update(status="SOURCE_EXACT_MODE_ADMITTED",
                   frequency_hz=m["frequency_hz"],
                   difference_from_2049_hz=m["frequency_hz"]-f,
                   relative_join_mismatch=m["relative_join_mismatch"],
                   original_source_cs2_python_callback_calls=m["source_callback_rhs_calls"],
                   native_source_knot_restarts=m["root_fortran_restart_count"])
            except (ValueError,RuntimeError) as exc:
              mode_rec.update(status="ROOT_OR_PROFILE_REJECTED",reason=str(exc))
            modes[label]=mode_rec
          rec["modes"]=modes
          results.append(rec)
          print("C17_TOV_GRID",model,factor,samples,
                [(k,v["status"],v.get("difference_from_2049_hz")) for k,v in modes.items()],
                flush=True)
        key=f"{model}@{factor:.3f}"
        cases[key]=dict(source_contract_sha256=hashlib.sha256(raw).hexdigest(),
          original_center_fm3=nc*factor,
          original_2049_reference=baseline,
          methods="Fortran native TOV Adams and exact-source Cowling Adams",
          results=results)
  dest=D/"TOV_GRID_CROSSCHECK.json"
  dest.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n")
  print("C17_TOV_CROSSCHECK_WRITTEN",dest,flush=True)

if __name__=="__main__":
  run()
