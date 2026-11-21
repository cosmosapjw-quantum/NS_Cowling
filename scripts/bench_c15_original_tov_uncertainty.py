#!/usr/bin/env python3
"""C15 original-source M2/M3 coupled TOV-to-Cowling numerical error-budget.

Python is strictly the frontend: load source ZIP, sample EOS once per
background grid, marshal Fortran phasewise TOV controls, serialize outputs.
TOV PCHIP slopes + integration, Cowling PCHIP + exact cs2,
strong-source-knot scheduling, DVODE and Brent all execute in Fortran.
The last Python mode profile node classification checks physical admission.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (
  solve_sharp_background_native_sweep, NativeSharpCowlingDVODE)
from compose_eos.stellar_cowling import CowlingSettings

SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
CASES={
 "M2":("SPG_M2","maxwell",.30385128926615257,"rapid",(5110.,5170.),1.3999999857354968,14.420419623651872),
 "M3":("SPG_M3","css",.192510767556007,"slow",(3710.,3760.),.7000000086373105,14.59794179069283)
}
BASE=dict(rtol=3e-10,atol=3e-11,max_step_h=.004)
TOV_RTOLS=(3e-8,1e-9,3e-10,1e-10)
TOV_MAXSTEPS=(.008,.004,.002,.001)
TOV_SAMPLES=(1025,2049,4097,8193)
METHODS=("adams","bdf")

def controls_for(grid):
    ret=[]
    keys=set()
    def add(axis,method,rtol,maxh):
        key=(method,rtol,maxh)
        if key not in keys:
            ret.append(dict(axis=axis,rtol=rtol,atol=BASE["atol"],
                max_step_h=maxh,method=method,samples=grid))
            keys.add(key)
    if grid==2049:
        for method in METHODS:
            for value in TOV_RTOLS:
                add("TOV_RTOL",method,value,BASE["max_step_h"])
            for value in TOV_MAXSTEPS:
                add("TOV_MAX_STEP_H",method,BASE["rtol"],value)
    else:
        for method in METHODS:
            add("TOV_THERMODYNAMIC_SOURCE_GRID",method,BASE["rtol"],
                BASE["max_step_h"])
    return ret

def main(case,out,quiet=False):
    raw=SOURCE.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SHA:raise RuntimeError("locked source EOS SHA mismatch")
    model,kind,nc,conversion,br,ref_m,ref_r=CASES[case]
    with zipfile.ZipFile(io.BytesIO(raw)) as a:
      with zipfile.ZipFile(io.BytesIO(a.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as b:
        blob=b.read(f"FRESH_5x2_WHOLE_STARS/{model}/{kind}/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(blob)))
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,
                  atol=1e-12,max_step_h=.002,root_xtol_hz=1e-5)
    records=[]
    # The same Cowling numerical controls and geometry resolution are
    # fixed across every TOV background perturbation.
    for samples in TOV_SAMPLES:
        controls=controls_for(samples)
        start=time.perf_counter()
        profiles=solve_sharp_background_native_sweep(eos,nc,controls,
                        surface_h=3e-8,native_samples=samples)
        batch_seconds=time.perf_counter()-start
        if len(profiles)!=len(controls):raise RuntimeError("TOV sweep row mismatch")
        for ctrl,profile in zip(controls,profiles):
            row=dict(axis=ctrl["axis"],TOV_method=ctrl["method"],
              TOV_rtol=ctrl["rtol"],TOV_atol=ctrl["atol"],
              TOV_max_step_h=ctrl["max_step_h"],TOV_source_samples=samples,
              TOV_effective_max_step_s=profile["effective_max_step_s"],
              TOV_status=profile["status"],
              TOV_fortran_rhs_evaluations=profile["tov_native_nfev"],
              background_phase_setup_group=samples,
              fixed_Cowling_rtol=settings.rtol,
              fixed_Cowling_max_step_h=settings.max_step_h,
              fixed_Cowling_samples=4097)
            if profile["status"]:
                row["mode_status"]="SKIPPED_TOV_FAILURE"
                records.append(row)
                continue
            star=profile["star"]
            row.update(TOV_mass_msun=float(star.mass_msun),
              TOV_radius_km=float(star.radius_km),
              TOV_interface_radius_km=float(star.interfaces[0]["r_km"]),
              TOV_interface_mass_geom_km=float(star.interfaces[0]["m_geom_km"]),
              TOV_compactness=2*star.mass_geom_km/star.radius_km,
              TOV_backend=star.config["backend"],
              TOV_thermodynamics_slopes_backend=star.config["tov_eos_slopes_backend"])
            try:
                mode_solver=NativeSharpCowlingDVODE(star,settings,samples=4097,
                   coefficient_backend="native_source_cs2",
                   source_restart_relative_jump=.01)
                mode=mode_solver.solve_mode(br,conversion=conversion,
                   profile_points=300,root_backend="fortran_source_exact")
                row.update(mode_status="PASS_P1",frequency_hz=float(mode["frequency_hz"]),
                    matching_residual=float(mode["matching_residual"]),
                    relative_join_mismatch=float(mode["relative_join_mismatch"]),
                    W_nodes=int(mode["W_nodes_per_phase"]),
                    V_nodes=int(mode["V_nodes_per_phase"]),
                    source_restart_count=int(mode["root_fortran_restart_count"]),
                    source_callback_rhs_calls=int(mode["source_callback_rhs_calls"]),
                    root_native_rhs_evaluations=int(mode["root_native_rhs_calls"]))
                if (row["W_nodes"],row["V_nodes"])!=(1,1) or row["relative_join_mismatch"]>=1e-5:
                    raise RuntimeError("physical p1 mode gate failed")
                if abs(row["frequency_hz"] - (5143.7197 if case=="M2" else 3736.2107))>.05:
                    raise RuntimeError("p1 branch anchor exceeded 0.05 Hz source-locked diagnostic band")
            except (RuntimeError,ValueError) as exc:
                row["mode_status"]="NUMERICALLY_OR_PHYSICALLY_REJECTED"
                row["reason"]=str(exc)
                for field in ("frequency_hz","matching_residual","relative_join_mismatch"):
                    row.pop(field,None)
            records.append(row)
            if not quiet:
                print("C15_TOV",case,ctrl["axis"],samples,ctrl["method"],
                  ctrl["rtol"],ctrl["max_step_h"],row["mode_status"],
                  row.get("frequency_hz"),flush=True)
    stable=[r for r in records if r["mode_status"]=="PASS_P1"]
    assert len(records)==20,(case,len(records))
    baseline=[r for r in stable if r["TOV_method"]=="adams" and
        r["TOV_rtol"]==BASE["rtol"] and r["TOV_max_step_h"]==BASE["max_step_h"]
        and r["TOV_source_samples"]==2049]
    if len(baseline)!=1:raise RuntimeError("missing or duplicate original TOV baseline")
    o=baseline[0]
    assert abs(o["TOV_mass_msun"]-ref_m)<.0001
    assert abs(o["TOV_radius_km"]-ref_r)<.01
    result=dict(schema="ns_cowling.c15.two_phase_native_TOV_uncertainty.v1",
        source_archive_sha256=SHA,source_contract_sha256=hashlib.sha256(blob).hexdigest(),
        baseline_commit="091c2f6f22d032379324aa900c910c6278a68437",
        model=model,kind=kind,star_center_n_fm3=nc,
        mode="p1",conversion=conversion,
        TOV_original_reference_M_msun=ref_m,TOV_original_reference_R_km=ref_r,
        baseline_frequency_hz=o["frequency_hz"],
        numerical_controls={"TOV_rtol":TOV_RTOLS,"TOV_max_step_h":TOV_MAXSTEPS,
            "TOV_source_samples":TOV_SAMPLES,"TOV_methods":METHODS,
            "fixed_Cowling_rtol":settings.rtol,
            "fixed_Cowling_samples":4097},
        samples_per_original_source_grid="source breakpoints are inserted; values are requested lower limits",
        all_rows=records,
        numerical_phase="TOV and Cowling equations computed in Fortran; Python source-input/output frontend",
        status="SCOPED_SOURCE_LOCKED_TOV_CONTROLS_MEASURED",
        scope="fixed source EOS, central density, Maxwell semantics, Cowling controls")
    out.parent.mkdir(exist_ok=True,parents=True)
    out.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("C15_FINAL",case,"CONTROLS",len(records),"MODES_ADMITTED",len(stable),
       "BASELINE_FREQ",o["frequency_hz"],"MODE_REJECTIONS",len(records)-len(stable),
       "INPUT_PREP_WALL_NOT_IN_SCIENCE",flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--case",choices=CASES,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--quiet",action="store_true")
    args=p.parse_args()
    main(args.case,args.output,args.quiet)
