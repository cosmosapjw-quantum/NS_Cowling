#!/usr/bin/env python3
"""Bounded native sharp-Maxwell/CSS convergence and independent Python oracle.

Source-locked snapshot from earlier exact code. Do not assume library first
build time or EOS source uncertainty is captured by 2049/4097/8193 knots.
"""
from __future__ import annotations
import hashlib,json,sys,time,traceback
from pathlib import Path
from math import isfinite
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import (SharpContractEnthalpyEOS,
                                      solve_sharp_background,
                                      solve_mode as solve_python_sharp)
from compose_eos.stellar_sharp_native import (solve_sharp_background_native,
                                              NativeSharpCowlingDVODE)
from compose_eos.stellar_cowling import CowlingSettings

BASE=ROOT/"data/cowling_maxwell_native_20261119"
DEST=BASE/"validated_convergence"
SAME=json.loads((BASE/"three_layer_original_fsu2/RESULT.json").read_text())
OFFICIAL=json.loads((BASE/"official_unified_css_sharp/RESULT.json").read_text())
BBP=json.loads((BASE/"bbp_bps/RESULT.json").read_text())
BBP_FSU2=next(x for x in BBP["combinations"] if x["branch"]=="pressure_primary_first_law_BBP"
              and x["core"]=="GPPVA_FSU2")

CASES=[
 ("same_FSU2_maxwell", BASE/"three_layer_original_fsu2/maxwell/REPRESENTATION.json",
   SAME["methods"]["maxwell"]["admitted"]["stellar"]),
 ("same_FSU2_css", BASE/"three_layer_original_fsu2/css/REPRESENTATION.json",
   SAME["methods"]["css"]["admitted"]["stellar"]),
 ("original_FSU2_css", BASE/"official_unified_css_sharp/FSU2/REPRESENTATION.json",
   OFFICIAL["cases"]["FSU2"]["native_fortran"]),
 ("original_DDME2_css", BASE/"official_unified_css_sharp/DDME2/REPRESENTATION.json",
   OFFICIAL["cases"]["DDME2"]["native_fortran"]),
 ("BBP_BPS_FSU2_css", BASE/"bbp_bps/whole_models/GPPVA_FSU2/css/REPRESENTATION.json",
   BBP_FSU2["css_three_layer"]["stellar"]),
]
MODE=CowlingSettings(surface_h=3e-8,rtol=3e-9,atol=3e-12,max_step_h=.001)
TOV=dict(surface_h=3e-8,rtol=3e-9,atol=3e-11,max_step_h=.004)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+"\n")

def nearest_bracket(freq,relative=.012):
    return [float(freq*(1-relative)),float(freq*(1+relative))]

def run_one(name,source,baseline):
    original_modes={
      mode:{v["label"]:float(v["frequency_hz"]) for v in baseline["modes"][mode]["modes"]}
      for mode in ("slow","rapid")}
    if any(set(original_modes[k])!={"f","p1"} for k in original_modes):
        raise RuntimeError(name+": original bounded mode list incomplete")
    whole=WholeCrustBarotrope.from_contract(json.loads(source.read_text()))
    eos=SharpContractEnthalpyEOS(whole)
    record=dict(name=name,contract=str(source.relative_to(ROOT)),
       sha256=sha(source),input_reference_frequencies_hz=original_modes,
       resolution={})
    for ns in (2049,4097,8193):
        start=time.perf_counter()
        b=solve_sharp_background_native(eos,.3,native_samples=ns,**TOV)
        builder=NativeSharpCowlingDVODE(b,MODE,samples=ns)
        xs=dict(background=dict(M_Msun=b.mass_msun,R_km=b.radius_km,
                  transition_radius_km=b.interfaces[0]["r_km"],
                  n_ode_evals=b.solution.nfev),
                modes={},wall_seconds=None)
        for conv in ("slow","rapid"):
            xs["modes"][conv]={}
            for label,f in original_modes[conv].items():
                try:
                    m=builder.solve_mode(nearest_bracket(f),conversion=conv,
                                         profile_points=180)
                    residual=float(m["relative_join_mismatch"])
                    jmax=max((z["traction_relative_residual"] for z in m["interface_junctions"]),default=0.)
                    xs["modes"][conv][label]=dict(status="CHECKED",
                        f_hz=float(m["frequency_hz"]),relative_join=residual,
                        interface_traction_residual=jmax,
                        W_nodes=m["W_nodes_per_phase"],V_nodes=m["V_nodes_per_phase"],
                        frequency_minus_4097_baseline_hz=float(m["frequency_hz"]-f))
                except Exception as e:
                    xs["modes"][conv][label]=dict(status="REJECTED",
                        type=type(e).__name__,error=str(e)[:280])
        xs["wall_seconds"]=time.perf_counter()-start
        record["resolution"][str(ns)]=xs
        print("CONVERGENCE",name,ns,"M",xs["background"]["M_Msun"],
              "R",xs["background"]["R_km"],
              "modes",[(k,l,v["status"]) for k,cs in xs["modes"].items()
                           for l,v in cs.items()],flush=True)
    # The original source check is independent from the native interpolating
    # table. It costs more, so limit to one source/control for full DOP853.
    if name=="same_FSU2_maxwell":
        start=time.perf_counter()
        py=solve_sharp_background(eos,.3,**TOV)
        ref=dict(M_Msun=py.mass_msun,R_km=py.radius_km,
                 native_DOP853_background_diff_msun=record["resolution"]["8193"]["background"]["M_Msun"]-py.mass_msun,
                 native_DOP853_background_diff_km=record["resolution"]["8193"]["background"]["R_km"]-py.radius_km,
                 modes={})
        for conv in ("slow","rapid"):
            ref["modes"][conv]={}
            for label,f in original_modes[conv].items():
                try:
                    result=solve_python_sharp(py,nearest_bracket(f),MODE,
                        conversion=conv,profile_points=180)
                    ref["modes"][conv][label]=dict(status="CHECKED",
                         frequency_hz=float(result["frequency_hz"]),
                         matched_residual=float(result["matching_residual"]),
                         native_minus_python_hz=record["resolution"]["8193"]["modes"][conv][label].get("f_hz",float("nan"))-result["frequency_hz"])
                except Exception as e:
                    ref["modes"][conv][label]=dict(status="REJECTED",type=type(e).__name__,
                                                   error=str(e)[:240])
        ref["wall_seconds"]=time.perf_counter()-start
        record["python_DOP853_oracle"]=ref
        print("PYTHON_SHARP_DOP853",name,ref,flush=True)
    return record

def main():
    record=dict(schema="ns_cowling.sharp_cowling_convergence.v1",
       claim_scope="fixed central density .3 fm-3, smooth one-sided EOS and actual sharp native jumps",
       source_approximations_not_tested_by_resolution="Original GPPVA source precision, BBP 3-sig-digit publication, physical conversion time, unknown phase competition",
       cases={},status="RUNNING")
    for name,path,baseline in CASES:
        try:record["cases"][name]=run_one(name,path,baseline)
        except Exception as ex:record["cases"][name]=dict(status="BLOCKED",
              type=type(ex).__name__,error=str(ex)[:450],
              traceback=traceback.format_exc()[-1800:])
        dump(DEST/"PARTIAL.json",record)
    record["status"]="FINISHED_WITH_CASE_STATUSES"
    dump(DEST/"RESULT.json",record)

if __name__=="__main__":main()
