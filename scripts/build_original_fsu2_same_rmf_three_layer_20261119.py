#!/usr/bin/env python3
"""Compose immutable GPPVA(FSU2) BSk22 outer, FSU2 inner and same-FSU2
homogeneous native source using minimal causal lower connector.

Test finite native Maxwell and CSS as distinct physical constructions;
approve three-layer stellar model only after strict WholeCrustBarotrope
contract, native Fortran TOV, Cowling slow/rapid modes and source readback.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path
import traceback

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos import read_compose, match_crust_core
from compose_eos.api import construct_maxwell_transition, _thermo_only_view
from compose_eos.maxwell_v2 import _segment_eosdata
from compose_eos.whole_crust import SurfacePolytrope, PiecewisePolytrope, CausalExponentialInterval, WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS, solve_sharp_background, solve_mode as solve_python_sharp_mode
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE, solve_sharp_background_native
from compose_eos.stellar_cowling import CowlingSettings

OUT=ROOT/"data/cowling_maxwell_native_20261119/three_layer_original_fsu2"
ORIGINAL=ROOT/"data/maxwell_three_families_20261117/intake/references/GPPVA_FSU2"
RECON=ROOT/"data/cowling_maxwell_native_20261119/same_rmf/generated_core/GPPVA_FSU2"
SURVEY=ROOT/"data/cowling_maxwell_native_20261119/same_rmf/SURVEY.json"

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n")

def mode_result(whole,source,mode_label,nc):
    """Run actual native Fortran phasewise TOV then two Cowling junction types."""
    eos=SharpContractEnthalpyEOS(whole)
    t0=time.perf_counter()
    opts=dict(surface_h=3e-8,rtol=3e-9,atol=3e-11,max_step_h=.004)
    star=solve_sharp_background_native(eos,nc,native_samples=4097,**opts)
    bgsec=time.perf_counter()-t0
    # Separate the source provenance from an independent Python DOP853 TOV cross-check.
    py=solve_sharp_background(eos,nc,**opts)
    cfg=CowlingSettings(surface_h=3e-8,rtol=3e-9,atol=3e-12,max_step_h=.001)
    solver=NativeSharpCowlingDVODE(star,cfg,samples=4097)
    records={}
    for conv in ("slow","rapid"):
        brackets,scan=solver.scan_brackets(np.linspace(1000.,11500.,85),conversion=conv)
        modes=[]
        for bracket in brackets[:6]:
            try:
                m=solver.solve_mode(bracket,conversion=conv,profile_points=160)
                if (m["W_nodes_per_phase"],m["V_nodes_per_phase"]) not in ((0,0),(1,1)):
                    continue
                modes.append(dict(label="f" if m["W_nodes_per_phase"]==0 else "p1",
                    frequency_hz=m["frequency_hz"],matched_relative=m["matching_residual"],
                    join_relative=m["relative_join_mismatch"],
                    surface_BC_relative=m["surface_BC_relative"],
                    interface_count=len(m["interface_junctions"]),
                    interface_max_traction=max((x["traction_relative_residual"]
                                               for x in m["interface_junctions"]),default=0),
                    W_nodes=m["W_nodes_per_phase"],V_nodes=m["V_nodes_per_phase"]))
                if len({x["label"] for x in modes})>=2:break
            except RuntimeError as e:
                # Retain numerical failures; do not admit fake sign jumps.
                print("MODE_REJECTED",source,mode_label,conv,bracket,str(e),flush=True)
        records[conv]=dict(status="NATIVE_MODES_VERIFIED" if len({x["label"] for x in modes})==2
                           else "PARTIAL_MODE_SET",modes=modes,detected_brackets=len(brackets))
    return dict(native_tov=dict(mass_msun=star.mass_msun,R_km=star.radius_km,
                 interfaces=star.interfaces,seconds=bgsec),
        independent_python_tov=dict(mass_msun=py.mass_msun,R_km=py.radius_km),
        mass_delta_msun=star.mass_msun-py.mass_msun,
        radius_delta_km=star.radius_km-py.radius_km,
        modes=records,source=source,matching=mode_label)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    raw=read_compose(ORIGINAL)
    high=read_compose(RECON)
    # exact original attachment: phase0 BSk22 outer and phase7 TF inner.
    edge=int(np.flatnonzero(raw.Iphase!=0)[0])
    assert np.all(raw.Iphase[:edge]==0) and np.all(raw.Iphase[edge:]==7)
    original_outer=_segment_eosdata(raw,(0,0,edge))
    original_inner=_segment_eosdata(raw,(0,edge,raw.n))
    outer_piece=PiecewisePolytrope(original_outer.nb,original_outer.p,
                                    original_outer.eps[0]/original_outer.nb[0])
    outer_validation=outer_piece.validate()
    assert outer_validation["status"]=="PASS"
    surface=SurfacePolytrope.from_anchor(original_outer.nb[0],
                                          original_outer.p[0],original_outer.eps[0])
    # Source n,P exactly unchanged in the outer layer. Only epsilon and mu
    # are re-derived by first-law integration anchored at the first original.
    deltaF=outer_piece.F-original_outer.eps/original_outer.nb
    original_report=dict(total_source_rows=raw.n,outer_rows=edge,inner_rows=original_inner.n,
        outer_last_density=float(original_outer.nb[-1]),
        inner_first_density=float(original_inner.nb[0]),
        unrepresented_density_gap=float(original_inner.nb[0]-original_outer.nb[-1]),
        surface_n0=float(surface.n0),surface_rest_energy_MeV=float(surface.rest_energy_MeV),
        max_outer_energy_per_baryon_correction_MeV=float(max(abs(deltaF))),
        source_original_input_hashes={f.name:sha(f) for f in ORIGINAL.glob("eos.*") if f.is_file()},
        generated_homogeneous_core_hashes={f.name:sha(f) for f in RECON.glob("eos.*") if f.is_file()},
        outer_source_preserved_fields="all n_B and P values; epsilon/mu thermodynamically corrected, zero external energy shift",
        original_lower_matched_source="original provider BSk22 outer + original FSU2 inner; a causal phenomenological connector closes their actual missing support",
        original_outer_density_firstlast=[float(original_outer.nb[0]),float(original_outer.nb[-1])])
    record=dict(schema="ns_cowling.original_fsu2_three_layer_same_rmf.v1",
          lower_layers=original_report,methods={},
          model_scope="three-layer locally neutral effective cold neutron star, no microscopic unification/global ground-state claim",
          mapping="P/epsilon continuity for lower connector, actual Maxwell/CSS jump; do not smooth jump",
          success_gate="actual native TOV and f/p1 Cowling both slow and rapid conversions")
    end=(float(outer_piece.nb[-1]),float(outer_piece.p[-1]),
         float(outer_piece.nb[-1]*outer_piece.F[-1]))
    outer_contract=outer_piece.contract()
    for method in ("maxwell","css"):
        entry=dict(status="RUNNING",matching=method,
           source_EOS_hashes=original_report["source_original_input_hashes"],
           same_functional_homogeneous_core=True,
           microscopic_unification_required=False)
        record["methods"][method]=entry
        dump(OUT/"PARTIAL.json",record)
        try:
            prepared_source=ROOT/"data/maxwell_three_families_20261117/prepared_v2/public_original/GPPVA_FSU2/compose"
            prepared_original=read_compose(prepared_source)
            for field in ("nb","p","eps","mub"):
                if not np.array_equal(getattr(prepared_original,field),getattr(original_inner,field)):
                    raise RuntimeError("prepared original source changed thermodynamic field: "+field)
            # The source-phase label and provider-specific extra fields remain
            # preserved in immutable original; the existing prepared source
            # provides the documented thermodynamics-only native phase view.
            low=_thermo_only_view(prepared_original)
            core=_thermo_only_view(high)
            if method=="maxwell":
                matched,info=construct_maxwell_transition(low,core,
                    assume_local_charge_neutrality=True,drop_composition=True,
                    return_report=True)
                j=int(info["low_phase_rows"])-1
                pt=float(info["pressure_transition"])
            else:
                matched,info=match_crust_core(low,core,"GPPVA_FSU2",
                    "GPPVA_FSU2","css",40,drop_composition=True)
                j=int(info["css_transition_low_index"])
                pt=float(info["css_transition_pressure"])
            n,p,e=matched.nb,matched.p,matched.eps
            assert info["consistency"]["overall"]["pass_"]
            assert 0<=j<len(n)-1
            assert n[j+1]>n[j]>0
            assert abs(p[j]-p[j+1])<2e-11*max(1,pt)
            entry.update(finite_status="PASS_NATIVE",matched_rows=matched.n,
                 transition_pressure=pt,transition_n=[float(n[j]),float(n[j+1])],
                 transition_jump_energy=float(e[j+1]-e[j]),
                 original_matching_low_pressure=float(original_inner.p[0]),
                 matching_report_kind=method)
            acceptance=None;rejections=[]
            # Iteratively omit only as many source inner points as necessary
            # to admit a thermodynamically causal lower connector.
            for k in range(min(50,j-1)):
                if n[k]<=end[0]:continue
                right=tuple(float(arr[k]) for arr in (n,p,e))
                try:
                    bridge=CausalExponentialInterval(*end,*right)
                    bpass=bridge.validate()
                    if bpass["status"]!="PASS":raise ValueError("lower causal bridge rejected")
                    contract=dict(surface=surface.contract(),
                       outer=outer_contract,lower_bridge=bridge.contract(),
                       upper=dict(n=list(map(float,n[k:])),
                           p=list(map(float,p[k:])),
                           epsilon=list(map(float,e[k:])),
                           sharp_interfaces=[j-k]))
                    whole=WholeCrustBarotrope.from_contract(contract)
                    # Explicitly verify 3 phases, same-func core and outer support.
                    ntest=float(whole.outer.nb[0])
                    if not whole.pressure_domain[1]>pt>0:raise ValueError("invalid upper P domain")
                    for pressure in (whole.outer.p[0],whole.outer.p[-1],bridge.p_right):
                        z=whole.state_at_pressure(float(pressure))
                        if not np.all(np.isfinite([z["p"],z["epsilon"],z["cs2_over_c2"]])):raise ValueError("nonfinite bound state")
                    acceptance=dict(join_index_k=k,
                        discarded_original_inner_rows=k,
                        join_n=float(n[k]),join_P=float(p[k]),bridge=bpass)
                    outdir=OUT/method
                    outdir.mkdir(exist_ok=True)
                    dump(outdir/"REPRESENTATION.json",contract)
                    # Choose a finite high-central-density point supported by the
                    # untouched generated homogeneous source.
                    result=mode_result(whole,"GPPVA_FSU2_original",method,.30)
                    acceptance["stellar"]=result
                    break
                except Exception as ex:
                    rejections.append(dict(k=k,type=type(ex).__name__,error=str(ex)[:200]))
            entry.update(status="THREE_LAYER_NATIVE_TOV_COWLING_CHECKED" if acceptance else
                "THREE_LAYER_NOT_ADMITTED",admitted=acceptance,
                first_rejections=rejections[:12],rejection_count=len(rejections))
            print("THREE_LAYER_METHOD",method,entry["status"],"first_k",
                  acceptance["join_index_k"] if acceptance else None,flush=True)
        except Exception as ex:
            entry.update(status="METHOD_FAILED",failure=dict(type=type(ex).__name__,error=str(ex)),
                traceback=traceback.format_exc()[-1300:])
            print("METHOD_ERROR",method,type(ex).__name__,str(ex),flush=True)
        dump(OUT/"PARTIAL.json",record)
    record["status"]="RUN_COMPLETED_WITH_SEPARATE_METHOD_AND_STELLAR_GATES"
    dump(OUT/"RESULT.json",record)
    print("THREE_LAYER_FINISH",[(k,v["status"]) for k,v in record["methods"].items()],flush=True)

if __name__=="__main__":
    main()
