#!/usr/bin/env python3
"""C15 evidence-based original-source TOV-to-Cowling sensitivity analysis.

Read only immutable original-source numeric outputs, preserve TOV
and Cowling root/mode failure statuses, and compare finite-control spans.
No Python physical/ODE solver here. Never extrapolate from non-asymptotic
finite differences into a certified frequency error.
"""
from __future__ import annotations
import hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"data/cowling_c15_tov_convergence_20261121"
SRC=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
EXPECTED={
 "M2":"59b7cb4671347c33c25133e6497e8be346b4e5b4f871eece9bcaf59155fd8765",
 "M3":"f23fb7020c9de1b65c99ac37690840f7a4d6611f3372409b015d8c3bd7e32f3b"}
C14_REF=ROOT/"data/cowling_c14_source_error_budget_20261121/C14_ERROR_BUDGET.json"

def axis_rows(rows,axis,method):
    base=next(r for r in rows if r["TOV_source_samples"]==2049 and
      r["axis"]=="TOV_RTOL" and r["TOV_method"]==method and
      r["TOV_rtol"]==3e-10 and r["TOV_max_step_h"]==.004)
    if axis=="TOV_THERMODYNAMIC_SOURCE_GRID":
        sub=[r for r in rows if r["axis"]==axis and r["TOV_method"]==method]
        sub.append(base)
        return sorted(sub,key=lambda r:r["TOV_source_samples"])
    if axis=="TOV_RTOL":
        return sorted([r for r in rows if r["axis"]==axis and r["TOV_method"]==method],
             key=lambda r:-r["TOV_rtol"])
    if axis=="TOV_MAX_STEP_H":
        sub=[r for r in rows if r["axis"]==axis and r["TOV_method"]==method]
        sub.append(base)
        return sorted(sub,key=lambda r:-r["TOV_max_step_h"])
    raise ValueError("unknown original TOV axis")

def span(values):
    return max(values)-min(values) if len(values)>1 else None

def oneaxis(rows,axis,method):
    ordered=axis_rows(rows,axis,method)
    if len(ordered)!=4:raise RuntimeError("incorrect dimensional control count")
    success=[r for r in ordered if r["mode_status"]=="PASS_P1"]
    stars=[r for r in ordered if r["TOV_status"]==0]
    fs=[r["frequency_hz"] for r in success]
    ms=[r["TOV_mass_msun"] for r in stars]
    rs=[r["TOV_radius_km"] for r in stars]
    steps=[r["TOV_source_samples"] if axis=="TOV_THERMODYNAMIC_SOURCE_GRID" else
           r["TOV_rtol"] if axis=="TOV_RTOL" else r["TOV_max_step_h"] for r in ordered]
    changes=[fs[i+1]-fs[i] for i in range(len(fs)-1)]
    return dict(axis=axis,method=method,
       controls=steps,
       statuses=[r["mode_status"] for r in ordered],
       tov_statuses=[r["TOV_status"] for r in ordered],
       mode_admissions=len(success),tov_valid=len(stars),
       frequency_span_admitted_hz=span(fs),
       tov_mass_span_valid_msun=span(ms),
       tov_radius_span_valid_km=span(rs),
       ordered_frequencies_hz=[r.get("frequency_hz") for r in ordered],
       ordered_masses_msun=[r.get("TOV_mass_msun") for r in ordered],
       ordered_radii_km=[r.get("TOV_radius_km") for r in ordered],
       signed_refinement_changes_admitted_hz=changes,
       signed_frequency_nonmonotone=(any(x*y<0 for x,y in zip(changes[:-1],changes[1:]))
          if len(success)==4 else None),
       max_join_among_admitted=max([r["relative_join_mismatch"] for r in success],default=None),
       source_knot_counts=sorted(set(r["source_restart_count"] for r in success)))

def run():
    if hashlib.sha256(SRC.read_bytes()).hexdigest()!=SHA:
        raise RuntimeError("original source EOS archive drift")
    a=(D/"REAL_P1_CONTIGUOUS/SCIENTIFIC.sha256").read_bytes()
    assert a==(D/"REAL_P2_CONTIGUOUS/SCIENTIFIC.sha256").read_bytes()
    assert a==(D/"REAL_P2/SCIENTIFIC.sha256").read_bytes()
    c14=json.loads(C14_REF.read_text())
    result=dict(schema="ns_cowling.c15.source_locked_TOV_to_Cowling_sensitivity.v1",
       prior_base_commit="091c2f6f22d032379324aa900c910c6278a68437",
       source_EOS_archive_sha256=SHA,
       numerical_controls="TOV grid, RTOL, max_step_h, two methods; Cowling identical",
       invariant="same original EOS, central density, Maxwell interface, Gamma1=Gamma_eq and Cowling controls",
       observations_are_not_bounds=True,models={})
    for key in ("M2","M3"):
        dat=json.loads((D/f"REAL_P1_CONTIGUOUS/{key}.json").read_text())
        assert dat["source_archive_sha256"]==SHA
        assert dat["source_contract_sha256"]==EXPECTED[key]
        assert len(dat["all_rows"])==20
        rows=dat["all_rows"]
        accepted=[r for r in rows if r["mode_status"]=="PASS_P1"]
        rejected=[r for r in rows if r["mode_status"]=="NUMERICALLY_OR_PHYSICALLY_REJECTED"]
        tov_failed=[r for r in rows if r["TOV_status"]!=0]
        assert all(r["W_nodes"]==r["V_nodes"]==1 and
                   r["relative_join_mismatch"]<1e-5 and
                   r["source_callback_rhs_calls"]==0 for r in accepted)
        assert all("frequency_hz" not in r for r in rejected+tov_failed)
        assert all(r["TOV_status"]==0 for r in accepted+rejected)
        groups={axis:{method:oneaxis(rows,axis,method) for method in ("adams","bdf")}
           for axis in ("TOV_THERMODYNAMIC_SOURCE_GRID","TOV_RTOL","TOV_MAX_STEP_H")}
        record=dict(original_EOS_SHA256=dat["source_contract_sha256"],
            original_center_n_fm3=dat["star_center_n_fm3"],
            source_model=dat["model"],
            conversion=dat["conversion"],
            original_reference_M_msun=dat["TOV_original_reference_M_msun"],
            original_reference_R_km=dat["TOV_original_reference_R_km"],
            admitted_p1=len(accepted),
            rejected_cowling_roots=len(rejected),
            failed_TOV_backgrounds=len(tov_failed),
            baseline_source_locked_frequency_hz=dat["baseline_frequency_hz"],
            axes=groups,
            rejection_codes=[dict(axis=r["axis"],rtol=r["TOV_rtol"],
              max_step_h=r["TOV_max_step_h"],method=r["TOV_method"],
              source_samples=r["TOV_source_samples"],status=r["mode_status"],
              TOV_status=r["TOV_status"],reason=r.get("reason")) for r in rejected+tov_failed],
            parent_C14_fixed_TOV_frequency_spans_hz=c14["models"][key]["summary"])
        result["models"][key]=record
        print("C15_NUMERIC",key,"PASSED",len(accepted),
          "TOV_FAILURE",len(tov_failed),"ROOT_REJECTION",len(rejected),flush=True)
        for axis,modes in groups.items():
            for method,g in modes.items():
                print("  AXIS",axis,method,"Admitted",g["mode_admissions"],
                  "f_span",g["frequency_span_admitted_hz"],
                  "M_span",g["tov_mass_span_valid_msun"],
                  "R_span",g["tov_radius_span_valid_km"],flush=True)
    result["status"]="PASS_SCOPED_TOV_BACKGROUND_SENSITIVITY_GLOBAL_FREQUENCY_CONVERGENCE_UNRESOLVED"
    out=D/"TOV_COWLING_SENSITIVITY_SUMMARY.json"
    out.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("C15_SUMMARY_SAVED",out,flush=True)

if __name__=="__main__":
    run()
