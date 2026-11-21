#!/usr/bin/env python3
"""Frontend-only parsing of C14 original-EOS native-Fortran error budget.

No Python eigenproblem, EOS physical evaluation, interpolation or roots.
Summarizes independently controlled Fortran calculations and reports
non-asymptotic behavior without claiming Richardson error bars.
"""
from __future__ import annotations
import gzip,hashlib,json,math,re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/cowling_c14_source_error_budget_20261121"
SOURCES=ROOT/"data/cowling_fortran_source_exact_root_c12_20261121"
ARCHIVE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SOURCE_SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
C12_INPUT_SHAS={
 "M2":"eee499eb50de4300f9987e05bf46a4fee3ddcbfb2e9b9d99deb9ac8a423c35e5",
 "M3":"f6a7abb566ae68b8ac5f838d410771823634a98024007ca02264bd6e44a0d3b1",
}

def decode_row(line):
    t=line.split()
    if not t or not t[0].startswith("C14_ROW"):
        return None
    axis=t[0][len("C14_ROW"):]
    if axis not in ("GEOMETRY","RTOL","MAXSTEP") or len(t)!=13:
        raise RuntimeError(f"unrecognized Fortran C14 scientific row: {line}")
    factor,nh,nl,mf=map(int,t[1:5])
    rtol,maxds,frequency,det=map(float,t[5:9])
    iters,nfev,ncuts,status=map(int,t[9:13])
    if not all(map(math.isfinite,(rtol,maxds,frequency,det))):
        raise RuntimeError("nonfinite Fortran C14 result")
    return dict(axis=axis,coarsen_factor=factor,n_hi=nh,n_lo=nl,
        mf=mf,relative_tolerance=rtol,max_step_s=maxds,
        frequency_hz=frequency,determinant=det,
        brent_iterations=iters,native_rhs_evaluations=nfev,
        source_knots=ncuts,status=status)

def summarize(rows,axis,mf):
    chosen=[r for r in rows if r["axis"]==axis and r["mf"]==mf]
    assert len(chosen)==4
    chosen=sorted(chosen,
      key=lambda r:r["coarsen_factor"] if axis=="GEOMETRY" else
       -r["relative_tolerance"] if axis=="RTOL" else
       -r["max_step_s"])
    vals=[r["frequency_hz"] for r in chosen if r["status"]==0]
    if len(vals)!=len(chosen):
        return dict(status="PARTIAL_FAILURE",attempted=len(chosen),
           accepted=len(vals),statuses=[r["status"] for r in chosen])
    diffs=[vals[i+1]-vals[i] for i in range(len(vals)-1)]
    # Nonmonotone means at least one signed consecutive refinement
    # increment changes sign. This is descriptive, not a convergence proof.
    nonmono=any(x*y<0 for x,y in zip(diffs[:-1],diffs[1:]))
    return dict(status="ALL_VALID",
      frequency_span_hz=max(vals)-min(vals),
      endpoint_difference_hz=vals[-1]-vals[0],
      max_change_per_neighbor_hz=max(abs(x) for x in diffs),
      frequencies_hz=vals,
      successive_deltas_hz=diffs,
      signed_step_nonmonotone=nonmono,
      max_native_rhs_evals=max(r["native_rhs_evaluations"] for r in chosen),
      min_native_rhs_evals=min(r["native_rhs_evaluations"] for r in chosen),
      all_source_knots_equal=len(set(r["source_knots"] for r in chosen))==1)

def run():
    if hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()!=SOURCE_SHA:
        raise RuntimeError("original EOS archive SHA mismatch")
    out=dict(schema="ns_cowling.c14.original_eos_native_separated_error_budget.v1",
      parent_commit="e0f6209f5c54b08b80959f9bf1dfe5742240bc2f",
      source_archive_sha256=SOURCE_SHA,
      physical_scope="cold equilibrium Gamma1=Gamma_eq, one Maxwell jump, original source cs2",
      fixed_parameters=["original EOS source","TOV background and phase endpoints",
        "Maxwell conversion","center regularity","surface DeltaP=0",
        "root_xtol_hz and absolute ODE tolerance"],
      independent_axes=["geometry coarsening with source knots pinned",
                        "DVODE relative_tolerance","DVODE maximum_step_s"],
      numerical_precision_claim="descriptive finite-control stability, not absolute error bound",
      models={})
    for model in ("M2","M3"):
        orig=gzip.decompress((SOURCES/f"ORIGINAL_FORTRAN_ROOT_{model}.dat.gz").read_bytes())
        if hashlib.sha256(orig).hexdigest()!=C12_INPUT_SHAS[model]:
            raise RuntimeError("original source Fortran contract does not match C12")
        data=(OUT/"OPT_1"/f"{model}.science").read_text().splitlines()
        rows=[r for r in map(decode_row,data) if r is not None]
        if len(rows)!=24:raise RuntimeError(f"expected 24 controls for {model}; got {len(rows)}")
        if len({(r["axis"],r["coarsen_factor"],r["mf"],
                 r["relative_tolerance"],r["max_step_s"]) for r in rows})!=24:
            raise RuntimeError("duplicate control detected")
        # A model row with status>0 remains an explicit failed scientific
        # observation and does not get a made-up frequency.
        for r in rows:
            if r["status"]==0 and not (
               r["brent_iterations"]>0 and r["native_rhs_evaluations"]>0
               and r["source_knots"]>0 and abs(r["determinant"])<=1e-5):
                raise RuntimeError("invalid certified Fortran root control")
            if r["status"]!=0 and r["frequency_hz"]!=0:
                raise RuntimeError("invalid control falsely has certified frequency")
        grouped={axis:{str(mf):summarize(rows,axis,mf) for mf in (10,22)}
            for axis in ("GEOMETRY","RTOL","MAXSTEP")}
        out["models"][model]=dict(
            original_source_input_uncompressed_sha256=C12_INPUT_SHAS[model],
            rows=rows,summary=grouped,
            controls=24,successes=sum(r["status"]==0 for r in rows))
        print("SOURCE_C14",model,"SUCCESS",out["models"][model]["successes"])
        for axis in grouped:
            for mf,s in grouped[axis].items():
                print("  ",axis,"mf",mf,"span_Hz",s.get("frequency_span_hz"),
                    "signed_nonmono",s.get("signed_step_nonmonotone"))
    out["total_fortran_original_source_roots"]=48
    out["status"]="SCOPED_NUMERICAL_STABILITY_CHARACTERIZED_NOT_ABSOLUTE_CONVERGENCE"
    path=OUT/"C14_ERROR_BUDGET.json"
    path.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("C14_ANALYSIS_WRITTEN",path)

if __name__=="__main__":
    run()
