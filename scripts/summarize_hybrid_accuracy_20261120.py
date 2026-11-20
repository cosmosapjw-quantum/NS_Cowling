#!/usr/bin/env python3
"""Select the fastest *measured* physically admissible Cowling backend.

This is a source-scoped comparison, not a global performance model.
Only a previously executed benchmark with actual phasewise native
and independently observed DOP853 frequencies may authorize a choice.
Failure at tighter tolerances remains explicit, never tuned away.
"""
from __future__ import annotations
import hashlib,json,math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data/cowling_hybrid_geometry_20261120"
INPUT=DATA/"RESULT.json"
LIMITS=(.05,.005,.001)


def choose_verified_backend(case, absolute_freq_tolerance_hz):
    if not (isinstance(absolute_freq_tolerance_hz,(int,float))
            and not isinstance(absolute_freq_tolerance_hz,bool)
            and math.isfinite(absolute_freq_tolerance_hz)
            and absolute_freq_tolerance_hz>0):
        raise ValueError("finite positive error threshold in Hz required")
    allowed=[]
    for name,r in case["backend_variants"].items():
        if r.get("status")!="MATCHED_REFERENCE_PASS":
            continue
        # This frozen comparison selects p1, never a faster f/g or
        # unclassified branch with a copied PASS label.
        if r.get("nodes")!=[1,1]:
            continue
        try:
            diff=float(r["reference_delta_hz"])
            t=float(r["total_seconds"])
            join=float(r["join_relative"])
            frequency=float(r["frequency_hz"])
            reference=float(case["reference_hz"])
        except (KeyError,TypeError,ValueError,OverflowError):
            continue
        if (not all(map(math.isfinite,(diff,t,join,frequency,reference))) or
            t<=0 or frequency<=0 or reference<=0 or not 0<=join<=1e-5):
            continue
        # A delta recorded against a different source/reference cannot
        # authorize an apparently accurate backend. Allow only subtraction
        # roundoff, not a physics-level numerical error budget.
        subtraction_roundoff=max(1e-9,32*math.ulp(max(frequency,reference)))
        if abs((frequency-reference)-diff)>subtraction_roundoff:
            continue
        if abs(diff)<=absolute_freq_tolerance_hz:
            allowed.append((t,abs(diff),name,r))
    if not allowed:
        return dict(status="NO_VERIFIED_CANDIDATE_AT_REQUESTED_FREQUENCY_TOLERANCE",
                    requested_abs_error_hz=float(absolute_freq_tolerance_hz),
                    selected=None)
    _,_,name,run=min(allowed,key=lambda x:x[:3])
    return dict(status="EMPIRICAL_CANDIDATE_SELECTED",
                selected=name,backend=run["backend"],
                sample_count=run["samples"],
                frequency_hz=run["frequency_hz"],
                reference_difference_hz=run["reference_delta_hz"],
                elapsed_seconds=run["total_seconds"],
                requested_abs_error_hz=float(absolute_freq_tolerance_hz),
                total_eligible_candidates=len(allowed))


def main():
    raw=INPUT.read_bytes()
    benchmark=json.loads(raw)
    result=dict(schema="ns_cowling.source_coefficient_hybrid_accuracy_selection.v1",
        input_benchmark_sha256=hashlib.sha256(raw).hexdigest(),
        original_source_zip_sha256=benchmark["frozen_EOS_bundle_sha256"],
        criteria="observed reference |df|, physical join <1e-5, strict W/V mode and measured wall time",
        method_selection_is_empirical_not_global_or_proven_optimal=True,
        target_errors_hz=list(LIMITS),cases={})
    for case,record in benchmark["cases"].items():
        values={f"{threshold:.3g} Hz":choose_verified_backend(
            record,threshold) for threshold in LIMITS}
        result["cases"][case]=values
        for tol,r in values.items():
            print("FASTEST_ACCURACY",case,tol,r["status"],
                  r.get("selected"),r.get("elapsed_seconds"),flush=True)
    (DATA/"DECISIONS.json").write_text(
        json.dumps(result,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
