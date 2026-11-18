"""Aggregate direct native EOS interpolation holdout and TF experiments.

Publication comparisons and figures exclude linear interpolation by user request.
Linear is retained only in original diagnostic JSON for provenance.
Methods with low coverage cannot win based on the successfully fitted subset.
"""
from __future__ import annotations
from collections import defaultdict,Counter
from pathlib import Path
import json,math,csv,statistics
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data/native_dense_interpolation_20261118"
OUT=DATA/"aggregated"
PLOT_METHODS=("hermite","exponential","rational","bernstein","power","pchip","akima","cubic")
DIAGNOSTIC_METHODS=("linear",)
assert "linear" not in PLOT_METHODS
MODELS={
 "rmf":("NL3_PROVIDER_PRECISION","PK1_CANONICAL","TM1_CANONICAL",
         "DDMEX_PROVIDER_ROUNDED_N0","DDMEX_ORIGINAL_TABLE_RHO0"),
 "core":("APR_APR","BBB_BHF2","GM_GM1","RG_KDE0v1","RG_SLY4"),
 "tf":("DDME2","DDMEX","NL3","NL3wrL55","TM1e")}
def source_results():
    records=[]
    for family,names in MODELS.items():
        for model in names:
            if family=="tf":
                x=json.loads((DATA/"tf_holdout"/f"{model}.json").read_text())
                rows=[dict(phase=None,train_stride="frozen_H_N128",heldout_size=x["heldout_distinct"],
                           native_physics=x.get("direct_native_fd_diagnostics"),
                           metrics=x["metrics"],case=x["case"])]
                if not rows[0]["metrics"]:raise RuntimeError("TF native holdout has no metrics")
            else:
                x=json.loads((DATA/"cv"/family/model/"RESULT.json").read_text())
                assert x["status"]=="COMPLETED_SCOPED_WITH_NATIVE_HOLDOUT"
                rows=[z for z in x["runs"] if z.get("status")=="COMPLETE_WITH_REJECTED_METHODS"]
                if not rows:raise RuntimeError(f"missing finite native holdout {family}/{model}")
            for row in rows:
                for method in PLOT_METHODS:
                    q=row["metrics"].get(method)
                    if q is None:raise RuntimeError(f"missing metric {method} on {family}/{model}")
                    if q["status"] not in ("COMPUTED_WITH_NATIVE_HOLDOUT",
                                           "EVALUATED_AGAINST_NEW_NATIVE_DROPLET_N128"):
                        raise RuntimeError(f"method crashed {family}/{model} {method}: {q}")
                    n=q["number_density_relative"]
                    cs=q["sound_speed_squared_absolute"]
                    rec=dict(family=family,model=model,phase=row.get("phase"),stride=row["train_stride"],
                        method=method,heldout_count=int(row["heldout_size"]),
                        compared_n_count=n["count"],
                        compared_cs_count=cs["count"],
                        native_causal_count=q["native_causal_holdout"],
                        coverage=float(q["finite_prediction_fraction"]),
                        physical=float(q["physical_among_native_causal"]),
                        baryon_relative_rms=n["rms"],
                        pressure_relative_rms=q["pressure_relative"]["rms"],
                        energy_per_baryon_rms_mev=q["energy_per_baryon_MeV"]["rms"],
                        cs2_rms=cs["rms"],
                        source_acausal_count=row.get("native_physics",{}).get("count_acausal",0),
                        source_cross_form_abs=row.get("native_physics",{}).get("max_abs_cross_form",None))
                    records.append(rec)
    return records
def aggregate(rows):
    group=defaultdict(list)
    for r in rows:group[(r["family"],r["method"])].append(r)
    result={}
    for (family,method),arr in sorted(group.items()):
        # RMS over predictions in all native holdout runs, not a mean of rounded RMS.
        def weighted_rms(key,count_key):
            keep=[x for x in arr if x[key] is not None and x[count_key]>0]
            n=sum(x[count_key] for x in keep)
            return math.sqrt(sum(x[count_key]*(x[key]**2) for x in keep)/n) if n else None
        hold=sum(x["heldout_count"] for x in arr)
        native_valid=sum(x["native_causal_count"] for x in arr)
        sum_cov=sum(x["heldout_count"]*x["coverage"] for x in arr)
        sum_phy=sum(x["native_causal_count"]*x["physical"] for x in arr)
        coverage=sum_cov/hold if hold else 0.
        phys=sum_phy/native_valid if native_valid else 0.
        result.setdefault(family,{})[method]=dict(
            runs=len(arr),heldout_count=hold,
            native_causal_count=native_valid,
            available_count=sum(x["compared_n_count"] for x in arr),
            coverage_fraction=coverage,physical_fraction=phys,
            rms_baryon_relative=weighted_rms("baryon_relative_rms","compared_n_count"),
            rms_pressure_relative=weighted_rms("pressure_relative_rms","compared_n_count"),
            rms_energy_per_baryon_mev=weighted_rms("energy_per_baryon_rms_mev","compared_n_count"),
            rms_cs2=weighted_rms("cs2_rms","compared_cs_count"),
            qualifies_for_global_accuracy_rank=bool(coverage>=.995 and phys>=.995))
    return result
def by_model(records):
    pairs=defaultdict(list)
    for r in records:pairs[(r["family"],r["model"],r["method"])].append(r)
    output=[]
    for (family,model,method),arr in sorted(pairs.items()):
        c=sum(x["compared_n_count"] for x in arr)
        total=sum(x["heldout_count"] for x in arr)
        n_causal=sum(x["native_causal_count"] for x in arr)
        if c==0:rms=None
        else:rms=math.sqrt(sum(x["compared_n_count"]*(x["baryon_relative_rms"]**2) for x in arr if x["baryon_relative_rms"] is not None)/c)
        physical=sum(x["native_causal_count"]*x["physical"] for x in arr)/n_causal if n_causal else 0.
        coverage=sum(x["heldout_count"]*x["coverage"] for x in arr)/total
        output.append(dict(family=family,model=model,method=method,rms_baryon_relative=rms,
                           coverage_fraction=coverage,physical_fraction=physical,
                           available_count=c,heldout_count=total,
                           qualifies=bool(coverage>=.995 and physical>=.995)))
    return output
def main():
    rows=source_results();summary=aggregate(rows);models=by_model(rows)
    tf_orig=json.loads((DATA/"tf"/"PARALLEL_RECEIPT.json").read_text())
    tf_reseed=json.loads((DATA/"tf"/"refined_source_bracket"/"PARALLEL_RECEIPT.json").read_text())
    tf_summary=json.loads((DATA/"tf_holdout"/"SUMMARY.json").read_text())
    rmf=json.loads((DATA/"rmf"/"PARALLEL_RECEIPT.json").read_text())
    core=json.loads((DATA/"core"/"PARALLEL_RECEIPT.json").read_text())
    assert rmf["all_admitted"] and core["success_count"]==5 and tf_summary["status"]=="COMPLETED_ALL_FIVE"
    all_rows=sum(x["native_state_count"] for x in rmf["models"])+sum(x["n_returned"] for x in core["results"])
    assert all_rows==25612
    tf_count=sum(x["new"] for x in tf_summary["cases"])
    assert tf_count==30
    orig_accepted=sum(x["accepted"] for x in tf_orig["results"])
    recovered=sum(x["recovered"] for x in tf_reseed["models"])
    audit=dict(schema="ns_cowling.native_dense_interp_summary.v1",
      scope="Frozen-model native high-resolution holdout: 5 RMF, 5 core, 5 TF source branches",
      original_h_commit="a0fa5fb93851a9647030f3589c24a62ddea527b1",
      generated_direct_native_state_count=all_rows,generated_uniform_rmf_models=5,
      generated_core_models=5,TF_requested_new_density=57,
      TF_direct_native_accepted=orig_accepted,TF_reseed_accepted=recovered,
      TF_distinct_new_withheld=30,
      TF_native_rejected_after_attempts=57-orig_accepted-recovered,
      interpolation_methods_plotted=list(PLOT_METHODS),
      excluded_plot_methods=list(DIAGNOSTIC_METHODS),
      accuracy_comparison="P(mu), n(mu), epsilon(mu) at withheld native chemical potentials; phase-separated with no APR phase mixing",
      physically_valid="Sampled withheld states pass c_s^2/c^2>0 and <=1, conditional on native source being causal; accepted shape interval validated continuously",
      source_limitations=[
          "Native uniform homogeneous branch is not proven ground-state below crust transition",
          "APR 6/8/7 phases are separately interpolated; numerically acausal native high-density states are flagged, not fitted as causal truth",
          "BBB_BHF2 underlying published energy PCHIP reconstruction is not a new microscopic BHF calculation",
          "TF finite-N128 droplet stationary points are not a global minimum across pasta geometries",
          "For first-order Maxwell interfaces, phase blending is prohibited; root-discovery sensitivity is a separate experiment",
          "Numerical cs2 reference is finite-difference based, with derivative error distinct from native thermodynamic values",
          "Low-coverage methods cannot win ranking based only on the successfully covered subset"],
      family_summary=summary,by_model=models,
      native_rmf_receipt="rmf/PARALLEL_RECEIPT.json",
      native_core_receipt="core/PARALLEL_RECEIPT.json",
      native_tf_receipt="tf/PARALLEL_RECEIPT.json",
      model_holdout_results="cv/{rmf|core}/{model}/RESULT.json and tf_holdout/{case}.json")
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"SUMMARY.json").write_text(json.dumps(audit,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    with (OUT/"BY_MODEL.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(models[0]))
        w.writeheader();w.writerows(models)
    with (OUT/"BY_RUN.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    for fam,methods in summary.items():
        ranking=sorted([(m,x) for m,x in methods.items() if x["qualifies_for_global_accuracy_rank"]],
                       key=lambda t:t[1]["rms_baryon_relative"])
        print(f"{fam}: ALL METHODS={len(methods)} eligible={len(ranking)}",
              "best",[(m,round(q["rms_baryon_relative"]*1e6,5)) for m,q in ranking],flush=True)
    print("NATIVE_STATE_COUNT",all_rows,"TF_HOLDOUT",tf_count,"PLOTTED_METHODS",list(PLOT_METHODS))
    return 0
if __name__=="__main__":
    raise SystemExit(main())
