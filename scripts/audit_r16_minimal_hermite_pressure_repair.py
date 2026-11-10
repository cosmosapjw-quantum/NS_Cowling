#!/usr/bin/env python3
"""Compare current q-linear projection with a diagnostic minimal Hermite pressure repair."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from r16_minimal_hermite_pressure_repair import (  # noqa: E402
    audit_hermite_cell,
    minimal_admissible_delta_p,
)

SUMMARY=ROOT/"oscillation_release_v1/data/spectrum_summary_234.csv"
SOURCE_DIR=ROOT/"oscillation_release_v1/data/eos_original_normalized"
PROJECTED_DIR=ROOT/"oscillation_release_v1/data/eos_projected"
SUMMARY_SHA256="ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965"
EXPECTED_TABLES=234
REPAIR_MARGIN_FRACTION=1.0e-10
REQUIRED=("source_row","mu_b_mev","pressure_mev_fm3","energy_mev_fm3","n_b_fm3")


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def read_rows(path:Path):
    with path.open(newline="",encoding="utf-8") as stream:
        reader=csv.DictReader(stream)
        missing=[x for x in REQUIRED if x not in (reader.fieldnames or ())]
        if missing: raise ValueError(f"{path.name}: missing columns {missing}")
        rows=list(reader)
    if len(rows)<3: raise ValueError(f"{path.name}: too few rows")
    for i,row in enumerate(rows,start=1):
        if int(row["source_row"])!=i: raise ValueError(f"{path.name}: source_row order mismatch")
    return rows


def relative_max(a,b):
    return max(abs(x-y)/max(abs(y),1e-300) for x,y in zip(a,b))


def audit_one(ref):
    name=ref["name"]
    src=SOURCE_DIR/f"{name}.csv"
    proj=PROJECTED_DIR/f"{name}.csv"
    source_hash=sha256(src)
    projected_hash=sha256(proj)
    if source_hash!=ref["raw_normalized_eos_sha256"]:
        raise ValueError(f"{name}: normalized source hash mismatch")
    if projected_hash!=ref["projected_eos_sha256"]:
        raise ValueError(f"{name}: projected source hash mismatch")
    sr=read_rows(src)
    pr=read_rows(proj)
    if len(sr)!=len(pr): raise ValueError(f"{name}: row count mismatch")

    mu=[float(x["mu_b_mev"]) for x in sr]
    n=[float(x["n_b_fm3"]) for x in sr]
    p=[float(x["pressure_mev_fm3"]) for x in sr]
    eps=[float(x["energy_mev_fm3"]) for x in sr]
    p_proj=[float(x["pressure_mev_fm3"]) for x in pr]
    n_proj=[float(x["n_b_fm3"]) for x in pr]
    mu_proj=[float(x["mu_b_mev"]) for x in pr]

    p_repaired=[p[0]]
    repair_records=[]
    source_hermite_failures=0
    source_status_counts={"REJECT_NONPOSITIVE_COMPRESSIBILITY":0,"REJECT_ACAUSAL_NATIVE_CELL":0,"PASS_NATIVE_CELL_ADMISSIBLE":0}
    min_d2=float("inf"); min_g=float("inf")
    max_abs_local_repair=0.0
    sum_abs_increment_repair=0.0
    adjusted=0

    for i in range(len(mu)-1):
        source_cell=audit_hermite_cell(
            mu0=mu[i],p0=p[i],n0=n[i],
            mu1=mu[i+1],p1=p[i+1],n1=n[i+1]
        )
        source_status_counts[source_cell["status"]]+=1
        if source_cell["status"]!="PASS_NATIVE_CELL_ADMISSIBLE":
            source_hermite_failures+=1

        delta_source=p[i+1]-p[i]
        repair=minimal_admissible_delta_p(
            mu0=mu[i],n0=n[i],mu1=mu[i+1],n1=n[i+1],
            delta_p_source=delta_source,
            interior_margin_fraction=REPAIR_MARGIN_FRACTION,
        )
        delta_star=repair["delta_p_repaired"]
        adjusted += int(repair["was_adjusted"])
        sum_abs_increment_repair += abs(repair["delta_p_adjustment"])
        max_abs_local_repair=max(max_abs_local_repair,abs(repair["delta_p_adjustment"])/max(abs(delta_source),1e-300))
        p_repaired.append(p_repaired[-1]+delta_star)

        cell=audit_hermite_cell(
            mu0=mu[i],p0=p_repaired[i],n0=n[i],
            mu1=mu[i+1],p1=p_repaired[i+1],n1=n[i+1]
        )
        if cell["status"]!="PASS_NATIVE_CELL_ADMISSIBLE":
            raise RuntimeError(f"{name}: repaired cell {i+1} not admitted: {cell}")
        min_d2=min(min_d2,cell["min_d2p"])
        min_g=min(min_g,cell["min_causal_margin"])
        repair_records.append((repair,cell))

    eps_repaired=[m*nn-pp for m,nn,pp in zip(mu,n,p_repaired)]
    if any(x<=0.0 for x in eps_repaired):
        raise RuntimeError(f"{name}: repaired energy density became nonpositive")
    if any(b<=a for a,b in zip(eps_repaired,eps_repaired[1:])):
        raise RuntimeError(f"{name}: repaired energy density is not monotone")
    if any(b<=a for a,b in zip(p_repaired,p_repaired[1:])):
        raise RuntimeError(f"{name}: repaired pressure is not monotone")

    eps_projected=[float(x["energy_mev_fm3"]) for x in pr]
    pmax=max(p)
    pdelta_source=[a-b for a,b in zip(p_repaired,p)]
    pdelta_projected=[a-b for a,b in zip(p_repaired,p_proj)]
    epsdelta_source=[a-b for a,b in zip(eps_repaired,eps)]
    epsdelta_projected=[a-b for a,b in zip(eps_repaired,eps_projected)]
    max_local_pressure_source=max(abs(d)/max(abs(x),1e-300) for d,x in zip(pdelta_source,p))
    max_local_pressure_projected=max(abs(d)/max(abs(x),1e-300) for d,x in zip(pdelta_projected,p_proj))
    source_unchanged=sha256(src)==source_hash
    projected_unchanged=sha256(proj)==projected_hash
    if not source_unchanged or not projected_unchanged:
        raise RuntimeError(f"{name}: source/projected bytes changed")

    return {
        "name":name,"crust":ref["crust"],"core":ref["core"],"matching_method":ref["matching_method"],
        "row_count":len(sr),"cell_count":len(sr)-1,
        "source_sha256":source_hash,"projected_sha256":projected_hash,
        "source_unchanged":source_unchanged,"projected_unchanged":projected_unchanged,
        "source_hermite_pass_cells":source_status_counts["PASS_NATIVE_CELL_ADMISSIBLE"],
        "source_hermite_nonpositive_curvature_cells":source_status_counts["REJECT_NONPOSITIVE_COMPRESSIBILITY"],
        "source_hermite_acausal_cells":source_status_counts["REJECT_ACAUSAL_NATIVE_CELL"],
        "source_hermite_reject_cells":source_hermite_failures,
        "minimal_repair_adjusted_cells":adjusted,
        "source_pressure_increments_preserved_cells":len(sr)-1-adjusted,
        "minimal_repair_max_relative_increment_change":max_abs_local_repair,
        "minimal_repair_sum_abs_increment_change":sum_abs_increment_repair,
        "minimal_repair_max_abs_pressure_delta_over_source_Pmax":max(abs(x) for x in pdelta_source)/pmax,
        "current_projector_max_abs_pressure_delta_over_source_Pmax":max(abs(a-b) for a,b in zip(p_proj,p))/pmax,
        "repair_vs_current_projector_max_abs_pressure_delta_over_source_Pmax":max(abs(x) for x in pdelta_projected)/pmax,
        "repair_vs_current_projector_rms_pressure_delta_over_source_Pmax":math.sqrt(sum(x*x for x in pdelta_projected)/len(pdelta_projected))/pmax,
        "minimal_repair_max_local_pressure_delta_pct":100*max_local_pressure_source,
        "repair_vs_projector_max_local_pressure_delta_pct":100*max_local_pressure_projected,
        "minimal_repair_max_abs_energy_delta_over_source_eps":max(abs(x)/max(abs(y),1e-300) for x,y in zip(epsdelta_source,eps)),
        "current_projector_max_abs_energy_delta_over_source_eps":relative_max(eps_projected,eps),
        "repair_vs_projector_max_abs_energy_delta_over_source_eps":max(abs(x)/max(abs(y),1e-300) for x,y in zip(epsdelta_projected,eps_projected)),
        "repaired_min_energy_mev_fm3":min(eps_repaired),
        "repaired_min_delta_energy_mev_fm3":min(b-a for a,b in zip(eps_repaired,eps_repaired[1:])),
        "repaired_min_cell_d2p":min_d2,
        "repaired_min_causal_margin":min_g,
        "source_row217_pressure":p[216],"source_row218_pressure":p[217],
        "repaired_row217_pressure":p_repaired[216],"repaired_row218_pressure":p_repaired[217],
        "projected_row217_pressure":p_proj[216],"projected_row218_pressure":p_proj[217],
        "repair_row218_pressure_delta_pct":100*(p_repaired[217]/p[217]-1),
        "projected_row218_pressure_delta_pct":100*(p_proj[217]/p[217]-1),
        "repair_row218_pressure_vs_projected_pct":100*(p_repaired[217]/p_proj[217]-1),
        "status":"PASS_DIAGNOSTIC_ONLY_NO_TABLE_WRITTEN",
    }


def write_csv(path,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r)) if rows else ["status"]
    with path.open("w",newline="",encoding="utf-8") as stream:
        w=csv.DictWriter(stream,fieldnames=fields,extrasaction="ignore",lineterminator="\n")
        w.writeheader();w.writerows(rows)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,default=ROOT/"dist/r16-minimal-hermite-repair")
    args=parser.parse_args(argv)
    if sha256(SUMMARY)!=SUMMARY_SHA256: raise RuntimeError("frozen summary hash mismatch")
    with SUMMARY.open(newline="",encoding="utf-8") as stream: refs=list(csv.DictReader(stream))
    if len(refs)!=EXPECTED_TABLES or len({x["name"] for x in refs})!=EXPECTED_TABLES:
        raise RuntimeError("frozen 234 table inventory mismatch")
    records=[];failures=[]
    for ref in refs:
        try: records.append(audit_one(ref))
        except Exception as exc: failures.append({"name":ref.get("name"),"error":f"{type(exc).__name__}: {exc}"[:2000]})

    def aggregate(key):
        rows=records if key is None else [x for x in records if x["core"]==key]
        if not rows:
            return {"tables":0,"cells":0,"summary_unavailable_due_to_table_failures":len(failures)}
        return {
            "tables":len(rows),
            "cells":sum(x["cell_count"] for x in rows),
            "source_hermite_reject_cells":sum(x["source_hermite_reject_cells"] for x in rows),
            "minimal_repair_adjusted_cells":sum(x["minimal_repair_adjusted_cells"] for x in rows),
            "source_pressure_increments_preserved_cells":sum(x["source_pressure_increments_preserved_cells"] for x in rows),
            "max_repair_pressure_delta_over_source_Pmax":max(x["minimal_repair_max_abs_pressure_delta_over_source_Pmax"] for x in rows),
            "max_current_projection_delta_over_source_Pmax":max(x["current_projector_max_abs_pressure_delta_over_source_Pmax"] for x in rows),
            "max_repair_vs_projection_delta_over_source_Pmax":max(x["repair_vs_current_projector_max_abs_pressure_delta_over_source_Pmax"] for x in rows),
            "max_repair_energy_delta_rel_source":max(x["minimal_repair_max_abs_energy_delta_over_source_eps"] for x in rows),
            "max_repaired_min_causal_margin":min(x["repaired_min_causal_margin"] for x in rows),
            "max_repaired_min_d2p":min(x["repaired_min_cell_d2p"] for x in rows),
        }
    summary={
        "all_cohort":aggregate(None),
        "by_core":{core:aggregate(core) for core in sorted({x["core"] for x in records})},
        "by_crust":{
            crust: {
                "table_count":sum(x["crust"]==crust for x in records),
                "source_hermite_reject_cells":sum(x["source_hermite_reject_cells"] for x in records if x["crust"]==crust),
                "minimal_repair_adjusted_cells":sum(x["minimal_repair_adjusted_cells"] for x in records if x["crust"]==crust),
                "max_repair_pressure_delta_over_source_Pmax":max(x["minimal_repair_max_abs_pressure_delta_over_source_Pmax"] for x in records if x["crust"]==crust),
                "max_current_projection_delta_over_source_Pmax":max(x["current_projector_max_abs_pressure_delta_over_source_Pmax"] for x in records if x["crust"]==crust),
                "max_repair_vs_projection_delta_over_source_Pmax":max(x["repair_vs_current_projector_max_abs_pressure_delta_over_source_Pmax"] for x in records if x["crust"]==crust),
            } for crust in sorted({x["crust"] for x in records})
        }
    }
    all_pass=(
        len(records)==EXPECTED_TABLES and not failures and
        all(x["source_hermite_reject_cells"]==x["minimal_repair_adjusted_cells"] for x in records) and
        all(x["source_hermite_pass_cells"]==x["source_pressure_increments_preserved_cells"] for x in records) and
        all(x["status"]=="PASS_DIAGNOSTIC_ONLY_NO_TABLE_WRITTEN" for x in records)
    )
    receipt={
        "schema":"ns-cowling-r16-minimal-hermite-pressure-repair-diagnostic-v1",
        "status":"PASS_R16_MINIMAL_HERMITE_PRESSURE_REPAIR_DIAGNOSTIC" if all_pass else "FAIL_R16_MINIMAL_HERMITE_PRESSURE_REPAIR_DIAGNOSTIC",
        "source_commit":os.environ.get("GITHUB_SHA","UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256":SUMMARY_SHA256,
        "repair_definition":"clip each raw DeltaP to exact Fortran cubic-Hermite admissible interval [L,U] with 1e-10 interior margin; preserve DeltaP exactly where already admitted; integrate repaired increments cumulatively from source P_1",
        "L":"h*(2*n0+n1)/3 + h^2*n0/(6*mu0)",
        "U":"h*(n0+2*n1)/3 - h^2*n1/(6*mu1)",
        "source_tables_modified":False,
        "projected_tables_modified":False,
        "synthetic_eos_tables_written":0,
        "tov_integrations_executed":0,
        "table_count":len(records),
        "failures":failures,
        "summary":summary,
        "tables":records,
    }
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    recpath=out/"r16_minimal_hermite_pressure_repair_receipt.json"
    recpath.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    write_csv(out/"r16_minimal_hermite_pressure_repair_tables.csv",records)
    manifest={
        "schema":"ns-cowling-r16-minimal-hermite-repair-manifest-v1",
        "receipt_sha256":sha256(recpath),
        "tables_csv_sha256":sha256(out/"r16_minimal_hermite_pressure_repair_tables.csv"),
        "source_summary_sha256":SUMMARY_SHA256,
        "table_count":len(records),
        "synthetic_eos_tables_written":0,
        "tov_integrations_executed":0,
    }
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":receipt["status"],"all_cohort":summary["all_cohort"],"failures":len(failures)},indent=2),flush=True)
    return 0 if all_pass else 1


if __name__=="__main__": raise SystemExit(main())
