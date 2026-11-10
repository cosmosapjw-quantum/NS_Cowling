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
REPAIR_MARGIN_FRACTION=1.0e-8
REQUIRED=("source_row","mu_b_mev","pressure_mev_fm3","energy_mev_fm3","n_b_fm3","mu_b_j","pressure_pa","n_b_m3")


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
    if len(sr)!=len(pr):
        raise ValueError(f"{name}: row count mismatch")

    # SI columns are authoritative because these are the values passed to
    # Fortran construct_bridge_endpoints. Repair is performed in SI directly.
    mu_si=[float(x["mu_b_j"]) for x in sr]
    n_si=[float(x["n_b_m3"]) for x in sr]
    p_source_si=[float(x["pressure_pa"]) for x in sr]
    eps_source_si=[float(x["energy_j_m3"]) for x in sr]
    p_projected_si=[float(x["pressure_pa"]) for x in pr]
    eps_projected_si=[float(x["energy_j_m3"]) for x in pr]

    pressure_unit_pa=1.602176634e32
    mu_mev=[float(x["mu_b_mev"]) for x in sr]
    n_mev=[float(x["n_b_fm3"]) for x in sr]

    p_repaired_si=[p_source_si[0]]
    source_status_counts={
        "REJECT_NONPOSITIVE_COMPRESSIBILITY":0,
        "REJECT_ACAUSAL_NATIVE_CELL":0,
        "PASS_NATIVE_CELL_ADMISSIBLE":0
    }
    min_d2_si=float("inf")
    min_g_si=float("inf")
    min_d2_mev=float("inf")
    min_g_mev=float("inf")
    max_abs_local_repair=0.0
    sum_abs_increment_repair=0.0
    adjusted=0
    si_pass_cells=0
    si_fail_cells=0
    min_si_margin_fraction=float("inf")
    max_si_margin_fraction=0.0

    for i in range(len(mu_si)-1):
        source_cell_si=audit_hermite_cell(
            mu0=mu_si[i],p0=p_source_si[i],n0=n_si[i],
            mu1=mu_si[i+1],p1=p_source_si[i+1],n1=n_si[i+1]
        )
        source_status_counts[source_cell_si["status"]]+=1

        delta_source_si=p_source_si[i+1]-p_source_si[i]
        repair=minimal_admissible_delta_p(
            mu0=mu_si[i],n0=n_si[i],mu1=mu_si[i+1],n1=n_si[i+1],
            delta_p_source=delta_source_si,
            interior_margin_fraction=REPAIR_MARGIN_FRACTION,
        )
        delta_star_si=repair["delta_p_repaired"]
        adjusted += int(repair["was_adjusted"])
        sum_abs_increment_repair += abs(repair["delta_p_adjustment"])
        max_abs_local_repair=max(
            max_abs_local_repair,
            abs(repair["delta_p_adjustment"])/max(abs(delta_source_si),1e-300)
        )
        width_si=repair["delta_p_width"]
        actual_margin_fraction=min(
            delta_star_si-repair["delta_p_lower"],
            repair["delta_p_upper"]-delta_star_si
        )/width_si
        min_si_margin_fraction=min(min_si_margin_fraction,actual_margin_fraction)
        max_si_margin_fraction=max(max_si_margin_fraction,actual_margin_fraction)
        p_repaired_si.append(p_repaired_si[-1]+delta_star_si)

        cell_si=audit_hermite_cell(
            mu0=mu_si[i],p0=p_repaired_si[i],n0=n_si[i],
            mu1=mu_si[i+1],p1=p_repaired_si[i+1],n1=n_si[i+1]
        )
        if cell_si["status"]=="PASS_NATIVE_CELL_ADMISSIBLE":
            si_pass_cells+=1
        else:
            si_fail_cells+=1
            raise RuntimeError(f"{name}: SI-native repaired cell {i+1} not admitted: {cell_si}")
        min_d2_si=min(min_d2_si,cell_si["min_d2p"])
        min_g_si=min(min_g_si,cell_si["min_causal_margin"])

        # Same repaired pressure, converted for a separate normalized-unit check.
        p0_mev=p_repaired_si[-2]/pressure_unit_pa
        p1_mev=p_repaired_si[-1]/pressure_unit_pa
        cell_mev=audit_hermite_cell(
            mu0=mu_mev[i],p0=p0_mev,n0=n_mev[i],
            mu1=mu_mev[i+1],p1=p1_mev,n1=n_mev[i+1]
        )
        if cell_mev["status"]=="PASS_NATIVE_CELL_ADMISSIBLE":
            min_d2_mev=min(min_d2_mev,cell_mev["min_d2p"])
            min_g_mev=min(min_g_mev,cell_mev["min_causal_margin"])

    eps_repaired_si=[
        m*nn-pp for m,nn,pp in zip(mu_si,n_si,p_repaired_si)
    ]
    if any(x<=0.0 for x in eps_repaired_si):
        raise RuntimeError(f"{name}: SI repaired energy density became nonpositive")
    if any(b<=a for a,b in zip(eps_repaired_si,eps_repaired_si[1:])):
        raise RuntimeError(f"{name}: SI repaired energy density is not monotone")
    if any(b<=a for a,b in zip(p_repaired_si,p_repaired_si[1:])):
        raise RuntimeError(f"{name}: SI repaired pressure is not monotone")

    p_repaired_mev=[x/pressure_unit_pa for x in p_repaired_si]
    eps_repaired_mev=[x/pressure_unit_pa for x in eps_repaired_si]
    pmax_si=max(p_source_si)
    pdelta_source_si=[a-b for a,b in zip(p_repaired_si,p_source_si)]
    pdelta_projected_si=[a-b for a,b in zip(p_repaired_si,p_projected_si)]
    epsdelta_source_si=[a-b for a,b in zip(eps_repaired_si,eps_source_si)]
    epsdelta_projected_si=[a-b for a,b in zip(eps_repaired_si,eps_projected_si)]
    max_local_pressure_source=max(
        abs(d)/max(abs(x),1e-300) for d,x in zip(pdelta_source_si,p_source_si)
    )
    max_local_pressure_projected=max(
        abs(d)/max(abs(x),1e-300) for d,x in zip(pdelta_projected_si,p_projected_si)
    )

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
        "source_hermite_reject_cells":(len(sr)-1)-source_status_counts["PASS_NATIVE_CELL_ADMISSIBLE"],
        "minimal_repair_adjusted_cells":adjusted,
        "source_pressure_increments_preserved_cells":len(sr)-1-adjusted,
        "minimal_repair_max_relative_increment_change":max_abs_local_repair,
        "minimal_repair_sum_abs_increment_change_pa":sum_abs_increment_repair,
        "minimum_actual_SI_interval_margin_fraction":min_si_margin_fraction,
        "maximum_actual_SI_interval_margin_fraction":max_si_margin_fraction,
        "minimal_repair_max_abs_pressure_delta_over_source_Pmax":max(abs(x) for x in pdelta_source_si)/pmax_si,
        "current_projector_max_abs_pressure_delta_over_source_Pmax":max(abs(a-b) for a,b in zip(p_projected_si,p_source_si))/pmax_si,
        "repair_vs_current_projector_max_abs_pressure_delta_over_source_Pmax":max(abs(x) for x in pdelta_projected_si)/pmax_si,
        "repair_vs_current_projector_rms_pressure_delta_over_source_Pmax":math.sqrt(sum(x*x for x in pdelta_projected_si)/len(pdelta_projected_si))/pmax_si,
        "minimal_repair_max_local_pressure_delta_pct":100*max_local_pressure_source,
        "repair_vs_projector_max_local_pressure_delta_pct":100*max_local_pressure_projected,
        "minimal_repair_max_abs_energy_delta_over_source_eps":max(
            abs(x)/max(abs(y),1e-300) for x,y in zip(epsdelta_source_si,eps_source_si)
        ),
        "current_projector_max_abs_energy_delta_over_source_eps":relative_max(eps_projected_si,eps_source_si),
        "repair_vs_projector_max_abs_energy_delta_over_source_eps":max(
            abs(x)/max(abs(y),1e-300) for x,y in zip(epsdelta_projected_si,eps_projected_si)
        ),
        "repaired_min_energy_mev_fm3":min(eps_repaired_mev),
        "repaired_min_delta_energy_mev_fm3":min(
            b-a for a,b in zip(eps_repaired_mev,eps_repaired_mev[1:])
        ),
        "repaired_min_cell_d2p_si":min_d2_si,
        "repaired_min_causal_margin_si":min_g_si,
        "repaired_min_cell_d2p_mev":min_d2_mev if min_d2_mev!=float("inf") else None,
        "repaired_min_causal_margin_mev":min_g_mev if min_g_mev!=float("inf") else None,
        "si_repaired_cell_pass_count":si_pass_cells,
        "si_repaired_cell_fail_count":si_fail_cells,
        "source_row217_pressure":p_source_si[216]/pressure_unit_pa,
        "source_row218_pressure":p_source_si[217]/pressure_unit_pa,
        "repaired_row217_pressure":p_repaired_mev[216],
        "repaired_row218_pressure":p_repaired_mev[217],
        "projected_row217_pressure":p_projected_si[216]/pressure_unit_pa,
        "projected_row218_pressure":p_projected_si[217]/pressure_unit_pa,
        "repair_row218_pressure_delta_pct":100*(p_repaired_si[217]/p_source_si[217]-1),
        "projected_row218_pressure_delta_pct":100*(p_projected_si[217]/p_source_si[217]-1),
        "repair_row218_pressure_vs_projected_pct":100*(p_repaired_si[217]/p_projected_si[217]-1),
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
            "minimum_repaired_causal_margin_si":min(x["repaired_min_causal_margin_si"] for x in rows),
            "minimum_repaired_d2p_si":min(x["repaired_min_cell_d2p_si"] for x in rows),
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
        all(x["status"]=="PASS_DIAGNOSTIC_ONLY_NO_TABLE_WRITTEN" for x in records) and
        all(x["si_repaired_cell_pass_count"]==x["cell_count"] and x["si_repaired_cell_fail_count"]==0 for x in records)
    )
    receipt={
        "schema":"ns-cowling-r16-minimal-hermite-pressure-repair-diagnostic-v1",
        "status":"PASS_R16_MINIMAL_HERMITE_PRESSURE_REPAIR_DIAGNOSTIC" if all_pass else "FAIL_R16_MINIMAL_HERMITE_PRESSURE_REPAIR_DIAGNOSTIC",
        "source_commit":os.environ.get("GITHUB_SHA","UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256":SUMMARY_SHA256,
        "repair_definition":"clip each raw SI DeltaP to the exact Fortran cubic-Hermite admissible interval [L_SI,U_SI] with 1e-8 interval-width interior margin; preserve DeltaP exactly where already admitted; integrate repaired increments cumulatively from source P_1 in Pa",
        "si_conversion": {"MeV_to_J":1.602176634e-13,"MeV_fm3_to_Pa":1.602176634e32,"fm3_inverse_to_m3_inverse":1.0e45},
        "si_unit_validation": "all repaired cells are rechecked using source mu_b_j and n_b_m3 columns and pressure_pa endpoint reconstruction",
        "repair_margin_fraction_of_SI_admissible_interval": REPAIR_MARGIN_FRACTION,
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
