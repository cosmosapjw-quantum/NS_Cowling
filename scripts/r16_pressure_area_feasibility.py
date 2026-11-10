"""Necessary/sufficient endpoint feasibility test for a causal q=n/mu interpolant."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import Counter
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "oscillation_release_v1/data/spectrum_summary_234.csv"
SOURCE_DIR = ROOT / "oscillation_release_v1/data/eos_original_normalized"
PROJECTED_DIR = ROOT / "oscillation_release_v1/data/eos_projected"
SUMMARY_SHA256 = "ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965"
EXPECTED_TABLES = 234
REQUIRED = ("source_row", "mu_b_mev", "pressure_mev_fm3", "n_b_fm3")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def audit_pressure_area_feasibility(
    *, mu0: float, p0: float, n0: float,
    mu1: float, p1: float, n1: float,
    relative_tolerance: float = 1.0e-12,
) -> dict[str, Any]:
    """Test q0 <= DeltaP / integral(mu dmu) <= q1 for q'=d(n/mu)/dmu >= 0.

    For a causal barotrope, chi=c_s^2/c^2=n/(mu dn/dmu)<=1 is equivalent
    to d(n/mu)/dmu >= 0. Therefore q(mu) is nondecreasing. The pressure
    increment is DeltaP=integral(mu*q(mu) dmu), yielding the exact endpoint
    bounds. If the strict bounds hold, a smooth monotone q interpolant can
    realize the integral; equality is a limiting/degenerate case.
    """
    vals = [float(mu0), float(p0), float(n0), float(mu1), float(p1), float(n1)]
    if any(not math.isfinite(x) or x <= 0.0 for x in vals):
        raise ValueError("endpoint values must be finite and positive")
    mu0, p0, n0, mu1, p1, n1 = vals
    if mu1 <= mu0 or p1 <= p0:
        raise ValueError("endpoint mu and pressure must be strictly increasing")

    q0 = n0 / mu0
    q1 = n1 / mu1
    integral_mu = 0.5 * (mu1 * mu1 - mu0 * mu0)
    qbar = (p1 - p0) / integral_mu
    lower_margin = qbar - q0
    upper_margin = q1 - qbar
    scale = max(abs(q0), abs(q1), abs(qbar), 1.0e-300)
    tol = relative_tolerance * scale

    if q1 < q0 - tol:
        status = "INFEASIBLE_NONMONOTONE_Q_ENDPOINTS"
    elif lower_margin < -tol:
        status = "INFEASIBLE_PRESSURE_INCREMENT_BELOW_Q_BOUNDS"
    elif upper_margin < -tol:
        status = "INFEASIBLE_PRESSURE_INCREMENT_ABOVE_Q_BOUNDS"
    elif abs(q1 - q0) <= tol and abs(qbar - q0) <= tol:
        status = "FEASIBLE_CONSTANT_Q_LIMIT"
    elif lower_margin > tol and upper_margin > tol:
        status = "FEASIBLE_MONOTONE_CAUSAL_INTERPOLANT"
    else:
        status = "FEASIBLE_ONLY_ASYMPTOTIC_ENDPOINT_LIMIT"

    span = q1 - q0
    fraction = (qbar - q0) / span if abs(span) > tol else None
    return {
        "status": status,
        "q_low": q0,
        "q_pressure_average": qbar,
        "q_high": q1,
        "integral_mu_dmu": integral_mu,
        "pressure_increment": p1 - p0,
        "lower_margin": lower_margin,
        "upper_margin": upper_margin,
        "q_average_fraction_of_endpoint_span": fraction,
        "relative_tolerance": relative_tolerance,
    }


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = [x for x in REQUIRED if x not in (reader.fieldnames or ())]
        if missing:
            raise ValueError(f"{path.name}: missing columns {missing}")
        rows = list(reader)
    if len(rows) < 2:
        raise ValueError(f"{path.name}: insufficient EOS rows")
    return rows


def audit_table(path: Path, expected_sha: str, *, name: str, cohort: str) -> tuple[dict, list[dict]]:
    before = sha256(path)
    if before != expected_sha:
        raise ValueError(f"{cohort} SHA-256 mismatch for {name}")
    rows = read_rows(path)
    arr = {
        k: [float(r[k]) for r in rows]
        for k in ("mu_b_mev", "pressure_mev_fm3", "n_b_fm3")
    }
    if any(not math.isfinite(x) for v in arr.values() for x in v):
        raise ValueError(f"{name}: nonfinite table values")
    cells = []
    counts = Counter()
    for i in range(len(rows)-1):
        if int(rows[i]["source_row"]) != i+1 or int(rows[i+1]["source_row"]) != i+2:
            raise ValueError(f"{name}: source_row sequence mismatch at cell {i+1}")
        result = audit_pressure_area_feasibility(
            mu0=arr["mu_b_mev"][i],
            p0=arr["pressure_mev_fm3"][i],
            n0=arr["n_b_fm3"][i],
            mu1=arr["mu_b_mev"][i+1],
            p1=arr["pressure_mev_fm3"][i+1],
            n1=arr["n_b_fm3"][i+1],
        )
        counts[result["status"]] += 1
        cells.append({
            "cohort": cohort,
            "name": name,
            "cell_index_one_based": i+1,
            "source_row_low": i+1,
            "source_row_high": i+2,
            **result,
        })
    after = sha256(path)
    if after != before:
        raise RuntimeError(f"{name}: source bytes changed during read-only audit")
    record = {
        "cohort": cohort,
        "name": name,
        "row_count": len(rows),
        "cell_count": len(cells),
        **{f"count_{k.lower()}": v for k,v in counts.items()},
        "feasible_cell_count": sum(v for k,v in counts.items() if k.startswith("FEASIBLE")),
        "infeasible_cell_count": sum(v for k,v in counts.items() if k.startswith("INFEASIBLE")),
        "strict_feasible_cell_count": counts["FEASIBLE_MONOTONE_CAUSAL_INTERPOLANT"],
        "source_unchanged": before == after,
        "sha256_expected": expected_sha,
        "sha256_before": before,
        "sha256_after": after,
    }
    return record, cells


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(k for r in rows for k in r)) if rows else ["status"]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT/"dist/r16-pressure-area-feasibility")
    args = parser.parse_args(argv)
    if sha256(SUMMARY) != SUMMARY_SHA256:
        raise RuntimeError("frozen cohort summary SHA mismatch")
    with SUMMARY.open(newline="",encoding="utf-8") as stream:
        refs=list(csv.DictReader(stream))
    if len(refs)!=EXPECTED_TABLES or len({r["name"] for r in refs})!=EXPECTED_TABLES:
        raise RuntimeError("frozen cohort summary inventory mismatch")

    tables,cells,failures=[],[],[]
    for ref in refs:
        name=ref["name"]
        for cohort,directory,hash_field in (
            ("normalized_source",SOURCE_DIR,"raw_normalized_eos_sha256"),
            ("projected_solver_input",PROJECTED_DIR,"projected_eos_sha256"),
        ):
            try:
                table,cell_rows=audit_table(
                    directory/f"{name}.csv",ref[hash_field],name=name,cohort=cohort
                )
                table.update({"crust":ref["crust"],"core":ref["core"],"matching_method":ref["matching_method"]})
                tables.append(table)
                for cell in cell_rows:
                    cell.update({"crust":ref["crust"],"core":ref["core"],"matching_method":ref["matching_method"]})
                cells.extend(cell_rows)
            except Exception as exc:
                failures.append({"cohort":cohort,"name":name,"error":f"{type(exc).__name__}: {exc}"[:2000]})

    cohort_summary={}
    for cohort in ("normalized_source","projected_solver_input"):
        rs=[r for r in tables if r["cohort"]==cohort]
        cs=[r for r in cells if r["cohort"]==cohort]
        status_counts=Counter(c["status"] for c in cs)
        cohort_summary[cohort]={
            "tables":len(rs),
            "cells":len(cs),
            "strict_feasible_cells":status_counts["FEASIBLE_MONOTONE_CAUSAL_INTERPOLANT"],
            "asymptotic_limit_cells":status_counts["FEASIBLE_ONLY_ASYMPTOTIC_ENDPOINT_LIMIT"]+status_counts["FEASIBLE_CONSTANT_Q_LIMIT"],
            "infeasible_low_pressure_cells":status_counts["INFEASIBLE_PRESSURE_INCREMENT_BELOW_Q_BOUNDS"],
            "infeasible_high_pressure_cells":status_counts["INFEASIBLE_PRESSURE_INCREMENT_ABOVE_Q_BOUNDS"],
            "infeasible_q_endpoint_cells":status_counts["INFEASIBLE_NONMONOTONE_Q_ENDPOINTS"],
            "all_source_hashes_unchanged":all(x["source_unchanged"] for x in rs),
            "strict_feasible_fraction":status_counts["FEASIBLE_MONOTONE_CAUSAL_INTERPOLANT"]/len(cs) if cs else None,
        }
    receipt={
        "schema":"ns-cowling-r16-pressure-area-feasibility-v1",
        "status":"PASS_R16_PRESSURE_AREA_FEASIBILITY_AUDIT" if len(tables)==2*EXPECTED_TABLES and not failures else "FAIL_R16_PRESSURE_AREA_FEASIBILITY_AUDIT",
        "source_commit":os.environ.get("GITHUB_SHA","UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256":SUMMARY_SHA256,
        "mathematical_test":"q0 <= DeltaP/[0.5*(mu1^2-mu0^2)] <= q1 for q=n/mu nondecreasing",
        "interpretation":"Strictly interior bounds admit a smooth monotone q interpolant preserving both endpoint n and P; outside bounds no causal convex barotrope can preserve all four endpoint values.",
        "cohort_summary":cohort_summary,
        "failures":failures,
        "table_records":tables,
        "cell_records":cells,
        "tov_integrations_executed":0,
        "synthetic_eos_tables_written":0,
    }
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    receipt_path=out/"r16_pressure_area_feasibility_receipt.json"
    receipt_path.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    write_csv(out/"r16_pressure_area_feasibility_tables.csv",tables)
    write_csv(out/"r16_pressure_area_feasibility_cells.csv",cells)
    manifest={
        "schema":"ns-cowling-r16-pressure-area-feasibility-manifest-v1",
        "receipt_sha256":sha256(receipt_path),
        "tables_csv_sha256":sha256(out/"r16_pressure_area_feasibility_tables.csv"),
        "cells_csv_sha256":sha256(out/"r16_pressure_area_feasibility_cells.csv"),
        "source_summary_sha256":SUMMARY_SHA256,
        "table_records":len(tables),
        "cell_records":len(cells),
    }
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":receipt["status"],"cohort_summary":cohort_summary,"failures":len(failures)},indent=2),flush=True)
    return 0 if receipt["status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
