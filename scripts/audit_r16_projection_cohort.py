#!/usr/bin/env python3
"""Full-cohort audit of the frozen normalized-to-projected EOS transformation."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import defaultdict
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from r16_projection_audit import reconstruct_projection  # noqa: E402

SUMMARY = ROOT / "oscillation_release_v1/data/spectrum_summary_234.csv"
SOURCE_DIR = ROOT / "oscillation_release_v1/data/eos_original_normalized"
PROJECTED_DIR = ROOT / "oscillation_release_v1/data/eos_projected"
SUMMARY_SHA256 = "ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965"
EXPECTED_TABLES = 234
METHODS = ("css", "cubic_F", "cubic_P", "maxwell", "quintic_F", "quintic_P")
REQUIRED_SOURCE_COLUMNS = (
    "source_row", "n_b_fm3", "pressure_mev_fm3", "energy_mev_fm3", "mu_b_mev",
)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = [x for x in REQUIRED_SOURCE_COLUMNS if x not in (reader.fieldnames or ())]
        if missing:
            raise ValueError(f"{path.name}: missing columns {missing}")
        rows = list(reader)
    if len(rows) < 3:
        raise ValueError(f"{path.name}: EOS table has fewer than three rows")
    return rows


def finite_column(rows: list[dict[str, str]], field: str) -> list[float]:
    values = [float(row[field]) for row in rows]
    if any(not math.isfinite(x) for x in values):
        raise ValueError(f"{field} contains nonfinite values")
    return values


def max_relative_mismatch(actual: list[float], expected: list[float]) -> float:
    return max(abs(a - e) / max(abs(e), 1.0e-300) for a, e in zip(actual, expected))


def verify_summary() -> list[dict[str, str]]:
    if digest(SUMMARY) != SUMMARY_SHA256:
        raise RuntimeError("frozen spectrum summary SHA-256 mismatch")
    with SUMMARY.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    names = [row.get("name", "") for row in rows]
    if len(rows) != EXPECTED_TABLES or len(set(names)) != EXPECTED_TABLES:
        raise RuntimeError("frozen 234-row spectrum summary is incomplete or non-unique")
    if {row.get("matching_method") for row in rows} != set(METHODS):
        raise RuntimeError("frozen matching-method inventory changed")
    if any(row.get("projected_native_cell_audit") != "PASS" for row in rows):
        raise RuntimeError("projected native-cell PASS missing from frozen reference cohort")
    return rows


def audit_one(reference: dict[str, str]) -> dict:
    name = reference["name"]
    source_path = SOURCE_DIR / f"{name}.csv"
    projected_path = PROJECTED_DIR / f"{name}.csv"
    expected_source_sha = str(reference.get("raw_normalized_eos_sha256", ""))
    expected_projected_sha = str(reference.get("projected_eos_sha256", ""))
    if not source_path.is_file() or not projected_path.is_file():
        raise FileNotFoundError(f"{name}: missing source or projected CSV")

    source_before = digest(source_path)
    projected_before = digest(projected_path)
    if not expected_source_sha or source_before != expected_source_sha:
        raise ValueError(f"{name}: normalized source hash mismatch")
    if not expected_projected_sha or projected_before != expected_projected_sha:
        raise ValueError(f"{name}: projected table hash mismatch")

    source_rows = read_table(source_path)
    projected_rows = read_table(projected_path)
    if len(source_rows) != len(projected_rows):
        raise ValueError(f"{name}: source/projected row count mismatch")

    for i, (s, p) in enumerate(zip(source_rows, projected_rows), start=1):
        if int(s["source_row"]) != i or int(p["source_row"]) != i:
            raise ValueError(f"{name}: source_row sequence/order mismatch at row {i}")
        for field in ("source_row", "i_t", "i_nb", "i_yq"):
            if field in s and field in p and s[field] != p[field]:
                raise ValueError(f"{name}: identity column {field} changed at row {i}")

    mu = finite_column(source_rows, "mu_b_mev")
    n_source = finite_column(source_rows, "n_b_fm3")
    p_source = finite_column(source_rows, "pressure_mev_fm3")
    eps_source = finite_column(source_rows, "energy_mev_fm3")
    result = reconstruct_projection(mu, n_source, p_source, eps_source)

    n_projected = finite_column(projected_rows, "n_b_fm3")
    p_projected = finite_column(projected_rows, "pressure_mev_fm3")
    eps_projected = finite_column(projected_rows, "energy_mev_fm3")
    mu_projected = finite_column(projected_rows, "mu_b_mev")

    mismatch_n = max_relative_mismatch(n_projected, result["projected_n"])
    mismatch_p = max_relative_mismatch(p_projected, result["projected_p"])
    mismatch_eps = max_relative_mismatch(eps_projected, result["projected_eps"])
    mismatch_mu = max_relative_mismatch(mu_projected, mu)
    max_mismatch = max(mismatch_n, mismatch_p, mismatch_eps, mismatch_mu)
    if max_mismatch > 5.0e-12:
        raise ValueError(
            f"{name}: independent reconstruction differs from frozen projected table "
            f"(max relative mismatch={max_mismatch:.6e})"
        )

    source_after = digest(source_path)
    projected_after = digest(projected_path)
    if source_after != source_before or projected_after != projected_before:
        raise RuntimeError(f"{name}: source/projected file bytes changed during audit")

    pooled_blocks = [b for b in result["pav_blocks"] if int(b["width"]) > 1]
    pooled_ranges = [
        f"{int(b['lo_zero_based']) + 1}-{int(b['hi_exclusive_zero_based'])}"
        for b in pooled_blocks
    ]
    return {
        "name": name,
        "crust": reference["crust"],
        "core": reference["core"],
        "matching_method": reference["matching_method"],
        "row_count": result["row_count"],
        "normalized_source_sha256_expected": expected_source_sha,
        "normalized_source_sha256_actual": source_before,
        "projected_sha256_expected": expected_projected_sha,
        "projected_sha256_actual": projected_before,
        "source_unchanged": source_after == source_before,
        "projected_unchanged": projected_after == projected_before,
        "independent_reconstruction_max_rel_mismatch": max_mismatch,
        "reconstruction_n_max_rel_mismatch": mismatch_n,
        "reconstruction_p_max_rel_mismatch": mismatch_p,
        "reconstruction_eps_max_rel_mismatch": mismatch_eps,
        "reconstruction_mu_max_rel_mismatch": mismatch_mu,
        "source_q_monotone": result["source_q_monotone"],
        "source_q_violation_count": result["source_q_violation_count"],
        "pav_block_count": result["pav_block_count"],
        "pav_pooled_block_count": result["pav_pooled_block_count"],
        "pav_pooled_knots": result["pav_pooled_knots"],
        "pav_knot_reduction": result["pav_knot_reduction"],
        "pav_row217_block_width": result["pav_row217_block_width"],
        "pav_row218_block_width": result["pav_row218_block_width"],
        "pav_pooled_ranges_one_based": ";".join(pooled_ranges),
        "max_abs_delta_n_over_n": result["max_abs_delta_n_over_n"],
        "max_abs_delta_p_over_p": result["max_abs_delta_p_over_p"],
        "max_abs_delta_p_over_pmax": result["max_abs_delta_p_over_pmax"],
        "max_abs_delta_eps_over_eps": result["max_abs_delta_eps_over_eps"],
        "max_abs_delta_eps_over_mu_n": result["max_abs_delta_eps_over_mu_n"],
        "source_first_law_max_rel": result["source_first_law_max_rel"],
        "projected_first_law_max_rel": result["projected_first_law_max_rel"],
        "source_cell_count": result["source_cell_count"],
        "source_nonstable_cell_count": result["source_nonstable_cell_count"],
        "source_noncausal_cell_count": result["source_noncausal_cell_count"],
        "projected_nonstable_cell_count": result["projected_nonstable_cell_count"],
        "projected_noncausal_cell_count": result["projected_noncausal_cell_count"],
        "source_min_d2p_all_cells": result["source_min_d2p_all_cells"],
        "source_min_causal_margin_all_cells": result["source_min_causal_margin_all_cells"],
        "projected_min_d2p_all_cells": result["projected_min_d2p_all_cells"],
        "projected_min_causal_margin_all_cells": result["projected_min_causal_margin_all_cells"],
        "source_row217_218_min_d2p": result["source_row217_218_cell"]["min_d2p"] if result["source_row217_218_cell"] else None,
        "source_row217_218_min_causal_margin": result["source_row217_218_cell"]["min_causal_margin"] if result["source_row217_218_cell"] else None,
        "source_row217_218_stable": result["source_row217_218_cell"]["stable"] if result["source_row217_218_cell"] else None,
        "source_row217_218_causal": result["source_row217_218_cell"]["causal"] if result["source_row217_218_cell"] else None,
        "projected_row217_218_min_d2p": result["projected_row217_218_cell"]["min_d2p"] if result["projected_row217_218_cell"] else None,
        "projected_row217_218_min_causal_margin": result["projected_row217_218_cell"]["min_causal_margin"] if result["projected_row217_218_cell"] else None,
        "projected_row217_218_stable": result["projected_row217_218_cell"]["stable"] if result["projected_row217_218_cell"] else None,
        "projected_row217_218_causal": result["projected_row217_218_cell"]["causal"] if result["projected_row217_218_cell"] else None,
        "row217_n_source": result["row217_n_source"],
        "row218_n_source": result["row218_n_source"],
        "row217_mu_source": result["row217_mu_source"],
        "row218_mu_source": result["row218_mu_source"],
        "row217_p_source": result["row217_p_source"],
        "row218_p_source": result["row218_p_source"],
        "row217_p_projected": result["row217_p_projected"],
        "row218_p_projected": result["row218_p_projected"],
        "row217_eps_source": result["row217_eps_source"],
        "row218_eps_source": result["row218_eps_source"],
        "row217_eps_projected": result["row217_eps_projected"],
        "row218_eps_projected": result["row218_eps_projected"],
        "row218_pressure_endpoint_delta_pct": result["row218_pressure_endpoint_delta_pct"],
        "row218_energy_endpoint_delta_pct": result["row218_energy_endpoint_delta_pct"],
        "row217_218_source_delta_p": result["row217_218_source_delta_p"],
        "row217_218_projected_delta_p": result["row217_218_projected_delta_p"],
        "row217_218_delta_p_increment_delta_pct": result["row217_218_delta_p_increment_delta_pct"],
        "status": "PASS_RECONSTRUCTION_MATCH",
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = tuple(rows[0].keys()) if rows else ("status",)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(records: list[dict]) -> list[dict]:
    groups: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in records:
        groups[("crust", row["crust"], row["core"])].append(row)
        groups[("core", row["core"], row["crust"])].append(row)
        groups[("method", row["matching_method"], row["core"])].append(row)
    summaries = []
    for (kind, key, context), rows in sorted(groups.items()):
        summaries.append({
            "group_kind": kind,
            "group": key,
            "context": context,
            "table_count": len(rows),
            "source_q_nonmonotone_table_count": sum(not x["source_q_monotone"] for x in rows),
            "pav_pooled_table_count": sum(x["pav_pooled_block_count"] > 0 for x in rows),
            "pav_pooled_block_count_total": sum(x["pav_pooled_block_count"] for x in rows),
            "pav_knot_reduction_total": sum(x["pav_knot_reduction"] for x in rows),
            "row218_pressure_delta_pct_min": min(x["row218_pressure_endpoint_delta_pct"] for x in rows),
            "row218_pressure_delta_pct_max": max(x["row218_pressure_endpoint_delta_pct"] for x in rows),
            "interval_delta_p_increment_pct_min": min(x["row217_218_delta_p_increment_delta_pct"] for x in rows),
            "interval_delta_p_increment_pct_max": max(x["row217_218_delta_p_increment_delta_pct"] for x in rows),
            "max_reconstruction_rel_mismatch": max(x["independent_reconstruction_max_rel_mismatch"] for x in rows),
            "max_source_first_law_rel": max(x["source_first_law_max_rel"] for x in rows),
            "source_row217_218_rejected_cell_count": sum(not x["source_row217_218_stable"] or not x["source_row217_218_causal"] for x in rows),
            "projected_row217_218_rejected_cell_count": sum(not x["projected_row217_218_stable"] or not x["projected_row217_218_causal"] for x in rows),
        })
    return summaries


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "dist/r16-full-cohort-projection-audit")
    args = parser.parse_args(argv)

    if digest(SUMMARY) != SUMMARY_SHA256:
        raise RuntimeError("frozen 234-row spectrum summary SHA-256 mismatch")
    with SUMMARY.open(newline="", encoding="utf-8") as stream:
        references = list(csv.DictReader(stream))
    names = {row.get("name", "") for row in references}
    if len(references) != EXPECTED_TABLES or len(names) != EXPECTED_TABLES:
        raise RuntimeError("frozen 234-row spectrum summary incomplete or non-unique")
    if {row.get("matching_method") for row in references} != set(METHODS):
        raise RuntimeError("frozen matching method inventory changed")

    source_names = {p.stem for p in SOURCE_DIR.glob("*.csv")}
    projected_names = {p.stem for p in PROJECTED_DIR.glob("*.csv")}
    if source_names != names or projected_names != names:
        raise RuntimeError("normalized/projected directory inventory differs from summary")

    by_name = {row["name"]: row for row in references}
    records, failures = [], []
    for name in sorted(names):
        try:
            records.append(audit_one(by_name[name]))
        except Exception as exc:
            failures.append({"name": name, "error": f"{type(exc).__name__}: {exc}"[:2000]})

    source_unchanged = all(x.get("source_unchanged") for x in records)
    projected_unchanged = all(x.get("projected_unchanged") for x in records)
    all_match = len(records) == EXPECTED_TABLES and all(x["status"] == "PASS_RECONSTRUCTION_MATCH" for x in records)
    projected_cells_admissible = (
        len(records) == EXPECTED_TABLES
        and all(
            x["projected_nonstable_cell_count"] == 0
            and x["projected_noncausal_cell_count"] == 0
            for x in records
        )
    )
    groups = summarize(records)
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    receipt = {
        "schema": "ns-cowling-r16-normalized-to-projected-audit-v1",
        "status": "PASS_R16_FULL_COHORT_PROJECTION_RECONSTRUCTION" if all_match and projected_cells_admissible and not failures and source_unchanged and projected_unchanged else "FAIL_R16_FULL_COHORT_PROJECTION_RECONSTRUCTION",
        "source_commit": os.environ.get("GITHUB_SHA", "UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256": SUMMARY_SHA256,
        "expected_table_count": EXPECTED_TABLES,
        "audited_table_count": len(records),
        "matching_methods": list(METHODS),
        "projection_definition": "weighted PAV of q=n_b/mu_b on trapezoidal mu weights; n=mu*q; exact interval integral of piecewise-linear q for P; epsilon=mu*n-P",
        "source_tables_unchanged": source_unchanged,
        "projected_tables_unchanged": projected_unchanged,
        "all_projector_reconstructions_match": all_match,
        "all_projected_hermite_cells_admissible": projected_cells_admissible,
        "source_q_nonmonotone_table_count": sum(not x["source_q_monotone"] for x in records),
        "pav_pooled_table_count": sum(x["pav_pooled_block_count"] > 0 for x in records),
        "pav_pooled_block_count_total": sum(x["pav_pooled_block_count"] for x in records),
        "pav_knot_reduction_total": sum(x["pav_knot_reduction"] for x in records),
        "source_row217_218_rejected_cell_count": sum(not x["source_row217_218_stable"] or not x["source_row217_218_causal"] for x in records),
        "projected_row217_218_rejected_cell_count": sum(not x["projected_row217_218_stable"] or not x["projected_row217_218_causal"] for x in records),
        "source_tables_with_any_nonstable_cell": sum(x["source_nonstable_cell_count"] > 0 for x in records),
        "source_tables_with_any_noncausal_cell": sum(x["source_noncausal_cell_count"] > 0 for x in records),
        "projected_tables_with_any_nonstable_cell": sum(x["projected_nonstable_cell_count"] > 0 for x in records),
        "projected_tables_with_any_noncausal_cell": sum(x["projected_noncausal_cell_count"] > 0 for x in records),
        "row218_pressure_endpoint_delta_pct_range": [min(x["row218_pressure_endpoint_delta_pct"] for x in records), max(x["row218_pressure_endpoint_delta_pct"] for x in records)] if records else None,
        "row217_218_delta_p_increment_delta_pct_range": [min(x["row217_218_delta_p_increment_delta_pct"] for x in records), max(x["row217_218_delta_p_increment_delta_pct"] for x in records)] if records else None,
        "max_reconstruction_rel_mismatch": max((x["independent_reconstruction_max_rel_mismatch"] for x in records), default=None),
        "max_source_first_law_rel": max((x["source_first_law_max_rel"] for x in records), default=None),
        "max_projected_first_law_rel": max((x["projected_first_law_max_rel"] for x in records), default=None),
        "tov_integrations_executed": 0,
        "synthetic_eos_tables_written": 0,
        "failures": failures,
        "tables": records,
        "group_summaries": groups,
    }

    receipt_path = out / "r16_projection_audit_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(out / "r16_projection_audit_tables.csv", records)
    write_csv(out / "r16_projection_audit_groups.csv", groups)
    manifest = {
        "schema": "ns-cowling-r16-projection-audit-manifest-v1",
        "receipt_sha256": digest(receipt_path),
        "tables_csv_sha256": digest(out / "r16_projection_audit_tables.csv"),
        "groups_csv_sha256": digest(out / "r16_projection_audit_groups.csv"),
        "source_summary_sha256": SUMMARY_SHA256,
        "audited_table_count": len(records),
        "tov_integrations_executed": 0,
        "synthetic_eos_table_count": 0,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": receipt["status"],
        "audited_tables": len(records),
        "failures": len(failures),
        "source_q_nonmonotone_tables": receipt["source_q_nonmonotone_table_count"],
        "pav_pooled_tables": receipt["pav_pooled_table_count"],
        "pav_knot_reduction": receipt["pav_knot_reduction_total"],
        "row218_pressure_delta_pct_range": receipt["row218_pressure_endpoint_delta_pct_range"],
        "row217_218_delta_p_increment_delta_pct_range": receipt["row217_218_delta_p_increment_delta_pct_range"],
        "max_reconstruction_rel_mismatch": receipt["max_reconstruction_rel_mismatch"],
        "source_row217_218_rejected_cells": receipt["source_row217_218_rejected_cell_count"],
        "projected_row217_218_rejected_cells": receipt["projected_row217_218_rejected_cell_count"],
        "source_tables_with_nonstable_cells": receipt["source_tables_with_any_nonstable_cell"],
        "projected_tables_with_nonstable_cells": receipt["projected_tables_with_any_nonstable_cell"],
        "source_tables_unchanged": source_unchanged,
        "projected_tables_unchanged": projected_unchanged,
    }, indent=2), flush=True)
    return 0 if receipt["status"].startswith("PASS_") else 1


if __name__ == "__main__":
    raise SystemExit(main())
