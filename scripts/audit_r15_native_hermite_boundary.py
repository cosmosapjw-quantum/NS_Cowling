#!/usr/bin/env python3
"""Audit the exact native cubic-Hermite P(mu_b) cell from source row 217 to 218."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import defaultdict
from pathlib import Path

from r15_native_hermite_audit import audit_native_interval

ROOT = Path(__file__).resolve().parents[1]
EOS_DIR = ROOT / "oscillation_release_v1/data/eos_projected"
SUMMARY = ROOT / "oscillation_release_v1/data/spectrum_summary_234.csv"
SUMMARY_SHA256 = "ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965"
EXPECTED_TABLES = 234
ROW_LOW = 217
ROW_HIGH = 218
REQUIRED_METHODS = {"css", "cubic_F", "cubic_P", "maxwell", "quintic_F", "quintic_P"}
REQUIRED_CORES = {"BBB_BHF2", "GM_GM1", "RG_SLY4"}
REQUIRED_COLUMNS = (
    "source_row", "pressure_pa", "energy_j_m3", "n_b_m3", "mu_b_j",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_references() -> list[dict[str, str]]:
    if sha256(SUMMARY) != SUMMARY_SHA256:
        raise RuntimeError("frozen 234-row spectrum summary SHA-256 mismatch")
    with SUMMARY.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    names = [row.get("name", "") for row in rows]
    if len(rows) != EXPECTED_TABLES or len(set(names)) != EXPECTED_TABLES:
        raise RuntimeError("frozen projected EOS reference cohort is incomplete")
    if {row.get("matching_method") for row in rows} != REQUIRED_METHODS:
        raise RuntimeError("frozen six-method inventory changed")
    if {row.get("core") for row in rows} != REQUIRED_CORES:
        raise RuntimeError("frozen three-core inventory changed")
    pairs = {(row.get("crust"), row.get("core")) for row in rows}
    if len(pairs) != 39:
        raise RuntimeError(f"expected 39 crust-core pairs, found {len(pairs)}")
    if any(row.get("projected_native_cell_audit") != "PASS" for row in rows):
        raise RuntimeError("one or more projected tables lack the frozen native-cell PASS")
    return rows


def read_boundary_rows(path: Path) -> tuple[list[dict[str, str]], int]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = [column for column in REQUIRED_COLUMNS if column not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"missing normalized columns: {missing}")
        rows = list(reader)
    if len(rows) < ROW_HIGH:
        raise ValueError(f"table has only {len(rows)} rows; row 218 is unavailable")
    low, high = rows[ROW_LOW - 1], rows[ROW_HIGH - 1]
    if int(low["source_row"]) != ROW_LOW or int(high["source_row"]) != ROW_HIGH:
        raise ValueError("source_row identity does not agree with file row order")
    return [low, high], len(rows)


def audit_table(path: Path, reference: dict[str, str]) -> dict:
    expected_sha = str(reference.get("projected_eos_sha256", ""))
    before = sha256(path)
    if not expected_sha or before != expected_sha:
        raise ValueError(f"projected EOS hash mismatch: expected {expected_sha}, got {before}")

    boundary_rows, row_count = read_boundary_rows(path)
    low, high = boundary_rows
    numeric = {}
    for label, row in (("low", low), ("high", high)):
        numeric[label] = {
            key: float(row[key])
            for key in ("pressure_pa", "energy_j_m3", "n_b_m3", "mu_b_j")
        }
        if any(not math.isfinite(value) or value <= 0.0 for value in numeric[label].values()):
            raise ValueError(f"nonpositive/nonfinite endpoint value at source row {ROW_LOW if label == 'low' else ROW_HIGH}")

    result = audit_native_interval(
        mu0_j=numeric["low"]["mu_b_j"],
        p0_pa=numeric["low"]["pressure_pa"],
        n0_m3=numeric["low"]["n_b_m3"],
        mu1_j=numeric["high"]["mu_b_j"],
        p1_pa=numeric["high"]["pressure_pa"],
        n1_m3=numeric["high"]["n_b_m3"],
    )

    epsilon_formula_low = (
        numeric["low"]["mu_b_j"] * numeric["low"]["n_b_m3"]
        - numeric["low"]["pressure_pa"]
    )
    epsilon_formula_high = (
        numeric["high"]["mu_b_j"] * numeric["high"]["n_b_m3"]
        - numeric["high"]["pressure_pa"]
    )
    closure_low = abs(epsilon_formula_low - numeric["low"]["energy_j_m3"]) / max(
        abs(numeric["low"]["energy_j_m3"]), 1.0
    )
    closure_high = abs(epsilon_formula_high - numeric["high"]["energy_j_m3"]) / max(
        abs(numeric["high"]["energy_j_m3"]), 1.0
    )
    after = sha256(path)
    if after != before:
        raise RuntimeError("source EOS changed during read-only Hermite audit")

    row = {
        "name": reference["name"],
        "crust": reference["crust"],
        "core": reference["core"],
        "matching_method": reference["matching_method"],
        "projected_eos_sha256_expected": expected_sha,
        "projected_eos_sha256_before": before,
        "projected_eos_sha256_after": after,
        "source_unchanged": before == after,
        "row_count": row_count,
        "low_source_row": ROW_LOW,
        "high_source_row": ROW_HIGH,
        "mu_low_j": numeric["low"]["mu_b_j"],
        "mu_high_j": numeric["high"]["mu_b_j"],
        "p_low_pa": numeric["low"]["pressure_pa"],
        "p_high_pa": numeric["high"]["pressure_pa"],
        "n_low_m3": numeric["low"]["n_b_m3"],
        "n_high_m3": numeric["high"]["n_b_m3"],
        "epsilon_low_pa": result["epsilon_low_pa"],
        "epsilon_high_pa": result["epsilon_high_pa"],
        "epsilon_closure_rel_low": closure_low,
        "epsilon_closure_rel_high": closure_high,
        **result,
    }
    return row


def summarize(records: list[dict]) -> tuple[dict, list[dict]]:
    groups: dict[str, list[dict]] = defaultdict(list)
    for row in records:
        groups[row["crust"]].append(row)
    family_rows = []
    for crust, rows in sorted(groups.items()):
        family_rows.append({
            "crust": crust,
            "table_count": len(rows),
            "density_ratio_min": min(x["density_ratio"] for x in rows),
            "density_ratio_max": max(x["density_ratio"] for x in rows),
            "pressure_ratio_min": min(x["pressure_ratio"] for x in rows),
            "pressure_ratio_max": max(x["pressure_ratio"] for x in rows),
            "chi_min_over_tables": min(x["chi_min"] for x in rows if x["chi_min"] is not None),
            "chi_max_over_tables": max(x["chi_max"] for x in rows if x["chi_max"] is not None),
            "max_stationary_points_per_cell": max(x["chi_stationary_count"] for x in rows),
            "min_d2p_over_tables_pa_per_j2": min(x["min_d2p_pa_per_j2"] for x in rows),
            "min_causal_margin_over_tables_m3": min(x["min_causal_margin_m3"] for x in rows),
            "max_endpoint_energy_closure_rel": max(
                max(x["epsilon_closure_rel_low"], x["epsilon_closure_rel_high"]) for x in rows
            ),
            "secant_dlnp_dlnn_min": min(x["secant_dlnp_dlnn"] for x in rows if x["secant_dlnp_dlnn"] is not None),
            "secant_dlnp_dlnn_max": max(x["secant_dlnp_dlnn"] for x in rows if x["secant_dlnp_dlnn"] is not None),
            "all_cells_native_pass": all(x["status"] == "PASS_NATIVE_CELL_ADMISSIBLE" for x in rows),
        })

    chi_values = [x["chi_min"] for x in records if x["chi_min"] is not None] + [
        x["chi_max"] for x in records if x["chi_max"] is not None
    ]
    secants = [x["secant_dlnp_dlnn"] for x in records if x["secant_dlnp_dlnn"] is not None]
    summary = {
        "table_count": len(records),
        "pass_cell_count": sum(x["status"] == "PASS_NATIVE_CELL_ADMISSIBLE" for x in records),
        "stable_interval_count": sum(bool(x["stable_interval"]) for x in records),
        "causal_interval_count": sum(bool(x["causal_interval"]) for x in records),
        "interior_chi_stationary_interval_count": sum(x["chi_stationary_count"] > 0 for x in records),
        "total_interior_chi_stationary_points": sum(x["chi_stationary_count"] for x in records),
        "density_ratio_range": [
            min(x["density_ratio"] for x in records),
            max(x["density_ratio"] for x in records),
        ],
        "pressure_ratio_range": [
            min(x["pressure_ratio"] for x in records),
            max(x["pressure_ratio"] for x in records),
        ],
        "chi_global_min": min(chi_values),
        "chi_global_max": max(chi_values),
        "secant_dlnp_dlnn_range": [min(secants), max(secants)],
        "minimum_d2p_pa_per_j2": min(x["min_d2p_pa_per_j2"] for x in records),
        "minimum_causal_margin_m3": min(x["min_causal_margin_m3"] for x in records),
        "maximum_endpoint_energy_closure_rel": max(
            max(x["epsilon_closure_rel_low"], x["epsilon_closure_rel_high"]) for x in records
        ),
        "source_tables_unchanged": all(x["source_unchanged"] for x in records),
        "synthetic_eos_tables_written": 0,
        "tov_integrations_executed": 0,
    }
    return summary, family_rows


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = tuple(rows[0].keys()) if rows else ("status",)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            clean = dict(row)
            if isinstance(clean.get("chi_stationary_t"), list):
                clean["chi_stationary_t"] = json.dumps(clean["chi_stationary_t"], separators=(",", ":"))
            writer.writerow(clean)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "dist/r15-native-hermite-row217-218")
    args = parser.parse_args(argv)

    references = load_references()
    expected_names = {row["name"] for row in references}
    csv_files = sorted(EOS_DIR.glob("*.csv"))
    actual_names = {path.stem for path in csv_files}
    if actual_names != expected_names or len(csv_files) != EXPECTED_TABLES:
        raise RuntimeError("projected EOS directory inventory differs from frozen 234-row summary")

    by_name = {row["name"]: row for row in references}
    records, failures = [], []
    for name in sorted(expected_names):
        try:
            records.append(audit_table(EOS_DIR / f"{name}.csv", by_name[name]))
        except Exception as exc:
            failures.append({"name": name, "error": f"{type(exc).__name__}: {exc}"[:2000]})

    summary, family_rows = summarize(records) if records else ({}, [])
    all_pass = (
        not failures
        and len(records) == EXPECTED_TABLES
        and summary.get("pass_cell_count") == EXPECTED_TABLES
        and summary.get("stable_interval_count") == EXPECTED_TABLES
        and summary.get("causal_interval_count") == EXPECTED_TABLES
        and summary.get("source_tables_unchanged") is True
    )

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    summary.update({
        "schema": "ns-cowling-r15-native-hermite-row217-218-audit-v1",
        "status": "PASS_R15_NATIVE_HERMITE_ROW217_218_AUDIT" if all_pass else "FAIL_R15_NATIVE_HERMITE_ROW217_218_AUDIT",
        "source_commit": os.environ.get("GITHUB_SHA", "UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256": SUMMARY_SHA256,
        "source_table_count_expected": EXPECTED_TABLES,
        "source_row_pair": [ROW_LOW, ROW_HIGH],
        "interpolation_id": "cubic_Hermite_P_of_mu_raw_n_v1",
        "fortran_acceptance_representation": (
            "min(P''(0),P''(1)) > 0 and min(mu*P''-P')(endpoints) >= 0; "
            "P'' is linear and d(mu*P''-P')/dmu=mu*P''' has fixed sign"
        ),
        "failures": failures,
        "tables": records,
        "crust_family_summary": family_rows,
    })

    receipt_path = output / "native_hermite_row217_218_receipt.json"
    receipt_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    table_fields = tuple(records[0].keys()) if records else ("status",)
    write_csv(output / "native_hermite_row217_218_tables.csv", records)
    write_csv(output / "native_hermite_row217_218_crust_summary.csv", family_rows)

    manifest = {
        "schema": "ns-cowling-r15-native-hermite-row217-218-manifest-v1",
        "receipt_sha256": sha256(receipt_path),
        "tables_csv_sha256": sha256(output / "native_hermite_row217_218_tables.csv"),
        "crust_summary_csv_sha256": sha256(output / "native_hermite_row217_218_crust_summary.csv"),
        "source_summary_sha256": SUMMARY_SHA256,
        "audited_projected_table_count": len(records),
        "source_tables_modified": False,
        "synthetic_eos_table_count": 0,
        "tov_integration_count": 0,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": summary["status"],
        "audited_tables": len(records),
        "pass_cells": summary.get("pass_cell_count"),
        "stable_intervals": summary.get("stable_interval_count"),
        "causal_intervals": summary.get("causal_interval_count"),
        "interior_chi_stationary_intervals": summary.get("interior_chi_stationary_interval_count"),
        "chi_global_min": summary.get("chi_global_min"),
        "chi_global_max": summary.get("chi_global_max"),
        "density_ratio_range": summary.get("density_ratio_range"),
        "pressure_ratio_range": summary.get("pressure_ratio_range"),
        "failure_count": len(failures),
    }, indent=2), flush=True)
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
