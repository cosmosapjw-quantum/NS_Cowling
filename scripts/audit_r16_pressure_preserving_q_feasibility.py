#!/usr/bin/env python3
"""Full-cohort feasibility audit for pressure-preserving monotone q=n/mu cells."""
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

from r16_pressure_preserving_q import (  # noqa: E402
    cell_pressure_increment_bounds,
    qexp_pressure_increment,
    qexp_shape_parameter,
)

SUMMARY = ROOT / "oscillation_release_v1/data/spectrum_summary_234.csv"
SOURCE_DIR = ROOT / "oscillation_release_v1/data/eos_original_normalized"
SUMMARY_SHA256 = "ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965"
EXPECTED_TABLES = 234
REGULARIZATION_FRACTION = 1.0e-3
EXACT_MATCH_RELATIVE_TOLERANCE = 5.0e-12
RELATIVE_TOLERANCE = 1.0e-12
METHODS = ("css", "cubic_F", "cubic_P", "maxwell", "quintic_F", "quintic_P")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def numeric(rows: list[dict[str, str]], key: str) -> list[float]:
    values = [float(row[key]) for row in rows]
    if any(not math.isfinite(x) for x in values):
        raise ValueError(f"{key}: nonfinite source values")
    return values


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(rows[0].keys()) if rows else ["status"]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def audit_table(reference: dict[str, str]) -> tuple[dict, list[dict]]:
    name = reference["name"]
    source_path = SOURCE_DIR / f"{name}.csv"
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    before = sha256(source_path)
    expected = str(reference.get("raw_normalized_eos_sha256", ""))
    if not expected or before != expected:
        raise ValueError(f"{name}: normalized source SHA-256 mismatch")

    rows = read_csv(source_path)
    mu = numeric(rows, "mu_b_mev")
    n = numeric(rows, "n_b_fm3")
    pressure = numeric(rows, "pressure_mev_fm3")
    if len(rows) < 3 or not (len(mu) == len(n) == len(pressure)):
        raise ValueError(f"{name}: invalid EOS array lengths")
    if any(int(row["source_row"]) != i for i, row in enumerate(rows, start=1)):
        raise ValueError(f"{name}: source-row identity/order mismatch")
    if any(b <= a for a, b in zip(mu, mu[1:])):
        raise ValueError(f"{name}: chemical potential is not strictly increasing")
    if any(b <= a for a, b in zip(n, n[1:])):
        raise ValueError(f"{name}: baryon density is not strictly increasing")
    if any(b <= a for a, b in zip(pressure, pressure[1:])):
        raise ValueError(f"{name}: source pressure is not strictly increasing")

    q = [ni / mui for ni, mui in zip(n, mu)]
    qbar = []
    cells = []
    previous_qbar = None
    drop_rows = []
    feasible_count = 0
    infeasible_below = 0
    infeasible_above = 0
    constant_q_infeasible = 0
    exact_shape_failures = 0
    regularized_shape_success = 0
    exact_k_values = []
    exact_mismatch_values = []
    exact_chi_values = []
    row217_218_record = None

    for i in range(len(mu) - 1):
        mu0, mu1 = mu[i], mu[i + 1]
        q0, q1 = q[i], q[i + 1]
        dp_source = pressure[i + 1] - pressure[i]
        area = 0.5 * (mu1 * mu1 - mu0 * mu0)
        lower, upper = cell_pressure_increment_bounds(mu0, mu1, q0, q1)
        qbar_i = dp_source / area
        qbar.append(qbar_i)

        span = upper - lower
        tolerance = RELATIVE_TOLERANCE * max(span, abs(dp_source), 1.0e-300)
        if q1 == q0:
            alpha = None
        else:
            alpha = (qbar_i - q0) / (q1 - q0)
        exact_result = qexp_shape_parameter(mu0, mu1, q0, q1, dp_source)
        exact_feasible = bool(exact_result.get("feasible", False))
        if exact_feasible:
            feasible_count += 1
            k_exact = float(exact_result["k"])
            reconstructed_dp = qexp_pressure_increment(mu0, mu1, q0, q1, k_exact)
            mismatch = abs(reconstructed_dp - dp_source) / max(abs(dp_source), 1.0e-300)
            exact_mismatch_values.append(mismatch)
            exact_k_values.append(k_exact)
            exact_chi_values.extend([
                float(exact_result["chi_lower"]),
                float(exact_result["chi_upper"]),
            ])
            regularized = None
            regularized_status = "NOT_NEEDED"
        else:
            k_exact = None
            mismatch = None
            exact_shape_failures += 1
            if exact_result["status"] == "BELOW_MONOTONE_Q_ENVELOPE":
                infeasible_below += 1
            elif exact_result["status"] == "ABOVE_MONOTONE_Q_ENVELOPE":
                infeasible_above += 1
            elif exact_result["status"] == "INFEASIBLE_CONSTANT_Q":
                constant_q_infeasible += 1
            regularized = None
            regularized_status = "NOT_ATTEMPTED"
            if q1 > q0:
                try:
                    regularized = qexp_shape_parameter(
                        mu0, mu1, q0, q1, dp_source,
                        regularization_fraction=REGULARIZATION_FRACTION,
                    )
                    if regularized.get("k") is not None:
                        regularized_shape_success += 1
                        regularized_status = "FINITE_SHAPE_AT_REGULARIZED_INCREMENT"
                    else:
                        regularized_status = str(regularized.get("status"))
                except Exception as exc:
                    regularized_status = f"FAILED:{type(exc).__name__}:{exc}"[:300]

        if previous_qbar is not None:
            if qbar_i < previous_qbar * (1.0 - RELATIVE_TOLERANCE):
                drop_rows.append(i + 1)
        previous_qbar = qbar_i

        cell = {
            "name": name,
            "crust": reference["crust"],
            "core": reference["core"],
            "matching_method": reference["matching_method"],
            "lower_source_row": i + 1,
            "upper_source_row": i + 2,
            "mu_lower_mev": mu0,
            "mu_upper_mev": mu1,
            "n_lower_fm3": n[i],
            "n_upper_fm3": n[i + 1],
            "q_lower_fm3_per_mev": q0,
            "q_upper_fm3_per_mev": q1,
            "source_delta_p_mev_fm3": dp_source,
            "monotone_q_delta_p_lower_bound": lower,
            "monotone_q_delta_p_upper_bound": upper,
            "source_qbar": qbar_i,
            "envelope_fraction_alpha": alpha,
            "exact_pressure_preserving_feasible": exact_feasible,
            "exact_shape_status": exact_result["status"],
            "exact_qexp_k": k_exact,
            "exact_reconstruction_relative_error": mismatch,
            "exact_chi_lower": exact_result.get("chi_lower"),
            "exact_chi_upper": exact_result.get("chi_upper"),
            "regularization_fraction": REGULARIZATION_FRACTION if not exact_feasible else None,
            "regularized_shape_status": regularized_status,
            "regularized_qexp_k": regularized.get("k") if regularized else None,
            "regularized_target_delta_p": regularized.get("target_pressure_increment") if regularized else None,
            "regularized_delta_p_change_relative": (
                (regularized["target_pressure_increment"] / dp_source - 1.0)
                if regularized and dp_source > 0.0 else None
            ),
        }
        cells.append(cell)
        if i == 216:
            row217_218_record = cell

    after = sha256(source_path)
    if before != after:
        raise RuntimeError(f"{name}: source file changed during read-only audit")

    qbar_drop_count = len(drop_rows)
    if not row217_218_record:
        raise RuntimeError(f"{name}: source rows 217 and 218 are not available")
    table = {
        "name": name,
        "crust": reference["crust"],
        "core": reference["core"],
        "matching_method": reference["matching_method"],
        "row_count": len(rows),
        "source_sha256_expected": expected,
        "source_sha256_before": before,
        "source_sha256_after": after,
        "source_unchanged": before == after,
        "native_cell_count": len(cells),
        "exact_pressure_preserving_cell_count": feasible_count,
        "exact_pressure_preserving_infeasible_cell_count": exact_shape_failures,
        "infeasible_below_envelope_count": infeasible_below,
        "infeasible_above_envelope_count": infeasible_above,
        "infeasible_constant_q_count": constant_q_infeasible,
        "regularized_finite_shape_success_count": regularized_shape_success,
        "qbar_decrease_count": qbar_drop_count,
        "qbar_decrease_lower_cell_indices": ";".join(str(x) for x in drop_rows),
        "qbar_monotone": qbar_drop_count == 0,
        "max_exact_pressure_reconstruction_relative_error": max(exact_mismatch_values, default=None),
        "exact_qexp_k_min": min(exact_k_values, default=None),
        "exact_qexp_k_max": max(exact_k_values, default=None),
        "exact_qexp_chi_endpoint_min": min(exact_chi_values, default=None),
        "exact_qexp_chi_endpoint_max": max(exact_chi_values, default=None),
        "row217_218_exact_feasible": row217_218_record["exact_pressure_preserving_feasible"],
        "row217_218_qexp_k": row217_218_record["exact_qexp_k"],
        "row217_218_alpha": row217_218_record["envelope_fraction_alpha"],
        "row217_218_chi_lower": row217_218_record["exact_chi_lower"],
        "row217_218_chi_upper": row217_218_record["exact_chi_upper"],
        "row217_218_delta_p_raw": row217_218_record["source_delta_p_mev_fm3"],
        "row217_218_delta_p_lower_bound": row217_218_record["monotone_q_delta_p_lower_bound"],
        "row217_218_delta_p_upper_bound": row217_218_record["monotone_q_delta_p_upper_bound"],
        "status": "PASS_SOURCE_INTEGRITY_AND_CELL_FEASIBILITY_AUDIT",
    }
    return table, cells


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path,
        default=ROOT / "dist/r16-pressure-preserving-q-feasibility",
    )
    args = parser.parse_args(argv)

    if sha256(SUMMARY) != SUMMARY_SHA256:
        raise RuntimeError("frozen spectrum summary SHA-256 mismatch")
    references = read_csv(SUMMARY)
    names = {x["name"] for x in references}
    if len(references) != EXPECTED_TABLES or len(names) != EXPECTED_TABLES:
        raise RuntimeError("frozen spectrum summary is not exactly 234 unique tables")
    if {x["matching_method"] for x in references} != set(METHODS):
        raise RuntimeError("matching method inventory changed")
    if {p.stem for p in SOURCE_DIR.glob("*.csv")} != names:
        raise RuntimeError("normalized source directory inventory differs from frozen cohort")

    table_rows, cell_rows, failures = [], [], []
    by_name = {x["name"]: x for x in references}
    for name in sorted(names):
        try:
            table, cells = audit_table(by_name[name])
            table_rows.append(table)
            cell_rows.extend(cells)
        except Exception as exc:
            failures.append({"name": name, "error": f"{type(exc).__name__}: {exc}"[:2000]})

    expected_cells = sum(x["row_count"] - 1 for x in table_rows)
    exact_feasible = sum(x["exact_pressure_preserving_cell_count"] for x in table_rows)
    exact_infeasible = sum(x["exact_pressure_preserving_infeasible_cell_count"] for x in table_rows)
    total_qbar_drops = sum(x["qbar_decrease_count"] for x in table_rows)
    tables_with_drops = sum(x["qbar_decrease_count"] > 0 for x in table_rows)
    boundary_feasible = sum(bool(x["row217_218_exact_feasible"]) for x in table_rows)
    all_sources_unchanged = (
        len(table_rows) == EXPECTED_TABLES
        and all(x["source_unchanged"] for x in table_rows)
    )
    coverage_ok = (
        len(cell_rows) == expected_cells
        and exact_feasible + exact_infeasible == len(cell_rows)
    )
    qbar_obstruction = total_qbar_drops > 0 and tables_with_drops == EXPECTED_TABLES
    max_exact_error = max(
        (x["max_exact_pressure_reconstruction_relative_error"] for x in table_rows),
        default=0.0,
    )
    exact_reconstruction_ok = max_exact_error <= EXACT_MATCH_RELATIVE_TOLERANCE
    overall_pass = (
        not failures and coverage_ok and all_sources_unchanged
        and boundary_feasible == EXPECTED_TABLES and qbar_obstruction
        and exact_reconstruction_ok
    )

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "r16_pressure_preserving_q_cells.csv", cell_rows)
    write_csv(out / "r16_pressure_preserving_q_tables.csv", table_rows)

    crust_groups: dict[str, list[dict]] = defaultdict(list)
    for row in table_rows:
        crust_groups[row["crust"]].append(row)
    crust_summary = []
    for crust, rows in sorted(crust_groups.items()):
        boundary = [x for x in rows if x["row217_218_exact_feasible"]]
        crust_summary.append({
            "crust": crust,
            "table_count": len(rows),
            "exact_pressure_preserving_cell_count": sum(x["exact_pressure_preserving_cell_count"] for x in rows),
            "exact_pressure_preserving_infeasible_cell_count": sum(x["exact_pressure_preserving_infeasible_cell_count"] for x in rows),
            "qbar_decrease_count": sum(x["qbar_decrease_count"] for x in rows),
            "row217_218_exact_feasible_count": len(boundary),
            "row217_218_qexp_k_min": min(x["row217_218_qexp_k"] for x in boundary),
            "row217_218_qexp_k_max": max(x["row217_218_qexp_k"] for x in boundary),
            "row217_218_alpha_min": min(x["row217_218_alpha"] for x in boundary),
            "row217_218_alpha_max": max(x["row217_218_alpha"] for x in boundary),
            "row217_218_chi_endpoint_min": min(
                min(x["row217_218_chi_lower"], x["row217_218_chi_upper"]) for x in boundary
            ),
            "row217_218_chi_endpoint_max": max(
                max(x["row217_218_chi_lower"], x["row217_218_chi_upper"]) for x in boundary
            ),
        })
    write_csv(out / "r16_pressure_preserving_q_crust_summary.csv", crust_summary)

    receipt = {
        "schema": "ns-cowling-r16-pressure-preserving-q-feasibility-v1",
        "status": "PASS_R16_PRESSURE_PRESERVING_Q_FEASIBILITY_AUDIT" if overall_pass else "FAIL_R16_PRESSURE_PRESERVING_Q_FEASIBILITY_AUDIT",
        "source_commit": os.environ.get("GITHUB_SHA", "UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256": SUMMARY_SHA256,
        "source_table_count_expected": EXPECTED_TABLES,
        "source_table_count_audited": len(table_rows),
        "native_cell_count_expected": expected_cells,
        "native_cell_count_audited": len(cell_rows),
        "exact_pressure_preserving_feasible_cell_count": exact_feasible,
        "exact_pressure_preserving_infeasible_cell_count": exact_infeasible,
        "exact_feasible_fraction": exact_feasible / len(cell_rows) if cell_rows else None,
        "infeasible_below_envelope_count": sum(x["infeasible_below_envelope_count"] for x in table_rows),
        "infeasible_above_envelope_count": sum(x["infeasible_above_envelope_count"] for x in table_rows),
        "infeasible_constant_q_count": sum(x["infeasible_constant_q_count"] for x in table_rows),
        "regularized_finite_shape_success_count": sum(x["regularized_finite_shape_success_count"] for x in table_rows),
        "qbar_decrease_total": total_qbar_drops,
        "tables_with_qbar_decrease": tables_with_drops,
        "qbar_sequence_global_monotone": total_qbar_drops == 0,
        "global_exact_pressure_preserving_monotone_q_possible": not qbar_obstruction,
        "row217_218_exact_pressure_preserving_feasible_table_count": boundary_feasible,
        "row217_218_qexp_k_range": [
            min(x["row217_218_qexp_k"] for x in table_rows),
            max(x["row217_218_qexp_k"] for x in table_rows),
        ] if table_rows else None,
        "row217_218_alpha_range": [
            min(x["row217_218_alpha"] for x in table_rows),
            max(x["row217_218_alpha"] for x in table_rows),
        ] if table_rows else None,
        "max_exact_pressure_reconstruction_relative_error": max(
            (x["max_exact_pressure_reconstruction_relative_error"] for x in table_rows),
            default=None,
        ),
        "source_tables_unchanged": all_sources_unchanged,
        "regularization_fraction_for_infeasible_cells": REGULARIZATION_FRACTION,
"exact_pressure_match_relative_tolerance": EXACT_MATCH_RELATIVE_TOLERANCE,
"exact_pressure_reconstruction_within_tolerance": exact_reconstruction_ok,
        "tov_integrations_executed": 0,
        "synthetic_eos_tables_written": 0,
        "failures": failures,
    }
    receipt_path = out / "r16_pressure_preserving_q_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "schema": "ns-cowling-r16-pressure-preserving-q-manifest-v1",
        "receipt_sha256": sha256(receipt_path),
        "cells_csv_sha256": sha256(out / "r16_pressure_preserving_q_cells.csv"),
        "tables_csv_sha256": sha256(out / "r16_pressure_preserving_q_tables.csv"),
        "crust_summary_csv_sha256": sha256(out / "r16_pressure_preserving_q_crust_summary.csv"),
        "source_summary_sha256": SUMMARY_SHA256,
        "audited_table_count": len(table_rows),
        "audited_cell_count": len(cell_rows),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2), flush=True)
    return 0 if overall_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
