#!/usr/bin/env python3
"""Read-only preflight for the full 234-table R15 native-density cutoff cohort."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_cowling.native_table_cutoff_scan import (  # noqa: E402
    load_native_eos_rows,
    select_native_density_knots,
)

EOS_DIR = ROOT / "oscillation_release_v1/data/eos_projected"
SUMMARY = ROOT / "oscillation_release_v1/data/spectrum_summary_234.csv"
SUMMARY_SHA256 = "ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965"
METHODS = ("css", "cubic_F", "cubic_P", "maxwell", "quintic_F", "quintic_P")
EXPECTED_TABLES = 234
CUTOFF_COUNT = 8
MAXIMUM_DENSITY_FM3 = 1.0e-3
REQUIRED_CRUST_COUNT = 13
REQUIRED_CORE_FAMILIES = {"BBB_BHF2", "GM_GM1", "RG_SLY4"}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_references() -> list[dict[str, str]]:
    if digest(SUMMARY) != SUMMARY_SHA256:
        raise RuntimeError("frozen 234-row spectrum summary SHA-256 mismatch")
    with SUMMARY.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    names = [row.get("name", "") for row in rows]
    if len(rows) != EXPECTED_TABLES or len(set(names)) != EXPECTED_TABLES:
        raise RuntimeError("frozen 234-row reference cohort is incomplete or non-unique")
    if {row.get("matching_method") for row in rows} != set(METHODS):
        raise RuntimeError("frozen matching-method set differs from R15")
    if any(row.get("projected_native_cell_audit") != "PASS" for row in rows):
        raise RuntimeError("one or more projected EOS tables lack native-cell PASS")
    pairs = Counter((row.get("crust"), row.get("core")) for row in rows)
    crusts = {row.get("crust") for row in rows}
    cores = {row.get("core") for row in rows}
    if len(pairs) != 39 or any(count != len(METHODS) for count in pairs.values()):
        raise RuntimeError("reference cohort is not 39 complete six-method pairs")
    if len(crusts) != REQUIRED_CRUST_COUNT or cores != REQUIRED_CORE_FAMILIES:
        raise RuntimeError("reference cohort crust/core family inventory changed")
    for row in rows:
        if row.get("name") != (
            f"{row.get('crust')}__{row.get('core')}__{row.get('matching_method')}"
        ):
            raise RuntimeError(f"reference name fields disagree: {row.get('name')}")
    return rows


def audit_native_table(
    path: Path,
    *,
    name: str,
    expected_sha256: str,
    central_pressure_pa: float,
    maximum_density_fm3: float = MAXIMUM_DENSITY_FM3,
    target_count: int = CUTOFF_COUNT,
) -> tuple[dict, list[dict]]:
    """Validate one projected CSV and select table-local native cutoff rows."""
    before = digest(path)
    if before != expected_sha256:
        raise ValueError(
            f"projected EOS SHA-256 mismatch: expected {expected_sha256}, got {before}"
        )
    rows = load_native_eos_rows(path)
    if "source_row" not in rows[0]:
        raise ValueError("projected EOS lacks source_row provenance column")
    for i, row in enumerate(rows, start=1):
        try:
            source_row = int(row["source_row"])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"source_row is not an integer at file row {i}") from exc
        if source_row != i:
            raise ValueError(
                f"source_row order mismatch at file row {i}: source_row={source_row}"
            )

    pc = float(central_pressure_pa)
    if not math.isfinite(pc) or pc <= 0.0:
        raise ValueError("central pressure must be finite and positive")

    knots = select_native_density_knots(
        rows,
        target_count=target_count,
        maximum_density_fm3=maximum_density_fm3,
    )
    if len(knots) != target_count:
        raise ValueError(f"selector returned {len(knots)} cutoffs, expected {target_count}")
    indices = [k.row_index for k in knots]
    if indices[0] != 0 or indices[-1] != max(
        i for i, row in enumerate(rows)
        if float(row["n_b_fm3"]) <= maximum_density_fm3
    ):
        raise ValueError("selected cutoffs do not include both eligible native-prefix endpoints")
    if len(set(indices)) != target_count or indices != sorted(indices):
        raise ValueError("selected native row indices are not unique and increasing")
    if any(k.pressure_pa >= pc for k in knots):
        raise ValueError("one or more selected cutoff pressures are not below central pressure")

    density = [float(row["n_b_fm3"]) for row in rows]
    pressure = [float(row["pressure_pa"]) for row in rows]
    last_eligible = max(i for i, value in enumerate(density) if value <= maximum_density)
    next_index = last_eligible + 1
    boundary = None
    if next_index < len(rows):
        boundary = {
            "last_eligible_source_row": last_eligible + 1,
            "first_ineligible_source_row": next_index + 1,
            "n_b_last_eligible_fm3": density[last_eligible],
            "n_b_first_ineligible_fm3": density[next_index],
            "density_jump_ratio": density[next_index] / density[last_eligible],
            "pressure_last_eligible_pa": pressure[last_eligible],
            "pressure_first_ineligible_pa": pressure[next_index],
            "pressure_jump_ratio": pressure[next_index] / pressure[last_eligible],
        }

    prefix_log_gaps = [
        math.log10(density[i + 1] / density[i]) for i in range(last_eligible)
    ]
    max_target_error = max(
        abs(math.log10(k.n_b_fm3) - k.target_log10_n_b_fm3) for k in knots
    )
    after = digest(path)
    if after != before:
        raise RuntimeError("projected EOS bytes changed during read-only preflight")

    table_record = {
        "name": name,
        "relative_path": path.relative_to(ROOT).as_posix(),
        "bytes": path.stat().st_size,
        "row_count": len(rows),
        "sha256_expected_from_frozen_summary": expected_sha256,
        "sha256_before": before,
        "sha256_after": after,
        "unchanged": before == after,
        "central_pressure_pa": pc,
        "low_density_prefix_count": last_eligible + 1,
        "last_eligible_source_row": last_eligible + 1,
        "selected_source_rows_one_based": indices_to_source_rows(knots),
        "selected_row_index_signature_zero_based": indices,
        "maximum_abs_target_quantization_error_dex": max_target_error,
        "maximum_adjacent_log10_density_gap_in_prefix_dex": max(prefix_log_gaps),
        "boundary_to_first_ineligible": boundary,
        "selected_native_knots": [
            {
                "source_row_one_based": int(k.source_row),
                "row_index_zero_based": k.row_index,
                "target_n_b_fm3": k.target_n_b_fm3,
                "target_log10_n_b_fm3": k.target_log10_n_b_fm3,
                "actual_offset_dex_from_table_edge": k.actual_offset_dex,
                "pressure_pa": k.pressure_pa,
                "pressure_mev_fm3": k.pressure_pa / 1.602176634e32,
                "energy_mev_fm3": k.energy_mev_fm3,
                "n_b_fm3": k.n_b_fm3,
                "mu_b_mev": k.mu_b_mev,
            }
            for k in knots
        ],
    }

    cutoff_rows = []
    for ordinal, k in enumerate(knots, start=1):
        cutoff_rows.append({
            "name": name,
            "cutoff_ordinal": ordinal,
            "source_row_one_based": int(k.source_row),
            "row_index_zero_based": k.row_index,
            "target_n_b_fm3": f"{k.target_n_b_fm3:.17g}",
            "target_log10_n_b_fm3": f"{k.target_log10_n_b_fm3:.17g}",
            "actual_offset_dex_from_table_edge": f"{k.actual_offset_dex:.17g}",
            "pressure_pa": f"{k.pressure_pa:.17g}",
            "pressure_mev_fm3": f"{k.pressure_pa / 1.602176634e32:.17g}",
            "energy_mev_fm3": f"{k.energy_mev_fm3:.17g}",
            "n_b_fm3": f"{k.n_b_fm3:.17g}",
            "mu_b_mev": f"{k.mu_b_mev:.17g}",
        })
    return table_record, cutoff_rows


def indices_to_source_rows(knots) -> list[int]:
    return [int(k.source_row) for k in knots]


def write_csv(path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path,
        default=ROOT / "dist/native-outer-crust-r15-preflight",
    )
    parser.add_argument("--target-count", type=int, default=CUTOFF_COUNT)
    parser.add_argument(
        "--maximum-density-fm3", type=float, default=MAXIMUM_DENSITY_FM3
    )
    args = parser.parse_args(argv)
    if args.target_count != CUTOFF_COUNT:
        raise SystemExit(f"R15 preflight is frozen to {CUTOFF_COUNT} cutoff targets")
    if args.maximum_density_fm3 != MAXIMUM_DENSITY_FM3:
        raise SystemExit(
            f"R15 preflight ceiling is frozen to {MAXIMUM_DENSITY_FM3:.8g} fm^-3"
        )

    references = load_references()
    expected_names = {row["name"] for row in references}
    csv_paths = sorted(EOS_DIR.glob("*.csv"))
    actual_names = {path.stem for path in csv_paths}
    sidecar_names = {path.name[:-len(".projection.json")] for path in EOS_DIR.glob("*.projection.json")}
    failures: list[dict] = []
    if len(csv_paths) != EXPECTED_TABLES or actual_names != expected_names:
        failures.append({
            "stage": "projected_directory_inventory",
            "expected_csv_count": EXPECTED_TABLES,
            "actual_csv_count": len(csv_paths),
            "missing_names": sorted(expected_names - actual_names),
            "unexpected_names": sorted(actual_names - expected_names),
        })
    if sidecar_names != expected_names:
        failures.append({
            "stage": "projection_sidecar_inventory",
            "expected_count": EXPECTED_TABLES,
            "actual_count": len(sidecar_names),
            "missing_names": sorted(expected_names - sidecar_names),
            "unexpected_names": sorted(sidecar_names - expected_names),
        })

    records: list[dict] = []
    selected_rows: list[dict] = []
    by_name = {row["name"]: row for row in references}
    for name in sorted(expected_names):
        ref = by_name[name]
        path = EOS_DIR / f"{name}.csv"
        if not path.is_file():
            continue
        try:
            record, knots = audit_native_table(
                path,
                name=name,
                expected_sha256=ref["projected_eos_sha256"],
                central_pressure_pa=float(ref["pc_pa"]),
                maximum_density_fm3=args.maximum_density_fm3,
                target_count=args.target_count,
            )
            record.update({
                "crust": ref["crust"],
                "core": ref["core"],
                "matching_method": ref["matching_method"],
                "projected_native_cell_audit": ref["projected_native_cell_audit"],
            })
            records.append(record)
            for knot in knots:
                knot.update({
                    "crust": ref["crust"],
                    "core": ref["core"],
                    "matching_method": ref["matching_method"],
                    "projected_eos_sha256": record["sha256_before"],
                })
                selected_rows.append(knot)
        except Exception as exc:
            failures.append({
                "name": name,
                "stage": "native_table_preflight",
                "error": f"{type(exc).__name__}: {exc}"[:2000],
            })

    signature_counts = Counter(
        tuple(record["selected_row_index_signature_zero_based"])
        for record in records
    )
    prefix_counts = [record["low_density_prefix_count"] for record in records]
    row_counts = [record["row_count"] for record in records]
    boundary_ratios = [
        record["boundary_to_first_ineligible"]["density_jump_ratio"]
        for record in records if record["boundary_to_first_ineligible"] is not None
    ]
    expected_selected_rows = EXPECTED_TABLES * CUTOFF_COUNT
    all_pass = (
        not failures
        and len(records) == EXPECTED_TABLES
        and len(selected_rows) == expected_selected_rows
        and all(record["unchanged"] for record in records)
        and all(record["projected_native_cell_audit"] == "PASS" for record in records)
    )

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    receipt = {
        "schema": "ns-cowling-r15-full-cohort-native-density-preflight-v1",
        "status": (
            "PASS_R15_FULL_COHORT_NATIVE_DENSITY_PREFLIGHT"
            if all_pass else "FAIL_R15_FULL_COHORT_NATIVE_DENSITY_PREFLIGHT"
        ),
        "source_commit": os.environ.get("GITHUB_SHA", "UNAVAILABLE_LOCAL_RUN"),
        "source_summary_path": SUMMARY.relative_to(ROOT).as_posix(),
        "source_summary_sha256": SUMMARY_SHA256,
        "authoritative_eos_directory": EOS_DIR.relative_to(ROOT).as_posix(),
        "expected_projected_table_count": EXPECTED_TABLES,
        "audited_projected_table_count": len(records),
        "expected_selected_native_knot_count": expected_selected_rows,
        "selected_native_knot_count": len(selected_rows),
        "cutoff_count_per_table": CUTOFF_COUNT,
        "maximum_density_fm3": args.maximum_density_fm3,
        "planned_full_cohort_tov_count_not_executed": expected_selected_rows,
        "tov_integrations_executed": 0,
        "synthetic_eos_tables_written": 0,
        "source_tables_modified": False,
        "source_tables_reordered": False,
        "row_index_signature_distribution_zero_based": [
            {"signature": list(signature), "table_count": count}
            for signature, count in sorted(signature_counts.items())
        ],
        "unique_row_index_signature_count": len(signature_counts),
        "row_count_range": (
            [min(row_counts), max(row_counts)] if row_counts else None
        ),
        "low_density_prefix_count_range": (
            [min(prefix_counts), max(prefix_counts)] if prefix_counts else None
        ),
        "boundary_density_jump_ratio_range": (
            [min(boundary_ratios), max(boundary_ratios)]
            if boundary_ratios else None
        ),
        "table_receipts": records,
        "failures": failures,
    }
    receipt_path = output / "native_outer_crust_r15_preflight.json"
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    table_fields = (
        "name", "crust", "core", "matching_method", "relative_path", "bytes",
        "row_count", "sha256_expected_from_frozen_summary", "sha256_before",
        "sha256_after", "unchanged", "central_pressure_pa",
        "low_density_prefix_count", "last_eligible_source_row",
        "selected_source_rows_one_based", "selected_row_index_signature_zero_based",
        "maximum_abs_target_quantization_error_dex",
        "maximum_adjacent_log10_density_gap_in_prefix_dex",
    )
    table_csv = []
    for record in records:
        item = dict(record)
        item["selected_source_rows_one_based"] = json.dumps(
            item["selected_source_rows_one_based"], separators=(",", ":")
        )
        item["selected_row_index_signature_zero_based"] = json.dumps(
            item["selected_row_index_signature_zero_based"], separators=(",", ":")
        )
        table_csv.append(item)
    write_csv(output / "native_outer_crust_r15_preflight_tables.csv", table_csv, table_fields)
    cutoff_fields = (
        "name", "crust", "core", "matching_method", "projected_eos_sha256",
        "cutoff_ordinal", "source_row_one_based", "row_index_zero_based",
        "target_n_b_fm3", "target_log10_n_b_fm3",
        "actual_offset_dex_from_table_edge", "pressure_pa",
        "pressure_mev_fm3", "energy_mev_fm3", "n_b_fm3", "mu_b_mev",
    )
    write_csv(
        output / "native_outer_crust_r15_preflight_selected_knots.csv",
        selected_rows, cutoff_fields,
    )

    manifest = {
        "schema": "ns-cowling-r15-full-cohort-native-density-preflight-manifest-v1",
        "receipt_sha256": digest(receipt_path),
        "tables_csv_sha256": digest(output / "native_outer_crust_r15_preflight_tables.csv"),
        "selected_knots_csv_sha256": digest(
            output / "native_outer_crust_r15_preflight_selected_knots.csv"
        ),
        "source_summary_sha256": SUMMARY_SHA256,
        "audited_table_count": len(records),
        "selected_native_knot_count": len(selected_rows),
        "source_tables_modified": False,
        "tov_integrations_executed": 0,
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print(json.dumps({
        "status": receipt["status"],
        "summary_sha256": SUMMARY_SHA256,
        "projected_tables_expected": EXPECTED_TABLES,
        "projected_tables_audited": len(records),
        "selected_native_knots_expected": expected_selected_rows,
        "selected_native_knots_recorded": len(selected_rows),
        "unique_source_row_index_signatures": len(signature_counts),
        "source_row_index_signatures": [
            {"signature": list(signature), "tables": count}
            for signature, count in sorted(signature_counts.items())
        ],
        "row_count_range": receipt["row_count_range"],
        "low_density_prefix_count_range": receipt["low_density_prefix_count_range"],
        "boundary_density_jump_ratio_range": receipt["boundary_density_jump_ratio_range"],
        "tov_integrations_executed": 0,
        "failure_count": len(failures),
    }, indent=2), flush=True)
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
