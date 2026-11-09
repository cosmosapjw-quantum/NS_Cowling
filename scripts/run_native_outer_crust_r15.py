#!/usr/bin/env python3
"""R15: native-knot low-density TOV cutoff pilot for three crossed EOS pairs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_cowling.native_table_cutoff_scan import (  # noqa: E402
    load_native_eos_rows,
    select_native_density_knots,
)
from compose_cowling.solver import solve_v5_tov_star  # noqa: E402

EOS_DIR = ROOT / "oscillation_release_v1/data/eos_projected"
SUMMARY = ROOT / "oscillation_release_v1/data/spectrum_summary_234.csv"
SUMMARY_SHA256 = "ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965"

METHODS = ("css", "cubic_F", "cubic_P", "maxwell", "quintic_F", "quintic_P")
PAIRS = (
    ("GPPVA_FSU2", "BBB_BHF2"),  # R14 sensitivity-tail stress pair
    ("SPG_M2", "BBB_BHF2"),     # same core, different crust family
    ("GPPVA_FSU2", "GM_GM1"),   # same crust, different core family
)
CUTOFF_COUNT = 8
MAXIMUM_DENSITY_FM3 = 1.0e-3
EXPECTED_SOURCE_ROWS = (1, 9, 17, 23, 52, 73, 116, 217)
SOLVER = {"r0_m": 10.0, "r_limit_m": 1.0e8, "rtol": 1.0e-9, "step_cap": 0.05}
SOLAR_MASS_KG = 1.988409870698051e30
MEV_FM3_TO_PA = 1.602176634e32


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path: Path, rows: list[dict], fields: tuple[str, ...] | None = None) -> None:
    if fields is None:
        fields = tuple(rows[0].keys()) if rows else ("status",)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def load_reference_rows() -> list[dict[str, str]]:
    if digest(SUMMARY) != SUMMARY_SHA256:
        raise RuntimeError("frozen 234-row spectrum summary SHA-256 mismatch")
    with SUMMARY.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 234 or len({row.get("name") for row in rows}) != 234:
        raise RuntimeError("frozen 234-row spectrum reference cohort is not intact")
    if {row.get("matching_method") for row in rows} != set(METHODS):
        raise RuntimeError("frozen six-method prescription set changed")
    if any(row.get("projected_native_cell_audit") != "PASS" for row in rows):
        raise RuntimeError("one or more projected EOS tables lack the frozen native-cell PASS")
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "dist/native-outer-crust-r15"
    )
    parser.add_argument("--baseline-relative-tolerance", type=float, default=2.0e-5)
    parser.add_argument("--target-count", type=int, default=CUTOFF_COUNT)
    parser.add_argument(
        "--maximum-density-fm3", type=float, default=MAXIMUM_DENSITY_FM3
    )
    args = parser.parse_args(argv)

    tolerance = float(args.baseline_relative_tolerance)
    maximum_density = float(args.maximum_density_fm3)
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise SystemExit("baseline-relative-tolerance must be finite and positive")
    if args.target_count != CUTOFF_COUNT:
        raise SystemExit(f"R15 pilot is frozen to {CUTOFF_COUNT} cutoff targets")
    if maximum_density != MAXIMUM_DENSITY_FM3:
        raise SystemExit(
            f"R15 pilot outer-segment ceiling is frozen to {MAXIMUM_DENSITY_FM3:.8g} fm^-3"
        )

    references = load_reference_rows()
    selected_pairs = []
    for crust, core in PAIRS:
        rows = [
            row for row in references
            if row.get("crust") == crust and row.get("core") == core
        ]
        if len(rows) != len(METHODS):
            raise RuntimeError(
                f"expected six methods for {crust} x {core}, found {len(rows)}"
            )
        if {row.get("matching_method") for row in rows} != set(METHODS):
            raise RuntimeError(f"incomplete method cell: {crust} x {core}")
        selected_pairs.append((crust, core, {r["matching_method"]: r for r in rows}))

    expected_cases = len(PAIRS) * len(METHODS) * CUTOFF_COUNT
    if expected_cases != 144:
        raise RuntimeError(f"internal R15 case-count contract changed: {expected_cases}")

    output = args.output_dir.resolve()
    backgrounds = output / "backgrounds"
    backgrounds.mkdir(parents=True, exist_ok=True)
    receipt_path = output / "native_outer_crust_r15.json"

    receipt = {
        "schema": "ns-cowling-native-outer-crust-r15-v1",
        "status": "RUNNING",
        "source_commit": os.environ.get("GITHUB_SHA", "UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256": SUMMARY_SHA256,
        "cohort": [
            {"crust": crust, "core": core, "methods": list(METHODS)}
            for crust, core in PAIRS
        ],
        "cohort_design": {
            "pair_count": len(PAIRS),
            "method_count_per_pair": len(METHODS),
            "projected_table_count": len(PAIRS) * len(METHODS),
            "cutoff_count_per_table": CUTOFF_COUNT,
            "expected_tov_case_count": expected_cases,
            "selection_contrasts": [
                "GPPVA_FSU2 x BBB_BHF2: R14 sensitivity-tail stress pair",
                "SPG_M2 x BBB_BHF2: crust-family contrast at fixed core",
                "GPPVA_FSU2 x GM_GM1: core-family contrast at fixed crust",
            ],
        },
        "source_policy": {
            "authoritative_solver_inputs": "oscillation_release_v1/data/eos_projected",
            "normalized_source_tables": "read-only provenance only; not solver inputs",
            "projected_native_cell_audit_required": "PASS",
            "source_rows_modified": False,
            "source_rows_reordered": False,
            "source_rows_appended": 0,
            "synthetic_eos_tables_written": 0,
            "cutoff_values_are_exact_native_rows": True,
            "cutoff_target_variable": "baryon number density n_b",
            "cutoff_selection": (
                "nearest existing rows to eight log-spaced n_b targets from the "
                "first source row through the last source row with n_b <= 1e-3 fm^-3"
            ),
            "maximum_density_fm3": maximum_density,
            "expected_selected_source_rows": list(EXPECTED_SOURCE_ROWS),
            "interpretation": (
                "finite-positive-pressure internal isobar endpoints; no P=0 surface, "
                "complete-star radius, crust mass/thickness or mode boundary is inferred"
            ),
        },
        "solver": {
            "backend": "frozen v5 native Fortran cell_solve task=star",
            **SOLVER,
        },
        "expected_tov_case_count": expected_cases,
        "completed_tov_case_count": 0,
        "baseline_relative_tolerance": tolerance,
        "source_table_receipts": [],
        "results": [],
        "pair_cutoff_summary": [],
        "failures": [],
    }
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    results: list[dict] = []
    source_receipts: list[dict] = []
    failures: list[dict] = []

    for crust, core, by_method in selected_pairs:
        for method in METHODS:
            ref = by_method[method]
            table = EOS_DIR / f"{ref['name']}.csv"
            if not table.is_file():
                failures.append({
                    "pair": f"{crust} x {core}", "method": method,
                    "stage": "input_missing", "path": str(table),
                })
                continue

            sha_before = digest(table)
            expected_sha = str(ref.get("projected_eos_sha256", ""))
            if not expected_sha or sha_before != expected_sha:
                failures.append({
                    "name": ref["name"], "stage": "projected_source_identity",
                    "expected_sha256": expected_sha, "actual_sha256": sha_before,
                })
                continue

            try:
                native_rows = load_native_eos_rows(table)
                knots = select_native_density_knots(
                    native_rows,
                    target_count=CUTOFF_COUNT,
                    maximum_density_fm3=maximum_density,
                )
                source_rows = tuple(int(k.source_row) for k in knots)
                if source_rows != EXPECTED_SOURCE_ROWS:
                    raise ValueError(
                        f"native outer-segment mesh changed: {source_rows}; "
                        f"expected {EXPECTED_SOURCE_ROWS}"
                    )
                pc = float(ref["pc_pa"])
                if any(k.pressure_pa >= pc for k in knots):
                    raise ValueError("a selected cutoff pressure is not below central pressure")
            except Exception as exc:
                failures.append({
                    "name": ref["name"], "stage": "native_table_intake",
                    "error": f"{type(exc).__name__}: {exc}"[:2000],
                })
                continue

            table_receipt = {
                "name": ref["name"],
                "crust": crust,
                "core": core,
                "method": method,
                "relative_path": table.relative_to(ROOT).as_posix(),
                "bytes": table.stat().st_size,
                "row_count": len(native_rows),
                "sha256_expected_from_frozen_summary": expected_sha,
                "sha256_before": sha_before,
                "projected_native_cell_audit": ref["projected_native_cell_audit"],
                "pc_pa": pc,
                "selected_native_rows": [
                    {
                        "row_index_zero_based": k.row_index,
                        "source_row_one_based": int(k.source_row),
                        "target_n_b_fm3": k.target_n_b_fm3,
                        "target_log10_n_b_fm3": k.target_log10_n_b_fm3,
                        "actual_offset_dex_from_table_edge": k.actual_offset_dex,
                        "pressure_pa": k.pressure_pa,
                        "energy_mev_fm3": k.energy_mev_fm3,
                        "n_b_fm3": k.n_b_fm3,
                        "mu_b_mev": k.mu_b_mev,
                    }
                    for k in knots
                ],
            }
            method_results: list[dict] = []
            ref_mass = float(ref["mass_kg"])
            ref_radius = float(ref["radius_m"])

            for knot in knots:
                row_number = int(knot.source_row)
                case_id = f"{crust}__{core}__{method}__row_{row_number}"
                profile_path = backgrounds / f"{case_id}.csv"
                try:
                    star = solve_v5_tov_star(
                        table,
                        pc,
                        p_cut_pa=knot.pressure_pa,
                        background_output_path=profile_path,
                        **SOLVER,
                    )
                    profile = star.get("background_profile")
                    if not isinstance(profile, dict) or not profile_path.is_file():
                        raise RuntimeError(
                            "native v5 TOV did not export the requested background profile"
                        )
                    profile_sha = digest(profile_path)
                    profile_rows = int(profile.get("row_count", 0))
                    if (
                        profile_sha != profile.get("sha256")
                        or profile_rows < 2
                    ):
                        raise RuntimeError(
                            "background profile failed SHA-256 or row-count validation"
                        )
                    if Path(str(profile.get("path", ""))) != profile_path:
                        raise RuntimeError("background-profile path receipt mismatch")
                    if digest(table) != sha_before:
                        raise RuntimeError("projected source EOS bytes changed during TOV")
                    mass = float(star["mass_kg"])
                    radius = float(star["radius_m"])
                    if not all(math.isfinite(v) and v > 0 for v in (mass, radius)):
                        raise ValueError("TOV returned nonfinite/nonpositive isobar M/R")
                    method_results.append({
                        "case_id": case_id,
                        "crust": crust,
                        "core": core,
                        "method": method,
                        "source_row": row_number,
                        "native_row_index_zero_based": knot.row_index,
                        "projected_eos_sha256": sha_before,
                        "target_n_b_fm3": knot.target_n_b_fm3,
                        "target_log10_n_b_fm3": knot.target_log10_n_b_fm3,
                        "actual_offset_dex_from_table_edge": knot.actual_offset_dex,
                        "p_cut_pa": knot.pressure_pa,
                        "p_cut_mev_fm3": knot.pressure_pa / MEV_FM3_TO_PA,
                        "epsilon_cut_mev_fm3": knot.energy_mev_fm3,
                        "n_b_cut_fm3": knot.n_b_fm3,
                        "mu_b_cut_mev": knot.mu_b_mev,
                        "mass_inside_isobar_kg": mass,
                        "mass_inside_isobar_msun": mass / SOLAR_MASS_KG,
                        "isobar_radius_m": radius,
                        "mass_rel_to_frozen_133_anchor": mass / ref_mass - 1.0,
                        "radius_rel_to_frozen_133_anchor": radius / ref_radius - 1.0,
                        "is_native_table_edge": knot.row_index == 0,
                        "is_complete_p_zero_star": False,
                        "background_profile_file": profile_path.relative_to(output).as_posix(),
                        "background_profile_sha256": profile_sha,
                        "background_profile_row_count": profile_rows,
                        "backend_build_key": star.get("backend_build_key"),
                        "status": "ok",
                    })
                except Exception as exc:
                    failures.append({
                        "case_id": case_id,
                        "stage": "native_tov_or_background_profile",
                        "error": f"{type(exc).__name__}: {exc}"[:2000],
                    })

                receipt["results"] = results + method_results
                receipt["failures"] = failures
                receipt["completed_tov_case_count"] = len(results) + len(method_results)
                receipt_path.write_text(
                    json.dumps(receipt, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )

            edge = next(
                (row for row in method_results if row["is_native_table_edge"]), None
            )
            if edge is not None:
                for row in method_results:
                    row["mass_rel_to_native_table_edge"] = (
                        row["mass_inside_isobar_kg"] / edge["mass_inside_isobar_kg"] - 1.0
                    )
                    row["radius_rel_to_native_table_edge"] = (
                        row["isobar_radius_m"] / edge["isobar_radius_m"] - 1.0
                    )
                edge["baseline_reproduction_ok"] = (
                    abs(edge["mass_rel_to_frozen_133_anchor"]) <= tolerance
                    and abs(edge["radius_rel_to_frozen_133_anchor"]) <= tolerance
                )
                if not edge["baseline_reproduction_ok"]:
                    failures.append({
                        "name": ref["name"], "stage": "native_edge_baseline",
                        "case_id": edge["case_id"],
                    })
            else:
                failures.append({
                    "name": ref["name"], "stage": "native_edge_baseline",
                    "error": "row-1 native-edge integration is missing",
                })

            results.extend(method_results)
            table_receipt["sha256_after"] = digest(table)
            table_receipt["unchanged"] = table_receipt["sha256_after"] == sha_before
            table_receipt["completed_cutoff_count"] = len(method_results)
            if not table_receipt["unchanged"]:
                failures.append({"name": ref["name"], "stage": "source_integrity"})
            if len(method_results) != CUTOFF_COUNT:
                failures.append({
                    "name": ref["name"], "stage": "case_coverage",
                    "expected": CUTOFF_COUNT, "actual": len(method_results),
                })
            source_receipts.append(table_receipt)

            receipt["results"] = results
            receipt["source_table_receipts"] = source_receipts
            receipt["failures"] = failures
            receipt["completed_tov_case_count"] = len(results)
            receipt_path.write_text(
                json.dumps(receipt, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )

    pair_cutoff_summary: list[dict] = []
    for crust, core in PAIRS:
        for source_row in EXPECTED_SOURCE_ROWS:
            cell = [
                row for row in results
                if row["crust"] == crust and row["core"] == core
                and row["source_row"] == source_row
            ]
            if len(cell) != len(METHODS):
                pair_cutoff_summary.append({
                    "crust": crust, "core": core, "source_row": source_row,
                    "complete_six_method_cell": False, "method_count": len(cell),
                })
                continue
            mass_response = [100.0 * row["mass_rel_to_native_table_edge"] for row in cell]
            radius_response = [100.0 * row["radius_rel_to_native_table_edge"] for row in cell]
            pair_cutoff_summary.append({
                "crust": crust, "core": core, "source_row": source_row,
                "complete_six_method_cell": True, "method_count": len(cell),
                "n_b_cut_fm3": cell[0]["n_b_cut_fm3"],
                "p_cut_mev_fm3": cell[0]["p_cut_mev_fm3"],
                "mean_enclosed_mass_change_pct": sum(mass_response) / len(mass_response),
                "max_abs_enclosed_mass_change_pct": max(abs(x) for x in mass_response),
                "mean_isobar_radius_change_pct": sum(radius_response) / len(radius_response),
                "max_abs_isobar_radius_change_pct": max(abs(x) for x in radius_response),
            })

    all_tables_unchanged = (
        len(source_receipts) == len(PAIRS) * len(METHODS)
        and all(item.get("unchanged") for item in source_receipts)
    )
    all_baselines_pass = all(
        any(
            row["crust"] == crust and row["core"] == core
            and row["method"] == method and row["is_native_table_edge"]
            and row.get("baseline_reproduction_ok") is True
            for row in results
        )
        for crust, core in PAIRS for method in METHODS
    )
    all_pair_cutoffs_complete = (
        len(pair_cutoff_summary) == len(PAIRS) * CUTOFF_COUNT
        and all(row.get("complete_six_method_cell") for row in pair_cutoff_summary)
    )
    profile_count = len(list(backgrounds.glob("*.csv")))
    receipt.update({
        "results": results,
        "source_table_receipts": source_receipts,
        "pair_cutoff_summary": pair_cutoff_summary,
        "failures": failures,
        "completed_tov_case_count": len(results),
        "source_tables_unchanged": all_tables_unchanged,
        "all_native_edge_baselines_reproduced": all_baselines_pass,
        "all_pair_cutoff_cells_complete": all_pair_cutoffs_complete,
        "background_profile_count": profile_count,
        "status": (
            "PASS_R15_NATIVE_OUTER_CRUST_CUTOFF_PILOT"
            if len(results) == expected_cases
            and not failures
            and all_tables_unchanged
            and all_baselines_pass
            and all_pair_cutoffs_complete
            and profile_count == expected_cases
            else "FAIL_R15_NATIVE_OUTER_CRUST_CUTOFF_PILOT"
        ),
    })
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    write_csv(output / "native_outer_crust_r15.csv", results)
    summary_fields = (
        "crust", "core", "source_row", "complete_six_method_cell", "method_count",
        "n_b_cut_fm3", "p_cut_mev_fm3", "mean_enclosed_mass_change_pct",
        "max_abs_enclosed_mass_change_pct", "mean_isobar_radius_change_pct",
        "max_abs_isobar_radius_change_pct",
    )
    write_csv(
        output / "native_outer_crust_r15_pair_cutoff_summary.csv",
        pair_cutoff_summary,
        summary_fields,
    )
    profiles = [
        {
            "path": row["background_profile_file"],
            "sha256": row["background_profile_sha256"],
            "row_count": row["background_profile_row_count"],
        }
        for row in results if row.get("background_profile_file")
    ]
    manifest = {
        "schema": "ns-cowling-native-outer-crust-r15-manifest-v1",
        "source_commit": receipt["source_commit"],
        "receipt_sha256": digest(receipt_path),
        "result_csv_sha256": digest(output / "native_outer_crust_r15.csv"),
        "pair_summary_csv_sha256": digest(
            output / "native_outer_crust_r15_pair_cutoff_summary.csv"
        ),
        "source_summary_sha256": SUMMARY_SHA256,
        "source_table_count": len(source_receipts),
        "source_tables_unchanged": all_tables_unchanged,
        "synthetic_eos_table_count": 0,
        "background_profile_count": len(profiles),
        "expected_tov_case_count": expected_cases,
        "completed_tov_case_count": len(results),
        "background_profiles": profiles,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print(json.dumps({
        "status": receipt["status"],
        "selected_pairs": len(PAIRS),
        "projected_source_tables": len(source_receipts),
        "expected_tov_cases": expected_cases,
        "completed_tov_cases": len(results),
        "source_tables_unchanged": all_tables_unchanged,
        "native_edge_baselines_reproduced": all_baselines_pass,
        "complete_pair_cutoff_cells": sum(
            bool(row.get("complete_six_method_cell")) for row in pair_cutoff_summary
        ),
        "background_profile_count": profile_count,
        "synthetic_eos_tables": 0,
        "failure_count": len(failures),
    }, indent=2), flush=True)
    return 0 if receipt["status"].startswith("PASS_") else 1


if __name__ == "__main__":
    raise SystemExit(main())
