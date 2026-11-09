#!/usr/bin/env python3
"""EOS-preserving TOV diagnostics at native finite-pressure table knots."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_cowling.native_table_cutoff_scan import (  # noqa: E402
    load_native_eos_rows,
    select_native_pressure_knots,
)
from compose_cowling.solver import solve_v5_tov_star  # noqa: E402

EOS_DIR = ROOT / "oscillation_release_v1/data/eos_projected"
SUMMARY = ROOT / "oscillation_release_v1/data/spectrum_summary_234.csv"
METHODS = ("css", "cubic_F", "cubic_P", "maxwell", "quintic_F", "quintic_P")
OFFSETS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0)
SUMMARY_SHA256 = "ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965"
BLOB_SHA1 = {
    "css": "5826c044e18054a8709dacde4bc648559cfb5bb9",
    "cubic_F": "d565b317d63b2b1b49876025d0c8c12520f71f7a",
    "cubic_P": "1269eeb036b5bc1a87a1b0ba0aa150eab644d798",
    "maxwell": "d02aeb26c822a138a34647360644d4daef5b031c",
    "quintic_F": "0f224fb2bfdbe7a81afefedb708464a9209b22e8",
    "quintic_P": "a5cb7b7e60895ee96c0183406c82034d1acf7d7c",
}
SOLVER = {"r0_m": 10.0, "r_limit_m": 1.0e8, "rtol": 1.0e-9, "step_cap": 0.05}
SOLAR_MASS_KG = 1.988409870698051e30


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def blob_sha1(path: Path) -> str:
    b = path.read_bytes()
    return hashlib.sha1(f"blob {len(b)}\0".encode() + b).hexdigest()


def write_csv(path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=ROOT / "dist/native-table-cutoff-scan")
    ap.add_argument("--baseline-relative-tolerance", type=float, default=2.0e-5)
    args = ap.parse_args(argv)
    tol = args.baseline_relative_tolerance
    if not math.isfinite(tol) or tol <= 0:
        raise SystemExit("baseline-relative-tolerance must be finite and positive")
    if digest(SUMMARY) != SUMMARY_SHA256:
        raise RuntimeError("frozen six-method reference summary SHA-256 mismatch")
    with SUMMARY.open(newline="", encoding="utf-8") as f:
        reference_rows = list(csv.DictReader(f))
    refs = [r for r in reference_rows if r.get("crust") == "GPPVA_FSU2" and r.get("core") == "BBB_BHF2"]
    if len(refs) != 6 or {r.get("matching_method") for r in refs} != set(METHODS):
        raise RuntimeError("frozen GPPVA_FSU2 x BBB_BHF2 six-method reference cell is not intact")
    ref_by_method = {r["matching_method"]: r for r in refs}

    out = args.output_dir.resolve()
    backgrounds = out / "backgrounds"
    backgrounds.mkdir(parents=True, exist_ok=True)
    receipt = {
        "schema": "ns-cowling-native-table-cutoff-scan-v1",
        "status": "RUNNING",
        "source_commit": os.environ.get("GITHUB_SHA", "UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256": SUMMARY_SHA256,
        "cohort": {"crust": "GPPVA_FSU2", "core": "BBB_BHF2", "methods": list(METHODS)},
        "eos_policy": {
            "source_rows_modified": False, "source_rows_reordered": False,
            "source_rows_appended": 0, "synthetic_eos_tables_written": 0,
            "cutoffs_are_exact_native_rows": True,
            "interpolation": "frozen v5 native-knot cubic Hermite P(mu_b)",
            "p_zero_surface_or_physical_crust_claimed": False,
            "interpretation": "finite-pressure endpoints are internal isobars; even the first table row has P>0",
        },
        "solver": {"backend": "frozen v5 calculation/methods/support cell_solve task=star", **SOLVER},
        "target_pressure_offsets_dex": list(OFFSETS),
        "expected_tov_case_count": 48, "completed_tov_case_count": 0,
        "baseline_relative_tolerance": tol,
        "results": [], "table_receipts": [], "failures": [],
    }
    receipt_path = out / "native_table_cutoff_scan.json"
    results: list[dict] = []
    failures: list[dict] = []
    table_receipts: list[dict] = []
    by_method: dict[str, list[dict]] = {}
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    for method in METHODS:
        ref = ref_by_method[method]
        eos = EOS_DIR / f"{ref['name']}.csv"
        before = digest(eos)
        actual_blob = blob_sha1(eos)
        if actual_blob != BLOB_SHA1[method]:
            failures.append({"method": method, "stage": "source_identity", "expected": BLOB_SHA1[method], "actual": actual_blob})
            continue
        try:
            rows = load_native_eos_rows(eos)
            if len(rows) != 323:
                raise ValueError(f"expected 323 native EOS rows, received {len(rows)}")
            knots = select_native_pressure_knots(rows, OFFSETS)
            pc = float(ref["pc_pa"])
            if any(k.pressure_pa >= pc for k in knots):
                raise ValueError("selected native cutoff must be below the frozen central pressure")
        except Exception as exc:
            failures.append({"method": method, "stage": "native_table_intake", "error": f"{type(exc).__name__}: {exc}"[:2000]})
            continue

        table_receipts.append({
            "method": method, "name": eos.name, "relative_path": eos.relative_to(ROOT).as_posix(),
            "bytes": eos.stat().st_size, "row_count": len(rows), "git_blob_sha1": actual_blob,
            "sha256_before": before,
            "selected_native_rows": [{
                "row_index": k.row_index, "source_row": k.source_row,
                "target_offset_dex": k.target_offset_dex, "actual_offset_dex": k.actual_offset_dex,
                "pressure_pa": k.pressure_pa, "energy_mev_fm3": k.energy_mev_fm3,
                "n_b_fm3": k.n_b_fm3, "mu_b_mev": k.mu_b_mev,
            } for k in knots],
        })
        method_rows: list[dict] = []
        by_method[method] = method_rows
        ref_mass, ref_radius = float(ref["mass_kg"]), float(ref["radius_m"])
        for k in knots:
            case = f"{method}__source_row_{k.source_row}"
            profile_path = backgrounds / f"{case}.csv"
            try:
                star = solve_v5_tov_star(eos, pc, p_cut_pa=k.pressure_pa,
                                          background_output_path=profile_path, **SOLVER)
                profile = star.get("background_profile")
                if not isinstance(profile, dict) or not profile_path.is_file():
                    raise RuntimeError("native v5 TOV did not export the requested background profile")
                profile_sha = digest(profile_path)
                if profile_sha != profile.get("sha256") or int(profile.get("row_count", 0)) < 2:
                    raise RuntimeError("native v5 background profile failed SHA/row-count validation")
                if Path(profile.get("path", "")) != profile_path:
                    raise RuntimeError("native v5 profile receipt path mismatch")
                if digest(eos) != before:
                    raise RuntimeError("source EOS bytes changed during TOV integration")
                mass, radius = float(star["mass_kg"]), float(star["radius_m"])
                if not all(math.isfinite(x) and x > 0 for x in (mass, radius)):
                    raise ValueError("native TOV returned nonfinite or nonpositive isobar M/R")
                method_rows.append({
                    "case_id": case, "method": method, "source_row": k.source_row,
                    "source_table_sha256": before, "native_row_index": k.row_index,
                    "target_offset_dex": k.target_offset_dex, "actual_offset_dex": k.actual_offset_dex,
                    "p_cut_pa": k.pressure_pa, "p_cut_mev_fm3": k.pressure_pa / 1.602176634e32,
                    "epsilon_cut_mev_fm3": k.energy_mev_fm3, "n_b_cut_fm3": k.n_b_fm3,
                    "mu_b_cut_mev": k.mu_b_mev, "mass_inside_isobar_kg": mass,
                    "mass_inside_isobar_msun": mass / SOLAR_MASS_KG, "isobar_radius_m": radius,
                    "mass_rel_to_frozen_ref": mass / ref_mass - 1,
                    "radius_rel_to_frozen_ref": radius / ref_radius - 1,
                    "is_native_table_edge": k.row_index == 0, "is_complete_p_zero_star": False,
                    "background_profile_file": profile_path.relative_to(out).as_posix(),
                    "background_profile_sha256": profile_sha,
                    "background_profile_row_count": int(profile["row_count"]),
                    "backend_build_key": star.get("backend_build_key"), "status": "ok",
                })
            except Exception as exc:
                failures.append({"case_id": case, "method": method, "stage": "tov_or_profile", "error": f"{type(exc).__name__}: {exc}"[:2000]})
            receipt["results"] = results + method_rows
            receipt["failures"] = failures
            receipt["completed_tov_case_count"] = len(results) + len(method_rows)
            receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        edge = next((r for r in method_rows if r["is_native_table_edge"]), None)
        for row in method_rows:
            if edge:
                row["mass_rel_to_native_edge"] = row["mass_inside_isobar_kg"] / edge["mass_inside_isobar_kg"] - 1
                row["radius_rel_to_native_edge"] = row["isobar_radius_m"] / edge["isobar_radius_m"] - 1
        if edge:
            edge["baseline_reproduction_ok"] = (
                abs(edge["mass_rel_to_frozen_ref"]) <= tol and abs(edge["radius_rel_to_frozen_ref"]) <= tol
            )
            if not edge["baseline_reproduction_ok"]:
                failures.append({"method": method, "stage": "baseline_reproduction", "case_id": edge["case_id"]})
        results.extend(method_rows)
        receipt["results"] = results
        receipt["failures"] = failures
        receipt["completed_tov_case_count"] = len(results)
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    for table in table_receipts:
        p = ROOT / table["relative_path"]
        table["sha256_after"] = digest(p)
        table["unchanged"] = table["sha256_after"] == table["sha256_before"]
        if not table["unchanged"]:
            failures.append({"method": table["method"], "stage": "source_integrity"})

    sensitivity = [{
        "method": r["method"], "source_row": r["source_row"], "p_cut_pa": r["p_cut_pa"],
        "actual_offset_dex": r["actual_offset_dex"],
        "mass_rel_to_native_edge": r.get("mass_rel_to_native_edge"),
        "radius_rel_to_native_edge": r.get("radius_rel_to_native_edge"),
    } for r in results]
    all_tables_unchanged = len(table_receipts) == 6 and all(t["unchanged"] for t in table_receipts)
    all_edges_pass = all(any(r.get("method") == m and r.get("is_native_table_edge") and r.get("baseline_reproduction_ok") for r in results) for m in METHODS)
    receipt.update({
        "results": results, "table_receipts": table_receipts, "failures": failures,
        "sensitivity": sensitivity, "completed_tov_case_count": len(results),
        "all_source_tables_unchanged": all_tables_unchanged,
        "all_native_edge_baselines_reproduced": all_edges_pass,
        "background_profile_count": len(list(backgrounds.glob("*.csv"))),
        "status": "PASS_EOS_PRESERVING_NATIVE_TABLE_CUTOFF_SCAN" if len(results) == 48 and not failures and all_tables_unchanged and all_edges_pass else "FAIL_EOS_PRESERVING_NATIVE_TABLE_CUTOFF_SCAN",
    })
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    fields = tuple(results[0].keys()) if results else ("case_id", "method", "status")
    write_csv(out / "native_table_cutoff_scan.csv", results, fields)
    sens_fields = ("method", "source_row", "p_cut_pa", "actual_offset_dex", "mass_rel_to_native_edge", "radius_rel_to_native_edge")
    write_csv(out / "native_table_cutoff_sensitivity.csv", sensitivity, sens_fields)
    profiles = [{"path": r["background_profile_file"], "sha256": r["background_profile_sha256"], "row_count": r["background_profile_row_count"]} for r in results if r.get("background_profile_file")]
    manifest = {
        "schema": "ns-cowling-native-table-cutoff-manifest-v1", "source_commit": receipt["source_commit"],
        "receipt_sha256": digest(receipt_path), "summary_sha256": SUMMARY_SHA256,
        "summary_csv_sha256": digest(out / "native_table_cutoff_scan.csv"),
        "sensitivity_csv_sha256": digest(out / "native_table_cutoff_sensitivity.csv"),
        "source_table_count": len(table_receipts), "source_tables_unchanged": all_tables_unchanged,
        "synthetic_eos_table_count": 0, "background_profile_count": len(profiles),
        "background_profiles": profiles, "expected_tov_case_count": 48, "completed_tov_case_count": len(results),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": receipt["status"], "completed_tov_case_count": len(results),
                      "source_tables_unchanged": all_tables_unchanged, "native_edge_baselines_reproduced": all_edges_pass,
                      "background_profile_count": len(profiles), "synthetic_eos_table_count": 0,
                      "failure_count": len(failures)}, indent=2))
    return 0 if receipt["status"].startswith("PASS_") else 1


if __name__ == "__main__":
    raise SystemExit(main())