#!/usr/bin/env python3
"""Measure six-method TOV sensitivity to thermodynamic P=0 surface continuations."""
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
# Direct script execution puts scripts/ ahead of PYTHONPATH.  This repository
# also contains scripts/compose_cowling.py, which otherwise shadows the
# python_frontend/compose_cowling package during imports.
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_cowling.solver import solve_v5_tov_star
from compose_cowling.surface_continuation import write_extended_projected_csv

SUMMARY = ROOT / "oscillation_release_v1" / "data" / "spectrum_summary_234.csv"
EOS_DIR = ROOT / "oscillation_release_v1" / "data" / "eos_projected"
METHODS = ("css", "cubic_F", "cubic_P", "maxwell", "quintic_F", "quintic_P")
SCENARIOS = (
    ("tangent_xmin_1e-8", "power_law", 0.0, "coarse"),
    ("tangent_xmin_1e-12", "power_law", 0.0, "fine"),
    ("exp_quadratic_alpha_m0p5", "exp_quadratic", -0.5, "fine"),
    ("exp_quadratic_alpha_p0p5", "exp_quadratic", 0.5, "fine"),
    ("poly_quadratic_alpha_m0p5", "poly_quadratic", -0.5, "fine"),
    ("poly_quadratic_alpha_p0p5", "poly_quadratic", 0.5, "fine"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def rel_diff(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(reference), 1.0)


def finite_cutoff_pa(eos_path: Path) -> float:
    with eos_path.open(newline="", encoding="utf-8") as stream:
        row = next(csv.DictReader(stream))
    value = float(row["pressure_pa"])
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"invalid source surface pressure in {eos_path}: {value}")
    return value


def solve_star(eos_path: Path, pc_pa: float, pcut_pa: float) -> dict:
    return solve_v5_tov_star(
        eos_path,
        pc_pa,
        p_cut_pa=pcut_pa,
        r0_m=10.0,
        r_limit_m=1.0e8,
        rtol=1.0e-9,
        step_cap=0.05,
    )


def persist_receipt(path: Path, receipt: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "dist" / "surface-continuation")
    parser.add_argument("--tail-points", type=int, default=512)
    parser.add_argument("--coarse-energy-ratio", type=float, default=1.0e-8)
    parser.add_argument("--fine-energy-ratio", type=float, default=1.0e-12)
    parser.add_argument("--baseline-relative-tolerance", type=float, default=2.0e-5)
    args = parser.parse_args(argv)
    if args.tail_points < 64:
        raise SystemExit("--tail-points must be >= 64")
    if not (0.0 < args.fine_energy_ratio < args.coarse_energy_ratio < 1.0):
        raise SystemExit("require 0 < fine-energy-ratio < coarse-energy-ratio < 1")
    if args.baseline_relative_tolerance <= 0.0 or not math.isfinite(args.baseline_relative_tolerance):
        raise SystemExit("baseline-relative-tolerance must be finite and positive")

    references = read_csv(SUMMARY)
    selected = [
        row for row in references
        if row.get("crust") == "GPPVA_FSU2" and row.get("core") == "BBB_BHF2"
    ]
    if len(selected) != 6 or {row["matching_method"] for row in selected} != set(METHODS):
        raise SystemExit("frozen GPPVA_FSU2 x BBB_BHF2 six-method reference cell is not intact")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = args.output_dir / "extended_tables"
    receipt_path = args.output_dir / "surface_continuation_summary.json"
    summary_rows: list[dict] = []
    table_receipts: list[dict] = []
    failures: list[dict] = []
    by_method: dict[str, dict[str, dict]] = {}
    receipt = {
        "schema": "ns-cowling-thermodynamic-surface-continuation-v1",
        "status": "RUNNING",
        "source_commit": os.environ.get("GITHUB_SHA", "UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256": sha256(SUMMARY),
        "cohort": {"crust": "GPPVA_FSU2", "core": "BBB_BHF2", "methods": list(METHODS)},
        "solver": {
            "backend": "frozen v5 calculation/methods/support cell_solve task=star",
            "r0_m": 10.0,
            "r_limit_m": 1.0e8,
            "rtol": 1.0e-9,
            "step_cap": 0.05,
        },
        "case_design": [
            {"id": "finite_cutoff_reference", "family": "source_table", "energy_ratio": 1.0},
            *[
                {"id": label, "family": family, "alpha": alpha,
                 "min_energy_ratio": args.coarse_energy_ratio if ratio_key == "coarse" else args.fine_energy_ratio}
                for label, family, alpha, ratio_key in SCENARIOS
            ],
        ],
        "results": [],
        "table_receipts": [],
        "failures": [],
        "acceptance": {
            "baseline_relative_tolerance": args.baseline_relative_tolerance,
            "baseline_reproduction_required": True,
            "surface_cutoff_and_family_sensitivity": "reported diagnostic; no physical-model selection threshold predeclared",
        },
        "scope": (
            "TOV-only equilibrium-barotrope continuation study. Continuations are mathematical C1 closures, "
            "not literature-sourced outer-crust EOS. No radial stability or f/p1 mode claim is made."
        ),
    }
    persist_receipt(receipt_path, receipt)

    solar_mass = 1.988409870698051e30
    for reference in selected:
        name = reference["name"]
        method = reference["matching_method"]
        eos_path = EOS_DIR / f"{name}.csv"
        method_cases: dict[str, dict] = {}
        by_method[method] = method_cases
        pc_pa = float(reference["pc_pa"])
        ref_mass = float(reference["mass_kg"])
        ref_radius = float(reference["radius_m"])
        source_sha = sha256(eos_path)
        try:
            source_pcut = finite_cutoff_pa(eos_path)
            base = solve_star(eos_path, pc_pa, source_pcut)
            baseline_mass = float(base["mass_kg"])
            baseline_radius = float(base["radius_m"])
            baseline_mass_ref = baseline_mass / ref_mass - 1.0
            baseline_radius_ref = baseline_radius / ref_radius - 1.0
            baseline_ok = (
                abs(baseline_mass_ref) <= args.baseline_relative_tolerance
                and abs(baseline_radius_ref) <= args.baseline_relative_tolerance
            )
            baseline_record = {
                "name": name,
                "method": method,
                "scenario": "finite_cutoff_reference",
                "status": "ok" if baseline_ok else "baseline_reference_mismatch",
                "pc_pa": pc_pa,
                "surface_pressure_pa": source_pcut,
                "mass_kg": baseline_mass,
                "radius_m": baseline_radius,
                "mass_msun": baseline_mass / solar_mass,
                "mass_relative_to_frozen_reference": baseline_mass_ref,
                "radius_relative_to_frozen_reference": baseline_radius_ref,
                "backend_build_key": base.get("backend_build_key"),
                "source_eos_sha256": source_sha,
            }
            summary_rows.append(baseline_record)
            method_cases["finite_cutoff_reference"] = baseline_record
            receipt["results"].append(baseline_record)
            if not baseline_ok:
                failures.append({"name": name, "scenario": "finite_cutoff_reference", "reason": "frozen baseline reproduction exceeds tolerance"})
        except Exception as exc:
            failure = {"name": name, "scenario": "finite_cutoff_reference", "error": f"{type(exc).__name__}: {exc}"[:2000]}
            failures.append(failure)
            receipt["failures"] = list(failures)
            persist_receipt(receipt_path, receipt)
            continue

        for label, family, alpha, ratio_key in SCENARIOS:
            energy_ratio = args.coarse_energy_ratio if ratio_key == "coarse" else args.fine_energy_ratio
            output_eos = tables_dir / f"{name}__{label}.csv"
            try:
                table_receipt = write_extended_projected_csv(
                    eos_path,
                    output_eos,
                    alpha=alpha,
                    family=family,
                    tail_points=args.tail_points,
                    min_energy_ratio=energy_ratio,
                )
                table_receipt.update({"name": name, "method": method, "scenario": label,
                                      "family": family, "alpha": alpha})
                table_receipts.append(table_receipt)
                extension_pcut = float(table_receipt["surface_pressure_pa"])
                star = solve_star(output_eos, pc_pa, extension_pcut)
                mass = float(star["mass_kg"])
                radius = float(star["radius_m"])
                record = {
                    "name": name,
                    "method": method,
                    "scenario": label,
                    "family": family,
                    "alpha": alpha,
                    "min_energy_ratio": energy_ratio,
                    "tail_points": args.tail_points,
                    "pc_pa": pc_pa,
                    "surface_pressure_pa": extension_pcut,
                    "mass_kg": mass,
                    "radius_m": radius,
                    "mass_msun": mass / solar_mass,
                    "mass_relative_to_finite_cutoff": mass / baseline_mass - 1.0,
                    "radius_relative_to_finite_cutoff": radius / baseline_radius - 1.0,
                    "mass_relative_to_frozen_reference": mass / ref_mass - 1.0,
                    "radius_relative_to_frozen_reference": radius / ref_radius - 1.0,
                    "surface_residual_pa": star.get("surface_residual_pa"),
                    "backend_build_key": star.get("backend_build_key"),
                    "extended_eos_sha256": table_receipt["output_sha256"],
                    "source_eos_sha256": source_sha,
                    "continuation_closure": table_receipt["diagnostics"]["closure"],
                    "continuation_cs2_range": table_receipt["diagnostics"]["tail"],
                }
                summary_rows.append(record)
                method_cases[label] = record
                receipt["results"].append(record)
            except Exception as exc:
                failure = {"name": name, "scenario": label, "error": f"{type(exc).__name__}: {exc}"[:2000]}
                failures.append(failure)
                receipt["failures"] = list(failures)
            receipt["table_receipts"] = list(table_receipts)
            receipt["failures"] = list(failures)
            persist_receipt(receipt_path, receipt)

    # Same-method paired sensitivities isolate the lower pressure cutoff and
    # two analytically distinct C1 curvature families at fixed central pressure.
    sensitivity = []
    for method, cases in by_method.items():
        coarse = cases.get("tangent_xmin_1e-8")
        fine = cases.get("tangent_xmin_1e-12")
        if coarse and fine:
            sensitivity.append({
                "method": method,
                "comparison": "tangent_surface_cutoff_1e-8_vs_1e-12",
                "mass_relative_change_fine_vs_coarse": fine["mass_kg"] / coarse["mass_kg"] - 1.0,
                "radius_relative_change_fine_vs_coarse": fine["radius_m"] / coarse["radius_m"] - 1.0,
            })
        for scenario in ("exp_quadratic_alpha_m0p5", "exp_quadratic_alpha_p0p5",
                         "poly_quadratic_alpha_m0p5", "poly_quadratic_alpha_p0p5"):
            alternative = cases.get(scenario)
            if fine and alternative:
                sensitivity.append({
                    "method": method,
                    "comparison": f"{scenario}_vs_tangent_xmin_1e-12",
                    "mass_relative_change_alternative_vs_tangent": alternative["mass_kg"] / fine["mass_kg"] - 1.0,
                    "radius_relative_change_alternative_vs_tangent": alternative["radius_m"] / fine["radius_m"] - 1.0,
                })
    receipt["sensitivity"] = sensitivity
    receipt["table_receipts"] = table_receipts
    receipt["failures"] = failures
    receipt["baseline_reproduction_failures"] = [x for x in failures if x.get("scenario") == "finite_cutoff_reference"]
    receipt["expected_tov_case_count"] = len(selected) * (1 + len(SCENARIOS))
    receipt["completed_tov_case_count"] = len(summary_rows)
    receipt["status"] = "PASS_TOV_CONTINUATION_EXECUTION" if not failures and len(summary_rows) == receipt["expected_tov_case_count"] else "FAIL_TOV_CONTINUATION_EXECUTION"
    persist_receipt(receipt_path, receipt)
    write_csv(args.output_dir / "surface_continuation_summary.csv", summary_rows)
    write_csv(args.output_dir / "surface_continuation_sensitivity.csv", sensitivity)
    manifest = {
        "schema": "ns-cowling-surface-continuation-manifest-v1",
        "receipt_sha256": sha256(receipt_path),
        "summary_csv_sha256": sha256(args.output_dir / "surface_continuation_summary.csv"),
        "sensitivity_csv_sha256": sha256(args.output_dir / "surface_continuation_sensitivity.csv"),
        "table_count": len(table_receipts),
        "source_summary_sha256": receipt["source_summary_sha256"],
        "source_commit": receipt["source_commit"],
    }
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": receipt["status"],
        "expected_tov_case_count": receipt["expected_tov_case_count"],
        "completed_tov_case_count": receipt["completed_tov_case_count"],
        "baseline_reproduction_failures": len(receipt["baseline_reproduction_failures"]),
        "continuation_failures": len(failures) - len(receipt["baseline_reproduction_failures"]),
        "sensitivity_pair_count": len(sensitivity),
        "artifact_manifest": str(args.output_dir / "manifest.json"),
    }, indent=2))
    return 0 if receipt["status"] == "PASS_TOV_CONTINUATION_EXECUTION" else 1


if __name__ == "__main__":
    raise SystemExit(main())
