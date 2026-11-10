#!/usr/bin/env python3
"""Run a matched-cutoff TOV pilot using SI-native cellwise minimal pressure repair."""
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
sys.path.insert(0, str(ROOT / "python_frontend"))
sys.path.insert(0, str(ROOT / "scripts"))

from compose_cowling.native_table_cutoff_scan import (  # noqa: E402
    load_native_eos_rows,
    select_native_density_knots,
)
from compose_cowling.solver import solve_v5_tov_star  # noqa: E402
from r16_minimal_hermite_pressure_repair import (  # noqa: E402
    audit_hermite_cell,
    minimal_admissible_delta_p,
)

SUMMARY = ROOT / "oscillation_release_v1/data/spectrum_summary_234.csv"
SOURCE_DIR = ROOT / "oscillation_release_v1/data/eos_original_normalized"
PROJECTED_DIR = ROOT / "oscillation_release_v1/data/eos_projected"
SUMMARY_SHA256 = "ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965"
METHODS = ("css", "cubic_F", "cubic_P", "maxwell", "quintic_F", "quintic_P")
PILOT_PAIRS = (
    ("GPPVA_FSU2", "BBB_BHF2"),
    ("SPG_M2", "BBB_BHF2"),
    ("GPPVA_FSU2", "GM_GM1"),
)
EXPECTED_ROWS = (1, 9, 17, 23, 52, 73, 116, 217)
CUTOFF_COUNT = 8
REPAIR_MARGIN_FRACTION = 1.0e-8
MEV_FM3_TO_PA = 1.602176634e32
SOLAR_MASS_KG = 1.988409870698051e30
SOLVER = {"r0_m": 10.0, "r_limit_m": 1.0e8, "rtol": 1.0e-9, "step_cap": 0.05}
REQUIRED = (
    "source_row", "pressure_pa", "energy_j_m3", "n_b_m3", "mu_b_j",
    "pressure_mev_fm3", "energy_mev_fm3", "n_b_fm3", "mu_b_mev",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        fields = list(reader.fieldnames or [])
        missing = [x for x in REQUIRED if x not in fields]
        if missing:
            raise ValueError(f"{path.name}: missing columns {missing}")
        rows = list(reader)
    if len(rows) < 3:
        raise ValueError(f"{path.name}: fewer than three EOS rows")
    for i, row in enumerate(rows, start=1):
        if int(row["source_row"]) != i:
            raise ValueError(f"{path.name}: source_row order mismatch at {i}")
    return fields, rows


def repair_normalized_rows(
    rows: list[dict[str, str]],
    *,
    margin_fraction: float = REPAIR_MARGIN_FRACTION,
) -> tuple[list[dict[str, str]], dict]:
    """Return a derived table whose pressure increments are clipped in SI units."""
    mu = [float(r["mu_b_j"]) for r in rows]
    n = [float(r["n_b_m3"]) for r in rows]
    p_source = [float(r["pressure_pa"]) for r in rows]
    eps_source = [float(r["energy_j_m3"]) for r in rows]
    p_repaired = [p_source[0]]
    adjusted = 0
    preserved = 0
    total_abs_increment_change = 0.0
    max_relative_increment_change = 0.0
    min_d2p = float("inf")
    min_margin = float("inf")
    cell_records = []

    for i in range(len(rows) - 1):
        source_cell = audit_hermite_cell(
            mu0=mu[i], p0=p_source[i], n0=n[i],
            mu1=mu[i + 1], p1=p_source[i + 1], n1=n[i + 1],
        )
        delta_source = p_source[i + 1] - p_source[i]
        repair = minimal_admissible_delta_p(
            mu0=mu[i], n0=n[i], mu1=mu[i + 1], n1=n[i + 1],
            delta_p_source=delta_source,
            interior_margin_fraction=margin_fraction,
        )
        p_repaired.append(p_repaired[-1] + repair["delta_p_repaired"])
        adjusted += int(repair["was_adjusted"])
        preserved += int(not repair["was_adjusted"])
        total_abs_increment_change += abs(repair["delta_p_adjustment"])
        max_relative_increment_change = max(
            max_relative_increment_change,
            abs(repair["delta_p_adjustment"]) / max(abs(delta_source), 1.0e-300),
        )

        repaired_cell = audit_hermite_cell(
            mu0=mu[i], p0=p_repaired[i], n0=n[i],
            mu1=mu[i + 1], p1=p_repaired[i + 1], n1=n[i + 1],
        )
        if repaired_cell["status"] != "PASS_NATIVE_CELL_ADMISSIBLE":
            raise RuntimeError(f"repaired table cell {i + 1} fails SI Hermite admission")
        min_d2p = min(min_d2p, repaired_cell["min_d2p"])
        min_margin = min(min_margin, repaired_cell["min_causal_margin"])
        cell_records.append({
            "source_cell": i + 1,
            "source_status": source_cell["status"],
            "was_adjusted": repair["was_adjusted"],
            "delta_p_source_pa": delta_source,
            "delta_p_repaired_pa": repair["delta_p_repaired"],
            "lower_bound_pa": repair["delta_p_lower"],
            "upper_bound_pa": repair["delta_p_upper"],
            "repaired_min_d2p_si": repaired_cell["min_d2p"],
            "repaired_min_causal_margin_si": repaired_cell["min_causal_margin"],
        })

    eps_repaired = [m * nb - pressure for m, nb, pressure in zip(mu, n, p_repaired)]
    if any(x <= 0.0 for x in p_repaired + eps_repaired):
        raise ValueError("repaired pressure/energy density became nonpositive")
    if any(b <= a for a, b in zip(p_repaired, p_repaired[1:])):
        raise ValueError("repaired pressure is not strictly increasing")
    if any(b <= a for a, b in zip(eps_repaired, eps_repaired[1:])):
        raise ValueError("repaired energy density is not strictly increasing")

    factor = MEV_FM3_TO_PA
    output_rows = []
    for row, p_si, eps_si in zip(rows, p_repaired, eps_repaired):
        out = dict(row)
        p_mev = p_si / factor
        eps_mev = eps_si / factor
        out["pressure_pa"] = format(p_si, ".17g")
        out["pressure_mev_fm3"] = format(p_mev, ".17g")
        out["energy_j_m3"] = format(eps_si, ".17g")
        out["energy_mev_fm3"] = format(eps_mev, ".17g")
        n_mev = float(row["n_b_fm3"])
        mu_mev = float(row["mu_b_mev"])
        out["q1_p_over_nb_mev"] = format(p_mev / n_mev, ".17g")
        if "q3_mu_b_over_mn_minus_1" in row:
            q3 = float(row["q3_mu_b_over_mn_minus_1"])
            mn_mev = mu_mev / (1.0 + q3)
            out["q7_e_over_nb_mn_minus_1"] = format(eps_mev / (n_mev * mn_mev) - 1.0, ".17g")
        output_rows.append(out)

    return output_rows, {
        "row_count": len(rows),
        "cell_count": len(rows) - 1,
        "adjusted_cells": adjusted,
        "preserved_cells": preserved,
        "max_relative_increment_change": max_relative_increment_change,
        "sum_abs_increment_change_pa": total_abs_increment_change,
        "minimum_repaired_d2p_si": min_d2p,
        "minimum_repaired_causal_margin_si": min_margin,
        "source_first_law_max_relative": max(
            abs(e + p - m * nb) / max(abs(e + p), 1.0e-300)
            for e, p, m, nb in zip(eps_source, p_source, mu, n)
        ),
        "repaired_first_law_max_relative": max(
            abs(e + p - m * nb) / max(abs(e + p), 1.0e-300)
            for e, p, m, nb in zip(eps_repaired, p_repaired, mu, n)
        ),
        "p_source_row217_pa": p_source[216],
        "p_source_row218_pa": p_source[217],
        "p_repaired_row217_pa": p_repaired[216],
        "p_repaired_row218_pa": p_repaired[217],
        "cell_records": cell_records,
    }


def write_csv(path: Path, rows: list[dict], fields: tuple[str, ...] | list[str] | None = None) -> None:
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row)) if rows else ["status"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_references() -> list[dict[str, str]]:
    if sha256(SUMMARY) != SUMMARY_SHA256:
        raise RuntimeError("frozen spectrum summary SHA-256 mismatch")
    with SUMMARY.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 234 or len({r.get("name") for r in rows}) != 234:
        raise RuntimeError("frozen spectrum reference cohort mismatch")
    if any(r.get("projected_native_cell_audit") != "PASS" for r in rows):
        raise RuntimeError("frozen projected native-cell PASS missing")
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r15-results-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "dist/r16-minimal-repair-tov-pilot")
    args = parser.parse_args(argv)

    refs = load_references()
    ref_by_name = {r["name"]: r for r in refs}
    with args.r15_results_csv.open(newline="", encoding="utf-8") as stream:
        r15_rows = list(csv.DictReader(stream))
    if len(r15_rows) != 1872 or len({r["case_id"] for r in r15_rows}) != 1872:
        raise RuntimeError("R15 comparator must contain 1,872 unique case IDs")
    r15_by_case = {r["case_id"]: r for r in r15_rows}

    out = args.output_dir.resolve()
    eos_out = out / "eos_repaired"
    profile_out = out / "backgrounds"
    eos_out.mkdir(parents=True, exist_ok=True)
    profile_out.mkdir(parents=True, exist_ok=True)

    selected_tables = []
    table_rows = []
    results = []
    profiles = []
    failures = []

    for crust, core in PILOT_PAIRS:
        for method in METHODS:
            name = f"{crust}__{core}__{method}"
            ref = ref_by_name[name]
            source_path = SOURCE_DIR / f"{name}.csv"
            projected_path = PROJECTED_DIR / f"{name}.csv"
            source_before = sha256(source_path)
            projected_before = sha256(projected_path)
            if source_before != ref["raw_normalized_eos_sha256"]:
                raise ValueError(f"{name}: normalized source hash mismatch")
            if projected_before != ref["projected_eos_sha256"]:
                raise ValueError(f"{name}: projected source hash mismatch")

            fields, source_rows = read_rows(source_path)
            repaired_rows, repair_summary = repair_normalized_rows(source_rows)
            derived_path = eos_out / f"{name}.csv"
            write_csv(derived_path, repaired_rows, fields)
            derived_sha = sha256(derived_path)

            # Native row selection and cutoff pressure are taken from the frozen
            # projected reference. This compares both EOS representations at
            # the same absolute P_cut and the same central pressure.
            projected_rows = load_native_eos_rows(projected_path)
            knots = select_native_density_knots(
                projected_rows, target_count=CUTOFF_COUNT, maximum_density_fm3=1.0e-3
            )
            if tuple(int(k.source_row) for k in knots) != EXPECTED_ROWS:
                raise RuntimeError(f"{name}: frozen cutoff row signature changed")
            pc_pa = float(ref["pc_pa"])

            table_receipt = {
                "name": name, "crust": crust, "core": core, "matching_method": method,
                "normalized_source_sha256": source_before,
                "projected_reference_sha256": projected_before,
                "repaired_derived_sha256": derived_sha,
                "row_count": len(repaired_rows),
                **{k: v for k, v in repair_summary.items() if k != "cell_records"},
                "solver_input_role": "temporary diagnostic derived table; not an authoritative source/projected table",
                "source_unchanged": sha256(source_path) == source_before,
                "projected_reference_unchanged": sha256(projected_path) == projected_before,
            }
            table_rows.append(table_receipt)
            selected_tables.append(name)

            for knot in knots:
                row_number = int(knot.source_row)
                case_id = f"{name}__row_{row_number}"
                comparator = r15_by_case.get(case_id)
                if comparator is None:
                    failures.append({"case_id": case_id, "stage": "reference_lookup", "error": "missing R15 comparator"})
                    continue
                p_cut_pa = float(comparator["p_cut_pa"])
                if not (0.0 < p_cut_pa < pc_pa):
                    failures.append({"case_id": case_id, "stage": "pressure_domain", "error": "P_cut not between 0 and P_c"})
                    continue

                profile_path = profile_out / f"{case_id}.csv"
                try:
                    star = solve_v5_tov_star(
                        derived_path,
                        pc_pa,
                        p_cut_pa=p_cut_pa,
                        background_output_path=profile_path,
                        **SOLVER,
                    )
                    profile = star.get("background_profile")
                    if not isinstance(profile, dict) or not profile_path.is_file():
                        raise RuntimeError("native solver did not export requested background")
                    profile_sha = sha256(profile_path)
                    profile_rows = int(profile.get("row_count", 0))
                    if profile_sha != profile.get("sha256") or profile_rows < 2:
                        raise RuntimeError("background profile SHA/row-count validation failed")
                    mass = float(star["mass_kg"])
                    radius = float(star["radius_m"])
                    if not all(math.isfinite(x) and x > 0.0 for x in (mass, radius)):
                        raise RuntimeError("TOV returned nonfinite/nonpositive M or R")

                    ref_mass = float(comparator["mass_inside_isobar_kg"])
                    ref_radius = float(comparator["isobar_radius_m"])
                    record = {
                        "case_id": case_id,
                        "crust": crust, "core": core, "matching_method": method,
                        "source_row": row_number,
                        "pc_pa": pc_pa,
                        "p_cut_pa_common_with_r15": p_cut_pa,
                        "projected_native_row_pressure_pa": knot.pressure_pa,
                        "repaired_native_row_pressure_pa": float(repaired_rows[knot.row_index]["pressure_pa"]),
                        "same_absolute_cutoff_pressure_as_reference": True,
                        "mass_repaired_kg": mass,
                        "radius_repaired_m": radius,
                        "mass_projected_r15_kg": ref_mass,
                        "radius_projected_r15_m": ref_radius,
                        "delta_mass_pct_vs_projected_r15": 100.0 * (mass / ref_mass - 1.0),
                        "delta_radius_pct_vs_projected_r15": 100.0 * (radius / ref_radius - 1.0),
                        "background_profile_path": profile_path.relative_to(out).as_posix(),
                        "background_profile_sha256": profile_sha,
                        "background_profile_row_count": profile_rows,
                        "backend_build_key": star.get("backend_build_key"),
                        "status": "ok",
                    }
                    results.append(record)
                    profiles.append({
                        "case_id": case_id, "path": record["background_profile_path"],
                        "sha256": profile_sha, "row_count": profile_rows,
                    })
                except Exception as exc:
                    failures.append({
                        "case_id": case_id,
                        "stage": "native_tov",
                        "error": f"{type(exc).__name__}: {exc}"[:2000],
                    })

            if sha256(source_path) != source_before or sha256(projected_path) != projected_before:
                raise RuntimeError(f"{name}: source/reference table mutated during pilot")

    if len(selected_tables) != 18:
        failures.append({"stage": "table_coverage", "error": f"expected 18 tables, got {len(selected_tables)}"})
    if len(results) != 144:
        failures.append({"stage": "case_coverage", "error": f"expected 144 TOV cases, got {len(results)}"})

    by_pair_cutoff = defaultdict(list)
    for row in results:
        by_pair_cutoff[(row["crust"], row["core"], row["source_row"])].append(row)
    summary = []
    for key, rows in sorted(by_pair_cutoff.items()):
        summary.append({
            "crust": key[0], "core": key[1], "source_row": key[2],
            "method_count": len(rows),
            "mean_delta_mass_pct": sum(r["delta_mass_pct_vs_projected_r15"] for r in rows) / len(rows),
            "mean_delta_radius_pct": sum(r["delta_radius_pct_vs_projected_r15"] for r in rows) / len(rows),
            "max_abs_delta_mass_pct": max(abs(r["delta_mass_pct_vs_projected_r15"]) for r in rows),
            "max_abs_delta_radius_pct": max(abs(r["delta_radius_pct_vs_projected_r15"]) for r in rows),
        })

    receipt = {
        "schema": "ns-cowling-r16-minimal-repair-matched-pressure-tov-pilot-v1",
        "status": "PASS_R16_MINIMAL_REPAIR_TOV_PILOT" if len(results) == 144 and len(table_rows) == 18 and not failures else "FAIL_R16_MINIMAL_REPAIR_TOV_PILOT",
        "source_commit": os.environ.get("GITHUB_SHA", "UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256": SUMMARY_SHA256,
        "r15_reference_results_sha256": sha256(args.r15_results_csv),
        "r15_reference_case_count": len(r15_rows),
        "pilot_pairs": [{"crust": c, "core": k} for c, k in PILOT_PAIRS],
        "method_count": len(METHODS),
        "native_cutoff_rows": list(EXPECTED_ROWS),
        "expected_tov_case_count": 144,
        "completed_tov_case_count": len(results),
        "comparison_contract": "same Pc and same absolute P_cut values as the corresponding frozen R15 projected-EOS case; source-row index is a density marker only",
        "repair_contract": "SI-native cellwise closest admissible DeltaP with 1e-8 interval-width interior margin; preserve originally admitted source increments exactly",
        "solver": {"backend": "v5 native Fortran cell_solve task=star", **SOLVER},
        "source_tables_modified": False,
        "projected_reference_tables_modified": False,
        "derived_repaired_table_count": len(table_rows),
        "synthetic_rows_added": 0,
        "radial_stability_or_cowling_modes_computed": False,
        "table_receipts": table_rows,
        "results": results,
        "pair_cutoff_summary": summary,
        "background_profiles": profiles,
        "failures": failures,
    }

    receipt_path = out / "r16_minimal_repair_tov_pilot.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(out / "r16_minimal_repair_tov_pilot_cases.csv", results)
    write_csv(out / "r16_minimal_repair_tov_pilot_tables.csv", table_rows)
    write_csv(out / "r16_minimal_repair_tov_pilot_summary.csv", summary)
    manifest = {
        "schema": "ns-cowling-r16-minimal-repair-tov-pilot-manifest-v1",
        "receipt_sha256": sha256(receipt_path),
        "cases_csv_sha256": sha256(out / "r16_minimal_repair_tov_pilot_cases.csv"),
        "tables_csv_sha256": sha256(out / "r16_minimal_repair_tov_pilot_tables.csv"),
        "summary_csv_sha256": sha256(out / "r16_minimal_repair_tov_pilot_summary.csv"),
        "source_summary_sha256": SUMMARY_SHA256,
        "r15_reference_results_sha256": receipt["r15_reference_results_sha256"],
        "table_count": len(table_rows),
        "case_count": len(results),
        "profile_count": len(profiles),
        "derived_eos_table_sha256": {
            row["name"]: row["repaired_derived_sha256"] for row in table_rows
        },
        "background_profiles": profiles,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": receipt["status"],
        "table_count": len(table_rows),
        "case_count": len(results),
        "profile_count": len(profiles),
        "failure_count": len(failures),
        "summary_by_cutoff": summary,
    }, indent=2), flush=True)
    return 0 if receipt["status"].startswith("PASS_") else 1


if __name__ == "__main__":
    raise SystemExit(main())
