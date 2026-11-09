#!/usr/bin/env python3
"""Run a bounded q-joint-fit sensitivity pilot on four frozen representative EOS tables."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from r16_projection_audit import hermite_cell_diagnostics  # noqa: E402
from r17_joint_q_fit import fit_monotone_q  # noqa: E402

SUMMARY = ROOT / "oscillation_release_v1/data/spectrum_summary_234.csv"
SOURCE_DIR = ROOT / "oscillation_release_v1/data/eos_original_normalized"
PROJECTED_DIR = ROOT / "oscillation_release_v1/data/eos_projected"
SUMMARY_SHA256 = "ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965"
SELECTED = (
    "GPPVA_DD2__BBB_BHF2__css",
    "GPPVA_TM1e__BBB_BHF2__css",
    "GPPVA_TM1e__BBB_BHF2__maxwell",
    "SPG_M1__BBB_BHF2__css",
)
LAMBDAS = (1.0, 0.75, 0.5, 0.25, 0.1)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(rows[0]) if rows else ["status"]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "dist/r17-joint-q-fit-pilot")
    parser.add_argument("--max-iterations", type=int, default=50000)
    parser.add_argument("--tolerance", type=float, default=1.0e-8)
    args = parser.parse_args(argv)

    if sha256(SUMMARY) != SUMMARY_SHA256:
        raise RuntimeError("frozen 234-row spectrum summary SHA-256 mismatch")
    summary_rows = read_csv(SUMMARY)
    by_name = {row["name"]: row for row in summary_rows}
    if len(by_name) != 234:
        raise RuntimeError("frozen cohort cardinality is not 234")
    missing = sorted(set(SELECTED) - set(by_name))
    if missing:
        raise RuntimeError(f"selected representative tables are absent: {missing}")

    summaries, knot_rows, failures = [], [], []
    source_hashes_unchanged = True
    projected_hashes_unchanged = True
    baseline_checks_pass = True
    all_candidates_admissible = True

    for name in SELECTED:
        ref = by_name[name]
        src = SOURCE_DIR / f"{name}.csv"
        proj = PROJECTED_DIR / f"{name}.csv"
        src_before, proj_before = sha256(src), sha256(proj)
        if src_before != ref["raw_normalized_eos_sha256"]:
            raise RuntimeError(f"{name}: normalized source hash mismatch")
        if proj_before != ref["projected_eos_sha256"]:
            raise RuntimeError(f"{name}: projected source hash mismatch")

        raw = read_csv(src)
        projected = read_csv(proj)
        mu = np.asarray([float(x["mu_b_mev"]) for x in raw])
        n = np.asarray([float(x["n_b_fm3"]) for x in raw])
        p = np.asarray([float(x["pressure_mev_fm3"]) for x in raw])
        eps = np.asarray([float(x["energy_mev_fm3"]) for x in raw])
        if len(raw) != len(projected) or len(raw) < 3:
            raise RuntimeError(f"{name}: raw/projected row counts disagree")
        if [int(x["source_row"]) for x in raw] != list(range(1, len(raw) + 1)):
            raise RuntimeError(f"{name}: source row order is invalid")

        for lam in LAMBDAS:
            result = fit_monotone_q(
                mu, n, p, lambda_density=lam,
                max_iterations=args.max_iterations, tolerance=args.tolerance,
            )
            cells = [
                hermite_cell_diagnostics(
                    float(mu[i]), float(result["p_fitted"][i]), float(result["n_fitted"][i]),
                    float(mu[i + 1]), float(result["p_fitted"][i + 1]), float(result["n_fitted"][i + 1]),
                )
                for i in range(len(mu) - 1)
            ]
            nonstable = sum(not x["stable"] for x in cells)
            noncausal = sum(not x["causal"] for x in cells)
            source_cell_count = len(cells)

            eps_fit = np.asarray(result["epsilon_fitted"])
            n_fit = np.asarray(result["n_fitted"])
            p_fit = np.asarray(result["p_fitted"])
            mu_fit = np.asarray(result["mu_b"])
            first_law = float(np.max(
                np.abs(eps_fit + p_fit - mu_fit * n_fit)
                / np.maximum(np.abs(eps_fit + p_fit), 1.0e-300)
            ))
            pressure_secants = np.diff(p_fit) / np.diff(mu_fit)
            secant_decreases = int(np.sum(
                np.diff(pressure_secants) < -1.0e-12 * np.maximum(np.abs(pressure_secants[:-1]), 1.0e-300)
            ))

            if not (
                result["q_monotone"] and result["positive_thermodynamics"]
                and result["pressure_strictly_increasing"]
                and result["energy_density_strictly_increasing"]
                and nonstable == 0 and noncausal == 0
                and first_law <= 2.0e-12 and secant_decreases == 0
            ):
                all_candidates_admissible = False
                raise RuntimeError(f"{name}, lambda={lam}: candidate violated acceptance")

            # Lambda=1 is baseline A and must reproduce the frozen projected solver table.
            baseline_n_mismatch = baseline_p_mismatch = baseline_eps_mismatch = None
            if lam == 1.0:
                p_proj = np.asarray([float(x["pressure_mev_fm3"]) for x in projected])
                n_proj = np.asarray([float(x["n_b_fm3"]) for x in projected])
                eps_proj = np.asarray([float(x["energy_mev_fm3"]) for x in projected])
                baseline_n_mismatch = float(np.max(np.abs(n_fit / n_proj - 1.0)))
                baseline_p_mismatch = float(np.max(np.abs(p_fit / p_proj - 1.0)))
                baseline_eps_mismatch = float(np.max(np.abs(eps_fit / eps_proj - 1.0)))
                if max(baseline_n_mismatch, baseline_p_mismatch, baseline_eps_mismatch) > 5.0e-12:
                    baseline_checks_pass = False
                    raise RuntimeError(f"{name}: lambda=1 does not reproduce frozen projected EOS")

            source_secants = np.diff(p) / np.diff(mu)
            source_secant_drop_count = int(np.sum(
                np.diff(source_secants) < -1.0e-12 * np.maximum(np.abs(source_secants[:-1]), 1.0e-300)
            ))
            summaries.append({
                "name": name,
                "crust": ref["crust"],
                "core": ref["core"],
                "matching_method": ref["matching_method"],
                "lambda_density": lam,
                "objective_before_source_q": result["objective_before_source_q"],
                "objective_after_fit": result["objective_after_fit"],
                "objective_reduction_fraction": result["objective_reduction_fraction"],
                "density_objective_after": result["density_objective_after"],
                "pressure_objective_after": result["pressure_objective_after"],
                "iterations": result["iterations"],
                "kkt_residual": result["projected_gradient_kkt_residual"],
                "density_rms_relative_change": result["density_rms_relative_change"],
                "density_max_abs_relative_change": result["density_max_abs_relative_change"],
                "pressure_rms_relative_change": result["pressure_rms_relative_change"],
                "pressure_max_abs_relative_change": result["pressure_max_abs_relative_change"],
                "energy_rms_relative_change": result["energy_rms_relative_change"],
                "source_pressure_secant_decrease_count": source_secant_drop_count,
                "candidate_pressure_secant_decrease_count": secant_decreases,
                "candidate_cell_count": source_cell_count,
                "candidate_nonstable_cell_count": nonstable,
                "candidate_noncausal_cell_count": noncausal,
                "first_law_max_relative_residual": first_law,
                "lambda_one_n_max_relative_mismatch": baseline_n_mismatch,
                "lambda_one_p_max_relative_mismatch": baseline_p_mismatch,
                "lambda_one_epsilon_max_relative_mismatch": baseline_eps_mismatch,
                "normalized_source_sha256": src_before,
                "projected_parent_sha256": proj_before,
                "status": "PASS_CANDIDATE_ADMISSIBLE",
            })

            for i in range(len(mu)):
                knot_rows.append({
                    "name": name,
                    "lambda_density": lam,
                    "source_row": i + 1,
                    "mu_b_mev": float(mu_fit[i]),
                    "n_source_fm3": float(n[i]),
                    "n_candidate_fm3": float(n_fit[i]),
                    "n_relative_change": float(n_fit[i] / n[i] - 1.0),
                    "p_source_mev_fm3": float(p[i]),
                    "p_candidate_mev_fm3": float(p_fit[i]),
                    "p_relative_change": float(p_fit[i] / p[i] - 1.0),
                    "epsilon_source_mev_fm3": float(eps[i]),
                    "epsilon_candidate_mev_fm3": float(eps_fit[i]),
                })

        source_hashes_unchanged &= sha256(src) == src_before
        projected_hashes_unchanged &= sha256(proj) == proj_before

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    receipt = {
        "schema": "ns-cowling-r17-joint-q-fit-pilot-v1",
        "status": "PASS_R17_JOINT_Q_FIT_PILOT" if (
            len(summaries) == len(SELECTED) * len(LAMBDAS)
            and baseline_checks_pass and all_candidates_admissible
            and source_hashes_unchanged and projected_hashes_unchanged
        ) else "FAIL_R17_JOINT_Q_FIT_PILOT",
        "source_commit": os.environ.get("GITHUB_SHA", "UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256": SUMMARY_SHA256,
        "selected_tables": list(SELECTED),
        "lambda_grid": list(LAMBDAS),
        "lambda_interpretation": "sensitivity parameter only; no statistical uncertainty model is assumed",
        "objective": (
            "lambda * mu-grid-weighted relative n_B residual^2 "
            "+ (1-lambda) * mu-grid-weighted relative P residual^2"
        ),
        "constraints": [
            "q=n_B/mu_B positive and nondecreasing",
            "P=P0+integral(mu_B*q_linear dmu_B)",
            "n_B=mu_B*q",
            "epsilon=mu_B*n_B-P",
            "whole-cell cubic-Hermite P''>0",
            "whole-cell mu_B*P''-P'>=0",
        ],
        "candidate_count": len(summaries),
        "candidate_knots_written": len(knot_rows),
        "baseline_lambda_one_reproduces_frozen_projected": baseline_checks_pass,
        "all_candidates_admissible": all_candidates_admissible,
        "source_tables_unchanged": source_hashes_unchanged,
        "projected_parent_tables_unchanged": projected_hashes_unchanged,
        "tov_integrations_executed": 0,
        "oscillation_modes_executed": 0,
        "synthetic_or_primary_eos_files_written": 0,
        "candidates": summaries,
        "failures": failures,
    }
    receipt_path = output / "r17_joint_q_fit_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(output / "r17_joint_q_fit_candidates.csv", summaries)
    write_csv(output / "r17_joint_q_fit_knots.csv", knot_rows)
    manifest = {
        "schema": "ns-cowling-r17-joint-q-fit-pilot-manifest-v1",
        "receipt_sha256": sha256(receipt_path),
        "candidates_csv_sha256": sha256(output / "r17_joint_q_fit_candidates.csv"),
        "knots_csv_sha256": sha256(output / "r17_joint_q_fit_knots.csv"),
        "source_summary_sha256": SUMMARY_SHA256,
        "selected_table_count": len(SELECTED),
        "candidate_count": len(summaries),
        "knot_row_count": len(knot_rows),
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": receipt["status"],
        "candidate_count": len(summaries),
        "candidate_knots": len(knot_rows),
        "lambda_one_baseline_pass": baseline_checks_pass,
        "all_candidates_admissible": all_candidates_admissible,
        "source_tables_unchanged": source_hashes_unchanged,
        "projected_parent_tables_unchanged": projected_hashes_unchanged,
    }, indent=2), flush=True)
    return 0 if receipt["status"].startswith("PASS_") else 1


if __name__ == "__main__":
    raise SystemExit(main())
