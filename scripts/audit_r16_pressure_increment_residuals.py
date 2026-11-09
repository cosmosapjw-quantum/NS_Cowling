#!/usr/bin/env python3
"""Decompose normalized-to-projected pressure changes into native-cell increments."""
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


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def values(rows: list[dict[str, str]], key: str) -> list[float]:
    out = [float(row[key]) for row in rows]
    if any(not math.isfinite(x) for x in out):
        raise ValueError(f"{key}: nonfinite values in source table")
    return out


def relative_mismatch(a: list[float], b: list[float]) -> float:
    return max(abs(x-y) / max(abs(y), 1.0e-300) for x, y in zip(a, b))


def audit_one(ref: dict[str, str]) -> tuple[dict, list[dict]]:
    name = ref["name"]
    source = SOURCE_DIR / f"{name}.csv"
    projected = PROJECTED_DIR / f"{name}.csv"
    if not source.is_file() or not projected.is_file():
        raise FileNotFoundError(f"{name}: normalized or projected table missing")

    source_sha_before, projected_sha_before = digest(source), digest(projected)
    if source_sha_before != ref.get("raw_normalized_eos_sha256"):
        raise ValueError(f"{name}: normalized source hash mismatch")
    if projected_sha_before != ref.get("projected_eos_sha256"):
        raise ValueError(f"{name}: projected hash mismatch")

    src, dst = load_csv(source), load_csv(projected)
    if len(src) != len(dst) or len(src) < 3:
        raise ValueError(f"{name}: source/projected row counts differ or are too small")
    for i, (a, b) in enumerate(zip(src, dst), start=1):
        if int(a["source_row"]) != i or int(b["source_row"]) != i:
            raise ValueError(f"{name}: source-row order mismatch at row {i}")
        for key in ("source_row", "i_t", "i_nb", "i_yq"):
            if key in a and key in b and a[key] != b[key]:
                raise ValueError(f"{name}: source identity field {key} changed at row {i}")

    mu = values(src, "mu_b_mev")
    n = values(src, "n_b_fm3")
    p = values(src, "pressure_mev_fm3")
    eps = values(src, "energy_mev_fm3")
    result = reconstruct_projection(mu, n, p, eps)

    n_dst = values(dst, "n_b_fm3")
    p_dst = values(dst, "pressure_mev_fm3")
    eps_dst = values(dst, "energy_mev_fm3")
    mu_dst = values(dst, "mu_b_mev")
    recon_mismatch = max(
        relative_mismatch(n_dst, result["n_projected"]),
        relative_mismatch(p_dst, result["p_projected"]),
        relative_mismatch(eps_dst, result["eps_projected"]),
        relative_mismatch(mu_dst, mu),
    )
    if recon_mismatch > 5.0e-12:
        raise ValueError(f"{name}: independent projection reconstruction mismatch {recon_mismatch:.6e}")

    source_sha_after, projected_sha_after = digest(source), digest(projected)
    if source_sha_before != source_sha_after or projected_sha_before != projected_sha_after:
        raise RuntimeError(f"{name}: source bytes changed during read-only audit")

    cells = []
    source_cells = result["source_cell_diagnostics"]
    projected_cells = result["projected_cell_diagnostics"]
    for i in range(len(mu)-1):
        dp_src = result["source_pressure_increments"][i]
        dp_proj = result["projected_pressure_increments"][i]
        residual = result["pressure_increment_residuals"][i]
        cumulative = result["cumulative_pressure_residuals"][i+1]
        sd, pd = source_cells[i], projected_cells[i]
        cells.append({
            "name": name,
            "crust": ref["crust"],
            "core": ref["core"],
            "matching_method": ref["matching_method"],
            "lower_source_row": i+1,
            "upper_source_row": i+2,
            "mu_lower_mev": mu[i],
            "mu_upper_mev": mu[i+1],
            "delta_mu_mev": mu[i+1]-mu[i],
            "n_lower_fm3": n[i],
            "n_upper_fm3": n[i+1],
            "p_source_lower_mev_fm3": p[i],
            "p_source_upper_mev_fm3": p[i+1],
            "p_projected_lower_mev_fm3": result["p_projected"][i],
            "p_projected_upper_mev_fm3": result["p_projected"][i+1],
            "delta_p_source_mev_fm3": dp_src,
            "delta_p_projected_mev_fm3": dp_proj,
            "delta_p_residual_mev_fm3": residual,
            "delta_p_residual_over_source_increment": residual/dp_src,
            "cumulative_p_residual_at_upper_mev_fm3": cumulative,
            "cumulative_p_residual_over_source_pressure": cumulative/p[i+1],
            "source_dP_dmu_secant": result["source_pressure_secants_dpdmu"][i],
            "projected_dP_dmu_secant": result["projected_pressure_secants_dpdmu"][i],
            "source_min_d2p": sd["min_d2p"],
            "source_min_causal_margin": sd["min_causal_margin"],
            "source_cell_stable": sd["stable"],
            "source_cell_causal": sd["causal"],
            "projected_min_d2p": pd["min_d2p"],
            "projected_min_causal_margin": pd["min_causal_margin"],
            "projected_cell_stable": pd["stable"],
            "projected_cell_causal": pd["causal"],
        })

    # The interval residuals must telescope to the endpoint pressure difference.
    residual_sum = sum(result["pressure_increment_residuals"])
    endpoint_difference = result["p_projected"][-1] - p[-1]
    telescope_error = abs(residual_sum-endpoint_difference) / max(abs(endpoint_difference), 1.0e-300)

    boundary_index = 216 if len(cells) > 216 else None
    boundary = cells[boundary_index] if boundary_index is not None else None
    endpoint217 = result["cumulative_pressure_residuals"][216] if len(mu) > 216 else None
    endpoint218 = result["cumulative_pressure_residuals"][217] if len(mu) > 217 else None

    row = {
        "name": name,
        "crust": ref["crust"],
        "core": ref["core"],
        "matching_method": ref["matching_method"],
        "row_count": len(src),
        "source_sha256": source_sha_before,
        "projected_sha256": projected_sha_before,
        "source_unchanged": source_sha_before == source_sha_after,
        "projected_unchanged": projected_sha_before == projected_sha_after,
        "reconstruction_max_relative_mismatch": recon_mismatch,
        "pressure_residual_telescope_relative_error": telescope_error,
        "source_nonstable_cell_count": result["source_nonstable_cell_count"],
        "source_noncausal_cell_count": result["source_noncausal_cell_count"],
        "projected_nonstable_cell_count": result["projected_nonstable_cell_count"],
        "projected_noncausal_cell_count": result["projected_noncausal_cell_count"],
        "source_cell_count": len(source_cells),
        "projected_cell_count": len(projected_cells),
        "source_min_d2p_all_cells": result["source_min_d2p_all_cells"],
        "projected_min_d2p_all_cells": result["projected_min_d2p_all_cells"],
        "max_abs_cell_pressure_increment_residual": max(abs(x) for x in result["pressure_increment_residuals"]),
        "max_abs_cell_pressure_increment_relative_residual": max(abs(x) for x in result["pressure_increment_relative_residuals"]),
        "max_abs_cumulative_pressure_residual": max(abs(x) for x in result["cumulative_pressure_residuals"]),
        "endpoint217_cumulative_pressure_residual": endpoint217,
        "endpoint218_cumulative_pressure_residual": endpoint218,
        "endpoint217_pressure_delta_pct": (
            100.0*endpoint217/p[216] if endpoint217 is not None else None
        ),
        "endpoint218_pressure_delta_pct": (
            100.0*endpoint218/p[217] if endpoint218 is not None else None
        ),
        "row217_218_increment_residual": (
            boundary["delta_p_residual_mev_fm3"] if boundary else None
        ),
        "row217_218_increment_residual_over_source_increment": (
            boundary["delta_p_residual_over_source_increment"] if boundary else None
        ),
        "row217_218_cell_source_stable": boundary["source_cell_stable"] if boundary else None,
        "row217_218_cell_source_causal": boundary["source_cell_causal"] if boundary else None,
        "row217_218_cell_projected_stable": boundary["projected_cell_stable"] if boundary else None,
        "row217_218_cell_projected_causal": boundary["projected_cell_causal"] if boundary else None,
    }
    return row, cells


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(rows[0].keys()) if rows else ["status"]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "dist/r16-pressure-increment-residuals")
    args = parser.parse_args(argv)

    if digest(SUMMARY) != SUMMARY_SHA256:
        raise RuntimeError("frozen spectrum summary SHA mismatch")
    refs = load_csv(SUMMARY)
    names = {x["name"] for x in refs}
    if len(refs) != EXPECTED_TABLES or len(names) != EXPECTED_TABLES:
        raise RuntimeError("frozen cohort is not exactly 234 unique tables")
    if {x["matching_method"] for x in refs} != set(METHODS):
        raise RuntimeError("matching method inventory changed")
    if {p.stem for p in SOURCE_DIR.glob("*.csv")} != names:
        raise RuntimeError("normalized source inventory differs from frozen cohort")
    if {p.stem for p in PROJECTED_DIR.glob("*.csv")} != names:
        raise RuntimeError("projected inventory differs from frozen cohort")

    table_rows, interval_rows, failures = [], [], []
    by_name = {x["name"]: x for x in refs}
    for name in sorted(names):
        try:
            table, cells = audit_one(by_name[name])
            table_rows.append(table)
            interval_rows.extend(cells)
        except Exception as exc:
            failures.append({"name": name, "error": f"{type(exc).__name__}: {exc}"[:2000]})

    source_cells = sum(x["source_cell_count"] for x in table_rows)
    projected_cells = sum(x["projected_cell_count"] for x in table_rows)
    expected_cells = sum(x["row_count"]-1 for x in table_rows)
    pass_integrity = (
        len(table_rows) == EXPECTED_TABLES
        and len(interval_rows) == expected_cells
        and not failures
        and all(x["source_unchanged"] and x["projected_unchanged"] for x in table_rows)
        and all(x["reconstruction_max_relative_mismatch"] <= 5.0e-12 for x in table_rows)
        and all(x["pressure_residual_telescope_relative_error"] <= 5.0e-12 for x in table_rows)
    )
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out/"r16_pressure_increment_residuals.csv", interval_rows)
    write_csv(out/"r16_pressure_increment_table_summary.csv", table_rows)

    crust_groups = defaultdict(list)
    for row in table_rows:
        crust_groups[row["crust"]].append(row)
    crust_summary=[]
    for crust, rows in sorted(crust_groups.items()):
        crust_summary.append({
            "crust":crust,
            "table_count":len(rows),
            "source_nonstable_cell_count":sum(x["source_nonstable_cell_count"] for x in rows),
            "source_noncausal_cell_count":sum(x["source_noncausal_cell_count"] for x in rows),
            "projected_nonstable_cell_count":sum(x["projected_nonstable_cell_count"] for x in rows),
            "projected_noncausal_cell_count":sum(x["projected_noncausal_cell_count"] for x in rows),
            "endpoint218_pressure_delta_pct_min":min(x["endpoint218_pressure_delta_pct"] for x in rows),
            "endpoint218_pressure_delta_pct_max":max(x["endpoint218_pressure_delta_pct"] for x in rows),
            "row217_218_increment_residual_over_source_increment_min":min(x["row217_218_increment_residual_over_source_increment"] for x in rows),
            "row217_218_increment_residual_over_source_increment_max":max(x["row217_218_increment_residual_over_source_increment"] for x in rows),
            "maximum_abs_cell_increment_relative_residual":max(x["max_abs_cell_pressure_increment_relative_residual"] for x in rows),
            "maximum_telescope_error":max(x["pressure_residual_telescope_relative_error"] for x in rows),
        })
    write_csv(out/"r16_pressure_increment_crust_summary.csv",crust_summary)

    receipt={
        "schema":"ns-cowling-r16-pressure-increment-residual-audit-v1",
        "status":"PASS_R16_CELLWISE_PRESSURE_INCREMENT_AUDIT" if pass_integrity else "FAIL_R16_CELLWISE_PRESSURE_INCREMENT_AUDIT",
        "source_commit":os.environ.get("GITHUB_SHA","UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256":SUMMARY_SHA256,
        "projector_identity":"weighted_PAV_q_linear_exact_integral_P_v1",
        "audited_table_count":len(table_rows),
        "expected_table_count":EXPECTED_TABLES,
        "audited_native_cell_count":len(interval_rows),
        "expected_native_cell_count":expected_cells,
        "source_nonstable_cell_count":sum(x["source_nonstable_cell_count"] for x in table_rows),
        "source_noncausal_cell_count":sum(x["source_noncausal_cell_count"] for x in table_rows),
        "projected_nonstable_cell_count":sum(x["projected_nonstable_cell_count"] for x in table_rows),
        "projected_noncausal_cell_count":sum(x["projected_noncausal_cell_count"] for x in table_rows),
        "max_reconstruction_relative_mismatch":max((x["reconstruction_max_relative_mismatch"] for x in table_rows),default=None),
        "max_pressure_residual_telescope_relative_error":max((x["pressure_residual_telescope_relative_error"] for x in table_rows),default=None),
        "max_abs_cell_increment_relative_residual":max((x["max_abs_cell_pressure_increment_relative_residual"] for x in table_rows),default=None),
        "endpoint218_pressure_delta_pct_range":[min(x["endpoint218_pressure_delta_pct"] for x in table_rows),max(x["endpoint218_pressure_delta_pct"] for x in table_rows)] if table_rows else None,
        "row217_218_increment_residual_over_source_increment_range":[min(x["row217_218_increment_residual_over_source_increment"] for x in table_rows),max(x["row217_218_increment_residual_over_source_increment"] for x in table_rows)] if table_rows else None,
        "source_tables_unchanged":all(x["source_unchanged"] for x in table_rows),
        "projected_tables_unchanged":all(x["projected_unchanged"] for x in table_rows),
        "tov_integrations_executed":0,
        "synthetic_eos_tables_written":0,
        "failures":failures,
    }
    receipt_path=out/"r16_pressure_increment_receipt.json"
    receipt_path.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    manifest={
        "schema":"ns-cowling-r16-pressure-increment-manifest-v1",
        "receipt_sha256":digest(receipt_path),
        "interval_csv_sha256":digest(out/"r16_pressure_increment_residuals.csv"),
        "table_summary_sha256":digest(out/"r16_pressure_increment_table_summary.csv"),
        "crust_summary_sha256":digest(out/"r16_pressure_increment_crust_summary.csv"),
        "source_summary_sha256":SUMMARY_SHA256,
        "audited_table_count":len(table_rows),
        "audited_native_cell_count":len(interval_rows),
    }
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2),flush=True)
    return 0 if pass_integrity else 1


if __name__=="__main__":
    raise SystemExit(main())
