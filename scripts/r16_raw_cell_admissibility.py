"""Apply the exact Fortran native cubic-Hermite cell admission test to both EOS cohorts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "oscillation_release_v1/data/spectrum_summary_234.csv"
SOURCE_DIR = ROOT / "oscillation_release_v1/data/eos_original_normalized"
PROJECTED_DIR = ROOT / "oscillation_release_v1/data/eos_projected"
SUMMARY_SHA256 = "ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965"
EXPECTED_TABLES = 234
REQUIRED_COLUMNS = (
    "source_row", "mu_b_mev", "pressure_mev_fm3", "n_b_fm3",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def audit_hermite_cell(
    *, mu0: float, p0: float, n0: float,
    mu1: float, p1: float, n1: float,
) -> dict[str, Any]:
    """Reproduce ns_mu_bridge.construct_bridge_endpoints exactly in one unit system."""
    values = [float(mu0), float(p0), float(n0), float(mu1), float(p1), float(n1)]
    if any(not math.isfinite(x) or x <= 0.0 for x in values):
        raise ValueError("cell endpoints must be finite and positive")
    mu0, p0, n0, mu1, p1, n1 = values
    if mu1 <= mu0 or p1 <= p0:
        raise ValueError("cell endpoints must be strictly ordered in mu and pressure")

    h = mu1 - mu0
    a0 = p0
    a1 = h * n0
    a2 = 3.0 * (p1 - p0) - h * (2.0 * n0 + n1)
    a3 = 2.0 * (p0 - p1) + h * (n0 + n1)
    d20 = 2.0 * a2 / (h * h)
    d21 = (2.0 * a2 + 6.0 * a3) / (h * h)
    g0 = mu0 * d20 - n0
    g1 = mu1 * d21 - n1
    min_d2 = min(d20, d21)
    min_g = min(g0, g1)
    stable = min_d2 > 0.0
    causal = stable and min_g >= 0.0
    if not stable:
        status = "REJECT_NONPOSITIVE_COMPRESSIBILITY"
    elif not causal:
        status = "REJECT_ACAUSAL_NATIVE_CELL"
    else:
        status = "PASS_NATIVE_CELL_ADMISSIBLE"

    chi0 = n0 / (mu0 * d20) if d20 > 0.0 else None
    chi1 = n1 / (mu1 * d21) if d21 > 0.0 else None
    return {
        "status": status,
        "stable_interval": stable,
        "causal_interval": causal,
        "h_mu": h,
        "a2": a2,
        "a3": a3,
        "d2p_low": d20,
        "d2p_high": d21,
        "min_d2p": min_d2,
        "causal_margin_low": g0,
        "causal_margin_high": g1,
        "min_causal_margin": min_g,
        "chi_low": chi0,
        "chi_high": chi1,
    }


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = [x for x in REQUIRED_COLUMNS if x not in (reader.fieldnames or ())]
        if missing:
            raise ValueError(f"{path.name}: missing columns {missing}")
        rows = list(reader)
    if len(rows) < 3:
        raise ValueError(f"{path.name}: fewer than three EOS rows")
    return rows


def audit_table(path: Path, expected_sha: str, *, name: str, cohort: str) -> tuple[dict, list[dict]]:
    before = sha256(path)
    if before != expected_sha:
        raise ValueError(f"{cohort} hash mismatch: expected {expected_sha}, got {before}")
    rows = read_table(path)
    arrays = {
        field: [float(row[field]) for row in rows]
        for field in ("mu_b_mev", "pressure_mev_fm3", "n_b_fm3")
    }
    if any(not math.isfinite(x) for values in arrays.values() for x in values):
        raise ValueError("nonfinite thermodynamic input")
    for i, row in enumerate(rows, start=1):
        if int(row["source_row"]) != i:
            raise ValueError(f"source_row order mismatch at row {i}")

    cell_records = []
    statuses = Counter()
    min_d2 = float("inf")
    min_g = float("inf")
    min_d2_cell = None
    min_g_cell = None
    for i in range(len(rows) - 1):
        try:
            result = audit_hermite_cell(
                mu0=arrays["mu_b_mev"][i],
                p0=arrays["pressure_mev_fm3"][i],
                n0=arrays["n_b_fm3"][i],
                mu1=arrays["mu_b_mev"][i + 1],
                p1=arrays["pressure_mev_fm3"][i + 1],
                n1=arrays["n_b_fm3"][i + 1],
            )
        except Exception as exc:
            result = {
                "status": "REJECT_INVALID_ENDPOINTS",
                "stable_interval": False,
                "causal_interval": False,
                "error": f"{type(exc).__name__}: {exc}",
            }
        statuses[result["status"]] += 1
        d2 = result.get("min_d2p")
        g = result.get("min_causal_margin")
        if d2 is not None and d2 < min_d2:
            min_d2, min_d2_cell = d2, i + 1
        if g is not None and g < min_g:
            min_g, min_g_cell = g, i + 1
        if result["status"] != "PASS_NATIVE_CELL_ADMISSIBLE":
            cell_records.append({
                "cohort": cohort,
                "name": name,
                "cell_index_one_based": i + 1,
                "source_row_low": i + 1,
                "source_row_high": i + 2,
                "mu_low_mev": arrays["mu_b_mev"][i],
                "mu_high_mev": arrays["mu_b_mev"][i + 1],
                "p_low_mev_fm3": arrays["pressure_mev_fm3"][i],
                "p_high_mev_fm3": arrays["pressure_mev_fm3"][i + 1],
                "n_low_fm3": arrays["n_b_fm3"][i],
                "n_high_fm3": arrays["n_b_fm3"][i + 1],
                **result,
            })

    after = sha256(path)
    if after != before:
        raise RuntimeError("source bytes changed during read-only cell audit")
    pass_count = statuses["PASS_NATIVE_CELL_ADMISSIBLE"]
    summary = {
        "cohort": cohort,
        "name": name,
        "row_count": len(rows),
        "cell_count": len(rows) - 1,
        "pass_cell_count": pass_count,
        "nonpositive_compressibility_cell_count": statuses["REJECT_NONPOSITIVE_COMPRESSIBILITY"],
        "acausal_cell_count": statuses["REJECT_ACAUSAL_NATIVE_CELL"],
        "invalid_endpoint_cell_count": statuses["REJECT_INVALID_ENDPOINTS"],
        "all_cells_admissible": pass_count == len(rows) - 1,
        "minimum_cell_d2p": min_d2 if min_d2 != float("inf") else None,
        "minimum_cell_d2p_location": min_d2_cell,
        "minimum_cell_causal_margin": min_g if min_g != float("inf") else None,
        "minimum_cell_causal_margin_location": min_g_cell,
        "sha256_expected": expected_sha,
        "sha256_before": before,
        "sha256_after": after,
        "unchanged": before == after,
    }
    return summary, cell_records


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("status\n", encoding="utf-8")
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "dist/r16-raw-cell-admissibility")
    args = parser.parse_args(argv)

    if sha256(SUMMARY) != SUMMARY_SHA256:
        raise RuntimeError("frozen 234-row summary SHA-256 mismatch")
    with SUMMARY.open(newline="", encoding="utf-8") as stream:
        refs = list(csv.DictReader(stream))
    names = {r["name"] for r in refs}
    if len(refs) != EXPECTED_TABLES or len(names) != EXPECTED_TABLES:
        raise RuntimeError("frozen cohort row count or names changed")

    table_records, bad_cells, failures = [], [], []
    for ref in refs:
        name = ref["name"]
        for cohort, directory, hash_field in (
            ("normalized_source", SOURCE_DIR, "raw_normalized_eos_sha256"),
            ("projected_solver_input", PROJECTED_DIR, "projected_eos_sha256"),
        ):
            path = directory / f"{name}.csv"
            try:
                record, failures_here = audit_table(
                    path, ref[hash_field], name=name, cohort=cohort
                )
                record.update({
                    "crust": ref["crust"],
                    "core": ref["core"],
                    "matching_method": ref["matching_method"],
                })
                table_records.append(record)
                bad_cells.extend(failures_here)
            except Exception as exc:
                failures.append({
                    "cohort": cohort,
                    "name": name,
                    "error": f"{type(exc).__name__}: {exc}"[:2000],
                })

    normalized = [r for r in table_records if r["cohort"] == "normalized_source"]
    projected = [r for r in table_records if r["cohort"] == "projected_solver_input"]
    def cohort_summary(rows: list[dict]) -> dict:
        return {
            "tables": len(rows),
            "total_cells": sum(x["cell_count"] for x in rows),
            "tables_all_cells_admissible": sum(x["all_cells_admissible"] for x in rows),
            "tables_with_any_nonpositive_compressibility": sum(x["nonpositive_compressibility_cell_count"] > 0 for x in rows),
            "tables_with_any_acausal_cell": sum(x["acausal_cell_count"] > 0 for x in rows),
            "pass_cells": sum(x["pass_cell_count"] for x in rows),
            "nonpositive_compressibility_cells": sum(x["nonpositive_compressibility_cell_count"] for x in rows),
            "acausal_cells": sum(x["acausal_cell_count"] for x in rows),
            "invalid_endpoint_cells": sum(x["invalid_endpoint_cell_count"] for x in rows),
            "all_hashes_unchanged": all(x["unchanged"] for x in rows),
            "minimum_d2p": min((x["minimum_cell_d2p"] for x in rows if x["minimum_cell_d2p"] is not None), default=None),
            "minimum_causal_margin": min((x["minimum_cell_causal_margin"] for x in rows if x["minimum_cell_causal_margin"] is not None), default=None),
        }
    ns, ps = cohort_summary(normalized), cohort_summary(projected)

    receipt = {
        "schema": "ns-cowling-r16-raw-normalized-native-cell-admissibility-v1",
        "status": (
            "PASS_R16_RAW_NORMALIZED_CELL_FAILURE_MAP"
            if len(normalized) == EXPECTED_TABLES and len(projected) == EXPECTED_TABLES and not failures
            else "FAIL_R16_RAW_NORMALIZED_CELL_FAILURE_MAP"
        ),
        "source_commit": os.environ.get("GITHUB_SHA", "UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256": SUMMARY_SHA256,
        "fortran_cell_acceptance": "min(d2P/dmu2 at endpoints)>0 and min(mu*d2P/dmu2-dP/dmu at endpoints)>=0",
        "normalized_source": ns,
        "projected_solver_input": ps,
        "failures": failures,
        "table_records": table_records,
        "first_failing_cells": bad_cells[:5000],
        "tov_integrations_executed": 0,
        "synthetic_eos_tables_written": 0,
    }

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    receipt_path = out / "r16_raw_normalized_cell_admissibility_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(out / "r16_raw_normalized_cell_admissibility_tables.csv", table_records)
    write_csv(out / "r16_raw_normalized_cell_failures.csv", bad_cells)
    manifest = {
        "schema": "ns-cowling-r16-raw-cell-admissibility-manifest-v1",
        "receipt_sha256": sha256(receipt_path),
        "table_csv_sha256": sha256(out / "r16_raw_normalized_cell_admissibility_tables.csv"),
        "failure_csv_sha256": sha256(out / "r16_raw_normalized_cell_failures.csv"),
        "source_summary_sha256": SUMMARY_SHA256,
        "normalized_table_count": len(normalized),
        "projected_table_count": len(projected),
        "failing_cell_records": len(bad_cells),
        "tov_integrations_executed": 0,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": receipt["status"],
        "normalized_source": ns,
        "projected_solver_input": ps,
        "cell_failure_record_count": len(bad_cells),
        "audit_failures": len(failures),
    }, indent=2), flush=True)
    return 0 if receipt["status"].startswith("PASS_") else 1


if __name__ == "__main__":
    raise SystemExit(main())
