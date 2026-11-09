#!/usr/bin/env python3
"""Repackage a completed R15 full-cohort artifact into sub-512-MiB shards."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

EXPECTED_CASES = 1872
EXPECTED_TABLES = 234
SHARD_COUNT = 4


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def verify_preflight(root: Path) -> dict:
    manifest_path = root / "manifest.json"
    receipt_path = root / "native_outer_crust_r15_preflight.json"
    tables_path = root / "native_outer_crust_r15_preflight_tables.csv"
    knots_path = root / "native_outer_crust_r15_preflight_selected_knots.csv"
    manifest = load_json(manifest_path)
    for key, path in (
        ("receipt_sha256", receipt_path),
        ("tables_csv_sha256", tables_path),
        ("selected_knots_csv_sha256", knots_path),
    ):
        if sha256(path) != manifest[key]:
            raise ValueError(f"preflight manifest mismatch: {key}")
    receipt = load_json(receipt_path)
    if receipt.get("status") != "PASS_R15_FULL_COHORT_NATIVE_DENSITY_PREFLIGHT":
        raise ValueError("full-cohort preflight receipt is not PASS")
    if receipt.get("audited_projected_table_count") != EXPECTED_TABLES:
        raise ValueError("preflight table count mismatch")
    if receipt.get("selected_native_knot_count") != EXPECTED_CASES:
        raise ValueError("preflight native-knot count mismatch")
    if receipt.get("tov_integrations_executed") != 0:
        raise ValueError("preflight unexpectedly reports TOV integrations")
    return {"manifest": manifest, "receipt": receipt}


def verify_tov(root: Path) -> tuple[dict, dict, list[dict]]:
    receipt_path = root / "native_outer_crust_r15.json"
    result_path = root / "native_outer_crust_r15.csv"
    summary_path = root / "native_outer_crust_r15_pair_cutoff_summary.csv"
    manifest_path = root / "manifest.json"

    receipt = load_json(receipt_path)
    manifest = load_json(manifest_path)

    if receipt.get("status") != "PASS_R15_FULL_COHORT_NATIVE_OUTER_CRUST_CUTOFF_SCAN":
        raise ValueError(f"TOV receipt is not full-cohort PASS: {receipt.get('status')}")
    if receipt.get("expected_tov_case_count") != EXPECTED_CASES:
        raise ValueError("TOV expected case count mismatch")
    if receipt.get("completed_tov_case_count") != EXPECTED_CASES:
        raise ValueError("TOV completed case count mismatch")
    if len(receipt.get("results", [])) != EXPECTED_CASES:
        raise ValueError("TOV receipt results length mismatch")
    if len(receipt.get("source_table_receipts", [])) != EXPECTED_TABLES:
        raise ValueError("TOV source-table receipt count mismatch")
    if len(receipt.get("failures", [])) != 0:
        raise ValueError("TOV receipt contains failures")
    for key in (
        "source_tables_unchanged",
        "all_native_edge_baselines_reproduced",
        "all_pair_cutoff_cells_complete",
    ):
        if receipt.get(key) is not True:
            raise ValueError(f"TOV acceptance flag is not true: {key}")
    if manifest.get("expected_tov_case_count") != EXPECTED_CASES:
        raise ValueError("TOV manifest expected case count mismatch")
    if manifest.get("completed_tov_case_count") != EXPECTED_CASES:
        raise ValueError("TOV manifest completed case count mismatch")
    if manifest.get("source_table_count") != EXPECTED_TABLES:
        raise ValueError("TOV manifest source table count mismatch")
    if manifest.get("background_profile_count") != EXPECTED_CASES:
        raise ValueError("TOV manifest profile count mismatch")
    if manifest.get("source_tables_unchanged") is not True:
        raise ValueError("TOV manifest says source tables changed")

    results = []
    with result_path.open(newline="", encoding="utf-8") as stream:
        results = list(csv.DictReader(stream))
    if len(results) != EXPECTED_CASES:
        raise ValueError("case CSV row count mismatch")
    if len({row.get("case_id") for row in results}) != EXPECTED_CASES:
        raise ValueError("case IDs are not unique")
    if any(row.get("status") != "ok" for row in results):
        raise ValueError("case CSV contains a non-ok case")
    return receipt, manifest, results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--source-run-id", required=True)
    args = parser.parse_args()

    source = args.source_dir.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)

    preflight_root = source / "r15-full-cohort-preflight"
    tov_root = source / "r15-full-cohort-tov"
    preflight = verify_preflight(preflight_root)
    receipt, manifest, results = verify_tov(tov_root)

    backgrounds = tov_root / "backgrounds"
    profiles = manifest.get("background_profiles", [])
    if len(profiles) != EXPECTED_CASES:
        raise ValueError("background profile manifest length mismatch")

    profile_records = []
    for item in profiles:
        relative = Path(item["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"unsafe profile path: {relative}")
        path = tov_root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = sha256(path)
        if actual != item["sha256"]:
            raise ValueError(f"profile SHA-256 mismatch: {relative}")
        profile_records.append({
            "path": relative.as_posix(),
            "sha256": actual,
            "row_count": int(item["row_count"]),
            "bytes": path.stat().st_size,
        })

    if len({record["path"] for record in profile_records}) != EXPECTED_CASES:
        raise ValueError("background profile paths are not unique")

    core = output / "core"
    core.mkdir()
    shutil.copytree(preflight_root, core / "r15-full-cohort-preflight")
    tov_core = core / "r15-full-cohort-tov"
    tov_core.mkdir()
    for name in (
        "native_outer_crust_r15.json",
        "native_outer_crust_r15.csv",
        "native_outer_crust_r15_pair_cutoff_summary.csv",
        "manifest.json",
    ):
        shutil.copy2(tov_root / name, tov_core / name)

    profile_records.sort(key=lambda item: item["path"])
    chunks = [
        profile_records[i::SHARD_COUNT]
        for i in range(SHARD_COUNT)
    ]
    shard_index = []
    for index, chunk in enumerate(chunks, start=1):
        shard_root = output / f"profiles_part_{index:02d}"
        (shard_root / "backgrounds").mkdir(parents=True)
        for record in chunk:
            shutil.copy2(
                tov_root / record["path"],
                shard_root / record["path"],
            )
        chunk_manifest = {
            "schema": "ns-cowling-r15-full-cohort-profile-shard-v1",
            "source_run_id": args.source_run_id,
            "shard_index": index,
            "shard_count": SHARD_COUNT,
            "profile_count": len(chunk),
            "profiles": chunk,
        }
        manifest_path = shard_root / "shard_manifest.json"
        manifest_path.write_text(
            json.dumps(chunk_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        shard_index.append({
            "artifact_name": f"ns-cowling-r15-full-cohort-profiles-{index:02d}-{args.source_run_id}",
            "profile_count": len(chunk),
            "profile_bytes": sum(x["bytes"] for x in chunk),
            "shard_manifest_sha256": sha256(manifest_path),
        })

    index = {
        "schema": "ns-cowling-r15-full-cohort-repacked-index-v1",
        "source_run_id": args.source_run_id,
        "source_artifact": f"ns-cowling-r15-full-cohort-{args.source_run_id}",
        "expected_tov_cases": EXPECTED_CASES,
        "verified_tov_cases": len(results),
        "verified_projected_tables": len(receipt["source_table_receipts"]),
        "verified_background_profiles": len(profile_records),
        "source_summary_sha256": receipt["source_summary_sha256"],
        "tov_result_status": receipt["status"],
        "preflight_status": preflight["receipt"]["status"],
        "profile_shards": shard_index,
    }
    (output / "index.json").write_text(
        json.dumps(index, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(index, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
