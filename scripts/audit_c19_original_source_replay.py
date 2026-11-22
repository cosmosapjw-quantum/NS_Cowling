#!/usr/bin/env python3
"""Frontend-only C17/C19 original-source semantic replay comparison.

All physical EOS/TOV, Cowling shooting, Maxwell matching, overlap arithmetic
and branch admission run in compiled Fortran. This parser only asserts output
identities and compares already-computed numerical receipts.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BASELINE = ROOT / "data/cowling_c17_extended_branch_20261122/REAL_P2_INITIAL"
SOURCE = ROOT / "data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SOURCE_SHA256 = "299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
OVERLAP_FIELDS = {"absolute_radial_W_overlap", "high_phase_overlap", "low_phase_overlap"}
OVERLAP_ABS_TOL = 5e-13  # Diagnostic roundoff only; NOT a physical mode-gate change.


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def compare_science(a, b, path: str, found: list[dict], errors: list[str]) -> None:
    """All source-root observables are exactly equal. Only overlap ULP drift allowed."""
    if type(a) is not type(b):
        errors.append(f"{path}: value types differ")
        return
    if isinstance(a, dict):
        if set(a) != set(b):
            errors.append(f"{path}: key sets differ")
        for key in a.keys() & b.keys():
            compare_science(a[key], b[key], f"{path}/{key}", found, errors)
    elif isinstance(a, list):
        if len(a) != len(b):
            errors.append(f"{path}: array lengths differ")
        for i, (x, y) in enumerate(zip(a, b)):
            compare_science(x, y, f"{path}[{i}]", found, errors)
    elif a != b:
        field = path.rsplit("/", 1)[-1]
        if field in OVERLAP_FIELDS and isinstance(a, float):
            difference = abs(a - b)
            found.append(dict(field=field, path=path, absolute_difference=difference))
            if difference > OVERLAP_ABS_TOL:
                errors.append(f"{path}: overlap delta {difference:.17g} > {OVERLAP_ABS_TOL:.1e}")
        else:
            errors.append(f"{path}: non-overlap science changed")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    p.add_argument("--replay", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        p.error("output must be a new file to preserve previous receipts")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    overlaps: list[dict] = []
    errors: list[str] = []
    archive_sha = sha(SOURCE)
    if archive_sha != SOURCE_SHA256:
        errors.append("original locked EOS ZIP SHA256 mismatch")
    record = dict(
        schema="ns_cowling.c19.c17_source_exact_semantic_replay.v1",
        source_sha256=archive_sha,
        baseline=str(args.baseline),
        replay=str(args.replay),
        overlap_absolute_roundoff_tolerance=OVERLAP_ABS_TOL,
        models={},
    )
    for model in ("M2", "M3"):
        a_path = args.baseline / f"{model}.json"
        b_path = args.replay / f"{model}.json"
        a = json.loads(a_path.read_text())
        b = json.loads(b_path.read_text())
        compare_science(a, b, model, overlaps, errors)
        if len(a["densities"]) != 20 or len(b["densities"]) != 20:
            errors.append(f"{model}: expected exact 20-density source cohort")
        original_scans = args.baseline / f"{model}_all_scan_points.json.gz"
        replay_scans = args.replay / f"{model}_all_scan_points.json.gz"
        ah, bh = sha(original_scans), sha(replay_scans)
        if ah != bh:
            errors.append(f"{model}: original Fortran determinant scan bytes changed")
        statuses = {}
        for kind in ("f", "p1"):
            a_status = [r["local_branches"][kind]["status"] for r in a["densities"]]
            b_status = [r["local_branches"][kind]["status"] for r in b["densities"]]
            if a_status != b_status:
                errors.append(f"{model}: {kind} branch admission statuses changed")
            statuses[kind] = dict(sorted(Counter(b_status).items()))
        record["models"][model] = dict(
            count_densities=len(b["densities"]),
            c17_json_sha256=sha(a_path),
            c19_json_sha256=sha(b_path),
            scan_bytes_identical=ah == bh,
            original_scans_sha256=ah,
            branch_status_counts=statuses,
        )
    record["overlap_values_differing_in_last_digits"] = len(overlaps)
    record["max_absolute_overlap_change"] = (
        max((x["absolute_difference"] for x in overlaps), default=0.0)
    )
    record["nonoverlap_science_identical"] = not errors
    record["errors"] = errors
    record["status"] = "PASS_SCOPED_C17_C19_REPLAY" if not errors else "FAIL_REPLAY_MISMATCH"
    args.output.write_text(json.dumps(record, sort_keys=True, indent=2) + "\n")
    print("C19_SOURCE_REPLAY", record["status"])
    print("C19_MAX_ABSOLUTE_OVERLAP_CHANGE", record["max_absolute_overlap_change"])
    print("C19_CHANGED_OVERLAP_VALUES", record["overlap_values_differing_in_last_digits"])
    for model, item in record["models"].items():
        print("C19_MODEL", model, "STARS", item["count_densities"],
              "SCANS_IDENTICAL", item["scan_bytes_identical"],
              "ADMISSION_COUNTS", item["branch_status_counts"])
    if errors:
        print("C19_FAILURES", errors[:12])
        raise SystemExit(1)


if __name__ == "__main__":
    main()
