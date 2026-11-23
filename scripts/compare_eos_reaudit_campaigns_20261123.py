#!/usr/bin/env python3
"""Compare repeated EOS campaign receipts and exact numerical file bytes."""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import tarfile

STABLE_FIELDS = (
    "status", "private", "crust_source_id", "core_source_id", "crust_model",
    "core_model", "method", "native_candidate_attempts", "historical_output_reused",
    "python_physical_prefilter", "source_coverage", "crust_view_id", "core_view_id",
    "verification", "native_report", "source_view_bounds", "failure_categories",
    "prior_or_total_failure_categories",
)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_campaign(root):
    summary_raw = (root / "SUMMARY.json").read_bytes()
    journal_raw = (root / "RESULTS.jsonl").read_bytes()
    summary = json.loads(summary_raw)
    if summary.get("status") != "COMPLETED" or summary.get("input_and_binary_hashes_unchanged") is not True:
        raise ValueError("comparison requires completed, immutable campaigns")
    records = {}
    for line in journal_raw.splitlines():
        record = json.loads(line)
        key = (record["crust_source_id"], record["core_source_id"], record["method"])
        if key in records:
            raise ValueError("duplicate source-pair/method")
        records[key] = record
    if len(records) != summary["requested_cells"]:
        raise ValueError("incomplete campaign")
    return summary, records, summary_raw, journal_raw


def numerical_files(record):
    artifact = record["artifact"]
    path = Path(artifact["path"])
    if sha256(path) != artifact["sha256"]:
        raise ValueError("archive identity changed")
    hashes = {}
    with tarfile.open(path, "r:gz") as archive:
        for member in archive:
            if member.isfile() and member.name.startswith(("internal/", "compose/")):
                if member.name in hashes:
                    raise ValueError("duplicate numerical archive member")
                stream = archive.extractfile(member)
                h = hashlib.sha256()
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    h.update(block)
                hashes[member.name] = h.hexdigest()
    if "internal/eos_data.csv" not in hashes or "internal/eos_data.nml" not in hashes:
        raise ValueError("missing native numerical tables")
    return hashes


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--previous", type=Path, required=True)
    p.add_argument("--expanded", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    previous, old, old_summary, old_journal = load_campaign(args.previous)
    expanded, new, new_summary, new_journal = load_campaign(args.expanded)
    missing = set(old) - set(new)
    if missing:
        raise ValueError("expanded campaign omitted previous cells")
    if previous["settings"]["binary_sha256"] != expanded["settings"]["binary_sha256"]:
        raise ValueError("native numerical binary differs")
    disagreements = []
    file_count = survivor_count = 0
    mapping = []
    for key in sorted(old):
        left, right = old[key], new[key]
        fields = [field for field in STABLE_FIELDS if left.get(field) != right.get(field)]
        if fields:
            disagreements.append(dict(source_pair_method=key, fields=fields))
        if left["status"] == right["status"] == "SURVIVED":
            lfiles, rfiles = numerical_files(left), numerical_files(right)
            if lfiles != rfiles:
                disagreements.append(dict(source_pair_method=key, numerical_file_mismatch=True))
            file_count += len(lfiles)
            survivor_count += 1
        mapping.append(dict(crust_source_id=key[0], core_source_id=key[1], method=key[2],
                            previous_cell_id=left["cell_id"], expanded_cell_id=right["cell_id"],
                            status=right["status"]))
    added = [new[k] for k in set(new) - set(old)]
    for root, summary_bytes, journal_bytes in ((args.previous, old_summary, old_journal),
                                               (args.expanded, new_summary, new_journal)):
        if ((root / "SUMMARY.json").read_bytes() != summary_bytes
                or (root / "RESULTS.jsonl").read_bytes() != journal_bytes):
            raise ValueError("campaign receipts changed during comparison")
    result = dict(
        schema="ns_cowling.reaudit_expansion_equivalence.v1",
        status="PASS" if not disagreements else "FAIL",
        previous_cells=len(old), expanded_cells=len(new), shared_cells=len(old),
        additional_cells=len(added), additional_statuses=dict(Counter(r["status"] for r in added)),
        shared_survivors=survivor_count, exact_numerical_files_compared=file_count,
        numerical_binary_sha256=expanded["settings"]["binary_sha256"],
        previous_journal_sha256=hashlib.sha256(old_journal).hexdigest(),
        expanded_journal_sha256=hashlib.sha256(new_journal).hexdigest(),
        previous_summary_sha256=hashlib.sha256(old_summary).hexdigest(),
        expanded_summary_sha256=hashlib.sha256(new_summary).hexdigest(),
        disagreement_count=len(disagreements), disagreements=disagreements,
        comparison_scope="same source-pair/method status, selected endpoints, native report, verification, failure categories, and all internal/CompOSE numerical file bytes; archive/provenance paths and elapsed times excluded",
    )
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out / "EQUIVALENCE.json").write_text(json.dumps(result, indent=2) + "\n")
    (args.out / "CELL_ID_MAPPING.json").write_text(json.dumps(mapping, indent=2) + "\n")
    print(json.dumps(result))
    if disagreements:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
