#!/usr/bin/env python3
"""Bundle original public failure receipts without changing their gzip bytes."""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import tarfile


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--campaign", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    root = args.campaign.resolve(strict=True)
    summary_raw = (root / "SUMMARY.json").read_bytes()
    summary = json.loads(summary_raw)
    if summary.get("status") != "COMPLETED" or summary.get("input_and_binary_hashes_unchanged") is not True:
        raise ValueError("a completed immutable campaign is required")
    journal = (root / "RESULTS.jsonl").read_bytes()
    records = [json.loads(line) for line in journal.splitlines()]
    if len(records) != summary["requested_cells"] or len(records) != summary["recorded_cells"]:
        raise ValueError("campaign count differs from completion receipt")
    if len({r["cell_id"] for r in records}) != len(records):
        raise ValueError("duplicate campaign cell identity")
    if dict(Counter(r["status"] for r in records)) != summary["statuses"]:
        raise ValueError("campaign outcome census differs from completion receipt")
    for record in records:
        if type(record.get("private")) is not bool:
            raise ValueError("explicit private flag is required")
        if any(record[k].startswith("user_zip:") for k in ("crust_source_id", "core_source_id")) and not record["private"]:
            raise ValueError("uploaded user source was not marked private")
    expected = {r["cell_id"]: r for r in records if not r["private"] and "failure_trace" in r}
    files = {p.name.removesuffix(".json.gz") for p in (root / "failure_traces").glob("*.json.gz")}
    if set(expected) != files:
        raise ValueError("public trace files differ from the completed journal")
    if args.out.exists():
        raise FileExistsError(args.out)
    manifest = {}
    with tarfile.open(args.out, "x:xz", format=tarfile.USTAR_FORMAT, preset=3) as archive:
        for cell, record in sorted(expected.items()):
            path = root / "failure_traces" / (cell + ".json.gz")
            if path.is_symlink() or not path.resolve(strict=True).is_relative_to(root):
                raise ValueError("trace leaves campaign directory")
            raw = path.read_bytes()
            trace = record["failure_trace"]
            if path.resolve() != Path(trace["path"]).resolve() or digest(raw) != trace["sha256"]:
                raise ValueError("public trace identity mismatch")
            name = "failure_traces/" + path.name
            manifest[name] = dict(sha256=digest(raw), bytes=len(raw))
            info = tarfile.TarInfo(name)
            info.size = len(raw)
            info.mode = 0o644
            info.mtime = 0
            archive.addfile(info, io.BytesIO(raw))
    verified = {}
    with tarfile.open(args.out, "r:xz") as archive:
        for member in archive:
            if not member.isfile() or member.name not in manifest or member.name in verified:
                raise ValueError("unexpected archive member")
            raw = archive.extractfile(member).read()
            verified[member.name] = dict(sha256=digest(raw), bytes=len(raw))
    if (verified != manifest or (root / "RESULTS.jsonl").read_bytes() != journal
            or (root / "SUMMARY.json").read_bytes() != summary_raw):
        raise ValueError("archive verification or campaign identity failed")
    receipt = dict(
        schema="ns_cowling.public_failure_receipt_bundle.v1",
        status="ALL_ORIGINAL_GZIP_BYTES_REVERIFIED", public_trace_files=len(manifest),
        original_gzip_bytes=sum(v["bytes"] for v in manifest.values()),
        archive_name=args.out.name, archive_sha256=digest(args.out.read_bytes()),
        archive_bytes=args.out.stat().st_size, journal_sha256=digest(journal),
        summary_sha256=digest(summary_raw),
        private_numeric_traces_included=False,
        extraction="Extract with tar -xJf into a fresh campaign copy; original failure_traces/ paths are preserved.",
        file_sha256_and_sizes=manifest,
    )
    receipt_path = args.out.with_name("FAILURE_TRACES_RECEIPT.json")
    with receipt_path.open("x") as stream:
        stream.write(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({k:v for k,v in receipt.items() if k != "file_sha256_and_sizes"}))


if __name__ == "__main__":
    main()
