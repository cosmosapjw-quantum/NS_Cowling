#!/usr/bin/env python3
"""Publish lossless solid-XZ shards of complete public EOS candidate archives.

Original gzip archives and campaign inputs are read-only. Private candidates
remain metadata-only. This frontend performs compression, hashing and I/O only.
"""
from __future__ import annotations
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import sys
import tarfile
import tempfile
import time

from extract_eos_reaudit_candidate_20261123 import (
    CELL, file_sha256, publish_fresh_directory, validated_members,
)

SCHEMA = "ns_cowling.eos_delivery.v1"
INPUT_FILES = ("SUMMARY.json", "SURVIVORS.json", "RESULTS.jsonl", "SETTINGS.json", "INPUT_IDENTITIES.json")
META_KEYS = ("cell_id", "status", "private", "crust_source_id", "core_source_id",
             "method", "source_coverage", "crust_view_id", "core_view_id")


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _metadata(record):
    return {key: record[key] for key in META_KEYS if key in record}


def _validate_records(records):
    seen = set()
    for r in records:
        cell = r.get("cell_id", "")
        if not CELL.fullmatch(cell) or cell in seen:
            raise ValueError("invalid or duplicate campaign cell identifier")
        seen.add(cell)
        if type(r.get("private")) is not bool:
            raise ValueError("every cell requires an explicit private flag")
        if r.get("status") == "RUNNING":
            raise ValueError("nonterminal campaign cell")
        if any(str(r.get(k, "")).startswith("user_zip:") for k in ("crust_source_id", "core_source_id")) and not r["private"]:
            raise ValueError("uploaded user source was not marked private")


def _completed_campaign(campaign):
    raw = {name: (campaign / name).read_bytes() for name in INPUT_FILES}
    summary = json.loads(raw["SUMMARY.json"])
    if summary.get("status") != "COMPLETED" or summary.get("input_and_binary_hashes_unchanged") is not True:
        raise ValueError("campaign must be COMPLETED with unchanged input/binary identities")
    records = [json.loads(line) for line in raw["RESULTS.jsonl"].splitlines() if line.strip()]
    _validate_records(records)
    if len(records) != summary.get("recorded_cells") or len(records) != summary.get("requested_cells"):
        raise ValueError("completed campaign cell count mismatch")
    if dict(Counter(r["status"] for r in records)) != summary.get("statuses"):
        raise ValueError("completed campaign status census mismatch")
    survivors = json.loads(raw["SURVIVORS.json"])
    _validate_records(survivors)
    admitted = {r["cell_id"]: r for r in records if r["status"] == "SURVIVED"}
    if len(survivors) != summary.get("survivors") or {r["cell_id"]: r for r in survivors} != admitted:
        raise ValueError("survivor inventory differs from completed campaign journal")
    settings = json.loads(raw["SETTINGS.json"])
    if settings != summary.get("settings"):
        raise ValueError("campaign settings differ from completion receipt")
    return records, summary, {name: _sha(data) for name, data in raw.items()}


def _candidate_bytes(campaign, record):
    cell = record["cell_id"]
    if record.get("private") is not False or record["status"] != "SURVIVED":
        raise ValueError("only public survivors may enter shards")
    source = campaign / "artifacts" / (cell + ".tar.gz")
    if source.is_symlink() or not source.resolve(strict=True).is_relative_to(campaign):
        raise ValueError("public archive must be inside the campaign directory")
    compressed = source.read_bytes()
    expected = record.get("artifact", {})
    if _sha(compressed) != expected.get("sha256"):
        raise ValueError("original gzip hash mismatch: " + cell)
    if len(compressed) != expected.get("bytes"):
        raise ValueError("original gzip byte count mismatch: " + cell)
    raw = gzip.decompress(compressed)
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:") as archive:
        members = validated_members(archive)
        matches = [m for m, _ in members if m.name == "PROVENANCE.json" and m.isfile()]
        if len(matches) != 1:
            raise ValueError("candidate provenance is absent or ambiguous")
        provenance = json.load(archive.extractfile(matches[0]))
        if provenance.get("private") is not False:
            raise ValueError("private candidate provenance cannot enter a public shard")
        for key in ("cell_id", "crust_source_id", "core_source_id", "method"):
            if provenance.get(key) != record.get(key):
                raise ValueError("candidate provenance identity mismatch: " + key)
    return compressed, raw


def _build_shard(number, records, campaign, stage, preset):
    relative = f"shards/shard_{number:05d}.tar.xz"
    path = stage / relative
    cells = {}
    with tarfile.open(path, "w:xz", format=tarfile.USTAR_FORMAT, preset=preset) as archive:
        for record in records:
            compressed, raw = _candidate_bytes(campaign, record)
            cell = record["cell_id"]
            member = cell + ".tar"
            info = tarfile.TarInfo(member)
            info.size = len(raw); info.mode = 0o644
            info.uid = info.gid = 0; info.mtime = 0
            info.uname = info.gname = ""
            archive.addfile(info, io.BytesIO(raw))
            cells[cell] = {**_metadata(record), "shard": relative, "member": member,
                           "original_gzip_sha256": _sha(compressed), "original_gzip_bytes": len(compressed),
                           "raw_tar_sha256": _sha(raw), "raw_tar_bytes": len(raw)}
    return dict(path=relative, sha256=file_sha256(path), bytes=path.stat().st_size,
                cells=len(records)), cells


def _build_delivery(campaign, output, records, summary, hashes, *, shard_cells, workers, xz_preset,
                    benchmark=False, progress=False):
    if shard_cells not in (128, 256) or not 1 <= workers <= 32 or not 0 <= xz_preset <= 9:
        raise ValueError("use 128/256 cells per shard, 1..32 workers, and XZ preset 0..9")
    output = Path(output).expanduser().absolute()
    if os.path.lexists(output):
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".eos-delivery-", dir=output.parent))
    start = time.perf_counter()
    public = sorted((r for r in records if r["status"] == "SURVIVED" and not r["private"]), key=lambda r: r["cell_id"])
    private = {r["cell_id"]: {**_metadata(r), "publication": "METADATA_ONLY_PRIVATE_BACKUP"}
               for r in records if r["private"]}
    try:
        (stage / "shards").mkdir()
        chunks = [public[i:i+shard_cells] for i in range(0, len(public), shard_cells)]
        shards, cells = [], {}
        with ThreadPoolExecutor(max_workers=workers) as pool:
            jobs = [pool.submit(_build_shard, i, group, campaign, stage, xz_preset)
                    for i, group in enumerate(chunks)]
            for job in as_completed(jobs):
                shard, group = job.result()
                shards.append(shard); cells.update(group)
                if progress:
                    print(json.dumps(dict(event="SHARD_COMPLETE", shards=len(shards), total=len(chunks),
                                          public_cells=len(cells))), file=sys.stderr, flush=True)
        if not benchmark:
            current = {name: file_sha256(campaign / name) for name in INPUT_FILES}
            if current != hashes:
                raise ValueError("campaign JSON/journal changed during packaging")
        original = sum(r["original_gzip_bytes"] for r in cells.values())
        published = sum(s["bytes"] for s in shards)
        index = dict(schema=SCHEMA,
            status="BENCHMARK_ONLY" if benchmark else "COMPLETED_LOSSLESS_PUBLIC_DELIVERY",
            campaign_file_sha256=hashes,
            source_catalog_sha256=summary.get("settings", {}).get("catalog_sha256"),
            native_binary_sha256=summary.get("settings", {}).get("binary_sha256"),
            campaign_requested_cells=summary.get("requested_cells"),
            campaign_recorded_cells=summary.get("recorded_cells"),
            public_survivors=len(cells), private_survivors=sum(r["status"] == "SURVIVED" for r in private.values()),
            private_cells_count=len(private),
            compression=dict(format="solid_tar_xz_of_original_tar_bytes", cells_per_shard=shard_cells,
                             xz_preset=xz_preset, workers=workers),
            original_public_gzip_bytes=original, published_shard_bytes=published,
            gzip_to_shard_ratio=original / published if published else None,
            saved_fraction=1 - published / original if original else None,
            packaging_seconds=time.perf_counter()-start,
            integrity_scope="each member equals gzip.decompress(original candidate archive) byte-for-byte",
            scientific_scope="complete original candidate content preserved; no EOS recalculation or new physical admission",
            private_policy="private numeric files and traces are absent; metadata only",
            shards=sorted(shards, key=lambda r: r["path"]),
            cells={k: cells[k] for k in sorted(cells)},
            private_cells={k: private[k] for k in sorted(private)})
        index_name = "BENCHMARK_INDEX.json" if benchmark else "DELIVERY_INDEX.json"
        (stage / index_name).write_text(json.dumps(index, indent=2, ensure_ascii=False, allow_nan=False)+"\n")
        publish_fresh_directory(stage, output)
        return index
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def package_campaign(campaign, output, *, shard_cells=128, workers=4, xz_preset=6, progress=False):
    campaign = Path(campaign).resolve(strict=True)
    if os.path.lexists(output):
        raise FileExistsError(output)
    records, summary, hashes = _completed_campaign(campaign)
    return _build_delivery(campaign, output, records, summary, hashes, shard_cells=shard_cells,
                           workers=workers, xz_preset=xz_preset, progress=progress)


def benchmark_campaign(campaign, output, *, limit=128, shard_cells=128, workers=1, xz_preset=6, progress=False):
    output = Path(output).resolve()
    if not output.is_relative_to(Path("/tmp").resolve()) or limit <= 0:
        raise ValueError("benchmark output must be under /tmp and limit must be positive")
    campaign = Path(campaign).resolve(strict=True)
    snapshot = (campaign / "RESULTS.jsonl").read_bytes()
    complete = snapshot[:snapshot.rfind(b"\n")+1]
    records = [json.loads(line) for line in complete.splitlines() if line.strip()]
    _validate_records(records)
    public = [r for r in records if r["status"] == "SURVIVED" and not r["private"]][:limit]
    if not public:
        raise ValueError("no completed public candidates in benchmark snapshot")
    summary = dict(recorded_cells=len(records), settings=json.loads((campaign / "SETTINGS.json").read_bytes()))
    return _build_delivery(campaign, output, public, summary, {"RESULTS.snapshot": _sha(complete)},
                           shard_cells=shard_cells, workers=workers, xz_preset=xz_preset,
                           benchmark=True, progress=progress)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("package", "benchmark"))
    parser.add_argument("--campaign", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--shard-cells", type=int, choices=(128, 256), default=128)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--xz-preset", type=int, default=6)
    parser.add_argument("--limit", type=int, default=128, help="benchmark only")
    args = parser.parse_args()
    kwargs = dict(shard_cells=args.shard_cells, workers=args.workers, xz_preset=args.xz_preset, progress=True)
    if args.mode == "benchmark":
        result = benchmark_campaign(args.campaign, args.out, limit=args.limit, **kwargs)
    else:
        result = package_campaign(args.campaign, args.out, **kwargs)
    print(json.dumps({k: result[k] for k in ("status", "public_survivors", "private_survivors",
                     "original_public_gzip_bytes", "published_shard_bytes", "gzip_to_shard_ratio", "packaging_seconds")}, sort_keys=True))


if __name__ == "__main__":
    main()
