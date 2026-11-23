#!/usr/bin/env python3
"""Verify a public delivery shard and restore one complete original candidate.

This frontend transports bytes only. It performs no EOS numerical calculation.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import re
import shutil
import tarfile
import tempfile

CELL = re.compile(r"cell_[0-9]{6}\Z")


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_parts(name):
    if not isinstance(name, str) or not name or "\\" in name or "\x00" in name:
        raise ValueError("unsafe archive path")
    parts = name.rstrip("/").split("/")
    if (PurePosixPath(name).is_absolute() or PureWindowsPath(name).drive
            or any(p in ("", ".", "..") for p in parts)):
        raise ValueError("unsafe archive path: " + name)
    return tuple(parts)


def validated_members(archive):
    """Allow only regular files/directories, without links or path aliases."""
    result, seen = [], {}
    for member in archive.getmembers():
        parts = safe_parts(member.name)
        if not (member.isfile() or member.isdir()) or member.issym() or member.islnk():
            raise ValueError("unsafe archive member type: " + member.name)
        if member.sparse is not None:
            raise ValueError("sparse archive members are unsupported")
        if parts in seen:
            raise ValueError("duplicate archive path: " + member.name)
        seen[parts] = member
        result.append((member, parts))
    for _, parts in result:
        for length in range(1, len(parts)):
            parent = seen.get(parts[:length])
            if parent is not None and not parent.isdir():
                raise ValueError("archive file is also a parent directory")
    return result


def publish_fresh_directory(stage, output):
    """Reserve a fresh target atomically, then replace only our empty reservation."""
    output.mkdir(exist_ok=False)
    try:
        os.replace(stage, output)
    except BaseException:
        try:
            output.rmdir()
        except OSError:
            pass
        raise


def extract_candidate(index_path, cell_id, output):
    if not CELL.fullmatch(cell_id):
        raise ValueError("invalid cell identifier")
    index_path = Path(index_path).resolve(strict=True)
    base = index_path.parent
    output = Path(output).expanduser().absolute()
    if os.path.lexists(output):
        raise FileExistsError(output)
    index = json.loads(index_path.read_bytes())
    if index.get("schema") != "ns_cowling.eos_delivery.v1":
        raise ValueError("unsupported delivery schema")
    if index.get("status") != "COMPLETED_LOSSLESS_PUBLIC_DELIVERY":
        raise ValueError("extraction requires a completed public delivery")
    if cell_id in index.get("private_cells", {}):
        raise ValueError("private candidate is not included in public delivery")
    entry = index.get("cells", {}).get(cell_id)
    if not isinstance(entry, dict) or entry.get("private") is not False:
        raise ValueError("public candidate not found")
    parts = safe_parts(entry["shard"])
    shard_path = base.joinpath(*parts)
    if shard_path.is_symlink() or not shard_path.resolve(strict=True).is_relative_to(base):
        raise ValueError("shard path leaves the delivery directory")
    shards = [s for s in index["shards"] if s["path"] == entry["shard"]]
    if len(shards) != 1:
        raise ValueError("ambiguous shard identity")
    shard = shards[0]
    actual_shard_hash = file_sha256(shard_path)
    if actual_shard_hash != shard["sha256"]:
        raise ValueError("shard hash mismatch")
    if shard_path.stat().st_size != shard["bytes"]:
        raise ValueError("shard byte count mismatch")
    expected_member = cell_id + ".tar"
    if entry["member"] != expected_member:
        raise ValueError("candidate member does not match cell identity")
    with tarfile.open(shard_path, "r:xz") as archive:
        members = validated_members(archive)
        if any(not m.isfile() or len(p) != 1 for m, p in members):
            raise ValueError("shard must contain only regular top-level candidate tars")
        matches = [m for m, _ in members if m.name == expected_member]
        if len(matches) != 1:
            raise ValueError("candidate tar absent or duplicated")
        member = matches[0]
        if member.size != entry["raw_tar_bytes"]:
            raise ValueError("candidate tar byte count mismatch")
        with archive.extractfile(member) as stream:
            raw = stream.read()
    raw_hash = hashlib.sha256(raw).hexdigest()
    if raw_hash != entry["raw_tar_sha256"]:
        raise ValueError("candidate tar hash mismatch")
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:") as archive:
        members = validated_members(archive)
        output.parent.mkdir(parents=True, exist_ok=True)
        stage = Path(tempfile.mkdtemp(prefix=".eos-extract-", dir=output.parent))
        try:
            count = 0
            for member, parts in members:
                target = stage.joinpath(*parts)
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                with archive.extractfile(member) as source, target.open("xb") as dest:
                    shutil.copyfileobj(source, dest)
                if target.stat().st_size != member.size:
                    raise ValueError("extracted file byte count mismatch")
                os.chmod(target, member.mode & 0o777)
                count += 1
            publish_fresh_directory(stage, output)
        finally:
            if stage.exists():
                shutil.rmtree(stage)
    return dict(status="EXTRACTED_HASH_VERIFIED", cell_id=cell_id,
                shard_sha256=actual_shard_hash, raw_tar_sha256=raw_hash,
                files=count, output=str(output),
                scope="complete original candidate file bytes; no physical recalculation")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", required=True, type=Path)
    parser.add_argument("--cell", required=True)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    result = extract_candidate(args.index, args.cell, args.out)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
